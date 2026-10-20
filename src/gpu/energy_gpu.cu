// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "energy_gpu.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <type_traits>
#include <vector>

#include "error.hpp"

namespace od {
namespace {

constexpr int kMaxThreads = 256;

void check(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw Error(std::string("CUDA: ") + what + ": " + cudaGetErrorString(e));
}

template <class T>
struct DevBuf {
  T* p = nullptr;
  std::size_t n = 0;
  DevBuf() = default;
  DevBuf(const DevBuf&) = delete;
  DevBuf& operator=(const DevBuf&) = delete;
  ~DevBuf() {
    if (p) cudaFree(p);
  }
  void reserve(std::size_t count) {
    if (count <= n) return;
    if (p) cudaFree(p);
    p = nullptr;
    n = 0;
    check(cudaMalloc(&p, count * sizeof(T)), "cudaMalloc");
    n = count;
  }
  void release() {
    if (p) cudaFree(p);
    p = nullptr;
    n = 0;
  }
  void upload(const std::vector<T>& v) {
    reserve(v.empty() ? 1 : v.size());
    if (!v.empty()) check(cudaMemcpy(p, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice), "upload");
  }
};

__device__ inline double warpSum(double v) {
  for (int o = 16; o > 0; o /= 2) v += __shfl_down_sync(0xffffffffu, v, o);
  return v;
}
__device__ inline double warpMax(double v) {
  for (int o = 16; o > 0; o /= 2) v = fmax(v, __shfl_down_sync(0xffffffffu, v, o));
  return v;
}

// Pair data of the kernel: pairs in the order of EnergyModel::pairs (tf: type << 2 | hbond << 1 |
// heavy), per-atom constants in T.
struct GpuPair {
  int i, j, tf;
};
template <class T>
struct GpuAtom {
  T cq, q, sv, vol;  // elecCoeff·q, q, desolvCoeff·S, V
};
template <class T>
struct GpuPairParams {
  T AC, BD, reqm;
  int m, n;
};
template <class T>
struct PairData {
  const GpuPair* pairs;
  const GpuAtom<T>* atoms;
  const GpuPairParams<T>* params;
  int count;
};

// One atom pair (the terms of flat::pairTerm with one reciprocal per function instead of
// divisions): adds to e[0..4] and clashPairs; d = xi − xj, f = dE/dr / r (0: no gradient).
// Returns false beyond the cutoff.
template <class T>
__device__ inline bool pairTermGpu(const EnergyModelView& m, const GpuPair& pr, const PairData<T>& pd, const T* x,
                                   double e[5], int& clashPairs, T d[3], T& f) {
  for (int a = 0; a < 3; ++a) d[a] = x[3 * pr.i + a] - x[3 * pr.j + a];
  const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (r2 >= static_cast<T>(kElecCutoff * kElecCutoff)) return false;
  const T r = sqrt(r2);
  T dedr = 0;
  if (r < static_cast<T>(kVdwCutoff)) {
    const GpuPairParams<T> q = pd.params[pr.tf >> 2];
    T ee, de;
    if (q.m == 12 && (q.n == 6 || q.n == 10)) {
      const T half = static_cast<T>(m.halfSmooth);
      T s = q.reqm, ds = 0;
      if (r <= q.reqm - half) {
        s = r + half;
        ds = 1;
      } else if (r >= q.reqm + half) {
        s = r - half;
        ds = 1;
      }
      const T inv = T(1) / s, inv2 = inv * inv, inv6 = inv2 * inv2 * inv2;
      const T rm = inv6 * inv6, rn = q.n == 6 ? inv6 : inv6 * inv2 * inv2;
      ee = q.AC * rm - q.BD * rn;
      de = ds * (static_cast<T>(q.n) * q.BD * rn - T(12) * q.AC * rm) * inv;
    } else {
      ad4VdwHb<T>(r, q.AC, q.BD, q.m, q.n, q.reqm, static_cast<T>(m.halfSmooth), ee, de);
    }
    e[pr.tf & 2 ? 1 : 0] += ee;
    dedr += de;
  }
  // Electrostatics and desolvation (ad4Elec, ad4Desolv): evaluated at rMin below rMin.
  const T rMin = static_cast<T>(m.rMin);
  const bool clamp = r < rMin;
  const T rr = clamp ? rMin : r;
  const T irr = T(1) / rr;
  const T sh = rr + T(1.26366), ish = T(1) / sh, ish2 = ish * ish;
  const T diel = T(1.10859) * ish2 + T(0.010358);
  const GpuAtom<T> ai = pd.atoms[pr.i], aj = pd.atoms[pr.j];
  const T qq = ai.cq * aj.q, sv = ai.sv * aj.vol + aj.sv * ai.vol;
  const T es = qq * irr;
  e[2] += es * diel;
  const T rr2 = rr * rr;
  const T num = T(12.96) - T(0.1063) * rr2 * (T(1) - T(0.001947) * rr2);
  const T den = T(12.96) + rr2 * (T(0.4137) + rr2 * (T(0.00357) + T(0.000112) * rr2));
  const T iden = T(1) / den;
  const T sn = sv * num * iden;
  e[3] += sn;
  if (!clamp) {
    const T dnum = -T(0.2126) * rr + T(4) * T(0.1063) * T(0.001947) * rr2 * rr;
    const T dden = rr * (T(2) * T(0.4137) + rr2 * (T(4) * T(0.00357) + T(6) * T(0.000112) * rr2));
    dedr += -es * irr * diel - es * T(2.21718) * ish2 * ish + (sv * dnum - sn * dden) * iden;
  }
  const T cd = static_cast<T>(m.clashDistance);
  if ((pr.tf & 1) && r2 < cd * cd) {
    T ee, de;
    clashPenalty<T>(r, cd, static_cast<T>(m.clashBase), static_cast<T>(m.clashK), ee, de);
    e[4] += ee;
    ++clashPairs;
    dedr += de;
  }
  f = r > 0 ? dedr * (clamp ? T(1) / r : irr) : T(0);
  return true;
}

// 9b probe atoms of one pose against field f (inlined: f may be a kernel parameter); adds to
// v[8..16] (inter, interClash, outside, translation and rotation gradient) and g.
template <int kThreads>
__device__ __forceinline__ void probeTerms(const DockView& dv, const FieldView& f, const float* x, const DockPose& ps,
                                           double* g, int cc, double* v) {
  for (int k = threadIdx.x; k < dv.probeCount; k += kThreads) {
    float gw[3], r[3], gy[3];
    if (!dock::probeTerm(dv, f, k, x, ps, v[8], v[9], g ? gw : nullptr, r)) {
      v[10] += 1;
      continue;
    }
    if (!g) continue;
    dock::probeGrad(ps, gw, r, v + 11, v + 14, gy);
    const int c = dv.probes[k].coord;
    for (int i = 0; i < 3; ++i) atomicAdd(g + i * cc + c, static_cast<double>(gy[i]));
  }
}

// Torsion transforms (spec 15.1.15). When the axis atoms of each torsion are moved only by
// torsions applied after it (lower index) and the torsions moving each coordinate are the chain
// of the last of them and its parents, coordinate a ends at M_deep[a] (x_ref) with M_t =
// M_parent[t] ∘ A_t, A_t the rotation of torsion t about its input axis; the products along the
// chains are formed in `rounds` pointer-jumping rounds. M_t is kept as R (3 × 3) and y:
// M_t(x) = y + R (x − o_t), o_t the input axis point of torsion t (local vectors only). Work memory: two sets of 12 T per torsion
// and two sets of parent indices, after the coordinates and torsion cos and sin in shared memory
// (shared) or per block in buf.
template <class T>
struct TorsionData {
  const T* axis;  // per torsion: o, unit k (k = 0: zero-length axis, A_t = identity)
  const int* parent;
  const int* deep;  // per coordinate; −1: not moved
  char* buf;
  int rounds, ok, shared;
};

// Blocks of NT threads; each block evaluates samples blockIdx.x, blockIdx.x + gridDim.x, … of the
// n samples. shared: the coordinates and the cos and sin of each torsion are kept in dynamic
// shared memory (the gradient stays in global memory: double atomics on shared memory are
// emulated on sm_70); xs is then only used to return the coordinates (keepX). The gradient (per
// block in gs) is stored by component (x of all coordinates, then y, then z): the atomics of a
// warp on consecutive atoms then share memory sectors. Tor: torsion transforms (TorsionData);
// otherwise the torsions are applied one by one. Dock: 0 energies only, 1 9b poses against
// dv.field, 2 9b poses against dv.fields[DockPose::field].
template <class T, int Dock, int NT, bool Tor>
__global__ void __launch_bounds__(NT) evaluateKernel(EnergyModelView m, PairData<T> pd, TorsionData<T> td,
                                                     const T* angles, T* xs, double* gs, FlatTerms* terms,
                                                     double* grad, DockView dv, const DockPose* poses, DockTerms* dterms,
                                                     double* grad6, int n, int shared, int keepX) {
  constexpr int kThreads = NT, kWarps = NT / 32;
  extern __shared__ __align__(16) unsigned char smem[];
  const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount);
  const int cc = m.coordCount;
  for (int s = blockIdx.x; s < n; s += gridDim.x) {
    const T* ang = angles + static_cast<std::size_t>(s) * m.torsionCount;
    T* x;
    double* g;
    x = shared ? reinterpret_cast<T*>(smem) : xs + s * nc;
    g = grad ? gs + blockIdx.x * nc : nullptr;

    T* tcs = shared ? x + nc : nullptr;
    for (std::size_t i = tid; i < nc; i += kThreads) {
      x[i] = static_cast<T>(m.ref[i]);
      if (g) g[i] = 0;
    }
    if (tcs)
      for (int t = tid; t < m.torsionCount; t += kThreads) flat::torsionCosSin(m, t, ang, tcs[2 * t], tcs[2 * t + 1]);
    __syncthreads();

    if constexpr (Tor) {
      const int nt = m.torsionCount;
      T* mb[2];
      int* pb[2];
      char* b = td.shared ? reinterpret_cast<char*>(x + nc + 2 * nt)
                          : td.buf + static_cast<std::size_t>(blockIdx.x) * (24 * sizeof(T) + 2 * sizeof(int)) * nt;
      mb[0] = reinterpret_cast<T*>(b);
      mb[1] = mb[0] + 12 * nt;
      pb[0] = reinterpret_cast<int*>(mb[1] + 12 * nt);
      pb[1] = pb[0] + nt;
      for (int u = tid; u < nt; u += kThreads) {
        T c, sn;
        if (tcs) {
          c = tcs[2 * u];
          sn = tcs[2 * u + 1];
        } else {
          flat::torsionCosSin(m, u, ang, c, sn);
        }
        const T* a = td.axis + 6 * u;
        const T o[3] = {a[0], a[1], a[2]}, k[3] = {a[3], a[4], a[5]};
        if (k[0] == 0 && k[1] == 0 && k[2] == 0) {
          c = 1;
          sn = 0;
        }
        // R v = c v + s (k × v) + (1 − c) k (k · v).
        const T R[3][3] = {{c + (1 - c) * k[0] * k[0], (1 - c) * k[0] * k[1] - sn * k[2], (1 - c) * k[0] * k[2] + sn * k[1]},
                           {(1 - c) * k[1] * k[0] + sn * k[2], c + (1 - c) * k[1] * k[1], (1 - c) * k[1] * k[2] - sn * k[0]},
                           {(1 - c) * k[2] * k[0] - sn * k[1], (1 - c) * k[2] * k[1] + sn * k[0], c + (1 - c) * k[2] * k[2]}};
        T* M = mb[0] + 12 * u;
        for (int r = 0; r < 3; ++r) {
          for (int q = 0; q < 3; ++q) M[3 * r + q] = R[r][q];
          M[9 + r] = o[r];
        }
        pb[0][u] = td.parent[u];
      }
      __syncthreads();
      int cur = 0;
      for (int round = 0; round < td.rounds; ++round, cur ^= 1) {
        const T* src = mb[cur];
        T* dst = mb[cur ^ 1];
        for (int u = tid; u < nt; u += kThreads) {
          const int pa = pb[cur][u];
          const T* A = src + 12 * u;
          T* D = dst + 12 * u;
          if (pa < 0) {
            for (int e = 0; e < 12; ++e) D[e] = A[e];
            pb[cur ^ 1][u] = -1;
            continue;
          }
          const T* P = src + 12 * pa;
          const T* op = td.axis + 6 * pa;
          const T w[3] = {A[9] - op[0], A[10] - op[1], A[11] - op[2]};
          for (int r = 0; r < 3; ++r) {
            for (int q = 0; q < 3; ++q) D[3 * r + q] = P[3 * r] * A[q] + P[3 * r + 1] * A[3 + q] + P[3 * r + 2] * A[6 + q];
            D[9 + r] = P[9 + r] + (P[3 * r] * w[0] + P[3 * r + 1] * w[1] + P[3 * r + 2] * w[2]);
          }
          pb[cur ^ 1][u] = pb[cur][pa];
        }
        __syncthreads();
      }
      const T* mf = mb[cur];
      for (int i = tid; i < cc; i += kThreads) {
        const int dd = td.deep[i];
        if (dd < 0) continue;
        const T* A = mf + 12 * dd;
        const T* od = td.axis + 6 * dd;
        const T p0 = x[3 * i] - od[0], p1 = x[3 * i + 1] - od[1], p2 = x[3 * i + 2] - od[2];
        for (int r = 0; r < 3; ++r) x[3 * i + r] = A[9 + r] + (A[3 * r] * p0 + A[3 * r + 1] * p1 + A[3 * r + 2] * p2);
      }
      __syncthreads();
    } else {
      __shared__ T so[3], sk[3], scs[2];
      __shared__ int sok;
      for (int t = m.torsionCount; t-- > 0;) {
        if (tid == 0) {
          T o[3], k[3], c = 0, sn = 0;
          if (tcs) {
            sok = flat::torsionAxisOnly(m, t, x, o, k);
            c = tcs[2 * t];
            sn = tcs[2 * t + 1];
          } else {
            sok = flat::torsionAxis(m, t, ang, x, o, k, c, sn);
          }
          for (int a = 0; a < 3; ++a) {
            so[a] = o[a];
            sk[a] = k[a];
          }
          scs[0] = c;
          scs[1] = sn;
        }
        __syncthreads();
        if (sok) {
          const FlatTorsion& tt = m.torsions[t];
          for (int q = tt.movingBegin + tid; q < tt.movingEnd; q += kThreads)
            flat::rotateAtom(x + 3 * m.moving[q], so, sk, scs[0], scs[1]);
        }
        __syncthreads();
      }
    }

    double e[5] = {0, 0, 0, 0, 0};
    int clashPairs = 0;
    // Warp w takes pairs base + lane. The pairs are sorted by i: the i contributions of a warp are
    // summed in T per i (segmented shuffle sum; one plain sum when all pairs share i) and added
    // once by the first lane of each segment. The pair of the next iteration is loaded ahead.
    GpuPair next{-1, 0, 0};
    if (warp * 32 + lane < pd.count) next = pd.pairs[warp * 32 + lane];
    for (int base = warp * 32; base < pd.count; base += kThreads) {
      const int p = base + lane;
      const GpuPair pr = next;
      if (p + kThreads < pd.count) next = pd.pairs[p + kThreads];
      T d[3], f;
      const bool ok = p < pd.count && pairTermGpu(m, pr, pd, x, e, clashPairs, d, f);
      if (!g) continue;
      T v[3] = {0, 0, 0};
      const int key = pr.i;
      if (ok && f != 0)
        for (int c = 0; c < 3; ++c) {
          v[c] = f * d[c];
          atomicAdd(g + c * cc + pr.j, -static_cast<double>(v[c]));
        }
      const int prev = __shfl_up_sync(0xffffffffu, key, 1);
      const unsigned heads = __ballot_sync(0xffffffffu, lane == 0 || prev != key);
      if (heads == 1u) {
        for (int o = 16; o > 0; o /= 2)
          for (int c = 0; c < 3; ++c) v[c] += __shfl_down_sync(0xffffffffu, v[c], o);
      } else {
        for (int o = 1; o < 32; o *= 2) {
          const int k = __shfl_down_sync(0xffffffffu, key, o);
          T nv[3];
          for (int c = 0; c < 3; ++c) nv[c] = __shfl_down_sync(0xffffffffu, v[c], o);
          if (lane + o < 32 && k == key)
            for (int c = 0; c < 3; ++c) v[c] += nv[c];
        }
      }
      if (key >= 0 && (heads >> lane & 1u) && (v[0] != 0 || v[1] != 0 || v[2] != 0))
        for (int c = 0; c < 3; ++c) atomicAdd(g + c * cc + key, static_cast<double>(v[c]));
    }
    double closure = 0, maxDev = 0;
    for (int c = tid; c < m.closureCount; c += kThreads) {
      T d[3], f;
      flat::closureTerm(m, c, x, closure, maxDev, d, f);
      if (!g) continue;
      const FlatClosure& cl = m.closures[c];
      for (int a = 0; a < 3; ++a) {
        const double v = f * d[a];
        atomicAdd(g + a * cc + cl.i, v);
        atomicAdd(g + a * cc + cl.j, -v);
      }
    }

    constexpr int kVals = Dock ? 17 : 8;
    double v[kVals] = {e[0], e[1], e[2], e[3], e[4], closure, maxDev, static_cast<double>(clashPairs)};
    if constexpr (Dock) {
      // 8: inter, 9: interClash, 10: outside, 11–13: translation, 14–16: rotation gradient.
      const DockPose ps = poses[s];
      // The field of this sample (spec 15.1.16): the kernel parameter, or of DockView::fields.
      if constexpr (Dock == 2)
        probeTerms<kThreads>(dv, dv.fields[ps.field], x, ps, g, cc, v);
      else
        probeTerms<kThreads>(dv, dv.field, x, ps, g, cc, v);
    }
    __shared__ double red[kWarps][kVals];
    for (int k = 0; k < kVals; ++k) {
      v[k] = k == 6 ? warpMax(v[k]) : warpSum(v[k]);
      if (lane == 0) red[warp][k] = v[k];
    }
    __syncthreads();
    if (tid == 0) {
      double r[kVals];
      for (int k = 0; k < kVals; ++k) {
        r[k] = red[0][k];
        for (int w = 1; w < kWarps; ++w) r[k] = k == 6 ? fmax(r[k], red[w][k]) : r[k] + red[w][k];
      }
      FlatTerms out;
      out.vdw = r[0];
      out.hbond = r[1];
      out.elec = r[2];
      out.desolv = r[3];
      out.clash = r[4];
      out.closure = r[5];
      out.closureMaxDeviation = r[6];
      out.clashPairs = static_cast<int>(r[7]);
      terms[s] = out;
      if constexpr (Dock) {
        dterms[s] = {r[8], r[9], static_cast<int>(r[10])};
        if (grad6)
          for (int i = 0; i < 6; ++i) grad6[static_cast<std::size_t>(s) * 6 + i] = r[11 + i];
      }
    }
    if (shared && keepX) {
      for (std::size_t i = tid; i < nc; i += kThreads) xs[s * nc + i] = x[i];
    }
    for (int t = warp; g && t < m.torsionCount; t += kWarps) {
      T a[3], u[3];
      double sum = 0;
      if (flat::torsionUnit(m, t, x, a, u)) {
        const FlatTorsion& tt = m.torsions[t];
        for (int q = tt.movingBegin + lane; q < tt.movingEnd; q += 32) {
          const int idx = m.moving[q];
          const T w[3] = {x[3 * idx] - a[0], x[3 * idx + 1] - a[1], x[3 * idx + 2] - a[2]};
          sum += g[idx] * (u[1] * w[2] - u[2] * w[1]) + g[cc + idx] * (u[2] * w[0] - u[0] * w[2]) +
                 g[2 * cc + idx] * (u[0] * w[1] - u[1] * w[0]);
        }
      }
      sum = warpSum(sum);
      if (lane == 0) grad[static_cast<std::size_t>(s) * m.torsionCount + t] = sum;
    }
    __syncthreads();
  }
}

}  // namespace

int gpuDeviceCount() {
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess) {
    cudaGetLastError();
    return 0;
  }
  return n;
}

GpuDeviceInfo gpuDeviceInfo(int device) {
  cudaDeviceProp p;
  check(cudaGetDeviceProperties(&p, device), "cudaGetDeviceProperties");
  GpuDeviceInfo info;
  info.name = p.name;
  info.major = p.major;
  info.minor = p.minor;
  int old = 0;
  check(cudaGetDevice(&old), "cudaGetDevice");
  check(cudaSetDevice(device), "cudaSetDevice");
  check(cudaMemGetInfo(&info.freeMem, &info.totalMem), "cudaMemGetInfo");
  check(cudaSetDevice(old), "cudaSetDevice");
  return info;
}

std::size_t gpuFreeMem(int device) {
  int old = 0;
  check(cudaGetDevice(&old), "cudaGetDevice");
  check(cudaSetDevice(device), "cudaSetDevice");
  std::size_t freeMem = 0, total = 0;
  check(cudaMemGetInfo(&freeMem, &total), "cudaMemGetInfo");
  check(cudaSetDevice(old), "cudaSetDevice");
  return freeMem;
}

void gpuDownload(const float* device, float* host, std::size_t n) {
  check(cudaMemcpy(host, device, n * sizeof(float), cudaMemcpyDeviceToHost), "gpuDownload");
}

// Reference maps of a field, shared by the slots with the same grids (FieldModel::base).
struct FieldMaps {
  std::shared_ptr<const void> base, baseDist;
  DevBuf<float> maps;
};

// One resident field.
struct FieldSlot {
  std::uint64_t id = 0, used = 0;
  std::shared_ptr<FieldMaps> maps;
  DevBuf<float> data;
  std::shared_ptr<const DeviceBlocks> deviceData;  // instead of data (spec 15.1.19)
  DevBuf<int> block, movedHeavy, cellStart, cellAtoms, heavyStart, heavyAtoms, pairM, pairN;
  DevBuf<FieldMoved> moved;
  DevBuf<double> solT, volT, pairAC, pairBD, pairReqm;
  FieldView view{};
  std::size_t bytes() const {
    return data.n * sizeof(float) + (deviceData ? deviceData->floats * sizeof(float) : 0) + block.n * sizeof(int) +
           moved.n * sizeof(FieldMoved);
  }
};

struct GpuField::Impl {
  int device = 0;
  std::vector<std::unique_ptr<FieldSlot>> slots;  // capacity: slots.size()
  int current = -1;
  std::uint64_t clock = 0;
  DevBuf<FieldView> views;  // view of each slot

  void refreshViews() {
    std::vector<FieldView> v(slots.size(), FieldView{});
    for (std::size_t i = 0; i < v.size(); ++i)
      if (slots[i]) v[i] = slots[i]->view;
    // Grown geometrically: cudaFree synchronizes the device (spec 15.1.20).
    if (v.size() > views.n) views.reserve(std::max(v.size(), 2 * views.n));
    views.upload(v);
  }
};

GpuField::GpuField(int device) : impl_(std::make_unique<Impl>()) {
  impl_->device = device;
  impl_->slots.resize(1);
}
GpuField::~GpuField() = default;
std::uint64_t GpuField::id() const {
  const Impl& d = *impl_;
  return d.current >= 0 ? d.slots[d.current]->id : 0;
}
const FieldView& GpuField::deviceView() const { return impl_->slots[impl_->current]->view; }

std::size_t GpuField::deviceBytes() const {
  const Impl& d = *impl_;
  std::size_t b = 0;
  std::vector<const FieldMaps*> seen;
  for (const auto& s : d.slots) {
    if (!s || !s->id) continue;
    b += s->bytes();
    if (s->maps && std::find(seen.begin(), seen.end(), s->maps.get()) == seen.end()) {
      seen.push_back(s->maps.get());
      b += s->maps->maps.n * sizeof(float);
    }
  }
  return b;
}

int GpuField::slots() const { return static_cast<int>(impl_->slots.size()); }

void GpuField::setSlots(int n) {
  Impl& d = *impl_;
  n = std::max(n, 1);
  if (n == static_cast<int>(d.slots.size())) return;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  // Keeps the most recently used fields.
  std::stable_sort(d.slots.begin(), d.slots.end(), [](const auto& a, const auto& b) {
    return (a ? a->used : 0) > (b ? b->used : 0);
  });
  d.slots.resize(static_cast<std::size_t>(n));
  d.current = -1;
  std::uint64_t best = 0;
  for (int i = 0; i < n; ++i)
    if (d.slots[i] && d.slots[i]->id && d.slots[i]->used > best) {
      best = d.slots[i]->used;
      d.current = i;
    }
  d.refreshViews();
}

int GpuField::slotOf(std::uint64_t id) const {
  const Impl& d = *impl_;
  if (id == 0) return -1;
  for (std::size_t i = 0; i < d.slots.size(); ++i)
    if (d.slots[i] && d.slots[i]->id == id) return static_cast<int>(i);
  return -1;
}

const FieldView* GpuField::deviceViews() const { return impl_->views.p; }

bool GpuField::touch(std::uint64_t id) {
  Impl& d = *impl_;
  const int have = slotOf(id);
  if (have < 0) return false;
  d.current = have;
  d.slots[have]->used = ++d.clock;
  return true;
}

void GpuField::set(const FieldModel& fm) {
  Impl& d = *impl_;
  if (touch(fm.id)) return;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  // An empty slot, otherwise the least recently used one.
  int k = 0;
  for (int i = 0; i < static_cast<int>(d.slots.size()); ++i) {
    if (!d.slots[i] || !d.slots[i]->id) {
      k = i;
      break;
    }
    if (d.slots[i]->used < d.slots[k]->used) k = i;
  }
  const int nm = fm.view.typeMapCount + 3;
  std::shared_ptr<FieldMaps> maps;
  for (const auto& s : d.slots)
    if (s && s->maps && s->maps->base == fm.base && s->maps->baseDist == fm.baseDist) maps = s->maps;
  d.slots[k].reset();
  d.current = -1;
  if (!maps) {
    maps = std::make_shared<FieldMaps>();
    maps->maps.reserve(static_cast<std::size_t>(nm) * fm.points);
    for (int m = 0; m < nm; ++m)
      check(cudaMemcpy(maps->maps.p + static_cast<std::size_t>(m) * fm.points, fm.view.maps[m], fm.points * sizeof(float),
                       cudaMemcpyHostToDevice),
            "upload map");
    maps->base = fm.base;
    maps->baseDist = fm.baseDist;
  }
  auto sp = std::make_unique<FieldSlot>();
  FieldSlot& s = *sp;
  s.maps = std::move(maps);
  // Delta blocks already in device memory (built there) are used in place.
  if (fm.deviceData && fm.deviceData->device == d.device) s.deviceData = fm.deviceData;
  else if (fm.data) s.data.upload(*fm.data);
  else if (fm.deviceData) throw Error("GpuField: delta blocks on another device");
  else s.data.upload(std::vector<float>{});
  s.block.upload(fm.block);
  s.moved.upload(fm.moved);
  s.movedHeavy.upload(fm.movedHeavy);
  s.cellStart.upload(fm.cellStart);
  s.cellAtoms.upload(fm.cellAtoms);
  s.heavyStart.upload(fm.heavyStart);
  s.heavyAtoms.upload(fm.heavyAtoms);
  s.solT.upload(fm.solT);
  s.volT.upload(fm.volT);
  s.pairAC.upload(fm.pairAC);
  s.pairBD.upload(fm.pairBD);
  s.pairReqm.upload(fm.pairReqm);
  s.pairM.upload(fm.pairM);
  s.pairN.upload(fm.pairN);
  FieldView v = fm.view;
  for (int m = 0; m < nm; ++m) v.maps[m] = s.maps->maps.p + static_cast<std::size_t>(m) * fm.points;
  v.data = s.deviceData ? s.deviceData->ptr : s.data.p;
  v.block = s.block.p;
  v.moved = s.moved.p;
  v.movedHeavy = s.movedHeavy.p;
  v.cells.start = s.cellStart.p;
  v.cells.atoms = s.cellAtoms.p;
  v.heavyCells.start = s.heavyStart.p;
  v.heavyCells.atoms = s.heavyAtoms.p;
  v.solT = s.solT.p;
  v.volT = s.volT.p;
  v.pairAC = s.pairAC.p;
  v.pairBD = s.pairBD.p;
  v.pairReqm = s.pairReqm.p;
  v.pairM = s.pairM.p;
  v.pairN = s.pairN.p;
  s.view = v;
  s.id = fm.id;
  s.used = ++d.clock;
  d.slots[k] = std::move(sp);
  d.current = k;
  d.refreshViews();
}

namespace {

// Smallest block (32 to kMaxThreads threads) with about 8 pairs per thread.
int blockThreads(const EnergyModelView& m, int probes) {
  int nt = 32;
  while (nt < kMaxThreads && (nt * 8 < m.pairCount || nt < probes)) nt *= 2;
  return nt;
}

// Largest dynamic shared memory per block (opt-in).
int sharedLimit() {
  static const int limit = [] {
    int dev = 0, v = 0;
    cudaGetDevice(&dev);
    if (cudaDeviceGetAttribute(&v, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev) != cudaSuccess) {
      cudaGetLastError();
      v = 48 * 1024;
    }
    return v;
  }();
  return limit;
}
// Shared memory per multiprocessor of the current device.
int sharedPerSm() {
  static const int v = [] {
    int dev = 0, s = 0;
    cudaGetDevice(&dev);
    if (cudaDeviceGetAttribute(&s, cudaDevAttrMaxSharedMemoryPerMultiprocessor, dev) != cudaSuccess) {
      cudaGetLastError();
      s = 48 * 1024;
    }
    return s;
  }();
  return v;
}
// Multiprocessor count and L2 cache size of the current device.
std::pair<int, std::size_t> smsAndL2() {
  static const std::pair<int, std::size_t> v = [] {
    int dev = 0, sms = 0, l2 = 0;
    cudaGetDevice(&dev);
    if (cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, dev) != cudaSuccess ||
        cudaDeviceGetAttribute(&l2, cudaDevAttrL2CacheSize, dev) != cudaSuccess) {
      cudaGetLastError();
      sms = 1;
      l2 = 0;
    }
    return std::pair<int, std::size_t>(std::max(sms, 1), static_cast<std::size_t>(std::max(l2, 0)));
  }();
  return v;
}
// Dynamic shared memory for the coordinates and torsion cos and sin of one sample; 0: they do not fit.
template <class T>
std::size_t sharedBytes(const EnergyModelView& m) {
  const std::size_t b = (3 * static_cast<std::size_t>(m.coordCount) + 2 * static_cast<std::size_t>(m.torsionCount)) * sizeof(T);
  return b <= static_cast<std::size_t>(sharedLimit()) - 4096 ? b : 0;  // the kernel's own shared memory
}

template <class T, int Dock, int NT>
void launchNT(cudaStream_t st, int n, std::size_t shm, bool keepX, const EnergyModelView& m, const PairData<T>& pd, const TorsionData<T>& td, const T* angles, T* xs,
              double* gs, FlatTerms* terms, double* grad, const DockView& dv, const DockPose* poses, DockTerms* dterms,
              double* grad6) {
  const auto kernelS = evaluateKernel<T, Dock, NT, false>;
  const auto kernelT = evaluateKernel<T, Dock, NT, true>;
  const auto [sms, l2] = smsAndL2();
  // Resident blocks per multiprocessor: their gradients (3 × coordCount doubles each) within the
  // L2 cache, at least 24 warps (spec 15.1.15); the blocks loop over the samples.
  constexpr int kMinBlocks = (24 * 32 + NT - 1) / NT;
  const auto perSm = [&](auto kernel, std::size_t bytes) {
    int occ = 0;
    if (bytes > 48 * 1024)
      check(cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(bytes)),
            "cudaFuncSetAttribute");
    check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&occ, kernel, NT, bytes), "occupancy");
    int r = std::max(occ, 1);
    if (gs) {
      const std::size_t g = 3 * static_cast<std::size_t>(m.coordCount) * sizeof(double) * sms;
      const int fit = static_cast<int>(std::min<std::size_t>(l2 / std::max<std::size_t>(g, 1), 1 << 20));
      r = std::clamp(std::max(fit, kMinBlocks), 1, r);
    }
    return r;
  };
  int per = perSm(kernelS, shm);
  // Torsion transforms (TorsionData): their work memory in shared memory when it does not lower
  // the resident blocks; with one wave of blocks otherwise in global memory. With more waves,
  // only in shared memory with at most half of the shared memory of a multiprocessor in use
  // (the rest stays L1 cache); otherwise the torsions are applied one by one.
  TorsionData<T> tdl = td;
  if (tdl.ok) {
    const std::size_t withM = shm + (24 * sizeof(T) + 2 * sizeof(int)) * static_cast<std::size_t>(m.torsionCount);
    const bool oneWave = n <= sms * per;
    const int perM = shm && withM <= static_cast<std::size_t>(sharedLimit()) - 4096 ? perSm(kernelT, withM) : 0;
    if (perM > 0 && (oneWave ? n <= sms * perM
                             : perM >= per && static_cast<std::size_t>(perM) * withM <= static_cast<std::size_t>(sharedPerSm()) / 2)) {
      shm = withM;
      per = perM;
      tdl.shared = 1;
    } else if (!oneWave) {
      tdl.ok = 0;
    } else {
      per = perSm(kernelT, shm);
    }
  }
  const int blocks = std::min(n, sms * per);
  (tdl.ok ? kernelT : kernelS)<<<blocks, NT, shm, st>>>(m, pd, tdl, angles, xs, gs, terms, grad, dv, poses, dterms, grad6,
                                                        n, shm ? 1 : 0, keepX ? 1 : 0);
}

// Samples evaluated at the same time by launchNT (resident blocks of all multiprocessors), without
// torsion transforms.
template <class T, int Dock, int NT>
int residentNT(const EnergyModelView& m, std::size_t shm) {
  const auto kernel = evaluateKernel<T, Dock, NT, false>;
  const auto [sms, l2] = smsAndL2();
  constexpr int kMinBlocks = (24 * 32 + NT - 1) / NT;
  int occ = 0;
  if (shm > 48 * 1024)
    check(cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(shm)),
          "cudaFuncSetAttribute");
  check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&occ, kernel, NT, shm), "occupancy");
  const std::size_t g = 3 * static_cast<std::size_t>(m.coordCount) * sizeof(double) * sms;
  const int fit = static_cast<int>(std::min<std::size_t>(l2 / std::max<std::size_t>(g, 1), 1 << 20));
  return sms * std::clamp(std::max(fit, kMinBlocks), 1, std::max(occ, 1));
}

template <class T, int Dock>
void launch(cudaStream_t st, int threads, int n, std::size_t shm, bool keepX, const EnergyModelView& m, const PairData<T>& pd, const TorsionData<T>& td, const T* angles,
            T* xs, double* gs, FlatTerms* terms, double* grad, const DockView& dv, const DockPose* poses,
            DockTerms* dterms, double* grad6) {
  switch (threads) {
    case 32: launchNT<T, Dock, 32>(st, n, shm, keepX, m, pd, td, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 64: launchNT<T, Dock, 64>(st, n, shm, keepX, m, pd, td, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 128: launchNT<T, Dock, 128>(st, n, shm, keepX, m, pd, td, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    default: launchNT<T, Dock, 256>(st, n, shm, keepX, m, pd, td, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
  }
}

// Page-locked host memory.
struct PinnedBuf {
  char* p = nullptr;
  std::size_t n = 0;
  PinnedBuf() = default;
  PinnedBuf(const PinnedBuf&) = delete;
  PinnedBuf& operator=(const PinnedBuf&) = delete;
  ~PinnedBuf() {
    if (p) cudaFreeHost(p);
  }
  void reserve(std::size_t bytes) {
    if (bytes <= n) return;
    release();
    check(cudaMallocHost(&p, bytes), "cudaMallocHost");
    n = bytes;
  }
  void release() {
    if (p) cudaFreeHost(p);
    p = nullptr;
    n = 0;
  }
};

constexpr std::size_t align8(std::size_t b) { return (b + 7) & ~std::size_t(7); }

// Work memory and stream of one lane (spec 15.1.8): inputs and outputs are each one contiguous
// block (one copy per call).
struct Lane {
  DevBuf<char> in, out, xs;
  DevBuf<double> gs;
  DevBuf<char> tbuf;
  PinnedBuf hIn, hOut;
  cudaStream_t stream = nullptr;
  // Pending dock evaluation.
  std::size_t n = 0, outBytes = 0, oDterms = 0, oGrad6 = 0, oGrad = 0;
  bool grad = false, pending = false;

  Lane() = default;
  Lane(const Lane&) = delete;
  Lane& operator=(const Lane&) = delete;
  ~Lane() {
    if (stream) cudaStreamDestroy(stream);
  }
  void release() {
    in.release();
    out.release();
    xs.release();
    gs.release();
    tbuf.release();
    hIn.release();
    hOut.release();
  }
};

}  // namespace

struct GpuEnergy::Impl {
  int device = 0;
  EnergyModelView view{};
  DevBuf<double> ref;
  DevBuf<FlatTorsion> torsions;
  DevBuf<int> moving;
  DevBuf<FlatClosure> closures;
  DevBuf<GpuPair> gpuPairs;
  DevBuf<GpuAtom<float>> atomsF;
  DevBuf<GpuAtom<double>> atomsD;
  DevBuf<GpuPairParams<float>> paramsF;
  DevBuf<GpuPairParams<double>> paramsD;
  Lane lanes[GpuEnergy::kLanes];
  DevBuf<float> axisF;
  DevBuf<double> axisD;
  DevBuf<int> tParent, tDeep;
  int tRounds = 0;
  bool tOk = false;

  template <class T>
  TorsionData<T> torsionData(Lane& L, std::size_t n) {
    if (!tOk) return {nullptr, nullptr, nullptr, nullptr, 0, 0, 0};
    L.tbuf.reserve(std::max<std::size_t>(1, n * (24 * sizeof(T) + 2 * sizeof(int)) * static_cast<std::size_t>(view.torsionCount)));
    if constexpr (std::is_same_v<T, float>)
      return {axisF.p, tParent.p, tDeep.p, L.tbuf.p, tRounds, 1, 0};
    else
      return {axisD.p, tParent.p, tDeep.p, L.tbuf.p, tRounds, 1, 0};
  }

  // Torsion transforms (TorsionData): tOk false when the conditions do not hold.
  void planTorsions(const EnergyModel& m) {
    const int nt = m.torsionCount, nc = m.coordCount;
    tOk = false;
    if (nt == 0) return;
    std::vector<std::vector<int>> by(nc);  // torsions moving each coordinate, ascending
    for (int u = 0; u < nt; ++u)
      for (int q = m.torsions[u].movingBegin; q < m.torsions[u].movingEnd; ++q) by[m.moving[q]].push_back(u);
    std::vector<int> parent(nt, -1), deep(nc, -1), depth(nt, 1);
    for (int u = 0; u < nt; ++u) {
      const FlatTorsion& tt = m.torsions[u];
      for (int a : {tt.fixedAtom, tt.movingAtom})
        for (int w : by[a])
          if (w > u) return;
      if (tt.movingEnd == tt.movingBegin) continue;
      const std::vector<int>& b = by[m.moving[tt.movingBegin]];
      const auto it = std::lower_bound(b.begin(), b.end(), u);
      if (it != b.begin()) parent[u] = *(it - 1);
      if (parent[u] >= 0) depth[u] = depth[parent[u]] + 1;
    }
    int maxDepth = 1;
    for (int a = 0; a < nc; ++a) {
      if (by[a].empty()) continue;
      deep[a] = by[a].back();
      int u = deep[a];
      for (auto it = by[a].rbegin(); it != by[a].rend(); ++it, u = parent[u])
        if (u != *it) return;
      if (u != -1) return;
      maxDepth = std::max(maxDepth, depth[deep[a]]);
    }
    std::vector<double> axis(6 * static_cast<std::size_t>(nt), 0.0);
    for (int u = 0; u < nt; ++u) {
      const FlatTorsion& tt = m.torsions[u];
      double o[3], k[3];
      for (int c = 0; c < 3; ++c) {
        o[c] = m.ref[3 * tt.fixedAtom + c];
        k[c] = m.ref[3 * tt.movingAtom + c] - o[c];
      }
      const double len = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
      for (int c = 0; c < 3; ++c) {
        axis[6 * u + c] = o[c];
        axis[6 * u + 3 + c] = len == 0 ? 0.0 : k[c] / len;
      }
    }
    axisD.upload(axis);
    axisF.upload(std::vector<float>(axis.begin(), axis.end()));
    tParent.upload(parent);
    tDeep.upload(deep);
    tRounds = 0;
    while ((1 << tRounds) < maxDepth) ++tRounds;
    tOk = true;
  }

  template <class T>
  PairData<T> pairData() const {
    if constexpr (std::is_same_v<T, float>)
      return {gpuPairs.p, atomsF.p, paramsF.p, view.pairCount};
    else
      return {gpuPairs.p, atomsD.p, paramsD.p, view.pairCount};
  }

  // desolvCoeff·(Si·Vj + Sj·Vi) as sv[i]·vol[j] + sv[j]·vol[i] (spec 15.1.15).
  template <class T>
  void uploadAtoms(const EnergyModel& m, DevBuf<GpuAtom<T>>& da, DevBuf<GpuPairParams<T>>& dpp) {
    const T elecCoeff = static_cast<T>(m.elecCoeff), desolvCoeff = static_cast<T>(m.desolvCoeff);
    std::vector<GpuAtom<T>> a(m.charge.size());
    for (std::size_t i = 0; i < a.size(); ++i)
      a[i] = {elecCoeff * static_cast<T>(m.charge[i]), static_cast<T>(m.charge[i]), desolvCoeff * static_cast<T>(m.sv[i]),
              static_cast<T>(m.vol[i])};
    da.upload(a);
    std::vector<GpuPairParams<T>> pp;
    for (const FlatPairParams& q : m.pairParams)
      pp.push_back({static_cast<T>(q.AC), static_cast<T>(q.BD), static_cast<T>(q.reqm), q.m, q.n});
    dpp.upload(pp);
  }
  DevBuf<DockProbe> probes;
  DockView dock{};

  template <class T>
  void run(const T* hAngles, int n, FlatTerms* hTerms, double* hGrad, T* hCoords) {
    if (n <= 0) return;
    check(cudaSetDevice(device), "cudaSetDevice");
    Lane& L = lanes[0];
    if (L.pending) throw Error("evaluate: lane 0 busy");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t inBytes = align8(nn * nt * sizeof(T));
    const std::size_t oTerms = 0, oGrad = align8(nn * sizeof(FlatTerms));
    const std::size_t outBytes = oGrad + (hGrad ? nn * nt * sizeof(double) : 0);
    const std::size_t shm = sharedBytes<T>(view);
    L.in.reserve(std::max<std::size_t>(8, inBytes));
    L.out.reserve(outBytes);
    if (!shm || hCoords) L.xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(T));
    if (hGrad) L.gs.reserve(std::max<std::size_t>(1, nn * nc));
    L.hIn.reserve(std::max<std::size_t>(8, inBytes));
    L.hOut.reserve(outBytes);
    if (nt) {
      std::copy_n(hAngles, nn * nt, reinterpret_cast<T*>(L.hIn.p));
      check(cudaMemcpyAsync(L.in.p, L.hIn.p, inBytes, cudaMemcpyHostToDevice, L.stream), "upload");
    }
    launch<T, 0>(L.stream, blockThreads(view, 0), n, shm, hCoords != nullptr, view, pairData<T>(), torsionData<T>(L, nn),
                     reinterpret_cast<const T*>(L.in.p), reinterpret_cast<T*>(L.xs.p), hGrad ? L.gs.p : nullptr,
                     reinterpret_cast<FlatTerms*>(L.out.p + oTerms),
                     hGrad ? reinterpret_cast<double*>(L.out.p + oGrad) : nullptr, DockView{}, nullptr, nullptr, nullptr);
    check(cudaGetLastError(), "evaluateKernel launch");
    check(cudaMemcpyAsync(L.hOut.p, L.out.p, outBytes, cudaMemcpyDeviceToHost, L.stream), "download");
    check(cudaStreamSynchronize(L.stream), "evaluateKernel");
    std::copy_n(reinterpret_cast<const FlatTerms*>(L.hOut.p + oTerms), nn, hTerms);
    if (hGrad) std::copy_n(reinterpret_cast<const double*>(L.hOut.p + oGrad), nn * nt, hGrad);
    if (hCoords) check(cudaMemcpy(hCoords, L.xs.p, nn * nc * sizeof(T), cudaMemcpyDeviceToHost), "download coordinates");
  }

  template <class T>
  void submitEval(int lane, const double* hAngles, int n, bool wantGrad) {
    Lane& L = lanes[lane];
    if (L.pending) throw Error("submit: lane busy");
    L.n = 0;
    if (n <= 0) return;
    check(cudaSetDevice(device), "cudaSetDevice");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t inBytes = align8(nn * nt * sizeof(T));
    L.oGrad = align8(nn * sizeof(FlatTerms));
    L.outBytes = L.oGrad + (wantGrad ? nn * nt * sizeof(double) : 0);
    const std::size_t shm = sharedBytes<T>(view);
    L.in.reserve(std::max<std::size_t>(8, inBytes));
    L.out.reserve(L.outBytes);
    if (!shm) L.xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(T));
    if (wantGrad) L.gs.reserve(std::max<std::size_t>(1, nn * nc));
    L.hIn.reserve(std::max<std::size_t>(8, inBytes));
    L.hOut.reserve(L.outBytes);
    if (nt) {
      T* ha = reinterpret_cast<T*>(L.hIn.p);
#pragma omp parallel for schedule(static)
      for (long long s = 0; s < static_cast<long long>(nn * nt); ++s) ha[s] = static_cast<T>(hAngles[s]);
      check(cudaMemcpyAsync(L.in.p, L.hIn.p, inBytes, cudaMemcpyHostToDevice, L.stream), "upload");
    }
    launch<T, 0>(L.stream, blockThreads(view, 0), n, shm, false, view, pairData<T>(), torsionData<T>(L, nn),
                          reinterpret_cast<const T*>(L.in.p), reinterpret_cast<T*>(L.xs.p),
                          wantGrad ? L.gs.p : nullptr, reinterpret_cast<FlatTerms*>(L.out.p),
                          wantGrad ? reinterpret_cast<double*>(L.out.p + L.oGrad) : nullptr, DockView{}, nullptr,
                          nullptr, nullptr);
    check(cudaGetLastError(), "evaluateKernel launch");
    check(cudaMemcpyAsync(L.hOut.p, L.out.p, L.outBytes, cudaMemcpyDeviceToHost, L.stream), "download");
    L.n = nn;
    L.grad = wantGrad;
    L.pending = true;
  }

  void completeEval(int lane, FlatTerms* hTerms, double* hGrad) {
    Lane& L = lanes[lane];
    if (!L.pending) return;
    L.pending = false;
    check(cudaStreamSynchronize(L.stream), "evaluateKernel");
    const std::size_t nn = L.n, nt = static_cast<std::size_t>(view.torsionCount);
    std::copy_n(reinterpret_cast<const FlatTerms*>(L.hOut.p), nn, hTerms);
    if (L.grad && hGrad) {
      const double* g = reinterpret_cast<const double*>(L.hOut.p + L.oGrad);
#pragma omp parallel for schedule(static)
      for (long long s = 0; s < static_cast<long long>(nn); ++s) std::copy_n(g + s * nt, nt, hGrad + s * nt);
    }
  }

  void submitDock(int lane, const GpuField::Impl& f, const float* hAngles, const DockPose* hPoses, int n, bool wantGrad) {
    Lane& L = lanes[lane];
    if (L.pending) throw Error("submitDock: lane busy");
    L.n = 0;
    if (n <= 0) return;
    if (f.current < 0) throw Error("evaluateDock: no field");
    check(cudaSetDevice(device), "cudaSetDevice");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t iPoses = align8(nn * nt * sizeof(float)), inBytes = iPoses + nn * sizeof(DockPose);
    L.oDterms = align8(nn * sizeof(FlatTerms));
    L.oGrad6 = L.oDterms + align8(nn * sizeof(DockTerms));
    L.oGrad = L.oGrad6 + nn * 6 * sizeof(double);
    L.outBytes = wantGrad ? L.oGrad + nn * nt * sizeof(double) : L.oGrad6;
    const std::size_t shm = sharedBytes<float>(view);
    L.in.reserve(inBytes);
    L.out.reserve(L.outBytes);
    if (!shm) L.xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(float));
    if (wantGrad) L.gs.reserve(std::max<std::size_t>(1, nn * nc));
    L.hIn.reserve(inBytes);
    L.hOut.reserve(L.outBytes);
    float* ha = reinterpret_cast<float*>(L.hIn.p);
    DockPose* hp = reinterpret_cast<DockPose*>(L.hIn.p + iPoses);
    // Poses of other resident fields than the current one: DockView::fields (spec 15.1.16).
    int other = 0;
#pragma omp parallel for schedule(static) reduction(| : other)
    for (long long s = 0; s < static_cast<long long>(nn); ++s) {
      std::copy_n(hAngles + s * nt, nt, ha + s * nt);
      hp[s] = hPoses[s];
      other |= hp[s].field != f.current;
    }
    check(cudaMemcpyAsync(L.in.p, L.hIn.p, inBytes, cudaMemcpyHostToDevice, L.stream), "upload");
    DockView dv = dock;
    dv.field = f.slots[f.current]->view;
    dv.fields = other ? f.views.p : nullptr;
    const auto run = [&]<int D>() {
      launch<float, D>(L.stream, blockThreads(view, dv.probeCount), n, shm, false, view, pairData<float>(), torsionData<float>(L, nn),
                       reinterpret_cast<const float*>(L.in.p), reinterpret_cast<float*>(L.xs.p),
                       wantGrad ? L.gs.p : nullptr, reinterpret_cast<FlatTerms*>(L.out.p),
                       wantGrad ? reinterpret_cast<double*>(L.out.p + L.oGrad) : nullptr, dv,
                       reinterpret_cast<const DockPose*>(L.in.p + iPoses),
                       reinterpret_cast<DockTerms*>(L.out.p + L.oDterms),
                       wantGrad ? reinterpret_cast<double*>(L.out.p + L.oGrad6) : nullptr);
    };
    // A separate kernel for several fields: the single-field kernel reads its field from the
    // kernel parameters (spec 15.1.16).
    if (other) run.template operator()<2>();
    else run.template operator()<1>();
    check(cudaGetLastError(), "evaluateKernel (dock) launch");
    check(cudaMemcpyAsync(L.hOut.p, L.out.p, L.outBytes, cudaMemcpyDeviceToHost, L.stream), "download");
    L.n = nn;
    L.grad = wantGrad;
    L.pending = true;
  }

  void completeDock(int lane, FlatTerms* hTerms, DockTerms* hDterms, double* hGrad) {
    Lane& L = lanes[lane];
    if (!L.pending) return;
    L.pending = false;
    check(cudaStreamSynchronize(L.stream), "evaluateKernel (dock)");
    const std::size_t nn = L.n, nt = static_cast<std::size_t>(view.torsionCount), ng = 6 + nt;
    const FlatTerms* ft = reinterpret_cast<const FlatTerms*>(L.hOut.p);
    const DockTerms* dt = reinterpret_cast<const DockTerms*>(L.hOut.p + L.oDterms);
    const double* g6 = reinterpret_cast<const double*>(L.hOut.p + L.oGrad6);
    const double* gt = reinterpret_cast<const double*>(L.hOut.p + L.oGrad);
    const bool grad = L.grad && hGrad;
#pragma omp parallel for schedule(static)
    for (long long s = 0; s < static_cast<long long>(nn); ++s) {
      hTerms[s] = ft[s];
      hDterms[s] = dt[s];
      if (!grad) continue;
      std::copy_n(g6 + s * 6, 6, hGrad + s * ng);
      std::copy_n(gt + s * nt, nt, hGrad + s * ng + 6);
    }
  }
};

GpuEnergy::GpuEnergy(const EnergyModel& m, int device) : impl_(std::make_unique<Impl>()) {
  Impl& d = *impl_;
  d.device = device;
  check(cudaSetDevice(device), "cudaSetDevice");
  d.ref.upload(m.ref);
  d.torsions.upload(m.torsions);
  d.moving.upload(m.moving);
  d.closures.upload(m.closures);
  std::vector<GpuPair> gp(m.pairs.size());
  for (std::size_t p = 0; p < gp.size(); ++p) {
    const FlatPair& pr = m.pairs[p];
    gp[p] = {pr.i, pr.j, pr.type << 2 | (m.pairParams[pr.type].n == 10 ? 2 : 0) | (pr.heavy ? 1 : 0)};
  }
  d.gpuPairs.upload(gp);
  d.uploadAtoms(m, d.atomsF, d.paramsF);
  d.uploadAtoms(m, d.atomsD, d.paramsD);
  d.planTorsions(m);
  // The pair terms use PairData (charge, sv, vol, pairParams and pairs are not uploaded).
  d.view = m.view();
  d.view.ref = d.ref.p;
  d.view.charge = d.view.sv = d.view.vol = nullptr;
  d.view.torsions = d.torsions.p;
  d.view.moving = d.moving.p;
  d.view.pairParams = nullptr;
  d.view.pairs = nullptr;
  d.view.closures = d.closures.p;
  // Highest priority: not delayed by the delta grids of the next StaM conformation (spec 15.1.20).
  int least = 0, greatest = 0;
  check(cudaDeviceGetStreamPriorityRange(&least, &greatest), "cudaDeviceGetStreamPriorityRange");
  for (Lane& L : d.lanes)
    check(cudaStreamCreateWithPriority(&L.stream, cudaStreamNonBlocking, greatest), "cudaStreamCreate");
}

GpuEnergy::~GpuEnergy() = default;

std::size_t GpuEnergy::bytesPerSample(bool fp64, bool grad) const {
  const std::size_t nt = static_cast<std::size_t>(impl_->view.torsionCount), nc = 3 * static_cast<std::size_t>(impl_->view.coordCount);
  const std::size_t r = fp64 ? sizeof(double) : sizeof(float);
  return (nt + nc) * r + sizeof(FlatTerms) + (grad ? (nc + nt) * sizeof(double) : 0) + nt * (24 * r + 2 * sizeof(int));
}

void GpuEnergy::setProbes(const std::vector<DockProbe>& probes, double clashDistance, double clashBase, double clashK) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  d.probes.upload(probes);
  d.dock = {static_cast<int>(probes.size()), d.probes.p, clashDistance, clashBase, clashK, FieldView{}};
}

void GpuEnergy::evaluateDock(const GpuField& field, const float* angles, const DockPose* poses, int n, FlatTerms* terms,
                             DockTerms* dterms, double* grad) {
  impl_->submitDock(0, *field.impl_, angles, poses, n, grad != nullptr);
  impl_->completeDock(0, terms, dterms, grad);
}

void GpuEnergy::submitDock(int lane, const GpuField& field, const float* angles, const DockPose* poses, int n, bool grad) {
  impl_->submitDock(lane, *field.impl_, angles, poses, n, grad);
}

void GpuEnergy::completeDock(int lane, FlatTerms* terms, DockTerms* dterms, double* grad) {
  impl_->completeDock(lane, terms, dterms, grad);
}

std::size_t GpuEnergy::bytesPerSampleDock() const {
  return bytesPerSample(false, true) + sizeof(DockPose) + sizeof(DockTerms) + 6 * sizeof(double);
}

int GpuEnergy::residentSamplesDock(int probeCount) const {
  check(cudaSetDevice(impl_->device), "cudaSetDevice");
  const EnergyModelView& m = impl_->view;
  const std::size_t shm = sharedBytes<float>(m);
  switch (blockThreads(m, probeCount)) {
    case 32: return residentNT<float, 2, 32>(m, shm);
    case 64: return residentNT<float, 2, 64>(m, shm);
    case 128: return residentNT<float, 2, 128>(m, shm);
    default: return residentNT<float, 2, 256>(m, shm);
  }
}

void GpuEnergy::releaseScratch() {
  for (Lane& L : impl_->lanes) {
    if (L.pending) {
      cudaStreamSynchronize(L.stream);
      L.pending = false;
    }
    L.release();
  }
}

void GpuEnergy::submit(int lane, const double* angles, int n, bool grad, bool fp32) {
  if (fp32)
    impl_->submitEval<float>(lane, angles, n, grad);
  else
    impl_->submitEval<double>(lane, angles, n, grad);
}

void GpuEnergy::complete(int lane, FlatTerms* terms, double* grad) { impl_->completeEval(lane, terms, grad); }

void GpuEnergy::evaluate(const double* angles, int n, FlatTerms* terms, double* grad, double* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

void GpuEnergy::evaluate(const float* angles, int n, FlatTerms* terms, double* grad, float* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

}  // namespace od
