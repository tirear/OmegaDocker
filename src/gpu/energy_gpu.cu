// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "energy_gpu.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <cstdint>
#include <cmath>
#include <string>
#include <utility>
#include <type_traits>
#include <vector>

#include "error.hpp"
#include "minimize.hpp"
#include "pool.hpp"

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

// Pair data of the kernel: pairs (tf: type << 2 | hbond << 1 | clash), per-atom constants in T.
struct GpuPair {
  int i, j, tf;
};
inline int pairFlags(const EnergyModel& m, const FlatPair& pr) {
  return pr.type << 2 | (m.pairParams[pr.type].n == 10 ? 2 : 0) | (pr.clash ? 1 : 0);
}
// Aligned to their size: one vector load each.
template <class T>
struct alignas(4 * sizeof(T)) GpuAtom {
  T cq, q, sv, vol;  // elecCoeff·q, q, desolvCoeff·S, V
};
// 12-6, or 12-10 when the pair flags have the hbond bit (pairFlags; uploadAtoms).
template <class T>
struct alignas(4 * sizeof(T)) GpuPairParams {
  T AC, BD, reqm, unused;
};
// Pairs in tiles of 32 × 32 atoms (spec 15.1.24). Tile t: atoms tileI[32t + l] and tileJ[32t + l]
// (−1: none); at step s, lane l takes atom tileI[32t + l] with tileJ[32t + (l + s) % 32], flags
// slots[1024t + 32s + l] (GpuPair::tf, | kSwap: i is the tileJ atom; −1: no pair). The steps of
// tile t are tileS0[t] … tileS0[t] + tileOff[t + 1] − tileOff[t] − 1; the warps of a block take
// equal ranges of the tileOff[tiles] steps of all tiles.
// pairs: in EnergyModel order (count), used for double and, with tiled 0, for float.
constexpr int kSwap = 1 << 30;
template <class T>
struct PairData {
  const GpuAtom<T>* atoms;
  const GpuPairParams<T>* params;
  const int *tileI, *tileJ, *tileS0, *tileOff, *slots;
  int tiles;
  const GpuPair* pairs;
  int count, tiled;
};

// Reciprocal of pairTermGpu: float with the hardware approximation (about 2 ulp; positive normal
// arguments only), double exact.
template <class T>
__device__ __forceinline__ T pairRcp(T v) {
  if constexpr (std::is_same_v<T, float>) return __fdividef(1.0f, v);
  else return T(1) / v;
}

// One atom pair (the terms of flat::pairTerm with one reciprocal per function instead of
// divisions; float: pairRcp and r, 1 / r from the hardware reciprocal square root): adds to
// e[0..4] and clashPairs; d = xi − xj, f = dE/dr / r (0: no gradient). Returns false beyond the
// cutoff.
template <class T>
__device__ inline bool pairTermGpu(const EnergyModelView& m, const GpuPair& pr, const PairData<T>& pd, const T* x,
                                   double e[5], int& clashPairs, T d[3], T& f) {
  for (int a = 0; a < 3; ++a) d[a] = x[3 * pr.i + a] - x[3 * pr.j + a];
  const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (r2 >= static_cast<T>(kElecCutoff * kElecCutoff)) return false;
  T r, ir;  // r, 1 / r (0 at r = 0)
  if constexpr (std::is_same_v<T, float>) {
    ir = r2 > 0 ? rsqrtf(r2) : 0.0f;
    r = r2 * ir;
  } else {
    r = sqrt(r2);
    ir = r > 0 ? T(1) / r : T(0);
  }
  T dedr = 0;
  if (r < static_cast<T>(kVdwCutoff)) {
    const GpuPairParams<T> q = pd.params[pr.tf >> 2];
    // 12-6 or 12-10 (uploadAtoms).
    const T half = static_cast<T>(m.halfSmooth);
    T s = q.reqm, ds = 0;
    if (r <= q.reqm - half) {
      s = r + half;
      ds = 1;
    } else if (r >= q.reqm + half) {
      s = r - half;
      ds = 1;
    }
    const T inv = pairRcp(s), inv2 = inv * inv, inv6 = inv2 * inv2 * inv2;
    const bool hb = pr.tf & 2;
    const T rm = inv6 * inv6, rn = hb ? inv6 * inv2 * inv2 : inv6;
    const T ee = q.AC * rm - q.BD * rn;
    const T de = ds * ((hb ? T(10) : T(6)) * q.BD * rn - T(12) * q.AC * rm) * inv;
    e[hb ? 1 : 0] += ee;
    dedr += de;
  }
  // Electrostatics and desolvation (ad4Elec, ad4Desolv): evaluated at rMin below rMin.
  const T rMin = static_cast<T>(m.rMin);
  const bool clamp = r < rMin;
  const T rr = clamp ? rMin : r;
  const T irr = clamp ? pairRcp(rMin) : ir;
  const T sh = rr + T(1.26366), ish = pairRcp(sh), ish2 = ish * ish;
  const T diel = T(1.10859) * ish2 + T(0.010358);
  const GpuAtom<T> ai = pd.atoms[pr.i], aj = pd.atoms[pr.j];
  const T qq = ai.cq * aj.q, sv = ai.sv * aj.vol + aj.sv * ai.vol;
  const T es = qq * irr;
  e[2] += es * diel;
  const T rr2 = rr * rr;
  const T num = T(12.96) - T(0.1063) * rr2 * (T(1) - T(0.001947) * rr2);
  const T den = T(12.96) + rr2 * (T(0.4137) + rr2 * (T(0.00357) + T(0.000112) * rr2));
  const T iden = pairRcp(den);
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
  f = dedr * ir;
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
// Torsion gradient (spec 15.1.29): own[ownOff[u], ownOff[u + 1]) the coordinates whose deepest
// torsion is u; pos[u] the position of u in a depth-first order of the torsions, whose subtree
// takes the sub[u] positions from pos[u].
template <class T>
struct TorsionData {
  const T* axis;  // per torsion: o, unit k (k = 0: zero-length axis, A_t = identity)
  const int* parent;
  const int* deep;  // per coordinate; −1: not moved
  char* buf;
  int rounds, ok, shared;
  const int *ownOff, *own, *pos, *sub;
};

// Energy terms and gradients of one sample by the block of NT threads (evaluateKernel,
// minimizeDockKernel). ang: torsion angles; x: coordinate memory (shared: dynamic shared memory,
// followed by the torsion cos and sin); g: gradient work memory of the block (nullptr: no
// gradient); pose: the 9b pose (Dock). On return (all threads): res[0..7] vdw, hbond, elec, desolv,
// clash, closure, closureMaxDeviation, clashPairs; Dock: res[8..10] inter, interClash, outside,
// res[11..16] translation and rotation gradient; with g, gradT the torsion gradient.
template <class T, int Dock, int NT, bool Tor>
__device__ __forceinline__ void evalOne(const EnergyModelView& m, const PairData<T>& pd, const TorsionData<T>& td,
                                        const T* ang, T* x, double* g, int shared, const DockView& dv,
                                        const DockPose* pose, double* res, double* gradT) {
  constexpr int kThreads = NT, kWarps = NT / 32;
  const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount);
  const int cc = m.coordCount;
  (void)pose;
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
    // Pair tiles (PairData): lane l keeps the gradient of its tileI atom and of the tileJ atom it
    // takes; the latter passes to lane l − 1 after each step. Both are summed in T and added once
    // per tile range of the warp.
    bool tiled = false;
    if constexpr (std::is_same_v<T, float>) tiled = pd.tiled;
    if (tiled) {
      const int steps = pd.tileOff[pd.tiles];
      int sb = static_cast<int>(static_cast<long long>(steps) * warp / kWarps);
      const int se = static_cast<int>(static_cast<long long>(steps) * (warp + 1) / kWarps);
      int tl = 0, th = pd.tiles;
      while (th - tl > 1) {
        const int mid = (tl + th) / 2;
        if (pd.tileOff[mid] <= sb) tl = mid; else th = mid;
      }
      for (int tt = tl; sb < se; ++tt) {
        const int te = min(se, pd.tileOff[tt + 1]);
        const int s0 = pd.tileS0[tt] + (sb - pd.tileOff[tt]), s1 = s0 + (te - sb);
        sb = te;
        const int ia = pd.tileI[32 * tt + lane], jr = pd.tileJ[32 * tt + lane];
        const int* sl = pd.slots + 1024 * static_cast<std::size_t>(tt) + lane;
        T gi[3] = {0, 0, 0}, gj[3] = {0, 0, 0};
        int tfNext = sl[32 * s0], jaNext = __shfl_sync(0xffffffffu, jr, (lane + s0) & 31);
#pragma unroll 2
        for (int s = s0; s < s1; ++s) {
          const int tf = tfNext, ja = jaNext;
          if (s + 1 < s1) tfNext = sl[32 * (s + 1)];
          jaNext = __shfl_sync(0xffffffffu, jr, (lane + s + 1) & 31);
          if (tf >= 0) {
            const bool sw = tf & kSwap;
            const GpuPair pr = sw ? GpuPair{ja, ia, tf & ~kSwap} : GpuPair{ia, ja, tf};
            T d[3], f;
            if (pairTermGpu(m, pr, pd, x, e, clashPairs, d, f) && g && f != 0)
              for (int c = 0; c < 3; ++c) {
                const T v = f * d[c];
                if (sw) {
                  gj[c] += v;
                  gi[c] -= v;
                } else {
                  gi[c] += v;
                  gj[c] -= v;
                }
              }
          }
          if (g && s + 1 < s1)
            for (int c = 0; c < 3; ++c) gj[c] = __shfl_sync(0xffffffffu, gj[c], (lane + 1) & 31);
        }
        if (!g) continue;
        const int jl = __shfl_sync(0xffffffffu, jr, (lane + s1 - 1) & 31);
        if (ia >= 0 && (gi[0] != 0 || gi[1] != 0 || gi[2] != 0))
          for (int c = 0; c < 3; ++c) atomicAdd(g + c * cc + ia, static_cast<double>(gi[c]));
        if (jl >= 0 && (gj[0] != 0 || gj[1] != 0 || gj[2] != 0))
          for (int c = 0; c < 3; ++c) atomicAdd(g + c * cc + jl, static_cast<double>(gj[c]));
      }
    } else {
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
    __shared__ double red[kWarps][kVals];
    {
      // Warp sums of the energy terms before the probe atoms (fewer registers there).
      double v[8] = {e[0], e[1], e[2], e[3], e[4], closure, maxDev, static_cast<double>(clashPairs)};
      for (int k = 0; k < 8; ++k) {
        v[k] = k == 6 ? warpMax(v[k]) : warpSum(v[k]);
        if (lane == 0) red[warp][k] = v[k];
      }
    }
    if constexpr (Dock) {
      // 8: inter, 9: interClash, 10: outside, 11–13: translation, 14–16: rotation gradient.
      double v[kVals] = {};
      const DockPose ps = *pose;
      // The field of this sample (spec 15.1.16): the kernel parameter, or of DockView::fields.
      if constexpr (Dock == 2)
        probeTerms<kThreads>(dv, dv.fields[ps.field], x, ps, g, cc, v);
      else
        probeTerms<kThreads>(dv, dv.field, x, ps, g, cc, v);
      for (int k = 8; k < kVals; ++k) {
        v[k] = warpSum(v[k]);
        if (lane == 0) red[warp][k] = v[k];
      }
    }
    __syncthreads();
    if (tid == 0) {
      for (int k = 0; k < kVals; ++k) {
        double r = red[0][k];
        for (int w = 1; w < kWarps; ++w) r = k == 6 ? fmax(r, red[w][k]) : r + red[w][k];
        res[k] = r;
      }
    }
    if constexpr (Tor) {
      // dE/dθ_u = k_u · Σ (x − a) × g = k_u · (Σ x × g − a × Σ g) over the coordinates moved by u:
      // the sums over the coordinates of each torsion, then over its subtree.
      const int nt = m.torsionCount;
      char* b = td.shared ? reinterpret_cast<char*>(x + nc + 2 * nt)
                          : td.buf + static_cast<std::size_t>(blockIdx.x) * (24 * sizeof(T) + 2 * sizeof(int)) * nt;
      double* S = reinterpret_cast<double*>((reinterpret_cast<std::uintptr_t>(b) + 7) & ~std::uintptr_t(7));
      for (int u = tid; g && u < nt; u += kThreads) {
        double s[6] = {0, 0, 0, 0, 0, 0};
        for (int q = td.ownOff[u]; q < td.ownOff[u + 1]; ++q) {
          const int i = td.own[q];
          const double f0 = g[i], f1 = g[cc + i], f2 = g[2 * cc + i];
          const double x0 = x[3 * i], x1 = x[3 * i + 1], x2 = x[3 * i + 2];
          s[0] += f0;
          s[1] += f1;
          s[2] += f2;
          s[3] += x1 * f2 - x2 * f1;
          s[4] += x2 * f0 - x0 * f2;
          s[5] += x0 * f1 - x1 * f0;
        }
        const int pu = td.pos[u];
        for (int k = 0; k < 6; ++k) S[k * nt + pu] = s[k];
      }
      __syncthreads();
      // Inclusive prefix sums of each component in the depth-first order (a warp per component).
      for (int c = warp; g && c < 6; c += kWarps) {
        double* Sc = S + c * nt;
        double carry = 0;
        for (int base = 0; base < nt; base += 32) {
          const int q = base + lane;
          double v = q < nt ? Sc[q] : 0.0;
          for (int o = 1; o < 32; o *= 2) {
            const double w = __shfl_up_sync(0xffffffffu, v, o);
            if (lane >= o) v += w;
          }
          v += carry;
          if (q < nt) Sc[q] = v;
          carry = __shfl_sync(0xffffffffu, v, 31);
        }
      }
      __syncthreads();
      for (int u = tid; g && u < nt; u += kThreads) {
        T a[3], k[3];
        double sum = 0;
        if (flat::torsionUnit(m, u, x, a, k)) {
          double s[6];
          const int p0 = td.pos[u], p1 = p0 + td.sub[u];
          for (int c = 0; c < 6; ++c) s[c] = S[c * nt + p1 - 1] - (p0 > 0 ? S[c * nt + p0 - 1] : 0.0);
          const double a0 = a[0], a1 = a[1], a2 = a[2];
          sum = k[0] * (s[3] - (a1 * s[2] - a2 * s[1])) + k[1] * (s[4] - (a2 * s[0] - a0 * s[2])) +
                k[2] * (s[5] - (a0 * s[1] - a1 * s[0]));
        }
        gradT[u] = sum;
      }
    }
    for (int t = warp; !Tor && g && t < m.torsionCount; t += kWarps) {
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
      if (lane == 0) gradT[t] = sum;
    }
    __syncthreads();
}

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
  constexpr int kThreads = NT;
  extern __shared__ __align__(16) unsigned char smem[];
  const int tid = threadIdx.x;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount);
  for (int s = blockIdx.x; s < n; s += gridDim.x) {
    const T* ang = angles + static_cast<std::size_t>(s) * m.torsionCount;
    T* x;
    double* g;
    x = shared ? reinterpret_cast<T*>(smem) : xs + s * nc;
    g = grad ? gs + blockIdx.x * nc : nullptr;

    constexpr int kVals = Dock ? 17 : 8;
    __shared__ double res[kVals];
    evalOne<T, Dock, NT, Tor>(m, pd, td, ang, x, g, shared, dv, Dock ? poses + s : nullptr, res,
                              grad ? grad + static_cast<std::size_t>(s) * m.torsionCount : nullptr);
    if (tid == 0) {
      FlatTerms out;
      out.vdw = res[0];
      out.hbond = res[1];
      out.elec = res[2];
      out.desolv = res[3];
      out.clash = res[4];
      out.closure = res[5];
      out.closureMaxDeviation = res[6];
      out.clashPairs = static_cast<int>(res[7]);
      terms[s] = out;
      if constexpr (Dock) {
        dterms[s] = {res[8], res[9], static_cast<int>(res[10])};
        if (grad6)
          for (int i = 0; i < 6; ++i) grad6[static_cast<std::size_t>(s) * 6 + i] = res[11 + i];
      }
    }
    if (shared && keepX) {
      for (std::size_t i = tid; i < nc; i += kThreads) xs[s * nc + i] = x[i];
    }
    __syncthreads();
  }
}


// Double-precision quaternion helpers of the device minimizer (as orient.cpp; w, x, y, z).
__device__ inline void quatMul(const double a[4], const double b[4], double o[4]) {
  o[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  o[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2];
  o[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1];
  o[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0];
}
__device__ inline void quatNormalize(double q[4]) {
  const double n = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  for (int i = 0; i < 4; ++i) q[i] = q[i] / n;
}
__device__ inline void quatFromRotationVector(const double v[3], double q[4]) {
  const double a = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (a < 1e-12) {
    q[0] = 1;
    for (int i = 0; i < 3; ++i) q[1 + i] = 0.5 * v[i];
    quatNormalize(q);
    return;
  }
  const double s = sin(0.5 * a) / a;
  q[0] = cos(0.5 * a);
  for (int i = 0; i < 3; ++i) q[1 + i] = s * v[i];
}
// wrapDeg (rng.hpp).
__device__ inline double wrapDegDev(double a) {
  a = a - 360.0 * static_cast<double>(static_cast<long long>(a / 360.0));
  if (a > 180.0) a -= 360.0;
  if (a <= -180.0) a += 360.0;
  return a;
}
// Sum over a warp, received by all lanes.
__device__ inline double warpAllSum(double v) {
  for (int o = 16; o > 0; o /= 2) v += __shfl_xor_sync(0xffffffffu, v, o);
  return v;
}

// Work memory of minimizeDockKernel per block (P = memory + 1 stored pairs, K = 2P + 1 basis
// vectors): doubles S and Y (P × ng each), the vectors (minVecDoubles: z, g, d, gt with ng each and
// the torsion gradient, nt) and the scalars (minScalarDoubles: the dot products of the basis vectors
// G (packed lower triangle, K (K + 1) / 2), rho (P), alpha (memory), the coefficients (K) and the
// order of the pairs (P, int)); floats (fwork): trial and best angles (nt each). The scalars are in
// shared memory when they fit (DockMinArgs::scalarOff). The other methods keep their vectors in
// the memory of S and Y (cg: g and d of the previous step; fire: v; adadelta: E[g²] and E[Δx²]).
__host__ __device__ inline std::size_t minVecDoubles(int ng, int nt) {
  return 4 * static_cast<std::size_t>(ng) + static_cast<std::size_t>(nt);
}
__host__ __device__ inline std::size_t minScalarDoubles(int mem) {
  const std::size_t P = static_cast<std::size_t>(mem) + 1, K = 2 * P + 1;
  return K * (K + 1) / 2 + P + static_cast<std::size_t>(mem) + K + (P + 1) / 2;
}
// Index of the basis vector pair (u, v) in the packed G.
__device__ inline int triIndex(int u, int v) { return u >= v ? u * (u + 1) / 2 + v : v * (v + 1) / 2 + u; }
__host__ __device__ inline std::size_t minWorkDoubles(int ng, int nt, int mem) {
  return 2 * (static_cast<std::size_t>(mem) + 1) * static_cast<std::size_t>(ng) + minVecDoubles(ng, nt) +
         minScalarDoubles(mem);
}
__host__ __device__ inline std::size_t minWorkWords(int nt) { return 2 * static_cast<std::size_t>(nt); }

// Requests of stealDock in mapped host memory (spec 15.1.35): req and the parameters written by the
// host, the reply (seq = req, newEnd, oldEnd) and taken by the kernel.
struct DockSteal {
  volatile int req, giveMin;
  volatile float frac;
  volatile int reply, newEnd, oldEnd, taken;
};

struct DockMinArgs {
  DockMinParams p;
  const DockMinStart* starts;
  const double* anglesIn;  // n × nt, degrees
  DockMinResult* out;
  float* anglesOut;  // n × nt: angles of the result
  double* work;
  float* fwork;
  // Sample counter: low 32 bits the next sample, high 32 bits the end (the samples from the end on
  // were given away, spec 15.1.35).
  unsigned long long* next;
  DockSteal* steal;  // mapped host memory (null: no stealing)
  int* stealSeen;    // the last request handled
  int n;
  int scalarOff;  // >= 0: the scalars (minScalarDoubles) in dynamic shared memory from this byte offset
};

// Whole 9b minimizations (spec 15.1.26–15.1.28): MobmDocker::minimize with the five methods of
// minimize.cpp (runSd, runCg, runFire, runLbfgs, runAdadelta) on the device. Each block takes the
// next sample until none is left; all threads evaluate the energy (evalOne) and do the per-element
// work (trial points, gradients, stored vectors); warp 0 computes each new step (lane-parallel).
// The L-BFGS direction is computed from the dot products of the basis vectors (the stored pairs s,
// y and the gradient; vector-free L-BFGS): after each accepted step the warps compute the new dot
// products together, then warp 0 runs the two-loop recursion on the coefficients of the basis (no
// vector memory) and forms the direction in one pass. Steps of the minimizers as on the CPU; the
// sums are in another order.
template <int Dock, int NT, bool Tor>
__global__ void __launch_bounds__(NT) minimizeDockKernel(EnergyModelView m, PairData<float> pd, TorsionData<float> td,
                                                         float* xs, double* gs, DockView dv, DockMinArgs a, int shared) {
  extern __shared__ __align__(16) unsigned char smem[];
  constexpr double kRad = 3.14159265358979323846 / 180.0, kArmijo = 1e-4;
  constexpr int kWarps = NT / 32, kC = (2 * (kDeviceLbfgsMaxMemory + 1) + 1 + 31) / 32, kD = 4;
  constexpr int kSd = static_cast<int>(MinMethod::Sd), kCg = static_cast<int>(MinMethod::Cg),
                kFire = static_cast<int>(MinMethod::Fire), kLbfgs = static_cast<int>(MinMethod::Lbfgs),
                kAda = static_cast<int>(MinMethod::Adadelta);
  const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
  const bool w0 = tid < 32;
  const int nt = m.torsionCount, ng = 6 + nt, M = a.p.memory, P = M + 1, K = 2 * P + 1, gI = 2 * P;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount), sng = static_cast<std::size_t>(ng);
  float* x = shared ? reinterpret_cast<float*>(smem) : xs + blockIdx.x * nc;
  double* g = gs + blockIdx.x * nc;
  const DockMinParams& p = a.p;
  const int meth = p.method;
  // The work memory pointers are read from shared memory where they are used (volatile: not kept
  // in registers across the energy evaluation, where they made the energy loops spill).
  __shared__ double *sWork, *sVec, *sScalar;
  __shared__ float *sAngT, *sAngB;
  // order[0, cnt): stored pairs, oldest first; order[cnt]: the next one.
#define OD_MINW                                                                                                       \
  double* const S = *reinterpret_cast<double* volatile*>(&sWork);                                                    \
  double* const Y = S + P * sng;                                                                                     \
  double* const z = *reinterpret_cast<double* volatile*>(&sVec);                                                     \
  double *const gv = z + ng, *const d = gv + ng, *const gt = d + ng, *const gradT = gt + ng;                         \
  double* const G = *reinterpret_cast<double* volatile*>(&sScalar);                                                  \
  double *const rho = G + static_cast<std::size_t>(K) * (K + 1) / 2, *const alpha = rho + P;                        \
  double* const coef = alpha + M;                                                                                    \
  int* const order = reinterpret_cast<int*>(coef + K);                                                               \
  float* const angT = *reinterpret_cast<float* volatile*>(&sAngT);                                                  \
  float* const angB = *reinterpret_cast<float* volatile*>(&sAngB);                                                  \
  const auto basis = [&](int l) -> const double* { return l < P ? S + l * sng : l < gI ? Y + (l - P) * sng : gv; }; \
  (void)gv, (void)d, (void)gt, (void)G, (void)rho, (void)alpha, (void)coef, (void)gradT, (void)angB,       \
      (void)order, (void)basis;

  __shared__ DockPose pose, bestPose;
  __shared__ double res[17], bestRes[17], q0[4];
  __shared__ double e, eT, bestE, slope, dScale, fTrial;
  __shared__ int sIdx, status, steps, evals, calm, cnt, accepted, pending;  // status: 0 running, 1 done, 2 left the grid
  // sd: step factor, rejected steps in a row; cg: alpha and slope of the previous step, slope of the
  // current direction, no previous step; fire: dt, alpha, steps with P > 0 in a row.
  __shared__ double sdA, cgAlpha, cgPrevSlope, cgSlope, fireDt, fireAlpha;
  __shared__ int sdFails, cgFirst, fireNpos;

  // Trial point v = z + f d (pose and angles; MobmDocker::minimize unpack); all threads.
  const auto trial = [&](double f) {
    OD_MINW
    if (tid == 0) {
      fTrial = f;
      double v[6];
      for (int i = 0; i < 6; ++i) v[i] = z[i] + f * d[i];
      const double rv[3] = {v[3] * p.sR, v[4] * p.sR, v[5] * p.sR};
      double r[4], q[4];
      quatFromRotationVector(rv, r);
      quatMul(r, q0, q);
      quatNormalize(q);
      for (int i = 0; i < 3; ++i) pose.t[i] = static_cast<float>(v[i] * p.sT);
      for (int i = 0; i < 4; ++i) pose.q[i] = static_cast<float>(q[i]);
    }
    for (int k = tid; k < nt; k += NT) {
      const double v = z[6 + k] + f * d[6 + k];
      angT[k] = static_cast<float>(wrapDegDev(v * p.sA / kRad));
    }
  };
  // Energy and scaled gradient (gout) of pose and angT; all threads; true unless outside the grid.
  const auto evaluate = [&](int which, double& eout) {
    __syncthreads();
    {
      OD_MINW
      evalOne<float, Dock, NT, Tor>(m, pd, td, angT, x, g, shared, dv, &pose, res, gradT);
    }
    OD_MINW
    double* const gout = which ? gt : gv;
    for (int i = tid; i < ng; i += NT) gout[i] = i < 3 ? res[11 + i] * p.sT : i < 6 ? res[11 + i] * p.sR : gradT[i - 6] * p.sA;
    if (tid == 0) {
      ++evals;
      eout = res[8] + res[9] + (((res[0] + res[1] + res[2] + res[3]) + res[4]) + res[5]);
    }
    return res[10] == 0;
  };
  // The current trial point becomes the best one: angles and terms (all threads), then pose and
  // energy (bestPoseE, thread 0, after a barrier when other threads read bestE).
  const auto keepBest = [&]() {
    OD_MINW
    if (tid < 17) bestRes[tid] = res[tid];
    for (int k = tid; k < nt; k += NT) angB[k] = angT[k];
  };
  // New dot products (all warps): with a new pair (the pair order[cnt] against the stored pairs,
  // itself and the gradient) and of the gradient against the stored pairs; only g.g otherwise.
  const auto dots = [&](bool withNew) {
    OD_MINW
    const int c = cnt, q = order[c], total = withNew ? 6 * c + 6 : 1 + 2 * c;
    for (int k = warp; k < total; k += kWarps) {
      int u, v;
      if (!withNew) {
        if (k == 0) u = v = gI;
        else {
          u = gI;
          v = (k - 1) % 2 ? P + order[(k - 1) / 2] : order[(k - 1) / 2];
        }
      } else if (k < 6 * c) {
        const int pp = order[k / 6], kind = k % 6;
        u = kind < 2 ? q : kind < 4 ? P + q : gI;
        v = kind % 2 ? P + pp : pp;
      } else {
        constexpr int ku[6] = {0, 0, 1, 0, 1, 2}, kv[6] = {0, 1, 1, 2, 2, 2};
        const int kk = k - 6 * c, id[3] = {q, P + q, gI};
        u = id[ku[kk]];
        v = id[kv[kk]];
      }
      const double* bu = basis(u);
      const double* bv = basis(v);
      double sum = 0;
      for (int i = lane; i < ng; i += 32) sum += bu[i] * bv[i];
      sum = warpAllSum(sum);
      if (lane == 0) {
        G[triIndex(u, v)] = sum;
      }
    }
  };
  // Sum over the basis of coefficient × G row (warp 0, lane-parallel).
  const auto rowDot = [&](const double (&cr)[kC], int row) {
    OD_MINW
    double s = 0;
#pragma unroll
    for (int r = 0; r < kC; ++r) {
      const int l = lane + 32 * r;
      if (l < K && cr[r] != 0) s += cr[r] * G[triIndex(row, l)];
    }
    return warpAllSum(s);
  };
  const auto addCoef = [&](double (&cr)[kC], int l, double v) {
#pragma unroll
    for (int r = 0; r < kC; ++r)
      if (l == lane + 32 * r) cr[r] += v;
  };

  if (tid == 0) {
    sWork = a.work + blockIdx.x * minWorkDoubles(ng, nt, M);
    sVec = sWork + 2 * P * sng;
    sScalar = a.scalarOff >= 0 ? reinterpret_cast<double*>(smem + a.scalarOff) : sVec + minVecDoubles(ng, nt);
    sAngB = a.fwork + blockIdx.x * minWorkWords(nt);
    sAngT = sAngB + nt;
  }
  for (;;) {
    if (tid == 0) {
      if (a.steal) {
        const int r = a.steal->req;
        const int seen = *reinterpret_cast<volatile int*>(a.stealSeen);
        if (r != seen && atomicCAS(a.stealSeen, seen, r) == seen) {
          // Gives away the last part of the untaken samples.
          unsigned long long old = *reinterpret_cast<volatile unsigned long long*>(a.next);
          int ne, oe;
          for (;;) {
            const int s0 = static_cast<int>(old & 0xffffffffULL);
            oe = static_cast<int>(old >> 32);
            const int give = oe > s0 ? static_cast<int>(static_cast<float>(oe - s0) * a.steal->frac) : 0;
            ne = give >= a.steal->giveMin && give > 0 ? oe - give : oe;
            if (ne == oe) break;
            const unsigned long long nw = (static_cast<unsigned long long>(ne) << 32) | (old & 0xffffffffULL);
            const unsigned long long got = atomicCAS(a.next, old, nw);
            if (got == old) break;
            old = got;
          }
          a.steal->newEnd = ne;
          a.steal->oldEnd = oe;
          __threadfence_system();
          a.steal->reply = r;
        }
      }
      const unsigned long long old = atomicAdd(a.next, 1ULL);
      const int s0 = static_cast<int>(old & 0xffffffffULL), e0 = static_cast<int>(old >> 32);
      sIdx = s0 < e0 ? s0 : a.n;
      if (a.steal && s0 < e0) a.steal->taken = s0 + 1;
    }
    __syncthreads();
    const int s = sIdx;
    if (s >= a.n) return;
    const DockMinStart& st = a.starts[s];
    {
    OD_MINW
    for (int i = tid; i < ng; i += NT)
      z[i] = i < 3 ? st.t[i] / p.sT : i < 6 ? 0.0 : a.anglesIn[static_cast<std::size_t>(s) * nt + (i - 6)] * kRad / p.sA;
    for (int i = tid; i < ng; i += NT) d[i] = 0;
    for (int k = tid; k < P; k += NT) order[k] = k;
    if (meth != kLbfgs)
      for (int i = tid; i < ng; i += NT) S[i] = Y[i] = 0;
    if (tid == 0) {
      sdA = 1;
      sdFails = 0;
      cgAlpha = 1;
      cgPrevSlope = 0;
      cgFirst = 1;
      fireDt = 0.1;
      fireAlpha = 0.1;
      fireNpos = 0;
      for (int i = 0; i < 4; ++i) q0[i] = st.q[i];
      for (int i = 0; i < 3; ++i) pose.center[i] = st.center[i];
      pose.field = st.field;
      status = steps = evals = calm = cnt = pending = 0;
    }
    }
    __syncthreads();
    // One evaluation per pass (a single copy of the energy code): the start point, then the trial
    // points of the line searches; a new direction after the start point and each accepted step.
    bool first = true;
    double al = 1, prevAl = 1;  // prevAl: factor of the last accepted step
    int bt = 0;
    for (;;) {
      trial(first ? 0.0 : al * dScale);
      const bool inside = evaluate(first ? 0 : 1, first ? e : eT);
      if (first) {
        if (!inside) {
          if (tid == 0) status = 2;
        } else {
          keepBest();
          if (tid == 0) {
            bestPose = pose;
            bestE = e;
          }
          if (meth == kLbfgs) dots(false);
        }
        __syncthreads();
        first = false;
      } else {
        // Armijo rule (sd: one trial per step; cg, lbfgs: backtracking); fire and adadelta accept every step.
        if (tid == 0) {
          if (!inside) status = 2;
          accepted = inside && (meth == kFire || meth == kAda || eT <= e + kArmijo * al * slope);
        }
        __syncthreads();
        if (status != 0) break;
        if (!accepted && meth != kSd) {
          if (++bt < p.maxBacktrack) {
            if (p.lineSearchReduce == static_cast<int>(LineSearchReduce::Interp)) {
              const double den = 2 * (eT - e - slope * al);
              const double an = den > 0 ? -slope * al * al / den : 0.5 * al;
              al = fmin(fmax(an, 0.1 * al), 0.5 * al);
            } else {
              al *= 0.5;
            }
            continue;
          }
          break;  // no decrease: converged
        }
        if (!accepted) {  // sd: half the step factor, the same point
          if (tid == 0) {
            sdA *= 0.5;
            if (++sdFails >= p.maxBacktrack) status = 1;
          }
          __syncthreads();
          if (status != 0) break;
        } else {
          // Accept (Run::accept, onAccept): lbfgs the pair order[cnt], cg the previous g and d;
          // re-base the rotation at the accepted pose.
          OD_MINW
          const bool better = eT < bestE;
          prevAl = al;
          const int q = meth == kLbfgs ? order[cnt] : 0;
          double* sq = S + q * sng;
          double* yq = Y + q * sng;
          for (int i = tid; i < ng; i += NT) {
            const double v = z[i] + fTrial * d[i];  // the accepted trial point
            if (meth == kLbfgs) {
              sq[i] = v - z[i];
              yq[i] = gt[i] - gv[i];
            } else if (meth == kCg) {
              S[i] = gv[i];
              Y[i] = d[i];
            }
            z[i] = i >= 3 && i < 6 ? 0.0 : v;
            gv[i] = gt[i];
          }
          if (better) keepBest();
          __syncthreads();
          if (tid == 0) {
            if (better) {
              bestPose = pose;
              bestE = eT;
            }
            ++steps;
            calm = fabs(eT - e) < p.deltaE ? calm + 1 : 0;
            if (calm >= p.convSteps) status = 1;
            e = eT;
            for (int i = 0; i < 4; ++i) q0[i] = pose.q[i];
            pending = meth == kLbfgs;
            if (meth == kSd) {
              sdFails = 0;
              sdA *= 1.2;
            } else if (meth == kCg) {
              cgAlpha = fTrial;
              cgPrevSlope = cgSlope;
              cgFirst = 0;
            }
          }
          if (meth == kLbfgs) dots(true);
        }
        __syncthreads();
      }
      if (status != 0 || steps >= p.maxSteps) break;
      if (w0 && meth != kLbfgs) {
        OD_MINW
        // Largest factor (<= 1) of the step f d (MobmDocker::minimize stepLimit).
        const auto limitOf = [&](double f) {
          double lim = 1, tq = 0, rq = 0;
          for (int i = lane; i < ng; i += 32) {
            const double v = f * d[i];
            if (i < 3) tq += v * v;
            else if (i < 6) rq += v * v;
            else {
              const double av = p.sA * fabs(v);
              if (av > p.maxAngle) lim = fmin(lim, p.maxAngle / av);
            }
          }
          for (int o = 16; o > 0; o /= 2) lim = fmin(lim, __shfl_xor_sync(0xffffffffu, lim, o));
          tq = warpAllSum(tq);
          rq = warpAllSum(rq);
          const double tn = p.sT * sqrt(tq);
          if (tn > p.maxTrans) lim = fmin(lim, p.maxTrans / tn);
          const double rn = p.sR * sqrt(rq);
          if (rn > p.maxAngle) lim = fmin(lim, p.maxAngle / rn);
          return lim;
        };
        const auto wdot = [&](const double* u, const double* v) {
          double s = 0;
          for (int i = lane; i < ng; i += 32) s += u[i] * v[i];
          return warpAllSum(s);
        };
        if (meth == kSd) {  // runSd
          const double a0 = sdA;
          for (int i = lane; i < ng; i += 32) d[i] = -a0 * gv[i];
          __syncwarp();
          const double s = limitOf(1.0);
          double sl = 0;
          for (int i = lane; i < ng; i += 32) {
            d[i] *= s;
            sl += gv[i] * d[i];
          }
          sl = warpAllSum(sl);
          if (lane == 0) {
            sdA = a0 * s;
            dScale = 1;
            slope = sl;
            if (!(sl < 0)) status = 1;  // zero gradient
          }
        } else if (meth == kCg) {  // runCg: S previous g, Y previous d
          const bool firstStep = cgFirst;
          if (firstStep) {
            for (int i = lane; i < ng; i += 32) d[i] = -gv[i];
          } else {
            double num = 0, den = 0;
            for (int i = lane; i < ng; i += 32) {
              num += gv[i] * (gv[i] - S[i]);
              den += S[i] * S[i];
            }
            const double beta = fmax(0.0, warpAllSum(num) / warpAllSum(den));
            for (int i = lane; i < ng; i += 32) d[i] = -gv[i] + beta * Y[i];
          }
          __syncwarp();
          double sl = wdot(gv, d);
          if (!(sl < 0)) {  // not a descent direction: restart
            for (int i = lane; i < ng; i += 32) d[i] = -gv[i];
            __syncwarp();
            sl = wdot(gv, d);
            if (!(sl < 0) && lane == 0) status = 1;
          }
          const double a0 = firstStep ? 1.0 : cgAlpha * cgPrevSlope / sl;
          const double s = limitOf(a0);
          if (lane == 0) {
            dScale = a0 * s;
            slope = a0 * s * sl;
            cgSlope = sl;
          }
        } else if (meth == kFire) {  // runFire: S velocity
          constexpr int kNmin = 5;
          constexpr double kFinc = 1.1, kFdec = 0.5, kAlpha0 = 0.1, kFalpha = 0.99, kDtMax = 1.0;
          double* v = S;
          double P = 0, vn = 0, fn = 0;
          for (int i = lane; i < ng; i += 32) {
            P -= gv[i] * v[i];
            vn += v[i] * v[i];
            fn += gv[i] * gv[i];
          }
          P = warpAllSum(P);
          vn = warpAllSum(vn);
          fn = warpAllSum(fn);
          double dt = fireDt, fa = fireAlpha;
          int np = fireNpos;
          if (P > 0) {
            vn = sqrt(vn);
            fn = sqrt(fn);
            if (fn > 0)
              for (int i = lane; i < ng; i += 32) v[i] = (1 - fa) * v[i] - fa * vn * gv[i] / fn;
            if (++np > kNmin) {
              dt = fmin(dt * kFinc, kDtMax);
              fa *= kFalpha;
            }
          } else {
            for (int i = lane; i < ng; i += 32) v[i] = 0;
            dt *= kFdec;
            fa = kAlpha0;
            np = 0;
          }
          for (int i = lane; i < ng; i += 32) {
            v[i] -= dt * gv[i];
            d[i] = dt * v[i];
          }
          __syncwarp();
          const double s = limitOf(1.0);
          if (lane == 0) {
            fireDt = dt;
            fireAlpha = fa;
            fireNpos = np;
            dScale = s;
          }
        } else {  // runAdadelta: S E[g²], Y E[Δx²]
          constexpr double kRho = 0.8, kEps = 1e-2;
          for (int i = lane; i < ng; i += 32) {
            S[i] = kRho * S[i] + (1 - kRho) * gv[i] * gv[i];
            d[i] = -gv[i] * sqrt((Y[i] + kEps) / (S[i] + kEps));
          }
          __syncwarp();
          const double s = limitOf(1.0);
          for (int i = lane; i < ng; i += 32) {
            const double st = d[i] * s;
            Y[i] = kRho * Y[i] + (1 - kRho) * st * st;
          }
          if (lane == 0) dScale = s;
        }
      } else if (w0) {
        OD_MINW
        // The new pair is stored when s.y > 1e-12 (the oldest one is dropped at memory pairs).
        if (pending) {
          const int c = cnt, q = order[c];
          const double sy = G[triIndex(q, P + q)];
          __syncwarp();
          if (sy > 1e-12 && lane == 0) {
            rho[q] = 1 / sy;
            if (c == M) {
              const int old = order[0];
              for (int k = 0; k + 1 < M; ++k) order[k] = order[k + 1];
              order[M - 1] = q;
              order[M] = old;
            } else {
              cnt = c + 1;
            }
          }
          __syncwarp();
        }
        // Direction: L-BFGS two-loop recursion on the coefficients cr of the basis, initial
        // Hessian (s.y / y.y) I.
        double cr[kC];
#pragma unroll
        for (int r = 0; r < kC; ++r) cr[r] = lane + 32 * r == gI ? 1.0 : 0.0;
        const int c = cnt;
        for (int k = c; k-- > 0;) {
          const int j = order[k];
          const double al = rho[j] * rowDot(cr, j);
          if (lane == 0) alpha[k] = al;
          addCoef(cr, P + j, -al);
        }
        if (c > 0) {
          const int j = order[c - 1];
          const double gamma = G[triIndex(j, P + j)] / G[triIndex(P + j, P + j)];
#pragma unroll
          for (int r = 0; r < kC; ++r) cr[r] *= gamma;
        }
        __syncwarp();
        for (int k = 0; k < c; ++k) {
          const int j = order[k];
          const double beta = rho[j] * rowDot(cr, P + j);
          addCoef(cr, j, alpha[k] - beta);
        }
        // d = -sum cr b; slope g.d from G.
        double sl = -rowDot(cr, gI);
        if (!(sl < 0)) {  // not a descent direction: restart from steepest descent
#pragma unroll
          for (int r = 0; r < kC; ++r) cr[r] = lane + 32 * r == gI ? 1.0 : 0.0;
          sl = -G[triIndex(gI, gI)];
          if (lane == 0) cnt = 0;
          if (!(sl < 0) && lane == 0) status = 1;  // zero gradient
        }
#pragma unroll
        for (int r = 0; r < kC; ++r)
          if (lane + 32 * r < K) coef[lane + 32 * r] = cr[r];
        __syncwarp();
        const int cu = cnt;
        // Step limit (MobmDocker::minimize stepLimit) while forming d.
        double lim = 1, tq = 0, rq = 0;
        for (int base = 0; base < ng; base += 32 * kD) {
          double acc[kD] = {};
          for (int b = -1; b < 2 * cu; ++b) {
            const int l = b < 0 ? gI : b % 2 ? P + order[b / 2] : order[b / 2];
            const double cl = coef[l];
            const double* bl = basis(l);
#pragma unroll
            for (int r = 0; r < kD; ++r) {
              const int i = base + lane + 32 * r;
              if (i < ng) acc[r] += cl * bl[i];
            }
          }
#pragma unroll
          for (int r = 0; r < kD; ++r) {
            const int i = base + lane + 32 * r;
            if (i >= ng) continue;
            const double v = -acc[r];
            d[i] = v;
            if (i < 3) tq += v * v;
            else if (i < 6) rq += v * v;
            else {
              const double av = p.sA * fabs(v);
              if (av > p.maxAngle) lim = fmin(lim, p.maxAngle / av);
            }
          }
        }
        for (int o = 16; o > 0; o /= 2) lim = fmin(lim, __shfl_xor_sync(0xffffffffu, lim, o));
        tq = warpAllSum(tq);
        rq = warpAllSum(rq);
        const double tn = p.sT * sqrt(tq);
        if (tn > p.maxTrans) lim = fmin(lim, p.maxTrans / tn);
        const double rn = p.sR * sqrt(rq);
        if (rn > p.maxAngle) lim = fmin(lim, p.maxAngle / rn);
        if (lane == 0) {
          dScale = lim;
          slope = sl * lim;
          pending = 0;
        }
      }
      __syncthreads();
      if (status != 0) break;
      al = meth == kLbfgs && p.lineSearchInit == static_cast<int>(LineSearchInit::Prev) ? fmin(1.0, 2 * prevAl) : 1.0;
      bt = 0;
    }
    if (tid == 0) {
      DockMinResult r{};
      r.steps = steps;
      r.evaluations = evals;
      r.ok = status != 2;
      if (r.ok) {
        for (int i = 0; i < 3; ++i) r.t[i] = bestPose.t[i];
        for (int i = 0; i < 4; ++i) r.q[i] = bestPose.q[i];
        r.terms = {bestRes[0], bestRes[1], bestRes[2], bestRes[3], bestRes[4], bestRes[5], bestRes[6],
                   static_cast<int>(bestRes[7])};
        r.dterms = {bestRes[8], bestRes[9], static_cast<int>(bestRes[10])};
      }
      a.out[s] = r;
    }
    if (status != 2) {
      OD_MINW
      for (int k = tid; k < nt; k += NT) a.anglesOut[static_cast<std::size_t>(s) * nt + k] = angB[k];
    }
    __syncthreads();
  }
#undef OD_MINW
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
  DevBuf<float> maps, packed;
  const float* ptr = nullptr;                 // maps.p, or the device maps of the field (adopted)
  std::shared_ptr<const DeviceBlocks> adopted;
};

// FieldView::packed from the maps (nm maps of n[0] × n[1] × n[2] points).
__global__ void packMapsKernel(const float* maps, float4* packed, int nm, int n0, int n1, int n2) {
  const std::size_t points = static_cast<std::size_t>(n0) * n1 * n2;
  for (std::size_t g = blockIdx.x * static_cast<std::size_t>(blockDim.x) + threadIdx.x; g < points * nm;
       g += static_cast<std::size_t>(gridDim.x) * blockDim.x) {
    const std::size_t m = g / points, pt = g % points;
    const int i = static_cast<int>(pt % n0), j = static_cast<int>(pt / n0 % n1), k = static_cast<int>(pt / n0 / n1);
    (void)i;
    const float* mp = maps + m * points;
    const std::size_t dy = n0, dz = static_cast<std::size_t>(n0) * n1;
    const bool y1 = j + 1 < n1, z1 = k + 1 < n2;
    packed[g] = make_float4(mp[pt], y1 ? mp[pt + dy] : 0.0f, z1 ? mp[pt + dz] : 0.0f, y1 && z1 ? mp[pt + dy + dz] : 0.0f);
  }
}

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
  Packed packed = Packed::Off;
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
void GpuField::setPackedMaps(Packed p) { impl_->packed = p; }
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
      b += (s->maps->maps.n + s->maps->packed.n + (s->maps->adopted ? s->maps->adopted->floats : 0)) * sizeof(float);
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
    if (fm.deviceMaps && fm.deviceMaps->device == d.device) {  // built here (spec 15.1.34)
      maps->adopted = fm.deviceMaps;
      maps->ptr = fm.deviceMaps->ptr;
    } else {
      maps->maps.reserve(static_cast<std::size_t>(nm) * fm.points);
      maps->ptr = maps->maps.p;
      if (fm.deviceMaps)
        check(cudaMemcpyPeer(maps->maps.p, d.device, fm.deviceMaps->ptr, fm.deviceMaps->device,
                             static_cast<std::size_t>(nm) * fm.points * sizeof(float)),
              "cudaMemcpyPeer maps");
      else
        for (int m = 0; m < nm; ++m)
          check(cudaMemcpy(maps->maps.p + static_cast<std::size_t>(m) * fm.points, fm.view.maps[m], fm.points * sizeof(float),
                           cudaMemcpyHostToDevice),
                "upload map");
    }
    maps->base = fm.base;
    maps->baseDist = fm.baseDist;
    const std::size_t packedBytes = 16 * static_cast<std::size_t>(nm) * fm.points;
    if (d.packed == Packed::On || (d.packed == Packed::IfFits && packedBytes <= gpuFreeMem(d.device) / 4)) {
      maps->packed.reserve(4 * static_cast<std::size_t>(nm) * fm.points);
      packMapsKernel<<<1024, 256>>>(maps->ptr, reinterpret_cast<float4*>(maps->packed.p), nm, fm.view.n[0], fm.view.n[1],
                                    fm.view.n[2]);
      check(cudaDeviceSynchronize(), "packMapsKernel");
    }
  }
  auto sp = std::make_unique<FieldSlot>();
  FieldSlot& s = *sp;
  s.maps = std::move(maps);
  // Delta blocks already in device memory (built there) are used in place.
  if (fm.deviceData && fm.deviceData->device == d.device) s.deviceData = fm.deviceData;
  else if (fm.deviceData && std::any_of(fm.deviceData->others.begin(), fm.deviceData->others.end(),
                                        [&](const auto& o) { return o->device == d.device; }))
    s.deviceData = *std::find_if(fm.deviceData->others.begin(), fm.deviceData->others.end(),
                                 [&](const auto& o) { return o->device == d.device; });
  else if (fm.data) s.data.upload(*fm.data);
  else if (fm.deviceData) {  // built on another device (spec 15.1.33)
    s.data.reserve(std::max<std::size_t>(1, fm.deviceData->floats));
    if (fm.deviceData->floats)
      check(cudaMemcpyPeer(s.data.p, d.device, fm.deviceData->ptr, fm.deviceData->device, fm.deviceData->floats * sizeof(float)),
            "cudaMemcpyPeer delta blocks");
  }
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
  for (int m = 0; m < nm; ++m) v.maps[m] = s.maps->ptr + static_cast<std::size_t>(m) * fm.points;
  v.packed = s.maps->packed.p;
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
  // Pair tiles for float with more samples than multiprocessors (spec 15.1.24).
  PairData<T> pdl = pd;
  pdl.tiled = std::is_same_v<T, float> && n > sms ? 1 : 0;
  (tdl.ok ? kernelT : kernelS)<<<blocks, NT, shm, st>>>(m, pdl, tdl, angles, xs, gs, terms, grad, dv, poses, dterms, grad6,
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
  DevBuf<int> tileI, tileJ, tileS0, tileOff, slots;
  DevBuf<GpuPair> gpuPairs;
  int tiles = 0;
  DevBuf<GpuAtom<float>> atomsF;
  DevBuf<GpuAtom<double>> atomsD;
  DevBuf<GpuPairParams<float>> paramsF;
  DevBuf<GpuPairParams<double>> paramsD;
  Lane lanes[GpuEnergy::kLanes];
  DevBuf<float> axisF;
  DevBuf<double> axisD;
  DevBuf<int> tParent, tDeep, tAux;
  int tRounds = 0, tAuxOff[4] = {};
  bool tOk = false;

  template <class T>
  TorsionData<T> torsionData(Lane& L, std::size_t n) {
    if (!tOk) return {nullptr, nullptr, nullptr, nullptr, 0, 0, 0};
    L.tbuf.reserve(std::max<std::size_t>(1, n * (24 * sizeof(T) + 2 * sizeof(int)) * static_cast<std::size_t>(view.torsionCount)));
    TorsionData<T> td{nullptr, tParent.p, tDeep.p, L.tbuf.p, tRounds, 1, 0};
    if constexpr (std::is_same_v<T, float>) td.axis = axisF.p;
    else td.axis = axisD.p;
    setGradientPlan(td);
    return td;
  }
  template <class T>
  void setGradientPlan(TorsionData<T>& td) const {
    td.ownOff = tAux.p + tAuxOff[0];
    td.own = tAux.p + tAuxOff[1];
    td.pos = tAux.p + tAuxOff[2];
    td.sub = tAux.p + tAuxOff[3];
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
    // Torsion gradient plan (TorsionData::own, pos, sub); parents have lower indices.
    {
      std::vector<std::vector<int>> own(nt), child(nt);
      for (int a = 0; a < nc; ++a)
        if (deep[a] >= 0) own[deep[a]].push_back(a);
      for (int u = 0; u < nt; ++u)
        if (parent[u] >= 0) child[parent[u]].push_back(u);
      std::vector<int> aux, pos(nt), sub(nt, 1);
      for (int u = nt; u-- > 0;)
        if (parent[u] >= 0) sub[parent[u]] += sub[u];
      int next = 0;
      std::vector<int> stack;
      for (int u = nt; u-- > 0;)
        if (parent[u] < 0) stack.push_back(u);
      while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        pos[u] = next++;
        for (auto it = child[u].rbegin(); it != child[u].rend(); ++it) stack.push_back(*it);
      }
      tAuxOff[0] = 0;
      int o = 0;
      for (int u = 0; u < nt; ++u) {
        aux.push_back(o);
        o += static_cast<int>(own[u].size());
      }
      aux.push_back(o);
      tAuxOff[1] = static_cast<int>(aux.size());
      for (const auto& e : own) aux.insert(aux.end(), e.begin(), e.end());
      tAuxOff[2] = static_cast<int>(aux.size());
      aux.insert(aux.end(), pos.begin(), pos.end());
      tAuxOff[3] = static_cast<int>(aux.size());
      aux.insert(aux.end(), sub.begin(), sub.end());
      tAux.upload(aux);
    }
    tRounds = 0;
    while ((1 << tRounds) < maxDepth) ++tRounds;
    tOk = true;
  }

  template <class T>
  PairData<T> pairData() const {
    if constexpr (std::is_same_v<T, float>)
      return {atomsF.p, paramsF.p, tileI.p, tileJ.p, tileS0.p, tileOff.p, slots.p, tiles, gpuPairs.p, view.pairCount, 0};
    else
      return {atomsD.p, paramsD.p, tileI.p, tileJ.p, tileS0.p, tileOff.p, slots.p, tiles, gpuPairs.p, view.pairCount, 0};
  }

  // Pair tiles (PairData): the atoms of the pairs, those moved by torsions first, in blocks of 32; a
  // tile per pair of blocks with pairs. Pair (a, b) of blocks A <= B is taken by the lane of its
  // block-A atom (of a when A = B and (lane b − lane a) mod 32 <= 16) at step (lane of the other −
  // lane) mod 32; the steps of a tile run from its first to its last used step.
  void uploadTiles(const EnergyModel& m) {
    const int nc = m.coordCount;
    std::vector<char> used(nc, 0), moved(nc, 0);
    for (const FlatPair& pr : m.pairs) used[pr.i] = used[pr.j] = 1;
    for (int q : m.moving) moved[q] = 1;
    std::vector<int> order;
    for (int pass = 0; pass < 2; ++pass)
      for (int a = 0; a < nc; ++a)
        if (used[a] && moved[a] == (pass == 0)) order.push_back(a);
    std::vector<int> pos(nc, -1);
    for (std::size_t k = 0; k < order.size(); ++k) pos[order[k]] = static_cast<int>(k);
    const int nb = static_cast<int>((order.size() + 31) / 32);
    std::vector<int> tileOf(static_cast<std::size_t>(nb) * nb, -1);
    std::vector<int> ti, tj, s0, s1, off{0}, sl;
    for (const FlatPair& pr : m.pairs) {
      const int pa = pos[pr.i], pb = pos[pr.j];
      int A = pa / 32, B = pb / 32, la = pa % 32, lb = pb % 32;
      bool swap = false;
      if (A > B || (A == B && ((la - lb) & 31) < ((lb - la) & 31)) || (A == B && ((lb - la) & 31) == 16 && la > lb)) {
        std::swap(A, B);
        std::swap(la, lb);
        swap = true;
      }
      int& id = tileOf[static_cast<std::size_t>(A) * nb + B];
      if (id < 0) {
        id = static_cast<int>(s0.size());
        for (int l = 0; l < 32; ++l) {
          const std::size_t ka = 32 * static_cast<std::size_t>(A) + l, kb = 32 * static_cast<std::size_t>(B) + l;
          ti.push_back(ka < order.size() ? order[ka] : -1);
          tj.push_back(kb < order.size() ? order[kb] : -1);
        }
        s0.push_back(32);
        s1.push_back(0);
        sl.resize(sl.size() + 1024, -1);
      }
      const int s = (lb - la) & 31;
      s0[id] = std::min(s0[id], s);
      s1[id] = std::max(s1[id], s + 1);
      sl[1024 * static_cast<std::size_t>(id) + 32 * s + la] =
          pairFlags(m, pr) | (swap ? kSwap : 0);
    }
    tiles = static_cast<int>(s0.size());
    for (int k = 0; k < tiles; ++k) off.push_back(off.back() + s1[k] - s0[k]);
    if (tiles == 0) {
      ti.assign(32, -1);
      tj.assign(32, -1);
      s0.push_back(0);
      sl.assign(1024, -1);
    }
    tileI.upload(ti);
    tileJ.upload(tj);
    tileS0.upload(s0);
    tileOff.upload(off);
    slots.upload(sl);
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
    for (const FlatPairParams& q : m.pairParams) {
      if (q.m != 12 || (q.n != 6 && q.n != 10)) throw Error("GPU energy: only 12-6 and 12-10 pair potentials");
      pp.push_back({static_cast<T>(q.AC), static_cast<T>(q.BD), static_cast<T>(q.reqm), T(0)});
    }
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
      parallelRange(static_cast<long long>(nn * nt), [&](long long s) { ha[s] = static_cast<T>(hAngles[s]); });
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
      parallelRange(static_cast<long long>(nn), [&](long long s) { std::copy_n(g + s * nt, nt, hGrad + s * nt); });
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
    std::atomic<int> otherA{0};
    parallelRange(static_cast<long long>(nn), [&](long long s) {
      std::copy_n(hAngles + s * nt, nt, ha + s * nt);
      hp[s] = hPoses[s];
      if (hp[s].field != f.current) otherA.store(1, std::memory_order_relaxed);
    });
    const int other = otherA.load();
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
    parallelRange(static_cast<long long>(nn), [&](long long s) {
      hTerms[s] = ft[s];
      hDterms[s] = dt[s];
      if (!grad) return;
      std::copy_n(g6 + s * 6, 6, hGrad + s * ng);
      std::copy_n(gt + s * nt, nt, hGrad + s * ng + 6);
    });
  }

  // Whole 9b minimizations (spec 15.1.26) on the stream of lane 0.
  DevBuf<double> minWork;
  DevBuf<float> minFwork;
  DevBuf<unsigned long long> minNext;
  DevBuf<int> stealSeen;
  PinnedBuf nextHost;  // [0] the initial counter, [1] the final one
  // Stealing (spec 15.1.35): stealHost in mapped host memory; stealMutex serializes the requests.
  DockSteal* stealHost = nullptr;
  int stealSeq = 0;
  std::mutex stealMutex;
  mutable std::mutex runMutex;
  bool running = false;  // runMutex
  std::chrono::steady_clock::time_point runStart;
  cudaEvent_t minDone = nullptr;
  ~Impl() {
    if (minDone) cudaEventDestroy(minDone);
    if (stealHost) cudaFreeHost(stealHost);
  }

  template <int Dock, int NT>
  void minimizeNT(Lane& L, const DockView& dv, DockMinArgs args, int n) {
    const auto kernelS = minimizeDockKernel<Dock, NT, false>;
    const auto kernelT = minimizeDockKernel<Dock, NT, true>;
    const auto [sms, l2] = smsAndL2();
    const int nt = view.torsionCount;
    const std::size_t nc = 3 * static_cast<std::size_t>(view.coordCount);
    // Resident blocks per multiprocessor as launchNT (gradients within the L2 cache, at least 24
    // warps); each block takes samples until none is left.
    constexpr int kMinBlocks = (24 * 32 + NT - 1) / NT;
    const auto perSm = [&](auto kernel, std::size_t bytes) {
      int occ = 0;
      if (bytes > 48 * 1024)
        check(cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(bytes)),
              "cudaFuncSetAttribute");
      check(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&occ, kernel, NT, bytes), "occupancy");
      const std::size_t gb = nc * sizeof(double) * sms;
      const int fit = static_cast<int>(std::min<std::size_t>(l2 / std::max<std::size_t>(gb, 1), 1 << 20));
      return std::clamp(std::max(fit, kMinBlocks), 1, std::max(occ, 1));
    };
    std::size_t shm = sharedBytes<float>(view);
    int per = perSm(kernelS, shm);
    // Torsion transforms only with their work memory in shared memory, as launchNT with several waves.
    bool tor = false;
    if (tOk && shm) {
      const std::size_t withM = shm + (24 * sizeof(float) + 2 * sizeof(int)) * static_cast<std::size_t>(nt);
      if (withM <= static_cast<std::size_t>(sharedLimit()) - 4096) {
        const int perM = perSm(kernelT, withM);
        if (perM >= per && static_cast<std::size_t>(perM) * withM <= static_cast<std::size_t>(sharedPerSm()) / 2) {
          tor = true;
          shm = withM;
          per = perM;
        }
      }
    }
    // The scalars of the minimizer in shared memory when the resident blocks and the shared memory
    // level of the multiprocessor (8 KiB × 2^k) stay the same: more shared memory leaves less L1
    // cache for the grids (spec 15.1.27).
    const bool coords = shm != 0;
    const auto level = [](std::size_t b) {
      std::size_t l = 8192;
      while (l < b) l *= 2;
      return l;
    };
    const std::size_t scalarOff = (shm + 15) & ~std::size_t(15),
                      withS = scalarOff + minScalarDoubles(args.p.memory) * sizeof(double);
    args.scalarOff = -1;
    if (withS <= static_cast<std::size_t>(sharedLimit()) - 4096 &&
        level(static_cast<std::size_t>(per) * withS) <= level(static_cast<std::size_t>(per) * shm) &&
        perSm(tor ? kernelT : kernelS, withS) >= per) {
      args.scalarOff = static_cast<int>(scalarOff);
      shm = withS;
    }
    const int blocks = std::min(n, sms * per);
    if (!coords) L.xs.reserve(static_cast<std::size_t>(blocks) * nc * sizeof(float));
    L.gs.reserve(static_cast<std::size_t>(blocks) * nc);
    minWork.reserve(static_cast<std::size_t>(blocks) * minWorkDoubles(6 + nt, nt, args.p.memory));
    minFwork.reserve(std::max<std::size_t>(1, static_cast<std::size_t>(blocks) * minWorkWords(nt)));
    minNext.reserve(1);
    args.work = minWork.p;
    args.fwork = minFwork.p;
    args.next = minNext.p;
    reinterpret_cast<unsigned long long*>(nextHost.p)[0] = static_cast<unsigned long long>(n) << 32;
    check(cudaMemcpyAsync(minNext.p, nextHost.p, sizeof(unsigned long long), cudaMemcpyHostToDevice, L.stream), "upload");
    TorsionData<float> td{axisF.p, tParent.p, tDeep.p, nullptr, tRounds, tor ? 1 : 0, tor ? 1 : 0};
    setGradientPlan(td);
    PairData<float> pdl = pairData<float>();
    pdl.tiled = n > sms ? 1 : 0;
    (tor ? kernelT : kernelS)<<<blocks, NT, shm, L.stream>>>(view, pdl, td, reinterpret_cast<float*>(L.xs.p), L.gs.p, dv,
                                                             args, coords ? 1 : 0);
    check(cudaGetLastError(), "minimizeDockKernel launch");
  }

  int minimizeDock(const GpuField::Impl& f, const DockMinParams& p, const DockMinStart* hStarts,
                   const double* hAngles, int n, DockMinResult* hOut, float* hAnglesOut, bool stealable) {
    if (n <= 0) return 0;
    nextHost.reserve(3 * sizeof(unsigned long long));
    if (f.current < 0) throw Error("minimizeDock: no field");
    if (p.memory < 1 || p.memory > kDeviceLbfgsMaxMemory) throw Error("minimizeDock: L-BFGS memory outside 1 to " + std::to_string(kDeviceLbfgsMaxMemory));
    check(cudaSetDevice(device), "cudaSetDevice");
    Lane& L = lanes[0];
    if (L.pending) throw Error("minimizeDock: lane 0 busy");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nn = static_cast<std::size_t>(n);
    const std::size_t iAngles = align8(nn * sizeof(DockMinStart)), inBytes = iAngles + nn * nt * sizeof(double);
    const std::size_t oAngles = align8(nn * sizeof(DockMinResult)), outBytes = oAngles + nn * nt * sizeof(float);
    L.in.reserve(inBytes);
    L.out.reserve(outBytes);
    L.hIn.reserve(inBytes);
    L.hOut.reserve(outBytes);
    DockMinStart* hs = reinterpret_cast<DockMinStart*>(L.hIn.p);
    double* ha = reinterpret_cast<double*>(L.hIn.p + iAngles);
    std::atomic<int> otherA{0};
    parallelRange(static_cast<long long>(nn), [&](long long s) {
      hs[s] = hStarts[s];
      std::copy_n(hAngles + s * nt, nt, ha + s * nt);
      if (hs[s].field != f.current) otherA.store(1, std::memory_order_relaxed);
    });
    const bool other = otherA.load();
    check(cudaMemcpyAsync(L.in.p, L.hIn.p, inBytes, cudaMemcpyHostToDevice, L.stream), "upload");
    DockView dv = dock;
    dv.field = f.slots[f.current]->view;
    dv.fields = other ? f.views.p : nullptr;
    DockMinArgs args{p, reinterpret_cast<const DockMinStart*>(L.in.p), reinterpret_cast<const double*>(L.in.p + iAngles),
                     reinterpret_cast<DockMinResult*>(L.out.p), reinterpret_cast<float*>(L.out.p + oAngles),
                     nullptr, nullptr, nullptr, nullptr, nullptr, n, -1};
    if (stealable) {
      if (!stealHost) {
        check(cudaHostAlloc(reinterpret_cast<void**>(&stealHost), sizeof(DockSteal), cudaHostAllocMapped | cudaHostAllocPortable),
              "cudaHostAlloc");
        stealHost->req = stealHost->reply = 0;
        stealSeen.reserve(1);
      }
      void* dp = nullptr;
      check(cudaHostGetDevicePointer(&dp, stealHost, 0), "cudaHostGetDevicePointer");
      args.steal = static_cast<DockSteal*>(dp);
      args.stealSeen = stealSeen.p;
      std::lock_guard lk(runMutex);
      stealHost->taken = 0;
      stealHost->oldEnd = stealHost->newEnd = n;
      // Requests before this call are not answered.
      int* const seen = reinterpret_cast<int*>(reinterpret_cast<unsigned long long*>(nextHost.p) + 2);
      *seen = stealSeq;
      check(cudaMemcpyAsync(stealSeen.p, seen, sizeof(int), cudaMemcpyHostToDevice, L.stream), "upload");
      running = true;
      runStart = std::chrono::steady_clock::now();
    }
    const auto run = [&]<int D>() {
      switch (blockThreads(view, dv.probeCount)) {
        case 32: minimizeNT<D, 32>(L, dv, args, n); break;
        case 64: minimizeNT<D, 64>(L, dv, args, n); break;
        case 128: minimizeNT<D, 128>(L, dv, args, n); break;
        default: minimizeNT<D, 256>(L, dv, args, n); break;
      }
    };
    if (other) run.template operator()<2>();
    else run.template operator()<1>();
    check(cudaMemcpyAsync(L.hOut.p, L.out.p, outBytes, cudaMemcpyDeviceToHost, L.stream), "download");
    unsigned long long* const nh = reinterpret_cast<unsigned long long*>(nextHost.p);
    check(cudaMemcpyAsync(nh + 1, minNext.p, sizeof(unsigned long long), cudaMemcpyDeviceToHost, L.stream), "download");
    // The calling thread sleeps until the kernel ends (no spinning on the CPU).
    if (!minDone) check(cudaEventCreateWithFlags(&minDone, cudaEventBlockingSync | cudaEventDisableTiming), "cudaEventCreate");
    check(cudaEventRecord(minDone, L.stream), "cudaEventRecord");
    const cudaError_t e = cudaEventSynchronize(minDone);
    if (stealable) {
      std::lock_guard lk(runMutex);
      running = false;
    }
    check(e, "minimizeDockKernel");
    const std::size_t end = std::min(nn, static_cast<std::size_t>(nh[1] >> 32));
    std::copy_n(reinterpret_cast<const DockMinResult*>(L.hOut.p), end, hOut);
    std::copy_n(reinterpret_cast<const float*>(L.hOut.p + oAngles), end * nt, hAnglesOut);
    return static_cast<int>(end);
  }

  std::pair<int, int> stealDock(double frac, int giveMin) {
    std::lock_guard sl(stealMutex);
    int seq;
    {
      std::lock_guard lk(runMutex);
      if (!running) return {0, 0};
      stealHost->frac = static_cast<float>(frac);
      stealHost->giveMin = std::max(1, giveMin);
      seq = ++stealSeq;
      std::atomic_thread_fence(std::memory_order_seq_cst);
      stealHost->req = seq;
    }
    for (;;) {
      bool run;
      {
        std::lock_guard lk(runMutex);
        run = running;
      }
      std::atomic_thread_fence(std::memory_order_seq_cst);
      if (stealHost->reply == seq) {
        std::atomic_thread_fence(std::memory_order_seq_cst);
        const int ne = stealHost->newEnd, oe = stealHost->oldEnd;
        return ne < oe ? std::pair{ne, oe} : std::pair{0, 0};
      }
      // The kernel ended before taking another sample: nothing given away.
      if (!run) return {0, 0};
      std::this_thread::sleep_for(std::chrono::microseconds(20));
    }
  }

  bool dockProgress(int& taken, int& end, double& seconds) const {
    std::lock_guard lk(runMutex);
    if (!running) return false;
    std::atomic_thread_fence(std::memory_order_seq_cst);
    taken = stealHost->taken;
    end = stealHost->newEnd;
    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - runStart).count();
    return true;
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
  for (std::size_t p = 0; p < gp.size(); ++p) gp[p] = {m.pairs[p].i, m.pairs[p].j, pairFlags(m, m.pairs[p])};
  d.gpuPairs.upload(gp);
  d.uploadTiles(m);
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

int GpuEnergy::minimizeDock(const GpuField& field, const DockMinParams& p, const DockMinStart* starts,
                            const double* anglesDeg, int n, DockMinResult* out, float* anglesOut, bool stealable) {
  return impl_->minimizeDock(*field.impl_, p, starts, anglesDeg, n, out, anglesOut, stealable);
}

std::pair<int, int> GpuEnergy::stealDock(double frac, int giveMin) { return impl_->stealDock(frac, giveMin); }

bool GpuEnergy::dockProgress(int& taken, int& end, double& seconds) const { return impl_->dockProgress(taken, end, seconds); }

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
