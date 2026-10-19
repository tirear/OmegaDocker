// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "shake_gpu.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <string>

#include "ad4_pair.hpp"
#include "error.hpp"

namespace od {
namespace {

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
  T* upload(const T* v, std::size_t count) {
    if (count > n) {
      if (p) cudaFree(p);
      p = nullptr;
      n = 0;
      check(cudaMalloc(&p, std::max<std::size_t>(count, 1) * sizeof(T)), "cudaMalloc");
      n = count;
    }
    if (count) check(cudaMemcpy(p, v, count * sizeof(T), cudaMemcpyHostToDevice), "upload");
    return p;
  }
};

constexpr int kThreads = 256;
constexpr int kMaxPointsPerThread = 16;  // sphere points per thread (radius / step <= 9)
constexpr int kTypes = static_cast<int>(Ad4Type::Count);

struct TypePair {
  float AC, BD, reqm;
  int m, n;
};

struct ShakeView {
  int na, np;
  const float* xyz;  // 3 × na
  const int* type;
  const float* q;   // charge
  const float* S;   // solpar + qasp |q|
  const float* V;   // volume
  const int* nearStart;  // na + 1
  const int* near;       // sorted per atom
  const TypePair* pairs;  // kTypes × kTypes
  const float* pts;       // 3 × np
  float elec, desolv, half, rMin, cut2, vdw;
  double reach2, rrt;
};

struct Partner {
  float x, y, z, AC, BD, reqm, q, sv;
  int mn;  // m * 256 + n
};

__global__ void __launch_bounds__(kThreads) shakeKernel(ShakeView v, const int* atoms, double* lnF) {
  __shared__ Partner tile[kThreads];
  __shared__ int count;
  __shared__ int warpSum[kThreads / 32];
  __shared__ double red[kThreads / 32];
  const int i = atoms[blockIdx.x];
  const int t = threadIdx.x, lane = t & 31, warp = t >> 5;
  const float xi = v.xyz[3 * i], yi = v.xyz[3 * i + 1], zi = v.xyz[3 * i + 2];
  const int ti = v.type[i];
  const float qi = v.q[i], Si = v.S[i], Vi = v.V[i];
  const int n0 = v.nearStart[i], n1 = v.nearStart[i + 1];
  float px[kMaxPointsPerThread], py[kMaxPointsPerThread], pz[kMaxPointsPerThread];
  double e[kMaxPointsPerThread];
  int nk = 0;
  for (int k = t; k < v.np; k += kThreads, ++nk) {
    px[nk] = xi + v.pts[3 * k];
    py[nk] = yi + v.pts[3 * k + 1];
    pz[nk] = zi + v.pts[3 * k + 2];
    e[nk] = 0;
  }
  for (int j0 = 0; j0 < v.na; j0 += kThreads) {
    // Partners of this tile in atom order (block-wide prefix sum of the flags).
    const int j = j0 + t;
    bool ok = j < v.na && j != i;
    float dx = 0, dy = 0, dz = 0;
    if (ok) {
      dx = v.xyz[3 * j] - xi;
      dy = v.xyz[3 * j + 1] - yi;
      dz = v.xyz[3 * j + 2] - zi;
      ok = static_cast<double>(dx) * dx + static_cast<double>(dy) * dy + static_cast<double>(dz) * dz < v.reach2;
    }
    if (ok) {
      int lo = n0, hi = n1;
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if (v.near[mid] < j) lo = mid + 1;
        else hi = mid;
      }
      ok = !(lo < n1 && v.near[lo] == j);
    }
    const unsigned ballot = __ballot_sync(0xffffffffu, ok);
    if (lane == 0) warpSum[warp] = __popc(ballot);
    __syncthreads();
    if (t == 0) {
      int s = 0;
      for (int w = 0; w < kThreads / 32; ++w) {
        const int c = warpSum[w];
        warpSum[w] = s;
        s += c;
      }
      count = s;
    }
    __syncthreads();
    if (ok) {
      const int slot = warpSum[warp] + __popc(ballot & ((1u << lane) - 1u));
      const TypePair& tp = v.pairs[ti * kTypes + v.type[j]];
      Partner& p = tile[slot];
      p.x = v.xyz[3 * j];
      p.y = v.xyz[3 * j + 1];
      p.z = v.xyz[3 * j + 2];
      p.AC = tp.AC;
      p.BD = tp.BD;
      p.reqm = tp.reqm;
      p.q = v.elec * qi * v.q[j];
      p.sv = v.desolv * (Si * v.V[j] + v.S[j] * Vi);
      p.mn = tp.m * 256 + tp.n;
    }
    __syncthreads();
    const int nb = count;
    for (int k = 0; k < nk; ++k) {
      float sum = 0;
      for (int b = 0; b < nb; ++b) {
        const Partner& p = tile[b];
        const float ddx = px[k] - p.x, ddy = py[k] - p.y, ddz = pz[k] - p.z;
        const float r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 >= v.cut2) continue;
        const float r = sqrtf(r2);
        float tt = 0, ee, de;
        if (r < v.vdw) {
          ad4VdwHb<float>(r, p.AC, p.BD, p.mn >> 8, p.mn & 255, p.reqm, v.half, ee, de);
          tt += ee;
        }
        ad4Elec<float>(r, p.q, v.rMin, ee, de);
        tt += ee;
        ad4Desolv<float>(r, p.sv, v.rMin, ee, de);
        sum += tt + ee;
      }
      e[k] += sum;
    }
    __syncthreads();
  }
  // Minimum over the points, then sum_k exp(-(E_k - E_min) * rrt).
  double m = INFINITY;
  for (int k = 0; k < nk; ++k) m = fmin(m, e[k]);
  for (int o = 16; o > 0; o >>= 1) m = fmin(m, __shfl_xor_sync(0xffffffffu, m, o));
  if (lane == 0) red[warp] = m;
  __syncthreads();
  if (t == 0) {
    for (int w = 1; w < kThreads / 32; ++w) red[0] = fmin(red[0], red[w]);
  }
  __syncthreads();
  const double emin = red[0];
  __syncthreads();
  double w = 0;
  for (int k = 0; k < nk; ++k) w += exp(-(e[k] - emin) * v.rrt);
  for (int o = 16; o > 0; o >>= 1) w += __shfl_xor_sync(0xffffffffu, w, o);
  if (lane == 0) red[warp] = w;
  __syncthreads();
  if (t == 0) {
    double s = 0;
    for (int x = 0; x < kThreads / 32; ++x) s += red[x];
    lnF[blockIdx.x] = log(s / static_cast<double>(v.np));
  }
}

}  // namespace

struct GpuShaker::Impl {
  int device = 0;
  ShakeView base{};
  std::vector<float> pts;
  DevBuf<float> xyz, q, S, V, dPts;
  DevBuf<int> type, nearStart, near, atoms;
  DevBuf<TypePair> pairs;
  DevBuf<double> out;
  double qasp = 0;
};

GpuShaker::GpuShaker(const EnergyParams& ep, const ShakeParams& p, int device) : impl_(std::make_unique<Impl>()) {
  Impl& d = *impl_;
  d.device = device;
  check(cudaSetDevice(device), "cudaSetDevice");
  const std::vector<Vec3f> pts = shakePoints(p.radius, p.step);
  if (pts.size() > static_cast<std::size_t>(kThreads * kMaxPointsPerThread))
    throw Error("GPU shaking: too many sphere points (od_0100.shaking.radius / step)");
  for (const Vec3f& x : pts) d.pts.insert(d.pts.end(), {x[0], x[1], x[2]});
  d.dPts.upload(d.pts.data(), d.pts.size());
  std::vector<TypePair> tp(kTypes * kTypes);
  for (int a = 0; a < kTypes; ++a)
    for (int b = 0; b < kTypes; ++b) {
      const Ad4PairParams q = ad4PairParams(static_cast<Ad4Type>(a), static_cast<Ad4Type>(b), ep.coeffs);
      tp[a * kTypes + b] = {static_cast<float>(q.AC), static_cast<float>(q.BD), static_cast<float>(q.reqm), q.m, q.n};
    }
  d.pairs.upload(tp.data(), tp.size());
  ShakeView& v = d.base;
  v.np = static_cast<int>(pts.size());
  v.pairs = d.pairs.p;
  v.pts = d.dPts.p;
  v.elec = static_cast<float>(ep.coeffs.elec);
  v.desolv = static_cast<float>(ep.coeffs.desolv);
  v.half = static_cast<float>(0.5 * ep.smooth);
  v.rMin = static_cast<float>(ep.elecMinDistance);
  const double reach = kElecCutoff + p.radius;
  v.reach2 = reach * reach;
  v.cut2 = static_cast<float>(kElecCutoff * kElecCutoff);
  v.vdw = static_cast<float>(kVdwCutoff);
  v.rrt = p.rrtScaling;
  d.qasp = ep.coeffs.qasp;
}

GpuShaker::~GpuShaker() = default;

std::vector<double> GpuShaker::lnF(const ShakeAtoms& atoms, std::span<const char> shake) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const std::size_t na = atoms.size();
  std::vector<double> lnF(na, 0.0);
  std::vector<int> list;
  for (std::size_t i = 0; i < na; ++i)
    if (atoms.heavy[i] && (shake.empty() || shake[i])) list.push_back(static_cast<int>(i));
  if (list.empty()) return lnF;
  std::vector<float> xyz(3 * na), q(na), S(na), V(na);
  std::vector<int> type(na), start(na + 1, 0), near;
  for (std::size_t a = 0; a < na; ++a) {
    for (int c = 0; c < 3; ++c) xyz[3 * a + c] = atoms.xyz[a][c];
    type[a] = static_cast<int>(atoms.type[a]);
    q[a] = atoms.charge[a];
    S[a] = static_cast<float>(ad4Params(atoms.type[a]).solpar + d.qasp * std::abs(atoms.charge[a]));
    V[a] = static_cast<float>(ad4Params(atoms.type[a]).volume);
    start[a] = static_cast<int>(near.size());
    near.insert(near.end(), atoms.near[a].begin(), atoms.near[a].end());
  }
  start[na] = static_cast<int>(near.size());
  ShakeView v = d.base;
  v.na = static_cast<int>(na);
  v.xyz = d.xyz.upload(xyz.data(), xyz.size());
  v.type = d.type.upload(type.data(), type.size());
  v.q = d.q.upload(q.data(), q.size());
  v.S = d.S.upload(S.data(), S.size());
  v.V = d.V.upload(V.data(), V.size());
  v.nearStart = d.nearStart.upload(start.data(), start.size());
  v.near = d.near.upload(near.data(), near.size());
  const int* dAtoms = d.atoms.upload(list.data(), list.size());
  std::vector<double> res(list.size());
  if (d.out.n < res.size()) d.out.upload(res.data(), res.size());
  shakeKernel<<<static_cast<unsigned>(list.size()), kThreads>>>(v, dAtoms, d.out.p);
  check(cudaGetLastError(), "shakeKernel launch");
  check(cudaMemcpy(res.data(), d.out.p, res.size() * sizeof(double), cudaMemcpyDeviceToHost), "shakeKernel");
  for (std::size_t k = 0; k < list.size(); ++k) lnF[list[k]] = res[k];
  return lnF;
}

}  // namespace od
