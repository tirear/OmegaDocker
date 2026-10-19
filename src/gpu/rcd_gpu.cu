// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "rcd_gpu.hpp"

#include <cooperative_groups.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <numbers>
#include <string>

#include "error.hpp"
#include "fiber_batch.hpp"

namespace od {
namespace {

namespace cg = cooperative_groups;

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
  T* reserve(std::size_t count) {
    if (count > n) {
      if (p) cudaFree(p);
      p = nullptr;
      n = 0;
      check(cudaMalloc(&p, std::max<std::size_t>(count, 1) * sizeof(T)), "cudaMalloc");
      n = count;
    }
    return p;
  }
  T* upload(const T* v, std::size_t count) {
    reserve(count);
    if (count) check(cudaMemcpy(p, v, count * sizeof(T), cudaMemcpyHostToDevice), "upload");
    return p;
  }
  void download(T* v, std::size_t count) const {
    if (count) check(cudaMemcpy(v, p, count * sizeof(T), cudaMemcpyDeviceToHost), "download");
  }
};

constexpr double kDeg = 180.0 / std::numbers::pi;
constexpr int kMaxSmem = 48 * 1024;
constexpr int kMaxWarps = 8;
// Updates per attempt and launch: a long attempt continues in the next launches.
constexpr int kChunk = 1024;
// Single-precision points drift slightly from the coordinates built from the angles: converge at
// 99% of rcdTolerance so that the check on the CPU passes.
constexpr double kTolFactor = 0.99;

struct FlatView {
  const int *pointNum, *dofOff, *dofNum, *pairOff, *pairNum, *copyOff, *copyNum;
  const int *movesOff, *moves, *axisFixed, *axisMoving, *pairPoint, *pairEnd;
  const unsigned *actMask, *actSecond;  // per dof: bit q: pair q active / its second point moves
};

// Per-slot state kept on the device between launches (an attempt runs over several launches).
struct SlotView {
  float* P;      // slots × 3 maxNp
  double* A;     // slots × maxNd
  int* tried;    // slots × maxNd
  int *link, *prev, *triedNum, *it, *maxIter, *conv, *done;
  unsigned long long* rng;
};
struct JobView {
  int n;
  const int* slot;
  const int* isNew;  // per job: start from the uploaded points, angles and random state
  const int *link, *pOff, *aOff, *maxIter;
  const double *P0, *A0;
  const unsigned long long* rng0;
};

struct V3 {
  float x, y, z;
};
__device__ inline V3 sub(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
__device__ inline V3 add(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
__device__ inline V3 mul(float s, V3 a) { return {s * a.x, s * a.y, s * a.z}; }
__device__ inline float dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__device__ inline V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
__device__ inline V3 unit(V3 a) {
  const float l = sqrtf(dot(a, a));
  return l > 0 ? mul(1.0f / l, a) : V3{0, 0, 0};
}
__device__ inline V3 perp(V3 v, V3 k) { return sub(v, mul(dot(v, k), k)); }
__device__ inline V3 rotateAbout(V3 p, V3 o, V3 k, float c, float s) {
  const V3 v = sub(p, o);
  return add(add(add(o, mul(c, v)), mul(s, cross(k, v))), mul(dot(k, v) * (1 - c), k));
}
__device__ inline V3 load(const float* P, int i) { return {P[3 * i], P[3 * i + 1], P[3 * i + 2]}; }

// Rng of rng.hpp.
__device__ inline double uniform(unsigned long long& s) {
  s += 0x9e3779b97f4a7c15ULL;
  unsigned long long z = s;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z = z ^ (z >> 31);
  return static_cast<double>(z >> 11) * 0x1.0p-53;
}
__device__ inline int below(unsigned long long& s, int n) {
  return static_cast<int>(static_cast<unsigned long long>(uniform(s) * static_cast<double>(n)));
}

template <class Tile>
__device__ inline float tileSum(const Tile& t, float v) {
  constexpr int G = Tile::num_threads();
#pragma unroll
  for (int o = G / 2; o > 0; o /= 2) v += t.shfl_xor(v, o);
  return v;
}

// One tile of G lanes per attempt, up to `chunk` updates per launch: lane q holds closure pair q;
// every lane draws the same random numbers; the points moved by an update are rotated by the
// lanes in turn. As LoopCloser::rcd without the clash term and metal distances, with the points
// in single precision (the angles are summed in double; the CPU checks the closure from them).
template <int G>
__global__ void rcdKernel(FlatView g, SlotView sv, JobView jv, float tol, float stallTol, int maxNp, int maxNd,
                          int strideFloats, int chunk) {
  extern __shared__ double smemD[];
  const auto tile = cg::tiled_partition<G>(cg::this_thread_block());
  const int lane = static_cast<int>(tile.thread_rank());
  const int grp = threadIdx.x / G, ng = blockDim.x / G;
  const int j = blockIdx.x * ng + grp;
  if (j >= jv.n) return;
  const int slot = jv.slot[j];
  double* A = smemD + static_cast<std::size_t>(grp) * strideFloats / 2;
  float* P = reinterpret_cast<float*>(A + maxNd);
  int* tried = reinterpret_cast<int*>(P + 3 * maxNp);
  float* Pg = sv.P + static_cast<std::size_t>(slot) * 3 * maxNp;
  double* Ag = sv.A + static_cast<std::size_t>(slot) * maxNd;
  int* Tg = sv.tried + static_cast<std::size_t>(slot) * maxNd;
  const bool fresh = jv.isNew[j] != 0;
  const int L = fresh ? jv.link[j] : sv.link[slot];
  const int np = __ldg(g.pointNum + L), nd = __ldg(g.dofNum + L), d0 = __ldg(g.dofOff + L);
  const int npair = __ldg(g.pairNum + L), nc = __ldg(g.copyNum + L);
  const int* pp = g.pairPoint + 2 * __ldg(g.pairOff + L);
  const int* pe = g.pairEnd + __ldg(g.copyOff + L);
  if (fresh) {
    const double* Pin = jv.P0 + jv.pOff[j];
    for (int i = lane; i < 3 * np; i += G) P[i] = static_cast<float>(Pin[i]);
    for (int i = lane; i < nd; i += G) {
      A[i] = jv.A0[jv.aOff[j] + i];
      tried[i] = 0;
    }
  } else {
    for (int i = lane; i < 3 * np; i += G) P[i] = Pg[i];
    for (int i = lane; i < nd; i += G) {
      A[i] = Ag[i];
      tried[i] = Tg[i];
    }
  }
  // Lane q: pair q (points i1, i2), its copy link and the size of the copy link.
  const bool hasPair = lane < npair;
  int i1 = 0, i2 = 0, myCopy = -1;
  if (hasPair) {
    i1 = __ldg(pp + 2 * lane);
    i2 = __ldg(pp + 2 * lane + 1);
    for (int c = 0, b = 0; c < nc; b = __ldg(pe + c++))
      if (lane >= b && lane < __ldg(pe + c)) myCopy = c;
  }
  tile.sync();
  float dd = 0;
  const auto pairDist = [&]() {
    if (hasPair) {
      const V3 v = sub(load(P, i1), load(P, i2));
      dd = dot(v, v);
    }
  };
  const auto dev = [&]() {
    float dm = 0;
    for (int c = 0, b = 0; c < nc; b = __ldg(pe + c++)) {
      const float s = tileSum(tile, myCopy == c ? dd : 0.0f);
      const float x = sqrtf(s / static_cast<float>(__ldg(pe + c) - b));
      dm = dm < x ? x : dm;
    }
    return dm;
  };
  pairDist();
  unsigned long long s;
  int prev, triedNum, it, maxIter;
  bool conv;
  if (fresh) {
    s = jv.rng0[j];
    prev = -1;
    triedNum = 0;
    it = 0;
    maxIter = jv.maxIter[j];
    conv = !(dev() > tol);
  } else {
    s = sv.rng[slot];
    prev = sv.prev[slot];
    triedNum = sv.triedNum[slot];
    it = sv.it[slot];
    maxIter = sv.maxIter[slot];
    conv = sv.conv[slot] != 0;
  }
  bool stalled = false;
  const int end = min(maxIter, it + chunk);
  for (; it < end && !conv; ++it) {
    int d = below(s, nd);
    if (nd > 1 && d == prev) d = (d + 1 + below(s, nd - 1)) % nd;
    prev = d;
    const int gd = d0 + d;
    const V3 o = load(P, __ldg(g.axisFixed + gd));
    const V3 k = unit(sub(load(P, __ldg(g.axisMoving + gd)), o));
    const unsigned mask = __ldg(g.actMask + gd);
    const bool active = hasPair && ((mask >> lane) & 1u);
    float nu = 0, de = 0;
    if (active) {
      const bool second = (__ldg(g.actSecond + gd) >> lane) & 1u;
      const V3 rp = perp(sub(load(P, second ? i2 : i1), o), k), tp = perp(sub(load(P, second ? i1 : i2), o), k);
      nu = dot(k, cross(rp, tp));
      de = dot(rp, tp);
    }
    const float num = tileSum(tile, nu), den = tileSum(tile, de);
    float delta = 0;
    if (mask != 0 && (num != 0 || den != 0)) delta = atan2f(num, den);
    if (delta != 0) {
      float sn, c;
      sincosf(delta, &sn, &c);
      tile.sync();
      const int mb = __ldg(g.movesOff + gd), me = __ldg(g.movesOff + gd + 1);
      for (int x = mb + lane; x < me; x += G) {
        const int i = __ldg(g.moves + x);
        const V3 r = rotateAbout(load(P, i), o, k, c, sn);
        P[3 * i] = r.x;
        P[3 * i + 1] = r.y;
        P[3 * i + 2] = r.z;
      }
      if (lane == 0) A[d] += static_cast<double>(delta) * kDeg;
      tile.sync();
      if (active) pairDist();
      conv = !(dev() > tol);
    }
    if (stallTol > 0) {
      if (fabsf(delta) > stallTol) {
        tile.sync();
        for (int i = lane; i < nd; i += G) tried[i] = 0;
        tile.sync();
        triedNum = 0;
      }
      const int was = tried[d];
      tile.sync();
      if (!was) {
        if (lane == 0) tried[d] = 1;
        ++triedNum;
      }
      tile.sync();
      if (triedNum == nd) {
        ++it;
        stalled = true;
        break;
      }
    }
  }
  const bool done = conv || stalled || it >= maxIter;
  tile.sync();
  for (int i = lane; i < nd; i += G) Ag[i] = A[i];
  if (!done) {
    for (int i = lane; i < 3 * np; i += G) Pg[i] = P[i];
    for (int i = lane; i < nd; i += G) Tg[i] = tried[i];
  }
  if (lane == 0) {
    sv.link[slot] = L;
    sv.rng[slot] = s;
    sv.prev[slot] = prev;
    sv.triedNum[slot] = triedNum;
    sv.it[slot] = it;
    sv.maxIter[slot] = maxIter;
    sv.conv[slot] = conv ? 1 : 0;
    sv.done[slot] = done ? 1 : 0;
  }
}

}  // namespace

struct GpuRcd::Impl {
  RcdFlat flat;
  int device = 0, nslots = 1;
  float tol = 0, stallTol = 0;
  int maxNp = 0, maxNd = 0, maxPairs = 0, strideFloats = 0, tileSize = 32, groups = 1;
  bool ok = false;
  DevBuf<int> pointNum, dofOff, dofNum, pairOff, pairNum, copyOff, copyNum, movesOff, moves, axisFixed, axisMoving,
      pairPoint, pairEnd;
  DevBuf<unsigned> actMask, actSecond;
  // Per slot.
  DevBuf<float> sP;
  DevBuf<double> sA;
  DevBuf<int> sTried, sLink, sPrev, sTriedNum, sIt, sMaxIter, sConv, sDone;
  DevBuf<unsigned long long> sRng;
  // Per launch.
  DevBuf<int> jSlot, jNew, jLink, jPOff, jAOff, jMaxIter;
  DevBuf<double> jP, jA;
  DevBuf<unsigned long long> jRng;
  std::vector<RcdJob*> pending;
  std::vector<char> started;  // per slot: the attempt of the slot has run in a launch
  std::vector<char> finished;
  long long kernels = 0, attempts = 0, updates = 0;
  double sec = 0;
};

GpuRcd::GpuRcd(const LoopCloser& closer, int device, int slots) : impl_(std::make_unique<Impl>()) {
  Impl& m = *impl_;
  m.flat = closer.flat();
  m.device = device;
  m.nslots = std::max(1, slots);
  m.tol = static_cast<float>(kTolFactor * closer.params().rcdTolerance);
  m.stallTol = static_cast<float>(closer.params().rcdStallTolerance);
  const RcdFlat& f = m.flat;
  for (std::size_t l = 0; l < f.eligible.size(); ++l)
    if (f.eligible[l]) {
      m.maxNp = std::max(m.maxNp, f.pointNum[l]);
      m.maxNd = std::max(m.maxNd, f.dofNum[l]);
      m.maxPairs = std::max(m.maxPairs, f.pairNum[l]);
    }
  // One lane per closure pair.
  if (m.maxNd == 0 || m.maxPairs > 32) return;
  m.tileSize = m.maxPairs <= 8 ? 8 : m.maxPairs <= 16 ? 16 : 32;
  m.strideFloats = 2 * m.maxNd + (3 * m.maxNp + m.maxNd + 1) / 2 * 2;  // angles (double), points, tried
  const int bytes = m.strideFloats * 4;
  const int perWarp = 32 / m.tileSize;
  if (bytes * perWarp > kMaxSmem) return;
  m.groups = perWarp * std::clamp(kMaxSmem / (bytes * perWarp), 1, kMaxWarps);
  // Per dof: active pairs as bit masks.
  std::vector<unsigned> mask(f.axisFixed.size(), 0), second(f.axisFixed.size(), 0);
  for (std::size_t g = 0; g < mask.size(); ++g)
    for (int a = f.activeOff[g]; a < f.activeOff[g + 1]; ++a) {
      const int q = f.active[a] / 2;
      mask[g] |= 1u << q;
      if (f.active[a] & 1) second[g] |= 1u << q;
    }
  check(cudaSetDevice(device), "cudaSetDevice");
  m.pointNum.upload(f.pointNum.data(), f.pointNum.size());
  m.dofOff.upload(f.dofOff.data(), f.dofOff.size());
  m.dofNum.upload(f.dofNum.data(), f.dofNum.size());
  m.pairOff.upload(f.pairOff.data(), f.pairOff.size());
  m.pairNum.upload(f.pairNum.data(), f.pairNum.size());
  m.copyOff.upload(f.copyOff.data(), f.copyOff.size());
  m.copyNum.upload(f.copyNum.data(), f.copyNum.size());
  m.movesOff.upload(f.movesOff.data(), f.movesOff.size());
  m.moves.upload(f.moves.data(), f.moves.size());
  m.axisFixed.upload(f.axisFixed.data(), f.axisFixed.size());
  m.axisMoving.upload(f.axisMoving.data(), f.axisMoving.size());
  m.pairPoint.upload(f.pairPoint.data(), f.pairPoint.size());
  m.pairEnd.upload(f.pairEnd.data(), f.pairEnd.size());
  m.actMask.upload(mask.data(), mask.size());
  m.actSecond.upload(second.data(), second.size());
  const std::size_t ns = static_cast<std::size_t>(m.nslots);
  m.sP.reserve(ns * 3 * m.maxNp);
  m.sA.reserve(ns * m.maxNd);
  m.sTried.reserve(ns * m.maxNd);
  for (auto* b : {&m.sLink, &m.sPrev, &m.sTriedNum, &m.sIt, &m.sMaxIter, &m.sConv, &m.sDone}) b->reserve(ns);
  m.sRng.reserve(ns);
  m.pending.assign(ns, nullptr);
  m.started.assign(ns, 0);
  m.finished.assign(ns, 0);
  m.ok = true;
}

GpuRcd::~GpuRcd() = default;

bool GpuRcd::active() const { return impl_->ok; }
int GpuRcd::slots() const { return impl_->nslots; }
long long GpuRcd::kernels() const { return impl_->kernels; }
long long GpuRcd::attempts() const { return impl_->attempts; }
long long GpuRcd::updates() const { return impl_->updates; }
double GpuRcd::kernelSec() const { return impl_->sec; }

void GpuRcd::run(int slot, RcdJob& job) {
  Impl& m = *impl_;
  m.started[slot] = 0;
  m.finished[slot] = 0;
  // An attempt runs over several launches of at most kChunk updates each.
  while (!m.finished[slot]) {
    m.pending[slot] = &job;
    batchWait();
  }
}

void GpuRcd::flush(std::span<const int> slots) {
  Impl& m = *impl_;
  std::vector<int> js;
  for (int s : slots)
    if (m.pending[s]) js.push_back(s);
  const int n = static_cast<int>(js.size());
  if (n == 0) return;
  // Jobs of the same link next to each other (the same control flow within a warp).
  std::stable_sort(js.begin(), js.end(), [&](int a, int b) { return m.pending[a]->link < m.pending[b]->link; });
  std::vector<int> isNew(n), link(n), pOff(n), aOff(n), maxIter(n);
  std::vector<unsigned long long> rng(n);
  std::vector<double> P, A;
  for (int i = 0; i < n; ++i) {
    const int s = js[i];
    const RcdJob& j = *m.pending[s];
    isNew[i] = m.started[s] ? 0 : 1;
    link[i] = j.link;
    pOff[i] = static_cast<int>(P.size());
    aOff[i] = static_cast<int>(A.size());
    maxIter[i] = j.maxIter;
    rng[i] = j.rng;
    if (!m.started[s]) {
      P.insert(P.end(), j.points.begin(), j.points.end());
      A.insert(A.end(), j.anglesDeg.begin(), j.anglesDeg.end());
    }
  }
  check(cudaSetDevice(m.device), "cudaSetDevice");
  const auto t0 = std::chrono::steady_clock::now();
  const JobView jv{n,
                   m.jSlot.upload(js.data(), n),
                   m.jNew.upload(isNew.data(), n),
                   m.jLink.upload(link.data(), n),
                   m.jPOff.upload(pOff.data(), n),
                   m.jAOff.upload(aOff.data(), n),
                   m.jMaxIter.upload(maxIter.data(), n),
                   m.jP.upload(P.data(), P.size()),
                   m.jA.upload(A.data(), A.size()),
                   m.jRng.upload(rng.data(), n)};
  const SlotView sv{m.sP.p, m.sA.p, m.sTried.p, m.sLink.p, m.sPrev.p, m.sTriedNum.p,
                    m.sIt.p, m.sMaxIter.p, m.sConv.p, m.sDone.p, m.sRng.p};
  const FlatView fv{m.pointNum.p, m.dofOff.p, m.dofNum.p, m.pairOff.p, m.pairNum.p, m.copyOff.p, m.copyNum.p,
                    m.movesOff.p, m.moves.p, m.axisFixed.p, m.axisMoving.p, m.pairPoint.p, m.pairEnd.p,
                    m.actMask.p, m.actSecond.p};
  const int blocks = (n + m.groups - 1) / m.groups;
  const int threads = m.groups * m.tileSize;
  const std::size_t smem = static_cast<std::size_t>(m.groups) * m.strideFloats * 4;
  if (m.tileSize == 8)
    rcdKernel<8><<<blocks, threads, smem>>>(fv, sv, jv, m.tol, m.stallTol, m.maxNp, m.maxNd, m.strideFloats, kChunk);
  else if (m.tileSize == 16)
    rcdKernel<16><<<blocks, threads, smem>>>(fv, sv, jv, m.tol, m.stallTol, m.maxNp, m.maxNd, m.strideFloats, kChunk);
  else
    rcdKernel<32><<<blocks, threads, smem>>>(fv, sv, jv, m.tol, m.stallTol, m.maxNp, m.maxNd, m.strideFloats, kChunk);
  check(cudaGetLastError(), "rcdKernel");
  check(cudaDeviceSynchronize(), "rcdKernel");
  // Results of the finished attempts.
  std::vector<int> done(m.nslots), conv(m.nslots), it(m.nslots);
  std::vector<unsigned long long> rs(m.nslots);
  std::vector<double> As(static_cast<std::size_t>(m.nslots) * m.maxNd);
  m.sDone.download(done.data(), m.nslots);
  m.sConv.download(conv.data(), m.nslots);
  m.sIt.download(it.data(), m.nslots);
  m.sRng.download(rs.data(), m.nslots);
  m.sA.download(As.data(), As.size());
  m.sec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  ++m.kernels;
  for (int s : js) {
    RcdJob& j = *m.pending[s];
    m.pending[s] = nullptr;
    m.started[s] = 1;
    if (!done[s]) continue;
    m.finished[s] = 1;
    std::copy_n(As.begin() + static_cast<std::ptrdiff_t>(s) * m.maxNd, j.anglesDeg.size(), j.anglesDeg.begin());
    j.rng = rs[s];
    j.converged = conv[s] != 0;
    j.iterations = it[s];
    ++m.attempts;
    m.updates += it[s];
  }
}

}  // namespace od
