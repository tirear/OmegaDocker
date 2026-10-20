// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "prep_gpu.hpp"

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

#include "dock_kernel.hpp"
#include "error.hpp"

namespace od {
namespace {

void check(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw Error(std::string("CUDA: ") + what + ": " + cudaGetErrorString(e));
}
void checkFft(cufftResult r, const char* what) {
  if (r != CUFFT_SUCCESS) throw Error(std::string("cuFFT: ") + what + ": error " + std::to_string(static_cast<int>(r)));
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
    release();
    check(cudaMalloc(&p, std::max<std::size_t>(count, 1) * sizeof(T)), "cudaMalloc");
    n = count;
  }
  void release() {
    if (p) cudaFree(p);
    p = nullptr;
    n = 0;
  }
  // On stream s; pageable v is copied to staging memory before the call returns.
  void upload(const T* v, std::size_t count, cudaStream_t s) {
    reserve(count);
    if (count) check(cudaMemcpyAsync(p, v, count * sizeof(T), cudaMemcpyHostToDevice, s), "upload");
  }
  void upload(const std::vector<T>& v, cudaStream_t s) { upload(v.data(), v.size(), s); }
};

constexpr int kFilterThreads = 128;
constexpr long long kPosm2Chunk = 1LL << 22;  // posm2 L1pos per launch
enum : int { kPass = 0, kOutside = 1, kOverlap = 2, kFar = 3 };

struct FilterArgs {
  FieldView field;
  const float* y;
  const int* heavy;
  int probes;
  float c[3], lo[3], hi[3], cd, nearD;
  int maxOverlap;
  int n;
  const GpuDockPrep::FilterPose* poses;  // explicit poses
  // posm2
  const float* quats;
  const int* range;
  const long long* offset;
  int orients;
  double base[3], step;
  long long k0;
  // posm3 (rotation r0 + j: quats[4j], base[6j], extent box[6j])
  const double* rotBase;
  const double* delta;
  const float* box;
  int nt;
  long long r0;
  float fitLo[3], fitHi[3];
  long long* pass;
  unsigned long long* passCount;
  unsigned long long* counts;  // outside, overlap, far
};

// MobmDocker::quickFilter; InBox (MobmDocker::quickFilterInBox): every probe atom is known to
// be inside the box, so the result is decided at the overlap limit.
template <bool InBox>
__device__ int classify(const FilterArgs& a, const float t[3], const float q[4]) {
  int overlap = 0;
  bool near = false;
  for (int k = 0; k < a.probes; ++k) {
    const float v[3] = {a.y[3 * k] - a.c[0], a.y[3 * k + 1] - a.c[1], a.y[3 * k + 2] - a.c[2]};
    float r[3];
    dock::rotate(q, v, r);
    const float w[3] = {r[0] + t[0], r[1] + t[1], r[2] + t[2]};
    if constexpr (!InBox)
      for (int d = 0; d < 3; ++d)
        if (w[d] < a.lo[d] || w[d] > a.hi[d]) return kOutside;
    if (!a.heavy[k]) continue;
    const float dd = field::distance<float>(a.field, w, nullptr);
    if (dd < a.cd && ++overlap >= a.maxOverlap && InBox) return kOverlap;
    if (dd <= a.nearD) near = true;
  }
  if (overlap >= a.maxOverlap) return kOverlap;
  if (!near) return kFar;
  return kPass;
}

// posm3: extent of the rotated probe atoms per rotation (MobmDocker::rotBox).
__global__ void rotBoxKernel(FilterArgs a, int rotations) {
  const int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= rotations) return;
  const float* q = a.quats + 4 * j;
  float lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
  for (int k = 0; k < a.probes; ++k) {
    const float v[3] = {a.y[3 * k] - a.c[0], a.y[3 * k + 1] - a.c[1], a.y[3 * k + 2] - a.c[2]};
    float r[3];
    dock::rotate(q, v, r);
    for (int d = 0; d < 3; ++d) {
      lo[d] = fminf(lo[d], r[d]);
      hi[d] = fmaxf(hi[d], r[d]);
    }
  }
  float* b = const_cast<float*>(a.box) + 6 * j;
  for (int d = 0; d < 3; ++d) {
    b[d] = lo[d];
    b[3 + d] = hi[d];
  }
}

enum : int { kExplicit = 0, kGenPosm2 = 1, kGenPosm3 = 2 };

template <int Mode>
__global__ void __launch_bounds__(kFilterThreads) filterKernel(FilterArgs a) {
  __shared__ unsigned int cnt[4];
  if (threadIdx.x < 4) cnt[threadIdx.x] = 0;
  __syncthreads();
  const long long i = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < a.n) {
    float t[3], q[4];
    long long k;
    int code = -1;
    if constexpr (Mode == kGenPosm3) {
      // MobmDocker::posm3Pose: t = anchorS + delta - rm, then fitsBox from the rotated extent.
      k = a.k0 + i;
      const int j = static_cast<int>(k / a.nt - a.r0), tr = static_cast<int>(k % a.nt);
      const double* bs = a.rotBase + 6 * j;
      const double* dl = a.delta + 3 * tr;
      for (int c = 0; c < 3; ++c) t[c] = static_cast<float>(bs[c] + dl[c] - bs[3 + c]);
      for (int c = 0; c < 4; ++c) q[c] = a.quats[4 * j + c];
      const float* b = a.box + 6 * j;
      for (int d = 0; d < 3; ++d)
        if (b[d] + t[d] < a.fitLo[d] || b[3 + d] + t[d] > a.fitHi[d]) code = kOutside;
    } else if constexpr (Mode == kGenPosm2) {
      k = a.k0 + i;
      // MobmDocker::l1Pose (posm2): orientation j covers [offset[j], offset[j + 1]).
      int lo = 0, hi = a.orients;  // first offset > k is at index hi
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if (a.offset[mid] <= k) lo = mid + 1;
        else hi = mid;
      }
      const int j = lo - 1;
      long long pos = k - a.offset[j];
      const int* rg = a.range + 6 * j;
      for (int c = 0; c < 3; ++c) {
        const long long w = rg[3 + c] - rg[c] + 1;
        const long long idx = rg[c] + pos % w;
        pos /= w;
        t[c] = static_cast<float>(a.base[c] + a.step * static_cast<double>(idx));
      }
      for (int c = 0; c < 4; ++c) q[c] = a.quats[4 * j + c];
    } else {
      const GpuDockPrep::FilterPose& p = a.poses[i];
      k = p.k;
      for (int c = 0; c < 3; ++c) t[c] = p.t[c];
      for (int c = 0; c < 4; ++c) q[c] = p.q[c];
    }
    if (code < 0) code = classify<Mode == kGenPosm3>(a, t, q);
    atomicAdd(&cnt[code], 1u);
    if (code == kPass) a.pass[atomicAdd(a.passCount, 1ull)] = k;
  }
  __syncthreads();
  if (threadIdx.x >= 1 && threadIdx.x < 4 && cnt[threadIdx.x]) atomicAdd(a.counts + threadIdx.x - 1, cnt[threadIdx.x]);
}

// posm2_fft kernels (fftCorrelate, fftSelect).
// zero: clears the touched points instead.
__global__ void spreadKernel(int n, int na, const int* atoms, const float* wts, int probes, const double* f, int px,
                             int py, int pz, float* real, std::size_t N, int zero) {
  const long long id = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (id >= static_cast<long long>(n) * na) return;
  const int b = static_cast<int>(id / na), j = static_cast<int>(id % na);
  const double* fa = f + (static_cast<std::size_t>(b) * probes + atoms[j]) * 3;
  const int p[3] = {px, py, pz};
  int base[3];
  float u[3];
  for (int d = 0; d < 3; ++d) {
    const double fl = floor(fa[d]);
    base[d] = ((static_cast<int>(fl) % p[d]) + p[d]) % p[d];
    u[d] = static_cast<float>(fa[d] - fl);
  }
  const float wt = wts[j];
  float* r = real + static_cast<std::size_t>(b) * N;
  for (int dz = 0; dz < 2; ++dz) {
    const int k = (base[2] + dz) % pz;
    const float wz = dz ? u[2] : 1.0f - u[2];
    for (int dy = 0; dy < 2; ++dy) {
      const int jj = (base[1] + dy) % py;
      const float wy = dy ? u[1] : 1.0f - u[1];
      for (int dx = 0; dx < 2; ++dx) {
        const int i = (base[0] + dx) % px;
        const float wx = dx ? u[0] : 1.0f - u[0];
        float* v = r + (static_cast<std::size_t>(k) * py + jj) * px + i;
        if (zero) *v = 0.0f;
        else atomicAdd(v, wt * wx * wy * wz);
      }
    }
  }
}

// acc (+)= conj(spec) S; blockIdx.y: orientation.
__global__ void macKernel(std::size_t Nc, const cufftComplex* spec, const cufftComplex* S, cufftComplex* acc, int first) {
  const std::size_t k = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (k >= Nc) return;
  const std::size_t id = static_cast<std::size_t>(blockIdx.y) * Nc + k;
  const cufftComplex a = spec[id], s = S[k];
  const float re = a.x * s.x + a.y * s.y, im = a.x * s.y - a.y * s.x;
  if (first) acc[id] = {re, im};
  else acc[id] = {acc[id].x + re, acc[id].y + im};
}

// Direct probe spectra (spec 15.1.5): the trilinear spreading factorizes per axis, so the
// transform of atom a is P_a(k) = prod_d [(1 - u_d) e(k_d b_d) + u_d e(k_d (b_d + 1))] with
// e(m) = exp(-2 pi i m / p_d). tables: per orientation and atom, hx + py + pz factors.
__global__ void tableKernel(int n, int probes, const double* f, int px, int py, int pz, cufftComplex* tables) {
  const int hx = px / 2 + 1, len = hx + py + pz;
  const long long id = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (id >= static_cast<long long>(n) * probes * len) return;
  const int k0 = static_cast<int>(id % len);
  const long long ba = id / len;
  const double* fa = f + ba * 3;
  int d, k;
  if (k0 < hx) d = 0, k = k0;
  else if (k0 < hx + py) d = 1, k = k0 - hx;
  else d = 2, k = k0 - hx - py;
  const int p = d == 0 ? px : d == 1 ? py : pz;
  const double fl = floor(fa[d]);
  const int b = ((static_cast<int>(fl) % p) + p) % p;
  const float u = static_cast<float>(fa[d] - fl);
  double s0, c0, s1, c1;
  sincospi(-2.0 * static_cast<double>((static_cast<long long>(k) * b) % p) / p, &s0, &c0);
  sincospi(-2.0 * static_cast<double>((static_cast<long long>(k) * ((b + 1) % p)) % p) / p, &s1, &c1);
  tables[id] = {static_cast<float>((1 - u) * c0 + u * c1), static_cast<float>((1 - u) * s0 + u * s1)};
}

// acc[b][k] = sum_a conj(P_a(k)) (S_map(a)(k) + q_a S_elec(k) + |q_a| S_desolv(k)); atoms grouped by map.
// Thread: one (kx, ky) and kDirectZ consecutive kz; blockIdx.y: kDirectOrients orientations.
constexpr int kDirectOrients = 4;
constexpr int kDirectZ = 4;
constexpr int kDirectAtomsPerFft = 40;
__global__ void __launch_bounds__(128) directKernel(int n, int probes, int px, int py, int pz, const cufftComplex* tables,
                                                    const cufftComplex* spectra, std::size_t Nc, int nt, int groups,
                                                    const int* groupMap, const int* groupStart, const int* atom,
                                                    const float* q, cufftComplex* acc) {
  const int hx = px / 2 + 1, len = hx + py + pz, zg = (pz + kDirectZ - 1) / kDirectZ;
  const long long id = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (id >= static_cast<long long>(hx) * py * zg) return;
  const int kx = static_cast<int>(id % hx), ky = static_cast<int>(id / hx % py), kz0 = static_cast<int>(id / hx / py) * kDirectZ;
  const int nz = min(kDirectZ, pz - kz0);
  std::size_t kk[kDirectZ];
  cufftComplex se[kDirectZ], sd[kDirectZ];
  for (int z = 0; z < kDirectZ; ++z) {
    kk[z] = (static_cast<std::size_t>(min(kz0 + z, pz - 1)) * py + ky) * hx + kx;
    se[z] = spectra[static_cast<std::size_t>(nt) * Nc + kk[z]];
    sd[z] = spectra[static_cast<std::size_t>(nt + 1) * Nc + kk[z]];
  }
  const int b0 = blockIdx.y * kDirectOrients, b1 = min(n, b0 + kDirectOrients);
  float ar[kDirectOrients][kDirectZ], ai[kDirectOrients][kDirectZ];
  for (int b = 0; b < kDirectOrients; ++b)
    for (int z = 0; z < kDirectZ; ++z) ar[b][z] = ai[b][z] = 0;
  for (int g = 0; g < groups; ++g) {
    cufftComplex sm[kDirectZ];
    for (int z = 0; z < kDirectZ; ++z) sm[z] = spectra[static_cast<std::size_t>(groupMap[g]) * Nc + kk[z]];
    for (int j = groupStart[g]; j < groupStart[g + 1]; ++j) {
      const float qa = q[j], aq = fabsf(qa);
      float cr[kDirectZ], ci[kDirectZ];
      for (int z = 0; z < kDirectZ; ++z) {
        cr[z] = sm[z].x + qa * se[z].x + aq * sd[z].x;
        ci[z] = sm[z].y + qa * se[z].y + aq * sd[z].y;
      }
      const int a = atom[j];
#pragma unroll
      for (int bb = 0; bb < kDirectOrients; ++bb) {
        const int b = b0 + bb;
        if (b >= b1) break;
        const cufftComplex* T = tables + (static_cast<std::size_t>(b) * probes + a) * len;
        const cufftComplex tx = T[kx], ty = T[hx + ky];
        const float xr = tx.x * ty.x - tx.y * ty.y, xi = tx.x * ty.y + tx.y * ty.x;
#pragma unroll
        for (int z = 0; z < kDirectZ; ++z) {
          const cufftComplex tz = T[hx + py + min(kz0 + z, pz - 1)];
          const float pr = xr * tz.x - xi * tz.y, pi = xr * tz.y + xi * tz.x;
          // conj(P) C
          ar[bb][z] += pr * cr[z] + pi * ci[z];
          ai[bb][z] += pr * ci[z] - pi * cr[z];
        }
      }
    }
  }
  for (int b = b0; b < b1; ++b)
    for (int z = 0; z < nz; ++z) acc[static_cast<std::size_t>(b) * Nc + kk[z]] = {ar[b - b0][z], ai[b - b0][z]};
}

// Several spectra sets (spec 15.1.22): acc[s][b][k] = sum_c conj(P[c][b][k]) S[s][chan[c]][k], P the probe
// spectra (channel-major, B orientations each). Thread: one k, kMacB orientations, kMacS sets at a time.
constexpr int kMacB = 8;
constexpr int kMacS = 4;
__global__ void macMultiKernel(std::size_t Nc, int n, int B, int nc, const int* chan, const cufftComplex* P,
                               const cufftComplex* S, int channels, int sets, cufftComplex* acc) {
  const std::size_t k = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (k >= Nc) return;
  const int b0 = blockIdx.y * kMacB, nb = min(kMacB, n - b0);
  for (int s0 = 0; s0 < sets; s0 += kMacS) {
    const int ns = min(kMacS, sets - s0);
    float ar[kMacS][kMacB], ai[kMacS][kMacB];
#pragma unroll
    for (int s = 0; s < kMacS; ++s)
#pragma unroll
      for (int b = 0; b < kMacB; ++b) ar[s][b] = ai[s][b] = 0;
    for (int c = 0; c < nc; ++c) {
      cufftComplex pv[kMacB];
#pragma unroll
      for (int b = 0; b < kMacB; ++b)
        pv[b] = b < nb ? P[(static_cast<std::size_t>(c) * B + b0 + b) * Nc + k] : cufftComplex{0, 0};
      const std::size_t m = static_cast<std::size_t>(chan[c]);
#pragma unroll
      for (int s = 0; s < kMacS; ++s) {
        if (s >= ns) break;
        const cufftComplex sv = S[(static_cast<std::size_t>(s0 + s) * channels + m) * Nc + k];
#pragma unroll
        for (int b = 0; b < kMacB; ++b) {
          ar[s][b] += pv[b].x * sv.x + pv[b].y * sv.y;
          ai[s][b] += pv[b].x * sv.y - pv[b].y * sv.x;
        }
      }
    }
    for (int s = 0; s < ns; ++s)
      for (int b = 0; b < nb; ++b)
        acc[(static_cast<std::size_t>(s0 + s) * B + b0 + b) * Nc + k] = {ar[s][b], ai[s][b]};
  }
}

struct Cand {
  float e;
  unsigned long long i;
};
__device__ inline bool less(const Cand& a, const Cand& b) { return a.e < b.e || (a.e == b.e && a.i < b.i); }

constexpr int kSelThreads = 256;
constexpr int kSelMax = 16;     // fftTopPerRot handled on the device
constexpr int kSelSlices = 16;  // blocks per orientation

// Block (slice, orientation): the top lowest candidates of its share of the translations.
__global__ void __launch_bounds__(kSelThreads) selectKernel(const float* real, std::size_t N, int px, int py,
                                                              const int* range, float scale, int top, Cand* part) {
  const int b = blockIdx.y, tid = threadIdx.x;
  const int* rg = range + 6 * b;
  const int wx = rg[3] - rg[0] + 1, wy = rg[4] - rg[1] + 1, wz = rg[5] - rg[2] + 1;
  const int count = wx * wy * wz;
  const float* r = real + static_cast<std::size_t>(b) * N;
  Cand loc[kSelMax];
  int m = 0;
  for (int t = blockIdx.x * kSelThreads + tid; t < count; t += kSelSlices * kSelThreads) {
    const int i = rg[0] + t % wx, j = rg[1] + t / wx % wy, k = rg[2] + t / wx / wy;
    const unsigned long long idx = (static_cast<unsigned long long>(k) * py + j) * px + i;
    const Cand c{r[idx] * scale, idx};
    if (m == top && !less(c, loc[m - 1])) continue;
    int p = m < top ? m++ : m - 1;
    while (p > 0 && less(c, loc[p - 1])) {
      loc[p] = loc[p - 1];
      --p;
    }
    loc[p] = c;
  }
  __shared__ Cand sc[kSelThreads];
  __shared__ int sw[kSelThreads];
  int head = 0;
  for (int round = 0; round < top; ++round) {
    sc[tid] = head < m ? loc[head] : Cand{INFINITY, ~0ull};
    sw[tid] = tid;
    __syncthreads();
    for (int s = kSelThreads / 2; s > 0; s /= 2) {
      if (tid < s && less(sc[tid + s], sc[tid])) {
        sc[tid] = sc[tid + s];
        sw[tid] = sw[tid + s];
      }
      __syncthreads();
    }
    if (tid == 0) part[(static_cast<std::size_t>(b) * kSelSlices + blockIdx.x) * top + round] = sc[0];
    if (tid == sw[0] && head < m) ++head;
    __syncthreads();
  }
}

// One thread per orientation: merges the slices' candidates.
__global__ void mergeKernel(int n, int top, const Cand* part, Cand* out) {
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= n) return;
  const Cand* p = part + static_cast<std::size_t>(b) * kSelSlices * top;
  int head[kSelSlices];
  for (int s = 0; s < kSelSlices; ++s) head[s] = 0;
  auto cand = [&](int s) { return head[s] < top ? p[s * top + head[s]] : Cand{INFINITY, ~0ull}; };
  for (int r = 0; r < top; ++r) {
    int best = 0;
    for (int s = 1; s < kSelSlices; ++s)
      if (less(cand(s), cand(best))) best = s;
    out[static_cast<std::size_t>(b) * top + r] = cand(best);
    if (head[best] < top) ++head[best];
  }
}

}  // namespace

struct GpuDockPrep::Impl {
  int device = 0;
  // Highest priority: not delayed by the delta grids of the next StaM conformation (spec 15.1.20).
  cudaStream_t stream = nullptr;
  // filter
  DevBuf<float> y, quats;
  DevBuf<int> heavy, range;
  DevBuf<long long> offset, pass;
  DevBuf<FilterPose> poses;
  DevBuf<unsigned long long> counters;  // passCount, outside, overlap, far
  DevBuf<double> base, delta;           // posm3
  DevBuf<float> box;
  // posm2_fft; several spectra sets (StaM conformations, spec 15.1.22)
  std::vector<std::uint64_t> fftIds;
  int sets = 0, planSets = 0;
  std::array<int, 3> pad{};
  std::size_t N = 0, Nc = 0;
  int channels = 0;
  DevBuf<cufftComplex> spectra;
  int planBatch = 0;
  cufftHandle r2c = 0, c2r = 0;
  std::size_t workBytes = 0;
  DevBuf<float> real;
  DevBuf<cufftComplex> spec, acc, tables, probeSpec;
  DevBuf<char> work;
  DevBuf<double> f;
  DevBuf<int> atoms, ranges, groupMap, groupStart;
  DevBuf<float> wts;
  DevBuf<Cand> part, outC;

  ~Impl() {
    destroyPlans();
    if (stream) {
      cudaStreamSynchronize(stream);
      cudaStreamDestroy(stream);
    }
  }
  // Before work on the stream: the copies of the legacy default stream (e.g. GpuField) are done.
  void begin() {
    check(cudaSetDevice(device), "cudaSetDevice");
    check(cudaStreamSynchronize(0), "cudaStreamSynchronize");
  }
  Spectra spectraMode = Spectra::Auto;
  // Direct spectra cost about probes × Nc per orientation, the FFTs channels × N log N.
  bool directSpectra(int probes, int channels) const {
    if (spectraMode != Spectra::Auto) return spectraMode == Spectra::Direct;
    return probes <= channels * kDirectAtomsPerFft;
  }
  void releaseFftWork() {
    real.release();
    spec.release();
    acc.release();
    probeSpec.release();
    tables.release();
    work.release();
  }
  void destroyPlans() {
    if (planBatch) {
      cufftDestroy(r2c);
      cufftDestroy(c2r);
    }
    planBatch = 0;
    planSets = 0;
    workBytes = 0;
  }
  void makePlans(int b) {
    destroyPlans();
    int dims[3] = {pad[2], pad[1], pad[0]};
    std::size_t w1 = 0, w2 = 0;
    checkFft(cufftCreate(&r2c), "cufftCreate");
    checkFft(cufftCreate(&c2r), "cufftCreate");
    planBatch = b;
    checkFft(cufftSetAutoAllocation(r2c, 0), "cufftSetAutoAllocation");
    checkFft(cufftSetAutoAllocation(c2r, 0), "cufftSetAutoAllocation");
    checkFft(cufftMakePlanMany(r2c, 3, dims, nullptr, 1, 0, nullptr, 1, 0, CUFFT_R2C, b, &w1), "plan R2C");
    checkFft(cufftMakePlanMany(c2r, 3, dims, nullptr, 1, 0, nullptr, 1, 0, CUFFT_C2R, b, &w2), "plan C2R");
    checkFft(cufftSetStream(r2c, stream), "cufftSetStream");
    checkFft(cufftSetStream(c2r, stream), "cufftSetStream");
    workBytes = std::max(w1, w2);
  }

  void runFilter(FilterArgs a, int mode, std::vector<long long>& out, FilterCounts& cnt) {
    check(cudaMemsetAsync(counters.p, 0, 4 * sizeof(unsigned long long), stream), "memset");
    a.passCount = counters.p;
    a.counts = counters.p + 1;
    a.pass = pass.p;
    const int blocks = (a.n + kFilterThreads - 1) / kFilterThreads;
    if (mode == kGenPosm3) filterKernel<kGenPosm3><<<blocks, kFilterThreads, 0, stream>>>(a);
    else if (mode == kGenPosm2) filterKernel<kGenPosm2><<<blocks, kFilterThreads, 0, stream>>>(a);
    else filterKernel<kExplicit><<<blocks, kFilterThreads, 0, stream>>>(a);
    check(cudaGetLastError(), "filterKernel launch");
    unsigned long long h[4];
    check(cudaMemcpyAsync(h, counters.p, sizeof h, cudaMemcpyDeviceToHost, stream), "download counters");
    const std::size_t old = out.size();
    out.resize(old + h[0]);
    if (h[0]) check(cudaMemcpyAsync(out.data() + old, pass.p, h[0] * sizeof(long long), cudaMemcpyDeviceToHost, stream), "download");
    cnt.outsideBox += static_cast<long long>(h[1]);
    cnt.overlap += static_cast<long long>(h[2]);
    cnt.far += static_cast<long long>(h[3]);
  }

  FilterArgs args(const GpuField& field, const FilterSetup& s) {
    y.upload(s.y, stream);
    heavy.upload(s.heavy, stream);
    counters.reserve(4);
    FilterArgs a{};
    a.field = field.deviceView();
    a.y = y.p;
    a.heavy = heavy.p;
    a.probes = static_cast<int>(s.heavy.size());
    for (int d = 0; d < 3; ++d) {
      a.c[d] = s.center[d];
      a.lo[d] = s.boxLo[d];
      a.hi[d] = s.boxHi[d];
    }
    a.cd = s.clashDistance;
    a.nearD = s.nearDistance;
    a.maxOverlap = s.maxOverlap;
    return a;
  }
};

GpuDockPrep::GpuDockPrep(int device) : impl_(std::make_unique<Impl>()) {
  impl_->device = device;
  check(cudaSetDevice(device), "cudaSetDevice");
  int least = 0, greatest = 0;
  check(cudaDeviceGetStreamPriorityRange(&least, &greatest), "cudaDeviceGetStreamPriorityRange");
  check(cudaStreamCreateWithPriority(&impl_->stream, cudaStreamNonBlocking, greatest), "cudaStreamCreate");
}
GpuDockPrep::~GpuDockPrep() = default;
void GpuDockPrep::setSpectra(Spectra s) { impl_->spectraMode = s; }

void GpuDockPrep::filter(const GpuField& field, const FilterSetup& s, int n, const FilterPose* poses,
                         std::vector<long long>& pass, FilterCounts& cnt) {
  if (n <= 0) return;
  Impl& d = *impl_;
  d.begin();
  FilterArgs a = d.args(field, s);
  d.poses.upload(poses, static_cast<std::size_t>(n), d.stream);
  d.pass.reserve(static_cast<std::size_t>(n));
  a.poses = d.poses.p;
  a.n = n;
  d.runFilter(a, kExplicit, pass, cnt);
}

void GpuDockPrep::filterPosm2(const GpuField& field, const FilterSetup& s, const Posm2Setup& p2, long long k0,
                              long long k1, std::vector<long long>& pass, FilterCounts& cnt) {
  if (k1 <= k0) return;
  Impl& d = *impl_;
  d.begin();
  FilterArgs a = d.args(field, s);
  d.quats.upload(p2.quats, d.stream);
  d.range.upload(p2.range, d.stream);
  d.offset.upload(p2.offset, d.stream);
  d.pass.reserve(static_cast<std::size_t>(std::min(kPosm2Chunk, k1 - k0)));
  a.quats = d.quats.p;
  a.range = d.range.p;
  a.offset = d.offset.p;
  a.orients = static_cast<int>(p2.quats.size() / 4);
  for (int c = 0; c < 3; ++c) a.base[c] = p2.base[c];
  a.step = p2.step;
  for (long long b = k0; b < k1; b += kPosm2Chunk) {
    a.k0 = b;
    a.n = static_cast<int>(std::min(kPosm2Chunk, k1 - b));
    d.runFilter(a, kGenPosm2, pass, cnt);
  }
}

void GpuDockPrep::filterPosm3(const GpuField& field, const FilterSetup& s, const Posm3Setup& p3, long long k0,
                              long long k1, std::vector<long long>& pass, FilterCounts& cnt) {
  if (k1 <= k0) return;
  Impl& d = *impl_;
  d.begin();
  FilterArgs a = d.args(field, s);
  const int rotations = static_cast<int>(p3.quats.size() / 4);
  d.quats.upload(p3.quats, d.stream);
  d.base.upload(p3.base, d.stream);
  d.delta.upload(p3.delta, d.stream);
  d.box.reserve(6 * static_cast<std::size_t>(rotations));
  d.pass.reserve(static_cast<std::size_t>(std::min(kPosm2Chunk, k1 - k0)));
  a.quats = d.quats.p;
  a.rotBase = d.base.p;
  a.delta = d.delta.p;
  a.box = d.box.p;
  a.nt = static_cast<int>(p3.delta.size() / 3);
  a.r0 = p3.r0;
  for (int c = 0; c < 3; ++c) {
    a.fitLo[c] = p3.fitLo[c];
    a.fitHi[c] = p3.fitHi[c];
  }
  rotBoxKernel<<<(rotations + kFilterThreads - 1) / kFilterThreads, kFilterThreads, 0, d.stream>>>(a, rotations);
  check(cudaGetLastError(), "rotBoxKernel launch");
  for (long long b = k0; b < k1; b += kPosm2Chunk) {
    a.k0 = b;
    a.n = static_cast<int>(std::min(kPosm2Chunk, k1 - b));
    d.runFilter(a, kGenPosm3, pass, cnt);
  }
}

void GpuDockPrep::setFft(const FftStam& s) {
  const FftStam* p = &s;
  setFft(std::span<const FftStam* const>(&p, 1));
}

void GpuDockPrep::setFft(std::span<const FftStam* const> s) {
  Impl& d = *impl_;
  if (s.empty()) throw Error("setFft: no spectra");
  std::vector<std::uint64_t> ids;
  for (const FftStam* f : s) ids.push_back(f->id);
  if (std::find(ids.begin(), ids.end(), std::uint64_t{0}) == ids.end() && ids == d.fftIds) return;
  const FftStam& s0 = *s[0];
  for (const FftStam* f : s)
    if (f->pad != s0.pad || f->spectra.size() != s0.spectra.size()) throw Error("setFft: spectra of different grid geometries");
  d.begin();
  if (s0.pad != d.pad) d.destroyPlans();
  d.pad = s0.pad;
  d.N = s0.realSize();
  d.Nc = s0.complexSize();
  d.channels = static_cast<int>(s0.spectra.size());
  d.sets = static_cast<int>(s.size());
  d.spectra.release();
  d.spectra.reserve(static_cast<std::size_t>(d.sets) * d.channels * d.Nc);
  for (int g = 0; g < d.sets; ++g)
    for (int m = 0; m < d.channels; ++m)
      check(cudaMemcpyAsync(d.spectra.p + (static_cast<std::size_t>(g) * d.channels + m) * d.Nc, s[g]->spectra[m].data(),
                            d.Nc * sizeof(cufftComplex), cudaMemcpyHostToDevice, d.stream),
            "upload spectrum");
  d.fftIds = std::move(ids);
}

int GpuDockPrep::fftBatch() {
  Impl& d = *impl_;
  if (d.fftIds.empty()) return 0;
  if (d.planBatch && d.planSets == d.sets) return d.planBatch;
  d.begin();
  d.destroyPlans();
  d.releaseFftWork();
  // One set: real, probe spectrum and acc per orientation; several: real, the probe spectra and an acc per set.
  const std::size_t per = d.sets == 1 ? d.N * sizeof(float) + 2 * d.Nc * sizeof(cufftComplex)
                                      : d.N * sizeof(float) + static_cast<std::size_t>(d.channels + d.sets) * d.Nc * sizeof(cufftComplex);
  const double avail = 0.5 * static_cast<double>(gpuFreeMem(d.device));
  int b = static_cast<int>(std::clamp(avail / static_cast<double>(per), 1.0, 64.0));
  for (;;) {
    d.makePlans(b);
    if (b == 1 || static_cast<double>(per * b + d.workBytes) <= avail) break;
    b /= 2;
  }
  d.planSets = d.sets;
  return b;
}

void GpuDockPrep::fftTop(int n, int probes, const double* f, const int* map, const double* q,
                         const std::array<int, 3>* lo, const std::array<int, 3>* hi, int top,
                         std::pair<float, std::size_t>* out) {
  if (n <= 0) return;
  Impl& d = *impl_;
  d.begin();
  if (d.fftIds.empty()) throw Error("fftTop: no spectra");
  if (n > fftBatch()) throw Error("fftTop: batch too large");
  const int B = d.planBatch, G = d.sets;
  d.real.reserve(static_cast<std::size_t>(B) * d.N);
  d.acc.reserve(static_cast<std::size_t>(G) * B * d.Nc);
  d.work.reserve(d.workBytes);
  checkFft(cufftSetWorkArea(d.r2c, d.work.p), "cufftSetWorkArea");
  checkFft(cufftSetWorkArea(d.c2r, d.work.p), "cufftSetWorkArea");
  d.f.upload(f, static_cast<std::size_t>(n) * probes * 3, d.stream);
  std::vector<int> rg(6 * static_cast<std::size_t>(n));
  for (int b = 0; b < n; ++b)
    for (int c = 0; c < 3; ++c) {
      rg[6 * b + c] = lo[b][c];
      rg[6 * b + 3 + c] = hi[b][c];
    }
  d.ranges.upload(rg, d.stream);
  // Channels as fftCorrelate: the type maps with probe atoms, then elec and desolv when charged.
  const int nt = d.channels - 2;
  std::vector<std::pair<int, std::vector<std::pair<int, float>>>> chans;
  for (int m = 0; m < nt; ++m) {
    std::vector<std::pair<int, float>> v;
    for (int a = 0; a < probes; ++a)
      if (map[a] == m) v.push_back({a, 1.0f});
    if (!v.empty()) chans.push_back({m, std::move(v)});
  }
  std::vector<std::pair<int, float>> ve, vd;
  for (int a = 0; a < probes; ++a)
    if (q[a] != 0) {
      ve.push_back({a, static_cast<float>(q[a])});
      vd.push_back({a, static_cast<float>(std::abs(q[a]))});
    }
  if (!ve.empty()) {
    chans.push_back({nt, std::move(ve)});
    chans.push_back({nt + 1, std::move(vd)});
  }
  const std::size_t outN = static_cast<std::size_t>(n) * top;
  if (chans.empty()) {
    for (std::size_t i = 0; i < G * outN; ++i) out[i] = {0.0f, SIZE_MAX};
    return;
  }
  const unsigned nbx = static_cast<unsigned>((d.Nc + 255) / 256);
  if (G > 1) {
    // Several sets (spec 15.1.22): the probe spectrum of each channel once (FFT of the spread
    // atoms), then the products with every set.
    const int nc = static_cast<int>(chans.size());
    d.probeSpec.reserve(static_cast<std::size_t>(nc) * B * d.Nc);
    check(cudaMemsetAsync(d.real.p, 0, static_cast<std::size_t>(B) * d.N * sizeof(float), d.stream), "memset");
    std::vector<int> cm;
    for (int c = 0; c < nc; ++c) {
      const auto& list = chans[c].second;
      cm.push_back(chans[c].first);
      std::vector<int> at;
      std::vector<float> wt;
      for (const auto& [a, w] : list) {
        at.push_back(a);
        wt.push_back(w);
      }
      d.atoms.upload(at, d.stream);
      d.wts.upload(wt, d.stream);
      const long long ns = static_cast<long long>(n) * static_cast<long long>(at.size());
      const unsigned nsb = static_cast<unsigned>((ns + 255) / 256);
      spreadKernel<<<nsb, 256, 0, d.stream>>>(n, static_cast<int>(at.size()), d.atoms.p, d.wts.p, probes, d.f.p, d.pad[0], d.pad[1],
                                 d.pad[2], d.real.p, d.N, 0);
      check(cudaGetLastError(), "spreadKernel launch");
      checkFft(cufftExecR2C(d.r2c, d.real.p, d.probeSpec.p + static_cast<std::size_t>(c) * B * d.Nc), "R2C");
      spreadKernel<<<nsb, 256, 0, d.stream>>>(n, static_cast<int>(at.size()), d.atoms.p, d.wts.p, probes, d.f.p, d.pad[0], d.pad[1],
                                 d.pad[2], d.real.p, d.N, 1);
      check(cudaGetLastError(), "spreadKernel (clear) launch");
    }
    d.groupMap.upload(cm, d.stream);
    macMultiKernel<<<dim3(nbx, static_cast<unsigned>((n + kMacB - 1) / kMacB)), 256, 0, d.stream>>>(
        d.Nc, n, B, nc, d.groupMap.p, d.probeSpec.p, d.spectra.p, d.channels, G, d.acc.p);
    check(cudaGetLastError(), "macMultiKernel launch");
  } else if (d.directSpectra(probes, static_cast<int>(chans.size()))) {
    // Direct probe spectra, all channels fused.
    const int len = d.pad[0] / 2 + 1 + d.pad[1] + d.pad[2];
    const std::size_t nTab = static_cast<std::size_t>(n) * probes * len;
    d.tables.reserve(nTab);
    tableKernel<<<static_cast<unsigned>((nTab + 255) / 256), 256, 0, d.stream>>>(n, probes, d.f.p, d.pad[0], d.pad[1], d.pad[2],
                                                                     d.tables.p);
    check(cudaGetLastError(), "tableKernel launch");
    std::vector<int> gm, gs{0}, at;
    std::vector<float> qs;
    for (int m = 0; m < nt; ++m) {
      for (int a = 0; a < probes; ++a)
        if (map[a] == m) {
          at.push_back(a);
          qs.push_back(static_cast<float>(q[a]));
        }
      if (static_cast<int>(at.size()) > gs.back()) {
        gm.push_back(m);
        gs.push_back(static_cast<int>(at.size()));
      }
    }
    d.groupMap.upload(gm, d.stream);
    d.groupStart.upload(gs, d.stream);
    d.atoms.upload(at, d.stream);
    d.wts.upload(qs, d.stream);
    const long long nthr = static_cast<long long>(d.pad[0] / 2 + 1) * d.pad[1] * ((d.pad[2] + kDirectZ - 1) / kDirectZ);
    const dim3 grid(static_cast<unsigned>((nthr + 127) / 128), static_cast<unsigned>((n + kDirectOrients - 1) / kDirectOrients));
    directKernel<<<grid, 128, 0, d.stream>>>(n, probes, d.pad[0], d.pad[1], d.pad[2], d.tables.p, d.spectra.p, d.Nc, nt,
                                static_cast<int>(gm.size()), d.groupMap.p, d.groupStart.p, d.atoms.p, d.wts.p, d.acc.p);
    check(cudaGetLastError(), "directKernel launch");
  } else {
    // FFT of the spread atoms per channel; the touched points are cleared after each one.
    d.spec.reserve(static_cast<std::size_t>(B) * d.Nc);
    check(cudaMemsetAsync(d.real.p, 0, static_cast<std::size_t>(B) * d.N * sizeof(float), d.stream), "memset");
    bool first = true;
    for (const auto& [m, list] : chans) {
      std::vector<int> at;
      std::vector<float> wt;
      for (const auto& [a, w] : list) {
        at.push_back(a);
        wt.push_back(w);
      }
      d.atoms.upload(at, d.stream);
      d.wts.upload(wt, d.stream);
      const long long ns = static_cast<long long>(n) * static_cast<long long>(at.size());
      const unsigned nsb = static_cast<unsigned>((ns + 255) / 256);
      spreadKernel<<<nsb, 256, 0, d.stream>>>(n, static_cast<int>(at.size()), d.atoms.p, d.wts.p, probes, d.f.p, d.pad[0], d.pad[1],
                                 d.pad[2], d.real.p, d.N, 0);
      check(cudaGetLastError(), "spreadKernel launch");
      checkFft(cufftExecR2C(d.r2c, d.real.p, d.spec.p), "R2C");
      spreadKernel<<<nsb, 256, 0, d.stream>>>(n, static_cast<int>(at.size()), d.atoms.p, d.wts.p, probes, d.f.p, d.pad[0], d.pad[1],
                                 d.pad[2], d.real.p, d.N, 1);
      check(cudaGetLastError(), "spreadKernel (clear) launch");
      macKernel<<<dim3(nbx, static_cast<unsigned>(n)), 256, 0, d.stream>>>(d.Nc, d.spec.p, d.spectra.p + static_cast<std::size_t>(m) * d.Nc,
                                                              d.acc.p, first ? 1 : 0);
      check(cudaGetLastError(), "macKernel launch");
      first = false;
    }
  }
  const float scale = 1.0f / static_cast<float>(d.N);
  std::vector<Cand> c;
  if (top <= kSelMax) {
    d.part.reserve(outN * kSelSlices);
    d.outC.reserve(static_cast<std::size_t>(G) * outN);
    c.resize(static_cast<std::size_t>(G) * outN);
  }
  std::vector<float> e;
  std::vector<std::pair<float, std::size_t>> best;
  for (int g = 0; g < G; ++g) {
    cufftComplex* acc = d.acc.p + static_cast<std::size_t>(g) * B * d.Nc;
    // Orientations beyond n: their spectra are cleared (the plan transforms B of them).
    if (n < B)
      check(cudaMemsetAsync(acc + static_cast<std::size_t>(n) * d.Nc, 0, static_cast<std::size_t>(B - n) * d.Nc * sizeof(cufftComplex),
                            d.stream),
            "memset");
    checkFft(cufftExecC2R(d.c2r, acc, d.real.p), "C2R");
    if (top <= kSelMax) {
      selectKernel<<<dim3(kSelSlices, static_cast<unsigned>(n)), kSelThreads, 0, d.stream>>>(d.real.p, d.N, d.pad[0], d.pad[1],
                                                                               d.ranges.p, scale, top, d.part.p);
      check(cudaGetLastError(), "selectKernel launch");
      mergeKernel<<<(n + 63) / 64, 64, 0, d.stream>>>(n, top, d.part.p, d.outC.p + static_cast<std::size_t>(g) * outN);
      check(cudaGetLastError(), "mergeKernel launch");
      continue;
    }
    // Larger fftTopPerRot: selection on the host (as fftSelect).
    e.resize(d.N);
    for (int b = 0; b < n; ++b) {
      check(cudaMemcpyAsync(e.data(), d.real.p + static_cast<std::size_t>(b) * d.N, d.N * sizeof(float), cudaMemcpyDeviceToHost,
                            d.stream),
            "download");
      check(cudaStreamSynchronize(d.stream), "cudaStreamSynchronize");
      best.clear();
      for (int k = lo[b][2]; k <= hi[b][2]; ++k)
        for (int j = lo[b][1]; j <= hi[b][1]; ++j)
          for (int i = lo[b][0]; i <= hi[b][0]; ++i) {
            const std::size_t idx = (static_cast<std::size_t>(k) * d.pad[1] + j) * d.pad[0] + i;
            const std::pair<float, std::size_t> cv{e[idx] * scale, idx};
            if (best.size() < static_cast<std::size_t>(top)) {
              best.push_back(cv);
              std::push_heap(best.begin(), best.end());
            } else if (cv < best.front()) {
              std::pop_heap(best.begin(), best.end());
              best.back() = cv;
              std::push_heap(best.begin(), best.end());
            }
          }
      std::sort_heap(best.begin(), best.end());
      for (int r = 0; r < top; ++r)
        out[static_cast<std::size_t>(g) * outN + static_cast<std::size_t>(b) * top + r] =
            r < static_cast<int>(best.size()) ? best[r] : std::pair<float, std::size_t>{0.0f, SIZE_MAX};
    }
  }
  if (top <= kSelMax) {
    check(cudaMemcpyAsync(c.data(), d.outC.p, c.size() * sizeof(Cand), cudaMemcpyDeviceToHost, d.stream), "download");
    check(cudaStreamSynchronize(d.stream), "cudaStreamSynchronize");
    for (std::size_t i = 0; i < c.size(); ++i) out[i] = {c[i].e, c[i].i == ~0ull ? SIZE_MAX : static_cast<std::size_t>(c[i].i)};
  }
}

void GpuDockPrep::releaseWork() {
  Impl& d = *impl_;
  d.releaseFftWork();
  d.spectra.release();
  d.fftIds.clear();
  d.f.release();
  d.part.release();
  d.outC.release();
  d.poses.release();
  d.pass.release();
  d.box.release();
}

}  // namespace od
