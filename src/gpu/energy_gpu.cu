// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "energy_gpu.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <string>
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

// One block of NT threads per sample. shared: the coordinates and the cos and sin of each torsion
// are kept in dynamic shared memory (the gradient stays in global memory: double atomics on shared
// memory are emulated on sm_70); xs is then only used to return the coordinates (keepX).
template <class T, bool Dock, int NT>
__global__ void __launch_bounds__(NT) evaluateKernel(EnergyModelView m, const T* angles, T* xs, double* gs,
                                                           FlatTerms* terms, double* grad, DockView dv,
                                                           const DockPose* poses, DockTerms* dterms, double* grad6,
                                                           int shared, int keepX) {
  constexpr int kThreads = NT, kWarps = NT / 32;
  extern __shared__ __align__(16) unsigned char smem[];
  const int s = blockIdx.x;
  const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount);
  const T* ang = angles + static_cast<std::size_t>(s) * m.torsionCount;
  T* x;
  double* g;
  x = shared ? reinterpret_cast<T*>(smem) : xs + s * nc;
  g = grad ? gs + s * nc : nullptr;

  T* tcs = shared ? x + nc : nullptr;
  for (std::size_t i = tid; i < nc; i += kThreads) {
    x[i] = static_cast<T>(m.ref[i]);
    if (g) g[i] = 0;
  }
  if (tcs)
    for (int t = tid; t < m.torsionCount; t += kThreads) flat::torsionCosSin(m, t, ang, tcs[2 * t], tcs[2 * t + 1]);
  __syncthreads();

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

  double e[5] = {0, 0, 0, 0, 0};
  int clashPairs = 0;
  // Warp w takes pairs base + lane. The pairs are sorted by i: the i contributions of a warp are
  // summed per i (segmented shuffle sum) and added once by the first lane of each segment.
  for (int base = warp * 32; base < m.pairCount; base += kThreads) {
    const int p = base + lane;
    T d[3], r, dedr;
    const bool ok = p < m.pairCount && flat::pairTerm(m, p, x, e, clashPairs, d, r, dedr);
    if (!g) continue;
    double v[3] = {0, 0, 0};
    const int key = p < m.pairCount ? m.pairs[p].i : -1;
    if (ok && r > 0) {
      const T f = dedr / r;
      const int j = m.pairs[p].j;
      for (int c = 0; c < 3; ++c) {
        v[c] = f * d[c];
        atomicAdd(g + 3 * j + c, -v[c]);
      }
    }
    for (int o = 1; o < 32; o *= 2) {
      const int k = __shfl_down_sync(0xffffffffu, key, o);
      double nv[3];
      for (int c = 0; c < 3; ++c) nv[c] = __shfl_down_sync(0xffffffffu, v[c], o);
      if (lane + o < 32 && k == key)
        for (int c = 0; c < 3; ++c) v[c] += nv[c];
    }
    const int prev = __shfl_up_sync(0xffffffffu, key, 1);
    if (key >= 0 && (lane == 0 || prev != key) && (v[0] != 0 || v[1] != 0 || v[2] != 0))
      for (int c = 0; c < 3; ++c) atomicAdd(g + 3 * key + c, v[c]);
  }
  double closure = 0, maxDev = 0;
  for (int c = tid; c < m.closureCount; c += kThreads) {
    T d[3], f;
    flat::closureTerm(m, c, x, closure, maxDev, d, f);
    if (!g) continue;
    const FlatClosure& cl = m.closures[c];
    for (int a = 0; a < 3; ++a) {
      const double v = f * d[a];
      atomicAdd(g + 3 * cl.i + a, v);
      atomicAdd(g + 3 * cl.j + a, -v);
    }
  }

  constexpr int kVals = Dock ? 17 : 8;
  double v[kVals] = {e[0], e[1], e[2], e[3], e[4], closure, maxDev, static_cast<double>(clashPairs)};
  if constexpr (Dock) {
    // 8: inter, 9: interClash, 10: outside, 11–13: translation, 14–16: rotation gradient.
    const DockPose ps = poses[s];
    for (int k = tid; k < dv.probeCount; k += kThreads) {
      float gw[3], r[3], gy[3];
      if (!dock::probeTerm(dv, k, x, ps, v[8], v[9], g ? gw : nullptr, r)) {
        v[10] += 1;
        continue;
      }
      if (!g) continue;
      dock::probeGrad(ps, gw, r, v + 11, v + 14, gy);
      const int c = dv.probes[k].coord;
      for (int i = 0; i < 3; ++i) atomicAdd(g + 3 * c + i, static_cast<double>(gy[i]));
    }
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
  if (!g) return;

  for (int t = warp; t < m.torsionCount; t += kWarps) {
    T a[3], u[3];
    double sum = 0;
    if (flat::torsionUnit(m, t, x, a, u)) {
      const FlatTorsion& tt = m.torsions[t];
      for (int q = tt.movingBegin + lane; q < tt.movingEnd; q += 32) sum += flat::torsionGradTerm(x, g, m.moving[q], a, u);
    }
    sum = warpSum(sum);
    if (lane == 0) grad[static_cast<std::size_t>(s) * m.torsionCount + t] = sum;
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

struct GpuField::Impl {
  int device = 0;
  std::uint64_t id = 0;
  std::shared_ptr<const void> base, baseDist;  // grids of the uploaded reference maps
  DevBuf<float> maps, data;
  DevBuf<int> block, movedHeavy, cellStart, cellAtoms, heavyStart, heavyAtoms, pairM, pairN;
  DevBuf<FieldMoved> moved;
  DevBuf<double> solT, volT, pairAC, pairBD, pairReqm;
  FieldView view{};
};

GpuField::GpuField(int device) : impl_(std::make_unique<Impl>()) { impl_->device = device; }
GpuField::~GpuField() = default;
std::uint64_t GpuField::id() const { return impl_->id; }
const FieldView& GpuField::deviceView() const { return impl_->view; }

std::size_t GpuField::deviceBytes() const {
  const Impl& d = *impl_;
  return d.maps.n * sizeof(float) + d.data.n * sizeof(float) + d.block.n * sizeof(int) + d.moved.n * sizeof(FieldMoved);
}

void GpuField::set(const FieldModel& fm) {
  Impl& d = *impl_;
  if (fm.id != 0 && fm.id == d.id) return;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const int nm = fm.view.typeMapCount + 3;
  if (fm.base != d.base || fm.baseDist != d.baseDist) {
    d.base.reset();
    d.baseDist.reset();
    d.maps.release();
    d.maps.reserve(static_cast<std::size_t>(nm) * fm.points);
    for (int m = 0; m < nm; ++m)
      check(cudaMemcpy(d.maps.p + static_cast<std::size_t>(m) * fm.points, fm.view.maps[m], fm.points * sizeof(float),
                       cudaMemcpyHostToDevice),
            "upload map");
    d.base = fm.base;
    d.baseDist = fm.baseDist;
  }
  d.data.release();
  d.data.upload(fm.data);
  d.block.upload(fm.block);
  d.moved.upload(fm.moved);
  d.movedHeavy.upload(fm.movedHeavy);
  d.cellStart.upload(fm.cellStart);
  d.cellAtoms.upload(fm.cellAtoms);
  d.heavyStart.upload(fm.heavyStart);
  d.heavyAtoms.upload(fm.heavyAtoms);
  d.solT.upload(fm.solT);
  d.volT.upload(fm.volT);
  d.pairAC.upload(fm.pairAC);
  d.pairBD.upload(fm.pairBD);
  d.pairReqm.upload(fm.pairReqm);
  d.pairM.upload(fm.pairM);
  d.pairN.upload(fm.pairN);
  FieldView v = fm.view;
  for (int m = 0; m < nm; ++m) v.maps[m] = d.maps.p + static_cast<std::size_t>(m) * fm.points;
  v.data = d.data.p;
  v.block = d.block.p;
  v.moved = d.moved.p;
  v.movedHeavy = d.movedHeavy.p;
  v.cells.start = d.cellStart.p;
  v.cells.atoms = d.cellAtoms.p;
  v.heavyCells.start = d.heavyStart.p;
  v.heavyCells.atoms = d.heavyAtoms.p;
  v.solT = d.solT.p;
  v.volT = d.volT.p;
  v.pairAC = d.pairAC.p;
  v.pairBD = d.pairBD.p;
  v.pairReqm = d.pairReqm.p;
  v.pairM = d.pairM.p;
  v.pairN = d.pairN.p;
  d.view = v;
  d.id = fm.id;
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
// Dynamic shared memory for the coordinates and torsion cos and sin of one sample; 0: they do not fit.
template <class T>
std::size_t sharedBytes(const EnergyModelView& m) {
  const std::size_t b = (3 * static_cast<std::size_t>(m.coordCount) + 2 * static_cast<std::size_t>(m.torsionCount)) * sizeof(T);
  return b <= static_cast<std::size_t>(sharedLimit()) - 4096 ? b : 0;  // the kernel's own shared memory
}

template <class T, bool Dock, int NT>
void launchNT(cudaStream_t st, int n, std::size_t shm, bool keepX, const EnergyModelView& m, const T* angles, T* xs,
              double* gs, FlatTerms* terms, double* grad, const DockView& dv, const DockPose* poses, DockTerms* dterms,
              double* grad6) {
  if (shm > 48 * 1024)
    check(cudaFuncSetAttribute(evaluateKernel<T, Dock, NT>, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(shm)),
          "cudaFuncSetAttribute");
  evaluateKernel<T, Dock, NT><<<n, NT, shm, st>>>(m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6, shm ? 1 : 0,
                                                  keepX ? 1 : 0);
}

template <class T, bool Dock>
void launch(cudaStream_t st, int threads, int n, std::size_t shm, bool keepX, const EnergyModelView& m, const T* angles,
            T* xs, double* gs, FlatTerms* terms, double* grad, const DockView& dv, const DockPose* poses,
            DockTerms* dterms, double* grad6) {
  switch (threads) {
    case 32: launchNT<T, Dock, 32>(st, n, shm, keepX, m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 64: launchNT<T, Dock, 64>(st, n, shm, keepX, m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 128: launchNT<T, Dock, 128>(st, n, shm, keepX, m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    default: launchNT<T, Dock, 256>(st, n, shm, keepX, m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
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
    hIn.release();
    hOut.release();
  }
};

}  // namespace

struct GpuEnergy::Impl {
  int device = 0;
  EnergyModelView view{};
  DevBuf<double> ref, charge, sv, vol;
  DevBuf<FlatTorsion> torsions;
  DevBuf<int> moving;
  DevBuf<FlatPairParams> pairParams;
  DevBuf<FlatPair> pairs;
  DevBuf<FlatClosure> closures;
  Lane lanes[GpuEnergy::kLanes];
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
    launch<T, false>(L.stream, blockThreads(view, 0), n, shm, hCoords != nullptr, view,
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
    launch<T, false>(L.stream, blockThreads(view, 0), n, shm, false, view,
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
    if (f.id == 0) throw Error("evaluateDock: no field");
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
#pragma omp parallel for schedule(static)
    for (long long s = 0; s < static_cast<long long>(nn); ++s) {
      std::copy_n(hAngles + s * nt, nt, ha + s * nt);
      hp[s] = hPoses[s];
    }
    check(cudaMemcpyAsync(L.in.p, L.hIn.p, inBytes, cudaMemcpyHostToDevice, L.stream), "upload");
    DockView dv = dock;
    dv.field = f.view;
    launch<float, true>(L.stream, blockThreads(view, dv.probeCount), n, shm, false, view,
                        reinterpret_cast<const float*>(L.in.p), reinterpret_cast<float*>(L.xs.p),
                        wantGrad ? L.gs.p : nullptr, reinterpret_cast<FlatTerms*>(L.out.p),
                        wantGrad ? reinterpret_cast<double*>(L.out.p + L.oGrad) : nullptr, dv,
                        reinterpret_cast<const DockPose*>(L.in.p + iPoses),
                        reinterpret_cast<DockTerms*>(L.out.p + L.oDterms),
                        wantGrad ? reinterpret_cast<double*>(L.out.p + L.oGrad6) : nullptr);
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
  d.charge.upload(m.charge);
  d.sv.upload(m.sv);
  d.vol.upload(m.vol);
  d.torsions.upload(m.torsions);
  d.moving.upload(m.moving);
  d.pairParams.upload(m.pairParams);
  d.pairs.upload(m.pairs);
  d.closures.upload(m.closures);
  d.view = m.view();
  d.view.ref = d.ref.p;
  d.view.charge = d.charge.p;
  d.view.sv = d.sv.p;
  d.view.vol = d.vol.p;
  d.view.torsions = d.torsions.p;
  d.view.moving = d.moving.p;
  d.view.pairParams = d.pairParams.p;
  d.view.pairs = d.pairs.p;
  d.view.closures = d.closures.p;
  for (Lane& L : d.lanes) check(cudaStreamCreateWithFlags(&L.stream, cudaStreamNonBlocking), "cudaStreamCreate");
}

GpuEnergy::~GpuEnergy() = default;

std::size_t GpuEnergy::bytesPerSample(bool fp64, bool grad) const {
  const std::size_t nt = static_cast<std::size_t>(impl_->view.torsionCount), nc = 3 * static_cast<std::size_t>(impl_->view.coordCount);
  const std::size_t r = fp64 ? sizeof(double) : sizeof(float);
  return (nt + nc) * r + sizeof(FlatTerms) + (grad ? (nc + nt) * sizeof(double) : 0);
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
