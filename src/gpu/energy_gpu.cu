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

// One block of NT threads per sample.
template <class T, bool Dock, int NT>
__global__ void __launch_bounds__(NT) evaluateKernel(EnergyModelView m, const T* angles, T* xs, double* gs,
                                                           FlatTerms* terms, double* grad, DockView dv,
                                                           const DockPose* poses, DockTerms* dterms, double* grad6) {
  constexpr int kThreads = NT, kWarps = NT / 32;
  const int s = blockIdx.x;
  const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
  const std::size_t nc = 3 * static_cast<std::size_t>(m.coordCount);
  const T* ang = angles + static_cast<std::size_t>(s) * m.torsionCount;
  T* x = xs + s * nc;
  double* g = grad ? gs + s * nc : nullptr;

  for (std::size_t i = tid; i < nc; i += kThreads) {
    x[i] = static_cast<T>(m.ref[i]);
    if (g) g[i] = 0;
  }
  __syncthreads();

  __shared__ T so[3], sk[3], scs[2];
  __shared__ int sok;
  for (int t = m.torsionCount; t-- > 0;) {
    if (tid == 0) {
      T o[3], k[3], c = 0, sn = 0;
      sok = flat::torsionAxis(m, t, ang, x, o, k, c, sn);
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
  for (int p = tid; p < m.pairCount; p += kThreads) {
    T d[3], r, dedr;
    if (!flat::pairTerm(m, p, x, e, clashPairs, d, r, dedr)) continue;
    if (!g || r <= 0) continue;
    const T f = dedr / r;
    const FlatPair& pr = m.pairs[p];
    for (int c = 0; c < 3; ++c) {
      const double v = f * d[c];
      atomicAdd(g + 3 * pr.i + c, v);
      atomicAdd(g + 3 * pr.j + c, -v);
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

template <class T, bool Dock>
void launch(int threads, int n, const EnergyModelView& m, const T* angles, T* xs, double* gs, FlatTerms* terms,
            double* grad, const DockView& dv, const DockPose* poses, DockTerms* dterms, double* grad6) {
  switch (threads) {
    case 32: evaluateKernel<T, Dock, 32><<<n, 32>>>(m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 64: evaluateKernel<T, Dock, 64><<<n, 64>>>(m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    case 128: evaluateKernel<T, Dock, 128><<<n, 128>>>(m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
    default: evaluateKernel<T, Dock, 256><<<n, 256>>>(m, angles, xs, gs, terms, grad, dv, poses, dterms, grad6); break;
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
  // Work memory: inputs and outputs are each one contiguous block (one copy per call).
  DevBuf<char> in, out, xs;
  DevBuf<double> gs;
  PinnedBuf hIn, hOut;
  DevBuf<DockProbe> probes;
  DockView dock{};

  template <class T>
  void run(const T* hAngles, int n, FlatTerms* hTerms, double* hGrad, T* hCoords) {
    if (n <= 0) return;
    check(cudaSetDevice(device), "cudaSetDevice");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t inBytes = align8(nn * nt * sizeof(T));
    const std::size_t oTerms = 0, oGrad = align8(nn * sizeof(FlatTerms));
    const std::size_t outBytes = oGrad + (hGrad ? nn * nt * sizeof(double) : 0);
    in.reserve(std::max<std::size_t>(8, inBytes));
    out.reserve(outBytes);
    xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(T));
    if (hGrad) gs.reserve(std::max<std::size_t>(1, nn * nc));
    hIn.reserve(std::max<std::size_t>(8, inBytes));
    hOut.reserve(outBytes);
    if (nt) {
      std::copy_n(hAngles, nn * nt, reinterpret_cast<T*>(hIn.p));
      check(cudaMemcpy(in.p, hIn.p, inBytes, cudaMemcpyHostToDevice), "upload");
    }
    launch<T, false>(blockThreads(view, 0), n, view, reinterpret_cast<const T*>(in.p), reinterpret_cast<T*>(xs.p),
                     hGrad ? gs.p : nullptr, reinterpret_cast<FlatTerms*>(out.p + oTerms),
                     hGrad ? reinterpret_cast<double*>(out.p + oGrad) : nullptr, DockView{}, nullptr, nullptr, nullptr);
    check(cudaGetLastError(), "evaluateKernel launch");
    check(cudaMemcpy(hOut.p, out.p, outBytes, cudaMemcpyDeviceToHost), "download");
    std::copy_n(reinterpret_cast<const FlatTerms*>(hOut.p + oTerms), nn, hTerms);
    if (hGrad) std::copy_n(reinterpret_cast<const double*>(hOut.p + oGrad), nn * nt, hGrad);
    if (hCoords) check(cudaMemcpy(hCoords, xs.p, nn * nc * sizeof(T), cudaMemcpyDeviceToHost), "download coordinates");
  }

  void runDock(const GpuField::Impl& f, const float* hAngles, const DockPose* hPoses, int n, FlatTerms* hTerms,
               DockTerms* hDterms, double* hGrad) {
    if (n <= 0) return;
    if (f.id == 0) throw Error("evaluateDock: no field");
    check(cudaSetDevice(device), "cudaSetDevice");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    const std::size_t iPoses = align8(nn * nt * sizeof(float)), inBytes = iPoses + nn * sizeof(DockPose);
    const std::size_t oTerms = 0, oDterms = align8(nn * sizeof(FlatTerms)), oGrad6 = oDterms + align8(nn * sizeof(DockTerms)),
                      oGrad = oGrad6 + nn * 6 * sizeof(double);
    const std::size_t outBytes = hGrad ? oGrad + nn * nt * sizeof(double) : oGrad6;
    in.reserve(inBytes);
    out.reserve(outBytes);
    xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(float));
    if (hGrad) gs.reserve(std::max<std::size_t>(1, nn * nc));
    hIn.reserve(inBytes);
    hOut.reserve(outBytes);
    std::copy_n(hAngles, nn * nt, reinterpret_cast<float*>(hIn.p));
    std::copy_n(hPoses, nn, reinterpret_cast<DockPose*>(hIn.p + iPoses));
    check(cudaMemcpy(in.p, hIn.p, inBytes, cudaMemcpyHostToDevice), "upload");
    DockView dv = dock;
    dv.field = f.view;
    launch<float, true>(blockThreads(view, dv.probeCount), n, view, reinterpret_cast<const float*>(in.p),
                        reinterpret_cast<float*>(xs.p), hGrad ? gs.p : nullptr, reinterpret_cast<FlatTerms*>(out.p + oTerms),
                        hGrad ? reinterpret_cast<double*>(out.p + oGrad) : nullptr, dv,
                        reinterpret_cast<const DockPose*>(in.p + iPoses), reinterpret_cast<DockTerms*>(out.p + oDterms),
                        hGrad ? reinterpret_cast<double*>(out.p + oGrad6) : nullptr);
    check(cudaGetLastError(), "evaluateKernel (dock) launch");
    check(cudaMemcpy(hOut.p, out.p, outBytes, cudaMemcpyDeviceToHost), "download");
    std::copy_n(reinterpret_cast<const FlatTerms*>(hOut.p + oTerms), nn, hTerms);
    std::copy_n(reinterpret_cast<const DockTerms*>(hOut.p + oDterms), nn, hDterms);
    if (!hGrad) return;
    const double* g6 = reinterpret_cast<const double*>(hOut.p + oGrad6);
    const double* gt = reinterpret_cast<const double*>(hOut.p + oGrad);
    const std::size_t ng = 6 + nt;
    for (std::size_t s = 0; s < nn; ++s) {
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
  impl_->runDock(*field.impl_, angles, poses, n, terms, dterms, grad);
}

std::size_t GpuEnergy::bytesPerSampleDock() const {
  return bytesPerSample(false, true) + sizeof(DockPose) + sizeof(DockTerms) + 6 * sizeof(double);
}

void GpuEnergy::releaseScratch() {
  impl_->in.release();
  impl_->out.release();
  impl_->xs.release();
  impl_->gs.release();
  impl_->hIn.release();
  impl_->hOut.release();
}

void GpuEnergy::evaluate(const double* angles, int n, FlatTerms* terms, double* grad, double* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

void GpuEnergy::evaluate(const float* angles, int n, FlatTerms* terms, double* grad, float* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

}  // namespace od
