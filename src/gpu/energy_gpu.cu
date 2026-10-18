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

constexpr int kThreads = 256;
constexpr int kWarps = kThreads / 32;

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

template <class T>
__global__ void __launch_bounds__(kThreads) evaluateKernel(EnergyModelView m, const T* angles, T* xs, double* gs,
                                                           FlatTerms* terms, double* grad) {
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

  __shared__ double red[kWarps][8];
  double v[8] = {e[0], e[1], e[2], e[3], e[4], closure, maxDev, static_cast<double>(clashPairs)};
  for (int k = 0; k < 8; ++k) {
    v[k] = k == 6 ? warpMax(v[k]) : warpSum(v[k]);
    if (lane == 0) red[warp][k] = v[k];
  }
  __syncthreads();
  if (tid == 0) {
    double r[8];
    for (int k = 0; k < 8; ++k) {
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

struct GpuEnergy::Impl {
  int device = 0;
  EnergyModelView view{};
  DevBuf<double> ref, charge, sv, vol;
  DevBuf<FlatTorsion> torsions;
  DevBuf<int> moving;
  DevBuf<FlatPairParams> pairParams;
  DevBuf<FlatPair> pairs;
  DevBuf<FlatClosure> closures;
  DevBuf<char> angles, xs;
  DevBuf<double> gs, grad;
  DevBuf<FlatTerms> terms;

  template <class T>
  void run(const T* hAngles, int n, FlatTerms* hTerms, double* hGrad, T* hCoords) {
    if (n <= 0) return;
    check(cudaSetDevice(device), "cudaSetDevice");
    const std::size_t nt = static_cast<std::size_t>(view.torsionCount), nc = 3 * static_cast<std::size_t>(view.coordCount);
    const std::size_t nn = static_cast<std::size_t>(n);
    angles.reserve(std::max<std::size_t>(1, nn * nt) * sizeof(T));
    xs.reserve(std::max<std::size_t>(1, nn * nc) * sizeof(T));
    terms.reserve(nn);
    if (hGrad) {
      gs.reserve(std::max<std::size_t>(1, nn * nc));
      grad.reserve(std::max<std::size_t>(1, nn * nt));
    }
    if (nt) check(cudaMemcpy(angles.p, hAngles, nn * nt * sizeof(T), cudaMemcpyHostToDevice), "upload angles");
    evaluateKernel<T><<<n, kThreads>>>(view, reinterpret_cast<const T*>(angles.p), reinterpret_cast<T*>(xs.p),
                                       hGrad ? gs.p : nullptr, terms.p, hGrad ? grad.p : nullptr);
    check(cudaGetLastError(), "evaluateKernel launch");
    check(cudaDeviceSynchronize(), "evaluateKernel");
    check(cudaMemcpy(hTerms, terms.p, nn * sizeof(FlatTerms), cudaMemcpyDeviceToHost), "download terms");
    if (hGrad && nt) check(cudaMemcpy(hGrad, grad.p, nn * nt * sizeof(double), cudaMemcpyDeviceToHost), "download gradient");
    if (hCoords) check(cudaMemcpy(hCoords, xs.p, nn * nc * sizeof(T), cudaMemcpyDeviceToHost), "download coordinates");
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

void GpuEnergy::releaseScratch() {
  impl_->angles.release();
  impl_->xs.release();
  impl_->gs.release();
  impl_->grad.release();
  impl_->terms.release();
}

void GpuEnergy::evaluate(const double* angles, int n, FlatTerms* terms, double* grad, double* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

void GpuEnergy::evaluate(const float* angles, int n, FlatTerms* terms, double* grad, float* coords) {
  impl_->run(angles, n, terms, grad, coords);
}

}  // namespace od
