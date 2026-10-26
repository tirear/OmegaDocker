// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "grid_gpu.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <string>
#include <vector>

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
  void reserve(std::size_t count) {
    if (count <= n) return;
    if (p) cudaFree(p);
    p = nullptr;
    n = 0;
    check(cudaMalloc(&p, std::max<std::size_t>(count, 1) * sizeof(T)), "cudaMalloc");
    n = count;
  }
  // Pageable v: copied to staging memory before the call returns (v may then be released).
  T* upload(const T* v, std::size_t count, cudaStream_t s = 0) {
    reserve(count);
    if (count) check(cudaMemcpyAsync(p, v, count * sizeof(T), cudaMemcpyHostToDevice, s), "upload");
    if (!s && count) check(cudaStreamSynchronize(0), "upload");
    return p;
  }
};

// Pinned host buffer.
template <class T>
struct HostBuf {
  T* p = nullptr;
  std::size_t n = 0;
  HostBuf() = default;
  HostBuf(const HostBuf&) = delete;
  HostBuf& operator=(const HostBuf&) = delete;
  ~HostBuf() {
    if (p) cudaFreeHost(p);
  }
  void reserve(std::size_t count) {
    if (count <= n) return;
    if (p) cudaFreeHost(p);
    p = nullptr;
    n = 0;
    check(cudaMallocHost(&p, std::max<std::size_t>(count, 1) * sizeof(T)), "cudaMallocHost");
    n = count;
  }
};

constexpr int kGridThreads = 128;
constexpr std::size_t kSlabValues = std::size_t(1) << 26;  // output floats per launch

// Points [p0, p0 + count) of the grid; out: (typeMapCount + 2) × count values.
__global__ void gridKernel(GridBuildView v, double clamp, std::size_t p0, std::size_t count, float* out) {
  const std::size_t s = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (s >= count) return;
  const std::size_t idx = p0 + s;
  const std::size_t row = idx / v.n[0];
  const int i = static_cast<int>(idx - row * v.n[0]);
  const int j = static_cast<int>(row % v.n[1]);
  const int k = static_cast<int>(row / v.n[1]);
  double e[kGridMaxTypeMaps];
  double el, ds;
  gridPoint(v, i, j, k, e, el, ds);
  const int nt = v.typeMapCount;
  for (int t = 0; t < nt; ++t) out[t * count + s] = static_cast<float>(clamp < e[t] ? clamp : e[t]);  // std::min(e, clamp)
  out[nt * count + s] = static_cast<float>(el);
  out[(nt + 1) * count + s] = static_cast<float>(ds);
}

__global__ void distKernel(DistBuildView v, std::size_t p0, std::size_t count, float* out) {
  const std::size_t s = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (s >= count) return;
  const std::size_t idx = p0 + s;
  const std::size_t row = idx / v.n[0];
  out[s] = distPoint(v, static_cast<int>(idx - row * v.n[0]), static_cast<int>(row % v.n[1]),
                     static_cast<int>(row / v.n[1]));
}

constexpr int kB = kDeltaBlock;
constexpr int kBlockPoints = kB * kB * kB;

struct DeltaArgs {
  int n[3], nb[3], maps, bz0;  // grid points, blocks per axis, maps, first block layer of the slab
  std::size_t points;          // grid points
  const float* ref;            // maps × points
  double threshold;
  double origin[3], spacing, reach2;
  const double* changed;  // changedCount positions; changedCount < 0: every block
  int changedCount;
};

// One thread block per grid block of the slab: all maps (type maps, elec, desolv, distance) of
// its points into scratch (per block maps × kBlockPoints values); flag[b] = 1 when a value
// differs from the reference. Blocks no changed position reaches: flag 0, nothing computed.
// MaxT >= type maps (e in registers).
template <int MaxT>
__global__ void __launch_bounds__(kBlockPoints, 2)
    deltaKernel(GridBuildView g, DistBuildView d, double clamp, DeltaArgs a, float* scratch, int* flag) {
  const int b = blockIdx.x;
  const int bx = b % a.nb[0], by = (b / a.nb[0]) % a.nb[1], bz = a.bz0 + b / (a.nb[0] * a.nb[1]);
  const int t = threadIdx.x;
  if (a.changedCount >= 0) {
    int near = 0;
    for (int c = t; c < a.changedCount && !near; c += blockDim.x)
      near = deltaBlockReached(a.origin, a.spacing, bx, by, bz, a.changed + 3 * c, a.reach2);
    if (!__syncthreads_or(near)) {
      if (t == 0) flag[b] = 0;
      return;
    }
  }
  const int x = bx * kB + t % kB, y = by * kB + (t / kB) % kB, z = bz * kB + t / (kB * kB);
  float* out = scratch + static_cast<std::size_t>(b) * a.maps * kBlockPoints;
  const int nt = g.typeMapCount;
  int differs = 0;
  if (x < a.n[0] && y < a.n[1] && z < a.n[2]) {
    double e[MaxT];
    double el, ds;
    gridPoint<MaxT>(g, x, y, z, e, el, ds);
    const std::size_t gl = (static_cast<std::size_t>(z) * a.n[1] + y) * a.n[0] + x;
    const auto put = [&](int m, float v) {
      out[m * kBlockPoints + t] = v;
      const double s = v, r = a.ref[m * a.points + gl];
      differs |= fabs(s - r) > a.threshold;
    };
#pragma unroll
    for (int m = 0; m < MaxT; ++m)
      if (m < nt) put(m, static_cast<float>(clamp < e[m] ? clamp : e[m]));
    put(nt, static_cast<float>(el));
    put(nt + 1, static_cast<float>(ds));
    put(nt + 2, distPoint(d, x, y, z));
  } else {
    for (int m = 0; m < a.maps; ++m) out[m * kBlockPoints + t] = 0.0f;
  }
  differs = __syncthreads_or(differs);
  if (t == 0) flag[b] = differs;
}

// Copies the blocks with slot[b] >= 0 to data + slot[b] × maps × kBlockPoints.
__global__ void gatherKernel(int maps, const float* scratch, const int* slot, float* data) {
  const int b = blockIdx.x;
  if (slot[b] < 0) return;
  const std::size_t per = static_cast<std::size_t>(maps) * kBlockPoints;
  const float* in = scratch + b * per;
  float* out = data + static_cast<std::size_t>(slot[b]) * per;
  for (std::size_t v = threadIdx.x; v < per; v += blockDim.x) out[v] = in[v];
}

// Delta blocks owned in device memory.
struct OwnedBlocks : DeviceBlocks {
  ~OwnedBlocks() override {
    if (ptr) cudaFree(const_cast<float*>(ptr));
  }
};

// Uploads the arrays of a view; returns it with device pointers.
template <class Bufs>
GridBuildView uploadGrid(Bufs& d, const GridBuildView& h, cudaStream_t s = 0) {
  const std::size_t na = static_cast<std::size_t>(h.cells.atomCount);
  const std::size_t nt = static_cast<std::size_t>(h.typeMapCount);
  const std::size_t npp = nt * static_cast<std::size_t>(h.nTypes);
  const std::size_t cellCount = static_cast<std::size_t>(h.cells.n[0]) * h.cells.n[1] * h.cells.n[2];
  GridBuildView v = h;
  v.xyz = d.xyz.upload(h.xyz, 3 * na, s);
  v.type = d.type.upload(h.type, na, s);
  v.elecQ = d.elecQ.upload(h.elecQ, na, s);
  v.S = d.S.upload(h.S, na, s);
  v.V = d.V.upload(h.V, na, s);
  v.solT = d.solT.upload(h.solT, nt, s);
  v.volT = d.volT.upload(h.volT, nt, s);
  v.desolvW = d.desolvW.upload(h.desolvW, nt * na, s);
  v.dsW = d.dsW.upload(h.dsW, na, s);
  v.AC = d.AC.upload(h.AC, npp, s);
  v.BD = d.BD.upload(h.BD, npp, s);
  v.reqm = d.reqm.upload(h.reqm, npp, s);
  v.pm = d.pm.upload(h.pm, npp, s);
  v.pn = d.pn.upload(h.pn, npp, s);
  if (na) {
    v.cells.start = d.start.upload(h.cells.start, cellCount + 1, s);
    v.cells.atoms = d.atoms.upload(h.cells.atoms, na, s);
  }
  return v;
}

template <class Bufs>
DistBuildView uploadDist(Bufs& d, const DistBuildView& h, cudaStream_t s = 0) {
  const std::size_t na = static_cast<std::size_t>(h.cells.atomCount);
  const std::size_t cellCount = static_cast<std::size_t>(h.cells.n[0]) * h.cells.n[1] * h.cells.n[2];
  DistBuildView v = h;
  v.xyz = d.dXyz.upload(h.xyz, 3 * na, s);
  if (na) {
    v.cells.start = d.dStart.upload(h.cells.start, cellCount + 1, s);
    v.cells.atoms = d.dAtoms.upload(h.cells.atoms, na, s);
  }
  return v;
}

}  // namespace

struct GpuGridBuilder::Impl {
  int device = 0;
  std::size_t slabValues = kSlabValues;
  DevBuf<double> xyz, elecQ, S, V, solT, volT, desolvW, dsW, AC, BD, reqm;
  DevBuf<int> type, pm, pn, start, atoms;
  DevBuf<double> dXyz;
  DevBuf<int> dStart, dAtoms;
  DevBuf<float> out;
  // delta
  std::uint64_t refId = 0;
  int refMaps = 0;
  std::size_t refPoints = 0;
  DevBuf<float> ref, blocks;
  DevBuf<int> flag, slot;
  DevBuf<double> changed;
  HostBuf<int> hFlag;
  cudaStream_t stream = nullptr;  // delta (non-blocking: overlaps the work of the other streams)
  // The submitted delta (submitDelta): views of the uploaded inputs and the launch arguments.
  bool pending = false;
  GridBuildView g{};
  DistBuildView dv{};
  DeltaArgs a{};
  double clamp = 0;
  int layers = 0;
  ~Impl() {
    if (stream) {
      cudaStreamSynchronize(stream);
      cudaStreamDestroy(stream);
    }
  }
  // Delta kernel of the block layers [bz0, bz0 + nl) and the copy of their flags into hFlag.
  void launch(int bz0, int nl) {
    a.bz0 = bz0;
    const int nblk = a.nb[0] * a.nb[1] * nl;
    if (g.typeMapCount <= 4) deltaKernel<4><<<nblk, kBlockPoints, 0, stream>>>(g, dv, clamp, a, out.p, flag.p);
    else if (g.typeMapCount <= 8) deltaKernel<8><<<nblk, kBlockPoints, 0, stream>>>(g, dv, clamp, a, out.p, flag.p);
    else if (g.typeMapCount <= 16) deltaKernel<16><<<nblk, kBlockPoints, 0, stream>>>(g, dv, clamp, a, out.p, flag.p);
    else deltaKernel<kGridMaxTypeMaps><<<nblk, kBlockPoints, 0, stream>>>(g, dv, clamp, a, out.p, flag.p);
    check(cudaGetLastError(), "deltaKernel launch");
    check(cudaMemcpyAsync(hFlag.p, flag.p, nblk * sizeof(int), cudaMemcpyDeviceToHost, stream), "deltaKernel");
  }
};

GpuGridBuilder::GpuGridBuilder(int device, std::size_t slabValues) : impl_(std::make_unique<Impl>()) {
  impl_->device = device;
  if (slabValues) impl_->slabValues = slabValues;
  check(cudaSetDevice(device), "cudaSetDevice");
  int least = 0, greatest = 0;
  check(cudaDeviceGetStreamPriorityRange(&least, &greatest), "cudaDeviceGetStreamPriorityRange");
  check(cudaStreamCreateWithPriority(&impl_->stream, cudaStreamNonBlocking, least), "cudaStreamCreate");
}

GpuGridBuilder::~GpuGridBuilder() = default;

void GpuGridBuilder::build(const GridBuildView& h, double clamp, float* const* typeMaps, float* elec, float* desolv) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const std::size_t nt = static_cast<std::size_t>(h.typeMapCount);
  const GridBuildView v = uploadGrid(d, h);
  const std::size_t total = static_cast<std::size_t>(h.n[0]) * h.n[1] * h.n[2];
  const std::size_t maps = nt + 2;
  const std::size_t slab = std::max<std::size_t>(1, std::min(total, d.slabValues / maps));
  d.out.reserve(slab * maps);
  for (std::size_t p0 = 0; p0 < total; p0 += slab) {
    const std::size_t count = std::min(slab, total - p0);
    const unsigned blocks = static_cast<unsigned>((count + kGridThreads - 1) / kGridThreads);
    gridKernel<<<blocks, kGridThreads>>>(v, clamp, p0, count, d.out.p);
    check(cudaGetLastError(), "gridKernel launch");
    for (std::size_t m = 0; m < maps; ++m) {
      float* dst = m < nt ? typeMaps[m] : m == nt ? elec : desolv;
      check(cudaMemcpy(dst + p0, d.out.p + m * count, count * sizeof(float), cudaMemcpyDeviceToHost), "gridKernel");
    }
  }
}

void GpuGridBuilder::buildDistance(const DistBuildView& h, float* out) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const DistBuildView v = uploadDist(d, h);
  const std::size_t total = static_cast<std::size_t>(h.n[0]) * h.n[1] * h.n[2];
  const std::size_t slab = std::max<std::size_t>(1, std::min(total, d.slabValues));
  d.out.reserve(slab);
  for (std::size_t p0 = 0; p0 < total; p0 += slab) {
    const std::size_t count = std::min(slab, total - p0);
    distKernel<<<static_cast<unsigned>((count + kGridThreads - 1) / kGridThreads), kGridThreads>>>(v, p0, count, d.out.p);
    check(cudaGetLastError(), "distKernel launch");
    check(cudaMemcpy(out + p0, d.out.p, count * sizeof(float), cudaMemcpyDeviceToHost), "distKernel");
  }
}

std::shared_ptr<const DeviceBlocks> GpuGridBuilder::buildDevice(const GridBuildView& gh, const DistBuildView& dh,
                                                                double clamp) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const std::size_t total = static_cast<std::size_t>(gh.n[0]) * gh.n[1] * gh.n[2];
  const std::size_t maps = static_cast<std::size_t>(gh.typeMapCount) + 3;
  auto out = std::make_shared<OwnedBlocks>();
  out->device = d.device;
  out->floats = maps * total;
  float* p = nullptr;
  check(cudaMalloc(&p, out->floats * sizeof(float)), "buildDevice: cudaMalloc");
  out->ptr = p;
  const unsigned blocks = static_cast<unsigned>((total + kGridThreads - 1) / kGridThreads);
  // gridKernel writes map m of point p0 + i at m × count + i: with p0 = 0 and count = total, the
  // layout of FieldView::maps.
  const GridBuildView v = uploadGrid(d, gh);
  gridKernel<<<blocks, kGridThreads>>>(v, clamp, 0, total, p);
  check(cudaGetLastError(), "gridKernel launch");
  const DistBuildView dv = uploadDist(d, dh);
  distKernel<<<blocks, kGridThreads>>>(dv, 0, total, p + (maps - 1) * total);
  check(cudaGetLastError(), "distKernel launch");
  check(cudaDeviceSynchronize(), "buildDevice");
  return out;
}

void GpuGridBuilder::setReference(std::uint64_t id, const float* const* maps, int mapCount, std::size_t points) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  d.refId = 0;
  d.ref.reserve(static_cast<std::size_t>(mapCount) * points);
  for (int m = 0; m < mapCount; ++m)
    check(cudaMemcpy(d.ref.p + m * points, maps[m], points * sizeof(float), cudaMemcpyHostToDevice), "upload reference");
  check(cudaDeviceSynchronize(), "upload reference");
  d.refMaps = mapCount;
  d.refPoints = points;
  d.refId = id;
}

std::uint64_t GpuGridBuilder::referenceId() const { return impl_->refId; }

void GpuGridBuilder::buildDelta(const GridBuildView& gh, const DistBuildView& dh, double clamp, double threshold,
                                std::vector<int>& block, std::vector<float>& data) {
  buildDelta(gh, dh, clamp, threshold, nullptr, -1, block, &data);
}

std::shared_ptr<const DeviceBlocks> GpuGridBuilder::buildDelta(const GridBuildView& gh, const DistBuildView& dh,
                                                               double clamp, double threshold, const double* changed,
                                                               int changedCount, std::vector<int>& block,
                                                               std::vector<float>* data) {
  submitDelta(gh, dh, clamp, threshold, changed, changedCount);
  return finishDelta(block, data);
}

void GpuGridBuilder::submitDelta(const GridBuildView& gh, const DistBuildView& dh, double clamp, double threshold,
                                 const double* changed, int changedCount) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  if (d.pending) throw Error("submitDelta: a delta is pending");
  const int maps = gh.typeMapCount + 3;
  const std::size_t points = static_cast<std::size_t>(gh.n[0]) * gh.n[1] * gh.n[2];
  if (d.refId == 0 || maps != d.refMaps || points != d.refPoints) throw Error("buildDelta: no matching reference");
  d.g = uploadGrid(d, gh, d.stream);
  d.dv = uploadDist(d, dh, d.stream);
  DeltaArgs& a = d.a;
  a = DeltaArgs{};
  for (int c = 0; c < 3; ++c) {
    a.n[c] = gh.n[c];
    a.nb[c] = (gh.n[c] + kB - 1) / kB;
    a.origin[c] = gh.origin[c];
  }
  a.maps = maps;
  a.points = points;
  a.ref = d.ref.p;
  a.threshold = threshold;
  a.spacing = gh.spacing;
  a.reach2 = deltaReach2(gh.spacing, std::max(kElecCutoff, dh.cap));
  a.changedCount = changed ? changedCount : -1;
  a.changed = changed && changedCount > 0
                  ? d.changed.upload(changed, 3 * static_cast<std::size_t>(changedCount), d.stream)
                  : nullptr;
  d.clamp = clamp;
  const std::size_t per = static_cast<std::size_t>(maps) * kBlockPoints;
  const int layerBlocks = a.nb[0] * a.nb[1];
  // Block layers per slab.
  d.layers = static_cast<int>(std::clamp<std::size_t>(d.slabValues / (per * layerBlocks), 1, a.nb[2]));
  d.out.reserve(per * layerBlocks * d.layers);
  d.flag.reserve(static_cast<std::size_t>(layerBlocks) * d.layers);
  d.slot.reserve(static_cast<std::size_t>(layerBlocks) * d.layers);
  d.hFlag.reserve(static_cast<std::size_t>(layerBlocks) * d.layers);
  d.launch(0, std::min(d.layers, a.nb[2]));
  d.pending = true;
}

std::shared_ptr<const DeviceBlocks> GpuGridBuilder::finishDelta(std::vector<int>& block, std::vector<float>* data) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  if (!d.pending) throw Error("finishDelta: no delta submitted");
  d.pending = false;
  const DeltaArgs& a = d.a;
  const int maps = a.maps, layers = d.layers;
  const std::size_t per = static_cast<std::size_t>(maps) * kBlockPoints;
  const int layerBlocks = a.nb[0] * a.nb[1];
  block.assign(static_cast<std::size_t>(layerBlocks) * a.nb[2], -1);
  int stored = 0;
  std::vector<int> slot;
  // Blocks of the slabs: the stored blocks of one slab go straight to the result when it is the
  // only slab, otherwise to d.blocks (grown keeping the earlier slabs).
  auto result = std::make_shared<OwnedBlocks>();
  result->device = d.device;
  for (int bz0 = 0; bz0 < a.nb[2]; bz0 += layers) {
    const int nl = std::min(layers, a.nb[2] - bz0);
    if (bz0 > 0) d.launch(bz0, nl);
    const int nblk = layerBlocks * nl;
    check(cudaStreamSynchronize(d.stream), "deltaKernel");
    slot.assign(nblk, -1);
    int ns = 0;
    for (int b = 0; b < nblk; ++b)
      if (d.hFlag.p[b]) {
        slot[b] = ns++;
        block[static_cast<std::size_t>(bz0) * layerBlocks + b] = stored + slot[b];
      }
    if (ns == 0) continue;
    d.slot.upload(slot.data(), nblk, d.stream);
    float* dst;
    if (layers >= a.nb[2]) {
      float* p = nullptr;
      check(cudaMalloc(&p, per * ns * sizeof(float)), "cudaMalloc");
      result->ptr = p;
      result->floats = per * ns;
      dst = p;
    } else {
      if (d.blocks.n < per * (stored + ns)) {
        DevBuf<float> grown;
        grown.reserve(std::max(per * (stored + ns), 2 * d.blocks.n));
        if (stored)
          check(cudaMemcpyAsync(grown.p, d.blocks.p, per * stored * sizeof(float), cudaMemcpyDeviceToDevice, d.stream),
                "grow blocks");
        check(cudaStreamSynchronize(d.stream), "grow blocks");
        std::swap(grown.p, d.blocks.p);
        std::swap(grown.n, d.blocks.n);
      }
      dst = d.blocks.p + per * stored;
    }
    gatherKernel<<<nblk, kBlockPoints, 0, d.stream>>>(maps, d.out.p, d.slot.p, dst);
    check(cudaGetLastError(), "gatherKernel launch");
    stored += ns;
  }
  if (layers < a.nb[2] && stored) {
    float* p = nullptr;
    check(cudaMalloc(&p, per * stored * sizeof(float)), "cudaMalloc");
    result->ptr = p;
    result->floats = per * stored;
    check(cudaMemcpyAsync(p, d.blocks.p, per * stored * sizeof(float), cudaMemcpyDeviceToDevice, d.stream), "blocks");
  }
  if (data) {
    data->resize(per * stored);
    if (stored)
      check(cudaMemcpyAsync(data->data(), result->ptr, per * stored * sizeof(float), cudaMemcpyDeviceToHost, d.stream),
            "blocks");
  }
  check(cudaStreamSynchronize(d.stream), "buildDelta");
  return result;
}

}  // namespace od
