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
  T* upload(const T* v, std::size_t count) {
    reserve(count);
    if (count) check(cudaMemcpy(p, v, count * sizeof(T), cudaMemcpyHostToDevice), "upload");
    return p;
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

// All maps (type maps, elec, desolv, distance) of the planes [z0, z0 + nz); out: maps × count.
__global__ void sampleKernel(GridBuildView g, DistBuildView d, double clamp, int z0, std::size_t count, float* out) {
  const std::size_t s = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (s >= count) return;
  const std::size_t row = s / g.n[0];
  const int i = static_cast<int>(s - row * g.n[0]);
  const int j = static_cast<int>(row % g.n[1]);
  const int k = z0 + static_cast<int>(row / g.n[1]);
  double e[kGridMaxTypeMaps];
  double el, ds;
  gridPoint(g, i, j, k, e, el, ds);
  const int nt = g.typeMapCount;
  for (int t = 0; t < nt; ++t) out[t * count + s] = static_cast<float>(clamp < e[t] ? clamp : e[t]);
  out[nt * count + s] = static_cast<float>(el);
  out[(nt + 1) * count + s] = static_cast<float>(ds);
  out[(nt + 2) * count + s] = distPoint(d, i, j, k);
}

struct DeltaArgs {
  int n[3], nb[3], maps, bz0;  // grid points, blocks per axis, maps, first block layer of the slab
  int z0;                      // first plane of the slab
  std::size_t count, points;   // slab points, grid points
  const float* cur;            // maps × count
  const float* ref;            // maps × points
  double threshold;
};

// One thread block per grid block of the slab: flag[b] = 1 when a value differs.
__global__ void compareKernel(DeltaArgs a, int* flag) {
  const int b = blockIdx.x;
  const int bx = b % a.nb[0], by = (b / a.nb[0]) % a.nb[1], bz = a.bz0 + b / (a.nb[0] * a.nb[1]);
  const int t = threadIdx.x;
  const int x = bx * kB + t % kB, y = by * kB + (t / kB) % kB, z = bz * kB + t / (kB * kB);
  int differs = 0;
  if (x < a.n[0] && y < a.n[1] && z < a.n[2]) {
    const std::size_t sl = (static_cast<std::size_t>(z - a.z0) * a.n[1] + y) * a.n[0] + x;
    const std::size_t gl = (static_cast<std::size_t>(z) * a.n[1] + y) * a.n[0] + x;
    for (int m = 0; m < a.maps && !differs; ++m) {
      const double s = a.cur[m * a.count + sl], r = a.ref[m * a.points + gl];
      differs = fabs(s - r) > a.threshold;
    }
  }
  differs = __syncthreads_or(differs);
  if (t == 0) flag[b] = differs;
}

// Copies the blocks with slot[b] >= 0 to data + slot[b] × maps × kBlockPoints.
__global__ void gatherKernel(DeltaArgs a, const int* slot, float* data) {
  const int b = blockIdx.x;
  if (slot[b] < 0) return;
  const int bx = b % a.nb[0], by = (b / a.nb[0]) % a.nb[1], bz = a.bz0 + b / (a.nb[0] * a.nb[1]);
  const int t = threadIdx.x;
  const int x = bx * kB + t % kB, y = by * kB + (t / kB) % kB, z = bz * kB + t / (kB * kB);
  const bool in = x < a.n[0] && y < a.n[1] && z < a.n[2];
  const std::size_t sl = in ? (static_cast<std::size_t>(z - a.z0) * a.n[1] + y) * a.n[0] + x : 0;
  float* out = data + static_cast<std::size_t>(slot[b]) * a.maps * kBlockPoints;
  for (int m = 0; m < a.maps; ++m) out[m * kBlockPoints + t] = in ? a.cur[m * a.count + sl] : 0.0f;
}

// Uploads the arrays of a view; returns it with device pointers.
template <class Bufs>
GridBuildView uploadGrid(Bufs& d, const GridBuildView& h) {
  const std::size_t na = static_cast<std::size_t>(h.cells.atomCount);
  const std::size_t nt = static_cast<std::size_t>(h.typeMapCount);
  const std::size_t npp = nt * static_cast<std::size_t>(h.nTypes);
  const std::size_t cellCount = static_cast<std::size_t>(h.cells.n[0]) * h.cells.n[1] * h.cells.n[2];
  GridBuildView v = h;
  v.xyz = d.xyz.upload(h.xyz, 3 * na);
  v.type = d.type.upload(h.type, na);
  v.elecQ = d.elecQ.upload(h.elecQ, na);
  v.S = d.S.upload(h.S, na);
  v.V = d.V.upload(h.V, na);
  v.solT = d.solT.upload(h.solT, nt);
  v.volT = d.volT.upload(h.volT, nt);
  v.AC = d.AC.upload(h.AC, npp);
  v.BD = d.BD.upload(h.BD, npp);
  v.reqm = d.reqm.upload(h.reqm, npp);
  v.pm = d.pm.upload(h.pm, npp);
  v.pn = d.pn.upload(h.pn, npp);
  if (na) {
    v.cells.start = d.start.upload(h.cells.start, cellCount + 1);
    v.cells.atoms = d.atoms.upload(h.cells.atoms, na);
  }
  return v;
}

template <class Bufs>
DistBuildView uploadDist(Bufs& d, const DistBuildView& h) {
  const std::size_t na = static_cast<std::size_t>(h.cells.atomCount);
  const std::size_t cellCount = static_cast<std::size_t>(h.cells.n[0]) * h.cells.n[1] * h.cells.n[2];
  DistBuildView v = h;
  v.xyz = d.dXyz.upload(h.xyz, 3 * na);
  if (na) {
    v.cells.start = d.dStart.upload(h.cells.start, cellCount + 1);
    v.cells.atoms = d.dAtoms.upload(h.cells.atoms, na);
  }
  return v;
}

}  // namespace

struct GpuGridBuilder::Impl {
  int device = 0;
  std::size_t slabValues = kSlabValues;
  DevBuf<double> xyz, elecQ, S, V, solT, volT, AC, BD, reqm;
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
};

GpuGridBuilder::GpuGridBuilder(int device, std::size_t slabValues) : impl_(std::make_unique<Impl>()) {
  impl_->device = device;
  if (slabValues) impl_->slabValues = slabValues;
  check(cudaSetDevice(device), "cudaSetDevice");
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

void GpuGridBuilder::setReference(std::uint64_t id, const float* const* maps, int mapCount, std::size_t points) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  d.refId = 0;
  d.ref.reserve(static_cast<std::size_t>(mapCount) * points);
  for (int m = 0; m < mapCount; ++m)
    check(cudaMemcpy(d.ref.p + m * points, maps[m], points * sizeof(float), cudaMemcpyHostToDevice), "upload reference");
  d.refMaps = mapCount;
  d.refPoints = points;
  d.refId = id;
}

std::uint64_t GpuGridBuilder::referenceId() const { return impl_->refId; }

void GpuGridBuilder::buildDelta(const GridBuildView& gh, const DistBuildView& dh, double clamp, double threshold,
                                std::vector<int>& block, std::vector<float>& data) {
  Impl& d = *impl_;
  check(cudaSetDevice(d.device), "cudaSetDevice");
  const int maps = gh.typeMapCount + 3;
  const std::size_t points = static_cast<std::size_t>(gh.n[0]) * gh.n[1] * gh.n[2];
  if (d.refId == 0 || maps != d.refMaps || points != d.refPoints) throw Error("buildDelta: no matching reference");
  const GridBuildView g = uploadGrid(d, gh);
  const DistBuildView dv = uploadDist(d, dh);
  DeltaArgs a{};
  for (int c = 0; c < 3; ++c) {
    a.n[c] = gh.n[c];
    a.nb[c] = (gh.n[c] + kB - 1) / kB;
  }
  a.maps = maps;
  a.points = points;
  a.ref = d.ref.p;
  a.threshold = threshold;
  const std::size_t plane = static_cast<std::size_t>(gh.n[0]) * gh.n[1];
  const int layerBlocks = a.nb[0] * a.nb[1];
  // Block layers per slab.
  const int layers = static_cast<int>(std::clamp<std::size_t>(d.slabValues / (maps * plane * kB), 1, a.nb[2]));
  d.out.reserve(static_cast<std::size_t>(maps) * plane * kB * layers);
  d.flag.reserve(static_cast<std::size_t>(layerBlocks) * layers);
  d.slot.reserve(static_cast<std::size_t>(layerBlocks) * layers);
  block.assign(static_cast<std::size_t>(layerBlocks) * a.nb[2], -1);
  data.clear();
  int stored = 0;
  std::vector<int> flag, slot;
  for (int bz0 = 0; bz0 < a.nb[2]; bz0 += layers) {
    const int nl = std::min(layers, a.nb[2] - bz0);
    const int z0 = bz0 * kB, nz = std::min(gh.n[2], z0 + nl * kB) - z0;
    a.bz0 = bz0;
    a.z0 = z0;
    a.count = plane * nz;
    a.cur = d.out.p;
    sampleKernel<<<static_cast<unsigned>((a.count + kGridThreads - 1) / kGridThreads), kGridThreads>>>(g, dv, clamp, z0,
                                                                                                       a.count, d.out.p);
    check(cudaGetLastError(), "sampleKernel launch");
    const int nblk = layerBlocks * nl;
    compareKernel<<<nblk, kBlockPoints>>>(a, d.flag.p);
    check(cudaGetLastError(), "compareKernel launch");
    flag.resize(nblk);
    check(cudaMemcpy(flag.data(), d.flag.p, nblk * sizeof(int), cudaMemcpyDeviceToHost), "compareKernel");
    slot.assign(nblk, -1);
    int ns = 0;
    for (int b = 0; b < nblk; ++b)
      if (flag[b]) {
        slot[b] = ns++;
        block[static_cast<std::size_t>(bz0) * layerBlocks + b] = stored + slot[b];
      }
    if (ns == 0) continue;
    check(cudaMemcpy(d.slot.p, slot.data(), nblk * sizeof(int), cudaMemcpyHostToDevice), "upload slots");
    const std::size_t per = static_cast<std::size_t>(maps) * kBlockPoints;
    d.blocks.reserve(per * ns);
    gatherKernel<<<nblk, kBlockPoints>>>(a, d.slot.p, d.blocks.p);
    check(cudaGetLastError(), "gatherKernel launch");
    data.resize(data.size() + per * ns);
    check(cudaMemcpy(data.data() + per * stored, d.blocks.p, per * ns * sizeof(float), cudaMemcpyDeviceToHost),
          "gatherKernel");
    stored += ns;
  }
}

}  // namespace od
