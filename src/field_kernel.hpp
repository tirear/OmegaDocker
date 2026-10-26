// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "ad4_pair.hpp"

namespace od {

// Flat form of a StamField (spec 8.1, 8.2) shared by the CPU and GPU back ends
// (spec 15.1). The functions below follow StamField::energy and StamField::distance
// operation by operation.

constexpr int kFieldMaxMaps = 40;  // type maps + elec + desolv + distance

struct FieldMoved {
  double xyz[3];
  int type;
  double elecQ, S, V;
  int heavy;
};

struct FieldCells {
  double cell, lo[3];
  int n[3];
  const int* start;
  const int* atoms;
  int atomCount;
};

struct FieldView {
  int mode;  // GridMode: 0 full, 1 flexres, 2 delta
  double origin[3], spacing;
  int n[3];
  int typeMapCount;
  const float* maps[kFieldMaxMaps];  // type maps, elec, desolv, distance (reference maps for delta)
  // GPU only (optional): the maps with the values of (i, j, k), (i, j + 1, k), (i, j, k + 1) and
  // (i, j + 1, k + 1) at each point (4 floats; map m from m × 4 × points).
  const float* packed;
  // delta
  int nb[3], mapCount;
  const int* block;
  const float* data;
  // flexres
  int movedCount, movedHeavyCount, nTypes;
  const FieldMoved* moved;
  const int* movedHeavy;
  FieldCells cells, heavyCells;
  const double *solT, *volT, *pairAC, *pairBD, *pairReqm;
  const int *pairM, *pairN;
  double cDesolv, qasp, halfSmooth, rMin, cap;
};

// Delta blocks in device memory (spec 15.1.19): blocks × maps × kDeltaBlock^3 floats in the
// layout of FieldModel::data; freed by the implementation's destructor.
struct DeviceBlocks {
  int device = 0;
  const float* ptr = nullptr;
  std::size_t floats = 0;
  std::vector<std::shared_ptr<const DeviceBlocks>> others;  // the same blocks on other devices (spec 15.1.33)
  virtual ~DeviceBlocks() = default;
};

// Host arrays of a FieldView. The maps are not copied: base keeps the grids that hold them.
struct FieldModel {
  std::uint64_t id = 0;               // StamField::id()
  std::shared_ptr<const void> base;   // owner of the reference maps (grids)
  std::shared_ptr<const void> baseDist;
  // full (optional): the maps in device memory, in the layout of FieldView::maps (spec 15.1.34);
  // view.maps are then null.
  std::shared_ptr<const DeviceBlocks> deviceMaps;
  std::size_t points = 0;             // values per map
  std::vector<int> block;
  std::shared_ptr<const std::vector<float>> data;  // delta blocks; null: only in device memory
  std::shared_ptr<const DeviceBlocks> deviceData;  // delta blocks in device memory (optional)
  std::vector<FieldMoved> moved;
  std::vector<int> movedHeavy, cellStart, cellAtoms, heavyStart, heavyAtoms;
  std::vector<double> solT, volT, pairAC, pairBD, pairReqm;
  std::vector<int> pairM, pairN;
  FieldView view{};  // arrays point into this model and the grids

  // Points the array members of view at the vectors above (after copy or move).
  void bind() {
    view.block = block.data();
    view.data = data ? data->data() : nullptr;
    view.moved = moved.data();
    view.movedHeavy = movedHeavy.data();
    view.cells.start = cellStart.data();
    view.cells.atoms = cellAtoms.data();
    view.heavyCells.start = heavyStart.data();
    view.heavyCells.atoms = heavyAtoms.data();
    view.solT = solT.data();
    view.volT = volT.data();
    view.pairAC = pairAC.data();
    view.pairBD = pairBD.data();
    view.pairReqm = pairReqm.data();
    view.pairM = pairM.data();
    view.pairN = pairN.data();
  }
};

namespace field {

OD_HD inline std::size_t pointIndex(const FieldView& f, int i, int j, int k) {
  return (static_cast<std::size_t>(k) * f.n[1] + j) * f.n[0] + i;
}

OD_HD inline float at(const FieldView& f, int m, int i, int j, int k) {
  if (f.mode != 2) return f.maps[m][pointIndex(f, i, j, k)];
  constexpr int B = 8;  // kDeltaBlock
  constexpr std::size_t bs = static_cast<std::size_t>(B) * B * B;
  const int b = f.block[(static_cast<std::size_t>(k / B) * f.nb[1] + j / B) * f.nb[0] + i / B];
  if (b < 0) return f.maps[m][pointIndex(f, i, j, k)];
  return f.data[(static_cast<std::size_t>(b) * f.mapCount + m) * bs + ((k % B) * B + j % B) * B + i % B];
}

// Trilinear interpolation of map m at p; grad (optional): dvalue/dp. Divisions by the spacing are
// multiplications by its reciprocal.
template <class T>
OD_HD inline T interpolate(const FieldView& f, int m, const T p[3], T* grad) {
  using std::floor;
  const T isp = T(1) / static_cast<T>(f.spacing);
  int i0[3];
  T fr[3];
  for (int c = 0; c < 3; ++c) {
    const T u = (p[c] - static_cast<T>(f.origin[c])) * isp;
    int i = static_cast<int>(floor(u));
    i = i < 0 ? 0 : i > f.n[c] - 2 ? f.n[c] - 2 : i;
    i0[c] = i;
    fr[c] = u - i;
  }
  T v[2][2][2];
  for (int z = 0; z < 2; ++z)
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 2; ++x) v[z][y][x] = at(f, m, i0[0] + x, i0[1] + y, i0[2] + z);
  const T fx = fr[0], fy = fr[1], fz = fr[2];
  auto lerp = [](T a, T b, T t) { return a + (b - a) * t; };
  const T c00 = lerp(v[0][0][0], v[0][0][1], fx), c10 = lerp(v[0][1][0], v[0][1][1], fx);
  const T c01 = lerp(v[1][0][0], v[1][0][1], fx), c11 = lerp(v[1][1][0], v[1][1][1], fx);
  const T c0 = lerp(c00, c10, fy), c1 = lerp(c01, c11, fy);
  if (grad) {
    const T dx0 = lerp(v[0][0][1] - v[0][0][0], v[0][1][1] - v[0][1][0], fy);
    const T dx1 = lerp(v[1][0][1] - v[1][0][0], v[1][1][1] - v[1][1][0], fy);
    grad[0] = lerp(dx0, dx1, fz) * isp;
    grad[1] = lerp(c10 - c00, c11 - c01, fz) * isp;
    grad[2] = (c1 - c0) * isp;
  }
  return lerp(c0, c1, fz);
}

// interpolate for the maps maps[0..count) at the same point: the cell and the corner addresses are
// found once (the values and gradients are those of interpolate). grad (optional): count × 3.
template <class T>
OD_HD inline void interpolateMaps(const FieldView& f, const int* maps, int count, const T p[3], T* val, T* grad) {
  using std::floor;
  const T isp = T(1) / static_cast<T>(f.spacing);
  int i0[3];
  T fr[3];
  for (int c = 0; c < 3; ++c) {
    const T u = (p[c] - static_cast<T>(f.origin[c])) * isp;
    int i = static_cast<int>(floor(u));
    i = i < 0 ? 0 : i > f.n[c] - 2 ? f.n[c] - 2 : i;
    i0[c] = i;
    fr[c] = u - i;
  }
  // Corner c = 4z + 2y + x: point index, or with delta blocks the block (−1: reference map) and the
  // offset in the block.
  constexpr int B = 8;  // kDeltaBlock
  constexpr std::size_t bs = static_cast<std::size_t>(B) * B * B;
  std::size_t at[8];
  int blk[8];
  for (int c = 0; c < 8; ++c) {
    const int i = i0[0] + (c & 1), j = i0[1] + (c >> 1 & 1), k = i0[2] + (c >> 2);
    blk[c] = f.mode == 2 ? f.block[(static_cast<std::size_t>(k / B) * f.nb[1] + j / B) * f.nb[0] + i / B] : -1;
    at[c] = blk[c] < 0 ? pointIndex(f, i, j, k)
                       : static_cast<std::size_t>(blk[c]) * f.mapCount * bs + ((k % B) * B + j % B) * B + i % B;
  }
  const T fx = fr[0], fy = fr[1], fz = fr[2];
  auto lerp = [](T a, T b, T t) { return a + (b - a) * t; };
#ifdef __CUDA_ARCH__
  // Packed maps when the corners are all of the reference maps.
  bool packed = f.packed != nullptr;
  for (int c = 0; packed && c < 8; ++c) packed = blk[c] < 0;
#endif
  for (int q = 0; q < count; ++q) {
    const int m = maps[q];
    T v[2][2][2];
#ifdef __CUDA_ARCH__
    if (packed) {
      // Two loads of 16 bytes: the corners x = 0 and x = 1.
      const float4* pk = reinterpret_cast<const float4*>(f.packed) + static_cast<std::size_t>(m) * f.n[0] * f.n[1] * f.n[2] + at[0];
      const float4 a = __ldg(pk), b = __ldg(pk + 1);
      v[0][0][0] = a.x, v[0][1][0] = a.y, v[1][0][0] = a.z, v[1][1][0] = a.w;
      v[0][0][1] = b.x, v[0][1][1] = b.y, v[1][0][1] = b.z, v[1][1][1] = b.w;
    } else
#endif
    for (int c = 0; c < 8; ++c)
      v[c >> 2][c >> 1 & 1][c & 1] = blk[c] < 0 ? f.maps[m][at[c]] : f.data[at[c] + static_cast<std::size_t>(m) * bs];
    const T c00 = lerp(v[0][0][0], v[0][0][1], fx), c10 = lerp(v[0][1][0], v[0][1][1], fx);
    const T c01 = lerp(v[1][0][0], v[1][0][1], fx), c11 = lerp(v[1][1][0], v[1][1][1], fx);
    const T c0 = lerp(c00, c10, fy), c1 = lerp(c01, c11, fy);
    if (grad) {
      const T dx0 = lerp(v[0][0][1] - v[0][0][0], v[0][1][1] - v[0][1][0], fy);
      const T dx1 = lerp(v[1][0][1] - v[1][0][0], v[1][1][1] - v[1][1][0], fy);
      grad[3 * q] = lerp(dx0, dx1, fz) * isp;
      grad[3 * q + 1] = lerp(c10 - c00, c11 - c01, fz) * isp;
      grad[3 * q + 2] = (c1 - c0) * isp;
    }
    val[q] = lerp(c0, c1, fz);
  }
}

// Calls fn(atom) for every atom of the cells overlapping the cube p ± r (detail::CellList::forNear).
template <class F>
OD_HD inline void forNear(const FieldCells& c, const double p[3], double r, F&& fn) {
  using std::floor;
  if (c.atomCount == 0) return;
  int a[3], b[3];
  for (int k = 0; k < 3; ++k) {
    const int lo = static_cast<int>(floor((p[k] - r - c.lo[k]) / c.cell));
    const int hi = static_cast<int>(floor((p[k] + r - c.lo[k]) / c.cell));
    a[k] = lo > 0 ? lo : 0;
    b[k] = hi < c.n[k] - 1 ? hi : c.n[k] - 1;
    if (a[k] > b[k]) return;
  }
  for (int z = a[2]; z <= b[2]; ++z)
    for (int y = a[1]; y <= b[1]; ++y) {
      const std::size_t row = (static_cast<std::size_t>(z) * c.n[1] + y) * c.n[0];
      for (int k = c.start[row + a[0]]; k < c.start[row + b[0] + 1]; ++k) fn(c.atoms[k]);
    }
}

// flexres: pair sum over the moved StaM atoms; adds dE/dp to grad.
template <class T>
OD_HD inline T movedSum(const FieldView& f, int m, T q, const T p[3], T* grad) {
  using std::abs;
  using std::sqrt;
  const T aq = abs(q);
  const T rMin = static_cast<T>(f.rMin), half = static_cast<T>(f.halfSmooth);
  const T cut2 = static_cast<T>(kElecCutoff * kElecCutoff), vdwCut = static_cast<T>(kVdwCutoff);
  const T solT = static_cast<T>(f.solT[m]), volT = static_cast<T>(f.volT[m]);
  const T cd = static_cast<T>(f.cDesolv), cdq = static_cast<T>(f.cDesolv * f.qasp);
  double sum = 0;
  double gs[3] = {0, 0, 0};
  const double pd[3] = {static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2])};
  forNear(f.cells, pd, kElecCutoff, [&](int j) {
    const FieldMoved& at = f.moved[j];
    const T d[3] = {p[0] - static_cast<T>(at.xyz[0]), p[1] - static_cast<T>(at.xyz[1]), p[2] - static_cast<T>(at.xyz[2])};
    const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (r2 >= cut2) return;
    const T r = sqrt(r2);
    T gauss, dgauss, el, del;
    ad4Desolv(r, T(1), rMin, gauss, dgauss);
    ad4Elec(r, static_cast<T>(at.elecQ), rMin, el, del);
    const T V_ = static_cast<T>(at.V), S_ = static_cast<T>(at.S);
    const T w = cd * (solT * V_ + S_ * volT) + aq * cdq * V_;
    T e = w * gauss + q * el, dedr = w * dgauss + q * del;
    if (r < vdwCut) {
      const std::size_t k = static_cast<std::size_t>(m) * f.nTypes + at.type;
      T ev, dv;
      ad4VdwHb(r, static_cast<T>(f.pairAC[k]), static_cast<T>(f.pairBD[k]), f.pairM[k], f.pairN[k],
               static_cast<T>(f.pairReqm[k]), half, ev, dv);
      e += ev;
      dedr += dv;
    }
    sum += e;
    if (grad && r > 0)
      for (int i = 0; i < 3; ++i) gs[i] += dedr / r * d[i];
  });
  if (grad)
    for (int i = 0; i < 3; ++i) grad[i] = static_cast<T>(grad[i] + gs[i]);
  return static_cast<T>(sum);
}

// Energy of a probe atom of map m and charge q at p (inside the grid); grad: dE/dp.
template <class T>
OD_HD inline T energy(const FieldView& f, int m, T q, const T p[3], T* grad) {
  using std::abs;
  const int nt = f.typeMapCount;
  T a[3] = {0, 0, 0}, b[3] = {0, 0, 0}, s[3] = {0, 0, 0};
  const T aq = abs(q);
  const T grid = interpolate(f, m, p, grad ? a : nullptr) + q * interpolate(f, nt, p, grad ? b : nullptr) +
                 aq * interpolate(f, nt + 1, p, grad ? s : nullptr);
  if (grad)
    for (int i = 0; i < 3; ++i) grad[i] = a[i] + q * b[i] + aq * s[i];
  if (f.mode != 1 || f.movedCount == 0) return grid;
  return static_cast<T>(grid + movedSum(f, m, q, p, grad));
}

// flexres part of distance: d, the map value, or the distance to the nearest moved StaM heavy atom
// when closer.
template <class T>
OD_HD inline T movedNearest(const FieldView& f, const T p[3], T d, T* grad);

// Distance to the nearest StaM heavy atom (at most the cap); grad: dd/dp.
template <class T>
OD_HD inline T distance(const FieldView& f, const T p[3], T* grad) {
  const T d = interpolate(f, f.typeMapCount + 2, p, grad);
  if (f.mode != 1 || f.movedHeavyCount == 0) return d;
  return movedNearest(f, p, d, grad);
}

// energy (map m, charge q) and, with dist, distance at p, the maps interpolated together
// (interpolateMaps; the same values as energy and distance).
template <class T>
OD_HD inline T energyDistance(const FieldView& f, int m, T q, const T p[3], T* grad, T* dist, T* dgrad) {
  using std::abs;
  const int nt = f.typeMapCount;
  const int maps[4] = {m, nt, nt + 1, nt + 2};
  T val[4], g[12];
  interpolateMaps(f, maps, dist ? 4 : 3, p, val, grad ? g : nullptr);
  const T aq = abs(q);
  const T grid = val[0] + q * val[1] + aq * val[2];
  if (grad)
    for (int i = 0; i < 3; ++i) grad[i] = g[i] + q * g[3 + i] + aq * g[6 + i];
  if (dist) {
    if (dgrad)
      for (int i = 0; i < 3; ++i) dgrad[i] = g[9 + i];
    *dist = f.mode != 1 || f.movedHeavyCount == 0 ? val[3] : movedNearest(f, p, val[3], dgrad);
  }
  if (f.mode != 1 || f.movedCount == 0) return grid;
  return static_cast<T>(grid + movedSum(f, m, q, p, grad));
}

template <class T>
OD_HD inline T movedNearest(const FieldView& f, const T p[3], T d, T* grad) {
  using std::sqrt;
  int best = -1;
  T best2 = d * d;
  const double pd[3] = {static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2])};
  forNear(f.heavyCells, pd, f.cap, [&](int h) {
    const double* x = f.moved[f.movedHeavy[h]].xyz;
    const T dx = p[0] - static_cast<T>(x[0]), dy = p[1] - static_cast<T>(x[1]), dz = p[2] - static_cast<T>(x[2]);
    const T r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < best2) {
      best2 = r2;
      best = f.movedHeavy[h];
    }
  });
  if (best < 0) return d;
  const T r = sqrt(best2);
  if (grad) {
    const double* x = f.moved[best].xyz;
    for (int i = 0; i < 3; ++i) grad[i] = r > 0 ? (p[i] - static_cast<T>(x[i])) / r : T(0);
  }
  return r;
}

}  // namespace field
}  // namespace od
