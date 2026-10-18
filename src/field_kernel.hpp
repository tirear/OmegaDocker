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

// Host arrays of a FieldView. The maps are not copied: base keeps the grids that hold them.
struct FieldModel {
  std::uint64_t id = 0;               // StamField::id()
  std::shared_ptr<const void> base;   // owner of the reference maps (grids)
  std::shared_ptr<const void> baseDist;
  std::size_t points = 0;             // values per map
  std::vector<int> block;
  std::vector<float> data;
  std::vector<FieldMoved> moved;
  std::vector<int> movedHeavy, cellStart, cellAtoms, heavyStart, heavyAtoms;
  std::vector<double> solT, volT, pairAC, pairBD, pairReqm;
  std::vector<int> pairM, pairN;
  FieldView view{};  // arrays point into this model and the grids

  // Points the array members of view at the vectors above (after copy or move).
  void bind() {
    view.block = block.data();
    view.data = data.data();
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

// Trilinear interpolation of map m at p; grad (optional): dvalue/dp.
template <class T>
OD_HD inline T interpolate(const FieldView& f, int m, const T p[3], T* grad) {
  using std::floor;
  const T spacing = static_cast<T>(f.spacing);
  int i0[3];
  T fr[3];
  for (int c = 0; c < 3; ++c) {
    const T u = (p[c] - static_cast<T>(f.origin[c])) / spacing;
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
    grad[0] = lerp(dx0, dx1, fz) / spacing;
    grad[1] = lerp(c10 - c00, c11 - c01, fz) / spacing;
    grad[2] = (c1 - c0) / spacing;
  }
  return lerp(c0, c1, fz);
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

// Distance to the nearest StaM heavy atom (at most the cap); grad: dd/dp.
template <class T>
OD_HD inline T distance(const FieldView& f, const T p[3], T* grad) {
  using std::sqrt;
  T d = interpolate(f, f.typeMapCount + 2, p, grad);
  if (f.mode != 1 || f.movedHeavyCount == 0) return d;
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
