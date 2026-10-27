// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cmath>

#include "ad4_pair.hpp"
#include "field_kernel.hpp"

namespace od {

// Energy grid values at one point (buildEnergyGrids, spec 8.1); shared by the CPU and the GPU
// (spec 15.1.10). Arrays per StaM atom (cell list order independent), per map, per (map, Ad4Type).
constexpr int kGridMaxTypeMaps = kFieldMaxMaps - 3;
// Blocks of kDeltaBlock^3 points (spec 8.1 delta).
constexpr int kDeltaBlock = 8;

// Delta blocks reached by a changed atom position (spec 15.1.19): squared distance from the centre
// of a block within which a position may contribute to one of its points (reach: the largest
// cutoff of the maps).
OD_HD inline double deltaReach2(double spacing, double reach) {
  using std::sqrt;
  const double h = reach + 0.5 * sqrt(3.0) * (kDeltaBlock - 1) * spacing;
  return h * h;
}
OD_HD inline bool deltaBlockReached(const double* origin, double spacing, int bx, int by, int bz, const double* p,
                                    double reach2) {
  const double c = 0.5 * (kDeltaBlock - 1);
  const double x = p[0] - (origin[0] + (bx * kDeltaBlock + c) * spacing);
  const double y = p[1] - (origin[1] + (by * kDeltaBlock + c) * spacing);
  const double z = p[2] - (origin[2] + (bz * kDeltaBlock + c) * spacing);
  return x * x + y * y + z * z <= reach2;
}

struct GridBuildView {
  double origin[3], spacing;
  int n[3];
  int typeMapCount, nTypes;
  const double* xyz;    // 3 per StaM atom
  const int* type;      // Ad4Type
  const double* elecQ;  // c_elec q
  const double* S;      // solpar + qasp |q|
  const double* V;
  const double* solT;   // per map
  const double* volT;
  const double* desolvW;  // per (StaM atom, map): cDesolv (solT V + S volT)
  const double* dsW;      // per StaM atom: cDesolvQasp V
  const double* AC;     // per (map, Ad4Type)
  const double* BD;
  const double* reqm;
  const int* pm;
  const int* pn;
  double cDesolv, cDesolvQasp, halfSmooth, rMin;
  FieldCells cells;     // over the StaM atoms, cell 4 Å
};

// e: typeMapCount values (unclamped); el, ds: elec and desolv map values. MaxT >= typeMapCount:
// with a small MaxT the GPU keeps e in registers (spec 15.1.19); the values do not depend on it.
#ifdef __CUDA_ARCH__
#define OD_UNROLL _Pragma("unroll")
#else
#define OD_UNROLL
#endif
template <int MaxT = kGridMaxTypeMaps>
OD_HD inline void gridPoint(const GridBuildView& v, int i, int j, int k, double* e, double& el, double& ds) {
  using std::sqrt;
  const double p[3] = {v.origin[0] + i * v.spacing, v.origin[1] + j * v.spacing, v.origin[2] + k * v.spacing};
  const int nt = v.typeMapCount;
  OD_UNROLL
  for (int t = 0; t < MaxT; ++t)
    if (t < nt) e[t] = 0.0;
  el = 0;
  ds = 0;
  field::forNear(v.cells, p, kElecCutoff, [&](int a) {
    const double x = p[0] - v.xyz[3 * a], y = p[1] - v.xyz[3 * a + 1], z = p[2] - v.xyz[3 * a + 2];
    const double r2 = x * x + y * y + z * z;
    if (r2 >= kElecCutoff * kElecCutoff) return;
    const double r = sqrt(r2);
    double ee, de;
    ad4Elec(r, v.elecQ[a], v.rMin, ee, de);
    el += ee;
    double gauss;
    ad4Desolv(r, 1.0, v.rMin, gauss, de);
    ds += v.dsW[a] * gauss;
    const int ta = v.type[a];
    const double* w = v.desolvW + static_cast<std::size_t>(a) * nt;
    OD_UNROLL
    for (int t = 0; t < MaxT; ++t) {
      if (t >= nt) break;
      e[t] += w[t] * gauss;
      if (r < kVdwCutoff) {
        const int q = t * v.nTypes + ta;
        ad4VdwHb(r, v.AC[q], v.BD[q], v.pm[q], v.pn[q], v.reqm[q], v.halfSmooth, ee, de);
        e[t] += ee;
      }
    }
  });
}

// Distance grid value at one point (buildDistanceGrid, spec 8.2): distance to the nearest
// StaM heavy atom other than a metal, at most cap.
struct DistBuildView {
  double origin[3], spacing;
  int n[3];
  double cap;
  const double* xyz;  // 3 per heavy atom
  FieldCells cells;   // over the heavy atoms, cell max(cap, 1) Å
};

OD_HD inline float distPoint(const DistBuildView& v, int i, int j, int k) {
  using std::sqrt;
  const double p[3] = {v.origin[0] + i * v.spacing, v.origin[1] + j * v.spacing, v.origin[2] + k * v.spacing};
  double best = v.cap * v.cap;
  field::forNear(v.cells, p, v.cap, [&](int a) {
    const double x = p[0] - v.xyz[3 * a], y = p[1] - v.xyz[3 * a + 1], z = p[2] - v.xyz[3 * a + 2];
    const double d2 = x * x + y * y + z * z;
    best = d2 < best ? d2 : best;  // std::min(best, d2)
  });
  return static_cast<float>(sqrt(best));
}

}  // namespace od
