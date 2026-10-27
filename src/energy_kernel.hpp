// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cmath>
#include <vector>

#include "ad4_pair.hpp"

namespace od {

// Flat form of an EnergySystem (spec 7.9) shared by the CPU and GPU back ends
// (spec 15.1). The per-item functions below follow EnergySystem::evaluate
// operation by operation; a serial loop over them gives the same result.

struct FlatPairParams {
  double AC, BD, reqm;
  int m, n;
};

struct FlatPair {
  int i, j;   // coordinate indices
  int type;   // index into pairParams
  int clash;  // both heavy atoms, neither a metal: clash penalty applies
};

// Torsions of all parts in part order; index = angle index. Applied last to first.
struct FlatTorsion {
  int fixedAtom, movingAtom;
  int movingBegin, movingEnd;  // range in moving
  double inputAngle;           // degrees
};

// Closure terms in the order of EnergySystem::evaluate: per part, virtual
// atoms (kind 0: G·|xi − xj|²) then distance restraints (kind 1: G·(|xi − xj| − r0)²).
struct FlatClosure {
  int i, j;
  int kind;
  double r0;
};

struct FlatTerms {
  double vdw, hbond, elec, desolv, clash, closure, closureMaxDeviation;
  int clashPairs;
};

struct EnergyModelView {
  int coordCount, torsionCount, pairCount, closureCount;
  const double* ref;  // 3 × coordCount: input coordinates (atoms, then virtual atoms)
  const FlatTorsion* torsions;
  const int* moving;
  const double* charge;
  const double* sv;
  const double* vol;
  const FlatPairParams* pairParams;
  const FlatPair* pairs;
  const FlatClosure* closures;
  double halfSmooth, clashDistance, clashBase, clashK, closureG, elecCoeff, desolvCoeff, rMin;
};

struct EnergyModel {
  int coordCount = 0, torsionCount = 0;
  std::vector<double> ref;
  std::vector<FlatTorsion> torsions;
  std::vector<int> moving;
  std::vector<double> charge, sv, vol;
  std::vector<FlatPairParams> pairParams;
  std::vector<FlatPair> pairs;
  std::vector<FlatClosure> closures;
  double halfSmooth = 0, clashDistance = 0, clashBase = 0, clashK = 0, closureG = 0, elecCoeff = 0, desolvCoeff = 0,
         rMin = 0;

  EnergyModelView view() const {
    return {coordCount,       torsionCount, static_cast<int>(pairs.size()), static_cast<int>(closures.size()),
            ref.data(),       torsions.data(), moving.data(), charge.data(), sv.data(), vol.data(),
            pairParams.data(), pairs.data(), closures.data(), halfSmooth, clashDistance, clashBase, clashK, closureG,
            elecCoeff,        desolvCoeff,  rMin};
  }
};

namespace flat {

template <class T>
OD_HD inline T piT() {
  return T(3.141592653589793);
}

// Cos and sin of the rotation of torsion t.
template <class T>
OD_HD inline void torsionCosSin(const EnergyModelView& m, int t, const T* angles, T& c, T& s) {
  using std::cos;
  using std::sin;
  const T delta = (angles[t] - static_cast<T>(m.torsions[t].inputAngle)) * piT<T>() / T(180);
  c = cos(delta);
  s = sin(delta);
}

// Axis point o and unit direction k of torsion t; false for a zero-length axis.
template <class T>
OD_HD inline bool torsionAxisOnly(const EnergyModelView& m, int t, const T* x, T o[3], T k[3]) {
  using std::sqrt;
  const FlatTorsion& tt = m.torsions[t];
  for (int a = 0; a < 3; ++a) {
    o[a] = x[3 * tt.fixedAtom + a];
    k[a] = x[3 * tt.movingAtom + a] - o[a];
  }
  const T len = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
  if (len == 0) return false;
  for (int a = 0; a < 3; ++a) k[a] /= len;
  return true;
}

// Axis point o, unit direction k, cos and sin of the rotation of torsion t.
// Returns false for a zero-length axis (the torsion is skipped).
template <class T>
OD_HD inline bool torsionAxis(const EnergyModelView& m, int t, const T* angles, const T* x, T o[3], T k[3], T& c, T& s) {
  if (!torsionAxisOnly(m, t, x, o, k)) return false;
  torsionCosSin(m, t, angles, c, s);
  return true;
}

template <class T>
OD_HD inline void rotateAtom(T* p, const T o[3], const T k[3], T c, T s) {
  const T v[3] = {p[0] - o[0], p[1] - o[1], p[2] - o[2]};
  const T dot = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
  const T cr[3] = {k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]};
  for (int i = 0; i < 3; ++i) p[i] = o[i] + v[i] * c + cr[i] * s + k[i] * dot * (1 - c);
}

// One atom pair: adds to e[0..4] (vdw, hbond, elec, desolv, clash) and
// clashPairs; d = xi − xj, r, dedr. Returns false beyond the cutoff.
template <class T>
OD_HD inline bool pairTerm(const EnergyModelView& m, int p, const T* x, double e[5], int& clashPairs, T d[3], T& r,
                           T& dedr) {
  using std::sqrt;
  const FlatPair& pr = m.pairs[p];
  for (int a = 0; a < 3; ++a) d[a] = x[3 * pr.i + a] - x[3 * pr.j + a];
  const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  const T elecCut2 = static_cast<T>(kElecCutoff * kElecCutoff);
  if (r2 >= elecCut2) return false;
  r = sqrt(r2);
  dedr = 0;
  T ee, de;
  if (r < static_cast<T>(kVdwCutoff)) {
    const FlatPairParams& pp = m.pairParams[pr.type];
    ad4VdwHb<T>(r, static_cast<T>(pp.AC), static_cast<T>(pp.BD), pp.m, pp.n, static_cast<T>(pp.reqm),
                static_cast<T>(m.halfSmooth), ee, de);
    e[pp.n == 10 ? 1 : 0] += ee;
    dedr += de;
  }
  const T rMin = static_cast<T>(m.rMin);
  ad4Elec<T>(r, static_cast<T>(m.elecCoeff) * static_cast<T>(m.charge[pr.i]) * static_cast<T>(m.charge[pr.j]), rMin, ee,
             de);
  e[2] += ee;
  dedr += de;
  ad4Desolv<T>(r,
               static_cast<T>(m.desolvCoeff) * (static_cast<T>(m.sv[pr.i]) * static_cast<T>(m.vol[pr.j]) +
                                                static_cast<T>(m.sv[pr.j]) * static_cast<T>(m.vol[pr.i])),
               rMin, ee, de);
  e[3] += ee;
  dedr += de;
  const T cd = static_cast<T>(m.clashDistance);
  if (pr.clash && r2 < cd * cd) {
    clashPenalty<T>(r, cd, static_cast<T>(m.clashBase), static_cast<T>(m.clashK), ee, de);
    e[4] += ee;
    ++clashPairs;
    dedr += de;
  }
  return true;
}

// One closure term: adds to energy and maxDev; the gradient on xi is f·d
// (−f·d on xj); f = 0 when no gradient applies.
template <class T>
OD_HD inline void closureTerm(const EnergyModelView& m, int c, const T* x, double& energy, double& maxDev, T d[3],
                              T& f) {
  using std::abs;
  using std::sqrt;
  const FlatClosure& cl = m.closures[c];
  for (int a = 0; a < 3; ++a) d[a] = x[3 * cl.i + a] - x[3 * cl.j + a];
  const T G = static_cast<T>(m.closureG);
  if (cl.kind == 0) {
    const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    energy += G * r2;
    const double dev = static_cast<double>(sqrt(r2));
    if (dev > maxDev) maxDev = dev;
    f = 2 * G;
  } else {
    const T r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    const T dd = r - static_cast<T>(cl.r0);
    energy += G * dd * dd;
    const double dev = static_cast<double>(abs(dd));
    if (dev > maxDev) maxDev = dev;
    f = r > 0 ? (2 * G * dd) / r : T(0);
  }
}

// Contribution of the moving atom idx to the gradient of torsion t: g · (u × v).
template <class T>
OD_HD inline double torsionGradTerm(const T* x, const double* g, int idx, const T a[3], const T u[3]) {
  const T v[3] = {x[3 * idx] - a[0], x[3 * idx + 1] - a[1], x[3 * idx + 2] - a[2]};
  const double* gi = g + 3 * idx;
  return gi[0] * (u[1] * v[2] - u[2] * v[1]) + gi[1] * (u[2] * v[0] - u[0] * v[2]) + gi[2] * (u[0] * v[1] - u[1] * v[0]);
}

// Unit vector from the fixed to the moving atom of torsion t; false for zero length.
template <class T>
OD_HD inline bool torsionUnit(const EnergyModelView& m, int t, const T* x, T a[3], T u[3]) {
  using std::sqrt;
  const FlatTorsion& tt = m.torsions[t];
  for (int c = 0; c < 3; ++c) {
    a[c] = x[3 * tt.fixedAtom + c];
    u[c] = x[3 * tt.movingAtom + c] - a[c];
  }
  const T len = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  if (len == 0) return false;
  for (int c = 0; c < 3; ++c) u[c] /= len;
  return true;
}

struct NoExternal {
  template <class T>
  void operator()(const T*, double*) const {}
};

// Serial evaluation (CPU check of the flat form). x: 3 × coordCount scratch;
// g: 3 × coordCount scratch when grad is requested; grad: torsionCount.
// ext(x, g or nullptr): external term after the closure terms (EnergySystem::evaluate).
template <class T, class Ext = NoExternal>
inline FlatTerms evaluate(const EnergyModelView& m, const T* angles, T* x, double* g, double* grad, Ext&& ext = {}) {
  for (int i = 0; i < 3 * m.coordCount; ++i) x[i] = static_cast<T>(m.ref[i]);
  for (int t = m.torsionCount; t-- > 0;) {
    T o[3], k[3], c, s;
    if (!torsionAxis(m, t, angles, x, o, k, c, s)) continue;
    for (int q = m.torsions[t].movingBegin; q < m.torsions[t].movingEnd; ++q) rotateAtom(x + 3 * m.moving[q], o, k, c, s);
  }
  if (grad)
    for (int i = 0; i < 3 * m.coordCount; ++i) g[i] = 0;
  double e[5] = {0, 0, 0, 0, 0};
  FlatTerms out{};
  for (int p = 0; p < m.pairCount; ++p) {
    T d[3], r, dedr;
    if (!pairTerm(m, p, x, e, out.clashPairs, d, r, dedr)) continue;
    if (!grad || r <= 0) continue;
    const T f = dedr / r;
    const FlatPair& pr = m.pairs[p];
    for (int c = 0; c < 3; ++c) {
      g[3 * pr.i + c] += f * d[c];
      g[3 * pr.j + c] -= f * d[c];
    }
  }
  for (int c = 0; c < m.closureCount; ++c) {
    T d[3], f;
    closureTerm(m, c, x, out.closure, out.closureMaxDeviation, d, f);
    if (!grad) continue;
    const FlatClosure& cl = m.closures[c];
    for (int a = 0; a < 3; ++a) {
      g[3 * cl.i + a] += f * d[a];
      g[3 * cl.j + a] -= f * d[a];
    }
  }
  ext(static_cast<const T*>(x), grad ? g : nullptr);
  out.vdw = e[0];
  out.hbond = e[1];
  out.elec = e[2];
  out.desolv = e[3];
  out.clash = e[4];
  if (grad)
    for (int t = 0; t < m.torsionCount; ++t) {
      T a[3], u[3];
      double s = 0;
      if (torsionUnit(m, t, x, a, u))
        for (int q = m.torsions[t].movingBegin; q < m.torsions[t].movingEnd; ++q) s += torsionGradTerm(x, g, m.moving[q], a, u);
      grad[t] = s;
    }
  return out;
}

}  // namespace flat
}  // namespace od
