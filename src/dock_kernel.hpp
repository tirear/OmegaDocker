// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include "energy_kernel.hpp"
#include "field_kernel.hpp"

namespace od {

// Flat form of the 9b pose energy (MobmDocker::evaluate, spec 8.6) shared by the CPU
// and GPU back ends (spec 15.1.4). Single precision per item; sums in double.

struct DockProbe {
  int coord;  // coordinate index in the MobM energy model
  int map;    // type map
  double charge;
  int heavy;
};

struct DockView {
  int probeCount;
  const DockProbe* probes;
  double clashDistance, clashBase, clashK;
  FieldView field;
  const FieldView* fields;  // several fields (device array, DockPose::field); nullptr: field
};

// World x = R(q) (y − center) + t.
struct DockPose {
  float t[3];
  float q[4];  // w, x, y, z
  float center[3];
  int field;  // index into DockView::fields
};

struct DockTerms {
  double inter, interClash;
  int outside;
};

// Whole 9b minimizations with L-BFGS on the device (MobmDocker::minimize, spec 15.1.26).
// Variables z = (t / sT, rotation vector / sR, torsions in rad / sA); steps limited to maxTrans
// (Å) and maxAngle (rad); VecMinParams maxSteps, convSteps, deltaE and memory.
struct DockMinParams {
  double sT, sR, sA, maxTrans, maxAngle, deltaE;
  int maxSteps, convSteps, memory;
};
// Start of one minimization: pose (double), centre (the Pose frame) and field as DockPose.
struct DockMinStart {
  double t[3], q[4];
  float center[3];
  int field;
};
// The lowest-energy accepted point (float pose) and its energies; ok 0: the pose left the grid
// (only steps and evaluations are set).
struct DockMinResult {
  float t[3], q[4];
  FlatTerms terms;
  DockTerms dterms;
  int steps, evaluations, ok;
};

namespace dock {

OD_HD inline void rotate(const float q[4], const float v[3], float out[3]) {
  const float tx = q[2] * v[2] - q[3] * v[1] + q[0] * v[0];
  const float ty = q[3] * v[0] - q[1] * v[2] + q[0] * v[1];
  const float tz = q[1] * v[1] - q[2] * v[0] + q[0] * v[2];
  out[0] = v[0] + 2 * (q[2] * tz - q[3] * ty);
  out[1] = v[1] + 2 * (q[3] * tx - q[1] * tz);
  out[2] = v[2] + 2 * (q[1] * ty - q[2] * tx);
}

// Probe atom k at the coordinates x (torsion-tree frame): adds its grid energy and clash
// penalty against the field f; gw (optional): dE/dw at the world position; r: R(q) (y − center).
// Returns false when the atom is outside the grid (nothing added).
OD_HD inline bool probeTerm(const DockView& d, const FieldView& f, int k, const float* x, const DockPose& ps,
                            double& inter, double& interClash, float gw[3], float r[3]) {
  const DockProbe& pr = d.probes[k];
  const float* y = x + 3 * pr.coord;
  const float v[3] = {y[0] - ps.center[0], y[1] - ps.center[1], y[2] - ps.center[2]};
  rotate(ps.q, v, r);
  const float w[3] = {r[0] + ps.t[0], r[1] + ps.t[1], r[2] + ps.t[2]};
  bool in = true;
  for (int c = 0; c < 3; ++c) {
    const float lo = static_cast<float>(f.origin[c]), hi = static_cast<float>(f.origin[c] + (f.n[c] - 1) * f.spacing);
    in &= w[c] >= lo && w[c] <= hi;
  }
  if (!in) return false;
  const float q = static_cast<float>(pr.charge);
  float g[3] = {0, 0, 0}, dg[3] = {0, 0, 0}, dd = 0;
  inter += field::energyDistance<float>(f, pr.map, q, w, gw ? g : nullptr, pr.heavy ? &dd : nullptr, gw ? dg : nullptr);
  if (pr.heavy) {
    const float cd = static_cast<float>(d.clashDistance), base = static_cast<float>(d.clashBase),
                kk = static_cast<float>(d.clashK);
    if (dd < cd) {
      interClash += base + kk * (cd - dd) * (cd - dd);
      for (int i = 0; i < 3; ++i) g[i] += -2 * kk * (cd - dd) * dg[i];
    }
  }
  if (gw)
    for (int i = 0; i < 3; ++i) gw[i] = g[i];
  return true;
}

// Gradient of probe k: adds to gt, gr (translation, rotation vector) and returns
// gy = R(q)^T gw, the gradient on its coordinates.
OD_HD inline void probeGrad(const DockPose& ps, const float gw[3], const float r[3], double gt[3], double gr[3],
                            float gy[3]) {
  for (int i = 0; i < 3; ++i) gt[i] += gw[i];
  gr[0] += r[1] * gw[2] - r[2] * gw[1];
  gr[1] += r[2] * gw[0] - r[0] * gw[2];
  gr[2] += r[0] * gw[1] - r[1] * gw[0];
  const float qc[4] = {ps.q[0], -ps.q[1], -ps.q[2], -ps.q[3]};
  rotate(qc, gw, gy);
}

// Serial evaluation (CPU check of the flat form). x: 3 × coordCount; g: 3 × coordCount
// when grad is requested; grad: 6 + torsionCount (translation, rotation vector, torsions).
inline FlatTerms evaluate(const EnergyModelView& m, const DockView& d, const float* angles, const DockPose& ps,
                          float* x, double* g, double* grad, DockTerms& dt) {
  dt = {0, 0, 0};
  double gt[3] = {0, 0, 0}, gr[3] = {0, 0, 0};
  const FlatTerms t = flat::evaluate(m, angles, x, g, grad ? grad + 6 : nullptr, [&](const float* xc, double* gc) {
    for (int k = 0; k < d.probeCount; ++k) {
      float gw[3], r[3], gy[3];
      if (!probeTerm(d, d.field, k, xc, ps, dt.inter, dt.interClash, gc ? gw : nullptr, r)) {
        ++dt.outside;
        continue;
      }
      if (!gc) continue;
      probeGrad(ps, gw, r, gt, gr, gy);
      for (int i = 0; i < 3; ++i) gc[3 * d.probes[k].coord + i] += gy[i];
    }
  });
  if (grad)
    for (int i = 0; i < 3; ++i) {
      grad[i] = gt[i];
      grad[3 + i] = gr[i];
    }
  return t;
}

}  // namespace dock
}  // namespace od
