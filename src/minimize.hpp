// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>
#include <span>
#include <vector>

namespace od {

// Torsion-space minimization (spec 12). Variables are angles in degrees; the
// objective returns the energy and its gradient per radian.
using TorsionObjective = std::function<double(std::span<const double> anglesDeg, std::vector<double>& gradRad)>;

struct MinParams {
  int maxSteps = 100;         // confMinMaxStep
  double deltaE = 10.0;       // confMinDeltaE: stop when one step changes E by less
  double maxStepDeg = 10.0;   // confMaxDihedralStep: largest change of one angle per step
};

struct MinResult {
  double energy = 0;  // energy at the returned angles
  int steps = 0;      // position updates
  int evaluations = 0;
};

// FIRE (Bitzek et al. 2006) with N_min = 5, f_inc = 1.1, f_dec = 0.5,
// alpha0 = 0.1, f_alpha = 0.99, dt0 = 0.1, dt_max = 1.0 (unit mass). A step whose
// largest angle change exceeds maxStepDeg is scaled down as a whole. The
// lowest-energy point visited is returned in anglesDeg, wrapped to (-180, 180].
MinResult minimizeFire(const TorsionObjective& f, std::span<double> anglesDeg, const MinParams& p);

}  // namespace od
