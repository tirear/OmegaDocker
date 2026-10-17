// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>
#include <span>
#include <string_view>
#include <vector>

namespace od {

// Minimization methods (spec 12).
enum class MinMethod { Sd, Cg, Fire, Lbfgs, Adadelta };
MinMethod parseMinMethod(std::string_view s);  // sd | cg | fire | lbfgs | adadelta
const char* minMethodName(MinMethod m);

// Objective of the general minimizer: energy and gradient at x; returns false to
// abort the minimization (e.g. a pose leaving the grid, spec 8.6).
using VecObjective = std::function<bool(std::span<const double> x, double& energy, std::vector<double>& grad)>;

struct VecMinParams {
  MinMethod method = MinMethod::Lbfgs;
  int maxSteps = 1000;    // accepted steps
  int convSteps = 10;     // stop after this many consecutive steps with |dE| < deltaE
  double deltaE = 1e-3;
  int memory = 7;         // lbfgsMemory
  // Largest factor (<= 1) applied to a step, enforcing per-step limits.
  std::function<double(std::span<const double> step)> stepLimit;
  // Called after each accepted step; best: the point is the lowest so far. May re-express x
  // (e.g. re-base a rotation).
  std::function<void(std::vector<double>& x, bool best)> onAccept;
};

struct VecMinResult {
  double energy = 0;     // lowest energy among the accepted points (and the start)
  int steps = 0;         // accepted steps
  int evaluations = 0;
  bool aborted = false;  // the objective returned false; x is then unspecified
};

// Minimizes from x; on return x holds the lowest-energy accepted point as it was when accepted.
// sd, cg and lbfgs only accept energy decreases, so that point is the last accepted one.
VecMinResult minimizeVec(const VecObjective& f, std::vector<double>& x, const VecMinParams& p);

// Torsion-space minimization (spec 7.9). Variables are angles in degrees; the
// objective returns the energy and its gradient per radian.
using TorsionObjective = std::function<double(std::span<const double> anglesDeg, std::vector<double>& gradRad)>;

struct MinParams {
  int maxSteps = 100;         // confMinMaxStep
  double deltaE = 10.0;       // confMinDeltaE: stop when one step changes E by less
  double maxStepDeg = 10.0;   // confMaxDihedralStep: largest change of one angle per step
  MinMethod method = MinMethod::Lbfgs;
};

struct MinResult {
  double energy = 0;  // energy at the returned angles
  int steps = 0;      // accepted steps
  int evaluations = 0;
};

// A step whose largest angle change exceeds maxStepDeg is scaled down as a whole. The
// lowest-energy accepted point is returned in anglesDeg, wrapped to (-180, 180].
MinResult minimizeTorsions(const TorsionObjective& f, std::span<double> anglesDeg, const MinParams& p);

}  // namespace od
