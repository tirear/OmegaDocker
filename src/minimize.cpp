// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "minimize.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "rng.hpp"

namespace od {

MinResult minimizeFire(const TorsionObjective& f, std::span<double> anglesDeg, const MinParams& p) {
  constexpr int kNmin = 5;
  constexpr double kFinc = 1.1, kFdec = 0.5, kAlpha0 = 0.1, kFalpha = 0.99, kDt0 = 0.1, kDtMax = 1.0;
  constexpr double kRad = std::numbers::pi / 180.0;
  const std::size_t n = anglesDeg.size();
  std::vector<double> x(anglesDeg.begin(), anglesDeg.end()), v(n, 0.0), g(n), best(x);
  const double maxStep = p.maxStepDeg * kRad;
  double dt = kDt0, alpha = kAlpha0;
  int npos = 0;
  MinResult r;
  double e = f(x, g);
  r.evaluations = 1;
  double bestE = e;
  while (r.steps < p.maxSteps) {
    double P = 0, vn = 0, fn = 0;
    for (std::size_t i = 0; i < n; ++i) {
      P -= g[i] * v[i];
      vn += v[i] * v[i];
      fn += g[i] * g[i];
    }
    if (P > 0) {
      vn = std::sqrt(vn);
      fn = std::sqrt(fn);
      if (fn > 0)
        for (std::size_t i = 0; i < n; ++i) v[i] = (1 - alpha) * v[i] - alpha * vn * g[i] / fn;
      if (++npos > kNmin) {
        dt = std::min(dt * kFinc, kDtMax);
        alpha *= kFalpha;
      }
    } else {
      std::fill(v.begin(), v.end(), 0.0);
      dt *= kFdec;
      alpha = kAlpha0;
      npos = 0;
    }
    double big = 0;
    for (std::size_t i = 0; i < n; ++i) {
      v[i] -= dt * g[i];
      big = std::max(big, std::abs(dt * v[i]));
    }
    const double scale = big > maxStep ? maxStep / big : 1.0;
    for (std::size_t i = 0; i < n; ++i) x[i] += scale * dt * v[i] / kRad;
    const double eNew = f(x, g);
    ++r.evaluations;
    ++r.steps;
    if (eNew < bestE) {
      bestE = eNew;
      best = x;
    }
    const double dE = std::abs(eNew - e);
    e = eNew;
    if (dE < p.deltaE) break;
  }
  for (std::size_t i = 0; i < n; ++i) anglesDeg[i] = wrapDeg(best[i]);
  r.energy = bestE;
  return r;
}

}  // namespace od
