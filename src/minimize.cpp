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

LbfgsResult minimizeLbfgs(const VecObjective& f, std::vector<double>& x, const LbfgsParams& p) {
  constexpr double kArmijo = 1e-4;
  constexpr int kMaxBacktrack = 20;
  const std::size_t n = x.size();
  LbfgsResult r;
  std::vector<double> g(n), gNew(n), d(n), xNew(n);
  double e;
  ++r.evaluations;
  if (!f(x, e, g)) {
    r.aborted = true;
    return r;
  }
  auto dot = [n](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
  };
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  int calm = 0;
  while (r.steps < p.maxSteps) {
    // Two-loop recursion: d = -H g.
    d = g;
    const std::size_t m = S.size();
    std::vector<double> alpha(m);
    for (std::size_t k = m; k-- > 0;) {
      alpha[k] = rho[k] * dot(S[k], d);
      for (std::size_t i = 0; i < n; ++i) d[i] -= alpha[k] * Y[k][i];
    }
    if (m > 0) {
      const double gamma = dot(S[m - 1], Y[m - 1]) / dot(Y[m - 1], Y[m - 1]);
      for (double& v : d) v *= gamma;
    }
    for (std::size_t k = 0; k < m; ++k) {
      const double beta = rho[k] * dot(Y[k], d);
      for (std::size_t i = 0; i < n; ++i) d[i] += (alpha[k] - beta) * S[k][i];
    }
    for (double& v : d) v = -v;
    double slope = dot(g, d);
    if (!(slope < 0)) {  // not a descent direction: restart from steepest descent
      S.clear();
      Y.clear();
      rho.clear();
      for (std::size_t i = 0; i < n; ++i) d[i] = -g[i];
      slope = dot(g, d);
      if (!(slope < 0)) break;  // zero gradient
    }
    if (p.stepLimit) {
      const double s = std::min(1.0, p.stepLimit(d));
      for (double& v : d) v *= s;
      slope *= s;
    }
    double a = 1, eNew = 0;
    bool accepted = false;
    for (int b = 0; b < kMaxBacktrack; ++b, a *= 0.5) {
      for (std::size_t i = 0; i < n; ++i) xNew[i] = x[i] + a * d[i];
      ++r.evaluations;
      if (!f(xNew, eNew, gNew)) {
        r.aborted = true;
        return r;
      }
      if (eNew <= e + kArmijo * a * slope) {
        accepted = true;
        break;
      }
    }
    if (!accepted) break;  // no decrease along the direction: converged
    ++r.steps;
    std::vector<double> s(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
      s[i] = xNew[i] - x[i];
      y[i] = gNew[i] - g[i];
    }
    const double sy = dot(s, y);
    if (sy > 1e-12) {
      if (static_cast<int>(S.size()) == p.memory) {
        S.erase(S.begin());
        Y.erase(Y.begin());
        rho.erase(rho.begin());
      }
      S.push_back(std::move(s));
      Y.push_back(std::move(y));
      rho.push_back(1 / sy);
    }
    const double dE = std::abs(eNew - e);
    x = xNew;
    g = gNew;
    e = eNew;
    if (p.onAccept) p.onAccept(x);
    calm = dE < p.deltaE ? calm + 1 : 0;
    if (calm >= p.convSteps) break;
  }
  r.energy = e;
  return r;
}

}  // namespace od
