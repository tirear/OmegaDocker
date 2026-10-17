// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "minimize.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "log.hpp"
#include "rng.hpp"

namespace od {

MinMethod parseMinMethod(std::string_view s) {
  if (s == "sd") return MinMethod::Sd;
  if (s == "cg") return MinMethod::Cg;
  if (s == "fire") return MinMethod::Fire;
  if (s == "lbfgs") return MinMethod::Lbfgs;
  if (s == "adadelta") return MinMethod::Adadelta;
  fail("expected sd, cg, fire, lbfgs or adadelta but got '{}'", s);
}

const char* minMethodName(MinMethod m) {
  switch (m) {
    case MinMethod::Sd: return "sd";
    case MinMethod::Cg: return "cg";
    case MinMethod::Fire: return "fire";
    case MinMethod::Lbfgs: return "lbfgs";
    case MinMethod::Adadelta: return "adadelta";
  }
  return "?";
}

namespace {

constexpr double kArmijo = 1e-4;
constexpr int kMaxBacktrack = 20;

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// State shared by the methods: the current point, its energy and gradient, the best point
// and the stop rules.
class Run {
 public:
  Run(const VecObjective& f, std::vector<double>& x, const VecMinParams& p) : f_(f), x_(x), p_(p) {
    g_.resize(x.size());
  }
  bool start() {
    if (!eval(x_, e_, g_)) return false;
    best_ = x_;
    bestE_ = e_;
    return true;
  }
  bool eval(std::span<const double> x, double& e, std::vector<double>& g) {
    ++r.evaluations;
    if (f_(x, e, g)) return true;
    r.aborted = true;
    return false;
  }
  double limit(std::span<const double> step) const { return p_.stepLimit ? std::min(1.0, p_.stepLimit(step)) : 1.0; }
  bool more() const { return r.steps < p_.maxSteps; }
  // Moves to (x, e, g); returns false when the energy-change rule stops the minimization.
  bool accept(std::vector<double>& x, double e, std::vector<double>& g) {
    ++r.steps;
    const double dE = std::abs(e - e_);
    std::swap(x_, x);
    std::swap(g_, g);
    e_ = e;
    const bool best = e_ < bestE_;
    if (p_.onAccept) p_.onAccept(x_, best);
    if (best) {
      best_ = x_;
      bestE_ = e_;
    }
    calm_ = dE < p_.deltaE ? calm_ + 1 : 0;
    return calm_ < p_.convSteps;
  }
  VecMinResult finish() {
    if (!r.aborted) {
      x_ = best_;
      r.energy = bestE_;
    }
    return r;
  }
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& g() const { return g_; }
  double e() const { return e_; }
  std::size_t n() const { return x_.size(); }
  const VecMinParams& p() const { return p_; }

  VecMinResult r;

 private:
  const VecObjective& f_;
  std::vector<double>& x_;
  const VecMinParams& p_;
  std::vector<double> g_, best_;
  double e_ = 0, bestE_ = 0;
  int calm_ = 0;
};

enum class Search { Accepted, Failed, Aborted };

// Armijo backtracking along d (already step-limited) from the current point; slope = g.d < 0.
Search backtrack(Run& run, const std::vector<double>& d, double slope, std::vector<double>& xNew, double& eNew,
                 std::vector<double>& gNew, double* used = nullptr) {
  double a = 1;
  for (int b = 0; b < kMaxBacktrack; ++b, a *= 0.5) {
    for (std::size_t i = 0; i < run.n(); ++i) xNew[i] = run.x()[i] + a * d[i];
    if (!run.eval(xNew, eNew, gNew)) return Search::Aborted;
    if (eNew <= run.e() + kArmijo * a * slope) {
      if (used) *used = a;
      return Search::Accepted;
    }
  }
  return Search::Failed;
}

// Steepest descent; the step factor grows by 1.2 after an accepted step and halves after a rejected one.
void runSd(Run& run) {
  const std::size_t n = run.n();
  std::vector<double> d(n), xNew(n), gNew(n);
  double a = 1, eNew;
  int fails = 0;
  while (run.more()) {
    for (std::size_t i = 0; i < n; ++i) d[i] = -a * run.g()[i];
    const double s = run.limit(d);
    a *= s;  // never larger than the step limit allows
    double slope = 0;
    for (std::size_t i = 0; i < n; ++i) {
      d[i] *= s;
      xNew[i] = run.x()[i] + d[i];
      slope += run.g()[i] * d[i];
    }
    if (!(slope < 0)) return;  // zero gradient
    if (!run.eval(xNew, eNew, gNew)) return;
    if (eNew <= run.e() + kArmijo * slope) {
      fails = 0;
      a *= 1.2;
      if (!run.accept(xNew, eNew, gNew)) return;
    } else {
      a *= 0.5;
      if (++fails >= kMaxBacktrack) return;
    }
  }
}

// Polak-Ribiere (PR+) conjugate gradient with Armijo backtracking. The trial step starts at
// alpha_prev (g_prev.d_prev) / (g.d) (Nocedal & Wright 3.60).
void runCg(Run& run) {
  const std::size_t n = run.n();
  std::vector<double> d(n), dPrev, gPrev, step(n), xNew(n), gNew(n);
  double alpha = 1, prevSlope = 0, eNew;
  while (run.more()) {
    double slope;
    if (dPrev.empty()) {
      for (std::size_t i = 0; i < n; ++i) d[i] = -run.g()[i];
    } else {
      double num = 0;
      for (std::size_t i = 0; i < n; ++i) num += run.g()[i] * (run.g()[i] - gPrev[i]);
      const double beta = std::max(0.0, num / dot(gPrev, gPrev));
      for (std::size_t i = 0; i < n; ++i) d[i] = -run.g()[i] + beta * dPrev[i];
    }
    slope = dot(run.g(), d);
    if (!(slope < 0)) {  // not a descent direction: restart
      for (std::size_t i = 0; i < n; ++i) d[i] = -run.g()[i];
      slope = dot(run.g(), d);
      if (!(slope < 0)) return;
    }
    const double a0 = dPrev.empty() ? 1.0 : alpha * prevSlope / slope;
    for (std::size_t i = 0; i < n; ++i) step[i] = a0 * d[i];
    const double s = run.limit(step);
    for (double& v : step) v *= s;
    double used = 1;
    const Search sr = backtrack(run, step, a0 * s * slope, xNew, eNew, gNew, &used);
    if (sr != Search::Accepted) return;
    alpha = a0 * s * used;
    prevSlope = slope;
    dPrev = d;
    gPrev = run.g();
    if (!run.accept(xNew, eNew, gNew)) return;
  }
}

// FIRE (Bitzek et al. 2006): N_min 5, f_inc 1.1, f_dec 0.5, alpha0 0.1, f_alpha 0.99, dt0 0.1,
// dt_max 1.0, unit mass. Every step is accepted.
void runFire(Run& run) {
  constexpr int kNmin = 5;
  constexpr double kFinc = 1.1, kFdec = 0.5, kAlpha0 = 0.1, kFalpha = 0.99, kDt0 = 0.1, kDtMax = 1.0;
  const std::size_t n = run.n();
  std::vector<double> v(n, 0.0), step(n), xNew(n), gNew(n);
  double dt = kDt0, alpha = kAlpha0, eNew;
  int npos = 0;
  while (run.more()) {
    const std::vector<double>& g = run.g();
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
    for (std::size_t i = 0; i < n; ++i) {
      v[i] -= dt * g[i];
      step[i] = dt * v[i];
    }
    const double s = run.limit(step);
    for (std::size_t i = 0; i < n; ++i) xNew[i] = run.x()[i] + s * step[i];
    if (!run.eval(xNew, eNew, gNew)) return;
    if (!run.accept(xNew, eNew, gNew)) return;
  }
}

// L-BFGS: two-loop recursion, initial Hessian (s.y / y.y) I, steps with s.y <= 1e-12 not stored.
void runLbfgs(Run& run) {
  const std::size_t n = run.n();
  std::vector<double> d(n), xNew(n), gNew(n), xOld, gOld;
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  double eNew;
  while (run.more()) {
    d = run.g();
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
    double slope = dot(run.g(), d);
    if (!(slope < 0)) {  // not a descent direction: restart from steepest descent
      S.clear();
      Y.clear();
      rho.clear();
      for (std::size_t i = 0; i < n; ++i) d[i] = -run.g()[i];
      slope = dot(run.g(), d);
      if (!(slope < 0)) return;  // zero gradient
    }
    const double s = run.limit(d);
    for (double& v : d) v *= s;
    slope *= s;
    if (backtrack(run, d, slope, xNew, eNew, gNew) != Search::Accepted) return;  // no decrease: converged
    std::vector<double> sv(n), yv(n);
    for (std::size_t i = 0; i < n; ++i) {
      sv[i] = xNew[i] - run.x()[i];
      yv[i] = gNew[i] - run.g()[i];
    }
    const double sy = dot(sv, yv);
    if (sy > 1e-12) {
      if (static_cast<int>(S.size()) == run.p().memory) {
        S.erase(S.begin());
        Y.erase(Y.begin());
        rho.erase(rho.begin());
      }
      S.push_back(std::move(sv));
      Y.push_back(std::move(yv));
      rho.push_back(1 / sy);
    }
    if (!run.accept(xNew, eNew, gNew)) return;
  }
}

// ADADELTA (Zeiler 2012) with rho 0.8 and epsilon 1e-2 as in AutoDock-GPU. Every step is accepted.
void runAdadelta(Run& run) {
  constexpr double kRho = 0.8, kEps = 1e-2;
  const std::size_t n = run.n();
  std::vector<double> eg2(n, 0.0), ex2(n, 0.0), step(n), xNew(n), gNew(n);
  double eNew;
  while (run.more()) {
    const std::vector<double>& g = run.g();
    for (std::size_t i = 0; i < n; ++i) {
      eg2[i] = kRho * eg2[i] + (1 - kRho) * g[i] * g[i];
      step[i] = -g[i] * std::sqrt((ex2[i] + kEps) / (eg2[i] + kEps));
    }
    const double s = run.limit(step);
    for (std::size_t i = 0; i < n; ++i) {
      step[i] *= s;
      ex2[i] = kRho * ex2[i] + (1 - kRho) * step[i] * step[i];
      xNew[i] = run.x()[i] + step[i];
    }
    if (!run.eval(xNew, eNew, gNew)) return;
    if (!run.accept(xNew, eNew, gNew)) return;
  }
}

}  // namespace

VecMinResult minimizeVec(const VecObjective& f, std::vector<double>& x, const VecMinParams& p) {
  Run run(f, x, p);
  if (!run.start()) return run.r;
  switch (p.method) {
    case MinMethod::Sd: runSd(run); break;
    case MinMethod::Cg: runCg(run); break;
    case MinMethod::Fire: runFire(run); break;
    case MinMethod::Lbfgs: runLbfgs(run); break;
    case MinMethod::Adadelta: runAdadelta(run); break;
  }
  return run.finish();
}

MinResult minimizeTorsions(const TorsionObjective& f, std::span<double> anglesDeg, const MinParams& p) {
  constexpr double kRad = std::numbers::pi / 180.0;
  // Variables in rad; adadelta in degrees (AutoDock-GPU units).
  const double u = p.method == MinMethod::Adadelta ? kRad : 1.0;
  const std::size_t n = anglesDeg.size();
  std::vector<double> x(n), deg(n);
  for (std::size_t i = 0; i < n; ++i) x[i] = anglesDeg[i] * kRad / u;
  const VecObjective vf = [&](std::span<const double> v, double& e, std::vector<double>& g) {
    for (std::size_t i = 0; i < n; ++i) deg[i] = v[i] * u / kRad;
    e = f(deg, g);
    for (double& gi : g) gi *= u;
    return true;
  };
  VecMinParams vp;
  vp.method = p.method;
  vp.maxSteps = p.maxSteps;
  vp.convSteps = 1;
  vp.deltaE = p.deltaE;
  const double maxStep = p.maxStepDeg * kRad / u;
  vp.stepLimit = [maxStep](std::span<const double> d) {
    double big = 0;
    for (double v : d) big = std::max(big, std::abs(v));
    return big > maxStep ? maxStep / big : 1.0;
  };
  const VecMinResult r = minimizeVec(vf, x, vp);
  for (std::size_t i = 0; i < n; ++i) anglesDeg[i] = wrapDeg(x[i] * u / kRad);
  return {r.energy, r.steps, r.evaluations};
}

}  // namespace od
