// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "csgm3.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "log.hpp"

namespace od {

namespace {

double norm(const std::vector<double>& x) {
  double s = 0;
  for (double v : x) s += v * v;
  return std::sqrt(s);
}

// Offset perpendicular to u (unit, or zero) with an integer length in [0, maxOffset].
std::vector<double> offset(Rng& rng, const std::vector<double>& u, int maxOffset) {
  const std::size_t d = u.size();
  const double len = static_cast<double>(rng.below(static_cast<std::uint64_t>(maxOffset) + 1));
  std::vector<double> v(d, 0.0);
  for (int attempt = 0; attempt < 64; ++attempt) {
    for (std::size_t i = 0; i < d; i += 2) {  // Box-Muller: normal components, uniform direction
      const double r = std::sqrt(-2.0 * std::log(1.0 - rng.uniform()));
      const double t = 2.0 * std::numbers::pi * rng.uniform();
      v[i] = r * std::cos(t);
      if (i + 1 < d) v[i + 1] = r * std::sin(t);
    }
    double dot = 0;
    for (std::size_t i = 0; i < d; ++i) dot += v[i] * u[i];
    for (std::size_t i = 0; i < d; ++i) v[i] -= dot * u[i];
    const double n = norm(v);
    if (n > 1e-9) {
      for (double& x : v) x *= len / n;
      return v;
    }
  }
  std::fill(v.begin(), v.end(), 0.0);  // no perpendicular direction (one torsion)
  return v;
}

}  // namespace

long long csgm3SampleNum(std::size_t k, const Csgm3Params& p) {
  return static_cast<long long>(k) * static_cast<long long>(k - 1) * p.perLine;
}

void csgm3Sample(Rng rng, std::span<const std::vector<float>* const> ends, long long line, int m,
                 const std::vector<bool>& cisTrans, const Csgm3Params& p, std::span<double> a) {
  const long long k = static_cast<long long>(ends.size());
  if (k < 2) fail("csgm3: at least 2 end points are needed");
  if (line < 0 || line >= k * (k - 1) || m < 0 || m >= p.perLine) fail("csgm3: sample index out of range");
  const std::size_t d = cisTrans.size();
  const long long i = line / (k - 1);
  long long j = line % (k - 1);
  if (j >= i) ++j;
  const std::vector<float>& pi = *ends[i];
  const std::vector<float>& pj = *ends[j];
  if (pi.size() != a.size() || pj.size() != a.size() || d > a.size()) fail("csgm3: end point size mismatch");
  for (std::size_t t = d; t < a.size(); ++t) a[t] = 2 * m < p.perLine ? pi[t] : pj[t];
  std::vector<double> u(d);
  for (std::size_t t = 0; t < d; ++t) u[t] = wrapDeg(static_cast<double>(pj[t]) - pi[t]);
  const double un = norm(u);
  for (double& x : u) x = un > 0 ? x / un : 0.0;
  const std::vector<double> v1 = offset(rng, u, p.maxOffset);
  const std::vector<double> v2 = offset(rng, u, p.maxOffset);
  std::vector<double> A(d), w(d);
  for (std::size_t t = 0; t < d; ++t) {
    A[t] = pi[t] + v1[t];
    w[t] = wrapDeg(pj[t] + v2[t] - A[t]);
  }
  const bool fromA = rng.uniform() < 0.5;
  const double wn = norm(w);
  double wmax = 0;
  for (double& x : w) {
    x = wn > 0 ? x / wn : 0.0;
    wmax = std::max(wmax, std::abs(x));
  }
  const double step = wmax > 0 ? p.coverAngle / wmax / p.perLine : 0.0;
  for (std::size_t t = 0; t < d; ++t) {
    const double s = fromA ? A[t] : pj[t] + v2[t];
    a[t] = wrapDeg(s + m * step * w[t]);
    if (cisTrans[t]) a[t] = std::abs(a[t]) < 90.0 ? 0.0 : 180.0;
  }
}

}  // namespace od
