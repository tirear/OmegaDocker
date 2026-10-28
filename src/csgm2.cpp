// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "csgm2.hpp"

#include <cmath>
#include <utility>

#include "log.hpp"

namespace od {

void csgm2Offspring(Rng& rng, std::span<const Csgm2Parent> parents, const std::vector<bool>& cisTrans,
                    const Csgm2Params& p, std::span<double> a, std::span<const int> ringK) {
  if (parents.empty()) fail("csgm2: no parents");
  const std::size_t k = parents.size();
  const int n = static_cast<int>(a.size());
  const int nt = static_cast<int>(cisTrans.size());
  if (nt + static_cast<int>(ringK.size()) != n) fail("csgm2: gene count mismatch");
  // Two binary tournaments; the four candidates need not be different.
  int cand[4];
  for (int& c : cand) c = static_cast<int>(rng.below(k));
  int par[2];
  for (int t = 0; t < 2; ++t) {
    const int x = cand[2 * t], y = cand[2 * t + 1];
    const bool xFitter = parents[x].fitness < parents[y].fitness;
    const bool fitterWins = rng.uniform() < p.tournamentRate;
    par[t] = xFitter == fitterWins ? x : y;
  }
  const std::vector<float>& g0 = *parents[par[0]].angles;
  const std::vector<float>& g1 = *parents[par[1]].angles;
  if (static_cast<int>(g0.size()) != n || static_cast<int>(g1.size()) != n) fail("csgm2: parent size mismatch");
  // Two-point crossover; single-point when both points coincide.
  if (rng.uniform() < p.crossoverRate) {
    int c0 = static_cast<int>((n - 1) * rng.uniform());
    int c1 = static_cast<int>((n - 1) * rng.uniform());
    if (c1 < c0) std::swap(c0, c1);
    for (int i = 0; i < n; ++i) {
      const bool first = c0 != c1 ? (i <= c0 || i > c1) : i <= c0;
      a[i] = first ? g0[i] : g1[i];
    }
  } else {
    for (int i = 0; i < n; ++i) a[i] = g0[i];
  }
  for (int i = 0; i < n; ++i) {
    if (rng.uniform() >= p.mutationRate) continue;
    if (i >= nt) a[i] = static_cast<double>(rng.below(static_cast<std::uint64_t>(ringK[i - nt])));
    else if (cisTrans[i]) a[i] = std::abs(wrapDeg(a[i])) < 90.0 ? 180.0 : 0.0;
    else a[i] = wrapDeg(a[i] + p.absMaxDang * (2.0 * rng.uniform() - 1.0));
  }
}

}  // namespace od
