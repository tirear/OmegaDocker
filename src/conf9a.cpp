// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "conf9a.hpp"

#include <algorithm>
#include <cmath>

#include "log.hpp"

namespace od {

EnergyTerms minimizeConf(const EnergySystem& sys, std::span<double> a, const MinParams& p, int* steps) {
  const TorsionObjective f = [&](std::span<const double> x, std::vector<double>& g) {
    return sys.evaluate(x, &g).total();
  };
  const MinResult r = minimizeFire(f, a, p);
  if (steps) *steps = r.steps;
  return sys.evaluate(a);
}

Fill9aStats fillL2Pool(const EnergySystem& sys, const ConfGenerator& gen, long long l1, const Fill9aParams& p,
                       std::uint64_t& counter, std::vector<ConfSample>& pool) {
  Fill9aStats st;
  const int d = sys.torsionCount();
  if (l1 < 1) l1 = 1;
  while (static_cast<long long>(pool.size()) < p.minL2) {
    long long added = 0;
    for (long long done = 0; done < l1;) {
      const long long n = std::min(l1 - done, p.maxL1Batch);
      std::vector<double> angles(static_cast<std::size_t>(n) * d);
      std::vector<EnergyTerms> terms(n);
      std::vector<int> steps(n);
#pragma omp parallel for schedule(dynamic, 1)
      for (long long i = 0; i < n; ++i) {
        Rng rng(p.seed, p.stream, counter + done + i);
        std::span<double> a(angles.data() + i * d, d);
        gen(rng, a);
        terms[i] = minimizeConf(sys, a, p.min, &steps[i]);
      }
      for (long long i = 0; i < n; ++i) {
        st.minSteps += steps[i];
        const EnergyTerms& e = terms[i];
        if (!(e.filterEnergy() < p.Emax)) {
          ++st.overEmax;
        } else if (e.closureMaxDeviation > p.maxBreakDist) {
          ++st.notClosed;
        } else {
          ConfSample s;
          s.angles.assign(angles.begin() + i * d, angles.begin() + (i + 1) * d);
          s.energy = static_cast<float>(e.filterEnergy());
          pool.push_back(std::move(s));
          ++added;
        }
      }
      ++st.batches;
      done += n;
    }
    counter += l1;
    st.l1 += l1;
    st.added += added;
    const long long nowL2 = static_cast<long long>(pool.size());
    long long next = added > 0 ? 2 * static_cast<long long>(std::ceil(static_cast<double>(p.minL2 - nowL2) * l1 / added))
                               : l1 * 1024;
    next = std::clamp(next, 1LL, p.maxL1Batch);
    st.nextL1 = next;
    if (nowL2 >= p.minL2) break;
    if (st.l1 >= p.maxL1Total)
      fail("9a: {} samples generated (--maxL1TotalNum) but only {} of {} passed the filters "
           "(energy >= --Emax: {}, closure deviation > --maxBreakDist: {})",
           st.l1, nowL2, p.minL2, st.overEmax, st.notClosed);
    l1 = std::min(next, p.maxL1Total - st.l1);
  }
  return st;
}

}  // namespace od
