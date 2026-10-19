// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "conf9a.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "fiber_batch.hpp"
#include "log.hpp"

namespace od {

EnergyTerms minimizeConf(const EnergySystem& sys, std::span<double> a, const MinParams& p, int* steps) {
  const TorsionObjective f = [&](std::span<const double> x, std::vector<double>& g) {
    return sys.evaluate(x, &g).total();
  };
  const MinResult r = minimizeTorsions(f, a, p);
  if (steps) *steps = r.steps;
  return sys.evaluate(a);
}

void minimizeConfBatch(ConfEvaluator& ev, int d, std::span<double> a, const MinParams& p, EnergyTerms* terms, int* steps) {
  const int count = d > 0 ? static_cast<int>(a.size() / d) : 0;
  if (count == 0) return;
  const int batch = std::max(1, ev.maxBatch());
  std::vector<MinResult> r;
  if (ev.lanes() > 1) {
    std::vector<std::vector<EnergyTerms>> buf(ev.lanes());
    std::vector<int> size(ev.lanes(), 0);
    BatchTorsionAsync f;
    f.lanes = ev.lanes();
    f.submit = [&](int l, const double* x, int n) {
      size[l] = n;
      ev.submit(l, x, n);
    };
    f.complete = [&](int l, double* e, double* g) {
      buf[l].resize(size[l]);
      ev.complete(l, buf[l].data(), g);
      for (int i = 0; i < size[l]; ++i) e[i] = buf[l][i].total();
    };
    r = minimizeTorsionsBatch(f, a, d, p, batch);
  } else {
    std::vector<EnergyTerms> buf;
    const BatchTorsionObjective f = [&](const double* x, int n, double* e, double* g) {
      buf.resize(n);
      ev.evaluate(x, n, buf.data(), g);
      for (int i = 0; i < n; ++i) e[i] = buf[i].total();
    };
    r = minimizeTorsionsBatch(f, a, d, p, batch);
  }
  for (int i = 0; i < count; ++i) steps[i] = r[i].steps;
  for (int done = 0; done < count; done += batch)
    ev.evaluate(a.data() + static_cast<std::size_t>(done) * d, std::min(batch, count - done), terms + done, nullptr);
}

Fill9aStats fillL2Pool(const EnergySystem& sys, const ConfGenerator& gen, long long l1, const Fill9aParams& p,
                       std::uint64_t& counter, std::vector<ConfSample>& pool) {
  Fill9aStats st;
  const int d = sys.torsionCount();
  if (p.fixedL1 > 0) l1 = p.fixedL1;
  if (l1 < 1) l1 = 1;
  const double preClose = p.preCloseFactor > 0 && sys.hasClosures() ? p.preCloseFactor * p.maxBreakDist : 0.0;
  while (static_cast<long long>(pool.size()) < p.minL2) {
    long long added = 0;
    for (long long done = 0; done < l1;) {
      const long long n = std::min(l1 - done, p.maxL1Batch);
      std::vector<double> angles(static_cast<std::size_t>(n) * d);
      std::vector<EnergyTerms> terms(n);
      std::vector<int> steps(n);
      std::vector<GenResult> built(n);
      const auto t0 = std::chrono::steady_clock::now();
      const auto one = [&](long long i, int slot) {
        Rng rng(p.seed, p.stream, counter + done + i);
        std::span<double> a(angles.data() + i * d, d);
        built[i] = gen(rng, GenIndex{counter + done + i, done + i, slot}, a);
        if (built[i] == GenResult::Built && preClose > 0 && sys.closureDeviation(a) > preClose)
          built[i] = GenResult::PreNotClosed;
        if (built[i] == GenResult::Built && !p.evaluator) terms[i] = minimizeConf(sys, a, p.min, &steps[i]);
      };
      if (p.genBatch && n >= p.genBatchMin) {
        runBatched(static_cast<int>(n), std::max(1, std::min(p.genBatch->slots(), static_cast<int>(n))),
                   [&](int i, int slot) { one(i, slot); }, [&](std::span<const int> s) { p.genBatch->flush(s); });
      } else {
#pragma omp parallel for schedule(dynamic, 1)
        for (long long i = 0; i < n; ++i) one(i, -1);
      }
      const auto t1 = std::chrono::steady_clock::now();
      st.genSec += std::chrono::duration<double>(t1 - t0).count();
      if (p.evaluator) {
        std::vector<long long> idx;
        for (long long i = 0; i < n; ++i)
          if (built[i] == GenResult::Built) idx.push_back(i);
        std::vector<double> a(idx.size() * d);
        for (std::size_t k = 0; k < idx.size(); ++k)
          std::copy_n(angles.begin() + idx[k] * d, d, a.begin() + static_cast<std::ptrdiff_t>(k * d));
        std::vector<EnergyTerms> t(idx.size());
        std::vector<int> s(idx.size());
        minimizeConfBatch(*p.evaluator, d, a, p.min, t.data(), s.data());
        for (std::size_t k = 0; k < idx.size(); ++k) {
          std::copy_n(a.begin() + static_cast<std::ptrdiff_t>(k * d), d, angles.begin() + idx[k] * d);
          terms[idx[k]] = t[k];
          steps[idx[k]] = s[k];
        }
        st.minSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
      }
      for (long long i = 0; i < n; ++i) {
        const EnergyTerms& e = terms[i];
        if (built[i] == GenResult::GrowFailed) {
          ++st.growFailed;
          continue;
        }
        if (built[i] == GenResult::CloseFailed) {
          ++st.closeFailed;
          continue;
        }
        if (built[i] == GenResult::SideGrowFailed) {
          ++st.sideGrowFailed;
          continue;
        }
        if (built[i] == GenResult::PreNotClosed) {
          ++st.preNotClosed;
          continue;
        }
        st.minSteps += steps[i];
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
      fail("9a: {} samples generated (stam.maxL1TotalNum / mobm.maxL1TotalNum) but only {} of {} passed the filters "
           "(discarded while growing: {}, loop closure failed: {}, discarded while growing after the closure: {}, "
           "closure deviation > conf9a.preCloseFactor x topology.maxBreakDist before minimization: {}, "
           "energy >= energy.Emax: {}, closure deviation > topology.maxBreakDist: {})",
           st.l1, nowL2, p.minL2, st.growFailed, st.closeFailed, st.sideGrowFailed, st.preNotClosed, st.overEmax, st.notClosed);
    l1 = std::min(p.fixedL1 > 0 ? p.fixedL1 : next, p.maxL1Total - st.l1);
  }
  if (p.evaluator) p.evaluator->idle();
  return st;
}

}  // namespace od
