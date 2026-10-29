// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "dock9c.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_set>

#include "ad4.hpp"
#include "orient.hpp"

namespace od {

namespace {

Vec3f centroid(std::span<const Vec3f> x) {
  double c[3] = {0, 0, 0};
  for (const auto& p : x)
    for (int d = 0; d < 3; ++d) c[d] += p[d];
  const double n = static_cast<double>(std::max<std::size_t>(x.size(), 1));
  return {static_cast<float>(c[0] / n), static_cast<float>(c[1] / n), static_cast<float>(c[2] / n)};
}

// Greedy filter over samples in order: a sample is skipped when dist() < rmsd to a
// kept sample with the same StaM conformation. Valid only for one-to-one pairings,
// for which RMSD >= centroid distance: kept samples are found through a hash of
// their centroids (cell = rmsd) and only those within rmsd are compared.
// Returns the indices of the kept samples (at most topNum).
template <class Dist>
std::vector<std::size_t> greedyKept(const std::vector<RoundSample>& samples, const MobmHeavyAtoms& atoms,
                                    std::size_t topNum, double rmsd, const Dist& dist) {
  std::vector<std::size_t> kept;
  std::vector<std::vector<Vec3f>> keptXyz;
  std::vector<Vec3f> keptC;
  std::unordered_map<std::uint64_t, std::vector<int>> cells;
  auto cellOf = [&](const Vec3f& c, int d) { return static_cast<std::int64_t>(std::floor(c[d] / rmsd)); };
  auto key = [](int stam, std::int64_t x, std::int64_t y, std::int64_t z) {
    std::uint64_t h = static_cast<std::uint64_t>(stam) * 0x9E3779B97F4A7C15ULL;
    for (std::int64_t v : {x, y, z}) h = (h ^ static_cast<std::uint64_t>(v)) * 0xBF58476D1CE4E5B9ULL + 0x94D049BB133111EBULL;
    return h;
  };
  for (std::size_t idx = 0; idx < samples.size(); ++idx) {
    if (kept.size() >= topNum) break;
    const RoundSample& m = samples[idx];
    auto xyz = atoms.world(m.s);
    const Vec3f c = centroid(xyz);
    const std::int64_t cx = cellOf(c, 0), cy = cellOf(c, 1), cz = cellOf(c, 2);
    bool dup = false;
    for (std::int64_t dx = -1; dx <= 1 && !dup; ++dx)
      for (std::int64_t dy = -1; dy <= 1 && !dup; ++dy)
        for (std::int64_t dz = -1; dz <= 1 && !dup; ++dz) {
          const auto it = cells.find(key(m.stam, cx + dx, cy + dy, cz + dz));
          if (it == cells.end()) continue;
          for (int k : it->second) {
            if (samples[kept[k]].stam != m.stam) continue;
            const double ex = c[0] - keptC[k][0], ey = c[1] - keptC[k][1], ez = c[2] - keptC[k][2];
            if (ex * ex + ey * ey + ez * ez >= rmsd * rmsd) continue;
            if (dist(xyz, keptXyz[k]) < rmsd) {
              dup = true;
              break;
            }
          }
        }
    if (dup) continue;
    cells[key(m.stam, cx, cy, cz)].push_back(static_cast<int>(kept.size()));
    keptXyz.push_back(std::move(xyz));
    keptC.push_back(c);
    kept.push_back(idx);
  }
  return kept;
}

std::vector<RoundSample> take(std::vector<RoundSample>& samples, const std::vector<std::size_t>& idx) {
  std::vector<RoundSample> out;
  out.reserve(idx.size());
  for (std::size_t i : idx) out.push_back(std::move(samples[i]));
  return out;
}

void sortByEnergy(std::vector<RoundSample>& v) {
  std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.rank() < b.rank(); });
}

double fixedDist(std::span<const Vec3f> a, std::span<const Vec3f> b) { return fixedRmsd(a, b); }

}  // namespace

MobmHeavyAtoms::MobmHeavyAtoms(const EnergySystem& mobSys) : sys_(mobSys) {
  const Ad4Typing& t = mobSys.typing(0);
  const int off = mobSys.coordOffset(0);
  for (std::size_t a = 0; a < t.type.size(); ++a) {
    const Ad4Type y = t.type[a];
    if (y == Ad4Type::Merged || y == Ad4Type::H || y == Ad4Type::HD || y == Ad4Type::HS) continue;
    index_.push_back(off + static_cast<int>(a));
    atom_.push_back(static_cast<int>(a));
    type_.push_back(static_cast<int>(y));
  }
}

std::vector<Vec3f> MobmHeavyAtoms::world(const DockSample& s) const {
  const auto x = sys_.coordinates(s.angles);
  const QuatF q = toFloat(s.pose.q);
  const Vec3f c{static_cast<float>(s.center[0]), static_cast<float>(s.center[1]), static_cast<float>(s.center[2])};
  const Vec3f t{static_cast<float>(s.pose.t[0]), static_cast<float>(s.pose.t[1]), static_cast<float>(s.pose.t[2])};
  std::vector<Vec3f> w;
  w.reserve(index_.size());
  for (int i : index_) {
    const Vec3f r = rotate(q, {x[i][0] - c[0], x[i][1] - c[1], x[i][2] - c[2]});
    w.push_back({r[0] + t[0], r[1] + t[1], r[2] + t[2]});
  }
  return w;
}

std::vector<RoundSample> selectTopSamples(std::vector<RoundSample> merged, const MobmHeavyAtoms& atoms, int topNum,
                                          double diversityRmsd) {
  sortByEnergy(merged);
  if (!(diversityRmsd > 0)) {
    if (merged.size() > static_cast<std::size_t>(topNum)) merged.resize(topNum);
    return merged;
  }
  return take(merged, greedyKept(merged, atoms, static_cast<std::size_t>(topNum), diversityRmsd, fixedDist));
}

PoolSelection selectPool(std::vector<RoundSample> merged, const MobmHeavyAtoms& atoms, std::size_t keepNum,
                         double diversityRmsd, std::size_t cap) {
  sortByEnergy(merged);
  const std::size_t n = merged.size();
  std::size_t end = n;
  std::vector<std::size_t> kept;
  if (!(diversityRmsd > 0)) {
    end = std::min(keepNum, n);
    for (std::size_t k = 0; k < end; ++k) kept.push_back(k);
  } else {
    kept = greedyKept(merged, atoms, keepNum, diversityRmsd, fixedDist);
    if (kept.size() >= keepNum) end = kept.back() + 1;
  }
  while (end > 0 && end < n && merged[end].rank() == merged[end - 1].rank()) ++end;
  PoolSelection r;
  if (end < n) r.droppedMin = merged[end].rank();
  if (end <= cap) {
    merged.resize(end);
    r.pool = std::move(merged);
    return r;
  }
  // Over the cap: all kept samples, then the lowest-energy skipped samples.
  std::vector<char> isKept(end, 0);
  for (std::size_t k : kept) isKept[k] = 1;
  std::size_t budget = cap > kept.size() ? cap - kept.size() : 0;
  r.pool.reserve(kept.size() + budget);
  for (std::size_t k = 0; k < end; ++k) {
    if (isKept[k]) {
      r.pool.push_back(std::move(merged[k]));
    } else if (budget > 0) {
      r.pool.push_back(std::move(merged[k]));
      --budget;
    } else {
      r.droppedMin = std::min(r.droppedMin, merged[k].rank());
    }
  }
  return r;
}

std::vector<RoundSample> symmetryFilter(std::vector<RoundSample> top, const MobmHeavyAtoms& atoms,
                                        const SymmetryRmsd& sym, double rmsd) {
  if (!(rmsd > 0)) return top;
  return take(top, greedyKept(top, atoms, top.size(), rmsd,
                              [&](std::span<const Vec3f> a, std::span<const Vec3f> b) { return sym.rmsd(a, b); }));
}

DcssGrowth growDcssSizes(const std::vector<int>& sizes, const std::vector<bool>& grow,
                         const std::function<std::size_t(const std::vector<int>&)>& estimate, std::size_t budget) {
  const auto scaled = [&](double t) {
    std::vector<int> s = sizes;
    for (std::size_t k = 0; k < s.size(); ++k)
      if (grow[k])
        s[k] = static_cast<int>(std::min(std::floor(static_cast<double>(sizes[k]) * t), static_cast<double>(std::numeric_limits<int>::max() / 4)));
    return s;
  };
  const auto fits = [&](const std::vector<int>& s) { return budget == 0 || estimate(s) <= budget; };
  DcssGrowth g;
  g.sizes = scaled(2.0);
  if (!fits(g.sizes)) {
    double lo = 1.0, hi = 2.0;
    for (int it = 0; it < 60; ++it) {
      const double mid = 0.5 * (lo + hi);
      (fits(scaled(mid)) ? lo : hi) = mid;
    }
    g.sizes = scaled(lo);
    if (!fits(g.sizes)) g.sizes = sizes;
  }
  g.estimate = estimate(g.sizes);
  for (std::size_t k = 0; k < sizes.size(); ++k)
    if (grow[k] && g.sizes[k] <= sizes[k]) g.atLimit = true;
  return g;
}

int filterDcss(std::vector<int>& dcss, const std::vector<RoundSample>& top, int RoundSample::* id) {
  const auto n = dcss.size();
  std::unordered_set<int> inTop;
  for (const auto& t : top) inTop.insert(t.*id);
  std::erase_if(dcss, [&](int s) { return !inTop.contains(s); });
  return static_cast<int>(n - dcss.size());
}

}  // namespace od
