// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "dock9c.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

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
template <class Dist>
std::vector<RoundSample> greedyFilter(std::vector<RoundSample>& samples, const MobmHeavyAtoms& atoms, std::size_t topNum,
                                      double rmsd, const Dist& dist) {
  std::vector<RoundSample> kept;
  std::vector<std::vector<Vec3f>> keptXyz;
  std::vector<Vec3f> keptC;
  std::unordered_map<std::uint64_t, std::vector<int>> cells;
  auto cellOf = [&](const Vec3f& c, int d) { return static_cast<std::int64_t>(std::floor(c[d] / rmsd)); };
  auto key = [](int stam, std::int64_t x, std::int64_t y, std::int64_t z) {
    std::uint64_t h = static_cast<std::uint64_t>(stam) * 0x9E3779B97F4A7C15ULL;
    for (std::int64_t v : {x, y, z}) h = (h ^ static_cast<std::uint64_t>(v)) * 0xBF58476D1CE4E5B9ULL + 0x94D049BB133111EBULL;
    return h;
  };
  for (auto& m : samples) {
    if (kept.size() >= topNum) break;
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
            if (kept[k].stam != m.stam) continue;
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
    kept.push_back(std::move(m));
  }
  return kept;
}

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
  std::stable_sort(merged.begin(), merged.end(), [](const auto& a, const auto& b) { return a.s.energy < b.s.energy; });
  if (!(diversityRmsd > 0)) {
    if (merged.size() > static_cast<std::size_t>(topNum)) merged.resize(topNum);
    return merged;
  }
  return greedyFilter(merged, atoms, static_cast<std::size_t>(topNum), diversityRmsd,
                      [](std::span<const Vec3f> a, std::span<const Vec3f> b) { return fixedRmsd(a, b); });
}

std::vector<RoundSample> symmetryFilter(std::vector<RoundSample> top, const MobmHeavyAtoms& atoms,
                                        const SymmetryRmsd& sym, double rmsd) {
  if (!(rmsd > 0)) return top;
  return greedyFilter(top, atoms, top.size(), rmsd,
                      [&](std::span<const Vec3f> a, std::span<const Vec3f> b) { return sym.rmsd(a, b); });
}

int filterDcss(std::vector<int>& dcss, const std::unordered_map<int, float>& lowest, float eTopMax) {
  const auto n = dcss.size();
  std::erase_if(dcss, [&](int id) {
    const auto it = lowest.find(id);
    return it == lowest.end() || it->second > eTopMax;
  });
  return static_cast<int>(n - dcss.size());
}

}  // namespace od
