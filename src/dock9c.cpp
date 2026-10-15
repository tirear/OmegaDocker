// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "dock9c.hpp"

#include <algorithm>
#include <cmath>

#include "ad4.hpp"
#include "orient.hpp"

namespace od {

MobmHeavyAtoms::MobmHeavyAtoms(const EnergySystem& mobSys) : sys_(mobSys) {
  const Ad4Typing& t = mobSys.typing(0);
  const int off = mobSys.coordOffset(0);
  for (std::size_t a = 0; a < t.type.size(); ++a) {
    const Ad4Type y = t.type[a];
    if (y == Ad4Type::Merged || y == Ad4Type::H || y == Ad4Type::HD || y == Ad4Type::HS) continue;
    index_.push_back(off + static_cast<int>(a));
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

double symmetricRmsd(std::span<const Vec3f> a, std::span<const Vec3f> b, std::span<const int> type) {
  if (a.empty()) return 0;
  double sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < b.size(); ++j) {
      if (type[j] != type[i]) continue;
      const double dx = a[i][0] - b[j][0], dy = a[i][1] - b[j][1], dz = a[i][2] - b[j][2];
      best = std::min(best, dx * dx + dy * dy + dz * dz);
    }
    sum += best;
  }
  return std::sqrt(sum / static_cast<double>(a.size()));
}

std::vector<RoundSample> selectTopSamples(std::vector<RoundSample> merged, const MobmHeavyAtoms& atoms, int topNum,
                                          double diversityRmsd) {
  std::stable_sort(merged.begin(), merged.end(), [](const auto& a, const auto& b) { return a.s.energy < b.s.energy; });
  std::vector<RoundSample> top;
  std::vector<std::vector<Vec3f>> keptXyz;
  for (auto& m : merged) {
    if (static_cast<int>(top.size()) >= topNum) break;
    if (diversityRmsd > 0) {
      auto xyz = atoms.world(m.s);
      bool dup = false;
      for (std::size_t k = 0; k < top.size() && !dup; ++k)
        dup = top[k].stam == m.stam && symmetricRmsd(xyz, keptXyz[k], atoms.types()) < diversityRmsd;
      if (dup) continue;
      keptXyz.push_back(std::move(xyz));
    }
    top.push_back(std::move(m));
  }
  return top;
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
