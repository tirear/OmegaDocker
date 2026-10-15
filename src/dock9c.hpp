// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <limits>
#include <span>
#include <unordered_map>
#include <vector>

#include "dock9b.hpp"
#include "energy.hpp"
#include "symrmsd.hpp"

namespace od {

// 9c: ranking, diversity filtering, topSampleSet and dcss filtering (spec 9).

// A docking sample with the conformation samples of its pair.
struct RoundSample {
  int stam = 0;  // StaM conformation sample id
  int mobm = 0;  // MobM conformation sample id
  DockSample s;
};

// MobM heavy atoms (probe atoms of 8.1 without H, HD, HS) and their world coordinates.
class MobmHeavyAtoms {
 public:
  explicit MobmHeavyAtoms(const EnergySystem& mobSys);
  std::size_t size() const { return index_.size(); }
  std::span<const int> types() const { return type_; }
  std::span<const int> atoms() const { return atom_; }  // atom indices in the MobM
  // World coordinates, single precision as in 9b: R(float q)(y - float center) + float t.
  std::vector<Vec3f> world(const DockSample& s) const;

 private:
  const EnergySystem& sys_;
  std::vector<int> index_;  // coordinate indices
  std::vector<int> atom_;   // atom indices
  std::vector<int> type_;   // Ad4Type
};

// Steps 1–4: merged = previous topSampleSet followed by the docking samples of
// the round. Sorted by energy (stable), greedy diversity filter (RMSD with atom i
// paired with atom i), first topNum. diversityRmsd <= 0 disables the filter.
std::vector<RoundSample> selectTopSamples(std::vector<RoundSample> merged, const MobmHeavyAtoms& atoms, int topNum,
                                          double diversityRmsd);

// 9d: in order, removes the samples closer than rmsd (symmetry-corrected) to a
// kept sample with the same StaM conformation sample.
std::vector<RoundSample> symmetryFilter(std::vector<RoundSample> top, const MobmHeavyAtoms& atoms,
                                        const SymmetryRmsd& sym, double rmsd);

// Step 5: removes the samples whose lowest docking energy of the round is
// above eTopMax, or that have none (absent from lowest). eTopMax: highest
// energy of the topSampleSet; -inf when it is empty. Returns the number removed.
int filterDcss(std::vector<int>& dcss, const std::unordered_map<int, float>& lowest, float eTopMax);

inline float topMaxEnergy(const std::vector<RoundSample>& top) {
  float e = -std::numeric_limits<float>::infinity();
  for (const auto& t : top) e = std::max(e, t.s.energy);
  return e;
}

}  // namespace od
