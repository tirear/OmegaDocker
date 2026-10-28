// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>
#include <span>
#include <utility>
#include <vector>

#include "energy.hpp"
#include "rng.hpp"
#include "topology.hpp"

namespace od {

// csgmBuild = grow / growClose (spec 7.8, 4.8): the tree is built from the
// root one segment at a time; after each segment the heavy-atom pairs completed
// by it are checked for clashes, and the distances to the closure targets of
// broken backbone links for reach. A failing segment is resampled; when its
// retries run out, its parent segment is resampled (backtracking).

struct GrowParams {
  int maxRetry = 16;           // growMaxRetry (csgm1; csgm2 builds without resampling)
  double clashDistance = 1.8;
  double reachProtein = 3.8;   // reachPerRes (Å per residue)
  double reachNucleic = 7.0;
  int maxBacktrack = 4;        // growMaxBacktrack: parent resamplings per segment
  bool reach = true;           // false with loop closure (spec 7.8)
  // growClose: per concatenated torsion, a torsion of a closure path. A segment
  // with such a torsion is applied as given: no checks, no resampling; the pairs
  // completed by it are left to the loop closure.
  std::vector<bool> fixed;
};

struct GrowResult {
  bool ok = false;
  int retries = 0;      // segment resamplings
  int backtracks = 0;   // parent segment resamplings
  int clashFails = 0;   // failed checks by cause
  int reachFails = 0;
};

using SegmentSampler = std::function<void(Rng& rng, int segment, std::span<double> anglesDeg)>;

class GrowBuilder {
 public:
  // parts as given to sys; topologies: Topology the trees were built from.
  // segments: indices into the concatenated angles, covering every torsion once.
  // Throws od::Error when the segments cannot be ordered from the root.
  GrowBuilder(const EnergySystem& sys, const std::vector<EnergyPart>& parts,
              const std::vector<const Topology*>& topologies, const std::vector<std::vector<int>>& segments,
              const GrowParams& p);

  // anglesDeg holds a complete sample; failing segments are resampled with
  // resample (at most maxRetry times each, then up to maxBacktrack resamplings
  // of the parent segment, one level up; never when resample is empty, as for
  // csgm2). ok = false: the sample is discarded.
  // coords (optional): coordinates of the built sample, as EnergySystem::evaluate.
  GrowResult build(Rng& rng, std::span<double> anglesDeg, const SegmentSampler& resample,
                   std::vector<Vec3>* coords = nullptr) const;

  int stepCount() const { return static_cast<int>(steps_.size()); }
  int fixedStepCount() const;
  std::size_t reachCheckCount() const;

 private:
  struct Tor {
    int parent = -1;  // global tree torsion index
    Vec3 o{}, k{};    // axis point and unit vector in input coordinates
    int fixedC = -1, movingC = -1;  // coordinate indices of the axis atoms
    double inputAngle = 0;
    std::vector<int> own;  // coordinate indices moved by this torsion but by none of its children
  };
  struct Reach {
    int current = -1, target = -1;  // coordinate indices
    double bound2 = 0;              // squared largest distance
  };
  struct Step {
    int segment = -1;
    bool fixed = false;
    int parentStep = -1;        // step of the parent torsion outside the segment
    std::vector<int> torsions;  // increasing (parents first)
    std::vector<std::pair<int, int>> pairs;
    std::vector<Reach> reach;
  };
  const EnergySystem* sys_ = nullptr;
  std::vector<Vec3> input_;
  std::vector<Tor> tor_;
  std::vector<Step> steps_;
  double clash2_ = 0;
  int maxRetry_ = 0;
  int maxBacktrack_ = 0;
};

}  // namespace od
