// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "energy.hpp"
#include "grow.hpp"
#include "rng.hpp"

namespace od {

// Loop closure of the broken bonds after csgm and before 9a minimization
// (spec 7.11): RCD (random coordinate descent, Chys & Chacón 2013) and KIC
// (kinematic closure, Coutsias et al. 2004; Mandell et al. 2009).

enum class ClosureMethod : std::uint8_t { None, Rcd, Kic };

struct ClosureParams {
  ClosureMethod method = ClosureMethod::Rcd;
  double rcdTolerance = 0.25;  // Å, RMSD of the closure atom pairs
  int rcdMaxIter = 1000;
  int rcdMaxRestart = 10;
  int kicMaxTrial = 100;
  double clashDistance = 1.8;
};

class LoopCloser {
 public:
  // freeMask: per concatenated torsion, a freely rotatable (not cisTrans) torsion.
  // segments: csgm1 segments (indices into the concatenated angles) used to resample.
  LoopCloser(const EnergySystem& sys, const std::vector<EnergyPart>& parts, const std::vector<bool>& freeMask,
             const std::vector<std::vector<int>>& segments, const ClosureParams& p);

  // Links with at least one torsion to adjust.
  bool active() const;
  // One line per link for the log.
  std::vector<std::string> describe() const;

  // Closes every link in turn (spec 7.11); false: a link could not be closed.
  // resample: csgm1 segment sampler (empty: uniform angles).
  bool close(Rng& rng, std::span<double> anglesDeg, const SegmentSampler& resample) const;

 private:
  struct Pair {
    int moving = -1, target = -1;  // part-local coordinate indices
  };
  struct Link {
    int part = 0;
    bool restraint = false;
    std::vector<Pair> pairs;        // copy links
    int ra = -1, rb = -1;           // restraint atoms
    double r0 = 0;
    std::vector<int> ring;          // copy links: path atoms a ... b (part-local)
    std::vector<int> bondTorsion;   // per ring bond (ring[j], ring[j+1]): global torsion index, -1 if not adjustable
    std::vector<int> dof;           // adjustable global torsion indices
    std::vector<int> pivots;        // ring indices usable as KIC pivots
    bool kic = false;
    bool partial = false;           // best effort: metal distance, or fewer than 6 adjustable torsions
    std::vector<int> resampleSegments;  // csgm1 segments made only of adjustable torsions
    std::vector<int> resampleCovered;   // their torsions, sorted
    // Tracked coordinates for RCD: points (part-local indices); per dof, the points it moves.
    std::vector<int> points;
    std::vector<std::vector<int>> moves;
    std::vector<int> axisFixed, axisMoving;  // per dof, index into points
    std::vector<int> pairPoint;              // per pair: moving, target as indices into points
    std::vector<std::pair<int, int>> clashPairs;  // system coordinate indices
  };
  double deviation(const Link& l, const std::vector<Vec3>& x) const;  // part-local coordinates
  bool clashFree(const Link& l, std::span<const double> anglesDeg) const;
  void resampleLink(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample) const;
  bool rcd(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample) const;
  bool kic(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample) const;
  std::vector<Vec3> partCoords(int part, std::span<const double> anglesDeg) const;

  const EnergySystem& sys_;
  std::vector<EnergyPart> parts_;
  ClosureParams p_;
  std::vector<int> partOfTorsion_, localTorsion_;
  std::vector<Link> links_;
};

}  // namespace od
