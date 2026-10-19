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
// Clash term in the RCD objective (spec 7.11): none; repair: from the first converged
// iterate that fails the clash check; full: from the first update.
enum class RcdClash : std::uint8_t { None, Repair, Full };

struct ClosureParams {
  ClosureMethod method = ClosureMethod::Rcd;
  double rcdTolerance = 0.25;  // Å, RMSD of the closure atom pairs
  int rcdMaxIter = 20000;
  // rad; stop an attempt once every torsion of the link was updated by at most this since the
  // last larger update (a stall); 0: run up to rcdMaxIter.
  double rcdStallTolerance = 1e-4;
  int rcdMaxRestart = 10;
  int kicMaxTrial = 100;
  double clashDistance = 1.8;
  RcdClash rcdClash = RcdClash::None;
  double rcdClashWeight = 1.0;  // weight of sum (clashDistance - d)^2 against the closure sum of squares (Å^2)
  // csgmBuild = growClose (spec 7.8): the clash check after closing covers every
  // atom moved by the closure except those built after it (side chains).
  bool growClose = false;
};

// Flat form of the RCD links for a GPU (spec 15.1.18). Indices are per link: points of the link,
// its pairs; per adjustable torsion (dof, flat over the links): the points it moves, its axis and
// the pairs with exactly one point moved (active: pair index × 2 + 1 when the second point moves).
struct RcdFlat {
  std::vector<char> eligible;  // per link: closed by RCD without the clash term or metal distances
  std::vector<int> pointNum, dofOff, dofNum, pairOff, pairNum, copyOff, copyNum, members;  // per link
  std::vector<int> movesOff, moves, axisFixed, axisMoving, activeOff, active;               // per dof
  std::vector<int> pairPoint;  // per pair: moving, target (points)
  std::vector<int> pairEnd;    // per copy link: end of its pairs within the link
};

// One RCD attempt of an eligible link run elsewhere (spec 15.1.18): the points and the adjustable
// angles (degrees) at the start; afterwards the angles, the random state and whether it converged.
struct RcdJob {
  int link = 0;
  std::vector<double> points;  // pointNum × 3
  std::vector<double> anglesDeg;  // dofNum, in the order of the link's dof
  std::uint64_t rng = 0;
  int maxIter = 0;
  bool converged = false;
  int iterations = 0;
};

// Runs a job of the coroutine `slot`: stores it and waits (batchWait) until it is done.
class RcdRunner {
 public:
  virtual ~RcdRunner() = default;
  virtual void run(int slot, RcdJob& job) = 0;
};
struct RcdSlot {
  RcdRunner* runner = nullptr;
  int slot = -1;
};

class LoopCloser {
 public:
  // freeMask: per concatenated torsion, a freely rotatable (not cisTrans) torsion.
  // segments: csgm1 segments (indices into the concatenated angles) used to resample.
  LoopCloser(const EnergySystem& sys, const std::vector<EnergyPart>& parts, const std::vector<bool>& freeMask,
             const std::vector<std::vector<int>>& segments, const ClosureParams& p);

  // Links with at least one torsion to adjust.
  bool active() const;
  // One line per group of links for the log.
  std::vector<std::string> describe() const;
  // Per concatenated torsion: a torsion of a closure path (ring path bonds; for a
  // metal distance, torsions moving one atom but not the other). Empty for none.
  const std::vector<bool>& pathTorsions() const { return path_; }

  // Closes every group of links in turn (spec 7.11); false: a group could not be closed.
  // resample: csgm1 segment sampler (empty: uniform angles).
  // rs: RCD attempts of eligible links run by rs.runner (spec 15.1.18).
  bool close(Rng& rng, std::span<double> anglesDeg, const SegmentSampler& resample, RcdSlot rs = {}) const;
  // The links in flat form for a GPU.
  RcdFlat flat() const;
  const ClosureParams& params() const { return p_; }

 private:
  struct Pair {
    int moving = -1, target = -1;  // part-local coordinate indices
  };
  struct Restraint {
    int a = -1, b = -1;  // part-local atom indices
    double r0 = 0;
  };
  // A group of links sharing adjustable torsions, closed together (spec 7.11).
  struct Link {
    int part = 0;
    std::vector<int> members;       // indices into TorsionTree::links
    std::vector<Pair> pairs;        // copy links: B', Y', A' pairs of each, in member order
    std::vector<int> pairEnd;       // per copy link: end of its pairs
    std::vector<Restraint> restraints;
    std::vector<int> pathAtoms;     // ring path atoms of the copy links, the atoms of the metal distances
    std::vector<int> ring;          // single copy link: path atoms a ... b (part-local)
    std::vector<int> bondTorsion;   // per ring bond (ring[j], ring[j+1]): global torsion index, -1 if not adjustable
    std::vector<int> dof;           // adjustable global torsion indices
    std::vector<int> pivots;        // ring indices usable as KIC pivots
    bool kic = false;
    bool partial = false;           // best effort: metal distances only, or fewer than 6 adjustable torsions per copy link
    std::vector<int> resampleSegments;  // csgm1 segments made only of adjustable torsions
    std::vector<int> resampleCovered;   // their torsions, sorted
    // Tracked coordinates for RCD: points (part-local indices); per dof, the points it moves.
    std::vector<int> points;
    std::vector<std::vector<int>> moves;
    std::vector<int> axisFixed, axisMoving;  // per dof, index into points
    std::vector<int> pairPoint;              // per pair: moving, target as indices into points
    std::vector<int> restraintPoint;         // per restraint: a, b as indices into points
    std::vector<std::pair<int, int>> clashPairs;  // system coordinate indices
    std::vector<std::vector<int>> movesSys;       // per dof, system coordinates it moves (sorted)
  };
  double deviation(const Link& l, const std::vector<Vec3>& x) const;  // part-local coordinates
  bool clashFree(const Link& l, std::span<const double> anglesDeg) const;
  void resampleLink(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample) const;
  bool rcd(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample, RcdSlot rs) const;
  bool rcdEligible(const Link& l) const;
  bool kic(Rng& rng, const Link& l, std::span<double> anglesDeg, const SegmentSampler& resample) const;
  std::vector<Vec3> partCoords(int part, std::span<const double> anglesDeg) const;

  const EnergySystem& sys_;
  std::vector<EnergyPart> parts_;
  ClosureParams p_;
  std::vector<int> partOfTorsion_, localTorsion_;
  std::vector<Link> links_;
  std::vector<bool> path_;
};

}  // namespace od
