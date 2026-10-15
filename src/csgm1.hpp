// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <memory>
#include <span>
#include <vector>

#include "molconfig.hpp"
#include "molecule.hpp"
#include "rama.hpp"
#include "rng.hpp"

namespace od {

// csgm1 (spec 7.5): conformations generated from nothing, csgm1a (uniform)
// or csgm1b (secondary-structure types of protein residues).

struct Csgm1Params {
  double angleStep = 1.0;  // csgm1AngleStep (degrees)
  double ramaFloor = 1e-4;
};

// Cumulative weights of a rama8000 grid, value + ramaFloor per cell.
class RamaSampler {
 public:
  RamaSampler(const RamaGrid& g, double floor);
  void sample(Rng& rng, double& phi, double& psi) const;
  double probability(double phi, double psi) const;  // of the cell containing (phi, psi)

 private:
  std::vector<double> cdf_;
};

class Csgm1Sampler {
 public:
  // One molecule; angles follow the order of settings.tree.torsions.
  Csgm1Sampler(const Molecule& m, const MolSettings& settings, const MolConfig& config, const Csgm1Params& params);

  int size() const { return static_cast<int>(input_.size()); }
  void sample(Rng& rng, std::span<double> anglesDeg) const;

  // Torsions sampled together (csgmBuild = grow, spec 7.8): a phi/psi pair of
  // csgm1b (with the cisTrans omega before it for RAMA_AUTO), otherwise one
  // torsion. Sorted by first index; members in increasing order.
  const std::vector<std::vector<int>>& segments() const { return segments_; }
  // Resamples the angles of one segment; the other angles are kept.
  void sampleSegment(Rng& rng, int segment, std::span<double> anglesDeg) const;

 private:
  enum class Kind : std::uint8_t { Uniform, CisTrans, PairPhi, PairPsi };
  enum class Rule : std::uint8_t { R, O, Custom, Rama, RamaAuto };
  struct Pair {
    int phi = -1, psi = -1;  // tree torsion indices
    Rule rule = Rule::R;
    int custom = -1;         // index into MolConfig::ssTypes
    RamaClass rama = RamaClass::General;
    // RamaAuto
    std::string resName, nextResName;
    int omegaBefore = -1;    // tree index of a cisTrans omega before the residue, -1 if none
    bool cisBefore = false;  // fixed omega before the residue is cis
  };
  void initPairs(const Molecule& m, const MolSettings& s, const std::vector<int>& treeOf);
  double uniformAngle(Rng& rng, int i) const;
  void samplePair(Rng& rng, const Pair& p, std::span<double> a) const;
  const RamaSampler& rama(RamaClass c) const;

  std::vector<double> input_;
  std::vector<Kind> kind_;
  std::vector<Pair> pairs_;
  std::vector<std::vector<int>> segments_;
  std::vector<int> segmentPair_;  // per segment: index into pairs_, -1 if none
  std::vector<SsType> ssTypes_;
  double halfWidth_ = 15.0, originalProb_ = 0.9;
  Csgm1Params params_;
  int stepMin_ = 0, stepCount_ = 0;  // Δ = (stepMin_ + k) * angleStep, k in [0, stepCount_)
  std::shared_ptr<std::vector<std::unique_ptr<RamaSampler>>> rama_;
};

// Several molecules (a StaM group): angles concatenated in part order.
class Csgm1GroupSampler {
 public:
  void add(Csgm1Sampler s) { parts_.push_back(std::move(s)); }
  int size() const;
  void sample(Rng& rng, std::span<double> anglesDeg) const;
  // Segments of all parts, with indices into the concatenated angles.
  std::vector<std::vector<int>> segments() const;
  void sampleSegment(Rng& rng, int segment, std::span<double> anglesDeg) const;

 private:
  std::vector<Csgm1Sampler> parts_;
};

}  // namespace od
