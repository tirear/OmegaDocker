// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#include "ad4.hpp"
#include "molecule.hpp"
#include "torsion_tree.hpp"

namespace od {

// 9a energy (spec 7.9): intramolecular AD4 terms, clash penalty and closure
// constraints, as a function of the tree torsions of one or more molecules
// (a StaM group or one MobM). CPU reference implementation.

struct EnergyParams {
  Ad4Coeffs coeffs;
  double smooth = 0.5;            // Å, vdW / H-bond smoothing (AutoDock-GPU default)
  double elecMinDistance = 0.01;  // Å, AutoDock-GPU default
  double clashDistance = 1.8;
  double clashBase = 10000.0;
  double clashK = 1000.0;
  double closureG = 200.0;
};

struct EnergyTerms {
  double vdw = 0, hbond = 0, elec = 0, desolv = 0;
  double clash = 0;
  double closure = 0;
  int clashPairs = 0;
  double closureMaxDeviation = 0;  // Å, compared with maxBreakDist (spec 4.8)
  double external = 0;             // ExternalTerm energy; not part of the sums below
  double intra() const { return vdw + hbond + elec + desolv; }
  // Compared with E_max and used for ranking; closure excluded (spec 4.5).
  double filterEnergy() const { return intra() + clash; }
  double total() const { return filterEnergy() + closure; }
};

struct EnergyPart {
  const Molecule* mol = nullptr;
  const TorsionTree* tree = nullptr;
};

// Extra energy term on the coordinates (e.g. the 9b grid energy): adds its
// dE/dx to grad (nullptr when no gradient is requested) and returns its energy.
template <class T>
using ExternalTermT = std::function<double(const std::vector<Vec3T<T>>& coords, std::vector<Vec3>* grad)>;
using ExternalTerm = ExternalTermT<double>;
using ExternalTermF = ExternalTermT<float>;

class EnergySystem {
 public:
  // Throws od::Error when an atom has no AD4 parameters, or when two heavy
  // atoms whose distance cannot change clash in the input structure.
  EnergySystem(std::vector<EnergyPart> parts, const EnergyParams& params);

  int torsionCount() const { return torsionCount_; }
  std::vector<double> inputAngles() const;  // degrees
  std::size_t pairCount() const { return pairs_.size(); }
  const Ad4Typing& typing(int part) const { return typing_[part]; }
  // Index of atom 0 of a part in the coordinate array; virtual atoms follow the atoms.
  int coordOffset(int part) const { return offset_[part]; }

  // anglesDeg: absolute angles of the tree torsions of all parts, in part
  // order. gradient (optional): dE/dangle of total() in kcal/mol/rad.
  // external (optional): its gradient enters the torsion gradient, its energy goes to EnergyTerms::external.
  EnergyTerms evaluate(std::span<const double> anglesDeg, std::vector<double>* gradient = nullptr,
                       std::vector<Vec3>* coords = nullptr, const ExternalTerm* external = nullptr) const;
  // Single precision (9b, spec 8.6): coordinates and per-pair terms in float;
  // energy sums and gradient accumulation in double.
  EnergyTerms evaluate(std::span<const float> anglesDeg, std::vector<double>* gradient, std::vector<Vec3f>* coords,
                       const ExternalTermF* external) const;

  struct Pair {
    int i, j;            // coordinate indices
    std::uint16_t type;  // index into the type-pair table
    bool heavy;          // both heavy atoms: clash penalty applies
  };
  const std::vector<Pair>& pairs() const { return pairs_; }

 private:
  template <class T>
  EnergyTerms evaluateT(std::span<const T> anglesDeg, std::vector<double>* gradient, std::vector<Vec3T<T>>* coords,
                        const ExternalTermT<T>* external) const;

  std::vector<EnergyPart> parts_;
  EnergyParams params_;
  std::vector<Ad4Typing> typing_;
  std::vector<int> offset_;
  std::vector<int> torsionOffset_;
  int torsionCount_ = 0;
  int coordCount_ = 0;
  std::vector<double> charge_;  // per coordinate index
  std::vector<double> sv_;      // per coordinate index: solpar + qasp*|q|
  std::vector<double> vol_;
  std::vector<Ad4PairParams> pairParams_;
  std::vector<Pair> pairs_;
};

}  // namespace od
