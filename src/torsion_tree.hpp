// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <span>
#include <vector>

#include "molecule.hpp"
#include "topology.hpp"

namespace od {

// Torsion trees, anchors, loops and closure constraints (spec 4.7, 4.8).

enum class TorsionMode : std::uint8_t { Fixed, Free, CisTrans };

// Consecutive non-fixed backbone units of a StaM chain (spec 4.7); only with
// stamFixedResMovable = false.
struct Loop {
  std::vector<int> units;  // chain order (indices into Topology::units)
  int before = -1;         // fixed unit before / after the loop, -1 at a chain end
  int after = -1;
  int breakBond = -1;      // first broken backbone link between before and after, -1 if none
};

struct TreeTorsion {
  int torsion = -1;                   // index into Topology::torsions
  int fixedAtom = -1, movingAtom = -1;  // bond atoms; movingAtom is on the moving side
  int parent = -1;                    // enclosing tree torsion, -1 at the root
  double inputAngle = 0.0;            // Torsion::angle (degrees)
  std::vector<int> moving;            // coordinate indices moved (>= atom count: virtual atoms)
};

// Copy of atom copyOf rigidly attached to parent; closure pulls it onto copyOf.
struct VirtualAtom {
  int parent = -1;
  int copyOf = -1;
  Vec3 xyz{};
};

// Metal bond broken with metalBond = "break": |a - b| is held at r0.
struct DistanceRestraint {
  int a = -1, b = -1;
  double r0 = 0.0;
};

// One broken bond a-b (spec 4.8): indices into virtuals (-1 if none) of the
// copies riding with a (B', Y') and with b (A'), or, for a metal bond, into restraints.
struct ClosureLink {
  int a = -1, b = -1, y = -1;  // y: the atom copied as Y'
  int copyB = -1, copyY = -1, copyA = -1;
  int restraint = -1;
};

struct TorsionTree {
  int atomCount = 0;
  std::vector<int> anchors;                // atoms that never move by torsions
  std::vector<TreeTorsion> torsions;       // rotatable torsions, parents before children
  std::vector<VirtualAtom> virtuals;       // coordinate index atomCount + i
  std::vector<DistanceRestraint> restraints;
  std::vector<ClosureLink> links;          // order of Topology::brokenBonds
  std::vector<Loop> loops;
};

struct TreeProblem {
  enum class Kind : std::uint8_t { BothSidesFixed, PartWithoutAnchor } kind;
  int torsion = -1;  // BothSidesFixed: index into Topology::torsions
  int atom = -1;     // PartWithoutAnchor: an atom of the part
};

// Units whose backbone angles are all fixed (spec 4.7).
std::vector<bool> fixedUnits(const Topology& t, const std::vector<TorsionMode>& mode);

// Anchor atoms (spec 4.7): MobM, or StaM with stamFixedResMovable or without
// fixed units: the heavy atom (with a backbone: the CA / C4') nearest the
// heavy-atom centroid; otherwise the backbone atoms of the fixed units.
std::vector<int> anchorAtoms(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                             bool stamFixedResMovable);

std::vector<Loop> findLoops(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                            bool stamFixedResMovable);

// Backbone link bond between unit and the next unit (protein C-N, nucleic O3'-P); -1 if none.
int backboneLink(const Molecule& m, const Topology& t, int unit, int next);

// Builds the tree from the torsion modes and Topology::brokenBonds.
TorsionTree buildTorsionTree(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                             bool stamFixedResMovable, std::vector<TreeProblem>& problems);

// Coordinates (atoms, then virtual atoms) for absolute angles in degrees, one
// per TorsionTree::torsions.
std::vector<Vec3> buildCoordinates(const Molecule& m, const TorsionTree& tree, std::span<const double> anglesDeg);
// Single precision (9b, spec 8.6): all arithmetic in float.
std::vector<Vec3f> buildCoordinates(const Molecule& m, const TorsionTree& tree, std::span<const float> anglesDeg);

struct ClosureState {
  double energy = 0.0;        // G * sum of squared deviations (spec 4.5)
  double maxDeviation = 0.0;  // compared with maxBreakDist (spec 4.8)
};
ClosureState evaluateClosure(const TorsionTree& tree, std::span<const Vec3> xyz, double G);

}  // namespace od
