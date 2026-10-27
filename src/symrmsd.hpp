// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstddef>
#include <span>
#include <vector>

#include "molecule.hpp"
#include "orient.hpp"

namespace od {

// RMSD without superposition between two coordinate sets of the same atoms.

// Atom i paired with atom i.
double fixedRmsd(std::span<const Vec3f> a, std::span<const Vec3f> b);

// Automorphisms of a labelled graph (adjacency lists, one label per vertex);
// stops after limit + 1 automorphisms. perm[i]: image of vertex i.
std::vector<std::vector<int>> graphAutomorphisms(const std::vector<std::vector<int>>& adj, std::span<const int> label,
                                                 std::size_t limit);

// One isomorphism from graph A onto graph B (labels must match); empty when there is none.
std::vector<int> graphIsomorphism(const std::vector<std::vector<int>>& adjA, std::span<const int> labelA,
                                  const std::vector<std::vector<int>>& adjB, std::span<const int> labelB);

// Reference structure (spec 10.3): for each of atoms (heavy atoms of mol), the index of the
// corresponding ref atom, by the heavy-atom bond graph labelled with the element (metal
// coordination bonds excluded); the same atom order when it fits. Throws od::Error when none.
std::vector<int> referenceAtomMap(const Molecule& mol, std::span<const int> atoms, const Molecule& ref);

// Smallest sum of squared distances over one-to-one pairings of a and b within
// each group (Hungarian method); groups: index lists covering all atoms.
double hungarianSquaredSum(std::span<const Vec3f> a, std::span<const Vec3f> b,
                           const std::vector<std::vector<int>>& groups);

// Symmetry-corrected RMSD (spec 10): the smallest RMSD over the automorphisms
// of the heavy-atom graph (labels: AD4 type) when there are at most maxPerm of
// them (C1); otherwise the best one-to-one pairing within each AD4 type (C2).
class SymmetryRmsd {
 public:
  // atoms: indices into mol of the compared atoms; type: label of each.
  SymmetryRmsd(const Molecule& mol, std::span<const int> atoms, std::span<const int> type, std::size_t maxPerm);
  bool usesPermutations() const { return !perms_.empty(); }
  std::size_t permutationCount() const { return perms_.size(); }  // 0 with C2
  double rmsd(std::span<const Vec3f> a, std::span<const Vec3f> b) const;

 private:
  std::vector<std::vector<int>> perms_;
  std::vector<std::vector<int>> groups_;  // per AD4 type, for C2
};

}  // namespace od
