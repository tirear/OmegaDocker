// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "molecule.hpp"
#include "ring_ops.hpp"
#include "topology.hpp"
#include "torsion_tree.hpp"

namespace od {

// Ring conformations (spec 4.9).

struct RingParams {
  int macroMinSize = 9;          // ring.macroMinSize: rings of at least this many atoms are macrocycles
  int confGenNum = 300;          // ring.confGenNum: ETKDG candidates per small ring system
  double confRmsd = 15.0;        // ring.confRmsd: degrees, ring dihedral RMSD of the same conformer
  double confEnergyWindow = 5.0; // ring.confEnergyWindow: kcal/mol above the lowest candidate
  int maxConfNum = 10;           // ring.maxConfNum: conformers kept per system, the input included
  std::uint64_t seed = 0;
};

// A small ring system of a molecule (indices into Molecule::atoms of the base molecule).
struct RingSystem {
  std::vector<int> ring;   // heavy atoms of the rings (sorted)
  std::vector<int> shell;  // heavy atoms bonded to a ring atom, outside the system (no metals; sorted)
  int largestRing = 0;     // atoms of the largest smallest-cycle of its bonds
  // Conformers: positions of the ring atoms, then of the shell atoms. Conformer 0 is the input;
  // the others are in their own frame (RDKit), aligned when the ring operations are built.
  std::vector<std::vector<Vec3>> conf;
  std::vector<double> energy;  // MMFF94 (UFF when MMFF has no parameters) of conformers 1.., kcal/mol, relative
  std::string note;            // why the system has only the input conformer (empty when it has more)
  int generated = 0;           // ETKDG conformers embedded
};

struct RingConformers {
  std::vector<RingSystem> systems;
  int conformerCount() const;  // over all systems, the input conformers included
};

// Smallest cycle (atom count) through each bond once the broken bonds and metal bonds are removed:
// the count when below limit, limit for longer cycles, 0 for bonds in no cycle. Indexed like Molecule::bonds.
std::vector<int> smallestCycles(const Molecule& m, const std::vector<int>& broken, int limit);

// Small ring systems: connected sets of the bonds whose smallest cycle is shorter than
// macroMinSize (no metal bonds; broken bonds removed). No conformers yet.
std::vector<RingSystem> findRingSystems(const Molecule& m, const std::vector<int>& broken, int macroMinSize);

// Bonds to break so that no macrocycle remains (spec 4.9): repeatedly, a bond of a cycle of
// at least macroMinSize atoms (not aromatic, not a metal bond, not in a small ring) is broken,
// preferring the torsion categories single, amide, resonance, symmetric, then any, and the lowest
// bond index. StaM: only macrocycles with a bond in the flexible range (Torsion::inFlexRange).
std::vector<int> macroRingBreaks(const Molecule& m, const Topology& t, int macroMinSize);

// Ring systems of a molecule with their conformers. A system gets conformers when any torsion
// (Topology::torsions, ring torsions included) with an atom in it is in the flexible range, not all
// of its ring bonds are aromatic and each bond from a ring atom to an atom outside it leaves the
// system (not part of a larger unbroken ring); otherwise only the input. Deterministic.
RingConformers generateRingConformers(const Molecule& m, const Topology& t, const RingParams& p,
                                      std::uint64_t stream);

// Ring operations of a molecule (possibly a protonation variant: toBase maps its atoms to the
// atoms of the molecule rc was made for, -1 for added atoms; empty: the same molecule) and its
// torsion tree. Systems with one conformer, or that cannot be placed, get no operation; notes
// (optional) receives the reasons.
RingOps buildRingOps(const Molecule& m, const TorsionTree& tree, const RingConformers& rc, std::span<const int> toBase,
                     std::vector<std::string>* notes = nullptr);

}  // namespace od
