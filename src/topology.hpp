// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

enum class BackboneType : std::uint8_t { None, Protein, Nucleic };

// Torsion classes of spec 4.4. Only Single is rotatable by default.
enum class TorsionCategory : std::uint8_t {
  Single,
  Amide,
  Resonance,
  Double,
  Triple,     // triple bond, or an end atom collinear with all its heavy neighbors
  Aromatic,
  Ring,
  Symmetric,
  Metal,
};

enum class BackboneAngle : std::uint8_t { None, Phi, Psi, Omega, Alpha, Beta, Gamma, Delta, Epsilon, Zeta };

// Backbone repeat unit (spec 4.1).
//   protein: atoms = {N, CA, C, O, -1, -1}
//   nucleic: atoms = {P (-1 at a 5' end without phosphate), O5', C5', C4', C3', O3'}
struct BackboneUnit {
  std::array<int, 6> atoms{-1, -1, -1, -1, -1, -1};
  int residue = -1;  // residue of CA (protein) or C4' (nucleic); -1 without residue information
  int prev = -1;     // neighboring units along the chain (indices into Topology::units)
  int next = -1;
};

// A dihedral defined by four heavy atoms around the bond atoms[1]-atoms[2].
// Bonds with an end atom that has no other heavy neighbor have no such
// dihedral and are not listed.
struct Torsion {
  std::array<int, 4> atoms{};
  int bond = -1;                 // index into Molecule::bonds
  TorsionCategory category = TorsionCategory::Single;
  BackboneAngle backbone = BackboneAngle::None;
  int unit = -1;                 // backbone unit for backbone angles
  int residue = -1;              // residue of the unit (backbone) or of atoms[1]
  double angle = 0.0;            // input value in degrees, (-180, 180]
  bool inFlexRange = true;       // StaM: within the docking box extended by stamFlexMargin
  bool broken = false;           // the bond is broken and closed by constraints (spec 4.8)
  bool rotatable = false;        // default of spec 4.2 / 4.4
};

struct Topology {
  BackboneType backboneType = BackboneType::None;
  std::vector<BackboneUnit> units;
  std::vector<Torsion> torsions;
  // Crosslinks to break automatically (spec 4.7): bond indices of inter-residue
  // bonds, other than backbone links and metal bonds, that close a cycle of
  // the bond graph without metal bonds.
  std::vector<int> crosslinkBreaks;
  // Metal bonds to break with metalBond = "break" (spec 4.7): bond indices of
  // metal bonds that close a cycle once crosslinkBreaks are broken.
  std::vector<int> metalBreaks;
  // Bonds broken so that no macrocycle remains (spec 4.9, macroRingBreaks); set by the caller.
  std::vector<int> macroBreaks;
  // Bonds treated as broken when classifying ring bonds (spec 4.5, 4.8).
  std::vector<int> brokenBonds;
  // DSSP code per residue (spec 4.6; set by assignDssp): H G I P E B T S, '.' for
  // none, ' ' for residues outside the protein backbone. Empty if not assigned.
  std::string dsspSS;
};

struct Box {
  Vec3 center{};
  Vec3 size{};
};

struct TopologyParams {
  // StaM only: torsions are rotatable by default only in residues (or, without
  // residue information, at atoms) within the box extended by stamFlexMargin.
  std::optional<Box> box;
  double stamFlexMargin = 4.0;
  // StaM (box.mode = fixAuto, spec 3.3): no torsion is rotatable by default.
  bool stamAllFixed = false;
  // StaM (spec 4.7): false keeps the backbone of fixed residues in place; with
  // false, torsions with backbone atoms on both sides are not rotatable by default.
  bool stamFixedResMovable = false;
  // Broken bonds (indices into Molecule::bonds); crosslinkBreaks when not given
  // (metal bonds are kept by default, spec 4.7).
  std::optional<std::vector<int>> brokenBonds;
};

Topology analyzeTopology(const Molecule& m, const TopologyParams& p);

// Indices of omega torsions whose input angle differs from both 0 and 180
// degrees by more than toleranceDeg (spec 4.3).
std::vector<int> nonPlanarOmegas(const Topology& t, double toleranceDeg);

// Bonds that lie in a ring once the broken bonds (indices into Molecule::bonds)
// are removed; broken bonds are not in a ring. Indexed like Molecule::bonds.
std::vector<bool> ringBondMask(const Molecule& m, const std::vector<int>& broken = {});

double dihedralDeg(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

const char* backboneTypeName(BackboneType t);       // "protein" | "nucleic" | "none"
const char* torsionCategoryName(TorsionCategory c);  // "single", "amide", ...
const char* backboneAngleName(BackboneAngle a);      // "phi", ..., "" for None

}  // namespace od
