// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <span>
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// Structure output of 9d (spec 10.3).

// One molecule with output coordinates (one per atom of mol).
struct StructPart {
  const Molecule* mol = nullptr;
  std::vector<Vec3f> xyz;
};

struct OutFormats {
  bool pdb = true, cif = false, sdf = false, mol2 = false;
};

// Writes base + ".pdb" / ".cif" / ".sdf" / ".mol2" for the enabled formats.
// PDB is replaced by mmCIF when the structure exceeds the PDB limits (only one
// mmCIF is written when both apply). Returns the file names written (without
// directory) and sets pdbFallback when PDB was replaced.
std::vector<std::string> writeStructure(const std::string& base, std::span<const StructPart> parts,
                                        const OutFormats& f, const std::string& title, bool& pdbFallback);

// True when the parts fit in the PDB format (atom serials, residue numbers,
// one-character chain IDs, residue names of at most 3 characters).
bool fitsPdb(std::span<const StructPart> parts);

// Writers (exposed for tests).
std::string pdbText(std::span<const StructPart> parts, const std::string& title);
std::string cifText(std::span<const StructPart> parts, const std::string& title);
std::string sdfText(std::span<const StructPart> parts, const std::string& title);
std::string mol2Text(std::span<const StructPart> parts, const std::string& title);
// SYBYL atom type of atom a (mol2).
std::string sybylType(const Molecule& m, int a);

}  // namespace od
