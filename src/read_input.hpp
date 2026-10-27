// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <string>
#include <vector>

#include "ccd.hpp"
#include "molecule.hpp"

namespace od {

// Reads a structure file by extension: .cif/.mmcif, .sdf/.sd/.mol, .mol2.
std::vector<RawStructure> readStructureFile(const std::string& path, const CcdLibrary& ccd);

// Reads a file, splits it into molecules, assigns ids and perceives chemistry.
//   id = "<sta|mob>:<file name>:<chains joined by '+'>" with residue information
//        (":<resName><seq><iCode>" of the first residue is appended when not unique),
//   id = "<sta|mob>:<file name>:<k>" without residue or chain information
//        (k = 1-based molecule index in the file).
std::vector<Molecule> loadMolecules(const std::string& path, Role role, const CcdLibrary& ccd);

// score / localMin (spec 6.1): every record (SDF / mol2 record, mmCIF model) of a MobM file is a
// pose of the same molecule. mol: the first record, as loadMolecules; xyz: the coordinates of
// every record (the first included). Records must have the same atoms (element, name) and bonds.
struct PoseRecords {
  Molecule mol;
  std::vector<std::vector<Vec3>> xyz;
};
PoseRecords loadPoseRecords(const std::string& path, Role role, const CcdLibrary& ccd);

// Reference structure (--refFiles, spec 10.3): one molecule, its first record only; atoms,
// coordinates and connectivity (hydrogens may be missing; no chemistry perception).
Molecule loadReference(const std::string& path, const CcdLibrary& ccd);

}  // namespace od
