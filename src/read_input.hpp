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

}  // namespace od
