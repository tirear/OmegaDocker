// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// SDF / MOL (V2000, V3000): one RawStructure per record. Hydrogens are kept as read.
std::vector<RawStructure> readSdf(const std::string& path);

// mol2: one RawStructure per @<TRIPOS>MOLECULE record; atom names, substructures
// (residues) and partial charges are taken from the file.
std::vector<RawStructure> readMol2(const std::string& path);

}  // namespace od
