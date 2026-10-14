// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <memory>

#include "molecule.hpp"

namespace RDKit {
class RWMol;
}

namespace od {

// Builds an RDKit molecule (explicit hydrogens only, Kekulé bonds, coordinates).
std::unique_ptr<RDKit::RWMol> toRdkit(const Molecule& m);

// Sanitizes the molecule, sets aromatic flags of atoms and bonds, and computes
// Gasteiger charges unless the charges were read from the file.
void perceiveChemistry(Molecule& m);

}  // namespace od
