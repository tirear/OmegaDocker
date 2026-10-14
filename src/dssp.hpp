// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <vector>

#include "molecule.hpp"
#include "topology.hpp"

namespace od {

// Secondary structure by the DSSP algorithm (Kabsch & Sander 1983; spec 4.6),
// including the kappa (PPII) helix of DSSP 4. All molecules are analyzed
// together, so bridges between molecules are found. Sets Topology::dsspSS of
// each molecule with a protein backbone.
void assignDssp(const std::vector<const Molecule*>& mols, const std::vector<Topology*>& topos);

}  // namespace od
