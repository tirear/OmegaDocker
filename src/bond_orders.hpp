// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// For structures whose bonds come from CCD definitions (mmCIF input).
//
// Non-aromatic bonds keep the CCD order; formal charges follow from the valence
// of each atom with the hydrogens actually present in the file.
// Bonds flagged aromatic in the CCD form conjugated systems whose single/double
// pattern is re-derived for the hydrogens present (e.g. HIS tautomers); a
// pyrrole-type N is made N+ only when no neutral pattern exists.
//
// `aromaticCandidate[i]` marks bonds flagged aromatic in the CCD.
// Problems are appended to `errors` as human-readable messages.
void assignChargesAndKekule(RawStructure& s, const std::vector<bool>& aromaticCandidate,
                            const std::vector<int>& ccdCharge,
                            std::vector<std::string>& errors);

// Describes an atom for messages, e.g. "A:12:HIS:NE2" or "#15".
std::string describeAtom(const RawStructure& s, int atom);

}  // namespace od
