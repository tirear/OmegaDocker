// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <string>
#include <vector>

#include "ccd.hpp"
#include "molecule.hpp"
#include "protonate.hpp"

namespace od {

// Reads the first model of an mmCIF file. Bonds come from CCD definitions,
// polymer links between consecutive residues, and _struct_conn (spec 3.1).
// connectivityOnly: no formal charges and bond orders (e.g. a reference structure without hydrogens).
// prot: hydrogens of the file are skipped and replaced (spec 3.4); planOut / planIn as protonate.
RawStructure readMmcif(const std::string& path, const CcdLibrary& ccd, bool connectivityOnly = false,
                       const ProtonationOptions* prot = nullptr, ProtonationReport* report = nullptr);
// Every model of an mmCIF file, each read as readMmcif (score and localMin poses, spec 6.1); with
// prot, the protonation of the first model is applied to the others.
std::vector<RawStructure> readMmcifModels(const std::string& path, const CcdLibrary& ccd,
                                          const ProtonationOptions* prot = nullptr, ProtonationReport* report = nullptr);

}  // namespace od
