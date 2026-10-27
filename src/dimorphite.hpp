// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// C++ port of the protonation-state rules of Dimorphite-DL 2.1.0 (Apache-2.0,
// https://github.com/durrantlab/dimorphite_dl, commit 69ec38e) for one pH (spec 3.4):
// neutralization reactions, site detection with atom protection, and the charge / hydrogen
// assignment of each site. The substructure file is embedded unchanged.

#pragma once
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
class RWMol;
}

namespace od::dimorphite {

struct Site {
  std::string name;  // substructure name (without the leading '*')
  int atom = -1;     // heavy-atom index
  double pka = 0;
  bool protonated = false;
  bool forced = false;  // metal-bonded site atom: deprotonated regardless of the pKa
};

struct Result {
  std::unique_ptr<RDKit::RWMol> mol;  // heavy atoms in the input order, sanitized, Kekulé bonds
  std::vector<Site> sites;
};

// heavy: heavy atoms only (implicit hydrogens), sanitized. Atoms flagged in `fixed` are left
// unchanged (no neutralization, no site); site atoms flagged in `metalBonded` are deprotonated.
// The dominant state at pH is used: protonated when pKa >= pH.
Result protonate(const RDKit::RWMol& heavy, double pH, const std::vector<bool>& fixed,
                 const std::vector<bool>& metalBonded);

}  // namespace od::dimorphite
