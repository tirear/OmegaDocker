// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Rule-based protonation states of small molecules, nucleic acids and non-standard residues
// (spec 3.4): functional groups are found with SMARTS and given literature pKa values of
// representative compounds; amines get beta-substituent corrections, sites that are charged
// earlier shift the pKa of later sites with the same charge sign.

#pragma once
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
class RWMol;
}

namespace od::pkarules {

struct Site {
  std::string name;  // rule name
  int atom = -1;     // heavy-atom index
  double pka = 0;    // pKa used for the decision (after corrections)
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

}  // namespace od::pkarules
