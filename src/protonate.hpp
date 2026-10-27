// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "ccd.hpp"
#include "molecule.hpp"

namespace od {

// [[protonation]] entry of the molecule configuration file (spec 3.4, 5).
struct ProtonationOverride {
  std::string file;            // input file as given on the command line
  std::string site;            // "<chain>:<seq><iCode>" (residue) or an atom reference
  std::string state;           // residue state (ASH, HIE, ...); empty for an atom entry
  std::optional<int> charge;   // atom entry: formal charge
  std::string origin;          // "<config file>:<line>" for messages
};

struct ProtonationOptions {
  double pH = 7.4;
  std::vector<ProtonationOverride> overrides;  // entries of this file only
};

struct ProtonationSite {
  std::string atom;    // describeAtom
  std::string method;  // "PROPKA", "Dimorphite-DL <substructure>", "override", "metal"
  double pka = std::numeric_limits<double>::quiet_NaN();
  std::string state;   // residue state (ASP, ASH, HIE, ...) or protonated / deprotonated
};

struct ProtonationReport {
  int removedH = 0;
  int addedH = 0;
  std::vector<std::string> addedOxt;    // residues
  std::vector<std::string> incomplete;  // residues with missing heavy atoms
  std::vector<ProtonationSite> sites;
  int rotors = 0, flipGroups = 0;
  std::vector<std::string> flipped;     // flipped Asn / Gln / His / carboxyl groups
  int components = 0, greedyComponents = 0;
};

// The hydrogens and charges chosen for one record, for the other records of the same molecule
// (score / localMin poses): the same atoms are created, only rotors and flips are optimized.
struct ProtonationPlan {
  std::vector<int> formalCharge;        // heavy atoms (after OXT completion)
  std::vector<BondOrder> bondOrder;     // heavy-atom bonds
  std::vector<std::pair<int, std::string>> hydrogens;  // (parent heavy atom, name)
  std::vector<std::string> addedOxt;
};

// Replaces the hydrogens of `s` by those of the dominant protonation state at opt.pH
// (spec 3.4). mmCIF structures (ccd != nullptr) carry per-atom CCD charges and per-bond CCD
// aromatic flags, which are kept in step with the atoms and bonds; their formal charges and
// Kekulé bonds are then assigned by assignChargesAndKekule (errors appended to `errors`).
// Other structures get formal charges and bond orders directly.
void protonate(RawStructure& s, const CcdLibrary* ccd, const ProtonationOptions& opt, ProtonationReport& rep,
               std::vector<int>* ccdCharge, std::vector<bool>* aromaticCandidate, std::vector<std::string>& errors,
               ProtonationPlan* planOut = nullptr, const ProtonationPlan* planIn = nullptr);

// Residue states accepted by [[protonation]] entries, per residue name.
const std::vector<std::string>& residueStates(const std::string& resName);

}  // namespace od
