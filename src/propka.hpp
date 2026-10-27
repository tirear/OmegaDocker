// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// C++ port of the pKa prediction of PROPKA 3 (LGPL-2.1, https://github.com/jensengroup/propka,
// commit 58ceb7b), version "VersionA" with the parameters of its propka.cfg (spec 3.4).

#pragma once
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od::propka {

// One heavy atom of the input (mmCIF atom_site order). Hydrogens are not passed.
struct InputAtom {
  std::string name;
  std::string resName;  // without padding
  std::string chain;
  int resNum = 0;
  std::string iCode;    // "" or one character
  std::string element;  // e.g. "C", "Zn"
  bool hetatm = false;  // HETATM record
  Vec3 xyz{};
};

struct GroupResult {
  int atom = -1;            // index into the input atoms
  std::string type;         // PROPKA group type: COO, HIS, CYS, TYR, LYS, ARG, N+, ...
  std::string residueType;  // ASP, GLU, C-, N+, HIS, ... (model pKa key)
  std::string label;        // e.g. "ASP  12 A"
  double pka = 0;
  double modelPka = 0;
  double buried = 0;        // fraction
  bool penalised = false;   // "NB: Discarded due to coupling"
};

// Titratable groups of the structure in PROPKA group order (one conformation).
std::vector<GroupResult> calculate(const std::vector<InputAtom>& atoms);

}  // namespace od::propka
