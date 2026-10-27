// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// AD4 atom types and parameters (spec 7.9), values of AutoDock-GPU 1.6
// (AD4.1_bound.dat, unbound model "bound").
enum class Ad4Type : std::int8_t {
  Merged = -1,  // non-polar hydrogen merged into its heavy atom (united atom)
  H, HD, HS, C, A, N, NA, NS, OA, OS, F, Mg, P, SA, S, Cl, Ca, Mn, Fe, Zn, Br, I, Si, B,
  Count
};

struct Ad4TypeParams {
  const char* name;
  double reqm;      // sum of vdW radii of two like atoms (Å)
  double eps;       // vdW well depth (kcal/mol)
  double reqmHb;    // H-bond radius sum contribution (Å)
  double epsHb;     // H-bond well depth (kcal/mol)
  double volume;    // atomic volume (Å³)
  double solpar;    // atomic solvation parameter
};

struct Ad4Coeffs {
  double vdW = 0.1662;
  double hb = 0.1209;
  double elec = 332.06363 * 0.1406;  // includes the Coulomb factor
  double desolv = 0.1322;
  double tors = 0.2983;
  double qasp = 0.01097;
};

const Ad4TypeParams& ad4Params(Ad4Type t);
const char* ad4TypeName(Ad4Type t);
bool ad4IsHbDonor(Ad4Type t);     // HD, HS
bool ad4IsHbAcceptor(Ad4Type t);  // NA, NS, OA, OS, SA
bool ad4IsMetal(Ad4Type t);       // Mg, Ca, Mn, Fe, Zn
// Heavy atom other than a metal: the atoms of clash checks (spec 7.9).
bool ad4Clashes(Ad4Type t);

struct Ad4Typing {
  std::vector<Ad4Type> type;    // per atom
  std::vector<double> charge;   // per atom; merged hydrogens carry 0, their charge moved to the heavy atom
  std::vector<int> unsupported; // atoms whose element has no AD4 parameters
};

// Types follow the AutoDock / Meeko rules; metal coordination bonds are not
// counted as neighbours.
Ad4Typing assignAd4Types(const Molecule& m);

// Pair parameters for the 12-6 (vdW) or 12-10 (H-bond) term:
// E = AC / s^m - BD / s^n with s the smoothed distance.
struct Ad4PairParams {
  double AC = 0, BD = 0, reqm = 0;
  int m = 12, n = 6;
};
Ad4PairParams ad4PairParams(Ad4Type a, Ad4Type b, const Ad4Coeffs& c);

}  // namespace od
