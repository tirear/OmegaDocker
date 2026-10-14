// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "molecule.hpp"

#include <algorithm>
#include <cctype>
#include <numeric>

namespace od {
namespace {

constexpr const char* kSymbols[] = {
    "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
constexpr int kNumElements = static_cast<int>(std::size(kSymbols)) - 1;

struct DisjointSet {
  std::vector<int> p;
  explicit DisjointSet(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
  int find(int x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  }
  void unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) p[std::max(a, b)] = std::min(a, b);
  }
};

}  // namespace

const char* elementSymbol(int z) { return (z >= 1 && z <= kNumElements) ? kSymbols[z] : "X"; }

int atomicNumber(const std::string& symbol) {
  if (symbol.empty() || symbol.size() > 3) return 0;
  if (symbol == "D" || symbol == "T" || symbol == "d" || symbol == "t") return 1;  // deuterium, tritium
  std::string s;
  s += static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
  for (std::size_t i = 1; i < symbol.size(); ++i)
    s += static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[i])));
  for (int z = 1; z <= kNumElements; ++z)
    if (s == kSymbols[z]) return z;
  return 0;
}

bool isMetal(int z) {
  switch (z) {
    case 3: case 4: case 11: case 12: case 13: case 19: case 20:
    case 31: case 37: case 38: case 49: case 50: case 55: case 56:
    case 81: case 82: case 83: case 84: case 87: case 88:
      return true;
    default:
      return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 57 && z <= 80) ||
             (z >= 89 && z <= 112);
  }
}

bool isMetalBond(const Molecule& m, const Bond& b) {
  return b.order == BondOrder::Metal || isMetal(m.atoms[b.a].atomicNum) || isMetal(m.atoms[b.b].atomicNum);
}

const char* roleTag(Role r) { return r == Role::StaM ? "sta" : "mob"; }

std::vector<Molecule> splitMolecules(const RawStructure& s) {
  const int n = static_cast<int>(s.atoms.size());
  DisjointSet ds(n);
  for (const Bond& b : s.bonds) ds.unite(b.a, b.b);

  std::vector<int> compOfRoot(n, -1), comp(n), local(n);
  std::vector<Molecule> mols;
  for (int i = 0; i < n; ++i) {
    int r = ds.find(i);
    if (compOfRoot[r] < 0) {
      compOfRoot[r] = static_cast<int>(mols.size());
      Molecule m;
      m.sourceFile = s.sourceFile;
      m.hasResidueInfo = s.hasResidueInfo;
      m.chargesFromFile = s.hasPartialCharges;
      mols.push_back(std::move(m));
    }
    comp[i] = compOfRoot[r];
    local[i] = static_cast<int>(mols[comp[i]].atoms.size());
    Atom a = s.atoms[i];
    a.residue = -1;
    mols[comp[i]].atoms.push_back(std::move(a));
  }
  for (const Bond& b : s.bonds) {
    Bond nb = b;
    nb.a = local[b.a];
    nb.b = local[b.b];
    mols[comp[b.a]].bonds.push_back(nb);
  }
  // A residue split over several components is kept as one residue per component.
  for (const Residue& res : s.residues) {
    std::vector<int> resIndexOfComp;
    std::vector<int> seenComps;
    for (int ai : res.atoms) {
      int c = comp[ai];
      auto it = std::find(seenComps.begin(), seenComps.end(), c);
      int ri;
      if (it == seenComps.end()) {
        seenComps.push_back(c);
        ri = static_cast<int>(mols[c].residues.size());
        resIndexOfComp.push_back(ri);
        mols[c].residues.push_back(Residue{res.info, {}});
      } else {
        ri = resIndexOfComp[it - seenComps.begin()];
      }
      mols[c].residues[ri].atoms.push_back(local[ai]);
      mols[c].atoms[local[ai]].residue = ri;
    }
  }
  return mols;
}

}  // namespace od
