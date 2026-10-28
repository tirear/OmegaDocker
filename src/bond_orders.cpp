// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "bond_orders.hpp"

#include <algorithm>
#include <format>
#include <functional>
#include <optional>

namespace od {
namespace {

// Formal charge implied by the bond-order sum `v`; nullopt if not supported.
std::optional<int> chargeFromValence(int z, int v) {
  auto pick = [v](std::initializer_list<std::pair<int, int>> table) -> std::optional<int> {
    for (auto [valence, charge] : table)
      if (valence == v) return charge;
    return std::nullopt;
  };
  switch (z) {
    case 1:  return pick({{1, 0}});
    case 5:  return pick({{3, 0}, {4, -1}});
    case 6:  return pick({{4, 0}});
    case 7:  return pick({{3, 0}, {4, 1}, {2, -1}});
    case 8:  return pick({{2, 0}, {1, -1}});
    case 14: return pick({{4, 0}});
    case 15:
    case 33: return pick({{3, 0}, {5, 0}, {4, 1}});
    case 16:
    case 34: return pick({{2, 0}, {4, 0}, {6, 0}, {1, -1}, {3, 1}});
    case 9:
    case 17:
    case 35:
    case 53: return pick({{1, 0}, {0, -1}});
    default: return std::nullopt;
  }
}

// Valences outside the table accepted when the charge is given (CCD or pKa rules):
// carbanion / carbocation (cyanide, isocyanide, carbon monoxide C) and oxonium (carbon monoxide O).
std::optional<int> chargeFromValenceGiven(int z, int v, int given) {
  if (z == 6 && v == 3 && (given == -1 || given == 1)) return given;
  if (z == 8 && v == 3 && given == 1) return given;
  return std::nullopt;
}

struct AromaticSystem {
  std::vector<int> atoms;
  std::vector<int> bonds;
};

// Perfect matching of `needy` atoms over the bonds of the system (backtracking).
bool matchSystem(const AromaticSystem& sys, const std::vector<Bond>& bonds,
                 const std::vector<int>& needy, std::vector<int>& doubleBonds) {
  std::vector<std::vector<std::pair<int, int>>> adj(needy.size());  // (other local, bond)
  auto localOf = [&](int atom) {
    auto it = std::find(needy.begin(), needy.end(), atom);
    return it == needy.end() ? -1 : static_cast<int>(it - needy.begin());
  };
  for (int bi : sys.bonds) {
    int la = localOf(bonds[bi].a), lb = localOf(bonds[bi].b);
    if (la >= 0 && lb >= 0) {
      adj[la].push_back({lb, bi});
      adj[lb].push_back({la, bi});
    }
  }
  std::vector<int> mate(needy.size(), -1);
  doubleBonds.clear();
  std::function<bool()> solve = [&]() -> bool {
    int best = -1;
    std::size_t bestDeg = SIZE_MAX;
    for (std::size_t i = 0; i < needy.size(); ++i) {
      if (mate[i] >= 0) continue;
      std::size_t deg = 0;
      for (auto [j, bi] : adj[i]) deg += (mate[j] < 0);
      if (deg < bestDeg) { bestDeg = deg; best = static_cast<int>(i); }
    }
    if (best < 0) return true;
    for (auto [j, bi] : adj[best]) {
      if (mate[j] >= 0) continue;
      mate[best] = j; mate[j] = best;
      doubleBonds.push_back(bi);
      if (solve()) return true;
      doubleBonds.pop_back();
      mate[best] = -1; mate[j] = -1;
    }
    return false;
  };
  return solve();
}

}  // namespace

std::string describeAtom(const RawStructure& s, int i) {
  const Atom& a = s.atoms[i];
  if (a.residue >= 0) {
    const ResidueInfo& r = s.residues[a.residue].info;
    return std::format("{}:{}{}:{}:{}", r.chain, r.seq, r.iCode, r.name, a.name);
  }
  return std::format("#{}", a.fileSerial);
}

void assignChargesAndKekule(RawStructure& s, const std::vector<bool>& aromaticCandidate,
                            const std::vector<int>& ccdCharge, std::vector<std::string>& errors) {
  const int n = static_cast<int>(s.atoms.size());
  std::vector<int> fixedValence(n, 0), aromaticDegree(n, 0), degree(n, 0);
  std::vector<std::vector<int>> aromaticBondsOf(n);
  for (std::size_t bi = 0; bi < s.bonds.size(); ++bi) {
    const Bond& b = s.bonds[bi];
    if (b.order == BondOrder::Metal) continue;
    ++degree[b.a];
    ++degree[b.b];
    if (aromaticCandidate[bi]) {
      ++aromaticDegree[b.a];
      ++aromaticDegree[b.b];
      aromaticBondsOf[b.a].push_back(static_cast<int>(bi));
      aromaticBondsOf[b.b].push_back(static_cast<int>(bi));
    } else {
      int o = static_cast<int>(b.order);
      fixedValence[b.a] += o;
      fixedValence[b.b] += o;
    }
  }

  // Atoms outside conjugated systems.
  for (int i = 0; i < n; ++i) {
    Atom& a = s.atoms[i];
    if (aromaticDegree[i] > 0) continue;
    if (isMetal(a.atomicNum) || degree[i] == 0) {
      a.formalCharge = ccdCharge[i];
      continue;
    }
    auto q = chargeFromValence(a.atomicNum, fixedValence[i]);
    if (!q) q = chargeFromValenceGiven(a.atomicNum, fixedValence[i], a.formalCharge != 0 ? a.formalCharge : ccdCharge[i]);
    if (!q) {
      errors.push_back(std::format("{}: element {} with bond-order sum {} is not supported "
                                   "(missing atoms or hydrogens? SDF/mol2 input can be used instead)",
                                   describeAtom(s, i), elementSymbol(a.atomicNum), fixedValence[i]));
      continue;
    }
    a.formalCharge = *q;
  }

  // Conjugated systems (connected through CCD aromatic bonds).
  std::vector<int> sysOf(n, -1);
  std::vector<AromaticSystem> systems;
  for (int i = 0; i < n; ++i) {
    if (aromaticDegree[i] == 0 || sysOf[i] >= 0) continue;
    AromaticSystem sys;
    std::vector<int> stack{i};
    sysOf[i] = static_cast<int>(systems.size());
    while (!stack.empty()) {
      int x = stack.back();
      stack.pop_back();
      sys.atoms.push_back(x);
      for (int bi : aromaticBondsOf[x]) {
        if (std::find(sys.bonds.begin(), sys.bonds.end(), bi) == sys.bonds.end()) sys.bonds.push_back(bi);
        int y = s.bonds[bi].a == x ? s.bonds[bi].b : s.bonds[bi].a;
        if (sysOf[y] < 0) {
          sysOf[y] = sysOf[i];
          stack.push_back(y);
        }
      }
    }
    systems.push_back(std::move(sys));
  }

  for (const AromaticSystem& sys : systems) {
    // need: number of extra bond orders (0 or 1) the atom must receive in the system.
    std::vector<int> needy, nPlusCandidates, nMinus;
    bool ok = true;
    for (int x : sys.atoms) {
      int base = fixedValence[x] + aromaticDegree[x];
      int z = s.atoms[x].atomicNum;
      const int given = s.atoms[x].formalCharge != 0 ? s.atoms[x].formalCharge : ccdCharge[x];
      int target = (z == 6) ? 4 : (z == 7) ? 3 : (z == 8 || z == 16 || z == 34) ? 2 : -1;
      if (target < 0 || base > target + (z == 7 ? 1 : 0) || target - base > 1) {
        errors.push_back(std::format("{}: unsupported conjugated-ring environment",
                                     describeAtom(s, x)));
        ok = false;
        continue;
      }
      s.atoms[x].formalCharge = 0;
      if (z == 7 && base == 4) {
        s.atoms[x].formalCharge = 1;  // e.g. N-substituted pyridinium with explicit H
      } else if (z == 7 && base == 2 && given == -1) {
        nMinus.push_back(x);  // deprotonated aromatic N (e.g. tetrazolate)
      } else if (target - base == 1) {
        needy.push_back(x);
      } else if (z == 7 && base == 3) {
        nPlusCandidates.push_back(x);
      }
    }
    if (!ok) continue;

    std::vector<int> doubles;
    bool solved = matchSystem(sys, s.bonds, needy, doubles);
    if (solved) {
      for (int x : nMinus) s.atoms[x].formalCharge = -1;
    } else if (!nMinus.empty()) {
      // The given anionic N is not compatible; treat it as pyridine-type N.
      needy.insert(needy.end(), nMinus.begin(), nMinus.end());
      solved = matchSystem(sys, s.bonds, needy, doubles);
    }
    // Promote pyrrole-type N to N+ (needing one double bond) until a pattern exists.
    for (std::size_t k = 1; !solved && k <= std::min<std::size_t>(3, nPlusCandidates.size()); ++k) {
      std::vector<bool> pick(nPlusCandidates.size(), false);
      std::fill(pick.begin(), pick.begin() + static_cast<long>(k), true);
      do {
        std::vector<int> trial = needy;
        for (std::size_t c = 0; c < pick.size(); ++c)
          if (pick[c]) trial.push_back(nPlusCandidates[c]);
        if (trial.size() % 2 == 0 && matchSystem(sys, s.bonds, trial, doubles)) {
          for (std::size_t c = 0; c < pick.size(); ++c)
            if (pick[c]) s.atoms[nPlusCandidates[c]].formalCharge = 1;
          solved = true;
        }
      } while (!solved && std::prev_permutation(pick.begin(), pick.end()));
    }
    if (!solved) {
      errors.push_back(std::format("{}: no valid single/double bond pattern for the conjugated "
                                   "system with the hydrogens present",
                                   describeAtom(s, sys.atoms.front())));
      continue;
    }
    for (int bi : sys.bonds) s.bonds[bi].order = BondOrder::Single;
    for (int bi : doubles) s.bonds[bi].order = BondOrder::Double;
  }
}

}  // namespace od
