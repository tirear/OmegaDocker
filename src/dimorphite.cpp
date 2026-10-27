// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// C++ port of Dimorphite-DL 2.1.0 (Apache-2.0, https://github.com/durrantlab/dimorphite_dl,
// commit 69ec38e): neutralize.py, protonate/detect.py, protonate/change.py.

#include "dimorphite.hpp"

#include <format>
#include <map>
#include <functional>
#include <sstream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "log.hpp"

extern "C" const char od_dimorphite_smarts_begin[];
extern "C" const char od_dimorphite_smarts_end[];

namespace od::dimorphite {
namespace {

struct PkaDatum {
  int site;
  double mean;
};
struct Substructure {
  std::string name;
  std::unique_ptr<RDKit::RWMol> query;
  std::vector<PkaDatum> pkas;
};

std::vector<Substructure> loadSubstructures() {
  std::vector<Substructure> v;
  std::istringstream in(std::string(od_dimorphite_smarts_begin, od_dimorphite_smarts_end));
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::vector<std::string> w;
    for (std::string x; ls >> x;) w.push_back(x);
    if (w.empty() || w[0][0] == '#') continue;
    if (w.size() < 3 || (w.size() - 2) % 3 != 0) fail("site_substructures.smarts: invalid line '{}'", line);
    Substructure s;
    s.name = w[0];
    s.query.reset(RDKit::SmartsToMol(w[1]));
    if (!s.query) fail("site_substructures.smarts: invalid SMARTS '{}'", w[1]);
    for (std::size_t i = 2; i < w.size(); i += 3) s.pkas.push_back({std::stoi(w[i]), std::stod(w[i + 1])});
    v.push_back(std::move(s));
  }
  return v;
}

const std::vector<Substructure>& substructures() {
  static const std::vector<Substructure> v = loadSubstructures();
  return v;
}

std::vector<std::vector<int>> matches(const RDKit::ROMol& mol, const RDKit::ROMol& q) {
  RDKit::SubstructMatchParameters p;
  p.maxMatches = 1u << 30;
  std::vector<std::vector<int>> out;
  for (const auto& m : RDKit::SubstructMatch(mol, q, p)) {
    std::vector<int> v(m.size());
    for (const auto& [qi, mi] : m) v[qi] = mi;
    out.push_back(std::move(v));
  }
  return out;
}

// ---------- neutralize.py ----------

struct Reaction {
  std::unique_ptr<RDKit::RWMol> query;
  std::function<void(RDKit::RWMol&, const std::vector<int>&)> apply;
};

void addH(RDKit::RWMol& m, int heavy) {
  const unsigned h = m.addAtom(new RDKit::Atom(1), false, true);
  m.addBond(static_cast<unsigned>(heavy), h, RDKit::Bond::SINGLE);
}
void setBond(RDKit::RWMol& m, int a, int b, RDKit::Bond::BondType t) {
  RDKit::Bond* bd = m.getBondBetweenAtoms(static_cast<unsigned>(a), static_cast<unsigned>(b));
  bd->setBondType(t);
  bd->setIsAromatic(false);
}
void charge(RDKit::RWMol& m, int a, int q) { m.getAtomWithIdx(static_cast<unsigned>(a))->setFormalCharge(q); }

std::vector<Reaction> makeReactions() {
  using B = RDKit::Bond;
  std::vector<std::pair<const char*, std::function<void(RDKit::RWMol&, const std::vector<int>&)>>> d{
      {"[#16+2:1]-[Ov1-1:2]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 1); charge(m, x[1], 0); setBond(m, x[0], x[1], B::DOUBLE); }},
      {"[#16+1:1]-[Ov1-1:2]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); charge(m, x[1], 0); setBond(m, x[0], x[1], B::DOUBLE); }},
      {"[#15+1:1]-[Ov1-1:2]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); charge(m, x[1], 0); setBond(m, x[0], x[1], B::DOUBLE); }},
      {"[Ov1-1;!$([O-]-[#7+;!$([#7+]=O)]);!$([O-]-[#7+](=O)-[#8]-[H]):1]",
       [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); addH(m, x[0]); }},
      {"[Sv1-1:1]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); addH(m, x[0]); }},
      {"[#7v4+1:1]-[H]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); m.removeAtom(static_cast<unsigned>(x[1])); }},
      {"[Ov2-:1]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); }},
      {"[#7v3+1:1]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); }},
      {"[#7v2-1;!$([#7-]=[#7+]=[#6]):1]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); addH(m, x[0]); }},
      {"[H]-[N:1]-[N:2]#[N:3]",
       [](RDKit::RWMol& m, const std::vector<int>& x) {
         setBond(m, x[1], x[2], B::DOUBLE);
         setBond(m, x[2], x[3], B::DOUBLE);
         charge(m, x[2], 1);
         addH(m, x[3]);
         m.removeAtom(static_cast<unsigned>(x[0]));
       }},
      {"[#7+1;H0;!$([#7+]~[O-]):1]=[#6:2]-[#7+0;!H0:3]",
       [](RDKit::RWMol& m, const std::vector<int>& x) {
         charge(m, x[0], 0);
         charge(m, x[2], 1);
         setBond(m, x[0], x[1], B::SINGLE);
         setBond(m, x[1], x[2], B::DOUBLE);
       }},
      {"[n+1;H0;!$([n+]~[O-]):1]:[c:2]:[n+0;!H0:3]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); charge(m, x[2], 1); }},
      {"[n+1;H0;!$([n+]~[O-]):1]:[n+0;!H0:2]", [](RDKit::RWMol& m, const std::vector<int>& x) { charge(m, x[0], 0); charge(m, x[1], 1); }},
  };
  std::vector<Reaction> r;
  for (auto& [s, f] : d) {
    Reaction x;
    x.query.reset(RDKit::SmartsToMol(s));
    x.apply = f;
    r.push_back(std::move(x));
  }
  return r;
}

const std::vector<Reaction>& reactions() {
  static const std::vector<Reaction> r = makeReactions();
  return r;
}

void neutralize(RDKit::RWMol& m, const std::vector<bool>& fixed) {
  m.updatePropertyCache(false);
  const std::size_t maxPasses = 10 * m.getNumAtoms();
  std::size_t passes = 0;
  for (bool changed = true; changed;) {
    changed = false;
    for (const Reaction& r : reactions()) {
      for (const auto& x : matches(m, *r.query)) {
        bool touchesFixed = false;
        for (int a : x)
          if (a < static_cast<int>(fixed.size()) && fixed[a]) touchesFixed = true;
        if (touchesFixed) continue;
        if (++passes > maxPasses) fail("Dimorphite-DL neutralization did not converge after {} passes", maxPasses);
        r.apply(m, x);
        m.updatePropertyCache(false);
        changed = true;
        break;
      }
      if (changed) break;
    }
  }
}

// ---------- change.py ----------

double bondOrderTotal(const RDKit::ROMol& m, const RDKit::Atom* a) {
  double t = 0;
  for (const auto* b : m.atomBonds(a)) t += b->getBondTypeAsDouble();
  return t;
}

enum class Applied { Ok, Invalid, Stop };

// _apply_charge_to_molecule. Invalid: no valid state (molecule unchanged); Stop: the molecule
// cannot be sanitized (RemoveHs fails), later sites are not applied.
Applied applyCharge(RDKit::RWMol& m, int idx, int q, int nq) {
  RDKit::RWMol copy(m);
  try {
    RDKit::MolOps::sanitizeMol(copy);
  } catch (const std::exception&) {
    return Applied::Stop;
  }
  const RDKit::RWMol parent(copy);
  RDKit::Atom* a = copy.getAtomWithIdx(static_cast<unsigned>(idx));
  const double bo = bondOrderTotal(copy, a);
  if (a->getAtomicNum() == 7) {
    if (!(nq == 1 && bo == 4.0 && a->getIsAromatic() && a->getDegree() == 3)) {
      a->setFormalCharge(nq);
      static const std::map<std::pair<int, int>, int> h{{{1, 1}, 3}, {{1, 2}, 2}, {{1, 3}, 1}, {{0, 1}, 2},
                                                        {{0, 2}, 1}, {{-1, 1}, 1}, {{-1, 2}, 0}};
      if (auto it = h.find({nq, static_cast<int>(bo)}); it != h.end()) a->setNumExplicitHs(static_cast<unsigned>(it->second));
    }
    if (a->getIsAromatic() && a->getFormalCharge() == -1 && a->getNumExplicitHs() > 0) a->setNumExplicitHs(0);
  } else {
    a->setFormalCharge(q);
    if (a->getAtomicNum() == 8 || a->getAtomicNum() == 16) {
      if (q == 0 && bo == 1) a->setNumExplicitHs(1);
      else if (q == -1 && bo == 1) a->setNumExplicitHs(0);
    }
  }
  copy.updatePropertyCache(false);
  RDKit::RWMol test(copy);
  unsigned failed = 0;
  try {
    RDKit::MolOps::sanitizeMol(test, failed);
  } catch (const std::exception&) {
    failed = 1;
  }
  if (failed != 0) {
    m = parent;
    return Applied::Invalid;
  }
  m = std::move(copy);
  return Applied::Ok;
}

}  // namespace

Result protonate(const RDKit::RWMol& heavy, double pH, const std::vector<bool>& fixed, const std::vector<bool>& metalBonded) {
  const int n = static_cast<int>(heavy.getNumAtoms());
  // MoleculeNeutralizer: add H, neutralize, sanitize, remove H.
  RDKit::RWMol work(heavy);
  RDKit::MolOps::addHs(work);
  neutralize(work, fixed);
  RDKit::MolOps::sanitizeMol(work);
  RDKit::MolOps::removeHs(work);
  if (static_cast<int>(work.getNumAtoms()) != n) fail("Dimorphite-DL: heavy atoms changed during neutralization");

  // ProtonationSiteDetector on the neutral molecule with explicit hydrogens.
  RDKit::RWMol det(work);
  RDKit::MolOps::addHs(det);
  std::vector<bool> prot(det.getNumAtoms(), false);
  for (int i = 0; i < n; ++i) prot[i] = fixed[i];
  struct Found {
    const Substructure* s;
    std::vector<int> match;
  };
  std::vector<Found> sites;
  for (const Substructure& s : substructures()) {
    std::vector<std::vector<int>> used;
    for (const auto& x : matches(det, *s.query)) {
      bool free = true;
      for (int a : x)
        if (prot[a]) free = false;
      if (!free) continue;
      sites.push_back({&s, x});
      used.push_back(x);
      for (const PkaDatum& p : s.pkas) {
        const int a = x[p.site];
        prot[a] = true;
        for (const auto* nb : det.atomNeighbors(det.getAtomWithIdx(static_cast<unsigned>(a))))
          if (nb->getAtomicNum() == 1) prot[nb->getIdx()] = true;
      }
    }
    for (const auto& x : used) {
      std::vector<int> siteAtoms;
      for (const PkaDatum& p : s.pkas) siteAtoms.push_back(x[p.site]);
      for (int a : x)
        if (std::find(siteAtoms.begin(), siteAtoms.end(), a) == siteAtoms.end() &&
            det.getAtomWithIdx(static_cast<unsigned>(a))->getAtomicNum() != 6)
          prot[a] = true;
    }
  }

  // protonate_site for each site in detection order (one state per pKa datum).
  Result res;
  RDKit::RWMol cur(work);
  bool stop = false;
  for (const Found& f : sites) {
    if (stop) break;
    const bool special = !f.s->name.empty() && f.s->name[0] == '*';
    std::vector<std::pair<int, bool>> states;  // (atom, protonated), unique in first-seen order
    std::vector<double> means;
    for (const PkaDatum& p : f.s->pkas) {
      const int a = f.match[p.site];
      if (a >= n) continue;
      bool protonated = p.mean >= pH;
      const bool forced = a < n && metalBonded[a];
      if (forced) protonated = false;
      if (std::find(states.begin(), states.end(), std::pair<int, bool>{a, protonated}) == states.end()) {
        states.push_back({a, protonated});
        res.sites.push_back({special ? f.s->name.substr(1) : f.s->name, a, p.mean, protonated, forced});
      }
    }
    for (const auto& [a, protonated] : states) {
      const int q = protonated ? 0 : -1;
      const int nq = special ? q : q + 1;
      const Applied r = applyCharge(cur, a, q, nq);
      if (r == Applied::Invalid)
        logWarn(std::format("Dimorphite-DL: no valid state with charge {} at atom {} ({}); the site is left unchanged", q, a + 1,
                            f.s->name));
      if (r == Applied::Stop) {
        logWarn("Dimorphite-DL: the molecule cannot be sanitized; the remaining sites are left unchanged");
        stop = true;
        break;
      }
    }
  }
  try {
    RDKit::MolOps::sanitizeMol(cur);
    RDKit::MolOps::Kekulize(cur, false);
  } catch (const std::exception& e) {
    fail("Dimorphite-DL: the protonated molecule is not valid: {}", e.what());
  }
  res.mol = std::make_unique<RDKit::RWMol>(cur);
  return res;
}

}  // namespace od::dimorphite
