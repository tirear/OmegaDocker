// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "pka_rules.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <format>
#include <limits>
#include <map>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "error.hpp"

namespace od::pkarules {
namespace {

enum class Kind { Acid, Base };  // Acid: neutral form loses H+ (charge -1); Base: neutral form gains H+ (+1)

struct Rule {
  const char* name;
  const char* smarts;  // :1 = site atom; :2 = atom shared by the sites of a polyprotic group
  Kind kind;
  std::vector<double> pka;  // one value, or the successive pKa values of a polyprotic group
  bool beta;                // aliphatic amine: beta-substituent corrections
};

// Priority order: an atom takes the first rule that matches it. pKa values: representative compounds
// (trailing comments), aqueous, 25 °C.
const std::vector<Rule>& rules() {
  static const std::vector<Rule> r{
      {"sulfonic acid", "[OX2H1:1][SX4](=O)=O", Kind::Acid, {-1.9}, false},  // methanesulfonic acid -1.9
      {"sulfinic acid", "[OX2H1:1][SX3]=O", Kind::Acid, {1.8}, false},  // benzenesulfinic acid 1.8
      {"diphosphoric acid", "[OX2H1:1][PX4](=O)([OX2H1])[OX2:2][PX4](=O)([OX2H1])[OX2H1]", Kind::Acid, {0.9, 2.0, 6.6, 9.4}, false},  // diphosphoric acid 0.9, 2.0, 6.6, 9.4
      {"phosphoric acid", "[OX2H1:1][PX4:2](=O)([OX2H1])[OX2H1]", Kind::Acid, {2.15, 7.2, 12.35}, false},  // phosphoric acid 2.15, 7.20, 12.35
      {"phosphate monoester", "[OX2H1:1][PX4:2](=O)([OX2H1])[OX2][#6,#15]", Kind::Acid, {1.5, 6.3}, false},  // methyl phosphate 1.5, 6.3
      {"phosphonic acid", "[OX2H1:1][PX4:2](=O)([OX2H1])[#6]", Kind::Acid, {2.4, 7.7}, false},  // methylphosphonic acid 2.4, 7.7
      {"phosphate diester", "[OX2H1:1][PX4](=O)([OX2][#6,#15])[OX2][#6,#15]", Kind::Acid, {1.3}, false},  // dimethyl phosphate 1.3
      {"phosphonate monoester", "[OX2H1:1][PX4](=O)([#6])[OX2][#6,#15]", Kind::Acid, {2.3}, false},  // methyl methylphosphonate 2.3
      {"phosphinic acid", "[OX2H1:1][PX4](=O)([#6])[#6]", Kind::Acid, {3.1}, false},  // dimethylphosphinic acid 3.1
      {"thiocarboxylic acid", "[SX2H1:1][CX3]=O", Kind::Acid, {3.4}, false},  // thioacetic acid 3.4
      {"acyl sulfonamide", "[NX3H1:1]([SX4](=O)=O)[CX3]=[O,S]", Kind::Acid, {4.5}, false},  // N-acetylbenzenesulfonamide 4.5
      {"tetrazole", "[nH1:1]1nnnc1", Kind::Acid, {4.9}, false},  // tetrazole 4.9
      {"tetrazole", "[nH1:1]1nncn1", Kind::Acid, {4.9}, false},  // tetrazole 4.9
      {"alpha-amino acid carboxyl", "[OX2H1:1][CX3](=O)[CX4][$([NX3;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(Na)]),$([NX4+])]", Kind::Acid, {2.3}, false},  // glycine 2.3
      {"aromatic carboxylic acid", "[OX2H1:1][CX3](=O)a", Kind::Acid, {4.2}, false},  // benzoic acid 4.2
      {"carboxylic acid", "[OX2H1:1][CX3]=O", Kind::Acid, {4.8}, false},  // acetic acid 4.8
      {"cyclic 1,3-dione enol", "[OX2H1:1][CX3;R]=[CX3;R][CX3;R]=O", Kind::Acid, {5.0}, false},  // 1,3-cyclohexanedione 5.3, ascorbic acid 4.1
      {"2,4-dinitrophenol", "[OX2H1:1]c1c([$([NX3+](=O)[O-]),$([NX3](=O)=O)])cc([$([NX3+](=O)[O-]),$([NX3](=O)=O)])cc1", Kind::Acid, {4.1}, false},  // 2,4-dinitrophenol 4.1
      {"2,6-dinitrophenol", "[OX2H1:1]c1c([$([NX3+](=O)[O-]),$([NX3](=O)=O)])cccc1[$([NX3+](=O)[O-]),$([NX3](=O)=O)]", Kind::Acid, {3.7}, false},  // 2,6-dinitrophenol 3.7
      {"4-nitrophenol", "[OX2H1:1]c1ccc([$([NX3+](=O)[O-]),$([NX3](=O)=O)])cc1", Kind::Acid, {7.2}, false},  // 4-nitrophenol 7.2
      {"2-nitrophenol", "[OX2H1:1]c1c([$([NX3+](=O)[O-]),$([NX3](=O)=O)])cccc1", Kind::Acid, {7.2}, false},  // 2-nitrophenol 7.2
      {"phenol", "[OX2H1:1]c1ccccc1", Kind::Acid, {10.0}, false},  // phenol 10.0
      {"thiophenol", "[SX2H1:1]c", Kind::Acid, {6.6}, false},  // thiophenol 6.6
      {"thiol", "[SX2H1:1][CX4]", Kind::Acid, {10.5}, false},  // ethanethiol 10.6
      {"N-heteroaryl sulfonamide", "[NX3H1:1]([SX4](=O)=O)[c;$(c:[n,o,s])]", Kind::Acid, {6.5}, false},  // sulfadiazine 6.5
      {"N-aryl sulfonamide", "[NX3H1:1]([SX4](=O)=O)c", Kind::Acid, {8.5}, false},  // N-phenylmethanesulfonamide 8.9
      {"sulfonamide", "[NX3H2:1][SX4](=O)=O", Kind::Acid, {10.1}, false},  // benzenesulfonamide 10.1
      {"hydroxamic acid", "[OX2H1:1][NX3][CX3]=O", Kind::Acid, {9.0}, false},  // acetohydroxamic acid 9.0
      {"thiazolidinedione", "[NX3H1:1]1C(=O)S[#6]C1=O", Kind::Acid, {6.5}, false},  // rosiglitazone 6.8
      {"acyl carbamate", "[NX3H1:1]([CX3](=O)[OX2])[CX3]=O", Kind::Acid, {6.1}, false},  // dimethadione 6.1
      {"barbiturate", "[NX3H1:1]1C(=O)[#7]C(=O)[#6]C1=O", Kind::Acid, {7.8}, false},  // barbital 7.9
      {"guanine N1", "[#7X3;H1:1]1~[#6](-[#7X3;H2])~[#7X2]~[#6]2~[#7X3]~[#6]~[#7X2]~[#6]~2~[#6]~1=O", Kind::Acid, {9.2}, false},  // guanosine 9.2
      {"thymine N3", "[#7X3;H1:1]1~[#6](=O)~[#7X3]~[#6]~[#6](-[CH3])~[#6]~1=O", Kind::Acid, {9.7}, false},  // thymidine 9.7
      {"uracil N3", "[#7X3;H1:1]1~[#6](=O)~[#7X3]~[#6]~[#6]~[#6]~1=O", Kind::Acid, {9.3}, false},  // uridine 9.3
      {"imide", "[NX3H1:1]([CX3]=O)[CX3]=O", Kind::Acid, {9.6}, false},  // succinimide 9.6
      {"acylguanidine", "[NX2;+0;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(NC#N);!$(N[N+](=O)[O-]);!$(N=O);!$(NO):1]=[CX3;!a](-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)])-[NX3][CX3]=O", Kind::Base, {8.3}, false},  // acetylguanidine 8.3
      {"acylguanidine", "[NX2;+0:1]([CX3]=O)=[CX3;!a](-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)])-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)]", Kind::Base, {8.3}, false},  // acetylguanidine 8.3
      {"guanidine", "[NX2;+0;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(NC#N);!$(N[N+](=O)[O-]);!$(N=O);!$(NO):1]=[CX3;!a](-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)])-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)]", Kind::Base, {13.6}, false},  // guanidine 13.6
      {"isourea", "[NX2;+0;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(NC#N);!$(N[N+](=O)[O-]);!$(N=O);!$(NO):1]=[CX3;!a](-[OX2,SX2])-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)]", Kind::Base, {9.7}, false},  // O-methylisourea 9.7, S-methylisothiourea 9.8
      {"amidine", "[NX2;+0;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(NC#N);!$(N[N+](=O)[O-]);!$(N=O);!$(NO):1]=[CX3;!a;!$(C(=N)(-N)-[N,O,S])]-[NX3;!$(N[C,S]=[O,S]);!$(N[SX4](=O)=O);!$(N[N+](=O)[O-]);!$(N[NX3](=O)=O);!$(NC#N)]", Kind::Base, {12.4}, false},  // acetamidine 12.4
      {"4-aminopyridine", "[nX2;+0:1]1ccc(-[NX3;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O)])cc1", Kind::Base, {9.2}, false},  // 4-aminopyridine 9.2
      {"2-aminopyridine", "[nX2;+0:1]1c(-[NX3;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O)])cccc1", Kind::Base, {6.9}, false},  // 2-aminopyridine 6.9
      {"cytosine N3", "[#7X2:1]1~[#6](=O)~[#7X3]~[#6]~[#6]~[#6]~1-[#7X3;H2]", Kind::Base, {4.2}, false},  // cytidine 4.2
      {"adenine N1", "[#7X2:1]1~[#6]~[#7X2]~[#6]2~[#7X3]~[#6]~[#7X2]~[#6]~2~[#6]~1-[#7X3;H2]", Kind::Base, {3.5}, false},  // adenosine 3.5
      {"guanine N7", "[#7X2:1]1~[#6]~[#7X3]~[#6]2~[#7X2]~[#6](-[#7X3;H2])~[#7X3;H1]~[#6](=O)~[#6]~1~2", Kind::Base, {2.2}, false},  // guanosine 2.2
      {"benzimidazole", "[nX2;+0:1]1c[nX3;+0]c2ccccc12", Kind::Base, {5.5}, false},  // benzimidazole 5.5
      {"imidazole", "[nX2;+0:1]1c[nX3;+0][c;!$(c(:n):n)][c;!$(c(:n):n)]1", Kind::Base, {7.0}, false},  // imidazole 7.0
      {"pyridine", "[nX2;+0:1]1ccccc1", Kind::Base, {5.2}, false},  // pyridine 5.2
      {"alpha-amino acid amine", "[NX3;H2,H1;+0;!$(N[!#6]);!$(N[#6]=[!#6]);!$(N[#6]#*);!$(Na):1][CX4][CX3](=O)[OX2H1]", Kind::Base, {9.6}, false},  // glycine 9.6
      {"primary amine", "[NX3;H2;+0:1][CX4]", Kind::Base, {10.6}, true},  // methylamine 10.6
      {"secondary amine", "[NX3;H1;+0:1]([CX4])[CX4]", Kind::Base, {11.0}, true},  // dimethylamine 10.7, piperidine 11.1
      {"tertiary amine", "[NX3;H0;+0:1]([CX4])([CX4])[CX4]", Kind::Base, {10.0}, true},  // trimethylamine 9.8, triethylamine 10.7
      {"hydrazine", "[NX3;H2;+0:1][NX3;!$(N[C,S]=[O,S,N]);!$(N[SX4](=O)=O);!$(Na)]", Kind::Base, {7.9}, false},  // methylhydrazine 7.9
      {"aniline", "[NX3;H2;+0:1]c1ccccc1", Kind::Base, {4.6}, false},  // aniline 4.6
      {"N-alkylaniline", "[NX3;H1;+0:1]([CX4])c1ccccc1", Kind::Base, {4.9}, false},  // N-methylaniline 4.9
      {"N,N-dialkylaniline", "[NX3;H0;+0:1]([CX4])([CX4])c1ccccc1", Kind::Base, {5.1}, false},  // N,N-dimethylaniline 5.1
  };
  return r;
}

struct CompiledRule {
  const Rule* rule;
  std::unique_ptr<RDKit::RWMol> query;
  int siteIdx = -1, groupIdx = -1;
};

const std::vector<CompiledRule>& compiled() {
  static const std::vector<CompiledRule> c = [] {
    std::vector<CompiledRule> v;
    for (const Rule& r : rules()) {
      CompiledRule x{&r, std::unique_ptr<RDKit::RWMol>(RDKit::SmartsToMol(r.smarts))};
      if (!x.query) throw Error(std::format("internal error: invalid SMARTS of protonation rule \"{}\"", r.name));
      for (const RDKit::Atom* a : x.query->atoms()) {
        if (a->getAtomMapNum() == 1) x.siteIdx = static_cast<int>(a->getIdx());
        if (a->getAtomMapNum() == 2) x.groupIdx = static_cast<int>(a->getIdx());
      }
      if (x.siteIdx < 0) throw Error(std::format("internal error: protonation rule \"{}\" has no site atom", r.name));
      v.push_back(std::move(x));
    }
    return v;
  }();
  return c;
}

bool isHalogen(int z) { return z == 9 || z == 17 || z == 35 || z == 53; }

bool hasDoubleTo(const RDKit::ROMol& m, const RDKit::Atom* a, std::initializer_list<int> zs) {
  for (const RDKit::Bond* b : m.atomBonds(a)) {
    if (b->getBondType() != RDKit::Bond::DOUBLE) continue;
    const int z = b->getOtherAtom(a)->getAtomicNum();
    if (std::find(zs.begin(), zs.end(), z) != zs.end()) return true;
  }
  return false;
}

// pKa decrease of an aliphatic amine (relative to methylamine 10.6 / piperidine 11.1 / trimethylamine 9.8)
// per path, by a heteroatom X in N-C-C-X (sp3 middle carbon) or N-C-X. Several decreases are applied
// in decreasing order, each weighted 0.75 times the previous one (tris 8.1, triethanolamine 7.8).
double heteroShift(const RDKit::ROMol& m, const RDKit::Atom* x) {
  const int z = x->getAtomicNum();
  if (z == 8) return 1.1;                                  // ethanolamine 9.5, tris 8.1
  if (z == 7) return 0.9;                                  // piperazine 9.7
  if (z == 16) return hasDoubleTo(m, x, {8}) ? 2.5 : 1.0;  // 2-(methylthio)ethylamine 9.2
  if (z == 15) return 2.5;
  if (isHalogen(z)) return 1.6;                            // 2,2,2-trifluoroethylamine 5.7
  return 0.0;
}

// pKa decrease by an unsaturated carbon Y in N-C-Y.
double unsaturatedShift(const RDKit::ROMol& m, const RDKit::Atom* y) {
  if (y->getAtomicNum() != 6) return 0.0;
  if (y->getIsAromatic()) return 1.2;  // benzylamine 9.3
  for (const RDKit::Bond* b : m.atomBonds(y))
    if (b->getBondType() == RDKit::Bond::TRIPLE && b->getOtherAtom(y)->getAtomicNum() == 7) return 5.3;  // aminoacetonitrile 5.3
  if (hasDoubleTo(m, y, {8, 16, 7})) return 2.6;  // glycinamide 8.0, glycine ethyl ester 7.7
  if (hasDoubleTo(m, y, {6})) return 1.0;         // allylamine 9.5, triallylamine 8.3
  return 0.0;
}

// pKa shift by an earlier charge of the same sign, by the number of bonds between the site atoms
// (second pKa of diamines and dicarboxylic acids, e.g. ethylenediamine 9.9 / 6.9, malonic acid 2.8 / 5.7);
// half of it, in the other direction, for an opposite charge (aminomethylphosphonic acid 5.4).
double interaction(int d) {
  static constexpr std::array<double, 9> t{0, 0, 4.0, 3.0, 2.3, 1.4, 1.0, 0.8, 0.6};
  return d >= 0 && d < static_cast<int>(t.size()) ? t[d] : 0.0;
}
constexpr int kMaxInteraction = 8;

// Bond distances from `start` up to kMaxInteraction.
std::map<int, int> nearby(const RDKit::ROMol& m, int start) {
  std::map<int, int> d{{start, 0}};
  std::deque<int> q{start};
  while (!q.empty()) {
    const int a = q.front();
    q.pop_front();
    if (d[a] >= kMaxInteraction) continue;
    for (const RDKit::Atom* nb : m.atomNeighbors(m.getAtomWithIdx(static_cast<unsigned>(a)))) {
      const int b = static_cast<int>(nb->getIdx());
      if (d.contains(b)) continue;
      d[b] = d[a] + 1;
      q.push_back(b);
    }
  }
  return d;
}

bool onPhosphorus(const RDKit::ROMol& m, int a) {
  for (const RDKit::Atom* nb : m.atomNeighbors(m.getAtomWithIdx(static_cast<unsigned>(a))))
    if (nb->getAtomicNum() == 15) return true;
  return false;
}

struct Unit {
  const Rule* rule;
  std::vector<int> atoms;   // ascending
  std::vector<double> pka;  // per step
  std::size_t next = 0;     // next undecided atom
  bool done = false;
};

}  // namespace

Result protonate(const RDKit::RWMol& heavy, double pH, const std::vector<bool>& fixed,
                 const std::vector<bool>& metalBonded) {
  const int n = static_cast<int>(heavy.getNumAtoms());
  auto work = std::make_unique<RDKit::RWMol>(heavy);
  std::vector<int> h(n), q(n);
  for (int i = 0; i < n; ++i) {
    const RDKit::Atom* a = work->getAtomWithIdx(static_cast<unsigned>(i));
    h[i] = static_cast<int>(a->getTotalNumHs());
    q[i] = a->getFormalCharge();
  }
  auto apply = [&](RDKit::RWMol& m) {
    for (int i = 0; i < n; ++i) {
      RDKit::Atom* a = m.getAtomWithIdx(static_cast<unsigned>(i));
      a->setNoImplicit(true);
      a->setNumExplicitHs(h[i]);
      a->setFormalCharge(q[i]);
    }
  };

  // 1. Neutral start: charges on N, O, S, P removed by adding or removing H; a cationic aromatic N
  //    without H takes the H of a neutral N of the same ring. Charges next to an opposite charge
  //    (nitro, N-oxide) and other charges (quaternary N) stay.
  const RDKit::RingInfo* rings = work->getRingInfo();
  for (int i = 0; i < n; ++i) {
    if (fixed[i] || q[i] == 0) continue;
    const RDKit::Atom* a = work->getAtomWithIdx(static_cast<unsigned>(i));
    const int z = a->getAtomicNum();
    bool paired = false;
    for (const RDKit::Atom* nb : work->atomNeighbors(a))
      if (nb->getFormalCharge() * q[i] < 0) paired = true;
    if (paired) continue;
    if (q[i] > 0 && (z == 7 || z == 8 || z == 15 || z == 16))
      while (q[i] > 0 && h[i] > 0) --q[i], --h[i];
    if (q[i] == 1 && z == 7 && a->getIsAromatic()) {
      for (const auto& ring : rings->atomRings()) {
        if (q[i] == 0 || std::find(ring.begin(), ring.end(), i) == ring.end()) continue;
        for (int b : ring) {
          const RDKit::Atom* ab = work->getAtomWithIdx(static_cast<unsigned>(b));
          if (b != i && !fixed[b] && ab->getAtomicNum() == 7 && ab->getIsAromatic() && q[b] == 0 && h[b] > 0) {
            --h[b], q[i] = 0;
            break;
          }
        }
      }
    }
    if (q[i] < 0 && (z == 7 || z == 8 || z == 16))
      while (q[i] < 0) ++q[i], ++h[i];
  }
  apply(*work);
  try {
    RDKit::MolOps::sanitizeMol(*work);
  } catch (const std::exception& e) {
    throw Error(std::format("protonation: the neutral form of the molecule is not valid: {}", e.what()));
  }

  // 2. Sites.
  std::vector<bool> taken(n, false);
  std::vector<Unit> units;
  std::map<std::pair<const Rule*, int>, std::size_t> groupUnit;
  RDKit::SubstructMatchParameters sp;
  sp.uniquify = false;
  sp.maxMatches = std::numeric_limits<unsigned>::max();
  for (const CompiledRule& cr : compiled()) {
    for (const auto& match : RDKit::SubstructMatch(*work, *cr.query, sp)) {
      int site = -1, group = -1;
      for (const auto& [qi, mi] : match) {
        if (qi == cr.siteIdx) site = mi;
        if (qi == cr.groupIdx) group = mi;
      }
      if (site < 0 || taken[site] || fixed[site]) continue;
      if (cr.rule->kind == Kind::Acid && h[site] == 0) continue;
      taken[site] = true;
      if (group >= 0) {
        auto [it, inserted] = groupUnit.try_emplace({cr.rule, group}, units.size());
        if (inserted) units.push_back({cr.rule, {}, cr.rule->pka});
        units[it->second].atoms.push_back(site);
        continue;
      }
      double pk = cr.rule->pka[0];
      if (cr.rule->beta) {
        const RDKit::Atom* sa = work->getAtomWithIdx(static_cast<unsigned>(site));
        std::vector<double> shifts;
        for (const RDKit::Atom* alpha : work->atomNeighbors(sa))
          for (const RDKit::Atom* y : work->atomNeighbors(alpha)) {
            if (y->getIdx() == sa->getIdx()) continue;
            shifts.push_back(heteroShift(*work, y) + unsaturatedShift(*work, y));
            if (y->getAtomicNum() == 6 && y->getHybridization() == RDKit::Atom::SP3)
              for (const RDKit::Atom* x : work->atomNeighbors(y))
                if (x->getIdx() != alpha->getIdx()) shifts.push_back(heteroShift(*work, x));
          }
        std::sort(shifts.rbegin(), shifts.rend());
        double w = 1.0;
        for (double s : shifts) {
          if (s <= 0) break;
          pk -= w * s;
          w *= 0.75;
        }
      }
      units.push_back({cr.rule, {site}, {pk}});
    }
  }
  for (Unit& u : units) std::sort(u.atoms.begin(), u.atoms.end());

  // 3. Decisions, one site at a time: the undecided site with the largest |pKa - pH| first (the
  //    steps of a polyprotic group in order), with shifts from the sites charged before it.
  std::vector<std::pair<int, int>> charged;  // (atom, sign)
  std::map<int, std::map<int, int>> near;    // site atom -> bond distances
  auto shiftOf = [&](int atom, int sign) {
    double s = 0;
    if (charged.empty()) return s;
    auto it = near.find(atom);
    if (it == near.end()) it = near.emplace(atom, nearby(*work, atom)).first;
    const bool pa = onPhosphorus(*work, atom);
    for (const auto& [b, sb] : charged) {
      auto d = it->second.find(b);
      if (d == it->second.end()) continue;
      if (sb == sign && pa && onPhosphorus(*work, b)) continue;  // phosphate pKa values already include it
      s += (sb == sign ? 1.0 : -0.5) * interaction(d->second);
    }
    return s;  // > 0: the charged form is less favourable
  };
  Result res;
  auto charge = [&](int a, int sign) {
    h[a] += sign;
    q[a] = sign;
    charged.push_back({a, sign});
  };
  for (Unit& u : units) {
    const int sign = u.rule->kind == Kind::Acid ? -1 : 1;
    std::vector<int> rest;
    for (int a : u.atoms) {
      if (!metalBonded[a]) {
        rest.push_back(a);
        continue;
      }
      res.sites.push_back({u.rule->name, a, std::numeric_limits<double>::quiet_NaN(), false, true});
      if (sign < 0) {
        charge(a, -1);
        u.pka.erase(u.pka.begin());
        if (u.pka.empty()) u.pka.push_back(u.rule->pka.back());
      }
    }
    u.atoms = rest;
    u.done = u.atoms.empty();
  }
  auto stepPka = [](const Unit& u) { return u.pka[std::min(u.next, u.pka.size() - 1)]; };
  for (;;) {
    Unit* best = nullptr;
    for (Unit& u : units)
      if (!u.done && (!best || std::abs(stepPka(u) - pH) > std::abs(stepPka(*best) - pH))) best = &u;
    if (!best) break;
    Unit& u = *best;
    const int sign = u.rule->kind == Kind::Acid ? -1 : 1;
    const int a = u.atoms[u.next];
    const double pk = stepPka(u) + sign * -shiftOf(a, sign);
    const bool protonated = pk >= pH;
    res.sites.push_back({u.rule->name, a, pk, protonated, false});
    const bool isCharged = sign < 0 ? !protonated : protonated;
    if (isCharged) charge(a, sign);
    ++u.next;
    if (u.next == u.atoms.size()) u.done = true;
    if (!isCharged && sign < 0) {
      // the remaining steps of the group have larger pKa values: protonated
      for (; u.next < u.atoms.size(); ++u.next) res.sites.push_back({u.rule->name, u.atoms[u.next], stepPka(u), true, false});
      u.done = true;
    }
  }
  std::sort(res.sites.begin(), res.sites.end(), [](const Site& x, const Site& y) { return x.atom < y.atom; });

  apply(*work);
  try {
    RDKit::MolOps::sanitizeMol(*work);
    RDKit::MolOps::Kekulize(*work, true);
  } catch (const std::exception& e) {
    throw Error(std::format("protonation: the protonated molecule is not valid: {}", e.what()));
  }
  res.mol = std::move(work);
  return res;
}

}  // namespace od::pkarules
