// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "protonate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>

#include "bond_orders.hpp"
#include "dimorphite.hpp"
#include "log.hpp"
#include "propka.hpp"

namespace od {
namespace {

// ---------- geometry ----------

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 unit(const Vec3& a) {
  const double n = norm(a);
  return n > 1e-12 ? mul(a, 1.0 / n) : Vec3{1, 0, 0};
}
double dist(const Vec3& a, const Vec3& b) { return norm(sub(a, b)); }
Vec3 perpendicular(const Vec3& u) {
  const Vec3 t = std::abs(u[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  return unit(cross(u, t));
}
// Position d with |cd| = len, angle bcd = ang, dihedral abcd = tor (degrees).
Vec3 placeTorsion(const Vec3& a, const Vec3& b, const Vec3& c, double len, double ang, double tor) {
  const double A = ang * M_PI / 180.0, T = tor * M_PI / 180.0;
  const Vec3 bc = unit(sub(c, b));
  Vec3 n = cross(sub(b, a), bc);
  n = norm(n) < 1e-9 ? perpendicular(bc) : unit(n);
  const Vec3 m = cross(n, bc);
  const Vec3 d2{-len * std::cos(A), len * std::sin(A) * std::cos(T), len * std::sin(A) * std::sin(T)};
  return add(c, add(mul(bc, d2[0]), add(mul(m, d2[1]), mul(n, d2[2]))));
}

const std::set<std::string> kStdAmino{"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                                      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

double bondLength(int z) {
  switch (z) {
    case 6: return 1.09;
    case 7: return 1.01;
    case 8: return 0.96;
    case 16: return 1.34;
    case 5: return 1.19;
    case 15: return 1.42;
    case 14: return 1.48;
    case 34: return 1.47;
    default: return 1.0;
  }
}

// ---------- structure editing ----------

// Rebuilds atoms in `order` (old indices); bonds whose atoms are kept are remapped.
void reorder(RawStructure& s, const std::vector<int>& order, std::vector<int>* q, std::vector<bool>* arom) {
  std::vector<int> map(s.atoms.size(), -1);
  for (std::size_t k = 0; k < order.size(); ++k) map[order[k]] = static_cast<int>(k);
  std::vector<Atom> atoms;
  std::vector<int> nq;
  for (int i : order) {
    atoms.push_back(s.atoms[i]);
    if (q) nq.push_back((*q)[i]);
  }
  std::vector<Bond> bonds;
  std::vector<bool> na;
  for (std::size_t b = 0; b < s.bonds.size(); ++b) {
    Bond x = s.bonds[b];
    if (map[x.a] < 0 || map[x.b] < 0) continue;
    x.a = map[x.a];
    x.b = map[x.b];
    bonds.push_back(x);
    if (arom) na.push_back((*arom)[b]);
  }
  for (Residue& r : s.residues) {
    std::vector<int> v;
    for (int i : r.atoms)
      if (map[i] >= 0) v.push_back(map[i]);
    r.atoms = std::move(v);
  }
  // Residue atom lists in structure order.
  for (Residue& r : s.residues) std::sort(r.atoms.begin(), r.atoms.end());
  s.atoms = std::move(atoms);
  s.bonds = std::move(bonds);
  if (q) *q = std::move(nq);
  if (arom) *arom = std::move(na);
}

std::vector<std::vector<int>> adjacency(const RawStructure& s) {
  std::vector<std::vector<int>> adj(s.atoms.size());
  for (const Bond& b : s.bonds) {
    adj[b.a].push_back(b.b);
    adj[b.b].push_back(b.a);
  }
  return adj;
}

int findInResidue(const RawStructure& s, int res, const std::string& name) {
  if (res < 0) return -1;
  for (int i : s.residues[res].atoms)
    if (s.atoms[i].name == name) return i;
  return -1;
}

bool isMetalBondOf(const RawStructure& s, const Bond& b) {
  return b.order == BondOrder::Metal || isMetal(s.atoms[b.a].atomicNum) || isMetal(s.atoms[b.b].atomicNum);
}

std::string resLabel(const RawStructure& s, int r) {
  const ResidueInfo& i = s.residues[r].info;
  return std::format("{}:{}{}:{}", i.chain, i.seq, i.iCode, i.name);
}

// Heavy-atom RDKit molecule: implicit hydrogens, no metal bonds.
std::unique_ptr<RDKit::RWMol> heavyMol(const RawStructure& s, const std::vector<int>& charge) {
  auto m = std::make_unique<RDKit::RWMol>();
  for (std::size_t i = 0; i < s.atoms.size(); ++i) {
    RDKit::Atom a(s.atoms[i].atomicNum);
    a.setFormalCharge(charge[i]);
    m->addAtom(&a, true, false);
  }
  for (const Bond& b : s.bonds) {
    if (isMetalBondOf(s, b)) continue;
    const auto t = b.order == BondOrder::Double ? RDKit::Bond::DOUBLE : b.order == BondOrder::Triple ? RDKit::Bond::TRIPLE : RDKit::Bond::SINGLE;
    m->addBond(static_cast<unsigned>(b.a), static_cast<unsigned>(b.b), t);
  }
  return m;
}

// ---------- hydrogen placement and optimization ----------

struct AtomView {
  int idx;
  Vec3 p;
  int parent;  // hydrogens
  bool polarH, acceptor, acceptorOnly;
  double radius;
  Vec3 pp{};   // hydrogens: parent position in the same state
  bool donorHeavy = false;  // N / O carrying hydrogens
};

struct GroupState {
  std::vector<std::pair<int, Vec3>> pos;    // moving atoms
  std::vector<std::pair<int, int>> parent;  // (hydrogen, parent) for movable hydrogens
  bool flipped = false;
  double bias = 0;
  double tor = 0;  // free rotors
};

struct OptGroup {
  std::string label;
  std::vector<int> moving;
  std::vector<GroupState> states;
  int chosen = 0;
  bool rotor = false;
  int rotorParent = -1;  // free rotors (refined torsion)
};

class Placer {
public:
  Placer(RawStructure& s, std::vector<int> hyb) : s_(s), hyb_(std::move(hyb)) {
    adj_ = adjacency(s);
    xyz_.resize(s.atoms.size());
    for (std::size_t i = 0; i < s.atoms.size(); ++i) xyz_[i] = s.atoms[i].xyz;
  }

  bool isH(int i) const { return s_.atoms[i].atomicNum == 1; }
  std::vector<int> heavyNbrs(int p) const {
    std::vector<int> v;
    for (int n : adj_[p])
      if (!isH(n)) v.push_back(n);
    return v;
  }
  std::vector<int> hydrogens(int p) const {
    std::vector<int> v;
    for (int n : adj_[p])
      if (isH(n)) v.push_back(n);
    return v;
  }
  int refAtom(int p, int n) const {
    for (int r : adj_[n])
      if (r != p && !isH(r)) return r;
    return -1;
  }
  double hAngle(int p) const {
    const int z = s_.atoms[p].atomicNum;
    if (z == 16 || z == 34) return 96.0;
    if (hyb_[p] == 2 && z != 8) return 120.0;
    return 109.47;
  }

  // Positions of h hydrogens on p; rotor torsion `tor` (k = 1), water orientation `ori` (k = 0).
  std::vector<Vec3> place(int p, int h, double tor, int ori = 0) const {
    std::vector<Vec3> out;
    const Vec3 c = xyz_[p];
    const double len = bondLength(s_.atoms[p].atomicNum);
    const std::vector<int> nb = heavyNbrs(p);
    const std::size_t k = nb.size();
    if (h <= 0) return out;
    if (k == 0) {
      static const std::vector<Vec3> ico = [] {
        const double g = (1 + std::sqrt(5.0)) / 2;
        std::vector<Vec3> v{{-1, g, 0}, {1, g, 0}, {-1, -g, 0}, {1, -g, 0}, {0, -1, g}, {0, 1, g},
                            {0, -1, -g}, {0, 1, -g}, {g, 0, -1}, {g, 0, 1}, {-g, 0, -1}, {-g, 0, 1}};
        for (Vec3& x : v) x = unit(x);
        return v;
      }();
      const Vec3 d1 = ico[(ori / 3) % 12];
      out.push_back(add(c, mul(d1, len)));
      if (h >= 2) {
        const Vec3 a = add(c, perpendicular(d1));
        out.push_back(placeTorsion(a, add(c, d1), c, len, 104.5, 120.0 * (ori % 3)));
      }
      for (int j = 2; j < h; ++j) out.push_back(placeTorsion(out[0], add(c, d1), c, len, 109.47, 120.0 * (j - 1)));
      return out;
    }
    if (k == 1) {
      const int n = nb[0];
      const int r = refAtom(p, n);
      const Vec3 a = r >= 0 ? xyz_[r] : add(xyz_[n], perpendicular(unit(sub(c, xyz_[n]))));
      if (hyb_[p] == 1) {
        out.push_back(add(c, mul(unit(sub(c, xyz_[n])), len)));
        return out;
      }
      const double ang = hAngle(p);
      if (hyb_[p] == 2 && h == 2) {
        out.push_back(placeTorsion(a, xyz_[n], c, len, ang, 0.0));
        out.push_back(placeTorsion(a, xyz_[n], c, len, ang, 180.0));
        return out;
      }
      for (int j = 0; j < h; ++j) out.push_back(placeTorsion(a, xyz_[n], c, len, ang, tor + 120.0 * j));
      return out;
    }
    std::vector<Vec3> u;
    for (int n : nb) u.push_back(unit(sub(xyz_[n], c)));
    Vec3 sum{0, 0, 0};
    for (const Vec3& x : u) sum = add(sum, x);
    const Vec3 b = norm(sum) < 1e-6 ? perpendicular(u[0]) : unit(mul(sum, -1.0));
    if (k == 2 && hyb_[p] == 3) {
      const Vec3 nn = unit(cross(u[0], u[1]));
      const double half = 54.75 * M_PI / 180.0;
      out.push_back(add(c, mul(add(mul(b, std::cos(half)), mul(nn, std::sin(half))), len)));
      if (h >= 2) out.push_back(add(c, mul(add(mul(b, std::cos(half)), mul(nn, -std::sin(half))), len)));
      return out;
    }
    out.push_back(add(c, mul(b, len)));
    return out;
  }

  // Places the hydrogens of every heavy atom (rotors at their first state).
  void placeAll() {
    for (std::size_t p = 0; p < s_.atoms.size(); ++p) {
      if (isH(static_cast<int>(p))) continue;
      const std::vector<int> hs = hydrogens(static_cast<int>(p));
      if (hs.empty()) continue;
      const std::vector<Vec3> v = place(static_cast<int>(p), static_cast<int>(hs.size()), 180.0);
      for (std::size_t j = 0; j < hs.size() && j < v.size(); ++j) xyz_[hs[j]] = v[j];
    }
  }

  std::vector<Vec3>& xyz() { return xyz_; }
  const std::vector<std::vector<int>>& adj() const { return adj_; }
  void rebuildAdjacency() { adj_ = adjacency(s_); }
  int hyb(int i) const { return hyb_[i]; }

private:
  RawStructure& s_;
  std::vector<int> hyb_;
  std::vector<std::vector<int>> adj_;
  std::vector<Vec3> xyz_;
};

// Score units: one ideal hydrogen bond = -1. Overlaps up to 0.4 Å cost kSmallOverlap per Å;
// larger overlaps are steric clashes. Flips must gain more than kFlipBias to be applied.
constexpr double kSmallOverlap = 2.0;
constexpr double kFlipBias = 1.0;

struct Scorer {
  const RawStructure& s;
  const std::vector<std::vector<int>>& adj;
  std::vector<int> hCountFixed;  // hydrogens not in any movable set, per heavy atom
  std::vector<std::vector<std::pair<int, int>>> near3;  // heavy atoms within 3 bonds: (atom, bonds)

  Scorer(const RawStructure& st, const std::vector<std::vector<int>>& a) : s(st), adj(a) {
    const int n = static_cast<int>(s.atoms.size());
    near3.resize(n);
    for (int i = 0; i < n; ++i) {
      if (s.atoms[i].atomicNum == 1) continue;
      std::vector<std::pair<int, int>> front{{i, 0}};
      std::set<int> seen{i};
      for (std::size_t f = 0; f < front.size(); ++f) {
        const auto [x, d] = front[f];
        near3[i].push_back({x, d});
        if (d == 3) continue;
        for (int y : adj[x])
          if (s.atoms[y].atomicNum != 1 && seen.insert(y).second) front.push_back({y, d + 1});
      }
    }
  }
  int heavyBonds(int a, int b) const {
    for (const auto& [x, d] : near3[a])
      if (x == b) return d;
    return 99;
  }
  // Bond separation (<= 3 excluded from scoring).
  bool excluded(const AtomView& a, const AtomView& b) const {
    const int pa = a.parent >= 0 ? a.parent : a.idx, pb = b.parent >= 0 ? b.parent : b.idx;
    const int extra = (a.parent >= 0) + (b.parent >= 0);
    return heavyBonds(pa, pb) + extra <= 3;
  }
  double pair(const AtomView& a, const AtomView& b) const {
    const double d = dist(a.p, b.p);
    if (d > 4.5 || excluded(a, b)) return 0;
    auto hb = [&](const AtomView& h, const AtomView& acc) {
      if (d < 1.5) return 10.0 * (1.5 - d) * (1.5 - d);
      const double f = d <= 2.2 ? 1.0 : d >= 2.8 ? 0.0 : (2.8 - d) / 0.6;
      const Vec3 dh = sub(h.pp, h.p), ah = sub(acc.p, h.p);
      const double c = dot(dh, ah) / (norm(dh) * norm(ah));
      const double g = c < 0 ? -c : 0.0;
      return -f * g;
    };
    if (a.polarH && b.acceptor) return hb(a, b);
    if (b.polarH && a.acceptor) return hb(b, a);
    double e = 0;
    if (a.polarH && b.polarH && d < 2.2) e += 2.0 * (2.2 - d);
    if ((a.polarH && b.donorHeavy && !b.acceptor && d < 2.6) || (b.polarH && a.donorHeavy && !a.acceptor && d < 2.6)) e += 2.6 - d;
    if (a.acceptorOnly && b.acceptorOnly && d < 3.2) e += 3.2 - d;
    const double ov = a.radius + b.radius - d;
    if (ov > 0) e += kSmallOverlap * std::min(ov, 0.4);
    if (ov > 0.4) e += 10.0 * (ov - 0.4) * (ov - 0.4);
    return e;
  }
};


double radiusOf(int z, bool polarH) {
  switch (z) {
    case 1: return polarH ? 1.0 : 1.17;
    case 6: return 1.75;
    case 7: return 1.55;
    case 8: return 1.40;
    case 16: return 1.80;
    case 15: return 1.80;
    case 9: return 1.30;
    case 17: return 1.77;
    case 35: return 1.95;
    case 53: return 2.10;
    default: return isMetal(z) ? 1.0 : 1.7;
  }
}

}  // namespace

const std::vector<std::string>& residueStates(const std::string& resName) {
  static const std::map<std::string, std::vector<std::string>> m{
      {"ASP", {"ASP", "ASH"}}, {"GLU", {"GLU", "GLH"}}, {"HIS", {"HID", "HIE", "HIP"}}, {"CYS", {"CYS", "CYM"}},
      {"TYR", {"TYR", "TYM"}}, {"LYS", {"LYS", "LYN"}}, {"ARG", {"ARG", "ARN"}}};
  static const std::vector<std::string> none;
  auto it = m.find(resName);
  return it == m.end() ? none : it->second;
}

void protonate(RawStructure& s, const CcdLibrary* ccd, const ProtonationOptions& opt, ProtonationReport& rep,
               std::vector<int>* ccdCharge, std::vector<bool>* aromaticCandidate, std::vector<std::string>& errors,
               ProtonationPlan* planOut, const ProtonationPlan* planIn) {
  const bool ccdMode = ccd != nullptr;
  // 1. Remove the hydrogens of the input.
  {
    std::vector<int> keep;
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
      if (s.atoms[i].atomicNum == 1) ++rep.removedH;
      else keep.push_back(static_cast<int>(i));
    }
    reorder(s, keep, ccdCharge, aromaticCandidate);
  }
  int maxSerial = 0;
  for (const Atom& a : s.atoms) maxSerial = std::max(maxSerial, a.fileSerial);

  // 2. C-terminal OXT of amino-acid residues whose C has no link (chain ends and gaps).
  auto residueComp = [&](int r) -> const CcdComp* { return ccdMode && r >= 0 ? ccd->find(s.residues[r].info.name) : nullptr; };
  if (ccdMode) {
    const auto adj = adjacency(s);
    const std::size_t nAtoms = s.atoms.size();
    for (std::size_t r = 0; r < s.residues.size(); ++r) {
      const CcdComp* comp = residueComp(static_cast<int>(r));
      if (!comp || comp->find("OXT") < 0) continue;
      const int c = findInResidue(s, static_cast<int>(r), "C"), o = findInResidue(s, static_cast<int>(r), "O"),
                ca = findInResidue(s, static_cast<int>(r), "CA");
      if (c < 0 || o < 0 || ca < 0 || findInResidue(s, static_cast<int>(r), "OXT") >= 0) continue;
      if (static_cast<std::size_t>(c) >= nAtoms || adj[c].size() != 2) continue;
      const bool wanted = planIn ? std::find(planIn->addedOxt.begin(), planIn->addedOxt.end(), resLabel(s, static_cast<int>(r))) != planIn->addedOxt.end()
                                 : true;
      if (!wanted) continue;
      Atom a;
      a.name = "OXT";
      a.atomicNum = 8;
      a.residue = static_cast<int>(r);
      a.hetatm = s.atoms[c].hetatm;
      const Vec3 u = add(unit(sub(s.atoms[ca].xyz, s.atoms[c].xyz)), unit(sub(s.atoms[o].xyz, s.atoms[c].xyz)));
      a.xyz = add(s.atoms[c].xyz, mul(unit(mul(u, -1.0)), 1.25));
      a.fileSerial = ++maxSerial;
      s.atoms.push_back(a);
      const int oi = static_cast<int>(s.atoms.size()) - 1;
      s.residues[r].atoms.push_back(oi);
      s.bonds.push_back({c, oi, BondOrder::Single, false});
      if (ccdCharge) ccdCharge->push_back(0);
      if (aromaticCandidate) aromaticCandidate->push_back(false);
      rep.addedOxt.push_back(resLabel(s, static_cast<int>(r)));
    }
  }
  const int n = static_cast<int>(s.atoms.size());
  if (n == 0) return;
  const auto adj0 = adjacency(s);
  std::vector<bool> metalBonded(n, false);
  for (const Bond& b : s.bonds)
    if (isMetalBondOf(s, b)) metalBonded[b.a] = metalBonded[b.b] = true;
  for (int i = 0; i < n; ++i)
    if (isMetal(s.atoms[i].atomicNum)) metalBonded[i] = false;

  // Residue completeness (heavy atoms of the CCD definition that are not leaving atoms).
  if (ccdMode && !planIn) {
    for (std::size_t r = 0; r < s.residues.size(); ++r) {
      const CcdComp* comp = residueComp(static_cast<int>(r));
      if (!comp || !kStdAmino.contains(comp->id)) continue;
      for (const CcdComp::Atom& ca : comp->atoms) {
        if (ca.atomicNum == 1 || ca.name == "OXT") continue;
        if (findInResidue(s, static_cast<int>(r), ca.name) < 0) {
          rep.incomplete.push_back(resLabel(s, static_cast<int>(r)));
          break;
        }
      }
    }
  }

  std::vector<int> charge(n), hCount(n, 0);
  for (int i = 0; i < n; ++i) charge[i] = ccdMode ? (*ccdCharge)[i] : s.atoms[i].formalCharge;
  std::vector<int> hisTautomer;  // HIS residues whose neutral tautomer is optimized
  std::vector<int> carboxylH;    // heavy O atoms carrying the H of a protonated carboxyl (CCD residues)

  if (planIn) {
    if (planIn->formalCharge.size() != static_cast<std::size_t>(n) || planIn->bondOrder.size() != s.bonds.size())
      fail("protonation: records differ in heavy atoms or bonds");
    charge = planIn->formalCharge;
    for (std::size_t b = 0; b < s.bonds.size(); ++b) s.bonds[b].order = planIn->bondOrder[b];
    for (const auto& [p, name] : planIn->hydrogens) ++hCount[p];
  } else {
    // 3. Overrides.
    std::vector<int> atomOverride(n, std::numeric_limits<int>::min());
    std::map<int, std::string> residueOverride;
    for (const ProtonationOverride& o : opt.overrides) {
      const std::string& st = o.site;
      int atom = -1, res = -1;
      if (!st.empty() && st[0] == '#') {
        int serial = 0;
        try {
          serial = std::stoi(st.substr(1));
        } catch (const std::exception&) {
        }
        for (int i = 0; i < n; ++i)
          if (s.atoms[i].fileSerial == serial && serial > 0) atom = i;
      } else {
        const auto c1 = st.find(':'), c2 = st.find(':', c1 == std::string::npos ? 0 : c1 + 1);
        if (c1 != std::string::npos) {
          const std::string chain = st.substr(0, c1);
          const std::string seqPart = st.substr(c1 + 1, c2 == std::string::npos ? std::string::npos : c2 - c1 - 1);
          std::size_t used = 0;
          int seq = 0;
          try {
            seq = std::stoi(seqPart, &used);
          } catch (const std::exception&) {
            used = 0;
          }
          const std::string icode = used > 0 ? seqPart.substr(used) : "";
          for (std::size_t r = 0; r < s.residues.size() && used > 0; ++r) {
            const ResidueInfo& ri = s.residues[r].info;
            if (ri.chain == chain && ri.seq == seq && ri.iCode == icode) res = static_cast<int>(r);
          }
          if (res >= 0 && c2 != std::string::npos) {
            atom = findInResidue(s, res, st.substr(c2 + 1));
            if (atom < 0) res = -2;
          }
        }
      }
      if (o.state.empty()) {
        if (atom < 0) fail("{}: [[protonation]] site \"{}\": atom not found in {}", o.origin, st, o.file);
        if (!o.charge) fail("{}: [[protonation]] site \"{}\": charge is required for an atom entry", o.origin, st);
        atomOverride[atom] = *o.charge;
      } else {
        if (res < 0 || atom >= 0) fail("{}: [[protonation]] site \"{}\": residue not found in {}", o.origin, st, o.file);
        const auto& ok = residueStates(s.residues[res].info.name);
        if (std::find(ok.begin(), ok.end(), o.state) == ok.end())
          fail("{}: [[protonation]] state \"{}\" is not valid for {} (valid: {})", o.origin, o.state, resLabel(s, res),
               [&] {
                 std::string v;
                 for (const auto& x : ok) v += (v.empty() ? "" : ", ") + x;
                 return v.empty() ? std::string("none") : v;
               }());
        if (!ccdMode) fail("{}: [[protonation]] residue states need mmCIF input", o.origin);
        residueOverride[res] = o.state;
      }
    }

    // 4. Heavy-atom molecule with the CCD / input charges; implicit hydrogens.
    std::vector<bool> fixed(n, false);
    for (int i = 0; i < n; ++i) {
      const int r = s.atoms[i].residue;
      if (ccdMode && r >= 0 && kStdAmino.contains(s.residues[r].info.name)) fixed[i] = true;
      if (isMetal(s.atoms[i].atomicNum)) fixed[i] = true;
      if (atomOverride[i] != std::numeric_limits<int>::min()) {
        fixed[i] = true;
        charge[i] = atomOverride[i];
      }
    }
    auto mol = heavyMol(s, charge);
    try {
      RDKit::MolOps::sanitizeMol(*mol);
    } catch (const std::exception& e) {
      fail("{}: chemistry check before protonation failed: {}", s.sourceFile, e.what());
    }
    std::vector<int> baseH(n);
    for (int i = 0; i < n; ++i) baseH[i] = static_cast<int>(mol->getAtomWithIdx(static_cast<unsigned>(i))->getTotalNumHs());

    // 5. Dimorphite-DL for all atoms outside standard amino-acid residues.
    dimorphite::Result dm = dimorphite::protonate(*mol, opt.pH, fixed, metalBonded);
    for (int i = 0; i < n; ++i) {
      const RDKit::Atom* a = dm.mol->getAtomWithIdx(static_cast<unsigned>(i));
      hCount[i] = fixed[i] ? baseH[i] : static_cast<int>(a->getTotalNumHs());
      if (!fixed[i]) charge[i] = a->getFormalCharge();
    }
    for (const auto& st : dm.sites)
      rep.sites.push_back({describeAtom(s, st.atom), st.forced ? "metal" : "Dimorphite-DL " + st.name, st.pka,
                           st.protonated ? "protonated" : "deprotonated"});
    for (std::size_t b = 0; b < s.bonds.size(); ++b) {
      const Bond& bd = s.bonds[b];
      if (isMetalBondOf(s, bd) || (aromaticCandidate && (*aromaticCandidate)[b])) continue;
      if (fixed[bd.a] && fixed[bd.b]) continue;
      const RDKit::Bond* rb = dm.mol->getBondBetweenAtoms(static_cast<unsigned>(bd.a), static_cast<unsigned>(bd.b));
      if (!rb) continue;
      const auto t = rb->getBondType();
      s.bonds[b].order = t == RDKit::Bond::DOUBLE ? BondOrder::Double : t == RDKit::Bond::TRIPLE ? BondOrder::Triple : BondOrder::Single;
    }
    for (int i = 0; i < n; ++i)
      if (atomOverride[i] != std::numeric_limits<int>::min())
        rep.sites.push_back({describeAtom(s, i), "override", std::numeric_limits<double>::quiet_NaN(),
                             std::format("charge {:+}", atomOverride[i])});

    // 6. Standard amino acids: PROPKA 3.
    bool anyAmino = false;
    for (std::size_t r = 0; r < s.residues.size() && ccdMode; ++r)
      if (kStdAmino.contains(s.residues[r].info.name)) anyAmino = true;
    if (anyAmino) {
      std::vector<propka::InputAtom> in;
      std::vector<int> inputOf;
      for (int i = 0; i < n; ++i) {
        const Atom& a = s.atoms[i];
        const ResidueInfo& ri = s.residues[a.residue].info;
        in.push_back({a.name, ri.name, ri.chain, ri.seq, ri.iCode, elementSymbol(a.atomicNum), a.hetatm, a.xyz});
        inputOf.push_back(i);
      }
      std::map<std::pair<int, std::string>, double> pka;  // (atom, residue type) -> pKa
      for (const auto& g : propka::calculate(in)) pka[{inputOf[g.atom], g.residueType}] = g.pka;
      auto pkaOf = [&](int atom, const std::string& type, double model) {
        auto it = pka.find({atom, type});
        return it == pka.end() ? std::pair<double, bool>{model, false} : std::pair<double, bool>{it->second, true};
      };
      for (std::size_t ri = 0; ri < s.residues.size(); ++ri) {
        const int r = static_cast<int>(ri);
        const std::string& rn = s.residues[r].info.name;
        if (!kStdAmino.contains(rn)) continue;
        auto at = [&](const char* nm) { return findInResidue(s, r, nm); };
        auto ovr = residueOverride.find(r);
        auto record = [&](int atom, double pk, bool fromPropka, const std::string& state, bool overridden, bool metal) {
          rep.sites.push_back({describeAtom(s, atom), overridden ? "override" : metal ? "metal" : fromPropka ? "PROPKA" : "PROPKA model pKa",
                               overridden ? std::numeric_limits<double>::quiet_NaN() : pk, state});
        };
        // Carboxyl: ASP, GLU side chains and C-terminal OXT.
        auto carboxyl = [&](int carbon, int oA, int oH, const std::string& type, const char* deprot, const char* prot, double model) {
          if (carbon < 0 || oA < 0 || oH < 0) return;
          const auto [pk, fp] = pkaOf(type == "C-" ? oH : carbon, type, model);
          bool protonated = pk >= opt.pH;
          const bool metal = metalBonded[oA] || metalBonded[oH];
          bool overridden = false;
          if (type != "C-" && ovr != residueOverride.end()) {
            protonated = ovr->second == prot;
            overridden = true;
          }
          if (metal && !overridden) protonated = false;
          if (atomOverride[oA] != std::numeric_limits<int>::min() || atomOverride[oH] != std::numeric_limits<int>::min()) return;
          hCount[oA] = 0;
          hCount[oH] = protonated ? 1 : 0;
          if (protonated) carboxylH.push_back(oH);
          record(oH, pk, fp, protonated ? prot : deprot, overridden, metal);
        };
        if (rn == "ASP") carboxyl(at("CG"), at("OD1"), at("OD2"), "ASP", "ASP", "ASH", 3.8);
        if (rn == "GLU") carboxyl(at("CD"), at("OE1"), at("OE2"), "GLU", "GLU", "GLH", 4.5);
        if (at("OXT") >= 0) carboxyl(at("C"), at("O"), at("OXT"), "C-", "COO-", "COOH", 3.2);
        // N terminus: N without a link to another residue.
        if (const int nn = at("N"); nn >= 0 && atomOverride[nn] == std::numeric_limits<int>::min()) {
          bool linked = false;
          int heavy = 0;
          for (int x : adj0[nn]) {
            if (s.atoms[x].residue != r && !isMetal(s.atoms[x].atomicNum)) linked = true;
            if (!isMetal(s.atoms[x].atomicNum)) ++heavy;
          }
          if (!linked) {
            const auto [pk, fp] = pkaOf(nn, "N+", 8.0);
            const bool metal = metalBonded[nn];
            const bool protonated = pk >= opt.pH && !metal;
            hCount[nn] = std::max(0, (protonated ? 4 : 3) - heavy);
            record(nn, pk, fp, protonated ? "NH3+" : "NH2", false, metal);
          }
        }
        auto simple = [&](int atom, int group, const std::string& type, int hProt, const std::string& deprot,
                          const std::string& prot, double model) {
          if (atom < 0 || atomOverride[atom] != std::numeric_limits<int>::min()) return;
          const auto [pk, fp] = pkaOf(group, type, model);
          bool protonated = pk >= opt.pH;
          bool overridden = false;
          if (ovr != residueOverride.end()) {
            protonated = ovr->second == prot;
            overridden = true;
          }
          const bool metal = metalBonded[atom];
          if (metal && !overridden) protonated = false;
          hCount[atom] = protonated ? hProt : hProt - 1;
          record(atom, pk, fp, protonated ? prot : deprot, overridden, metal);
        };
        if (rn == "CYS") {
          const int sg = at("SG");
          bool bridged = false;
          if (sg >= 0)
            for (int x : adj0[sg])
              if (s.atoms[x].atomicNum == 16) bridged = true;
          if (!bridged) simple(sg, sg, "CYS", 1, "CYM", "CYS", 9.0);
        }
        if (rn == "TYR") simple(at("OH"), at("OH"), "TYR", 1, "TYM", "TYR", 10.0);
        if (rn == "LYS") simple(at("NZ"), at("NZ"), "LYS", 3, "LYN", "LYS", 10.5);
        if (rn == "ARG") {
          const int nh2 = at("NH2");
          if (nh2 >= 0 && at("CZ") >= 0) simple(nh2, at("CZ"), "ARG", 2, "ARN", "ARG", 12.5);
        }
        if (rn == "HIS") {
          const int cg = at("CG"), nd1 = at("ND1"), ne2 = at("NE2");
          if (cg >= 0 && nd1 >= 0 && ne2 >= 0 && atomOverride[nd1] == std::numeric_limits<int>::min() &&
              atomOverride[ne2] == std::numeric_limits<int>::min()) {
            const auto [pk, fp] = pkaOf(cg, "HIS", 6.5);
            std::string state = pk >= opt.pH ? "HIP" : "HIE";
            bool overridden = false, metal = false;
            if (ovr != residueOverride.end()) {
              state = ovr->second;
              overridden = true;
            } else if (metalBonded[nd1] || metalBonded[ne2]) {
              metal = true;
              state = metalBonded[nd1] && metalBonded[ne2] ? "HIM" : metalBonded[nd1] ? "HIE" : "HID";
            } else if (state == "HIE") {
              hisTautomer.push_back(r);
            }
            hCount[nd1] = (state == "HID" || state == "HIP") ? 1 : 0;
            hCount[ne2] = (state == "HIE" || state == "HIP") ? 1 : 0;
            record(cg, pk, fp, state == "HIM" ? "HIS (both N metal-bonded)" : (state == "HIE" && !overridden && !metal) ? "HID/HIE" : state,
                   overridden, metal);
          }
        }
      }
    }
  }

  // 7. Create the hydrogens (positions follow), ordered after the heavy atoms of their residue.
  const int nHeavy = n;
  std::vector<int> parentOfH;
  for (int p = 0; p < nHeavy; ++p)
    for (int k = 0; k < hCount[p]; ++k) {
      Atom h;
      h.atomicNum = 1;
      h.residue = s.atoms[p].residue;
      h.hetatm = s.atoms[p].hetatm;
      h.xyz = s.atoms[p].xyz;
      h.fileSerial = ++maxSerial;
      s.atoms.push_back(h);
      const int hi = static_cast<int>(s.atoms.size()) - 1;
      if (h.residue >= 0) s.residues[h.residue].atoms.push_back(hi);
      s.bonds.push_back({p, hi, BondOrder::Single, false});
      if (ccdCharge) ccdCharge->push_back(0);
      if (aromaticCandidate) aromaticCandidate->push_back(false);
      parentOfH.push_back(p);
    }
  if (planIn) {
    std::vector<int> next(nHeavy, 0);
    std::vector<std::vector<std::string>> names(nHeavy);
    for (const auto& [p, name] : planIn->hydrogens) names[p].push_back(name);
    for (std::size_t k = 0; k < parentOfH.size(); ++k) {
      const int p = parentOfH[k];
      s.atoms[nHeavy + static_cast<int>(k)].name = names[p][next[p]++];
    }
  }
  rep.addedH = static_cast<int>(parentOfH.size());
  for (int i = 0; i < nHeavy; ++i) s.atoms[i].formalCharge = charge[i];

  // 8. Charges and Kekulé bonds with the hydrogens present; hybridization.
  auto assign = [&]() {
    if (ccdMode) {
      std::vector<std::string> err;
      assignChargesAndKekule(s, *aromaticCandidate, *ccdCharge, err);
      return err;
    }
    return std::vector<std::string>{};
  };
  if (auto err = assign(); !err.empty()) {
    errors.insert(errors.end(), err.begin(), err.end());
    return;
  }
  std::vector<int> hyb(s.atoms.size(), 3);
  {
    Molecule tmp;
    tmp.atoms = s.atoms;
    tmp.bonds = s.bonds;
    auto full = std::make_unique<RDKit::RWMol>();
    for (const Atom& a : tmp.atoms) {
      RDKit::Atom x(a.atomicNum);
      x.setFormalCharge(a.formalCharge);
      x.setNoImplicit(true);
      full->addAtom(&x, true, false);
    }
    for (const Bond& b : tmp.bonds) {
      const auto t = b.order == BondOrder::Metal || isMetalBondOf(s, b) ? RDKit::Bond::ZERO
                     : b.order == BondOrder::Double ? RDKit::Bond::DOUBLE
                     : b.order == BondOrder::Triple ? RDKit::Bond::TRIPLE : RDKit::Bond::SINGLE;
      full->addBond(static_cast<unsigned>(b.a), static_cast<unsigned>(b.b), t);
    }
    try {
      RDKit::MolOps::sanitizeMol(*full);
    } catch (const std::exception& e) {
      errors.push_back(std::format("chemistry check after protonation failed: {}", e.what()));
      return;
    }
    for (std::size_t i = 0; i < tmp.atoms.size(); ++i) {
      const auto h = full->getAtomWithIdx(static_cast<unsigned>(i))->getHybridization();
      hyb[i] = h == RDKit::Atom::SP ? 1 : h == RDKit::Atom::SP2 ? 2 : 3;
    }
  }

  // 9. Placement and optimization of rotatable hydrogens and flips.
  Placer pl(s, hyb);
  pl.placeAll();
  std::vector<OptGroup> groups;
  std::vector<bool> inGroup(s.atoms.size(), false);
  auto stateFromXyz = [&](const std::vector<int>& moving) {
    GroupState st;
    for (int a : moving) st.pos.push_back({a, pl.xyz()[a]});
    return st;
  };
  // Flip groups of CCD residues.
  if (ccdMode) {
    for (std::size_t ri = 0; ri < s.residues.size(); ++ri) {
      const int r = static_cast<int>(ri);
      const std::string& rn = s.residues[r].info.name;
      auto at = [&](const char* nm) { return findInResidue(s, r, nm); };
      auto swapStates = [&](const std::vector<std::pair<int, int>>& pairs, const std::vector<int>& parents) {
        std::vector<int> moving;
        for (const auto& [a, b] : pairs) {
          moving.push_back(a);
          moving.push_back(b);
        }
        for (int p : parents)
          for (int h : pl.hydrogens(p)) moving.push_back(h);
        OptGroup g;
        g.label = resLabel(s, r);
        g.moving = moving;
        for (int f = 0; f < 2; ++f) {
          std::vector<Vec3> saved = pl.xyz();
          if (f == 1)
            for (const auto& [a, b] : pairs) std::swap(pl.xyz()[a], pl.xyz()[b]);
          for (int p : parents) {
            const auto hs = pl.hydrogens(p);
            const auto v = pl.place(p, static_cast<int>(hs.size()), 180.0);
            for (std::size_t j = 0; j < hs.size(); ++j) pl.xyz()[hs[j]] = v[j];
          }
          GroupState st = stateFromXyz(moving);
          st.flipped = f == 1;
          st.bias = f == 1 ? kFlipBias : 0.0;
          g.states.push_back(std::move(st));
          pl.xyz() = saved;
        }
        return g;
      };
      if (rn == "ASN" || rn == "GLN") {
        const int o = at(rn == "ASN" ? "OD1" : "OE1"), nn = at(rn == "ASN" ? "ND2" : "NE2");
        if (o < 0 || nn < 0 || metalBonded[o] || metalBonded[nn]) continue;
        OptGroup g = swapStates({{o, nn}}, {nn});
        for (int a : g.moving) inGroup[a] = true;
        groups.push_back(std::move(g));
      }
      if (rn == "HIS") {
        const int nd1 = at("ND1"), cd2 = at("CD2"), ce1 = at("CE1"), ne2 = at("NE2");
        if (nd1 < 0 || cd2 < 0 || ce1 < 0 || ne2 < 0 || metalBonded[nd1] || metalBonded[ne2]) continue;
        const bool tautomer = std::find(hisTautomer.begin(), hisTautomer.end(), r) != hisTautomer.end();
        if (!tautomer) {
          OptGroup g = swapStates({{nd1, cd2}, {ce1, ne2}}, {nd1, cd2, ce1, ne2});
          for (int a : g.moving) inGroup[a] = true;
          groups.push_back(std::move(g));
          continue;
        }
        // Neutral His: the N-H hydrogen on ND1 (HID) or NE2 (HIE), with and without the flip.
        int hN = -1;
        for (int h : pl.hydrogens(nd1)) hN = h;
        for (int h : pl.hydrogens(ne2)) hN = h;
        if (hN < 0) continue;
        OptGroup g;
        g.label = resLabel(s, r);
        g.moving = {nd1, cd2, ce1, ne2, hN};
        for (int p : {cd2, ce1})
          for (int h : pl.hydrogens(p)) g.moving.push_back(h);
        for (int f = 0; f < 2; ++f)
          for (int t = 0; t < 2; ++t) {
            const int np = t == 0 ? ne2 : nd1;  // HIE first
            std::vector<Vec3> saved = pl.xyz();
            if (f == 1) {
              std::swap(pl.xyz()[nd1], pl.xyz()[cd2]);
              std::swap(pl.xyz()[ce1], pl.xyz()[ne2]);
            }
            GroupState st;
            for (int p : {cd2, ce1}) {
              const auto hs = pl.hydrogens(p);
              const auto v = pl.place(p, static_cast<int>(hs.size()), 180.0);
              for (std::size_t j = 0; j < hs.size(); ++j) pl.xyz()[hs[j]] = v[j];
            }
            // Bisector placement on the chosen N (two ring neighbours).
            {
              const auto nbs = pl.heavyNbrs(np);
              Vec3 sum{0, 0, 0};
              for (int x : nbs) sum = add(sum, unit(sub(pl.xyz()[x], pl.xyz()[np])));
              pl.xyz()[hN] = add(pl.xyz()[np], mul(unit(mul(sum, -1.0)), bondLength(7)));
            }
            for (int a : g.moving) st.pos.push_back({a, pl.xyz()[a]});
            st.parent.push_back({hN, np});
            st.flipped = f == 1;
            st.bias = f == 1 ? kFlipBias : 0.0;
            g.states.push_back(std::move(st));
            pl.xyz() = saved;
          }
        for (int a : g.moving) inGroup[a] = true;
        groups.push_back(std::move(g));
      }
    }
    // Protonated carboxyl groups: which oxygen carries the H, and syn / anti.
    for (int oh : carboxylH) {
      int c = -1;
      for (int x : pl.adj()[oh])
        if (s.atoms[x].atomicNum == 6) c = x;
      if (c < 0) continue;
      int oo = -1;
      for (int x : pl.adj()[c])
        if (x != oh && s.atoms[x].atomicNum == 8 && pl.heavyNbrs(x).size() == 1 && pl.hydrogens(x).empty()) oo = x;
      const auto hs = pl.hydrogens(oh);
      if (oo < 0 || hs.size() != 1 || inGroup[oh] || inGroup[oo]) continue;
      OptGroup g;
      g.label = describeAtom(s, oh);
      g.moving = {oh, oo, hs[0]};
      for (int f = 0; f < 2; ++f)
        for (int t = 0; t < 2; ++t) {
          std::vector<Vec3> saved = pl.xyz();
          if (f == 1) std::swap(pl.xyz()[oh], pl.xyz()[oo]);
          pl.xyz()[hs[0]] = placeTorsion(pl.xyz()[oo], pl.xyz()[c], pl.xyz()[oh], bondLength(8), 109.47, t == 0 ? 0.0 : 180.0);
          GroupState st = stateFromXyz(g.moving);
          st.flipped = f == 1;
          st.bias = (f == 1 ? kFlipBias : 0.0) + (t == 1 ? 0.05 : 0.0);
          g.states.push_back(std::move(st));
          pl.xyz() = saved;
        }
      for (int a : g.moving) inGroup[a] = true;
      groups.push_back(std::move(g));
    }
  }
  // Rotors: OH, SH, NH3+, sp3 NH2, waters.
  for (int p = 0; p < nHeavy; ++p) {
    const int z = s.atoms[p].atomicNum;
    if ((z != 7 && z != 8 && z != 16) || inGroup[p]) continue;
    const auto hs = pl.hydrogens(p);
    if (hs.empty()) continue;
    bool any = false;
    for (int h : hs) any = any || inGroup[h];
    if (any) continue;
    const auto nb = pl.heavyNbrs(p);
    OptGroup g;
    g.label = describeAtom(s, p);
    g.moving = hs;
    g.rotor = true;
    if (nb.empty()) {
      if (hs.size() > 2) continue;
      for (int o = 0; o < 36; ++o) {
        const auto v = pl.place(p, static_cast<int>(hs.size()), 0.0, o);
        GroupState st;
        for (std::size_t j = 0; j < hs.size(); ++j) st.pos.push_back({hs[j], v[j]});
        g.states.push_back(std::move(st));
      }
    } else if (nb.size() == 1) {
      bool single = true;
      for (const Bond& b : s.bonds)
        if ((b.a == p && b.b == nb[0]) || (b.b == p && b.a == nb[0])) single = b.order == BondOrder::Single || b.order == BondOrder::Metal;
      if (!single || pl.hyb(p) == 1) continue;
      if (pl.hyb(p) == 2) {
        if (hs.size() != 1) continue;
        for (double t : {180.0, 0.0}) {
          const auto v = pl.place(p, 1, t);
          GroupState st;
          st.pos.push_back({hs[0], v[0]});
          g.states.push_back(std::move(st));
        }
      } else {
        // 30-degree steps (NH3+: one 120-degree period), refined after the search.
        const int nStates = hs.size() == 3 ? 4 : 12;
        g.rotorParent = p;
        for (int k = 0; k < nStates; ++k) {
          const auto v = pl.place(p, static_cast<int>(hs.size()), 180.0 + 30.0 * k);
          GroupState st;
          st.tor = 180.0 + 30.0 * k;
          for (std::size_t j = 0; j < hs.size(); ++j) st.pos.push_back({hs[j], v[j]});
          g.states.push_back(std::move(st));
        }
      }
    } else {
      continue;
    }
    for (int a : g.moving) inGroup[a] = true;
    groups.push_back(std::move(g));
  }

  // Scoring.
  const int nAll = static_cast<int>(s.atoms.size());
  Scorer sc(s, pl.adj());
  std::vector<int> parentOf(nAll, -1);
  for (const Bond& b : s.bonds) {
    if (s.atoms[b.a].atomicNum == 1) parentOf[b.a] = b.b;
    if (s.atoms[b.b].atomicNum == 1) parentOf[b.b] = b.a;
  }
  std::vector<int> hOn(nAll, 0);  // hydrogens per heavy atom (current assignment)
  for (int i = 0; i < nAll; ++i)
    if (parentOf[i] >= 0) ++hOn[parentOf[i]];
  auto view = [&](int idx, const Vec3& p, int parent, int hCnt) {
    const int z = s.atoms[idx].atomicNum;
    AtomView v{idx, p, parent, false, false, false, 0};
    if (z == 1) {
      const int pz = parent >= 0 ? s.atoms[parent].atomicNum : 6;
      v.polarH = pz == 7 || pz == 8 || pz == 16;
    } else if (z == 8 || z == 16) {
      v.acceptor = s.atoms[idx].formalCharge <= 0;
    } else if (z == 7) {
      v.acceptor = hCnt == 0 && s.atoms[idx].formalCharge <= 0 && pl.heavyNbrs(idx).size() < 3;
    }
    v.acceptorOnly = v.acceptor && hCnt == 0;
    v.donorHeavy = (z == 7 || z == 8 || z == 16) && hCnt > 0;
    if (parent >= 0) v.pp = pl.xyz()[parent];
    v.radius = radiusOf(z, v.polarH);
    return v;
  };
  // Hydrogen counts with movable N-H hydrogens removed (state dependent).
  std::vector<int> hStatic = hOn;
  for (const OptGroup& g : groups)
    for (const GroupState& st : g.states)
      for (const auto& [h, p] : st.parent) {
        (void)p;
        if (parentOf[h] >= 0) {
          hStatic[parentOf[h]] = hOn[parentOf[h]] - 1;
        }
      }
  std::vector<AtomView> fixedViews;
  for (int i = 0; i < nAll; ++i)
    if (!inGroup[i]) fixedViews.push_back(view(i, pl.xyz()[i], parentOf[i], hStatic[i]));
  // Cell list of fixed atoms.
  const double cell = 4.5;
  std::unordered_map<long long, std::vector<int>> cells;
  auto key = [&](const Vec3& p, int dx, int dy, int dz) {
    const long long x = static_cast<long long>(std::floor(p[0] / cell)) + dx, y = static_cast<long long>(std::floor(p[1] / cell)) + dy,
                    z = static_cast<long long>(std::floor(p[2] / cell)) + dz;
    return (x * 73856093LL) ^ (y * 19349663LL) ^ (z * 83492791LL);
  };
  for (std::size_t k = 0; k < fixedViews.size(); ++k) cells[key(fixedViews[k].p, 0, 0, 0)].push_back(static_cast<int>(k));
  auto stateViews = [&](const GroupState& st) {
    std::vector<AtomView> v;
    std::map<int, int> par;
    for (const auto& [h, p] : st.parent) par[h] = p;
    std::map<int, int> extraH;
    for (const auto& [h, p] : st.parent) ++extraH[p];
    std::map<int, Vec3> inState;
    for (const auto& [a, p] : st.pos) inState[a] = p;
    for (const auto& [a, p] : st.pos) {
      const int parent = par.contains(a) ? par[a] : parentOf[a];
      const int hc = hStatic[a] + (extraH.contains(a) ? extraH[a] : 0);
      AtomView x = view(a, p, parent, hc);
      if (parent >= 0 && inState.contains(parent)) x.pp = inState[parent];
      v.push_back(x);
    }
    return v;
  };
  const int G = static_cast<int>(groups.size());
  std::vector<std::vector<std::vector<AtomView>>> views(G);
  std::vector<std::vector<double>> self(G);
  for (int g = 0; g < G; ++g) {
    for (const GroupState& st : groups[g].states) {
      auto v = stateViews(st);
      double e = st.bias;
      for (const AtomView& a : v) {
        std::set<int> done;
        for (int dx = -1; dx <= 1; ++dx)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
              auto it = cells.find(key(a.p, dx, dy, dz));
              if (it == cells.end()) continue;
              for (int k : it->second)
                if (done.insert(k).second) e += sc.pair(a, fixedViews[k]);
            }
      }
      self[g].push_back(e);
      views[g].push_back(std::move(v));
    }
  }
  // Interacting group pairs.
  std::vector<std::vector<int>> nbrG(G);
  std::map<std::pair<int, int>, std::vector<std::vector<double>>> pairE;
  {
    std::vector<Vec3> lo(G, Vec3{1e30, 1e30, 1e30}), hi(G, Vec3{-1e30, -1e30, -1e30});
    for (int g = 0; g < G; ++g)
      for (const auto& v : views[g])
        for (const AtomView& a : v)
          for (int d = 0; d < 3; ++d) {
            lo[g][d] = std::min(lo[g][d], a.p[d]);
            hi[g][d] = std::max(hi[g][d], a.p[d]);
          }
    for (int a = 0; a < G; ++a)
      for (int b = a + 1; b < G; ++b) {
        bool close = true;
        for (int d = 0; d < 3; ++d)
          if (lo[a][d] > hi[b][d] + 4.5 || lo[b][d] > hi[a][d] + 4.5) close = false;
        if (!close) continue;
        std::vector<std::vector<double>> t(groups[a].states.size(), std::vector<double>(groups[b].states.size(), 0));
        bool any = false;
        for (std::size_t i = 0; i < t.size(); ++i)
          for (std::size_t j = 0; j < t[i].size(); ++j) {
            double e = 0;
            for (const AtomView& x : views[a][i])
              for (const AtomView& y : views[b][j]) e += sc.pair(x, y);
            t[i][j] = e;
            any = any || e != 0;
          }
        if (!any) continue;
        nbrG[a].push_back(b);
        nbrG[b].push_back(a);
        pairE[{a, b}] = std::move(t);
      }
  }
  auto pairVal = [&](int a, int sa, int b, int sb) {
    if (a < b) return pairE.at({a, b})[sa][sb];
    return pairE.at({b, a})[sb][sa];
  };
  // Components.
  std::vector<int> comp(G, -1);
  int nComp = 0;
  for (int g = 0; g < G; ++g) {
    if (comp[g] >= 0) continue;
    std::vector<int> st{g};
    comp[g] = nComp;
    while (!st.empty()) {
      const int x = st.back();
      st.pop_back();
      for (int y : nbrG[x])
        if (comp[y] < 0) {
          comp[y] = nComp;
          st.push_back(y);
        }
    }
    ++nComp;
  }
  rep.components = nComp;
  constexpr double kExhaustiveLimit = 1e6;
  for (int c = 0; c < nComp; ++c) {
    std::vector<int> mem;
    for (int g = 0; g < G; ++g)
      if (comp[g] == c) mem.push_back(g);
    std::sort(mem.begin(), mem.end(), [&](int a, int b) { return nbrG[a].size() != nbrG[b].size() ? nbrG[a].size() > nbrG[b].size() : a < b; });
    // Dead-end elimination (Goldstein): drop states that some other state of the group always beats.
    std::vector<std::vector<int>> alive(G);
    for (int g : mem)
      for (std::size_t si = 0; si < groups[g].states.size(); ++si) alive[g].push_back(static_cast<int>(si));
    for (bool changed = true; changed;) {
      changed = false;
      for (int g : mem) {
        if (alive[g].size() <= 1) continue;
        std::vector<int> keep;
        for (int r : alive[g]) {
          bool dead = false;
          for (int u : alive[g]) {
            if (u == r) continue;
            double d = self[g][r] - self[g][u];
            for (int y : nbrG[g]) {
              double m = std::numeric_limits<double>::infinity();
              for (int sy : alive[y]) m = std::min(m, pairVal(g, r, y, sy) - pairVal(g, u, y, sy));
              d += m;
            }
            if (d > 1e-9) {
              dead = true;
              break;
            }
          }
          if (!dead) keep.push_back(r);
        }
        if (keep.size() < alive[g].size()) {
          alive[g] = std::move(keep);
          changed = true;
        }
      }
    }
    double product = 1;
    for (int g : mem) product *= static_cast<double>(alive[g].size());
    auto bestSelf = [&](int g) {
      int b = alive[g][0];
      for (int si : alive[g])
        if (self[g][si] < self[g][b]) b = si;
      return b;
    };
    if (product <= kExhaustiveLimit) {
      std::vector<int> cur(G, -1), best(G, -1);
      double bestE = std::numeric_limits<double>::infinity();
      std::vector<double> minSelf(mem.size());
      for (std::size_t k = 0; k < mem.size(); ++k) {
        minSelf[k] = std::numeric_limits<double>::infinity();
        for (int si : alive[mem[k]]) minSelf[k] = std::min(minSelf[k], self[mem[k]][si]);
      }
      std::vector<double> suffix(mem.size() + 1, 0);
      for (int k = static_cast<int>(mem.size()) - 1; k >= 0; --k) suffix[k] = suffix[k + 1] + minSelf[k];
      // Lower bound of pair energies: negative parts only.
      double pairLow = 0;
      for (std::size_t k = 0; k < mem.size(); ++k)
        for (int y : nbrG[mem[k]])
          if (y > mem[k]) {
            double m = 0;
            for (int si : alive[mem[k]])
              for (int sy : alive[y]) m = std::min(m, pairVal(mem[k], si, y, sy));
            pairLow += m;
          }
      std::function<void(std::size_t, double)> dfs = [&](std::size_t k, double e) {
        if (k == mem.size()) {
          if (e < bestE - 1e-12) {
            bestE = e;
            best = cur;
          }
          return;
        }
        if (e + suffix[k] + pairLow >= bestE - 1e-12) return;
        const int g = mem[k];
        for (int si : alive[g]) {
          double de = self[g][si];
          for (int y : nbrG[g])
            if (cur[y] >= 0) de += pairVal(g, si, y, cur[y]);
          cur[g] = si;
          dfs(k + 1, e + de);
          cur[g] = -1;
        }
      };
      dfs(0, 0.0);
      for (int g : mem) groups[g].chosen = best[g] >= 0 ? best[g] : 0;
    } else {
      ++rep.greedyComponents;
      for (int g : mem) groups[g].chosen = bestSelf(g);
      for (int sweep = 0; sweep < 50; ++sweep) {
        bool changed = false;
        for (int g : mem) {
          int bs = groups[g].chosen;
          double be = std::numeric_limits<double>::infinity();
          for (int si : alive[g]) {
            double e = self[g][si];
            for (int y : nbrG[g]) e += pairVal(g, si, y, groups[y].chosen);
            if (e < be - 1e-12) {
              be = e;
              bs = si;
            }
          }
          if (bs != groups[g].chosen) {
            groups[g].chosen = bs;
            changed = true;
          }
        }
        if (!changed) break;
      }
    }
  }
  // Refinement of free rotors in 5-degree steps (+-15 degrees), the other groups at their states.
  auto energyOf = [&](int g, const std::vector<AtomView>& v) {
    double e = 0;
    for (const AtomView& a : v) {
      std::set<int> done;
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz) {
            auto it = cells.find(key(a.p, dx, dy, dz));
            if (it == cells.end()) continue;
            for (int k : it->second)
              if (done.insert(k).second) e += sc.pair(a, fixedViews[k]);
          }
      for (int y : nbrG[g])
        for (const AtomView& b : views[y][groups[y].chosen]) e += sc.pair(a, b);
    }
    return e;
  };
  for (int sweep = 0; sweep < 2; ++sweep)
    for (int g = 0; g < G; ++g) {
      OptGroup& og = groups[g];
      if (og.rotorParent < 0) continue;
      const GroupState base = og.states[og.chosen];
      double bestE = energyOf(g, views[g][og.chosen]);
      GroupState best = base;
      std::vector<AtomView> bestV = views[g][og.chosen];
      for (int dk = -3; dk <= 3; ++dk) {
        if (dk == 0) continue;
        GroupState st;
        st.tor = base.tor + 5.0 * dk;
        const auto v = pl.place(og.rotorParent, static_cast<int>(og.moving.size()), st.tor);
        for (std::size_t j = 0; j < og.moving.size(); ++j) st.pos.push_back({og.moving[j], v[j]});
        auto sv = stateViews(st);
        const double e = energyOf(g, sv);
        if (e < bestE - 1e-9) {
          bestE = e;
          best = st;
          bestV = std::move(sv);
        }
      }
      og.states[og.chosen] = best;
      views[g][og.chosen] = std::move(bestV);
    }

  // Apply.
  bool tautomerChanged = false;
  for (OptGroup& g : groups) {
    const GroupState& st = g.states[g.chosen];
    for (const auto& [a, p] : st.pos) pl.xyz()[a] = p;
    for (const auto& [h, p] : st.parent) {
      for (Bond& b : s.bonds) {
        if (b.a == h) b.b = p;
        else if (b.b == h) b.a = p;
      }
      if (parentOf[h] != p) tautomerChanged = true;
    }
    if (g.rotor) ++rep.rotors;
    else ++rep.flipGroups;
    if (st.flipped) rep.flipped.push_back(g.label);
  }
  for (int i = 0; i < nAll; ++i) s.atoms[i].xyz = pl.xyz()[i];
  (void)tautomerChanged;

  // 10. Hydrogen names: CCD names of the parent (matched by the CCD coordinates), then extras.
  if (!planIn) {
    const auto adj = adjacency(s);
    for (int p = 0; p < nHeavy; ++p) {
      std::vector<int> hs;
      for (int x : adj[p])
        if (s.atoms[x].atomicNum == 1) hs.push_back(x);
      if (hs.empty()) continue;
      const int r = s.atoms[p].residue;
      std::vector<std::string> names;
      const CcdComp* comp = residueComp(r);
      std::vector<Vec3> ideal;
      bool haveIdeal = false;
      if (comp) {
        const int cp = comp->find(s.atoms[p].name);
        if (cp >= 0) {
          std::vector<int> ch;
          for (const CcdComp::Bond& b : comp->bonds) {
            const int o = b.a == cp ? b.b : b.b == cp ? b.a : -1;
            if (o >= 0 && comp->atoms[o].atomicNum == 1) ch.push_back(o);
          }
          for (int o : ch) names.push_back(comp->atoms[o].name);
          // Superpose the CCD parent and its heavy neighbours present in the structure.
          std::vector<Vec3> from, to;
          bool ok = comp->atoms[cp].hasXyz;
          if (ok) {
            from.push_back({comp->atoms[cp].xyz[0], comp->atoms[cp].xyz[1], comp->atoms[cp].xyz[2]});
            to.push_back(s.atoms[p].xyz);
          }
          for (int x : adj[p]) {
            if (s.atoms[x].atomicNum == 1 || s.atoms[x].residue != r) continue;
            const int cx = comp->find(s.atoms[x].name);
            if (cx < 0 || !comp->atoms[cx].hasXyz) continue;
            from.push_back({comp->atoms[cx].xyz[0], comp->atoms[cx].xyz[1], comp->atoms[cx].xyz[2]});
            to.push_back(s.atoms[x].xyz);
          }
          // One-neighbour parents (CH3, NH3+, ...): the next atom along the chain fixes the frame.
          if (from.size() == 2) {
            for (int x : adj[p]) {
              if (s.atoms[x].atomicNum == 1 || s.atoms[x].residue != r) continue;
              for (int y : adj[x]) {
                if (y == p || s.atoms[y].atomicNum == 1 || s.atoms[y].residue != r || from.size() >= 3) continue;
                const int cy = comp->find(s.atoms[y].name);
                if (cy < 0 || !comp->atoms[cy].hasXyz) continue;
                from.push_back({comp->atoms[cy].xyz[0], comp->atoms[cy].xyz[1], comp->atoms[cy].xyz[2]});
                to.push_back(s.atoms[y].xyz);
              }
            }
          }
          bool allH = true;
          for (int o : ch) allH = allH && comp->atoms[o].hasXyz;
          if (ok && allH && from.size() >= 3 && hs.size() >= 2 && ch.size() >= 2) {
            // Kabsch via the covariance and an orthonormal frame (Gram-Schmidt of the first two arms).
            auto frame = [](const std::vector<Vec3>& v) {
              const Vec3 e1 = unit(sub(v[1], v[0]));
              Vec3 t = sub(v[2], v[0]);
              const Vec3 e2 = unit(sub(t, mul(e1, dot(t, e1))));
              return std::array<Vec3, 3>{e1, e2, cross(e1, e2)};
            };
            const auto F = frame(from), T = frame(to);
            for (int o : ch) {
              const Vec3 d = sub({comp->atoms[o].xyz[0], comp->atoms[o].xyz[1], comp->atoms[o].xyz[2]}, from[0]);
              const double c0 = dot(d, F[0]), c1 = dot(d, F[1]), c2 = dot(d, F[2]);
              ideal.push_back(add(to[0], add(mul(T[0], c0), add(mul(T[1], c1), mul(T[2], c2)))));
            }
            haveIdeal = true;
          }
        }
      }
      std::vector<std::string> assigned(hs.size());
      std::vector<bool> usedName(names.size(), false);
      if (haveIdeal) {
        // Greedy nearest assignment.
        std::vector<std::tuple<double, std::size_t, std::size_t>> cand;
        for (std::size_t i = 0; i < hs.size(); ++i)
          for (std::size_t j = 0; j < ideal.size(); ++j) cand.push_back({dist(s.atoms[hs[i]].xyz, ideal[j]), i, j});
        std::sort(cand.begin(), cand.end());
        for (const auto& [d, i, j] : cand) {
          (void)d;
          if (!assigned[i].empty() || usedName[j]) continue;
          assigned[i] = names[j];
          usedName[j] = true;
        }
      } else {
        for (std::size_t i = 0, j = 0; i < hs.size() && j < names.size(); ++i, ++j) {
          assigned[i] = names[j];
          usedName[j] = true;
        }
      }
      // Remaining hydrogens: names not used in the residue.
      std::set<std::string> taken;
      if (r >= 0)
        for (int x : s.residues[r].atoms) taken.insert(s.atoms[x].name);
      for (const auto& a : assigned) taken.insert(a);
      const bool named = !s.atoms[p].name.empty();
      // Charged amino-acid N terminus: H1, H2, H3 (PDB convention).
      if (comp && kStdAmino.contains(comp->id) && s.atoms[p].name == "N" && hs.size() == 3) {
        assigned = {"H1", "H2", "H3"};
      }
      for (std::size_t i = 0; i < hs.size(); ++i) {
        if (!assigned[i].empty()) continue;
        if (!named) continue;
        const std::string base = "H" + s.atoms[p].name.substr(std::min<std::size_t>(1, s.atoms[p].name.size()));
        std::string nm;
        for (int k = 1;; ++k) {
          nm = (k == 1 && !taken.contains(base)) ? base : base + std::to_string(k);
          if (!taken.contains(nm)) break;
        }
        assigned[i] = nm;
        taken.insert(nm);
      }
      for (std::size_t i = 0; i < hs.size(); ++i) s.atoms[hs[i]].name = assigned[i];
    }
  }

  // 11. Final charges and Kekulé bonds (the His tautomers may have changed).
  if (auto err = assign(); !err.empty()) {
    errors.insert(errors.end(), err.begin(), err.end());
    return;
  }
  if (planOut) {
    planOut->formalCharge.assign(charge.begin(), charge.end());
    for (int i = 0; i < nHeavy; ++i) planOut->formalCharge[i] = s.atoms[i].formalCharge;
    planOut->bondOrder.clear();
    planOut->hydrogens.clear();
    for (const Bond& b : s.bonds)
      if (b.a < nHeavy && b.b < nHeavy) planOut->bondOrder.push_back(b.order);
    for (int i = nHeavy; i < nAll; ++i)
      for (const Bond& b : s.bonds)
        if (b.a == i || b.b == i) planOut->hydrogens.push_back({b.a == i ? b.b : b.a, s.atoms[i].name});
    planOut->addedOxt = rep.addedOxt;
  }

  // 12. Order: hydrogens after the heavy atoms of their residue (atoms without residue: at the end).
  {
    std::vector<int> order;
    std::vector<std::vector<int>> hOfRes(s.residues.size());
    std::vector<int> loose;
    std::vector<int> par(nAll, -1);
    for (const Bond& b : s.bonds) {
      if (b.b >= nHeavy) par[b.b] = b.a;
      if (b.a >= nHeavy) par[b.a] = b.b;
    }
    std::vector<std::pair<int, int>> hs;
    for (int i = nHeavy; i < nAll; ++i) hs.push_back({par[i], i});
    std::stable_sort(hs.begin(), hs.end());
    for (const auto& [p, h] : hs) {
      if (s.atoms[h].residue >= 0) hOfRes[s.atoms[h].residue].push_back(h);
      else loose.push_back(h);
    }
    std::vector<int> lastOfRes(s.residues.size(), -1);
    for (int i = 0; i < nHeavy; ++i)
      if (s.atoms[i].residue >= 0) lastOfRes[s.atoms[i].residue] = i;
    for (int i = 0; i < nHeavy; ++i) {
      order.push_back(i);
      const int r = s.atoms[i].residue;
      if (r >= 0 && lastOfRes[r] == i) order.insert(order.end(), hOfRes[r].begin(), hOfRes[r].end());
    }
    order.insert(order.end(), loose.begin(), loose.end());
    reorder(s, order, ccdCharge, aromaticCandidate);
  }
}

}  // namespace od
