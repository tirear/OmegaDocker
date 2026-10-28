// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "reprotonate.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <set>

#include "chem.hpp"
#include "log.hpp"
#include "propka.hpp"

namespace od {
namespace {

const std::set<std::string> kStdAmino{"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                                      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec3 unit(const Vec3& a) {
  const double n = std::sqrt(dot(a, a));
  return n > 0 ? mul(a, 1.0 / n) : Vec3{1, 0, 0};
}
// Position d with |cd| = len, angle bcd = ang, dihedral abcd = tor (degrees).
Vec3 placeTorsion(const Vec3& a, const Vec3& b, const Vec3& c, double len, double ang, double tor) {
  const double th = ang * M_PI / 180.0, ph = tor * M_PI / 180.0;
  const Vec3 bc = unit(sub(c, b));
  Vec3 n = cross(sub(b, a), bc);
  if (dot(n, n) < 1e-12) n = cross(std::abs(bc[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}, bc);
  n = unit(n);
  const Vec3 m = cross(n, bc);
  const Vec3 d2{-len * std::cos(th), len * std::sin(th) * std::cos(ph), len * std::sin(th) * std::sin(ph)};
  return add(c, add(add(mul(bc, d2[0]), mul(m, d2[1])), mul(n, d2[2])));
}
double hBondLength(int z) { return z == 7 ? 1.01 : z == 8 ? 0.96 : z == 16 ? 1.34 : 1.0; }

std::vector<std::vector<int>> adjacency(const Molecule& m) {
  std::vector<std::vector<int>> adj(m.atoms.size());
  for (const Bond& b : m.bonds) {
    adj[b.a].push_back(b.b);
    adj[b.b].push_back(b.a);
  }
  return adj;
}
int findInResidue(const Molecule& m, int r, const std::string& name) {
  for (int a : m.residues[r].atoms)
    if (m.atoms[a].name == name) return a;
  return -1;
}
bool metalBonded(const Molecule& m, int a) {
  for (const Bond& b : m.bonds)
    if ((b.a == a || b.b == a) && isMetalBond(m, b)) return true;
  return false;
}
int hydrogensOn(const Molecule& m, const std::vector<std::vector<int>>& adj, int a) {
  int n = 0;
  for (int x : adj[a]) n += m.atoms[x].atomicNum == 1;
  return n;
}
Bond* bondBetween(Molecule& m, int a, int b) {
  for (Bond& x : m.bonds)
    if ((x.a == a && x.b == b) || (x.a == b && x.b == a)) return &x;
  return nullptr;
}
const char* propkaType(SiteType t) {
  switch (t) {
    case SiteType::Asp: return "ASP";
    case SiteType::Glu: return "GLU";
    case SiteType::Cterm: return "C-";
    case SiteType::Nterm: return "N+";
    case SiteType::Cys: return "CYS";
    case SiteType::Tyr: return "TYR";
    case SiteType::Lys: return "LYS";
    case SiteType::Arg: return "ARG";
    case SiteType::His: return "HIS";
  }
  return "";
}
std::string newState(const TitratableSite& s, bool protonated) {
  switch (s.type) {
    case SiteType::Asp: return protonated ? "ASH" : "ASP";
    case SiteType::Glu: return protonated ? "GLH" : "GLU";
    case SiteType::Cterm: return protonated ? "COOH" : "COO-";
    case SiteType::Nterm: return protonated ? "NH3+" : "NH2";
    case SiteType::Cys: return protonated ? "CYS" : "CYM";
    case SiteType::Tyr: return protonated ? "TYR" : "TYM";
    case SiteType::Lys: return protonated ? "LYS" : "LYN";
    case SiteType::Arg: return protonated ? "ARG" : "ARN";
    case SiteType::His: return protonated ? "HIP" : (s.state == "HIP" ? "HIE" : s.state);
  }
  return s.state;
}

}  // namespace

bool protonatedState(const std::string& s) {
  return s == "ASH" || s == "GLH" || s == "COOH" || s == "NH3+" || s == "CYS" || s == "TYR" || s == "LYS" || s == "ARG" ||
         s == "HIP";
}

std::vector<TitratableSite> titratableSites(const Molecule& m, const std::vector<ProtonationOverride>& overrides) {
  std::vector<TitratableSite> out;
  if (!m.hasResidueInfo) return out;
  const auto adj = adjacency(m);
  for (std::size_t ri = 0; ri < m.residues.size(); ++ri) {
    const int r = static_cast<int>(ri);
    const ResidueInfo& info = m.residues[r].info;
    if (!kStdAmino.contains(info.name)) continue;
    const std::string res = std::format("{}:{}{}", info.chain, info.seq, info.iCode);
    bool overridden = false;
    for (const auto& o : overrides)
      if (o.site == res || o.site.starts_with(res + ":")) overridden = true;
    auto at = [&](const char* n) { return findInResidue(m, r, n); };
    auto add = [&](SiteType t, int key, std::vector<int> atoms, std::string state) {
      TitratableSite s;
      s.residue = r;
      s.type = t;
      s.key = key;
      s.state = std::move(state);
      s.label = res + ":" + info.name;
      s.locked = overridden;
      for (int a : atoms) s.locked = s.locked || metalBonded(m, a);
      out.push_back(std::move(s));
    };
    auto carboxyl = [&](SiteType t, int c, int o1, int o2, const char* deprot, const char* prot) {
      if (c < 0 || o1 < 0 || o2 < 0) return;
      add(t, t == SiteType::Cterm ? o2 : c, {o1, o2}, hydrogensOn(m, adj, o1) + hydrogensOn(m, adj, o2) > 0 ? prot : deprot);
    };
    if (info.name == "ASP") carboxyl(SiteType::Asp, at("CG"), at("OD1"), at("OD2"), "ASP", "ASH");
    if (info.name == "GLU") carboxyl(SiteType::Glu, at("CD"), at("OE1"), at("OE2"), "GLU", "GLH");
    if (at("OXT") >= 0) carboxyl(SiteType::Cterm, at("C"), at("O"), at("OXT"), "COO-", "COOH");
    if (const int n = at("N"); n >= 0) {
      bool linked = false;
      int heavy = 0;
      for (int x : adj[n]) {
        if (m.atoms[x].atomicNum == 1 || isMetal(m.atoms[x].atomicNum)) continue;
        ++heavy;
        if (m.atoms[x].residue != r) linked = true;
      }
      if (!linked) add(SiteType::Nterm, n, {n}, hydrogensOn(m, adj, n) + heavy == 4 ? "NH3+" : "NH2");
    }
    if (info.name == "CYS")
      if (const int sg = at("SG"); sg >= 0) {
        bool bridged = false;
        for (int x : adj[sg]) bridged = bridged || m.atoms[x].atomicNum == 16;
        if (!bridged) add(SiteType::Cys, sg, {sg}, hydrogensOn(m, adj, sg) > 0 ? "CYS" : "CYM");
      }
    if (info.name == "TYR")
      if (const int oh = at("OH"); oh >= 0) add(SiteType::Tyr, oh, {oh}, hydrogensOn(m, adj, oh) > 0 ? "TYR" : "TYM");
    if (info.name == "LYS")
      if (const int nz = at("NZ"); nz >= 0) add(SiteType::Lys, nz, {nz}, hydrogensOn(m, adj, nz) == 3 ? "LYS" : "LYN");
    if (info.name == "ARG") {
      const int cz = at("CZ"), ne = at("NE"), n1 = at("NH1"), n2 = at("NH2");
      if (cz >= 0 && ne >= 0 && n1 >= 0 && n2 >= 0) {
        int q = 0;
        for (int a : {ne, n1, n2}) q += m.atoms[a].formalCharge;
        add(SiteType::Arg, cz, {ne, n1, n2}, q > 0 ? "ARG" : "ARN");
      }
    }
    if (info.name == "HIS") {
      const int cg = at("CG"), nd1 = at("ND1"), ne2 = at("NE2"), cd2 = at("CD2"), ce1 = at("CE1");
      if (cg >= 0 && nd1 >= 0 && ne2 >= 0 && cd2 >= 0 && ce1 >= 0) {
        const bool h1 = hydrogensOn(m, adj, nd1) > 0, h2 = hydrogensOn(m, adj, ne2) > 0;
        if (h1 || h2) add(SiteType::His, cg, {nd1, ne2}, h1 && h2 ? "HIP" : h1 ? "HID" : "HIE");
      }
    }
  }
  return out;
}

std::vector<PkaMap> calculatePka(const std::vector<PkaMolecule>& mols) {
  std::vector<propka::InputAtom> in;
  std::vector<std::pair<int, int>> origin;  // (molecule, atom)
  for (std::size_t k = 0; k < mols.size(); ++k) {
    const Molecule& m = *mols[k].mol;
    for (std::size_t i = 0; i < m.atoms.size(); ++i) {
      const Atom& a = m.atoms[i];
      if (a.atomicNum == 1) continue;
      propka::InputAtom x;
      x.name = a.name.empty() ? std::format("{}{}", elementSymbol(a.atomicNum), i + 1) : a.name;
      if (a.residue >= 0 && m.hasResidueInfo) {
        const ResidueInfo& ri = m.residues[a.residue].info;
        x.resName = ri.name;
        x.chain = ri.chain;
        x.resNum = ri.seq;
        x.iCode = ri.iCode;
        x.hetatm = a.hetatm;
      } else {
        x.resName = "LIG";
        x.chain = std::format("L{}", k);
        x.resNum = 1;
        x.hetatm = true;
      }
      x.element = elementSymbol(a.atomicNum);
      x.xyz = mols[k].xyz[i];
      in.push_back(std::move(x));
      origin.push_back({static_cast<int>(k), static_cast<int>(i)});
    }
  }
  std::vector<PkaMap> out(mols.size());
  for (const auto& g : propka::calculate(in)) {
    const auto [k, i] = origin[g.atom];
    out[k][{i, g.residueType}] = g.pka;
  }
  return out;
}

std::vector<StateChange> stateChanges(const std::vector<TitratableSite>& sites, const PkaMap& pka, double pH) {
  std::vector<StateChange> out;
  for (std::size_t k = 0; k < sites.size(); ++k) {
    const TitratableSite& s = sites[k];
    if (s.locked) continue;
    auto it = pka.find({s.key, propkaType(s.type)});
    if (it == pka.end()) continue;
    const std::string st = newState(s, it->second >= pH);
    if (st != s.state) out.push_back({static_cast<int>(k), st, it->second});
  }
  return out;
}

double protonationFreeEnergy(const std::vector<TitratableSite>& sites, const std::vector<StateChange>& changes, double pH,
                             double temperature) {
  const double c = std::log(10.0) * 0.0019872036 * temperature;
  double g = 0;
  for (const auto& ch : changes) {
    const bool before = protonatedState(sites[ch.site].state), after = protonatedState(ch.state);
    if (after && !before) g += c * (pH - ch.pka);
    else if (before && !after) g += c * (ch.pka - pH);
  }
  return g;
}

std::string describeChanges(const std::vector<TitratableSite>& sites, const std::vector<StateChange>& changes) {
  std::string s;
  for (const auto& c : changes) s += std::format("{}{} {}->{}", s.empty() ? "" : " ", sites[c.site].label, sites[c.site].state, c.state);
  return s;
}

MoleculeVariant makeVariant(const Molecule& base, const TorsionTree& tree, const std::vector<TitratableSite>& sites,
                            const std::vector<StateChange>& changes) {
  Molecule m = base;
  const auto adj = adjacency(base);
  std::vector<bool> removed(base.atoms.size(), false);
  struct NewH {
    int parent;
    Vec3 xyz;
    std::string name;
  };
  std::vector<NewH> added;
  auto hsOf = [&](int a) {
    std::vector<int> v;
    for (int x : adj[a])
      if (base.atoms[x].atomicNum == 1 && !removed[x]) v.push_back(x);
    return v;
  };
  // Name of a hydrogen added to atom a: "H" + name without its first character (+ index when several).
  auto hName = [&](int a, int count, int index) {
    const std::string& n = base.atoms[a].name;
    const std::string stem = "H" + (n.size() > 1 ? n.substr(1) : std::string());
    return count > 1 ? stem + std::to_string(index) : stem;
  };
  auto addH = [&](int a, const std::string& name) {
    std::vector<Vec3> nb;
    int heavyNb = -1;
    for (int x : adj[a]) {
      if (removed[x] || isMetal(base.atoms[x].atomicNum)) continue;
      nb.push_back(base.atoms[x].xyz);
      if (base.atoms[x].atomicNum != 1) heavyNb = x;
    }
    for (const auto& h : added)
      if (h.parent == a) nb.push_back(h.xyz);
    const Vec3 c = base.atoms[a].xyz;
    const int z = base.atoms[a].atomicNum;
    const double len = hBondLength(z);
    Vec3 p;
    if (nb.size() >= 2 || heavyNb < 0) {
      Vec3 s{0, 0, 0};
      for (const Vec3& x : nb) s = add(s, unit(sub(x, c)));
      p = add(c, mul(unit(mul(s, -1.0)), len));
    } else {
      // One neighbour: torsion placement; sp2 neighbour (carboxyl, Tyr): syn / in plane (0 degrees) to a
      // double-bonded or ring neighbour, otherwise 180 degrees to a heavy neighbour of the neighbour.
      int ref = -1;
      bool sp2 = false;
      for (int y : adj[heavyNb]) {
        if (y == a || base.atoms[y].atomicNum == 1) continue;
        const Bond* b = nullptr;
        for (const Bond& x : base.bonds)
          if ((x.a == heavyNb && x.b == y) || (x.b == heavyNb && x.a == y)) b = &x;
        const bool dbl = b && (b->order == BondOrder::Double || b->aromatic);
        if (ref < 0 || (dbl && !sp2)) {
          ref = y;
          sp2 = dbl;
        }
      }
      const Vec3 r = ref >= 0 ? base.atoms[ref].xyz : add(base.atoms[heavyNb].xyz, Vec3{1, 0, 0});
      p = placeTorsion(r, base.atoms[heavyNb].xyz, c, len, z == 16 ? 96.0 : 109.47, sp2 ? 0.0 : 180.0);
    }
    added.push_back({a, p, name});
  };
  auto removeOneH = [&](int a) {
    auto hs = hsOf(a);
    if (hs.empty()) return;
    std::sort(hs.begin(), hs.end(), [&](int x, int y) { return base.atoms[x].name < base.atoms[y].name; });
    removed[hs.back()] = true;
  };

  for (const auto& ch : changes) {
    const TitratableSite& s = sites[ch.site];
    const int r = s.residue;
    auto at = [&](const char* n) { return findInResidue(base, r, n); };
    const bool gain = protonatedState(ch.state) && !protonatedState(s.state);
    switch (s.type) {
      case SiteType::Asp:
      case SiteType::Glu:
      case SiteType::Cterm: {
        const int o1 = s.type == SiteType::Asp ? at("OD1") : s.type == SiteType::Glu ? at("OE1") : at("O");
        const int o2 = s.type == SiteType::Asp ? at("OD2") : s.type == SiteType::Glu ? at("OE2") : at("OXT");
        if (gain) {
          const int o = base.atoms[o2].formalCharge < 0 ? o2 : base.atoms[o1].formalCharge < 0 ? o1 : o2;
          m.atoms[o].formalCharge = 0;
          addH(o, hName(o, 1, 1));
        } else {
          for (int o : {o1, o2})
            if (!hsOf(o).empty()) {
              removeOneH(o);
              m.atoms[o].formalCharge = -1;
            }
        }
        break;
      }
      case SiteType::Nterm:
      case SiteType::Lys: {
        const int n = s.key;
        if (gain) {
          m.atoms[n].formalCharge = 1;
          addH(n, "");
        } else {
          removeOneH(n);
          m.atoms[n].formalCharge = 0;
        }
        // Names: N terminus H1, H2, H3 (charged) or H, H2 (neutral); Lys HZ1, HZ2, HZ3.
        break;
      }
      case SiteType::Cys:
      case SiteType::Tyr: {
        const int a = s.key;
        if (gain) {
          m.atoms[a].formalCharge = 0;
          addH(a, hName(a, 1, 1));
        } else {
          removeOneH(a);
          m.atoms[a].formalCharge = -1;
        }
        break;
      }
      case SiteType::Arg: {
        const int cz = at("CZ");
        if (gain) {
          int n = -1;
          for (const char* nm : {"NH2", "NH1", "NE"}) {
            const int x = at(nm);
            const Bond* b = nullptr;
            for (const Bond& y : base.bonds)
              if ((y.a == x && y.b == cz) || (y.b == x && y.a == cz)) b = &y;
            if (n < 0 && b && b->order == BondOrder::Double) n = x;
          }
          if (n < 0) fail("{}: guanidine without a C=N double bond", s.label);
          m.atoms[n].formalCharge = 1;
          addH(n, "");
        } else {
          for (const char* nm : {"NE", "NH1", "NH2"})
            if (const int x = at(nm); base.atoms[x].formalCharge > 0) {
              removeOneH(x);
              m.atoms[x].formalCharge = 0;
              break;
            }
        }
        break;
      }
      case SiteType::His: {
        const int cg = at("CG"), nd1 = at("ND1"), cd2 = at("CD2"), ce1 = at("CE1"), ne2 = at("NE2");
        const bool hd = ch.state == "HID" || ch.state == "HIP", he = ch.state == "HIE" || ch.state == "HIP";
        if (!hd) {
          if (!hsOf(nd1).empty()) removeOneH(nd1);
        }
        if (!he) {
          if (!hsOf(ne2).empty()) removeOneH(ne2);
        }
        if (hd && hsOf(nd1).empty()) addH(nd1, "HD1");
        if (he && hsOf(ne2).empty()) addH(ne2, "HE2");
        auto setOrder = [&](int a, int b, BondOrder o) {
          if (Bond* x = bondBetween(m, a, b)) x->order = o;
        };
        setOrder(cg, cd2, BondOrder::Double);
        setOrder(cg, nd1, BondOrder::Single);
        setOrder(cd2, ne2, BondOrder::Single);
        const bool ndDouble = ch.state != "HID";
        setOrder(nd1, ce1, ndDouble ? BondOrder::Double : BondOrder::Single);
        setOrder(ce1, ne2, ndDouble ? BondOrder::Single : BondOrder::Double);
        m.atoms[nd1].formalCharge = ch.state == "HIP" ? 1 : 0;
        m.atoms[ne2].formalCharge = 0;
        break;
      }
    }
  }

  // Hydrogen names of the changed N atoms (N terminus, Lys, Arg).
  std::map<int, std::vector<int>> nameGroups;  // parent -> kept base hydrogens
  for (const auto& ch : changes) {
    const TitratableSite& s = sites[ch.site];
    std::vector<int> parents;
    if (s.type == SiteType::Nterm || s.type == SiteType::Lys) parents = {s.key};
    if (s.type == SiteType::Arg)
      for (const char* nm : {"NE", "NH1", "NH2"}) parents.push_back(findInResidue(base, s.residue, nm));
    for (int p : parents) {
      bool any = false;
      for (const auto& h : added) any = any || h.parent == p;
      if (!any && !std::ranges::any_of(adj[p], [&](int x) { return removed[x]; })) continue;
      const auto hs = hsOf(p);
      int n = static_cast<int>(hs.size());
      for (const auto& h : added) n += h.parent == p;
      const bool nterm = s.type == SiteType::Nterm && p == s.key;
      const bool charged = m.atoms[p].formalCharge > 0;
      std::vector<std::string> names;
      for (int k = 1; k <= n; ++k)
        names.push_back(nterm ? (charged ? std::format("H{}", k) : (k == 1 ? std::string("H") : std::format("H{}", k)))
                              : hName(p, n, k));
      int k = 0;
      for (int h : hs) m.atoms[h].name = names[k++];
      for (auto& h : added)
        if (h.parent == p) h.name = names[k++];
    }
  }

  // Atom order: kept atoms; added hydrogens after the last atom of their residue.
  std::vector<int> lastOfResidue(base.residues.size(), -1);
  for (std::size_t i = 0; i < base.atoms.size(); ++i)
    if (base.atoms[i].residue >= 0) lastOfResidue[base.atoms[i].residue] = static_cast<int>(i);
  MoleculeVariant v;
  Molecule& out = v.mol;
  out.id = base.id;
  out.role = base.role;
  out.sourceFile = base.sourceFile;
  out.hasResidueInfo = base.hasResidueInfo;
  out.chargesFromFile = false;
  std::vector<int> map(base.atoms.size(), -1);
  std::vector<std::vector<int>> addedAfter(base.atoms.size());
  for (std::size_t k = 0; k < added.size(); ++k) {
    const int r = base.atoms[added[k].parent].residue;
    addedAfter[r >= 0 ? lastOfResidue[r] : added[k].parent].push_back(static_cast<int>(k));
  }
  std::vector<int> addedIndex(added.size(), -1);
  for (std::size_t i = 0; i < base.atoms.size(); ++i) {
    if (!removed[i]) {
      map[i] = static_cast<int>(out.atoms.size());
      out.atoms.push_back(m.atoms[i]);
      v.toBase.push_back(static_cast<int>(i));
    }
    for (int k : addedAfter[i]) {
      Atom h;
      h.name = added[k].name;
      h.atomicNum = 1;
      h.xyz = added[k].xyz;
      h.residue = base.atoms[added[k].parent].residue;
      h.hetatm = base.atoms[added[k].parent].hetatm;
      addedIndex[k] = static_cast<int>(out.atoms.size());
      out.atoms.push_back(h);
      v.toBase.push_back(-1);
    }
  }
  for (const Bond& b : m.bonds)
    if (map[b.a] >= 0 && map[b.b] >= 0) out.bonds.push_back({map[b.a], map[b.b], b.order, false});
  for (std::size_t k = 0; k < added.size(); ++k) out.bonds.push_back({map[added[k].parent], addedIndex[k], BondOrder::Single, false});
  out.residues = base.residues;
  for (auto& r : out.residues) r.atoms.clear();
  for (std::size_t i = 0; i < out.atoms.size(); ++i)
    if (out.atoms[i].residue >= 0) out.residues[out.atoms[i].residue].atoms.push_back(static_cast<int>(i));
  perceiveChemistry(out);

  // Torsion tree: coordinate indices remapped; virtual atoms after the atoms.
  const int nBase = tree.atomCount, nVar = static_cast<int>(out.atoms.size());
  auto coord = [&](int c) { return c < 0 ? -1 : c >= nBase ? nVar + (c - nBase) : map[c]; };
  TorsionTree& t = v.tree;
  t = tree;
  t.atomCount = nVar;
  t.anchors.clear();
  for (int a : tree.anchors)
    if (map[a] >= 0) t.anchors.push_back(map[a]);
  for (auto& tt : t.torsions) {
    tt.fixedAtom = coord(tt.fixedAtom);
    tt.movingAtom = coord(tt.movingAtom);
    std::vector<int> mv;
    std::set<int> movedBase(tt.moving.begin(), tt.moving.end());
    for (int c : tt.moving)
      if (const int x = coord(c); x >= 0) mv.push_back(x);
    for (std::size_t k = 0; k < added.size(); ++k)
      if (movedBase.contains(added[k].parent)) mv.push_back(addedIndex[k]);
    std::sort(mv.begin(), mv.end());
    tt.moving = std::move(mv);
  }
  for (auto& va : t.virtuals) {
    va.parent = coord(va.parent);
    va.copyOf = coord(va.copyOf);
  }
  for (auto& r : t.restraints) {
    r.a = coord(r.a);
    r.b = coord(r.b);
  }
  for (auto& l : t.links) {
    l.a = coord(l.a);
    l.b = coord(l.b);
    l.y = coord(l.y);
    l.copyB = coord(l.copyB);
    l.copyY = coord(l.copyY);
    l.copyA = coord(l.copyA);
  }
  return v;
}

ProtonationVariants::ProtonationVariants(const EnergySystem& base, std::vector<EnergyPart> parts,
                                         const std::vector<std::vector<ProtonationOverride>>& overrides, const EnergyParams& ep)
    : base_(base), parts_(std::move(parts)), ep_(ep) {
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    sites_.push_back(titratableSites(*parts_[k].mol, k < overrides.size() ? overrides[k] : std::vector<ProtonationOverride>{}));
    for (const auto& s : sites_.back()) active_ = active_ || !s.locked;
  }
}

std::vector<PkaMolecule> ProtonationVariants::pkaMolecules(std::span<const Vec3f> x) const {
  std::vector<PkaMolecule> out;
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    PkaMolecule p;
    p.mol = parts_[k].mol;
    const int off = base_.coordOffset(static_cast<int>(k));
    for (std::size_t i = 0; i < p.mol->atoms.size(); ++i) {
      const Vec3f& q = x[off + i];
      p.xyz.push_back({q[0], q[1], q[2]});
    }
    out.push_back(std::move(p));
  }
  return out;
}

ProtonationVariants::Key ProtonationVariants::keyOf(const std::vector<std::vector<StateChange>>& changes) const {
  Key key(parts_.size());
  for (std::size_t k = 0; k < changes.size() && k < parts_.size(); ++k)
    for (const auto& c : changes[k]) key[k].push_back({c.site, c.state});
  for (auto& v : key) std::sort(v.begin(), v.end());
  return key;
}

const ProtonationVariants::Variant& ProtonationVariants::entry(const Key& key) {
  std::lock_guard lock(mutex_);
  if (auto it = cache_.find(key); it != cache_.end()) return *it->second;
  auto var = std::make_unique<Variant>();
  var->changes.resize(parts_.size());
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    for (const auto& [site, state] : key[k]) var->changes[k].push_back({site, state, 0.0});
    if (var->changes[k].empty()) {
      var->mols.push_back(nullptr);
      var->parts.push_back(parts_[k]);
      continue;
    }
    var->mols.push_back(std::make_unique<MoleculeVariant>(makeVariant(*parts_[k].mol, *parts_[k].tree, sites_[k], var->changes[k])));
    var->parts.push_back({&var->mols.back()->mol, &var->mols.back()->tree, parts_[k].rings, &var->mols.back()->toBase});
  }
  var->sys = std::make_unique<EnergySystem>(var->parts, ep_);
  if (var->sys->torsionCount() != base_.torsionCount())
    fail("protonation variant: {} torsions, base system {}", var->sys->torsionCount(), base_.torsionCount());
  return *cache_.emplace(key, std::move(var)).first->second;
}

int ProtonationVariants::variantOf(const std::vector<std::vector<StateChange>>& changes) {
  const Key key = keyOf(changes);
  if (std::all_of(key.begin(), key.end(), [](const auto& v) { return v.empty(); })) return 0;
  const Variant& var = entry(key);
  std::lock_guard lock(mutex_);
  if (auto it = index_.find(key); it != index_.end()) return it->second;
  variants_.push_back(&var);
  const int id = static_cast<int>(variants_.size());
  index_[key] = id;
  return id;
}

const EnergySystem& ProtonationVariants::sysOf(const std::vector<std::vector<StateChange>>& changes) {
  const Key key = keyOf(changes);
  if (std::all_of(key.begin(), key.end(), [](const auto& v) { return v.empty(); })) return base_;
  return *entry(key).sys;
}

const std::vector<EnergyPart>& ProtonationVariants::partsOf(const std::vector<std::vector<StateChange>>& changes) {
  const Key key = keyOf(changes);
  if (std::all_of(key.begin(), key.end(), [](const auto& v) { return v.empty(); })) return parts_;
  return entry(key).parts;
}

int ProtonationVariants::chooseTautomers(std::vector<std::vector<StateChange>>& changes,
                                         const TautomerEnergy& energy) {
  changes.resize(parts_.size());
  // State of site k of part p in changes, and setting it (no entry when it is the state of the molecule;
  // an entry keeps its pKa).
  auto stateOf = [&](std::size_t p, int k) {
    for (const auto& c : changes[p])
      if (c.site == k) return c.state;
    return sites_[p][k].state;
  };
  auto setState = [&](std::size_t p, int k, const std::string& st) {
    auto& v = changes[p];
    auto it = std::find_if(v.begin(), v.end(), [&](const StateChange& c) { return c.site == k; });
    if (st == sites_[p][k].state) {
      if (it != v.end()) v.erase(it);
      return;
    }
    if (it != v.end()) {
      it->state = st;
      return;
    }
    v.insert(std::lower_bound(v.begin(), v.end(), k, [](const StateChange& c, int s) { return c.site < s; }), {k, st, 0.0});
  };
  std::vector<std::pair<std::size_t, int>> his;
  for (std::size_t p = 0; p < parts_.size(); ++p)
    for (std::size_t k = 0; k < sites_[p].size(); ++k) {
      const TitratableSite& s = sites_[p][k];
      const std::string st = stateOf(p, static_cast<int>(k));
      if (s.type == SiteType::His && !s.locked && (st == "HID" || st == "HIE")) his.emplace_back(p, static_cast<int>(k));
    }
  if (his.empty()) return 0;
  int flips = 0;
  std::map<Key, double> seen;
  auto eval = [&] {
    const Key key = keyOf(changes);
    if (auto it = seen.find(key); it != seen.end()) return it->second;
    return seen[key] = energy(sysOf(changes), partsOf(changes));
  };
  double best = eval();
  for (int pass = 0; pass < kTautomerPasses; ++pass) {
    bool changed = false;
    for (const auto& [p, k] : his) {
      const std::string cur = stateOf(p, k), alt = cur == "HIE" ? "HID" : "HIE";
      setState(p, k, alt);
      const double e = eval();
      if (e < best) {
        best = e;
        changed = true;
        ++flips;
      } else {
        setState(p, k, cur);
      }
    }
    if (!changed) break;
  }
  return flips;
}

std::size_t ProtonationVariants::count() const {
  std::lock_guard lock(mutex_);
  return variants_.size() + 1;
}

const EnergySystem& ProtonationVariants::sys(int id) const {
  if (id == 0) return base_;
  std::lock_guard lock(mutex_);
  return *variants_[id - 1]->sys;
}

std::vector<EnergyPart> ProtonationVariants::parts(int id) const {
  if (id == 0) return parts_;
  std::lock_guard lock(mutex_);
  return variants_[id - 1]->parts;
}

std::vector<std::vector<StateChange>> ProtonationVariants::changes(int id) const {
  if (id == 0) return std::vector<std::vector<StateChange>>(parts_.size());
  std::lock_guard lock(mutex_);
  return variants_[id - 1]->changes;
}

std::string ProtonationVariants::describe(int id) const {
  const auto ch = changes(id);
  std::string s;
  for (std::size_t k = 0; k < ch.size(); ++k) {
    const std::string d = describeChanges(sites_[k], ch[k]);
    if (!d.empty()) s += (s.empty() ? "" : " ") + d;
  }
  return s;
}

}  // namespace od
