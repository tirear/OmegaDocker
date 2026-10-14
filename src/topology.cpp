// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "topology.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <numbers>

#include "log.hpp"

namespace od {
namespace {

constexpr double kCollinearDeg = 170.0;

struct Edge {
  int to;
  int bond;
};

struct Graph {
  const Molecule& m;
  std::vector<std::vector<Edge>> adj;
  std::vector<bool> ringBond;  // ring bonds of the whole bond graph (backbone detection)
  std::vector<bool> ringAtom;
  std::vector<bool> classRing;  // ring bonds after removing broken bonds (classification)
  std::map<std::pair<int, int>, int> bondOf;

  explicit Graph(const Molecule& mol) : m(mol), adj(mol.atoms.size()), ringAtom(mol.atoms.size(), false) {
    for (std::size_t i = 0; i < m.bonds.size(); ++i) {
      const Bond& b = m.bonds[i];
      adj[b.a].push_back({b.b, static_cast<int>(i)});
      adj[b.b].push_back({b.a, static_cast<int>(i)});
      bondOf[{std::min(b.a, b.b), std::max(b.a, b.b)}] = static_cast<int>(i);
    }
    ringBond = markBridges({});
    for (std::size_t i = 0; i < m.bonds.size(); ++i)
      if (ringBond[i]) ringAtom[m.bonds[i].a] = ringAtom[m.bonds[i].b] = true;
    classRing = ringBond;
  }

  void setBroken(const std::vector<int>& broken) {
    std::vector<bool> skip(m.bonds.size(), false);
    for (int b : broken) skip[b] = true;
    classRing = markBridges(skip);
    for (int b : broken) classRing[b] = false;
  }

  // A bond lies in a ring if and only if it is not a bridge (iterative Tarjan).
  // Bonds marked in skip are treated as absent.
  std::vector<bool> markBridges(const std::vector<bool>& skip) const {
    std::vector<bool> ring(m.bonds.size(), true);
    const int n = static_cast<int>(adj.size());
    std::vector<int> disc(n, -1), low(n, 0);
    int t = 0;
    struct Frame {
      int v, parentBond;
      std::size_t next;
    };
    for (int s = 0; s < n; ++s) {
      if (disc[s] >= 0) continue;
      std::vector<Frame> st{{s, -1, 0}};
      disc[s] = low[s] = t++;
      while (!st.empty()) {
        Frame& f = st.back();
        if (f.next < adj[f.v].size()) {
          const Edge e = adj[f.v][f.next++];
          if (e.bond == f.parentBond || (!skip.empty() && skip[e.bond])) continue;
          if (disc[e.to] < 0) {
            disc[e.to] = low[e.to] = t++;
            st.push_back({e.to, e.bond, 0});
          } else {
            low[f.v] = std::min(low[f.v], disc[e.to]);
          }
        } else {
          const Frame done = f;
          st.pop_back();
          if (!st.empty()) {
            low[st.back().v] = std::min(low[st.back().v], low[done.v]);
            if (low[done.v] > disc[st.back().v]) ring[done.parentBond] = false;
          }
        }
      }
    }
    return ring;
  }

  int z(int i) const { return m.atoms[i].atomicNum; }
  bool heavy(int i) const { return z(i) > 1; }
  int bond(int a, int b) const {
    auto it = bondOf.find({std::min(a, b), std::max(a, b)});
    return it == bondOf.end() ? -1 : it->second;
  }
  int heavyDegree(int i) const {
    int d = 0;
    for (const Edge& e : adj[i]) d += heavy(e.to);
    return d;
  }
  int hydrogenCount(int i) const {
    int h = 0;
    for (const Edge& e : adj[i]) h += !heavy(e.to);
    return h;
  }
  bool isDouble(int bondIdx) const {
    return m.bonds[bondIdx].order == BondOrder::Double && !m.bonds[bondIdx].aromatic;
  }
};

double angleDeg(const Vec3& a, const Vec3& b, const Vec3& c) {
  double u[3], v[3], uu = 0, vv = 0, uv = 0;
  for (int k = 0; k < 3; ++k) {
    u[k] = a[k] - b[k];
    v[k] = c[k] - b[k];
    uu += u[k] * u[k];
    vv += v[k] * v[k];
    uv += u[k] * v[k];
  }
  if (uu == 0 || vv == 0) return 0;
  return std::acos(std::clamp(uv / std::sqrt(uu * vv), -1.0, 1.0)) * 180.0 / std::numbers::pi;
}

// ---------- backbone detection (spec 4.1) ----------

std::vector<BackboneUnit> linkUnits(std::vector<BackboneUnit> cand, int linkFrom, int linkTo, const Graph& g) {
  // Keep candidates whose atoms are not used by an earlier candidate.
  std::vector<bool> used(g.adj.size(), false);
  std::vector<BackboneUnit> units;
  for (const BackboneUnit& u : cand) {
    bool clash = false;
    for (int a : u.atoms) clash |= a >= 0 && used[a];
    if (clash) continue;
    for (int a : u.atoms)
      if (a >= 0) used[a] = true;
    units.push_back(u);
  }
  std::map<int, int> unitOfLinkTo;
  for (std::size_t i = 0; i < units.size(); ++i)
    if (units[i].atoms[linkTo] >= 0) unitOfLinkTo[units[i].atoms[linkTo]] = static_cast<int>(i);
  for (std::size_t i = 0; i < units.size(); ++i)
    for (const Edge& e : g.adj[units[i].atoms[linkFrom]]) {
      auto it = unitOfLinkTo.find(e.to);
      if (it != unitOfLinkTo.end() && it->second != static_cast<int>(i)) {
        units[i].next = it->second;
        units[it->second].prev = static_cast<int>(i);
      }
    }
  // A backbone needs at least two linked units.
  std::vector<int> newIndex(units.size(), -1);
  std::vector<BackboneUnit> kept;
  for (std::size_t i = 0; i < units.size(); ++i)
    if (units[i].prev >= 0 || units[i].next >= 0) {
      newIndex[i] = static_cast<int>(kept.size());
      kept.push_back(units[i]);
    }
  for (BackboneUnit& u : kept) {
    if (u.prev >= 0) u.prev = newIndex[u.prev];
    if (u.next >= 0) u.next = newIndex[u.next];
  }
  return kept;
}

std::vector<BackboneUnit> findProteinUnits(const Graph& g) {
  std::vector<BackboneUnit> cand;
  for (int c = 0; c < static_cast<int>(g.adj.size()); ++c) {
    if (g.z(c) != 6) continue;
    int o = -1;
    for (const Edge& e : g.adj[c])
      if (g.z(e.to) == 8 && g.heavyDegree(e.to) == 1 && g.isDouble(e.bond)) o = e.to;
    if (o < 0) continue;
    for (const Edge& e1 : g.adj[c]) {
      if (g.z(e1.to) != 6) continue;
      for (const Edge& e2 : g.adj[e1.to])
        if (g.z(e2.to) == 7) {
          BackboneUnit u;
          u.atoms = {e2.to, e1.to, c, o, -1, -1};
          cand.push_back(u);
        }
    }
  }
  auto units = linkUnits(std::move(cand), 2, 0, g);  // C(i) - N(i+1)
  for (BackboneUnit& u : units) u.residue = g.m.atoms[u.atoms[1]].residue;
  return units;
}

std::vector<BackboneUnit> findNucleicUnits(const Graph& g) {
  std::vector<BackboneUnit> cand;
  const int n = static_cast<int>(g.adj.size());
  for (int c3 = 0; c3 < n; ++c3) {
    if (g.z(c3) != 6 || !g.ringAtom[c3]) continue;
    for (const Edge& eo3 : g.adj[c3]) {
      const int o3 = eo3.to;
      if (g.z(o3) != 8 || g.ringAtom[o3]) continue;
      for (const Edge& e4 : g.adj[c3]) {
        const int c4 = e4.to;
        if (g.z(c4) != 6 || !g.ringBond[e4.bond]) continue;
        for (const Edge& e5 : g.adj[c4]) {
          const int c5 = e5.to;
          if (g.z(c5) != 6 || g.ringAtom[c5]) continue;
          for (const Edge& eo5 : g.adj[c5]) {
            const int o5 = eo5.to;
            if (g.z(o5) != 8 || g.ringAtom[o5]) continue;
            int p = -1;
            for (const Edge& ep : g.adj[o5])
              if (g.z(ep.to) == 15) p = ep.to;
            BackboneUnit u;
            u.atoms = {p, o5, c5, c4, c3, o3};
            cand.push_back(u);
          }
        }
      }
    }
  }
  auto units = linkUnits(std::move(cand), 5, 0, g);  // O3'(i) - P(i+1)
  for (BackboneUnit& u : units) u.residue = g.m.atoms[u.atoms[3]].residue;
  return units;
}

// ---------- torsion classification (spec 4.4) ----------

bool collinearEnd(const Graph& g, int end, int other) {
  bool any = false;
  for (const Edge& e : g.adj[end]) {
    if (e.to == other || !g.heavy(e.to)) continue;
    any = true;
    if (angleDeg(g.m.atoms[e.to].xyz, g.m.atoms[end].xyz, g.m.atoms[other].xyz) <= kCollinearDeg) return false;
  }
  return any;
}

// end is a carbon with a double bond to O, S or N (not other); returns that atom or -1.
int carbonylLike(const Graph& g, int end, int other) {
  if (g.z(end) != 6) return -1;
  for (const Edge& e : g.adj[end]) {
    if (e.to == other || !g.isDouble(e.bond)) continue;
    const int zx = g.z(e.to);
    if (zx == 8 || zx == 16 || zx == 7) return e.to;
  }
  return -1;
}

// All neighbors of end except other are equivalent terminal atoms (at least two).
bool symmetricEnd(const Graph& g, int end, int other) {
  int count = 0, zKey = -1, hKey = -1;
  for (const Edge& e : g.adj[end]) {
    if (e.to == other) continue;
    if (g.heavy(e.to) && g.heavyDegree(e.to) != 1) return false;
    const int zz = g.z(e.to), hh = g.heavy(e.to) ? g.hydrogenCount(e.to) : 0;
    if (count > 0 && (zz != zKey || hh != hKey)) return false;
    zKey = zz;
    hKey = hh;
    ++count;
  }
  return count >= 2 && zKey > 1;
}

TorsionCategory classify(const Graph& g, int bondIdx) {
  const Bond& b = g.m.bonds[bondIdx];
  if (b.order == BondOrder::Metal || isMetal(g.z(b.a)) || isMetal(g.z(b.b))) return TorsionCategory::Metal;
  if (b.aromatic) return TorsionCategory::Aromatic;
  if (g.classRing[bondIdx]) return TorsionCategory::Ring;
  if (b.order == BondOrder::Triple || collinearEnd(g, b.a, b.b) || collinearEnd(g, b.b, b.a))
    return TorsionCategory::Triple;
  if (b.order == BondOrder::Double) return TorsionCategory::Double;
  for (int k = 0; k < 2; ++k) {
    const int c = k == 0 ? b.a : b.b, y = k == 0 ? b.b : b.a;
    const int x = carbonylLike(g, c, y);
    if (x < 0 || (g.z(y) != 7 && g.z(y) != 8)) continue;
    return g.z(x) == 8 && g.z(y) == 7 ? TorsionCategory::Amide : TorsionCategory::Resonance;
  }
  if (symmetricEnd(g, b.a, b.b) || symmetricEnd(g, b.b, b.a)) return TorsionCategory::Symmetric;
  return TorsionCategory::Single;
}

// Lowest-index heavy neighbor of end other than other, preferring non-metals.
int dihedralNeighbor(const Graph& g, int end, int other) {
  int best = -1, bestMetal = -1;
  for (const Edge& e : g.adj[end]) {
    if (e.to == other || !g.heavy(e.to)) continue;
    int& slot = isMetal(g.z(e.to)) ? bestMetal : best;
    if (slot < 0 || e.to < slot) slot = e.to;
  }
  return best >= 0 ? best : bestMetal;
}

// Bonds of the given class that close a cycle. Bonds outside skip are added
// to a spanning forest in index order with the candidates last, so the
// molecule stays connected after the returned bonds are broken.
std::vector<int> cycleClosing(const Molecule& m, const std::vector<bool>& candidate, const std::vector<bool>& skip) {
  std::vector<int> parent(m.atoms.size());
  for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
  auto find = [&](int x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  std::vector<int> out;
  for (int pass = 0; pass < 2; ++pass)
    for (int bi = 0; bi < static_cast<int>(m.bonds.size()); ++bi) {
      if (skip[bi] || candidate[bi] != (pass == 1)) continue;
      const int a = find(m.bonds[bi].a), b = find(m.bonds[bi].b);
      if (a == b && pass == 1)
        out.push_back(bi);
      else
        parent[a] = b;
    }
  return out;
}

// Inter-residue bonds other than backbone links and metal bonds that close a
// cycle of the bond graph without metal bonds.
std::vector<int> findCrosslinkBreaks(const Graph& g, const Topology& t) {
  const Molecule& m = g.m;
  std::vector<bool> link(m.bonds.size(), false);
  for (const BackboneUnit& u : t.units) {
    if (u.next < 0) continue;
    const int from = t.backboneType == BackboneType::Protein ? u.atoms[2] : u.atoms[5];
    const int bi = g.bond(from, t.units[u.next].atoms[0]);
    if (bi >= 0) link[bi] = true;
  }
  std::vector<bool> candidate(m.bonds.size(), false), metal(m.bonds.size(), false);
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) {
    const Bond& b = m.bonds[bi];
    const int ra = m.atoms[b.a].residue, rb = m.atoms[b.b].residue;
    metal[bi] = isMetalBond(m, b);
    candidate[bi] = ra >= 0 && rb >= 0 && ra != rb && !link[bi] && !metal[bi];
  }
  return cycleClosing(m, candidate, metal);
}

// Metal bonds that close a cycle once crosslinks are broken.
std::vector<int> findMetalBreaks(const Molecule& m, const std::vector<int>& crosslinks) {
  std::vector<bool> metal(m.bonds.size(), false), skip(m.bonds.size(), false);
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) metal[bi] = isMetalBond(m, m.bonds[bi]);
  for (int bi : crosslinks) skip[bi] = true;
  return cycleClosing(m, metal, skip);
}

bool insideBox(const Vec3& x, const Box& box, double margin) {
  for (int k = 0; k < 3; ++k)
    if (std::abs(x[k] - box.center[k]) > box.size[k] / 2 + margin) return false;
  return true;
}

}  // namespace

Topology analyzeTopology(const Molecule& m, const TopologyParams& p) {
  Graph g(m);
  Topology t;
  auto prot = findProteinUnits(g);
  auto nuc = findNucleicUnits(g);
  if (!prot.empty() && !nuc.empty())
    logWarn(std::format("{}: both protein ({} units) and nucleic-acid ({} units) backbones found; using the "
                        "larger one",
                        m.id, prot.size(), nuc.size()));
  if (!prot.empty() && prot.size() >= nuc.size()) {
    t.backboneType = BackboneType::Protein;
    t.units = std::move(prot);
  } else if (!nuc.empty()) {
    t.backboneType = BackboneType::Nucleic;
    t.units = std::move(nuc);
  }

  t.crosslinkBreaks = findCrosslinkBreaks(g, t);
  t.metalBreaks = findMetalBreaks(m, t.crosslinkBreaks);
  if (p.brokenBonds) {
    t.brokenBonds = *p.brokenBonds;
  } else {
    t.brokenBonds = t.crosslinkBreaks;
    std::sort(t.brokenBonds.begin(), t.brokenBonds.end());
  }
  g.setBroken(t.brokenBonds);
  std::vector<bool> broken(m.bonds.size(), false);
  for (int b : t.brokenBonds) broken[b] = true;

  // Backbone angles keyed by bond index.
  std::map<int, Torsion> bb;
  auto addBackbone = [&](BackboneAngle a, int unit, std::array<int, 4> atoms) {
    const int bi = g.bond(atoms[1], atoms[2]);
    if (bi < 0) return;
    Torsion tor;
    tor.atoms = atoms;
    tor.bond = bi;
    tor.backbone = a;
    tor.unit = unit;
    tor.residue = t.units[unit].residue;
    bb[bi] = tor;
  };
  for (int i = 0; i < static_cast<int>(t.units.size()); ++i) {
    const auto& u = t.units[i].atoms;
    const int pv = t.units[i].prev, nx = t.units[i].next;
    if (t.backboneType == BackboneType::Protein) {
      if (pv >= 0) addBackbone(BackboneAngle::Phi, i, {t.units[pv].atoms[2], u[0], u[1], u[2]});
      if (nx >= 0) {
        const auto& v = t.units[nx].atoms;
        addBackbone(BackboneAngle::Psi, i, {u[0], u[1], u[2], v[0]});
        addBackbone(BackboneAngle::Omega, i, {u[1], u[2], v[0], v[1]});
      }
    } else {
      if (pv >= 0 && u[0] >= 0) addBackbone(BackboneAngle::Alpha, i, {t.units[pv].atoms[5], u[0], u[1], u[2]});
      if (u[0] >= 0) addBackbone(BackboneAngle::Beta, i, {u[0], u[1], u[2], u[3]});
      addBackbone(BackboneAngle::Gamma, i, {u[1], u[2], u[3], u[4]});
      addBackbone(BackboneAngle::Delta, i, {u[2], u[3], u[4], u[5]});
      if (nx >= 0) {
        const auto& v = t.units[nx].atoms;
        addBackbone(BackboneAngle::Epsilon, i, {u[3], u[4], u[5], v[0]});
        addBackbone(BackboneAngle::Zeta, i, {u[4], u[5], v[0], v[1]});
      }
    }
  }

  // StaM flexible range (spec 4.4).
  const bool restrict = m.role == Role::StaM && p.box.has_value();
  std::vector<bool> atomIn(m.atoms.size(), true), resIn(m.residues.size(), !restrict);
  if (restrict) {
    for (std::size_t i = 0; i < m.atoms.size(); ++i) {
      atomIn[i] = g.heavy(static_cast<int>(i)) && insideBox(m.atoms[i].xyz, *p.box, p.staFlexMargin);
      if (atomIn[i] && m.atoms[i].residue >= 0) resIn[m.atoms[i].residue] = true;
    }
  }
  auto inRange = [&](int a) { return m.atoms[a].residue >= 0 ? resIn[m.atoms[a].residue] : atomIn[a]; };

  for (int bi = 0; bi < static_cast<int>(m.bonds.size()); ++bi) {
    const Bond& b = m.bonds[bi];
    if (!g.heavy(b.a) || !g.heavy(b.b)) continue;
    Torsion tor;
    if (auto it = bb.find(bi); it != bb.end()) {
      tor = it->second;
    } else {
      const int a0 = dihedralNeighbor(g, b.a, b.b), d0 = dihedralNeighbor(g, b.b, b.a);
      if (a0 < 0 || d0 < 0) continue;
      tor.atoms = {a0, b.a, b.b, d0};
      tor.bond = bi;
      tor.residue = m.atoms[b.a].residue;
    }
    tor.category = classify(g, bi);
    const auto& x = tor.atoms;
    tor.angle = dihedralDeg(m.atoms[x[0]].xyz, m.atoms[x[1]].xyz, m.atoms[x[2]].xyz, m.atoms[x[3]].xyz);
    tor.inFlexRange = !restrict || inRange(b.a) || inRange(b.b);
    tor.broken = broken[bi];
    if (tor.backbone == BackboneAngle::Omega || tor.broken)
      tor.rotatable = false;
    else if (tor.backbone != BackboneAngle::None)
      tor.rotatable = tor.category == TorsionCategory::Single && m.role == Role::MobM;
    else
      tor.rotatable = tor.category == TorsionCategory::Single && tor.inFlexRange;
    t.torsions.push_back(tor);
  }
  return t;
}

std::vector<bool> ringBondMask(const Molecule& m, const std::vector<int>& broken) {
  Graph g(m);
  if (!broken.empty()) g.setBroken(broken);
  return g.classRing;
}

std::vector<int> nonPlanarOmegas(const Topology& t, double toleranceDeg) {
  std::vector<int> out;
  for (std::size_t i = 0; i < t.torsions.size(); ++i) {
    const Torsion& tor = t.torsions[i];
    if (tor.backbone != BackboneAngle::Omega) continue;
    const double d0 = std::abs(tor.angle), d180 = 180.0 - d0;
    if (d0 > toleranceDeg && d180 > toleranceDeg) out.push_back(static_cast<int>(i));
  }
  return out;
}

double dihedralDeg(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  auto sub = [](const Vec3& u, const Vec3& v) { return Vec3{u[0] - v[0], u[1] - v[1], u[2] - v[2]}; };
  auto cross = [](const Vec3& u, const Vec3& v) {
    return Vec3{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  };
  auto dot = [](const Vec3& u, const Vec3& v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
  const Vec3 b1 = sub(b, a), b2 = sub(c, b), b3 = sub(d, c);
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  double deg = std::atan2(std::sqrt(dot(b2, b2)) * dot(b1, n2), dot(n1, n2)) * 180.0 / std::numbers::pi;
  if (deg <= -180.0) deg += 360.0;
  return deg;
}

const char* backboneTypeName(BackboneType t) {
  switch (t) {
    case BackboneType::Protein: return "protein";
    case BackboneType::Nucleic: return "nucleic";
    default: return "none";
  }
}

const char* torsionCategoryName(TorsionCategory c) {
  switch (c) {
    case TorsionCategory::Single: return "single";
    case TorsionCategory::Amide: return "amide";
    case TorsionCategory::Resonance: return "resonance";
    case TorsionCategory::Double: return "double";
    case TorsionCategory::Triple: return "triple";
    case TorsionCategory::Aromatic: return "aromatic";
    case TorsionCategory::Ring: return "ring";
    case TorsionCategory::Symmetric: return "symmetric";
    case TorsionCategory::Metal: return "metal";
  }
  return "";
}

const char* backboneAngleName(BackboneAngle a) {
  switch (a) {
    case BackboneAngle::Phi: return "phi";
    case BackboneAngle::Psi: return "psi";
    case BackboneAngle::Omega: return "omega";
    case BackboneAngle::Alpha: return "alpha";
    case BackboneAngle::Beta: return "beta";
    case BackboneAngle::Gamma: return "gamma";
    case BackboneAngle::Delta: return "delta";
    case BackboneAngle::Epsilon: return "epsilon";
    case BackboneAngle::Zeta: return "zeta";
    default: return "";
  }
}

}  // namespace od
