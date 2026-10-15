// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "torsion_tree.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace od {
namespace {

struct Edge {
  int to;
  int bond;
};

std::vector<std::vector<Edge>> adjacency(const Molecule& m, const std::vector<int>& broken) {
  std::vector<bool> skip(m.bonds.size(), false);
  for (int b : broken) skip[b] = true;
  std::vector<std::vector<Edge>> adj(m.atoms.size());
  for (std::size_t i = 0; i < m.bonds.size(); ++i) {
    if (skip[i]) continue;
    adj[m.bonds[i].a].push_back({m.bonds[i].b, static_cast<int>(i)});
    adj[m.bonds[i].b].push_back({m.bonds[i].a, static_cast<int>(i)});
  }
  return adj;
}

// Atoms reached from start without crossing bond; stops early at an anchor
// when stopAtAnchor. Returns whether an anchor was reached.
bool side(const std::vector<std::vector<Edge>>& adj, int start, int bond, const std::vector<bool>& anchor,
          bool stopAtAnchor, std::vector<int>& out) {
  out.clear();
  std::vector<bool> seen(adj.size(), false);
  seen[start] = true;
  out.push_back(start);
  bool found = anchor[start];
  for (std::size_t i = 0; i < out.size() && !(found && stopAtAnchor); ++i)
    for (const Edge& e : adj[out[i]]) {
      if (e.bond == bond || seen[e.to]) continue;
      seen[e.to] = true;
      out.push_back(e.to);
      if (anchor[e.to]) found = true;
    }
  return found;
}

bool heavy(const Molecule& m, int a) { return m.atoms[a].atomicNum > 1; }

Vec3 heavyCentroid(const Molecule& m) {
  Vec3 c{};
  int n = 0;
  for (const Atom& a : m.atoms)
    if (a.atomicNum > 1) {
      for (int k = 0; k < 3; ++k) c[k] += a.xyz[k];
      ++n;
    }
  if (n > 0)
    for (double& x : c) x /= n;
  return c;
}

double dist2(const Vec3& a, const Vec3& b) {
  double s = 0;
  for (int k = 0; k < 3; ++k) s += (a[k] - b[k]) * (a[k] - b[k]);
  return s;
}

// Chains as unit lists in chain order.
std::vector<std::vector<int>> chains(const Topology& t) {
  std::vector<std::vector<int>> out;
  std::vector<bool> seen(t.units.size(), false);
  auto walk = [&](int s) {
    out.emplace_back();
    for (int u = s; u >= 0 && !seen[u]; u = t.units[u].next) {
      seen[u] = true;
      out.back().push_back(u);
    }
  };
  for (std::size_t i = 0; i < t.units.size(); ++i)
    if (t.units[i].prev < 0) walk(static_cast<int>(i));
  for (std::size_t i = 0; i < t.units.size(); ++i)
    if (!seen[i]) walk(static_cast<int>(i));
  return out;
}

// Rotates p around the axis through o with unit direction k by angle (radians).
Vec3 rotate(const Vec3& p, const Vec3& o, const Vec3& k, double c, double s) {
  const Vec3 v{p[0] - o[0], p[1] - o[1], p[2] - o[2]};
  const double dot = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
  const Vec3 cr{k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]};
  Vec3 r;
  for (int i = 0; i < 3; ++i) r[i] = o[i] + v[i] * c + cr[i] * s + k[i] * dot * (1 - c);
  return r;
}

}  // namespace

std::vector<bool> fixedUnits(const Topology& t, const std::vector<TorsionMode>& mode) {
  std::vector<bool> fixed(t.units.size(), true);
  for (std::size_t i = 0; i < t.torsions.size(); ++i)
    if (t.torsions[i].backbone != BackboneAngle::None && mode[i] != TorsionMode::Fixed) fixed[t.torsions[i].unit] = false;
  return fixed;
}

std::vector<int> anchorAtoms(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                             bool stamFixedResMovable) {
  if (m.role == Role::StaM && !stamFixedResMovable && t.backboneType != BackboneType::None) {
    const auto fixed = fixedUnits(t, mode);
    std::vector<int> out;
    for (std::size_t u = 0; u < t.units.size(); ++u)
      if (fixed[u])
        for (int a : t.units[u].atoms)
          if (a >= 0) out.push_back(a);
    if (!out.empty()) {
      std::sort(out.begin(), out.end());
      return out;
    }
  }
  const Vec3 c = heavyCentroid(m);
  int best = -1;
  double bestD = std::numeric_limits<double>::max();
  auto consider = [&](int a) {
    if (a < 0) return;
    if (const double d = dist2(m.atoms[a].xyz, c); d < bestD) {
      bestD = d;
      best = a;
    }
  };
  if (m.role == Role::StaM && t.backboneType != BackboneType::None) {
    for (const auto& u : t.units) consider(u.atoms[t.backboneType == BackboneType::Protein ? 1 : 3]);
  } else {
    for (int a = 0; a < static_cast<int>(m.atoms.size()); ++a)
      if (heavy(m, a)) consider(a);
  }
  if (best < 0 && !m.atoms.empty()) best = 0;
  return best < 0 ? std::vector<int>{} : std::vector<int>{best};
}

int backboneLink(const Molecule& m, const Topology& t, int unit, int next) {
  if (unit < 0 || next < 0) return -1;
  const bool protein = t.backboneType == BackboneType::Protein;
  const int a = t.units[unit].atoms[protein ? 2 : 5], b = t.units[next].atoms[0];
  if (a < 0 || b < 0) return -1;
  for (std::size_t i = 0; i < m.bonds.size(); ++i)
    if ((m.bonds[i].a == a && m.bonds[i].b == b) || (m.bonds[i].a == b && m.bonds[i].b == a)) return static_cast<int>(i);
  return -1;
}

std::vector<Loop> findLoops(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                            bool stamFixedResMovable) {
  std::vector<Loop> loops;
  if (m.role != Role::StaM || stamFixedResMovable || t.backboneType == BackboneType::None) return loops;
  const auto fixed = fixedUnits(t, mode);
  std::vector<bool> broken(m.bonds.size(), false);
  for (int b : t.brokenBonds) broken[b] = true;
  for (const auto& ch : chains(t)) {
    for (std::size_t i = 0; i < ch.size();) {
      if (fixed[ch[i]]) {
        ++i;
        continue;
      }
      Loop l;
      if (i > 0) l.before = ch[i - 1];
      while (i < ch.size() && !fixed[ch[i]]) l.units.push_back(ch[i++]);
      if (i < ch.size()) l.after = ch[i];
      if (l.before >= 0 && l.after >= 0) {
        std::vector<int> seq{l.before};
        seq.insert(seq.end(), l.units.begin(), l.units.end());
        seq.push_back(l.after);
        for (std::size_t k = 0; k + 1 < seq.size() && l.breakBond < 0; ++k)
          if (const int bi = backboneLink(m, t, seq[k], seq[k + 1]); bi >= 0 && broken[bi]) l.breakBond = bi;
      }
      loops.push_back(std::move(l));
    }
  }
  return loops;
}

TorsionTree buildTorsionTree(const Molecule& m, const Topology& t, const std::vector<TorsionMode>& mode,
                             bool stamFixedResMovable, std::vector<TreeProblem>& problems) {
  TorsionTree tree;
  const int n = static_cast<int>(m.atoms.size());
  tree.atomCount = n;
  tree.anchors = anchorAtoms(m, t, mode, stamFixedResMovable);
  tree.loops = findLoops(m, t, mode, stamFixedResMovable);
  std::vector<bool> anchor(n, false);
  for (int a : tree.anchors) anchor[a] = true;
  const auto adj = adjacency(m, t.brokenBonds);

  // Every part left by the broken bonds needs an anchor.
  {
    std::vector<bool> seen(n, false);
    std::vector<int> part;
    for (int s = 0; s < n; ++s) {
      if (seen[s]) continue;
      if (!side(adj, s, -1, anchor, false, part)) problems.push_back({TreeProblem::Kind::PartWithoutAnchor, -1, s});
      for (int a : part) seen[a] = true;
    }
  }

  // Closure: virtual copies for broken bonds, distance restraints for metal bonds (spec 4.8).
  std::vector<int> torsionOfBond(m.bonds.size(), -1);
  for (std::size_t i = 0; i < t.torsions.size(); ++i) torsionOfBond[t.torsions[i].bond] = static_cast<int>(i);
  for (int bi : t.brokenBonds) {
    const Bond& b = m.bonds[bi];
    if (isMetalBond(m, b)) {
      tree.restraints.push_back({b.a, b.b, std::sqrt(dist2(m.atoms[b.a].xyz, m.atoms[b.b].xyz))});
      continue;
    }
    int A = b.a, B = b.b, Y = -1;
    if (const int ti = torsionOfBond[bi]; ti >= 0) {
      const auto& x = t.torsions[ti].atoms;
      A = x[1];
      B = x[2];
      Y = x[3];
    } else {
      for (const Bond& o : m.bonds) {
        const int other = o.a == B ? o.b : o.b == B ? o.a : -1;
        if (other >= 0 && other != A && heavy(m, other) && (Y < 0 || other < Y)) Y = other;
      }
    }
    tree.virtuals.push_back({A, B, m.atoms[B].xyz});
    if (Y >= 0) tree.virtuals.push_back({A, Y, m.atoms[Y].xyz});
    tree.virtuals.push_back({B, A, m.atoms[A].xyz});
  }

  // Moving sides.
  std::vector<int> sideB, sideA;
  for (std::size_t ti = 0; ti < t.torsions.size(); ++ti) {
    if (mode[ti] == TorsionMode::Fixed) continue;
    const Torsion& tor = t.torsions[ti];
    const int a = tor.atoms[1], b = tor.atoms[2];
    TreeTorsion tt;
    tt.torsion = static_cast<int>(ti);
    tt.inputAngle = tor.angle;
    if (!side(adj, b, tor.bond, anchor, true, sideB)) {
      tt.fixedAtom = a;
      tt.movingAtom = b;
      side(adj, b, tor.bond, anchor, false, tt.moving);
    } else if (!side(adj, a, tor.bond, anchor, true, sideA)) {
      tt.fixedAtom = b;
      tt.movingAtom = a;
      tt.moving = sideA;
    } else {
      problems.push_back({TreeProblem::Kind::BothSidesFixed, static_cast<int>(ti), -1});
      continue;
    }
    std::sort(tt.moving.begin(), tt.moving.end());
    tree.torsions.push_back(std::move(tt));
  }
  // Virtual atoms move with their parents.
  {
    std::vector<bool> in(n, false);
    for (auto& tt : tree.torsions) {
      for (int a : tt.moving) in[a] = true;
      for (std::size_t v = 0; v < tree.virtuals.size(); ++v)
        if (in[tree.virtuals[v].parent]) tt.moving.push_back(n + static_cast<int>(v));
      for (int a : tt.moving)
        if (a < n) in[a] = false;
    }
  }
  // Parents before children: moving sets are nested or disjoint.
  std::stable_sort(tree.torsions.begin(), tree.torsions.end(),
                   [](const TreeTorsion& x, const TreeTorsion& y) { return x.moving.size() > y.moving.size(); });
  std::vector<int> owner(n, -1);
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    auto& tt = tree.torsions[i];
    tt.parent = owner[tt.movingAtom];
    for (int a : tt.moving)
      if (a < n) owner[a] = static_cast<int>(i);
  }
  return tree;
}

std::vector<Vec3> buildCoordinates(const Molecule& m, const TorsionTree& tree, std::span<const double> anglesDeg) {
  std::vector<Vec3> xyz;
  xyz.reserve(m.atoms.size() + tree.virtuals.size());
  for (const Atom& a : m.atoms) xyz.push_back(a.xyz);
  for (const VirtualAtom& v : tree.virtuals) xyz.push_back(v.xyz);
  // Children first; the axis atoms of a torsion are not moved by its descendants.
  for (std::size_t i = tree.torsions.size(); i-- > 0;) {
    const TreeTorsion& tt = tree.torsions[i];
    const Vec3 o = xyz[tt.fixedAtom], p = xyz[tt.movingAtom];
    Vec3 k{p[0] - o[0], p[1] - o[1], p[2] - o[2]};
    const double len = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
    if (len == 0) continue;
    for (double& x : k) x /= len;
    const double delta = (anglesDeg[i] - tt.inputAngle) * std::numbers::pi / 180.0;
    const double c = std::cos(delta), s = std::sin(delta);
    for (int a : tt.moving) xyz[a] = rotate(xyz[a], o, k, c, s);
  }
  return xyz;
}

ClosureState evaluateClosure(const TorsionTree& tree, std::span<const Vec3> xyz, double G) {
  ClosureState st;
  for (std::size_t v = 0; v < tree.virtuals.size(); ++v) {
    const double d2 = dist2(xyz[tree.atomCount + v], xyz[tree.virtuals[v].copyOf]);
    st.energy += G * d2;
    st.maxDeviation = std::max(st.maxDeviation, std::sqrt(d2));
  }
  for (const auto& r : tree.restraints) {
    const double dd = std::abs(std::sqrt(dist2(xyz[r.a], xyz[r.b])) - r.r0);
    st.energy += G * dd * dd;
    st.maxDeviation = std::max(st.maxDeviation, dd);
  }
  return st;
}

}  // namespace od
