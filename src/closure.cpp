// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "closure.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <queue>

#include "topology.hpp"

namespace od {

namespace {

constexpr double kDeg = 180.0 / std::numbers::pi;
// Torsions needed to close the three atom pairs of a broken bond exactly (spec 4.8).
constexpr int kFullDof = 6;

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 unit(const Vec3& a) {
  const double l = norm(a);
  return l > 0 ? (1.0 / l) * a : Vec3{};
}
// Component of v perpendicular to the unit vector k.
Vec3 perp(const Vec3& v, const Vec3& k) { return v - dot(v, k) * k; }
Vec3 rotateAbout(const Vec3& p, const Vec3& o, const Vec3& k, double c, double s) {
  const Vec3 v = p - o;
  return o + c * v + s * cross(k, v) + (dot(k, v) * (1 - c)) * k;
}
double angleCos(const Vec3& a, const Vec3& pivot, const Vec3& b) { return dot(unit(a - pivot), unit(b - pivot)); }

// Orthonormal frame (columns) at origin o from three points: e1 along p - o, e2 in the plane of q.
struct Frame {
  Vec3 o, e1, e2, e3;
};
Frame frameOf(const Vec3& o, const Vec3& p, const Vec3& q) {
  Frame f;
  f.o = o;
  f.e1 = unit(p - o);
  f.e2 = unit(perp(q - o, f.e1));
  f.e3 = cross(f.e1, f.e2);
  return f;
}
// Rigid transform taking frame `from` onto frame `to`.
Vec3 mapFrame(const Frame& from, const Frame& to, const Vec3& x) {
  const Vec3 v = x - from.o;
  return to.o + dot(v, from.e1) * to.e1 + dot(v, from.e2) * to.e2 + dot(v, from.e3) * to.e3;
}

bool contains(const std::vector<int>& sorted, int v) { return std::binary_search(sorted.begin(), sorted.end(), v); }

}  // namespace

LoopCloser::LoopCloser(const EnergySystem& sys, const std::vector<EnergyPart>& parts, const std::vector<bool>& freeMask,
                       const std::vector<std::vector<int>>& segments, const ClosureParams& p)
    : sys_(sys), parts_(parts), p_(p) {
  for (std::size_t k = 0; k < parts_.size(); ++k)
    for (std::size_t t = 0; t < parts_[k].tree->torsions.size(); ++t) {
      partOfTorsion_.push_back(static_cast<int>(k));
      localTorsion_.push_back(static_cast<int>(t));
    }
  if (p_.method == ClosureMethod::None) return;
  path_.assign(partOfTorsion_.size(), false);
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    const Molecule& m = *parts_[k].mol;
    const TorsionTree& tree = *parts_[k].tree;
    const int off = sys_.torsionOffset(static_cast<int>(k));
    const int n = tree.atomCount;
    auto adjustable = [&](int local) { return static_cast<bool>(freeMask[off + local]); };
    // Graph without the broken bonds; node n joins the anchor atoms, which never move.
    std::vector<std::vector<int>> adj(n + 1);
    if (tree.anchors.size() > 1)
      for (int x : tree.anchors) {
        adj[n].push_back(x);
        adj[x].push_back(n);
      }
    for (const Bond& b : m.bonds) {
      bool broken = false;
      for (const auto& l : tree.links) broken |= (l.a == b.a && l.b == b.b) || (l.a == b.b && l.b == b.a);
      if (broken) continue;
      adj[b.a].push_back(b.b);
      adj[b.b].push_back(b.a);
    }
    // Each link alone: its pairs or restraint, ring path and adjustable torsions.
    std::vector<Link> single;
    for (const ClosureLink& cl : tree.links) {
      Link l;
      l.part = static_cast<int>(k);
      l.members.push_back(static_cast<int>(single.size()));
      if (cl.restraint >= 0) {
        const auto& r = tree.restraints[cl.restraint];
        l.restraints.push_back({r.a, r.b, r.r0});
        l.pathAtoms = {r.a, r.b};
        for (std::size_t t = 0; t < tree.torsions.size(); ++t) {
          const auto& tt = tree.torsions[t];
          const auto& mv = tt.moving;
          // An atom on the axis does not move.
          auto moves = [&](int x) { return contains(mv, x) && x != tt.fixedAtom && x != tt.movingAtom; };
          auto still = [&](int x) { return !contains(mv, x) || x == tt.fixedAtom || x == tt.movingAtom; };
          if (!((moves(r.a) && still(r.b)) || (moves(r.b) && still(r.a)))) continue;
          path_[off + t] = true;
          if (adjustable(static_cast<int>(t))) l.dof.push_back(off + static_cast<int>(t));
        }
      } else {
        l.pairs.push_back({n + cl.copyB, cl.b});
        if (cl.copyY >= 0) l.pairs.push_back({n + cl.copyY, cl.y});
        l.pairs.push_back({n + cl.copyA, cl.a});
        l.pairEnd.push_back(static_cast<int>(l.pairs.size()));
        // Ring path a ... b without the broken bond.
        std::vector<int> prev(n + 1, -2);
        std::queue<int> q;
        prev[cl.a] = -1;
        q.push(cl.a);
        while (!q.empty() && prev[cl.b] == -2) {
          const int u = q.front();
          q.pop();
          for (int v : adj[u])
            if (prev[v] == -2) {
              prev[v] = u;
              q.push(v);
            }
        }
        if (prev[cl.b] != -2) {
          for (int v = cl.b; v >= 0; v = prev[v])
            if (v < n) l.ring.push_back(v);
          std::reverse(l.ring.begin(), l.ring.end());
        }
        for (std::size_t j = 0; j + 1 < l.ring.size(); ++j) {
          int g = -1;
          for (std::size_t t = 0; t < tree.torsions.size(); ++t) {
            const auto& tt = tree.torsions[t];
            if ((tt.fixedAtom == l.ring[j] && tt.movingAtom == l.ring[j + 1]) ||
                (tt.fixedAtom == l.ring[j + 1] && tt.movingAtom == l.ring[j])) {
              path_[off + t] = true;
              if (adjustable(static_cast<int>(t))) g = off + static_cast<int>(t);
            }
          }
          l.bondTorsion.push_back(g);
          if (g >= 0) l.dof.push_back(g);
        }
        for (std::size_t i = 1; i + 1 < l.ring.size(); ++i)
          if (l.bondTorsion[i - 1] >= 0 && l.bondTorsion[i] >= 0) l.pivots.push_back(static_cast<int>(i));
        // KIC needs Y' and three pivots i1 < i2 < i3 at least two bonds apart.
        int found = 0, last = -10;
        for (int i : l.pivots)
          if (i >= last + 2) {
            ++found;
            last = i;
          }
        l.kic = cl.copyY >= 0 && found >= 3;
        l.pathAtoms = l.ring;
      }
      single.push_back(std::move(l));
    }
    // Groups: links sharing an adjustable torsion (union-find), in order of their first link.
    const int L = static_cast<int>(single.size());
    std::vector<int> root(L);
    for (int i = 0; i < L; ++i) root[i] = i;
    auto find = [&](int i) {
      while (root[i] != i) i = root[i] = root[root[i]];
      return i;
    };
    for (int i = 0; i < L; ++i)
      for (int j = i + 1; j < L; ++j) {
        const bool shared = std::any_of(single[i].dof.begin(), single[i].dof.end(), [&](int g) {
          return std::find(single[j].dof.begin(), single[j].dof.end(), g) != single[j].dof.end();
        });
        if (shared) root[find(j)] = find(i);
      }
    std::vector<char> seen(L, 0);
    for (int i = 0; i < L; ++i) {
      if (seen[i]) continue;
      std::vector<int> mem;
      for (int j = i; j < L; ++j)
        if (find(j) == find(i)) {
          mem.push_back(j);
          seen[j] = 1;
        }
      Link l;
      if (mem.size() == 1) {
        l = std::move(single[i]);
      } else {
        l.part = static_cast<int>(k);
        for (int j : mem) {
          Link& s = single[j];
          l.members.push_back(j);
          if (!s.pairs.empty()) {
            l.pairs.insert(l.pairs.end(), s.pairs.begin(), s.pairs.end());
            l.pairEnd.push_back(static_cast<int>(l.pairs.size()));
          }
          l.restraints.insert(l.restraints.end(), s.restraints.begin(), s.restraints.end());
          l.pathAtoms.insert(l.pathAtoms.end(), s.pathAtoms.begin(), s.pathAtoms.end());
          for (int g : s.dof)
            if (std::find(l.dof.begin(), l.dof.end(), g) == l.dof.end()) l.dof.push_back(g);
        }
      }
      const int copies = static_cast<int>(l.pairEnd.size());
      l.partial = copies == 0 || static_cast<int>(l.dof.size()) < kFullDof * copies;
      std::vector<int> sortedDof = l.dof;
      std::sort(sortedDof.begin(), sortedDof.end());
      for (std::size_t s = 0; s < segments.size(); ++s)
        if (!segments[s].empty() && std::all_of(segments[s].begin(), segments[s].end(),
                                                [&](int g) { return contains(sortedDof, g); }))
        {
          l.resampleSegments.push_back(static_cast<int>(s));
          l.resampleCovered.insert(l.resampleCovered.end(), segments[s].begin(), segments[s].end());
        }
      std::sort(l.resampleCovered.begin(), l.resampleCovered.end());
      // RCD tracking: closure points and the axis atoms of the adjustable torsions.
      for (const auto& pr : l.pairs) {
        l.points.push_back(pr.moving);
        l.points.push_back(pr.target);
      }
      for (const auto& r : l.restraints) {
        l.points.push_back(r.a);
        l.points.push_back(r.b);
      }
      for (int g : l.dof) {
        const auto& tt = tree.torsions[localTorsion_[g]];
        l.points.push_back(tt.fixedAtom);
        l.points.push_back(tt.movingAtom);
      }
      std::vector<int> pts = l.points;
      std::sort(pts.begin(), pts.end());
      pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
      auto at = [&](int c) { return static_cast<int>(std::lower_bound(pts.begin(), pts.end(), c) - pts.begin()); };
      for (const auto& pr : l.pairs) {
        l.pairPoint.push_back(at(pr.moving));
        l.pairPoint.push_back(at(pr.target));
      }
      for (const auto& r : l.restraints) {
        l.restraintPoint.push_back(at(r.a));
        l.restraintPoint.push_back(at(r.b));
      }
      for (int g : l.dof) {
        const auto& tt = tree.torsions[localTorsion_[g]];
        l.axisFixed.push_back(at(tt.fixedAtom));
        l.axisMoving.push_back(at(tt.movingAtom));
        std::vector<int> mv;
        for (std::size_t i = 0; i < pts.size(); ++i)
          if (contains(tt.moving, pts[i])) mv.push_back(static_cast<int>(i));
        l.moves.push_back(std::move(mv));
      }
      l.points = std::move(pts);
      links_.push_back(std::move(l));
    }
  }
  // Clash pairs (heavy atoms). Default: pairs with a ring path atom (metal
  // distance: its two atoms) moved by the adjustable torsions. growClose:
  // pairs with any atom moved by them, neither atom built after the closure.
  const int nc = sys_.coordOffset(static_cast<int>(parts_.size()) - 1) + parts_.back().tree->atomCount +
                 static_cast<int>(parts_.back().tree->virtuals.size());
  std::vector<char> later(nc, 0);
  if (p_.growClose)
    for (std::size_t k = 0; k < parts_.size(); ++k) {
      // Deepest torsion moving each atom; the moving bond atom goes with the parent (as GrowBuilder).
      const auto& tt = parts_[k].tree->torsions;
      const int off = sys_.torsionOffset(static_cast<int>(k)), co = sys_.coordOffset(static_cast<int>(k));
      std::vector<int> owner(parts_[k].tree->atomCount, -1);
      for (std::size_t i = 0; i < tt.size(); ++i) {
        for (int c : tt[i].moving)
          if (c < static_cast<int>(owner.size())) owner[c] = static_cast<int>(i);
        owner[tt[i].movingAtom] = tt[i].parent;
      }
      for (std::size_t c = 0; c < owner.size(); ++c)
        if (owner[c] >= 0 && !path_[off + owner[c]]) later[co + c] = 1;
    }
  for (Link& l : links_) {
    const TorsionTree& tree = *parts_[l.part].tree;
    const int n = tree.atomCount, co = sys_.coordOffset(l.part);
    std::vector<char> moved(n, 0), checked(nc, 0);
    for (int g : l.dof)
      for (int c : tree.torsions[localTorsion_[g]].moving)
        if (c < n) moved[c] = 1;
    if (p_.growClose) {
      for (int c = 0; c < n; ++c)
        if (moved[c] && !later[co + c]) checked[co + c] = 1;
    } else {
      for (int c : l.pathAtoms)
        if (moved[c]) checked[co + c] = 1;
    }
    for (const auto& pr : sys_.pairs())
      if (pr.heavy && (checked[pr.i] || checked[pr.j]) && !later[pr.i] && !later[pr.j])
        l.clashPairs.push_back({pr.i, pr.j});
    for (int g : l.dof) {
      std::vector<int> mv;
      for (int c : tree.torsions[localTorsion_[g]].moving) mv.push_back(co + c);
      std::sort(mv.begin(), mv.end());
      l.movesSys.push_back(std::move(mv));
    }
  }
}

bool LoopCloser::active() const {
  return std::any_of(links_.begin(), links_.end(), [](const Link& l) { return !l.dof.empty(); });
}

std::vector<std::string> LoopCloser::describe() const {
  std::vector<std::string> out;
  for (const Link& l : links_) {
    const Molecule& m = *parts_[l.part].mol;
    const TorsionTree& tree = *parts_[l.part].tree;
    auto name = [&](int a) {
      const Atom& x = m.atoms[a];
      if (x.residue < 0) return std::format("#{}", x.fileSerial);
      const ResidueInfo& r = m.residues[x.residue].info;
      return std::format("{}:{}{}:{}", r.chain, r.seq, r.iCode, x.name);
    };
    std::string label;
    for (int j : l.members) {
      const ClosureLink& cl = tree.links[j];
      if (!label.empty()) label += " + ";
      if (cl.restraint >= 0)
        label += std::format("{} - {}", name(tree.restraints[cl.restraint].a), name(tree.restraints[cl.restraint].b));
      else
        label += std::format("{} - {}", name(cl.a), name(cl.b));
    }
    const int copies = static_cast<int>(l.pairEnd.size());
    const bool kicAsked = p_.method == ClosureMethod::Kic;
    std::string method;
    if (l.dof.empty()) method = "none (no adjustable torsion; constraint only)";
    else if (copies == 0) method = "rcd, best effort (metal distance)";
    else if (l.partial) method = std::format("rcd, best effort (fewer than {} torsions)", kFullDof * copies);
    else if (l.members.size() > 1) method = kicAsked ? "rcd (kic not applicable to coupled links)" : "rcd";
    else if (kicAsked && l.kic) method = "kic";
    else if (kicAsked) method = "rcd (kic not applicable)";
    else method = "rcd";
    if (l.members.size() > 1)
      out.push_back(std::format("{}: coupled, {} links, {} torsions, {}", label, l.members.size(), l.dof.size(), method));
    else if (copies == 0)
      out.push_back(std::format("{} (metal distance): {} torsions, {}", label, l.dof.size(), method));
    else
      out.push_back(std::format("{}: ring path {} atoms, {} torsions, {} KIC pivots, {}", label, l.ring.size(),
                                l.dof.size(), l.pivots.size(), method));
  }
  return out;
}

std::vector<Vec3> LoopCloser::partCoords(int part, std::span<const double> a) const {
  const auto& tree = *parts_[part].tree;
  return buildCoordinates(*parts_[part].mol, tree, a.subspan(sys_.torsionOffset(part), tree.torsions.size()));
}

double LoopCloser::deviation(const Link& l, const std::vector<Vec3>& x) const {
  // Copy links: the largest RMSD of the pairs of one link; metal distances only: the largest |d - r0|.
  double d = 0;
  if (l.pairEnd.empty()) {
    for (const auto& r : l.restraints) d = std::max(d, std::abs(norm(x[r.a] - x[r.b]) - r.r0));
    return d;
  }
  for (std::size_t c = 0, b = 0; c < l.pairEnd.size(); b = l.pairEnd[c++]) {
    double s = 0;
    for (auto q = b; q < static_cast<std::size_t>(l.pairEnd[c]); ++q)
      s += dot(x[l.pairs[q].moving] - x[l.pairs[q].target], x[l.pairs[q].moving] - x[l.pairs[q].target]);
    d = std::max(d, std::sqrt(s / static_cast<double>(l.pairEnd[c] - b)));
  }
  return d;
}

bool LoopCloser::clashFree(const Link& l, std::span<const double> a) const {
  const auto x = sys_.coordinates(a);
  const double c2 = p_.clashDistance * p_.clashDistance;
  for (const auto& [i, j] : l.clashPairs)
    if (dot(x[i] - x[j], x[i] - x[j]) < c2) return false;
  return true;
}

void LoopCloser::resampleLink(Rng& rng, const Link& l, std::span<double> a, const SegmentSampler& resample) const {
  std::vector<int> covered;
  if (resample)
    for (int s : l.resampleSegments) resample(rng, s, a);
  if (resample) covered = l.resampleCovered;
  // Adjustable torsions outside every resampled segment: uniform.
  for (int g : l.dof)
    if (!contains(covered, g)) a[g] = rng.uniform(-180.0, 180.0);
}

bool LoopCloser::rcd(Rng& rng, const Link& l, std::span<double> a, const SegmentSampler& resample) const {
  // Best effort: one attempt; the deviation never grows, the result is kept.
  const int restarts = l.partial ? 0 : p_.rcdMaxRestart;
  // A group of coupled links converges more slowly: rcdMaxIter per link.
  const int maxIter = p_.rcdMaxIter * static_cast<int>(l.members.size());
  const std::vector<double> orig(a.begin(), a.end());
  const int nd = static_cast<int>(l.dof.size());
  const double tol = p_.rcdTolerance;
  std::vector<Vec3> P(l.points.size());
  std::vector<char> mv(l.points.size(), 0);
  auto dev = [&]() {
    double d = 0;
    if (l.pairEnd.empty()) {
      for (std::size_t q = 0; q < l.restraints.size(); ++q)
        d = std::max(d, std::abs(norm(P[l.restraintPoint[2 * q]] - P[l.restraintPoint[2 * q + 1]]) - l.restraints[q].r0));
      return d;
    }
    for (std::size_t c = 0, b = 0; c < l.pairEnd.size(); b = l.pairEnd[c++]) {
      double s = 0;
      for (auto q = b; q < static_cast<std::size_t>(l.pairEnd[c]); ++q) {
        const Vec3 v = P[l.pairPoint[2 * q]] - P[l.pairPoint[2 * q + 1]];
        s += dot(v, v);
      }
      d = std::max(d, std::sqrt(s / static_cast<double>(l.pairEnd[c] - b)));
    }
    return d;
  };
  // Clash term: sum over clash pairs of w (c - d)^2 for d < c, on a neighbour list of the pairs
  // closer than c + 4 Å, rebuilt every 20 updates; the final clash check is exact.
  const double cd = p_.clashDistance, w = p_.rcdClashWeight;
  const bool clashTerm = p_.rcdClash != RcdClash::None && cd > 0 && w > 0 && !l.partial && !l.clashPairs.empty();
  std::vector<Vec3> X;
  std::vector<std::pair<int, int>> nl;
  int nlAge = 0;
  auto buildNl = [&]() {
    nl.clear();
    const double r2 = (cd + 4) * (cd + 4);
    for (const auto& [i, j] : l.clashPairs)
      if (dot(X[i] - X[j], X[i] - X[j]) < r2) nl.push_back({i, j});
    nlAge = 0;
  };
  auto clashes = [&]() {
    const double c2 = cd * cd;
    for (const auto& [i, j] : l.clashPairs)
      if (dot(X[i] - X[j], X[i] - X[j]) < c2) return true;
    return false;
  };
  const double stallTol = p_.rcdStallTolerance;
  std::vector<char> tried(nd, 0);
  for (int att = 0; att <= restarts; ++att) {
    if (att > 0) {
      std::copy(orig.begin(), orig.end(), a.begin());
      resampleLink(rng, l, a, resample);
    }
    const auto x = partCoords(l.part, a);
    for (std::size_t i = 0; i < l.points.size(); ++i) P[i] = x[l.points[i]];
    bool useClash = clashTerm && p_.rcdClash == RcdClash::Full;
    if (clashTerm) {
      X = sys_.coordinates(std::span<const double>(a.data(), a.size()));
      buildNl();
    }
    // Converged: closed, and with the clash term also clash free.
    auto converged = [&]() {
      if (dev() > tol) return false;
      if (!clashTerm) return true;
      if (!clashes()) return true;
      if (!useClash) {
        useClash = true;
        buildNl();
      }
      return false;
    };
    bool conv = converged();
    int prev = -1;
    // Stall: every torsion updated by at most rcdStallTolerance since the last larger update.
    std::fill(tried.begin(), tried.end(), 0);
    int triedNum = 0;
    for (int it = 0; it < maxIter && !conv; ++it) {
      int d = static_cast<int>(rng.below(nd));
      if (nd > 1 && d == prev) d = (d + 1 + static_cast<int>(rng.below(nd - 1))) % nd;
      prev = d;
      for (int i : l.moves[d]) mv[i] = 1;
      const Vec3 o = P[l.axisFixed[d]];
      const Vec3 k = unit(P[l.axisMoving[d]] - o);
      // Copy pairs: sum |m - s|^2 = const - 2 (den cos d + num sin d). Metal distance:
      // |m - s|^2 = A - 2 B cos(d - d0).
      double num = 0, den = 0;
      bool copyTerm = false;
      struct RTerm {
        double A, B, d0, r0;
      };
      std::vector<RTerm> rt;
      for (std::size_t q = 0; q < l.pairPoint.size(); q += 2) {
        const int i1 = l.pairPoint[q], i2 = l.pairPoint[q + 1];
        if (mv[i1] == mv[i2]) continue;
        const int m = mv[i1] ? i1 : i2, s = mv[i1] ? i2 : i1;
        const Vec3 rp = perp(P[m] - o, k), tp = perp(P[s] - o, k);
        num += dot(k, cross(rp, tp));
        den += dot(rp, tp);
        copyTerm = true;
      }
      for (std::size_t q = 0; q < l.restraints.size(); ++q) {
        const int im = l.restraintPoint[2 * q], is = l.restraintPoint[2 * q + 1];
        if (mv[im] == mv[is]) continue;
        const int m = mv[im] ? im : is, s = mv[im] ? is : im;
        const Vec3 r = P[m] - o, tv = P[s] - o;
        const Vec3 rp = perp(r, k), tp = perp(tv, k);
        const double dn = dot(rp, tp), nm = dot(k, cross(rp, tp));
        const double B = std::hypot(dn, nm);
        if (B <= 1e-12) continue;
        rt.push_back({dot(rp, rp) + dot(tp, tp) + std::pow(dot(r, k) - dot(tv, k), 2), B, std::atan2(nm, dn),
                      l.restraints[q].r0});
      }
      // Clash candidates: pairs with one atom moved whose circle comes closer than c to the other atom.
      struct CTerm {
        Vec3 m, s;
      };
      std::vector<CTerm> ct;
      if (useClash) {
        if (++nlAge >= 20) buildNl();
        const auto& ms = l.movesSys[d];
        for (const auto& [i, j] : nl) {
          const bool mi = std::binary_search(ms.begin(), ms.end(), i), mj = std::binary_search(ms.begin(), ms.end(), j);
          if (mi == mj) continue;
          const Vec3& m = X[mi ? i : j];
          const Vec3& s = X[mi ? j : i];
          const Vec3 rm = m - o, rs = s - o;
          const double h = dot(rs - rm, k), rho = norm(perp(rm, k)), rad = norm(perp(rs, k));
          if (h * h + (rad - rho) * (rad - rho) < cd * cd) ct.push_back({m, s});
        }
      }
      double delta = 0;
      if (!ct.empty()) {
        // Closure terms plus w (c - d)^2: 1-degree scan, golden-section refinement; kept only when lower than at 0.
        auto f = [&](double x) {
          double v = -2 * (den * std::cos(x) + num * std::sin(x));
          for (const RTerm& r : rt) {
            const double e = std::sqrt(std::max(0.0, r.A - 2 * r.B * std::cos(x - r.d0))) - r.r0;
            v += e * e;
          }
          const double c = std::cos(x), sn = std::sin(x);
          for (const CTerm& q : ct) {
            const double dd = norm(rotateAbout(q.m, o, k, c, sn) - q.s);
            if (dd < cd) v += w * (cd - dd) * (cd - dd);
          }
          return v;
        };
        const double f0 = f(0), st = std::numbers::pi / 180;
        double best = 0, fb = f0;
        for (int s = -179; s <= 180; ++s)
          if (const double v = f(s * st); v < fb) {
            fb = v;
            best = s * st;
          }
        double lo = best - st, hi = best + st;
        constexpr double g = 0.6180339887498949;
        for (int i2 = 0; i2 < 40; ++i2) {
          const double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
          if (f(x1) < f(x2)) hi = x2;
          else lo = x1;
        }
        if (const double xx = 0.5 * (lo + hi); f(xx) < fb) best = xx;
        delta = f(best) < f0 ? best : 0;
      } else if (rt.empty()) {
        if (copyTerm && (num != 0 || den != 0)) delta = std::atan2(num, den);
      } else if (!copyTerm && rt.size() == 1) {
        const RTerm& r = rt[0];
        const double c = (r.A - r.r0 * r.r0) / (2 * r.B);
        if (c >= 1) delta = r.d0;
        else if (c <= -1) delta = r.d0 + std::numbers::pi;
        else {
          const double w = std::acos(c);
          const double x1 = std::remainder(r.d0 + w, 2 * std::numbers::pi), x2 = std::remainder(r.d0 - w, 2 * std::numbers::pi);
          delta = std::abs(x1) <= std::abs(x2) ? x1 : x2;
        }
      } else {
        // Several terms: 1-degree scan, then golden-section refinement; kept only when lower than at 0.
        auto f = [&](double x) {
          double v = -2 * (den * std::cos(x) + num * std::sin(x));
          for (const RTerm& r : rt) {
            const double e = std::sqrt(std::max(0.0, r.A - 2 * r.B * std::cos(x - r.d0))) - r.r0;
            v += e * e;
          }
          return v;
        };
        const double f0 = f(0), st = std::numbers::pi / 180;
        double best = 0, fb = f0;
        for (int s = -179; s <= 180; ++s)
          if (const double v = f(s * st); v < fb) {
            fb = v;
            best = s * st;
          }
        double lo = best - st, hi = best + st;
        constexpr double g = 0.6180339887498949;
        for (int it = 0; it < 40; ++it) {
          const double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
          if (f(x1) < f(x2)) hi = x2;
          else lo = x1;
        }
        if (const double x = 0.5 * (lo + hi); f(x) < fb) best = x;
        delta = f(best) < f0 ? best : 0;
      }
      if (delta != 0) {
        const double c = std::cos(delta), s = std::sin(delta);
        for (int i : l.moves[d]) P[i] = rotateAbout(P[i], o, k, c, s);
        if (clashTerm)
          for (int i : l.movesSys[d]) X[i] = rotateAbout(X[i], o, k, c, s);
        a[l.dof[d]] += delta * kDeg;
      }
      for (int i : l.moves[d]) mv[i] = 0;
      conv = converged();
      if (stallTol > 0) {
        if (std::abs(delta) > stallTol) {
          std::fill(tried.begin(), tried.end(), 0);
          triedNum = 0;
        }
        if (!tried[d]) {
          tried[d] = 1;
          ++triedNum;
        }
        if (triedNum == nd) break;
      }
    }
    for (int g : l.dof) a[g] = wrapDeg(a[g]);
    if (l.partial) return true;
    if (!conv) continue;
    if (deviation(l, partCoords(l.part, a)) > tol + 1e-6) continue;
    if (!clashFree(l, a)) continue;
    return true;
  }
  std::copy(orig.begin(), orig.end(), a.begin());
  return false;
}

bool LoopCloser::kic(Rng& rng, const Link& l, std::span<double> a, const SegmentSampler& resample) const {
  const std::vector<double> orig(a.begin(), a.end());
  const auto& tree = *parts_[l.part].tree;
  const int n = static_cast<int>(l.ring.size()) - 1;  // ring indices 0 (a) ... n (b)
  // Pivot candidates: f3 / f2 = a valid later pivot / pair of later pivots exists.
  const auto& cand = l.pivots;
  const int nc = static_cast<int>(cand.size());
  std::vector<char> f3(nc, 0), f2(nc, 0);
  for (int i = 0; i < nc; ++i)
    for (int j = i + 1; j < nc; ++j)
      if (cand[j] >= cand[i] + 2) f3[i] = 1;
  for (int i = 0; i < nc; ++i)
    for (int j = i + 1; j < nc; ++j)
      if (cand[j] >= cand[i] + 2 && f3[j]) f2[i] = 1;
  auto pick = [&](auto ok) {
    std::vector<int> v;
    for (int i = 0; i < nc; ++i)
      if (ok(i)) v.push_back(i);
    return v.empty() ? -1 : v[rng.below(v.size())];
  };
  const int vA = l.pairs.back().moving, vB = l.pairs[0].moving, vY = l.pairs[1].moving;
  const int A = l.pairs.back().target, B = l.pairs[0].target, Y = l.pairs[1].target;
  const double tol = p_.rcdTolerance;

  for (int trial = 0; trial < p_.kicMaxTrial; ++trial) {
    if (trial > 0) {
      std::copy(orig.begin(), orig.end(), a.begin());
      resampleLink(rng, l, a, resample);
    }
    const int c1 = pick([&](int i) { return f2[i] != 0; });
    if (c1 < 0) break;
    const int c2 = pick([&](int i) { return cand[i] >= cand[c1] + 2 && f3[i]; });
    const int c3 = pick([&](int i) { return cand[i] >= cand[c2] + 2; });
    const int i1 = cand[c1], i2 = cand[c2], i3 = cand[c3];
    const auto x = partCoords(l.part, a);
    auto ringX = [&](int idx) { return x[l.ring[idx]]; };
    // B-side atoms into the frame of the a side: {A', B, Y} -> {A, B', Y'}.
    const Frame src = frameOf(x[B], x[vA], x[Y]), dst = frameOf(x[vB], x[A], x[vY]);
    // Traversal from pivot i1 down to a, across the break, from b down to i1 again.
    std::vector<int> seq;
    for (int i = i1; i >= 0; --i) seq.push_back(i);
    for (int i = n; i >= i1; --i) seq.push_back(i);
    std::vector<Vec3> R(seq.size());
    for (std::size_t s = 0; s < seq.size(); ++s)
      R[s] = static_cast<int>(s) <= i1 ? ringX(seq[s]) : mapFrame(src, dst, ringX(seq[s]));
    const int sb = i1 + 1 + (n - i3), sc = i1 + 1 + (n - i2), sd = n + 1;
    auto pivotCos = [&](int i) { return angleCos(ringX(i - 1), ringX(i), ringX(i + 1)); };
    const double cos1 = pivotCos(i1), cos2 = pivotCos(i3), cos3 = pivotCos(i2);
    const Vec3 P1 = R[0], P2 = R[sb];
    const Vec3 out1 = unit(R[1] - P1), in2 = unit(R[sb - 1] - P2);
    const Vec3 p2l = R[sb], p3l = R[sc], out2l = unit(R[sb + 1] - p2l), in3l = unit(R[sc - 1] - p3l);
    const Vec3 p3m = R[sc], p1m = R[sd], out3l = unit(R[sc + 1] - p3m), in1l = unit(R[sd - 1] - p1m);
    const double d12 = norm(P2 - P1), d23 = norm(p3l - p2l), d31 = norm(p1m - p3m);
    if (d12 < 1e-9) continue;
    const double xx = (d12 * d12 + d31 * d31 - d23 * d23) / (2 * d12);
    const double rho2 = d31 * d31 - xx * xx;
    if (rho2 <= 1e-12) continue;
    const double rho = std::sqrt(rho2);
    const Vec3 ax = unit(P2 - P1);
    Vec3 mref = perp(out1, ax);
    if (norm(mref) < 1e-6) mref = perp(Vec3{1, 0, 0}, ax);
    if (norm(mref) < 1e-6) mref = perp(Vec3{0, 1, 0}, ax);
    const Vec3 mm = unit(mref), nn = cross(ax, mm);
    const Vec3 e2l = unit(p3l - p2l), f2l = unit(perp(out2l, e2l)), h2l = cross(e2l, f2l);
    const double cc2 = dot(out2l, e2l), ss2 = norm(perp(out2l, e2l));
    const Vec3 e3l = unit(p1m - p3m), f3l = unit(perp(in1l, e3l)), h3l = cross(e3l, f3l);
    const double cc3 = dot(in1l, e3l), ss3 = norm(perp(in1l, e3l));
    if (ss2 < 1e-6 || ss3 < 1e-6) continue;

    struct Place {
      Vec3 P3, e2, f2, g2, e3, f3, g3;
    };
    // Arc placement for tau and the branch signs; false where the pivot equations have no solution.
    auto solveAt = [&](double tau, int sg2, int sg3, double& F, Place& pl) {
      pl.P3 = P1 + xx * ax + rho * (std::cos(tau) * mm + std::sin(tau) * nn);
      pl.e2 = unit(pl.P3 - P2);
      const Vec3 f2w = unit(perp(P1 - P2, pl.e2)), g2w = cross(pl.e2, f2w);
      const double A2 = ss2 * dot(in2, f2w), B2 = ss2 * dot(in2, g2w), C2 = cos2 - cc2 * dot(in2, pl.e2);
      const double R2 = std::hypot(A2, B2);
      if (R2 < 1e-12 || std::abs(C2) > R2) return false;
      const double s2 = std::atan2(B2, A2) + sg2 * std::acos(C2 / R2);
      pl.f2 = std::cos(s2) * f2w + std::sin(s2) * g2w;
      pl.g2 = cross(pl.e2, pl.f2);
      pl.e3 = unit(P1 - pl.P3);
      const Vec3 f3w = unit(perp(P2 - pl.P3, pl.e3)), g3w = cross(pl.e3, f3w);
      const double A3 = ss3 * dot(out1, f3w), B3 = ss3 * dot(out1, g3w), C3 = cos1 - cc3 * dot(out1, pl.e3);
      const double R3 = std::hypot(A3, B3);
      if (R3 < 1e-12 || std::abs(C3) > R3) return false;
      const double s3 = std::atan2(B3, A3) + sg3 * std::acos(C3 / R3);
      pl.f3 = std::cos(s3) * f3w + std::sin(s3) * g3w;
      pl.g3 = cross(pl.e3, pl.f3);
      const Vec3 in3w = dot(e2l, in3l) * pl.e2 + dot(f2l, in3l) * pl.f2 + dot(h2l, in3l) * pl.g2;
      const Vec3 out3w = dot(e3l, out3l) * pl.e3 + dot(f3l, out3l) * pl.f3 + dot(h3l, out3l) * pl.g3;
      F = dot(in3w, out3w) - cos3;
      return true;
    };

    constexpr int kScan = 720;
    const double step = 2 * std::numbers::pi / kScan;
    std::vector<Place> sols;
    for (int sg2 : {1, -1})
      for (int sg3 : {1, -1}) {
        std::vector<double> F(kScan);
        std::vector<char> ok(kScan);
        Place pl;
        for (int s = 0; s < kScan; ++s) ok[s] = solveAt(s * step, sg2, sg3, F[s], pl);
        for (int s = 0; s < kScan; ++s) {
          const int t = (s + 1) % kScan;
          if (!ok[s] || !ok[t] || (F[s] < 0) == (F[t] < 0)) continue;
          double lo = s * step, hi = lo + step, flo = F[s];
          bool good = true;
          for (int it = 0; it < 60 && good; ++it) {
            const double mid = 0.5 * (lo + hi);
            double fm;
            if (!solveAt(mid, sg2, sg3, fm, pl)) {
              good = false;
              break;
            }
            if ((fm < 0) == (flo < 0)) {
              lo = mid;
              flo = fm;
            } else {
              hi = mid;
            }
          }
          double fr;
          if (good && solveAt(0.5 * (lo + hi), sg2, sg3, fr, pl) && std::abs(fr) < 1e-8) sols.push_back(pl);
        }
      }
    if (sols.empty()) continue;
    const Place& pl = sols[rng.below(sols.size())];
    // Closed ring positions by ring index.
    std::vector<Vec3> ring(n + 1);
    for (int s = 0; s <= sd; ++s) {
      Vec3 y = R[s];
      if (s > sb && s <= sc) {
        const Vec3 v = R[s] - p2l;
        y = P2 + dot(e2l, v) * pl.e2 + dot(f2l, v) * pl.f2 + dot(h2l, v) * pl.g2;
      } else if (s > sc) {
        const Vec3 v = R[s] - p3m;
        y = pl.P3 + dot(e3l, v) * pl.e3 + dot(f3l, v) * pl.f3 + dot(h3l, v) * pl.g3;
      }
      if (s < sd) ring[seq[s]] = y;
    }
    auto solved = [&](int idx) { return idx < 0 ? ring[n] : idx > n ? ring[0] : ring[idx]; };
    auto current = [&](int idx) { return idx < 0 ? x[vB] : idx > n ? x[vA] : ringX(idx); };
    auto currentIndex = [&](int idx) { return idx < 0 ? vB : idx > n ? vA : l.ring[idx]; };
    for (int i : {i1, i2, i3})
      for (int j : {i - 1, i}) {
        const int g = l.bondTorsion[j];
        const auto& tt = tree.torsions[localTorsion_[g]];
        const double dNew = dihedralDeg(solved(j - 1), solved(j), solved(j + 1), solved(j + 2));
        const double dOld = dihedralDeg(current(j - 1), current(j), current(j + 1), current(j + 2));
        // Sign of the dihedral change when the moving side turns by +1 degree.
        const Vec3 o = x[tt.fixedAtom], k = unit(x[tt.movingAtom] - o);
        const double c = std::cos(1.0 / kDeg), s = std::sin(1.0 / kDeg);
        Vec3 q0 = current(j - 1), q3 = current(j + 2);
        if (contains(tt.moving, currentIndex(j + 2))) q3 = rotateAbout(q3, o, k, c, s);
        else q0 = rotateAbout(q0, o, k, c, s);
        const double sign = wrapDeg(dihedralDeg(q0, current(j), current(j + 1), q3) - dOld) > 0 ? 1.0 : -1.0;
        a[g] = wrapDeg(a[g] + sign * wrapDeg(dNew - dOld));
      }
    if (deviation(l, partCoords(l.part, a)) > tol) continue;
    if (!clashFree(l, a)) continue;
    return true;
  }
  std::copy(orig.begin(), orig.end(), a.begin());
  return false;
}

bool LoopCloser::close(Rng& rng, std::span<double> a, const SegmentSampler& resample) const {
  for (const Link& l : links_) {
    if (l.dof.empty()) continue;
    if (deviation(l, partCoords(l.part, a)) <= p_.rcdTolerance) continue;
    const bool ok = p_.method == ClosureMethod::Kic && l.kic && !l.partial ? kic(rng, l, a, resample) : rcd(rng, l, a, resample);
    if (!ok) return false;
  }
  return true;
}

}  // namespace od
