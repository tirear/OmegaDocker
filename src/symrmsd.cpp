// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "symrmsd.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_map>

#include "log.hpp"

namespace od {

namespace {

double squared(const Vec3f& p, const Vec3f& q) {
  const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
  return dx * dx + dy * dy + dz * dz;
}

// Colour refinement: vertices with different colours are never mapped onto each other.
std::vector<int> refineColours(const std::vector<std::vector<int>>& adj, std::span<const int> label) {
  const std::size_t n = adj.size();
  std::vector<int> colour(n);
  {
    std::map<std::pair<int, std::size_t>, int> id;
    for (std::size_t v = 0; v < n; ++v) colour[v] = id.try_emplace({label[v], adj[v].size()}, static_cast<int>(id.size())).first->second;
  }
  for (std::size_t classes = 0;;) {
    std::map<std::vector<int>, int> id;
    std::vector<int> next(n);
    for (std::size_t v = 0; v < n; ++v) {
      std::vector<int> key{colour[v]};
      for (int w : adj[v]) key.push_back(colour[w]);
      std::sort(key.begin() + 1, key.end());
      next[v] = id.try_emplace(std::move(key), static_cast<int>(id.size())).first->second;
    }
    colour = std::move(next);
    if (id.size() == classes) break;
    classes = id.size();
  }
  return colour;
}

}  // namespace

double fixedRmsd(std::span<const Vec3f> a, std::span<const Vec3f> b) {
  if (a.empty()) return 0;
  double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i) s += squared(a[i], b[i]);
  return std::sqrt(s / static_cast<double>(a.size()));
}

namespace {

// Isomorphisms from graph A onto graph B (at most limit + 1); map[i]: image in B of vertex i of A.
std::vector<std::vector<int>> isomorphisms(const std::vector<std::vector<int>>& adjA, std::span<const int> labelA,
                                           const std::vector<std::vector<int>>& adjB, std::span<const int> labelB,
                                           std::size_t limit) {
  const int n = static_cast<int>(adjA.size());
  std::vector<std::vector<int>> out;
  if (static_cast<int>(adjB.size()) != n) return out;
  if (n == 0) return {{}};
  // Colours refined on the disjoint union, so that they are comparable between A and B.
  std::vector<std::vector<int>> adjU(adjA);
  std::vector<int> labelU(labelA.begin(), labelA.end());
  for (const auto& nb : adjB) {
    adjU.emplace_back();
    for (int w : nb) adjU.back().push_back(w + n);
  }
  labelU.insert(labelU.end(), labelB.begin(), labelB.end());
  const std::vector<int> colourU = refineColours(adjU, labelU);
  const std::vector<int> colourA(colourU.begin(), colourU.begin() + n), colourB(colourU.begin() + n, colourU.end());
  {
    std::vector<int> a = colourA, b = colourB;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    if (a != b) return out;
  }
  // Breadth-first order of A; parent[k]: an earlier neighbour of order[k], -1 for component roots.
  std::vector<int> order, parent;
  std::vector<char> seen(n, 0);
  for (int r = 0; r < n; ++r) {
    if (seen[r]) continue;
    seen[r] = 1;
    order.push_back(r);
    parent.push_back(-1);
    for (std::size_t h = order.size() - 1; h < order.size(); ++h)
      for (int w : adjA[order[h]])
        if (!seen[w]) {
          seen[w] = 1;
          order.push_back(w);
          parent.push_back(order[h]);
        }
  }
  std::vector<std::vector<char>> isAdjA(n, std::vector<char>(n, 0)), isAdjB(n, std::vector<char>(n, 0));
  for (int v = 0; v < n; ++v) {
    for (int w : adjA[v]) isAdjA[v][w] = 1;
    for (int w : adjB[v]) isAdjB[v][w] = 1;
  }
  std::vector<int> map(n, -1);
  std::vector<char> used(n, 0);
  auto fits = [&](int k, int u, int v) {
    if (used[v] || colourB[v] != colourA[u]) return false;
    for (int j = 0; j < k; ++j)
      if (isAdjA[u][order[j]] != isAdjB[v][map[order[j]]]) return false;
    return true;
  };
  auto rec = [&](auto&& self, int k) -> void {
    if (out.size() > limit) return;
    if (k == n) {
      out.push_back(map);
      return;
    }
    const int u = order[k];
    auto tryV = [&](int v) {
      if (!fits(k, u, v)) return;
      map[u] = v;
      used[v] = 1;
      self(self, k + 1);
      used[v] = 0;
      map[u] = -1;
    };
    if (parent[k] >= 0)
      for (int v : adjB[map[parent[k]]]) tryV(v);
    else
      for (int v = 0; v < n; ++v) tryV(v);
  };
  rec(rec, 0);
  return out;
}

// Adjacency of the atoms (positions in atoms), metal coordination bonds excluded.
std::vector<std::vector<int>> atomGraph(const Molecule& mol, std::span<const int> atoms) {
  std::unordered_map<int, int> pos;
  for (std::size_t k = 0; k < atoms.size(); ++k) pos[atoms[k]] = static_cast<int>(k);
  std::vector<std::vector<int>> adj(atoms.size());
  for (const Bond& b : mol.bonds) {
    if (isMetalBond(mol, b)) continue;
    const auto ia = pos.find(b.a), ib = pos.find(b.b);
    if (ia == pos.end() || ib == pos.end()) continue;
    adj[ia->second].push_back(ib->second);
    adj[ib->second].push_back(ia->second);
  }
  return adj;
}

}  // namespace

std::vector<std::vector<int>> graphAutomorphisms(const std::vector<std::vector<int>>& adj, std::span<const int> label,
                                                 std::size_t limit) {
  return isomorphisms(adj, label, adj, label, limit);
}

std::vector<int> graphIsomorphism(const std::vector<std::vector<int>>& adjA, std::span<const int> labelA,
                                  const std::vector<std::vector<int>>& adjB, std::span<const int> labelB) {
  auto v = isomorphisms(adjA, labelA, adjB, labelB, 0);
  return v.empty() ? std::vector<int>{} : std::move(v.front());
}

std::vector<int> referenceAtomMap(const Molecule& mol, std::span<const int> atoms, const Molecule& ref) {
  std::vector<int> refAtoms;
  for (std::size_t i = 0; i < ref.atoms.size(); ++i)
    if (ref.atoms[i].atomicNum != 1) refAtoms.push_back(static_cast<int>(i));
  if (refAtoms.size() != atoms.size())
    fail("{}: {} heavy atoms, {} has {}", ref.id, refAtoms.size(), mol.id, atoms.size());
  std::vector<int> la, lb;
  for (int a : atoms) la.push_back(mol.atoms[a].atomicNum);
  for (int a : refAtoms) lb.push_back(ref.atoms[a].atomicNum);
  const auto adjA = atomGraph(mol, atoms), adjB = atomGraph(ref, refAtoms);
  // The same atom order first (the usual case), then a search.
  bool same = la == lb;
  for (std::size_t k = 0; same && k < adjA.size(); ++k) {
    auto x = adjA[k], y = adjB[k];
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    same = x == y;
  }
  std::vector<int> map;
  if (same) {
    map.resize(atoms.size());
    std::iota(map.begin(), map.end(), 0);
  } else {
    map = graphIsomorphism(adjA, la, adjB, lb);
    if (map.empty()) fail("{}: the heavy-atom bond graph differs from {}", ref.id, mol.id);
  }
  for (int& m : map) m = refAtoms[m];
  return map;
}

double hungarianSquaredSum(std::span<const Vec3f> a, std::span<const Vec3f> b,
                           const std::vector<std::vector<int>>& groups) {
  double total = 0;
  for (const auto& g : groups) {
    const int n = static_cast<int>(g.size());
    if (n == 1) {
      total += squared(a[g[0]], b[g[0]]);
      continue;
    }
    std::vector<double> cost(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) cost[i * n + j] = squared(a[g[i]], b[g[j]]);
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> u(n + 1, 0), v(n + 1, 0), minv(n + 1);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);
    std::vector<char> usedCol(n + 1);
    for (int i = 1; i <= n; ++i) {
      p[0] = i;
      int j0 = 0;
      std::fill(minv.begin(), minv.end(), inf);
      std::fill(usedCol.begin(), usedCol.end(), 0);
      do {
        usedCol[j0] = 1;
        const int i0 = p[j0];
        double delta = inf;
        int j1 = 0;
        for (int j = 1; j <= n; ++j) {
          if (usedCol[j]) continue;
          const double cur = cost[(i0 - 1) * n + j - 1] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1 = j;
          }
        }
        for (int j = 0; j <= n; ++j) {
          if (usedCol[j]) {
            u[p[j]] += delta;
            v[j] -= delta;
          } else {
            minv[j] -= delta;
          }
        }
        j0 = j1;
      } while (p[j0] != 0);
      do {
        const int j1 = way[j0];
        p[j0] = p[j1];
        j0 = j1;
      } while (j0);
    }
    for (int j = 1; j <= n; ++j) total += cost[(p[j] - 1) * n + j - 1];
  }
  return total;
}

SymmetryRmsd::SymmetryRmsd(const Molecule& mol, std::span<const int> atoms, std::span<const int> type,
                           std::size_t maxPerm) {
  const std::size_t n = atoms.size();
  const std::vector<std::vector<int>> adj = atomGraph(mol, atoms);
  if (maxPerm >= 1) {
    auto perms = graphAutomorphisms(adj, type, maxPerm);
    if (perms.size() <= maxPerm) perms_ = std::move(perms);
  }
  std::map<int, std::vector<int>> byType;
  for (std::size_t k = 0; k < n; ++k) byType[type[k]].push_back(static_cast<int>(k));
  for (auto& [t, g] : byType) groups_.push_back(std::move(g));
}

double SymmetryRmsd::rmsd(std::span<const Vec3f> a, std::span<const Vec3f> b) const {
  if (a.empty()) return 0;
  double best;
  if (perms_.empty()) {
    best = hungarianSquaredSum(a, b, groups_);
  } else {
    best = std::numeric_limits<double>::infinity();
    for (const auto& p : perms_) {
      double s = 0;
      for (std::size_t i = 0; i < a.size() && s < best; ++i) s += squared(a[i], b[p[i]]);
      best = std::min(best, s);
    }
  }
  return std::sqrt(best / static_cast<double>(a.size()));
}

}  // namespace od
