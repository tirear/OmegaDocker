// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "energy.hpp"

#include <algorithm>
#include <cmath>
#include <format>

#include "ad4_pair.hpp"
#include "log.hpp"

namespace od {
namespace {

std::string atomRef(const Molecule& m, int a) {
  const Atom& x = m.atoms[a];
  if (x.residue < 0 || x.name.empty()) return std::format("#{} {}", x.fileSerial, elementSymbol(x.atomicNum));
  const auto& r = m.residues[x.residue].info;
  return std::format("{}:{}{}:{}", r.chain, r.seq, r.iCode, x.name);
}

// Rigid groups an atom keeps a constant position in: its own group, plus the
// group across a tree torsion whose axis it lies on (AutoDock-GPU
// get_intraE_contributors, cases 1, 3 and 4).
using Groups = std::vector<std::vector<int>>;

Groups rigidGroups(const Molecule& m, const TorsionTree& tree, int rootGroup, int& nextGroup) {
  const int n = static_cast<int>(m.atoms.size());
  std::vector<int> owner(n, -1);
  for (std::size_t i = 0; i < tree.torsions.size(); ++i)
    for (int a : tree.torsions[i].moving)
      if (a < n) owner[a] = static_cast<int>(i);
  const int base = nextGroup;
  nextGroup += static_cast<int>(tree.torsions.size());
  auto group = [&](int t) { return t < 0 ? rootGroup : base + t; };
  Groups g(n);
  for (int a = 0; a < n; ++a) g[a].push_back(group(owner[a]));
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    const auto& tt = tree.torsions[i];
    g[tt.movingAtom].push_back(group(owner[tt.fixedAtom]));
    g[tt.fixedAtom].push_back(group(static_cast<int>(i)));
  }
  for (auto& v : g) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  return g;
}

bool intersect(const std::vector<int>& a, const std::vector<int>& b) {
  for (std::size_t i = 0, j = 0; i < a.size() && j < b.size();) {
    if (a[i] == b[j]) return true;
    if (a[i] < b[j]) ++i;
    else ++j;
  }
  return false;
}

}  // namespace

EnergySystem::EnergySystem(std::vector<EnergyPart> parts, const EnergyParams& params)
    : parts_(std::move(parts)), params_(params) {
  constexpr int kTypes = static_cast<int>(Ad4Type::Count);
  pairParams_.resize(kTypes * kTypes);
  for (int a = 0; a < kTypes; ++a)
    for (int b = 0; b < kTypes; ++b)
      pairParams_[a * kTypes + b] = ad4PairParams(static_cast<Ad4Type>(a), static_cast<Ad4Type>(b), params_.coeffs);

  std::string unsupported;
  for (const EnergyPart& p : parts_) {
    typing_.push_back(assignAd4Types(*p.mol));
    for (int a : typing_.back().unsupported)
      unsupported += std::format("\n  {} {}", p.mol->id, atomRef(*p.mol, a));
    offset_.push_back(coordCount_);
    torsionOffset_.push_back(torsionCount_);
    coordCount_ += static_cast<int>(p.mol->atoms.size() + p.tree->virtuals.size());
    torsionCount_ += static_cast<int>(p.tree->torsions.size());
  }
  if (!unsupported.empty()) fail("no AD4 parameters for the element of:{}", unsupported);

  charge_.assign(coordCount_, 0.0);
  sv_.assign(coordCount_, 0.0);
  vol_.assign(coordCount_, 0.0);
  struct EAtom {
    int part, atom, index;
    Ad4Type type;
    bool heavy;
  };
  std::vector<EAtom> atoms;
  std::vector<Groups> groups;
  int nextGroup = 1;  // 0: the fixed frame shared by all StaM roots
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    const Molecule& m = *parts_[k].mol;
    const int root = m.role == Role::StaM ? 0 : nextGroup++;
    groups.push_back(rigidGroups(m, *parts_[k].tree, root, nextGroup));
    for (int a = 0; a < static_cast<int>(m.atoms.size()); ++a) {
      const Ad4Type t = typing_[k].type[a];
      if (t == Ad4Type::Merged) continue;
      const int idx = offset_[k] + a;
      const double q = typing_[k].charge[a];
      charge_[idx] = q;
      sv_[idx] = ad4Params(t).solpar + params_.coeffs.qasp * std::abs(q);
      vol_[idx] = ad4Params(t).volume;
      atoms.push_back({static_cast<int>(k), a, idx, t, m.atoms[a].atomicNum > 1});
    }
  }

  // Atoms within three bonds (all bonds, including broken and metal bonds) are excluded.
  std::vector<std::vector<std::vector<int>>> near(parts_.size());
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    const Molecule& m = *parts_[k].mol;
    const int n = static_cast<int>(m.atoms.size());
    std::vector<std::vector<int>> adj(n);
    for (const Bond& b : m.bonds) {
      adj[b.a].push_back(b.b);
      adj[b.b].push_back(b.a);
    }
    near[k].resize(n);
    std::vector<int> depth(n, -1);
    for (int s = 0; s < n; ++s) {
      std::vector<int> q{s};
      depth[s] = 0;
      for (std::size_t i = 0; i < q.size(); ++i)
        if (depth[q[i]] < 3)
          for (int o : adj[q[i]])
            if (depth[o] < 0) {
              depth[o] = depth[q[i]] + 1;
              q.push_back(o);
            }
      for (int a : q) {
        if (a > s) near[k][s].push_back(a);
        depth[a] = -1;
      }
      std::sort(near[k][s].begin(), near[k][s].end());
    }
  }

  std::vector<std::string> fixedClashes;
  const double cd2 = params_.clashDistance * params_.clashDistance;
  for (std::size_t x = 0; x < atoms.size(); ++x) {
    const EAtom& p = atoms[x];
    const auto& gp = groups[p.part][p.atom];
    for (std::size_t y = x + 1; y < atoms.size(); ++y) {
      const EAtom& q = atoms[y];
      if (p.part == q.part && std::binary_search(near[p.part][p.atom].begin(), near[p.part][p.atom].end(), q.atom))
        continue;
      if (intersect(gp, groups[q.part][q.atom])) {
        if (!p.heavy || !q.heavy) continue;
        const Molecule& mp = *parts_[p.part].mol;
        const Molecule& mq = *parts_[q.part].mol;
        const Vec3& a = mp.atoms[p.atom].xyz;
        const Vec3& b = mq.atoms[q.atom].xyz;
        const double r2 = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
        if (r2 < cd2)
          fixedClashes.push_back(std::format("\n  {} {} - {} {}: {:.3f} Å", mp.id, atomRef(mp, p.atom), mq.id,
                                             atomRef(mq, q.atom), std::sqrt(r2)));
        continue;
      }
      pairs_.push_back({p.index, q.index, static_cast<std::uint16_t>(static_cast<int>(p.type) * kTypes + static_cast<int>(q.type)),
                        p.heavy && q.heavy});
    }
  }
  if (!fixedClashes.empty()) {
    std::string list;
    for (std::size_t i = 0; i < fixedClashes.size() && i < 50; ++i) list += fixedClashes[i];
    if (fixedClashes.size() > 50) list += std::format("\n  ... ({} more)", fixedClashes.size() - 50);
    fail("input structure: {} pair(s) of heavy atoms closer than --clashDistance {} Å whose distance cannot change "
         "(make torsions between them rotatable or fix the input):{}",
         fixedClashes.size(), params_.clashDistance, list);
  }
}

std::vector<double> EnergySystem::inputAngles() const {
  std::vector<double> v;
  v.reserve(torsionCount_);
  for (const EnergyPart& p : parts_)
    for (const auto& tt : p.tree->torsions) v.push_back(tt.inputAngle);
  return v;
}

template <class T>
std::vector<Vec3T<T>> EnergySystem::coordinatesT(std::span<const T> anglesDeg) const {
  std::vector<Vec3T<T>> x(coordCount_);
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    const auto& tree = *parts_[k].tree;
    const auto xyz = buildCoordinates(*parts_[k].mol, tree, anglesDeg.subspan(torsionOffset_[k], tree.torsions.size()));
    std::copy(xyz.begin(), xyz.end(), x.begin() + offset_[k]);
  }
  return x;
}

std::vector<Vec3f> EnergySystem::coordinates(std::span<const float> anglesDeg) const { return coordinatesT(anglesDeg); }

template <class T>
EnergyTerms EnergySystem::evaluateT(std::span<const T> anglesDeg, std::vector<double>* gradient,
                                    std::vector<Vec3T<T>>* coords, const ExternalTermT<T>* external) const {
  using V = Vec3T<T>;
  std::vector<V> x = coordinatesT(anglesDeg);
  const bool grad = gradient != nullptr;
  // Accumulators (energy sums, per-coordinate gradients) are double.
  std::vector<Vec3> g(grad ? coordCount_ : 0, Vec3{0, 0, 0});
  auto addPairGrad = [&](int i, int j, const V& d, T r, T dedr) {
    if (!grad || r <= 0) return;
    const T f = dedr / r;
    for (int c = 0; c < 3; ++c) {
      g[i][c] += f * d[c];
      g[j][c] -= f * d[c];
    }
  };

  EnergyTerms e;
  const T half = static_cast<T>(0.5 * params_.smooth);
  const T cd = static_cast<T>(params_.clashDistance), cd2 = cd * cd;
  const T elecCoeff = static_cast<T>(params_.coeffs.elec), desolvCoeff = static_cast<T>(params_.coeffs.desolv);
  const T rMin = static_cast<T>(params_.elecMinDistance);
  const T clashBase = static_cast<T>(params_.clashBase), clashK = static_cast<T>(params_.clashK);
  constexpr T kVdwCut = static_cast<T>(kVdwCutoff), kElecCut2 = static_cast<T>(kElecCutoff * kElecCutoff);
  for (const Pair& p : pairs_) {
    const V d{x[p.i][0] - x[p.j][0], x[p.i][1] - x[p.j][1], x[p.i][2] - x[p.j][2]};
    const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (r2 >= kElecCut2) continue;
    const T r = std::sqrt(r2);
    T dedr = 0, ee, de;
    if (r < kVdwCut) {
      const Ad4PairParams& pp = pairParams_[p.type];
      ad4VdwHb<T>(r, static_cast<T>(pp.AC), static_cast<T>(pp.BD), pp.m, pp.n, static_cast<T>(pp.reqm), half, ee, de);
      (pp.n == 10 ? e.hbond : e.vdw) += ee;
      dedr += de;
    }
    ad4Elec<T>(r, elecCoeff * static_cast<T>(charge_[p.i]) * static_cast<T>(charge_[p.j]), rMin, ee, de);
    e.elec += ee;
    dedr += de;
    ad4Desolv<T>(r,
                 desolvCoeff * (static_cast<T>(sv_[p.i]) * static_cast<T>(vol_[p.j]) +
                                static_cast<T>(sv_[p.j]) * static_cast<T>(vol_[p.i])),
                 rMin, ee, de);
    e.desolv += ee;
    dedr += de;
    if (p.heavy && r2 < cd2) {
      clashPenalty<T>(r, cd, clashBase, clashK, ee, de);
      e.clash += ee;
      ++e.clashPairs;
      dedr += de;
    }
    addPairGrad(p.i, p.j, d, r, dedr);
  }

  const T G = static_cast<T>(params_.closureG);
  for (std::size_t k = 0; k < parts_.size(); ++k) {
    const auto& tree = *parts_[k].tree;
    const int o = offset_[k];
    for (std::size_t v = 0; v < tree.virtuals.size(); ++v) {
      const int i = o + tree.atomCount + static_cast<int>(v), j = o + tree.virtuals[v].copyOf;
      const V d{x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]};
      const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      e.closure += G * r2;
      e.closureMaxDeviation = std::max(e.closureMaxDeviation, static_cast<double>(std::sqrt(r2)));
      if (grad)
        for (int c = 0; c < 3; ++c) {
          g[i][c] += 2 * G * d[c];
          g[j][c] -= 2 * G * d[c];
        }
    }
    for (const auto& rs : tree.restraints) {
      const int i = o + rs.a, j = o + rs.b;
      const V d{x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]};
      const T r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      const T dd = r - static_cast<T>(rs.r0);
      e.closure += G * dd * dd;
      e.closureMaxDeviation = std::max(e.closureMaxDeviation, static_cast<double>(std::abs(dd)));
      addPairGrad(i, j, d, r, 2 * G * dd);
    }
  }

  if (external) e.external = (*external)(x, grad ? &g : nullptr);

  if (grad) {
    gradient->assign(torsionCount_, 0.0);
    for (std::size_t k = 0; k < parts_.size(); ++k) {
      const auto& tree = *parts_[k].tree;
      const int o = offset_[k];
      for (std::size_t t = 0; t < tree.torsions.size(); ++t) {
        const auto& tt = tree.torsions[t];
        const V& a = x[o + tt.fixedAtom];
        V u{x[o + tt.movingAtom][0] - a[0], x[o + tt.movingAtom][1] - a[1], x[o + tt.movingAtom][2] - a[2]};
        const T len = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (len == 0) continue;
        for (T& c : u) c /= len;
        double s = 0;
        for (int idx : tt.moving) {
          const V& p = x[o + idx];
          const V v{p[0] - a[0], p[1] - a[1], p[2] - a[2]};
          const Vec3& gi = g[o + idx];
          // g . (u x v)
          s += gi[0] * (u[1] * v[2] - u[2] * v[1]) + gi[1] * (u[2] * v[0] - u[0] * v[2]) + gi[2] * (u[0] * v[1] - u[1] * v[0]);
        }
        (*gradient)[torsionOffset_[k] + t] = s;
      }
    }
  }
  if (coords) *coords = std::move(x);
  return e;
}

EnergyTerms EnergySystem::evaluate(std::span<const double> anglesDeg, std::vector<double>* gradient,
                                   std::vector<Vec3>* coords, const ExternalTerm* external) const {
  return evaluateT(anglesDeg, gradient, coords, external);
}

EnergyTerms EnergySystem::evaluate(std::span<const float> anglesDeg, std::vector<double>* gradient,
                                   std::vector<Vec3f>* coords, const ExternalTermF* external) const {
  return evaluateT(anglesDeg, gradient, coords, external);
}

}  // namespace od
