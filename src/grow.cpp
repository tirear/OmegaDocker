// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "grow.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <set>

#include "log.hpp"
#include "torsion_tree.hpp"

namespace od {

namespace {

using Mat3 = std::array<double, 9>;

Mat3 mul(const Mat3& a, const Mat3& b) {
  Mat3 r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) r[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
  return r;
}

Vec3 matVec(const Mat3& m, const Vec3& v) {
  return {m[0] * v[0] + m[1] * v[1] + m[2] * v[2], m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
          m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
}

// Rotation by angle (radians) about unit axis k.
Mat3 rotation(const Vec3& k, double angle) {
  const double c = std::cos(angle), s = std::sin(angle), t = 1 - c;
  return {t * k[0] * k[0] + c,        t * k[0] * k[1] - s * k[2], t * k[0] * k[2] + s * k[1],
          t * k[0] * k[1] + s * k[2], t * k[1] * k[1] + c,        t * k[1] * k[2] - s * k[0],
          t * k[0] * k[2] - s * k[1], t * k[1] * k[2] + s * k[0], t * k[2] * k[2] + c};
}

double dist2(const Vec3& a, const Vec3& b) {
  const double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return x * x + y * y + z * z;
}

}  // namespace

GrowBuilder::GrowBuilder(const EnergySystem& sys, const std::vector<EnergyPart>& parts,
                         const std::vector<const Topology*>& topologies, const std::vector<std::vector<int>>& segments,
                         const GrowParams& p)
    : clash2_(p.clashDistance * p.clashDistance), maxRetry_(p.maxRetry), maxBacktrack_(p.maxBacktrack) {
  // Input coordinates and global torsions.
  std::vector<int> torOffset;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const Molecule& m = *parts[k].mol;
    const TorsionTree& tree = *parts[k].tree;
    const int co = sys.coordOffset(static_cast<int>(k)), to = static_cast<int>(tor_.size());
    torOffset.push_back(to);
    input_.resize(co + m.atoms.size() + tree.virtuals.size());
    for (std::size_t a = 0; a < m.atoms.size(); ++a) input_[co + a] = m.atoms[a].xyz;
    for (std::size_t v = 0; v < tree.virtuals.size(); ++v) input_[co + m.atoms.size() + v] = tree.virtuals[v].xyz;
    for (const TreeTorsion& tt : tree.torsions) {
      Tor t;
      t.parent = tt.parent < 0 ? -1 : tt.parent + to;
      t.o = input_[co + tt.fixedAtom];
      const Vec3 q = input_[co + tt.movingAtom];
      Vec3 kv{q[0] - t.o[0], q[1] - t.o[1], q[2] - t.o[2]};
      const double len = std::sqrt(kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2]);
      t.k = len > 0 ? Vec3{kv[0] / len, kv[1] / len, kv[2] / len} : Vec3{0, 0, 0};
      t.inputAngle = tt.inputAngle;
      tor_.push_back(std::move(t));
    }
  }
  const int T = static_cast<int>(tor_.size());
  if (T != sys.torsionCount()) fail("grow: torsion count mismatch");
  if (!p.fixed.empty() && static_cast<int>(p.fixed.size()) != T) fail("grow: fixed torsion mask size mismatch");

  // Deepest torsion moving each coordinate (parents precede children). The
  // moving bond atom lies on the axis and is placed with the parent torsion.
  std::vector<int> owner(input_.size(), -1);
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const int co = sys.coordOffset(static_cast<int>(k));
    const auto& tt = parts[k].tree->torsions;
    for (std::size_t i = 0; i < tt.size(); ++i) {
      for (int c : tt[i].moving) owner[co + c] = torOffset[k] + static_cast<int>(i);
      owner[co + tt[i].movingAtom] = tor_[torOffset[k] + i].parent;
    }
  }
  for (std::size_t c = 0; c < owner.size(); ++c)
    if (owner[c] >= 0) tor_[owner[c]].own.push_back(static_cast<int>(c));

  // Steps: segments in order of first index, each once all outside parents are built.
  std::vector<int> segOf(T, -1);
  for (std::size_t s = 0; s < segments.size(); ++s)
    for (int i : segments[s]) {
      if (i < 0 || i >= T || segOf[i] >= 0) fail("grow: invalid torsion segments");
      segOf[i] = static_cast<int>(s);
    }
  if (std::count(segOf.begin(), segOf.end(), -1) > 0) fail("grow: invalid torsion segments");
  std::vector<int> stepOf(T, -1);
  std::vector<bool> done(segments.size(), false);
  for (std::size_t n = 0; n < segments.size(); ++n) {
    std::size_t s = 0;
    for (; s < segments.size(); ++s) {
      if (done[s]) continue;
      bool ready = true;
      for (int i : segments[s]) {
        const int par = tor_[i].parent;
        ready = ready && (par < 0 || stepOf[par] >= 0 || segOf[par] == static_cast<int>(s));
      }
      if (ready) break;
    }
    if (s == segments.size()) fail("grow: the torsion segments cannot be built from the root");
    done[s] = true;
    Step st;
    st.segment = static_cast<int>(s);
    st.torsions = segments[s];
    std::sort(st.torsions.begin(), st.torsions.end());
    for (int i : st.torsions) {
      if (!p.fixed.empty() && p.fixed[i]) st.fixed = true;
      const int par = tor_[i].parent;
      if (par >= 0 && segOf[par] != static_cast<int>(s)) st.parentStep = std::max(st.parentStep, stepOf[par]);
    }
    for (int i : st.torsions) stepOf[i] = static_cast<int>(steps_.size());
    steps_.push_back(std::move(st));
  }
  auto stepOfCoord = [&](int c) { return owner[c] < 0 ? -1 : stepOf[owner[c]]; };

  // Clash pairs (heavy atoms, no metal), checked at the step that places the later atom.
  for (const auto& pr : sys.pairs()) {
    if (!pr.clash) continue;
    const int s = std::max(stepOfCoord(pr.i), stepOfCoord(pr.j));
    if (s >= 0 && !steps_[s].fixed) steps_[s].pairs.emplace_back(pr.i, pr.j);
  }

  // Reach (spec 4.8): broken backbone links. For each end of the link and
  // each torsion on the tree path to that end, the moving atom of the torsion
  // must be within (residues from its residue to the other end, inclusive) x
  // reachPerRes of the other end, once the other end is built.
  for (std::size_t k = 0; k < parts.size() && p.reach; ++k) {
    const Molecule& m = *parts[k].mol;
    const Topology& t = *topologies[k];
    const int co = sys.coordOffset(static_cast<int>(k));
    const double perRes = t.backboneType == BackboneType::Protein   ? p.reachProtein
                          : t.backboneType == BackboneType::Nucleic ? p.reachNucleic
                                                                    : 0.0;
    if (perRes <= 0) continue;
    // Chain and position of the unit of each residue.
    std::vector<int> chainOf(m.residues.size(), -1), posOf(m.residues.size(), 0);
    for (std::size_t u = 0; u < t.units.size(); ++u) {
      if (t.units[u].prev >= 0) continue;
      int pos = 0;
      for (int v = static_cast<int>(u); v >= 0; v = t.units[v].next, ++pos)
        if (const int r = t.units[v].residue; r >= 0) {
          chainOf[r] = static_cast<int>(u);
          posOf[r] = pos;
        }
    }
    const std::set<int> broken(t.brokenBonds.begin(), t.brokenBonds.end());
    for (std::size_t u = 0; u < t.units.size(); ++u) {
      if (t.units[u].next < 0) continue;
      const int bi = backboneLink(m, t, static_cast<int>(u), t.units[u].next);
      if (bi < 0 || !broken.count(bi)) continue;
      const Bond& b = m.bonds[bi];
      for (auto [end, target] : {std::pair{b.a, b.b}, std::pair{b.b, b.a}}) {
        const int targetRes = m.atoms[target].residue;
        if (targetRes < 0 || chainOf[targetRes] < 0) continue;
        const int targetStep = stepOfCoord(co + target);
        std::set<int> seen;
        for (int ti = owner[co + end]; ti >= 0; ti = tor_[ti].parent) {
          const int s = stepOf[ti];
          if (!seen.insert(s).second || targetStep > s || steps_[s].fixed) continue;
          const int local = ti - torOffset[k];
          const int moving = parts[k].tree->torsions[local].movingAtom;
          const int r = m.atoms[moving].residue;
          if (r < 0 || chainOf[r] != chainOf[targetRes]) continue;
          const double bound = (std::abs(posOf[r] - posOf[targetRes]) + 1) * perRes;
          steps_[s].reach.push_back({co + moving, co + target, bound * bound});
        }
      }
    }
  }
}

int GrowBuilder::fixedStepCount() const {
  return static_cast<int>(std::count_if(steps_.begin(), steps_.end(), [](const Step& s) { return s.fixed; }));
}

std::size_t GrowBuilder::reachCheckCount() const {
  std::size_t n = 0;
  for (const auto& s : steps_) n += s.reach.size();
  return n;
}

GrowResult GrowBuilder::build(Rng& rng, std::span<double> a, const SegmentSampler& resample,
                              std::vector<Vec3>* coords) const {
  constexpr double kRad = std::numbers::pi / 180.0;
  GrowResult res;
  std::vector<Vec3> x(input_);
  std::vector<Mat3> R(tor_.size());
  std::vector<Vec3> tr(tor_.size());
  std::vector<int> backtracks(steps_.size(), 0);
  for (std::size_t si = 0; si < steps_.size();) {
    const Step& st = steps_[si];
    bool ok = true;
    for (int attempt = 0;; ++attempt) {
      for (int i : st.torsions) {
        const Tor& t = tor_[i];
        const Mat3 rk = rotation(t.k, (a[i] - t.inputAngle) * kRad);
        const Vec3 ro = matVec(rk, t.o);
        const Vec3 d{t.o[0] - ro[0], t.o[1] - ro[1], t.o[2] - ro[2]};
        if (t.parent < 0) {
          R[i] = rk;
          tr[i] = d;
        } else {
          R[i] = mul(R[t.parent], rk);
          const Vec3 pd = matVec(R[t.parent], d);
          const Vec3& pt = tr[t.parent];
          tr[i] = {pd[0] + pt[0], pd[1] + pt[1], pd[2] + pt[2]};
        }
        for (int c : t.own) {
          const Vec3 v = matVec(R[i], input_[c]);
          x[c] = {v[0] + tr[i][0], v[1] + tr[i][1], v[2] + tr[i][2]};
        }
      }
      ok = true;
      for (const auto& [i, j] : st.pairs)
        if (dist2(x[i], x[j]) < clash2_) {
          ok = false;
          ++res.clashFails;
          break;
        }
      if (ok)
        for (const Reach& r : st.reach)
          if (dist2(x[r.current], x[r.target]) > r.bound2) {
            ok = false;
            ++res.reachFails;
            break;
          }
      if (ok || attempt >= maxRetry_ || !resample) break;
      resample(rng, st.segment, a);
      ++res.retries;
    }
    if (ok) {
      ++si;
      continue;
    }
    // Backtrack one level: resample the parent segment and rebuild from it.
    const int ps = st.parentStep;
    if (!resample || ps < 0 || steps_[ps].fixed || backtracks[si] >= maxBacktrack_) return res;
    ++backtracks[si];
    ++res.backtracks;
    resample(rng, steps_[ps].segment, a);
    si = static_cast<std::size_t>(ps);
  }
  res.ok = true;
  if (coords) *coords = std::move(x);
  return res;
}

}  // namespace od
