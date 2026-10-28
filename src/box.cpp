// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "box.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "closure.hpp"
#include "csgm1.hpp"
#include "energy.hpp"
#include "log.hpp"
#include "ring_ops.hpp"
#include "rng.hpp"
#include "torsion_tree.hpp"

namespace od {
namespace {

using Axes = std::array<Vec3, 3>;

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

double dist(const Vec3& a, const Vec3& b) {
  const double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return std::sqrt(x * x + y * y + z * z);
}

Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Rows of e times rows of a: the box axes a turned by e.
Axes mul(const Axes& e, const Axes& a) {
  Axes r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) r[i][j] = e[i][0] * a[0][j] + e[i][1] * a[1][j] + e[i][2] * a[2][j];
  return r;
}

// Rz(g) Ry(b) Rx(a), angles in degrees.
Axes turn(double a, double b, double g) {
  const double k = std::numbers::pi / 180;
  const double ca = std::cos(a * k), sa = std::sin(a * k), cb = std::cos(b * k), sb = std::sin(b * k);
  const double cg = std::cos(g * k), sg = std::sin(g * k);
  return {{{cg * cb, cg * sb * sa - sg * ca, cg * sb * ca + sg * sa},
           {sg * cb, sg * sb * sa + cg * ca, sg * sb * ca - cg * sa},
           {-sb, cb * sa, cb * ca}}};
}

struct Extent {
  Vec3 lo{}, hi{};
};

Extent extentOf(std::span<const Vec3> x, const Axes& a) {
  Extent e;
  for (int k = 0; k < 3; ++k) {
    e.lo[k] = std::numeric_limits<double>::infinity();
    e.hi[k] = -std::numeric_limits<double>::infinity();
  }
  for (const Vec3& p : x)
    for (int k = 0; k < 3; ++k) {
      const double v = dot(a[k], p);
      e.lo[k] = std::min(e.lo[k], v);
      e.hi[k] = std::max(e.hi[k], v);
    }
  return e;
}

// Per axis: the extent, or with a center c (input coordinates) twice the largest extent from it.
Vec3 spanOf(const Extent& e, const Axes& a, const std::optional<Vec3>& c) {
  Vec3 s{};
  for (int k = 0; k < 3; ++k) {
    if (c) {
      const double v = dot(a[k], *c);
      s[k] = 2 * std::max(e.hi[k] - v, v - e.lo[k]);
    } else {
      s[k] = e.hi[k] - e.lo[k];
    }
  }
  return s;
}

double volumeOf(std::span<const Vec3> x, const Axes& a, double add, const std::optional<Vec3>& c) {
  const Vec3 s = spanOf(extentOf(x, a), a, c);
  return (s[0] + add) * (s[1] + add) * (s[2] + add);
}

}  // namespace

Vec3 BoxFrame::toBox(const Vec3& x) const { return {dot(axes[0], x), dot(axes[1], x), dot(axes[2], x)}; }

Vec3 BoxFrame::toInput(const Vec3& x) const {
  Vec3 r{};
  for (int j = 0; j < 3; ++j) r[j] = axes[0][j] * x[0] + axes[1][j] * x[1] + axes[2][j] * x[2];
  return r;
}

std::vector<Vec3> heavyCoords(const Molecule& m) {
  std::vector<Vec3> v;
  for (const Atom& a : m.atoms)
    if (a.atomicNum > 1) v.push_back(a.xyz);
  return v;
}

double radiusOf(std::span<const Vec3> x) {
  if (x.empty()) return 0;
  Vec3 c{};
  for (const Vec3& p : x)
    for (int k = 0; k < 3; ++k) c[k] += p[k];
  for (int k = 0; k < 3; ++k) c[k] /= static_cast<double>(x.size());
  double r = 0;
  for (const Vec3& p : x) r = std::max(r, dist(p, c));
  return r;
}

double mobmRadiusMax(const Molecule& m, const TorsionTree& tree, const RingOps* rings, int* overflow,
                     long long maxCombos) {
  const int na = static_cast<int>(m.atoms.size());
  // Rigid part of an atom: 1 + the innermost tree torsion moving it, 0 for the root part.
  std::vector<int> part(na, 0);
  for (std::size_t t = 0; t < tree.torsions.size(); ++t)
    for (int a : tree.torsions[t].moving)
      if (a < na) part[a] = static_cast<int>(t) + 1;
  const int np = static_cast<int>(tree.torsions.size()) + 1;
  // Heavy atoms rigid with each part: its own and the bond atoms of the torsions between it and its neighbours.
  std::vector<int> heavy;
  std::vector<int> local(na, -1);
  for (int a = 0; a < na; ++a)
    if (m.atoms[a].atomicNum > 1) {
      local[a] = static_cast<int>(heavy.size());
      heavy.push_back(a);
    }
  const int n = static_cast<int>(heavy.size());
  if (n == 0) return 0;
  std::vector<std::vector<int>> rigid(np);
  for (int i = 0; i < n; ++i) rigid[part[heavy[i]]].push_back(i);
  struct Edge {
    int to;
    std::array<int, 2> ports;
  };
  std::vector<std::vector<Edge>> adj(np);
  for (std::size_t t = 0; t < tree.torsions.size(); ++t) {
    const TreeTorsion& tt = tree.torsions[t];
    const int c = static_cast<int>(t) + 1, p = part[tt.fixedAtom];
    const std::array<int, 2> ports{local[tt.fixedAtom], local[tt.movingAtom]};
    if (ports[0] < 0 || ports[1] < 0) fail("box.mobmRadius = max: torsion without heavy bond atoms");
    for (int q : {c, p})
      for (int x : ports)
        if (std::find(rigid[q].begin(), rigid[q].end(), x) == rigid[q].end()) rigid[q].push_back(x);
    adj[c].push_back({p, ports});
    adj[p].push_back({c, ports});
  }
  std::vector<Vec3> x(n);
  for (int i = 0; i < n; ++i) x[i] = m.atoms[heavy[i]].xyz;
  // Distances in each rigid part (its local order): fixed, or the largest over the ring conformers.
  std::vector<std::vector<double>> dmat(np);
  std::vector<std::vector<int>> where(n);  // per heavy atom: (part, slot) pairs
  for (int p = 0; p < np; ++p) {
    const auto& R = rigid[p];
    const std::size_t k = R.size();
    dmat[p].assign(k * k, 0.0);
    for (std::size_t a = 0; a < k; ++a) {
      where[R[a]].push_back(p);
      where[R[a]].push_back(static_cast<int>(a));
      for (std::size_t b = 0; b < k; ++b) dmat[p][a * k + b] = dist(x[R[a]], x[R[b]]);
    }
  }
  if (overflow) *overflow = 0;
  if (rings && !rings->ops.empty()) {
    std::size_t nc = static_cast<std::size_t>(na);
    for (const RingOp& op : rings->ops) {
      for (int c : op.move) nc = std::max(nc, static_cast<std::size_t>(c) + 1);
      for (int c : op.exitAtom) nc = std::max(nc, static_cast<std::size_t>(c) + 1);
    }
    std::vector<Vec3> base(nc, Vec3{});
    for (int a = 0; a < na; ++a) base[a] = m.atoms[a].xyz;
    const int nr = static_cast<int>(rings->ops.size());
    auto tokenOf = [&](const RingOp& op, int atom) {
      return static_cast<std::size_t>(atom) < op.token.size() ? op.token[atom] : 0;
    };
#pragma omp parallel for schedule(dynamic, 1)
    for (int p = 0; p < np; ++p) {
      const auto& R = rigid[p];
      const std::size_t k = R.size();
      if (k < 2) continue;
      // Ring operations that change distances in this part: its atoms are not all one rigid piece.
      std::vector<int> ops;
      for (int r = 0; r < nr; ++r) {
        const RingOp& op = rings->ops[r];
        const int t0 = tokenOf(op, heavy[R[0]]);
        for (std::size_t a = 1; a < k; ++a)
          if (tokenOf(op, heavy[R[a]]) != t0) {
            ops.push_back(r);
            break;
          }
      }
      if (ops.empty()) continue;
      long long combos = 1;
      bool over = false;
      for (int r : ops) {
        combos *= rings->ops[r].confCount;
        if (combos > maxCombos) {
          over = true;
          break;
        }
      }
      std::vector<double> sel(nr, 0.0);
      std::vector<Vec3> xyz;
      auto take = [&] {
        xyz = base;
        applyRingOps(*rings, std::span<const double>(sel), xyz);
        for (std::size_t a = 0; a < k; ++a)
          for (std::size_t b = a + 1; b < k; ++b) {
            const double d = dist(xyz[heavy[R[a]]], xyz[heavy[R[b]]]);
            double& v = dmat[p][a * k + b];
            if (d > v) v = dmat[p][b * k + a] = d;
          }
      };
      if (!over) {
        for (long long c = 1; c < combos; ++c) {
          long long v = c;
          for (int r : ops) {
            sel[r] = static_cast<double>(v % rings->ops[r].confCount);
            v /= rings->ops[r].confCount;
          }
          take();
        }
      } else {
        for (int r : ops) {
          for (int c = 1; c < rings->ops[r].confCount; ++c) {
            sel[r] = c;
            take();
          }
          sel[r] = 0;
        }
#pragma omp critical
        if (overflow) ++*overflow;
      }
    }
  }
  auto partDist = [&](int p, int i, int j) {
    const auto& w = where[i];
    int si = -1;
    for (std::size_t q = 0; q < w.size(); q += 2)
      if (w[q] == p) si = w[q + 1];
    const auto& v = where[j];
    int sj = -1;
    for (std::size_t q = 0; q < v.size(); q += 2)
      if (v[q] == p) sj = v[q + 1];
    return dmat[p][static_cast<std::size_t>(si) * rigid[p].size() + static_cast<std::size_t>(sj)];
  };
  double best = 0;
#pragma omp parallel
  {
    std::vector<double> d(n);
    std::vector<int> stack;
    std::vector<char> seen(np);
    double mine = 0;
#pragma omp for schedule(dynamic, 16)
    for (int i = 0; i < n; ++i) {
      std::fill(d.begin(), d.end(), std::numeric_limits<double>::infinity());
      std::fill(seen.begin(), seen.end(), 0);
      const int p0 = part[heavy[i]];
      for (int j : rigid[p0]) d[j] = partDist(p0, i, j);
      seen[p0] = 1;
      stack.assign(1, p0);
      while (!stack.empty()) {
        const int p = stack.back();
        stack.pop_back();
        for (const Edge& e : adj[p]) {
          if (seen[e.to]) continue;
          seen[e.to] = 1;
          for (int j : rigid[e.to])
            for (int q : e.ports) d[j] = std::min(d[j], d[q] + partDist(e.to, q, j));
          stack.push_back(e.to);
        }
      }
      double s = 0;
      for (int j = 0; j < n; ++j) s += d[j];
      mine = std::max(mine, s / n);
    }
#pragma omp critical
    best = std::max(best, mine);
  }
  return best;
}

SampledRadius mobmRadiusSampled(const Molecule& m, const EnergySystem& sys, const Csgm1GroupSampler& sampler,
                                int sampleNum, std::uint64_t seed, std::uint64_t stream, const LoopCloser* closer,
                                double maxDeviation) {
  std::vector<int> heavy;
  for (int a = 0; a < static_cast<int>(m.atoms.size()); ++a)
    if (m.atoms[a].atomicNum > 1) heavy.push_back(a);
  SampledRadius out;
  {
    const std::vector<Vec3> h = heavyCoords(m);
    out.radius = radiusOf(h);
  }
  if (sys.torsionCount() == 0) return out;
  const int off = sys.coordOffset(0);
  std::vector<double> r(sampleNum, -1.0);
  std::vector<char> notClosed(sampleNum, 0);
  const bool closing = closer && closer->active();
  const SegmentSampler resample = [&](Rng& g, int seg, std::span<double> a) { sampler.sampleSegment(g, seg, a); };
#pragma omp parallel for schedule(dynamic, 4)
  for (int k = 0; k < sampleNum; ++k) {
    Rng rng(seed, stream, static_cast<std::uint64_t>(k));
    std::vector<double> a = sys.inputAngles();
    sampler.sample(rng, a);
    if (sys.ringSlotCount() > 0) {
      Rng rr(seed, (1ULL << 61) + stream, static_cast<std::uint64_t>(k));
      const int t0 = sys.treeTorsionCount();
      for (int s = 0; s < sys.ringSlotCount(); ++s)
        a[t0 + s] = static_cast<double>(rr.below(static_cast<std::uint64_t>(sys.ringConfCounts()[s])));
    }
    if (closing && (!closer->close(rng, a, resample) || (maxDeviation > 0 && sys.closureDeviation(a) > maxDeviation))) {
      notClosed[k] = 1;
      continue;
    }
    std::vector<Vec3> c;
    if (sys.evaluate(a, nullptr, &c).clashPairs > 0) continue;
    std::vector<Vec3> h;
    h.reserve(heavy.size());
    for (int i : heavy) h.push_back(c[off + i]);
    r[k] = radiusOf(h);
  }
  for (char v : notClosed) out.notClosed += v;
  for (double v : r)
    if (v >= 0) {
      ++out.accepted;
      out.radius = std::max(out.radius, v);
    }
  return out;
}

std::array<Vec3, 3> principalAxes(std::span<const Vec3> x) {
  Axes a{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (x.size() < 2) return a;
  Eigen::Vector3d c = Eigen::Vector3d::Zero();
  for (const Vec3& p : x) c += Eigen::Vector3d(p[0], p[1], p[2]);
  c /= static_cast<double>(x.size());
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const Vec3& p : x) {
    const Eigen::Vector3d d = Eigen::Vector3d(p[0], p[1], p[2]) - c;
    cov += d * d.transpose();
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  for (int k = 0; k < 2; ++k) {
    const Eigen::Vector3d v = es.eigenvectors().col(2 - k);
    int big = 0;
    for (int j = 1; j < 3; ++j)
      if (std::abs(v[j]) > std::abs(v[big])) big = j;
    const double s = v[big] < 0 ? -1.0 : 1.0;
    a[k] = {s * v[0], s * v[1], s * v[2]};
  }
  a[2] = cross(a[0], a[1]);
  return a;
}

BoxFrame autoBox(std::span<const Vec3> stamHeavy, const AutoBoxParams& p) {
  if (stamHeavy.empty()) fail("automatic box: the StaM group has no heavy atoms");
  const double add = p.stamMargin + 2 * p.mobmRadius;
  BoxFrame b;
  if (p.autoRotate) {
    // Three passes over a 7 x 7 x 7 grid of turns about the box axes, each around the best so far.
    Axes best = principalAxes(stamHeavy);
    double vBest = volumeOf(stamHeavy, best, add, p.center);
    for (double step : {10.0, 2.0, 0.4}) {
      const Axes base = best;
      for (int i = -3; i <= 3; ++i)
        for (int j = -3; j <= 3; ++j)
          for (int k = -3; k <= 3; ++k) {
            if (i == 0 && j == 0 && k == 0) continue;
            const Axes a = mul(turn(i * step, j * step, k * step), base);
            const double v = volumeOf(stamHeavy, a, add, p.center);
            if (v < vBest * (1 - 1e-12)) {
              vBest = v;
              best = a;
            }
          }
    }
    b.axes = best;
    b.rotated = true;
  }
  const Extent e = extentOf(stamHeavy, b.axes);
  const Vec3 s = spanOf(e, b.axes, p.center);
  const Vec3 c = p.center ? b.toBox(*p.center) : Vec3{};
  for (int k = 0; k < 3; ++k) {
    b.center[k] = p.center ? c[k] : (e.lo[k] + e.hi[k]) / 2;
    b.size[k] = p.size ? (*p.size)[k] : s[k] + add;
  }
  return b;
}

}  // namespace od
