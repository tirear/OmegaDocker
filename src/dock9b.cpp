// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "dock9b.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <type_traits>

#include "log.hpp"
#include "rng.hpp"

namespace od {

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

Quat conj(const Quat& q) { return {q.w, -q.x, -q.y, -q.z}; }
QuatF conj(const QuatF& q) { return {q.w, -q.x, -q.y, -q.z}; }

template <class T>
Vec3T<T> castVec(const Vec3& v) {
  return {static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2])};
}

std::uint64_t splitmix(std::uint64_t z) {
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

}  // namespace

Dock9bStats& Dock9bStats::operator+=(const Dock9bStats& o) {
  l1 += o.l1;
  overlap += o.overlap;
  far += o.far;
  outsideBox += o.outsideBox;
  l2 += o.l2;
  leftGrid += o.leftGrid;
  steps += o.steps;
  evaluations += o.evaluations;
  return *this;
}

std::uint64_t dockStream(std::uint64_t round, std::uint64_t mobm, std::uint64_t stamSample, std::uint64_t mobmSample) {
  std::uint64_t s = splitmix(0x9b9b9b9bULL ^ round);
  s = splitmix(s ^ mobm);
  s = splitmix(s ^ stamSample);
  return splitmix(s ^ mobmSample);
}

MobmDocker::MobmDocker(const EnergySystem& sys, const EnergyGrids& grids, const DistanceGrid& dist,
                       const Dock9bParams& p)
    : sys_(sys), grids_(grids), dist_(dist), p_(p) {
  const Ad4Typing& t = sys.typing(0);
  const int off = sys.coordOffset(0);
  for (std::size_t a = 0; a < t.type.size(); ++a) {
    if (t.type[a] == Ad4Type::Merged) continue;
    const int m = grids.mapOf(t.type[a]);
    if (m < 0) fail("no grid map for AD4 type {}", ad4TypeName(t.type[a]));
    probe_.push_back(off + static_cast<int>(a));
    map_.push_back(m);
    charge_.push_back(t.charge[a]);
    heavy_.push_back(t.type[a] != Ad4Type::H && t.type[a] != Ad4Type::HD && t.type[a] != Ad4Type::HS);
  }
  for (int c = 0; c < 3; ++c) posN_[c] = static_cast<int>(std::floor(p.boxSize[c] / p.posStep + 1e-9)) + 1;
  if (p.posm == Posm::Posm2) hopf_ = hopfGrid(p.hopfLevel);
}

long long MobmDocker::PosRange::count() const {
  long long n = 1;
  for (int d = 0; d < 3; ++d) n *= std::max(0, hi[d] - lo[d] + 1);
  return n;
}

long long MobmDocker::l1Count(const Conf& c) const {
  if (p_.posm == Posm::Posm1) return p_.posm1SampleNum;
  return c.hopfOffset.back();
}

MobmDocker::PosRange MobmDocker::posRange(const Conf& c, const Quat& q) const {
  const QuatF qf = toFloat(q);
  const Vec3f cf = castVec<float>(c.center);
  Vec3f lo{1e30f, 1e30f, 1e30f}, hi{-1e30f, -1e30f, -1e30f};
  for (int k : probe_) {
    const Vec3f& y = c.coords[k];
    const Vec3f r = rotate(qf, {y[0] - cf[0], y[1] - cf[1], y[2] - cf[2]});
    for (int d = 0; d < 3; ++d) {
      lo[d] = std::min(lo[d], r[d]);
      hi[d] = std::max(hi[d], r[d]);
    }
  }
  // Centre t = boxMin + i * posStep with boxMin + margin <= t + lo and t + hi <= boxMax - margin.
  PosRange pr;
  for (int d = 0; d < 3; ++d) {
    pr.lo[d] = std::max(0, static_cast<int>(std::ceil((kPosMargin - lo[d]) / p_.posStep - 1e-9)));
    pr.hi[d] = std::min(posN_[d] - 1,
                        static_cast<int>(std::floor((p_.boxSize[d] - kPosMargin - hi[d]) / p_.posStep + 1e-9)));
  }
  return pr;
}

MobmDocker::Conf MobmDocker::conformation(std::span<const double> anglesDeg) const {
  Conf c;
  c.anglesDeg.assign(anglesDeg.begin(), anglesDeg.end());
  const std::vector<float> af(anglesDeg.begin(), anglesDeg.end());
  sys_.evaluate(std::span<const float>(af), nullptr, &c.coords, nullptr);
  int n = 0;
  for (std::size_t k = 0; k < probe_.size(); ++k)
    if (heavy_[k]) {
      for (int d = 0; d < 3; ++d) c.center[d] += c.coords[probe_[k]][d];
      ++n;
    }
  for (double& v : c.center) v /= n;
  if (p_.posm == Posm::Posm2) {
    c.hopfRange.resize(hopf_.size());
    c.hopfOffset.assign(hopf_.size() + 1, 0);
#pragma omp parallel for schedule(static)
    for (std::size_t j = 0; j < hopf_.size(); ++j) c.hopfRange[j] = posRange(c, hopf_[j]);
    for (std::size_t j = 0; j < hopf_.size(); ++j) c.hopfOffset[j + 1] = c.hopfOffset[j] + c.hopfRange[j].count();
  }
  return c;
}

std::optional<Pose> MobmDocker::l1Pose(const Conf& conf, std::uint64_t stream, long long k) const {
  Pose pose;
  std::array<long long, 3> idx{};
  if (p_.posm == Posm::Posm1) {
    Rng rng(p_.seed, stream, static_cast<std::uint64_t>(k));
    pose.q = randomQuat(rng);
    const PosRange pr = posRange(conf, pose.q);
    if (pr.count() == 0) return std::nullopt;
    for (int c = 0; c < 3; ++c)
      idx[c] = pr.lo[c] + static_cast<long long>(rng.below(static_cast<std::uint64_t>(pr.hi[c] - pr.lo[c] + 1)));
  } else {
    // Orientation-major order: orientation j covers [hopfOffset[j], hopfOffset[j + 1]).
    const auto it = std::upper_bound(conf.hopfOffset.begin(), conf.hopfOffset.end(), k);
    const std::size_t j = static_cast<std::size_t>(it - conf.hopfOffset.begin()) - 1;
    pose.q = hopf_[j];
    const PosRange& pr = conf.hopfRange[j];
    long long pos = k - conf.hopfOffset[j];
    for (int c = 0; c < 3; ++c) {
      const long long w = pr.hi[c] - pr.lo[c] + 1;
      idx[c] = pr.lo[c] + pos % w;
      pos /= w;
    }
  }
  for (int c = 0; c < 3; ++c)
    pose.t[c] = p_.boxCenter[c] - 0.5 * p_.boxSize[c] + p_.posStep * static_cast<double>(idx[c]);
  return pose;
}

bool MobmDocker::quickFilter(const Conf& c, const Pose& pose, Dock9bStats& stats) const {
  int overlap = 0;
  bool near = false;
  const QuatF qf = toFloat(pose.q);
  const Vec3f cf = castVec<float>(c.center), tf = castVec<float>(pose.t);
  Vec3f boxLo, boxHi;
  for (int d = 0; d < 3; ++d) {
    boxLo[d] = static_cast<float>(p_.boxCenter[d] - 0.5 * p_.boxSize[d]);
    boxHi[d] = static_cast<float>(p_.boxCenter[d] + 0.5 * p_.boxSize[d]);
  }
  const float cd = static_cast<float>(p_.clashDistance), nearD = static_cast<float>(p_.maxNearDistance);
  for (std::size_t k = 0; k < probe_.size(); ++k) {
    const Vec3f& y = c.coords[probe_[k]];
    const Vec3f r = rotate(qf, {y[0] - cf[0], y[1] - cf[1], y[2] - cf[2]});
    const Vec3f w{r[0] + tf[0], r[1] + tf[1], r[2] + tf[2]};
    for (int d = 0; d < 3; ++d)
      if (w[d] < boxLo[d] || w[d] > boxHi[d]) {
        ++stats.outsideBox;
        return false;
      }
    if (!heavy_[k]) continue;
    const float dd = interpolate(dist_.geom, dist_.d, w);
    if (dd < cd) ++overlap;
    if (dd <= nearD) near = true;
  }
  if (overlap >= p_.maxOverlapAtomNum) {
    ++stats.overlap;
    return false;
  }
  if (!near) {
    ++stats.far;
    return false;
  }
  return true;
}

template <class T>
MobmDocker::PoseEnergy MobmDocker::evaluateT(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                                             std::vector<double>* grad) const {
  using V = Vec3T<T>;
  PoseEnergy pe;
  Vec3 gt{}, gr{};  // double accumulators
  const QuatF qf = toFloat(pose.q);
  auto rot = [&](const V& v) -> V {
    if constexpr (std::is_same_v<T, float>) return rotate(qf, v);
    else return rotate(pose.q, v);
  };
  auto rotBack = [&](const V& v) -> V {
    if constexpr (std::is_same_v<T, float>) return rotate(conj(qf), v);
    else return rotate(conj(pose.q), v);
  };
  const V ctr = castVec<T>(c.center), tt = castVec<T>(pose.t);
  const V gLo = castVec<T>(grids_.geom.origin), gHi = castVec<T>(grids_.geom.max());
  const T cd = static_cast<T>(p_.clashDistance), clashBase = static_cast<T>(p_.clashBase),
          clashK = static_cast<T>(p_.clashK);
  const ExternalTermT<T> ext = [&](const std::vector<V>& x, std::vector<Vec3>* g) {
    for (std::size_t k = 0; k < probe_.size(); ++k) {
      const V& y = x[probe_[k]];
      const V r = rot({y[0] - ctr[0], y[1] - ctr[1], y[2] - ctr[2]});
      const V w{r[0] + tt[0], r[1] + tt[1], r[2] + tt[2]};
      bool in = true;
      for (int d = 0; d < 3; ++d) in &= w[d] >= gLo[d] && w[d] <= gHi[d];
      if (!in) {
        ++pe.outside;
        continue;
      }
      const T q = static_cast<T>(charge_[k]), aq = std::abs(q);
      V a{}, b{}, s{}, dg{};
      const T ew = interpolate(grids_.geom, grids_.typeMaps[map_[k]], w, g ? &a : nullptr) +
                   q * interpolate(grids_.geom, grids_.elecMap, w, g ? &b : nullptr) +
                   aq * interpolate(grids_.geom, grids_.desolvMap, w, g ? &s : nullptr);
      pe.inter += ew;
      V gw{a[0] + q * b[0] + aq * s[0], a[1] + q * b[1] + aq * s[1], a[2] + q * b[2] + aq * s[2]};
      if (heavy_[k]) {
        const T d = interpolate(dist_.geom, dist_.d, w, g ? &dg : nullptr);
        if (d < cd) {
          pe.interClash += clashBase + clashK * (cd - d) * (cd - d);
          for (int i = 0; i < 3; ++i) gw[i] += -2 * clashK * (cd - d) * dg[i];
        }
      }
      if (!g) continue;
      for (int i = 0; i < 3; ++i) gt[i] += gw[i];
      gr[0] += r[1] * gw[2] - r[2] * gw[1];
      gr[1] += r[2] * gw[0] - r[0] * gw[2];
      gr[2] += r[0] * gw[1] - r[1] * gw[0];
      const V gy = rotBack(gw);
      for (int i = 0; i < 3; ++i) (*g)[probe_[k]][i] += gy[i];
    }
    return pe.inter + pe.interClash;
  };
  const std::vector<T> ang(anglesDeg.begin(), anglesDeg.end());
  std::vector<double> gTors;
  pe.intra = sys_.evaluate(std::span<const T>(ang), grad ? &gTors : nullptr, nullptr, &ext);
  if (grad) {
    grad->assign(6 + gTors.size(), 0.0);
    for (int i = 0; i < 3; ++i) {
      (*grad)[i] = gt[i];
      (*grad)[3 + i] = gr[i];
    }
    std::copy(gTors.begin(), gTors.end(), grad->begin() + 6);
  }
  return pe;
}

MobmDocker::PoseEnergy MobmDocker::evaluate(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                                            std::vector<double>* grad) const {
  return evaluateT<float>(c, pose, anglesDeg, grad);
}

MobmDocker::PoseEnergy MobmDocker::evaluateDouble(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                                                  std::vector<double>* grad) const {
  return evaluateT<double>(c, pose, anglesDeg, grad);
}

bool MobmDocker::minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats) const {
  const std::size_t nt = c.anglesDeg.size();
  const double sT = p_.scaleTrans, sR = p_.scaleRot, sA = p_.scaleTors;
  Quat q0 = start.q;
  // z = (t / sT, rotation vector / sR, torsions in rad / sA); the rotation is re-based after each step.
  std::vector<double> z(6 + nt);
  for (int i = 0; i < 3; ++i) z[i] = start.t[i] / sT;
  for (std::size_t k = 0; k < nt; ++k) z[6 + k] = c.anglesDeg[k] * kDeg / sA;
  // A point of the energy surface: the float pose and angles the energy is evaluated at.
  struct Point {
    Pose pose;
    std::vector<double> angles;  // float values, wrapped to (-180, 180]
    PoseEnergy pe;
  };
  Point trial, accepted;
  trial.angles.resize(nt);
  bool started = false;
  std::vector<double> g;
  auto unpack = [&](std::span<const double> v, Point& pt) {
    for (int i = 0; i < 3; ++i) pt.pose.t[i] = static_cast<float>(v[i] * sT);
    const Quat q = normalized(fromRotationVector({v[3] * sR, v[4] * sR, v[5] * sR}) * q0);
    const QuatF qf = toFloat(q);
    pt.pose.q = {qf.w, qf.x, qf.y, qf.z};
    for (std::size_t k = 0; k < nt; ++k) pt.angles[k] = static_cast<float>(wrapDeg(v[6 + k] * sA / kDeg));
  };
  const VecObjective f = [&](std::span<const double> v, double& e, std::vector<double>& grad) {
    unpack(v, trial);
    trial.pe = evaluate(c, trial.pose, trial.angles, &g);
    if (trial.pe.outside > 0) return false;
    if (!started) {
      accepted = trial;
      started = true;
    }
    e = trial.pe.objective();
    grad.resize(g.size());
    for (std::size_t i = 0; i < g.size(); ++i) grad[i] = g[i] * (i < 3 ? sT : i < 6 ? sR : sA);
    return true;
  };
  LbfgsParams lp = p_.lbfgs;
  const double maxA = p_.maxAngleStep * kDeg;
  lp.stepLimit = [&](std::span<const double> d) {
    double s = 1;
    const double tn = sT * std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (tn > p_.maxTransStep) s = std::min(s, p_.maxTransStep / tn);
    const double rn = sR * std::sqrt(d[3] * d[3] + d[4] * d[4] + d[5] * d[5]);
    if (rn > maxA) s = std::min(s, maxA / rn);
    for (std::size_t k = 6; k < d.size(); ++k)
      if (sA * std::abs(d[k]) > maxA) s = std::min(s, maxA / (sA * std::abs(d[k])));
    return s;
  };
  // The accepted point is the last one evaluated (the line search stops at it).
  lp.onAccept = [&](std::vector<double>& v) {
    accepted = trial;
    q0 = accepted.pose.q;
    v[3] = v[4] = v[5] = 0;
  };
  const LbfgsResult r = minimizeLbfgs(f, z, lp);
  stats.steps += r.steps;
  stats.evaluations += r.evaluations;
  if (r.aborted) {
    ++stats.leftGrid;
    return false;
  }
  const PoseEnergy& pe = accepted.pe;
  out.pose = accepted.pose;
  out.angles.assign(accepted.angles.begin(), accepted.angles.end());
  out.energy = static_cast<float>(pe.ranking());
  out.inter = static_cast<float>(pe.inter);
  out.interClash = static_cast<float>(pe.interClash);
  out.intra = static_cast<float>(pe.intra.filterEnergy());
  out.closureMaxDeviation = static_cast<float>(pe.intra.closureMaxDeviation);
  out.steps = r.steps;
  out.center = c.center;
  return true;
}

Dock9bStats MobmDocker::dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream,
                                 std::vector<DockSample>& out) const {
  const Conf c = conformation(mobAnglesDeg);
  const long long n = l1Count(c);
  Dock9bStats total;
  total.l1 = n;
  std::vector<long long> l2;
#pragma omp parallel
  {
    Dock9bStats st;
    std::vector<long long> mine;
#pragma omp for schedule(static) nowait
    for (long long k = 0; k < n; ++k) {
      const auto pose = l1Pose(c, stream, k);
      if (!pose) ++st.outsideBox;
      else if (quickFilter(c, *pose, st)) mine.push_back(k);
    }
#pragma omp critical
    {
      total += st;
      l2.insert(l2.end(), mine.begin(), mine.end());
    }
  }
  std::sort(l2.begin(), l2.end());
  total.l2 = static_cast<long long>(l2.size());
  std::vector<DockSample> res(l2.size());
  std::vector<char> ok(l2.size(), 0);
#pragma omp parallel
  {
    Dock9bStats st;
#pragma omp for schedule(dynamic, 1) nowait
    for (std::size_t i = 0; i < l2.size(); ++i) ok[i] = minimize(c, *l1Pose(c, stream, l2[i]), res[i], st);
#pragma omp critical
    {
      total.leftGrid += st.leftGrid;
      total.steps += st.steps;
      total.evaluations += st.evaluations;
    }
  }
  for (std::size_t i = 0; i < res.size(); ++i)
    if (ok[i]) out.push_back(std::move(res[i]));
  return total;
}

}  // namespace od
