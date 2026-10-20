// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "dock9b.hpp"

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <numbers>
#include <type_traits>

#include "error.hpp"
#include "fiber_batch.hpp"
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

double dist2(const Vec3& a, const Vec3& b) {
  const double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return x * x + y * y + z * z;
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
  sepaPairs += o.sepaPairs;
  sepaSingles += o.sepaSingles;
  sepaTriples += o.sepaTriples;
  return *this;
}

void MobmDocker::minimizeBatched(PoseEvaluator& pev, std::span<const MobmDocker* const> dks,
                                 std::span<const BatchItem> items, std::vector<DockSample>& res, std::vector<char>& ok,
                                 std::vector<Dock9bStats>& stats) {
  if (dks.size() == 1) pev.bind(*dks[0]);
  else pev.bindAll(dks);
  const int m = static_cast<int>(items.size());
  const int nf = std::max(1, std::min(pev.maxBatch(), m));
  const std::size_t ng = 6 + static_cast<std::size_t>(dks[0]->sys_.torsionCount());
  std::vector<PoseRequest> req(nf);
  std::vector<PoseEnergy> pe(nf);
  std::vector<double> grad(nf * ng);
  const auto task = [&](int i, int slot) {
    const Conf& c = *items[i].conf;
    const MobmDocker& dk = *items[i].dk;
    const int field = items[i].field;
    const PoseEval ev = [&, slot, field](const Pose& pose, std::span<const double> a, std::vector<double>* g) {
      req[slot] = {pose, c.center, a.data(), field};
      batchWait();
      if (g) {
        const auto b = grad.begin() + static_cast<std::ptrdiff_t>(slot * ng);
        g->assign(b, b + static_cast<std::ptrdiff_t>(ng));
      }
      return pe[slot];
    };
    ok[i] = dk.minimize(c, *dk.l1Pose(c, items[i].stream, items[i].k), res[i], stats[i], ev);
  };
  const int lanes = std::max(1, pev.lanes());
  std::vector<std::vector<PoseRequest>> breq(lanes);
  std::vector<std::vector<PoseEnergy>> bout(lanes);
  std::vector<std::vector<double>> bgrad(lanes);
  std::vector<std::vector<int>> bslots(lanes);
  const auto gather = [&](int l, std::span<const int> slots) {
    const std::size_t n = slots.size();
    breq[l].resize(n);
    bout[l].resize(n);
    bgrad[l].resize(n * ng);
    bslots[l].assign(slots.begin(), slots.end());
    for (std::size_t i = 0; i < n; ++i) breq[l][i] = req[slots[i]];
  };
  const auto scatter = [&](int l) {
    const std::vector<int>& slots = bslots[l];
    const long long n = static_cast<long long>(slots.size());
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < n; ++i) {
      pe[slots[i]] = bout[l][i];
      std::copy_n(bgrad[l].begin() + static_cast<std::ptrdiff_t>(i * ng), ng,
                  grad.begin() + static_cast<std::ptrdiff_t>(slots[i] * ng));
    }
  };
  if (lanes == 1) {
    const auto flush = [&](std::span<const int> slots) {
      gather(0, slots);
      pev.evaluate(static_cast<int>(slots.size()), breq[0].data(), bout[0].data(), bgrad[0].data());
      scatter(0);
    };
    runBatched(m, nf, task, flush);
    return;
  }
  const auto submit = [&](int l, std::span<const int> slots) {
    gather(l, slots);
    pev.submit(l, static_cast<int>(slots.size()), breq[l].data());
  };
  const auto complete = [&](int l) {
    pev.complete(l, bout[l].data(), bgrad[l].data());
    scatter(l);
  };
  runBatchedAsync(m, nf, lanes, task, submit, complete);
}

void MobmDocker::dockPairs(std::vector<PairJob>& jobs, std::size_t chunkSize) const {
  const PairGroup g{this, &jobs};
  dockPairs(std::span<const PairGroup>(&g, 1), chunkSize);
}

void MobmDocker::prepareJobs(std::vector<PairJob>& jobs) const {
  for (PairJob& j : jobs) {
    j.prepared = prepare(j.mobAnglesDeg, j.stream, j.posm1Only);
    j.stats = j.prepared->stats;
  }
}

void MobmDocker::dockPairs(std::span<const PairGroup> groups, std::size_t chunkSize) {
  chunkSize = std::max<std::size_t>(chunkSize, 1);
  if (groups.empty()) return;
  PoseEvaluator* pev = groups[0].dk->p_.evaluator;
  if (!pev) {
    for (const PairGroup& g : groups)
      for (PairJob& j : *g.jobs) {
        Prepared pr = j.prepared ? std::move(*j.prepared) : g.dk->prepare(j.mobAnglesDeg, j.stream, j.posm1Only);
        j.prepared.reset();
        j.stats = g.dk->dockPrepared(std::move(pr), j.stream, chunkSize, j.sink);
      }
    return;
  }
  if (static_cast<int>(groups.size()) > std::max(1, pev->fields()))
    throw Error("dockPairs: more StaM conformations than the evaluator binds together");
  std::vector<const MobmDocker*> dks;
  for (const PairGroup& g : groups) {
    if (g.dk->p_.evaluator != pev || &g.dk->sys_ != &groups[0].dk->sys_)
      throw Error("dockPairs: groups of different MobMs or evaluators");
    dks.push_back(g.dk);
  }
  // The L2poss of all jobs in group and job order, minimized in global chunks; each job's sink
  // receives the same chunks as from dockPair, in group and job order.
  struct State {
    int group = 0;
    PairJob* job = nullptr;
    Prepared pr;
    std::vector<DockSample> res;
    std::vector<char> ok;
    std::size_t done = 0, delivered = 0;
  };
  std::vector<State> st;
  for (std::size_t g = 0; g < groups.size(); ++g)
    for (PairJob& j : *groups[g].jobs) {
      State& s = st.emplace_back();
      s.group = static_cast<int>(g);
      s.job = &j;
    }
  std::vector<std::pair<std::size_t, std::size_t>> order;  // (job, L2pos)
  for (std::size_t j = 0; j < st.size(); ++j) {
    PairJob& job = *st[j].job;
    st[j].pr = job.prepared ? std::move(*job.prepared) : dks[st[j].group]->prepare(job.mobAnglesDeg, job.stream, job.posm1Only);
    job.prepared.reset();
    job.stats = st[j].pr.stats;
    const std::size_t n = st[j].pr.l2.size();
    st[j].res.resize(n);
    st[j].ok.assign(n, 0);
    for (std::size_t k = 0; k < n; ++k) order.emplace_back(j, k);
  }
  const std::size_t global = std::max<std::size_t>(chunkSize, static_cast<std::size_t>(std::max(1, pev->maxBatch())));
  std::size_t next = 0;  // first job not fully delivered
  for (std::size_t begin = 0; begin < order.size(); begin += global) {
    const std::size_t m = std::min(global, order.size() - begin);
    std::vector<BatchItem> items(m);
    for (std::size_t i = 0; i < m; ++i) {
      const auto [j, k] = order[begin + i];
      items[i] = {dks[st[j].group], st[j].group, &st[j].pr.c, st[j].job->stream, st[j].pr.l2[k]};
    }
    std::vector<DockSample> res(m);
    std::vector<char> ok(m, 0);
    std::vector<Dock9bStats> ist(m);
    minimizeBatched(*pev, dks, items, res, ok, ist);
    for (std::size_t i = 0; i < m; ++i) {
      const auto [j, k] = order[begin + i];
      st[j].res[k] = std::move(res[i]);
      st[j].ok[k] = ok[i];
      ++st[j].done;
      PairJob& job = *st[j].job;
      job.stats.leftGrid += ist[i].leftGrid;
      job.stats.steps += ist[i].steps;
      job.stats.evaluations += ist[i].evaluations;
    }
    for (; next < st.size(); ++next) {
      State& s = st[next];
      const std::size_t n = s.pr.l2.size();
      while (s.delivered < n && (s.done == n || s.done - s.delivered >= chunkSize)) {
        const std::size_t e = std::min(n, s.delivered + chunkSize);
        std::vector<DockSample> chunk;
        for (std::size_t k = s.delivered; k < e; ++k)
          if (s.ok[k]) chunk.push_back(std::move(s.res[k]));
        s.job->sink(chunk);
        s.delivered = e;
      }
      if (s.delivered < n) break;
      s.res = {};
      s.ok = {};
    }
  }
}

std::uint64_t dockStream(std::uint64_t round, std::uint64_t mobm, std::uint64_t stamSample, std::uint64_t mobmSample) {
  std::uint64_t s = splitmix(0x9b9b9b9bULL ^ round);
  s = splitmix(s ^ mobm);
  s = splitmix(s ^ stamSample);
  return splitmix(s ^ mobmSample);
}

std::vector<StamSepa> stamSepas(const StamField& field, const Dock9bParams& p, StamSepaInfo* info, int neighbours) {
  const GridGeometry& g = field.geom();
  const std::size_t n = static_cast<std::size_t>(g.n[0]) * g.n[1] * g.n[2];
  std::vector<float> phi(n, 0.0f);
  std::vector<char> shell(n, 0);
  const int m0 = 0;
#pragma omp parallel for schedule(static)
  for (int z = 0; z < g.n[2]; ++z)
    for (int y = 0; y < g.n[1]; ++y)
      for (int x = 0; x < g.n[0]; ++x) {
        const Vec3 pt = g.point(x, y, z);
        const double d = field.distance<double>(pt, nullptr);
        if (d < p.sepaShellMin || d > p.sepaShellMax) continue;
        const std::size_t i = g.index(x, y, z);
        shell[i] = 1;
        // Potential = energy per unit charge: (E(q = 1) - E(q = -1)) / 2 (spec 8.1).
        phi[i] = static_cast<float>(0.5 * (field.energy<double>(m0, 1.0, pt, nullptr) - field.energy<double>(m0, -1.0, pt, nullptr)));
      }
  std::vector<float> mags;
  for (std::size_t i = 0; i < n; ++i)
    if (shell[i]) mags.push_back(std::abs(phi[i]));
  StamSepaInfo inf;
  inf.shellPoints = static_cast<long long>(mags.size());
  double cut;
  if (p.sepaPotThreshold) {
    cut = *p.sepaPotThreshold;
  } else if (mags.empty()) {
    cut = std::numeric_limits<double>::infinity();
  } else {
    // The points ranked in the top sepaPotPercentile % of |potential| (at least one).
    const std::size_t keep = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::ceil(p.sepaPotPercentile / 100.0 * static_cast<double>(mags.size()) - 1e-9)));
    std::nth_element(mags.begin(), mags.begin() + static_cast<std::ptrdiff_t>(keep - 1), mags.end(), std::greater<float>());
    cut = mags[keep - 1];
  }
  inf.potCut = cut;
  std::vector<char> sel(n, 0);
  for (std::size_t i = 0; i < n; ++i)
    if (shell[i] && phi[i] != 0.0f && std::abs(phi[i]) >= cut) {
      sel[i] = phi[i] > 0 ? 1 : -1;
      ++inf.selected;
    }
  // Connected regions of the same sign.
  std::vector<std::array<int, 3>> offs;
  for (int dz = -1; dz <= 1; ++dz)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx) {
        const int k = std::abs(dx) + std::abs(dy) + std::abs(dz);
        if (k == 0 || (neighbours == 6 && k > 1)) continue;
        offs.push_back({dx, dy, dz});
      }
  std::vector<StamSepa> out;
  std::vector<char> seen(n, 0);
  std::vector<std::array<int, 3>> stack;
  for (int z = 0; z < g.n[2]; ++z)
    for (int y = 0; y < g.n[1]; ++y)
      for (int x = 0; x < g.n[0]; ++x) {
        const std::size_t i0 = g.index(x, y, z);
        if (!sel[i0] || seen[i0]) continue;
        StamSepa s;
        s.sign = sel[i0];
        seen[i0] = 1;
        stack.assign(1, {x, y, z});
        while (!stack.empty()) {
          const auto c = stack.back();
          stack.pop_back();
          const Vec3 pt = g.point(c[0], c[1], c[2]);
          for (int k = 0; k < 3; ++k) s.xyz[k] += pt[k];
          ++s.points;
          for (const auto& o : offs) {
            const int a = c[0] + o[0], b = c[1] + o[1], e = c[2] + o[2];
            if (a < 0 || b < 0 || e < 0 || a >= g.n[0] || b >= g.n[1] || e >= g.n[2]) continue;
            const std::size_t j = g.index(a, b, e);
            if (sel[j] != s.sign || seen[j]) continue;
            seen[j] = 1;
            stack.push_back({a, b, e});
          }
        }
        for (double& v : s.xyz) v /= s.points;
        out.push_back(s);
      }
  if (info) *info = inf;
  return out;
}

namespace {

// Rotation taking the unit vector a to the unit vector b.
Quat alignQuat(const Vec3& a, const Vec3& b) {
  const double d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  if (d < -1 + 1e-12) {
    // Antiparallel: half turn about any axis perpendicular to a.
    Vec3 p = std::abs(a[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    Vec3 ax{a[1] * p[2] - a[2] * p[1], a[2] * p[0] - a[0] * p[2], a[0] * p[1] - a[1] * p[0]};
    const double l = std::sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
    return {0, ax[0] / l, ax[1] / l, ax[2] / l};
  }
  return normalized(Quat{1 + d, a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]});
}

Vec3 unit(const Vec3& v) {
  const double l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  return {v[0] / l, v[1] / l, v[2] / l};
}

// Height of the triangle over its longest side.
double triangleHeight(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  const Vec3 x{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
  const double area2 = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double longest = std::sqrt(std::max({dist2(a, b), dist2(b, c), dist2(a, c)}));
  return longest > 0 ? area2 / longest : 0.0;
}

// Proper rotation minimising sum |R (from_i - cf) - (to_i - ct)|^2 (Kabsch), as a quaternion; rmsd receives the
// superposition RMSD.
Quat kabschQuat(const std::array<Vec3, 3>& from, const std::array<Vec3, 3>& to, double& rmsd) {
  Eigen::Vector3d cf = Eigen::Vector3d::Zero(), ct = Eigen::Vector3d::Zero();
  for (int i = 0; i < 3; ++i) {
    cf += Eigen::Vector3d(from[i][0], from[i][1], from[i][2]) / 3.0;
    ct += Eigen::Vector3d(to[i][0], to[i][1], to[i][2]) / 3.0;
  }
  Eigen::Matrix3d h = Eigen::Matrix3d::Zero();
  for (int i = 0; i < 3; ++i)
    h += (Eigen::Vector3d(from[i][0], from[i][1], from[i][2]) - cf) *
         (Eigen::Vector3d(to[i][0], to[i][1], to[i][2]) - ct).transpose();
  const Eigen::JacobiSVD<Eigen::Matrix3d> svd(h, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d d = Eigen::Matrix3d::Identity();
  d(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant() < 0 ? -1.0 : 1.0;
  const Eigen::Matrix3d r = svd.matrixV() * d * svd.matrixU().transpose();
  double sq = 0;
  for (int i = 0; i < 3; ++i)
    sq += (r * (Eigen::Vector3d(from[i][0], from[i][1], from[i][2]) - cf) -
           (Eigen::Vector3d(to[i][0], to[i][1], to[i][2]) - ct)).squaredNorm();
  rmsd = std::sqrt(sq / 3.0);
  const Eigen::Quaterniond q(r);
  return normalized(Quat{q.w(), q.x(), q.y(), q.z()});
}

}  // namespace

MobmDocker::MobmDocker(const EnergySystem& sys, const StamField& field, const Dock9bParams& p,
                       const std::vector<StamSepa>* sepas, const FftStam* fft)
    : sys_(sys), field_(field), p_(p), sepas_(sepas), fft_(fft) {
  if (p.posm == Posm::Posm2Fft && !fft) fail("posm2_fft: the StaM spectra are missing");
  const Ad4Typing& t = sys.typing(0);
  const int off = sys.coordOffset(0);
  for (std::size_t a = 0; a < t.type.size(); ++a) {
    if (t.type[a] == Ad4Type::Merged) continue;
    const int m = field.mapOf(t.type[a]);
    if (m < 0) fail("no grid map for AD4 type {}", ad4TypeName(t.type[a]));
    probe_.push_back(off + static_cast<int>(a));
    map_.push_back(m);
    charge_.push_back(t.charge[a]);
    heavy_.push_back(t.type[a] != Ad4Type::H && t.type[a] != Ad4Type::HD && t.type[a] != Ad4Type::HS);
  }
  for (int c = 0; c < 3; ++c) posN_[c] = static_cast<int>(std::floor(p.boxSize[c] / p.posStep + 1e-9)) + 1;
  if (p.posm != Posm::Posm1) hopf_ = hopfGrid(p.hopfLevel);
  if (p.posm == Posm::Posm3 || p.posm == Posm::Posm4) {
    if (!sepas_) fail("posm3: the StaM regions are missing");
    axisAngles_ = std::max(1, static_cast<int>(std::lround(360.0 / p.rotStep)));
    transN_ = 2 * static_cast<int>(std::floor(p.sepaPosRange / p.sepaPosStep + 1e-9)) + 1;
  }
}

long long MobmDocker::PosRange::count() const {
  long long n = 1;
  for (int d = 0; d < 3; ++d) n *= std::max(0, hi[d] - lo[d] + 1);
  return n;
}

long long MobmDocker::l1Count(const Conf& c) const {
  switch (p_.posm) {
    case Posm::Posm1: return p_.posm1SampleNum;
    case Posm::Posm2: return c.hopfOffset.back();
    case Posm::Posm3: return c.sepaPoses();
    case Posm::Posm4: return c.sepaPoses() + p_.posm1SampleNum;
    case Posm::Posm2Fft: return static_cast<long long>(c.fftPoses.size());
  }
  return 0;
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
  if (p_.posm == Posm::Posm2Fft) fftPoses(c);
  if (p_.posm == Posm::Posm3 || p_.posm == Posm::Posm4) {
    // MobM side: probe atoms with |q| >= sepaChargeThreshold (spec 8.4.4).
    std::vector<int> ms;
    for (std::size_t k = 0; k < probe_.size(); ++k)
      if (std::abs(charge_[k]) >= p_.sepaChargeThreshold && charge_[k] != 0) ms.push_back(static_cast<int>(k));
    auto sgn = [&](int k) { return charge_[k] > 0 ? 1 : -1; };
    auto yAt = [&](int k) {
      const Vec3f& y = c.coords[probe_[k]];
      return Vec3{y[0], y[1], y[2]};
    };
    const auto& ss = *sepas_;
    std::vector<Conf::PairMatch> pm;
    if (p_.sepaPair || p_.sepaTriple)
      for (std::size_t i = 0; i < ss.size(); ++i)
        for (std::size_t j = i + 1; j < ss.size(); ++j) {
          const double dS = std::sqrt(dist2(ss[i].xyz, ss[j].xyz));
          for (int a : ms)
            for (int b : ms) {
              if (a == b || sgn(a) != -ss[i].sign || sgn(b) != -ss[j].sign) continue;
              if (std::abs(dS - std::sqrt(dist2(yAt(a), yAt(b)))) > p_.sepaDistTolerance) continue;
              pm.push_back({static_cast<int>(i), static_cast<int>(j), a, b});
            }
        }
    // Triples extend a matched pair (s1 < s2) by a region s3 > s2 and a third probe atom (spec 8.4.4).
    if (p_.sepaTriple)
      for (const auto& q : pm)
        for (std::size_t k = static_cast<std::size_t>(q.s2) + 1; k < ss.size(); ++k) {
          const double d13 = std::sqrt(dist2(ss[q.s1].xyz, ss[k].xyz)), d23 = std::sqrt(dist2(ss[q.s2].xyz, ss[k].xyz));
          for (int e : ms) {
            if (e == q.m1 || e == q.m2 || sgn(e) != -ss[k].sign) continue;
            if (std::abs(d13 - std::sqrt(dist2(yAt(q.m1), yAt(e)))) > p_.sepaDistTolerance ||
                std::abs(d23 - std::sqrt(dist2(yAt(q.m2), yAt(e)))) > p_.sepaDistTolerance)
              continue;
            const std::array<Vec3, 3> sv{ss[q.s1].xyz, ss[q.s2].xyz, ss[k].xyz}, mv{yAt(q.m1), yAt(q.m2), yAt(e)};
            if (triangleHeight(sv[0], sv[1], sv[2]) < kTripleMinHeight ||
                triangleHeight(mv[0], mv[1], mv[2]) < kTripleMinHeight)
              continue;
            double rmsd;
            const Quat rq = kabschQuat(mv, sv, rmsd);
            if (rmsd > p_.sepaDistTolerance) continue;
            c.triples.push_back({{q.s1, q.s2, static_cast<int>(k)}, {q.m1, q.m2, e}, rq});
          }
        }
    if (p_.sepaPair) c.pairs = std::move(pm);
    if (p_.sepaSingle)
      for (std::size_t i = 0; i < ss.size(); ++i)
        for (int a : ms)
          if (sgn(a) == -ss[i].sign) c.singles.push_back({static_cast<int>(i), a});
    const long long nt = static_cast<long long>(transN_) * transN_ * transN_;
    c.pairPoses = static_cast<long long>(c.pairs.size()) * axisAngles_ * nt;
    c.singlePoses = static_cast<long long>(c.singles.size()) * static_cast<long long>(hopf_.size()) * nt;
    c.triplePoses = static_cast<long long>(c.triples.size()) * nt;
  }
  return c;
}

// posm2_fft (spec 8.4.3): per Hopf orientation, the fftTopPerRot lowest-energy grid translations
// keeping every probe atom kPosMargin inside the box.
void MobmDocker::fftPoses(Conf& c) const {
  const FftStam& s = *fft_;
  const GridGeometry& g = s.geom;
  const std::size_t no = hopf_.size(), np = probe_.size(), top = static_cast<std::size_t>(std::max(1, p_.fftTopPerRot));
  std::vector<std::vector<std::pair<float, std::size_t>>> kept(no);
  const Vec3f cf = castVec<float>(c.center);
  // Grid coordinates f (np × 3) of the probe atoms in orientation o and its translation range; false when empty.
  auto geometry = [&](std::size_t o, Vec3* f, std::array<int, 3>& lo, std::array<int, 3>& hi) {
    const QuatF qf = toFloat(hopf_[o]);
    bool empty = false;
    for (int d = 0; d < 3; ++d) {
      float rlo = 1e30f, rhi = -1e30f;
      for (std::size_t a = 0; a < np; ++a) {
        const Vec3f& y = c.coords[probe_[a]];
        const float r = rotate(qf, {y[0] - cf[0], y[1] - cf[1], y[2] - cf[2]})[d];
        f[a][d] = r / g.spacing;
        rlo = std::min(rlo, r);
        rhi = std::max(rhi, r);
      }
      // Centre t = origin + T spacing with boxMin + margin <= t + rlo and t + rhi <= boxMax - margin.
      const double bmin = p_.boxCenter[d] - 0.5 * p_.boxSize[d] + kPosMargin - g.origin[d];
      const double bmax = p_.boxCenter[d] + 0.5 * p_.boxSize[d] - kPosMargin - g.origin[d];
      lo[d] = std::max(0, static_cast<int>(std::ceil((bmin - rlo) / g.spacing - 1e-9)));
      hi[d] = std::min(g.n[d] - 1, static_cast<int>(std::floor((bmax - rhi) / g.spacing + 1e-9)));
      if (lo[d] > hi[d]) empty = true;
    }
    return !empty;
  };
  const int batch = p_.evaluator ? p_.evaluator->fftBatch(*this) : 0;
  if (batch > 0) {
    // Device correlation (spec 15.1.5): the non-empty orientations in batches.
    std::vector<std::array<int, 3>> lo(no), hi(no);
    std::vector<char> use(no);
    {
      std::vector<Vec3> f(np);
#pragma omp parallel for schedule(static) firstprivate(f)
      for (std::size_t o = 0; o < no; ++o) use[o] = geometry(o, f.data(), lo[o], hi[o]);
    }
    std::vector<std::size_t> which;
    for (std::size_t o = 0; o < no; ++o)
      if (use[o]) which.push_back(o);
    const std::size_t nb = static_cast<std::size_t>(batch);
    std::vector<Vec3> f(nb * np);
    std::vector<std::array<int, 3>> blo(nb), bhi(nb);
    std::vector<std::pair<float, std::size_t>> out(nb * top);
    for (std::size_t b0 = 0; b0 < which.size(); b0 += nb) {
      const std::size_t m = std::min(nb, which.size() - b0);
#pragma omp parallel for schedule(static)
      for (std::size_t b = 0; b < m; ++b) geometry(which[b0 + b], f.data() + b * np, blo[b], bhi[b]);
      p_.evaluator->fftTop(static_cast<int>(m), static_cast<int>(np), f[0].data(), map_.data(), charge_.data(), blo.data(),
                           bhi.data(), static_cast<int>(top), out.data());
      for (std::size_t b = 0; b < m; ++b)
        for (std::size_t r = 0; r < top; ++r)
          if (out[b * top + r].second != SIZE_MAX) kept[which[b0 + b]].push_back(out[b * top + r]);
    }
  } else {
#pragma omp parallel
    {
      FftWork w;
      std::vector<Vec3> f(np);
#pragma omp for schedule(dynamic, 1)
      for (std::size_t o = 0; o < no; ++o) {
        std::array<int, 3> lo{}, hi{};
        if (!geometry(o, f.data(), lo, hi)) continue;
        fftCorrelate(s, f, map_, charge_, w);
        fftSelect(s, w.energy.data(), lo, hi, top, kept[o]);
      }
    }
  }
  for (std::size_t o = 0; o < no; ++o)
    for (const auto& [e, idx] : kept[o]) {
      const int i = static_cast<int>(idx % s.pad[0]), j = static_cast<int>(idx / s.pad[0] % s.pad[1]),
                k = static_cast<int>(idx / s.pad[0] / s.pad[1]);
      c.fftPoses.push_back({g.point(i, j, k), hopf_[o]});
      c.fftEnergy.push_back(e);
    }
}

void MobmDocker::fitBounds(float lo[3], float hi[3]) const {
  const float m = static_cast<float>(kPosMargin);
  for (int d = 0; d < 3; ++d) {
    lo[d] = static_cast<float>(p_.boxCenter[d] - 0.5 * p_.boxSize[d]) + m;
    hi[d] = static_cast<float>(p_.boxCenter[d] + 0.5 * p_.boxSize[d]) - m;
  }
}

bool MobmDocker::fitsBox(const Conf& c, const Pose& pose) const {
  const QuatF qf = toFloat(pose.q);
  const Vec3f cf = castVec<float>(c.center), tf = castVec<float>(pose.t);
  float lo[3], hi[3];
  fitBounds(lo, hi);
  for (int k : probe_) {
    const Vec3f& y = c.coords[k];
    const Vec3f r = rotate(qf, {y[0] - cf[0], y[1] - cf[1], y[2] - cf[2]});
    for (int d = 0; d < 3; ++d) {
      const float w = r[d] + tf[d];
      if (w < lo[d] || w > hi[d]) return false;
    }
  }
  return true;
}

MobmDocker::RotBox MobmDocker::rotBox(const Conf& c, const Quat& q) const {
  const QuatF qf = toFloat(q);
  const Vec3f cf = castVec<float>(c.center);
  RotBox b;
  for (int d = 0; d < 3; ++d) {
    b.lo[d] = std::numeric_limits<float>::infinity();
    b.hi[d] = -std::numeric_limits<float>::infinity();
  }
  for (int k : probe_) {
    const Vec3f& y = c.coords[k];
    const Vec3f r = rotate(qf, {y[0] - cf[0], y[1] - cf[1], y[2] - cf[2]});
    for (int d = 0; d < 3; ++d) {
      b.lo[d] = std::min(b.lo[d], r[d]);
      b.hi[d] = std::max(b.hi[d], r[d]);
    }
  }
  return b;
}

bool MobmDocker::fitsBox(const RotBox& b, const Pose& pose) const {
  const Vec3f tf = castVec<float>(pose.t);
  float lo[3], hi[3];
  fitBounds(lo, hi);
  for (int d = 0; d < 3; ++d)
    if (b.lo[d] + tf[d] < lo[d] || b.hi[d] + tf[d] > hi[d]) return false;
  return true;
}

Vec3 MobmDocker::posm3Delta(int tr) const {
  const int half = transN_ / 2;
  Vec3 delta;
  for (int d = 0; d < 3; ++d) {
    delta[d] = p_.sepaPosStep * static_cast<double>(tr % transN_ - half);
    tr /= transN_;
  }
  return delta;
}

MobmDocker::Posm3Rot MobmDocker::posm3Rot(const Conf& c, long long r) const {
  const long long nt = posm3Translations();
  const auto& ss = *sepas_;
  auto yAt = [&](int pk) {
    const Vec3f& y = c.coords[probe_[pk]];
    return Vec3{y[0], y[1], y[2]};
  };
  Posm3Rot rot;
  Vec3 anchorM;  // conformation point made to coincide with rot.anchorS (before delta)
  Vec3& anchorS = rot.anchorS;
  if (r < c.pairPoses / nt) {
    // Pair: match-major, then the rotation about the axis, then the translation.
    const auto& pm = c.pairs[static_cast<std::size_t>(r / axisAngles_)];
    const int a = static_cast<int>(r % axisAngles_);
    const Vec3 s1 = ss[pm.s1].xyz, s2 = ss[pm.s2].xyz, m1 = yAt(pm.m1), m2 = yAt(pm.m2);
    const Vec3 us = unit({s2[0] - s1[0], s2[1] - s1[1], s2[2] - s1[2]});
    const Quat q0 = alignQuat(unit({m2[0] - m1[0], m2[1] - m1[1], m2[2] - m1[2]}), us);
    const double th = 2 * std::numbers::pi * a / axisAngles_;
    rot.q = normalized(fromRotationVector({us[0] * th, us[1] * th, us[2] * th}) * q0);
    for (int d = 0; d < 3; ++d) {
      anchorS[d] = 0.5 * (s1[d] + s2[d]);
      anchorM[d] = 0.5 * (m1[d] + m2[d]);
    }
  } else if (r >= (c.pairPoses + c.singlePoses) / nt) {
    // Triple: match-major, then the translation; the triangle centroids coincide.
    const auto& tm = c.triples[static_cast<std::size_t>(r - (c.pairPoses + c.singlePoses) / nt)];
    rot.q = tm.q;
    anchorS = anchorM = Vec3{0, 0, 0};
    for (int v = 0; v < 3; ++v) {
      const Vec3 sv = ss[tm.s[v]].xyz, mv = yAt(tm.m[v]);
      for (int d = 0; d < 3; ++d) {
        anchorS[d] += sv[d] / 3.0;
        anchorM[d] += mv[d] / 3.0;
      }
    }
  } else {
    // Single: placement-major, then the Hopf orientation, then the translation.
    const long long rs = r - c.pairPoses / nt;
    const long long nh = static_cast<long long>(hopf_.size());
    const auto& sp = c.singles[static_cast<std::size_t>(rs / nh)];
    rot.q = hopf_[static_cast<std::size_t>(rs % nh)];
    anchorS = ss[sp[0]].xyz;
    anchorM = yAt(sp[1]);
  }
  rot.rm = rotate(rot.q, {anchorM[0] - c.center[0], anchorM[1] - c.center[1], anchorM[2] - c.center[2]});
  return rot;
}

Pose MobmDocker::posm3Pose(const Posm3Rot& rot, int tr) const {
  const Vec3 delta = posm3Delta(tr);
  Pose pose;
  pose.q = rot.q;
  for (int d = 0; d < 3; ++d) pose.t[d] = rot.anchorS[d] + delta[d] - rot.rm[d];
  return pose;
}

std::optional<Pose> MobmDocker::posm3Pose(const Conf& c, long long k) const {
  const long long nt = posm3Translations();
  const Pose pose = posm3Pose(posm3Rot(c, k / nt), static_cast<int>(k % nt));
  if (!fitsBox(c, pose)) return std::nullopt;
  return pose;
}

std::optional<Pose> MobmDocker::l1Pose(const Conf& conf, std::uint64_t stream, long long k) const {
  if (p_.posm == Posm::Posm2Fft) return conf.fftPoses[static_cast<std::size_t>(k)];
  if (p_.posm == Posm::Posm3 || p_.posm == Posm::Posm4) {
    const long long n3 = conf.sepaPoses();
    if (k < n3) return posm3Pose(conf, k);
    k -= n3;  // posm4: the posm1 samples follow
  }
  Pose pose;
  std::array<long long, 3> idx{};
  if (p_.posm != Posm::Posm2) {
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
  return quickFilterT<false>(c, pose, stats);
}

bool MobmDocker::quickFilterInBox(const Conf& c, const Pose& pose, Dock9bStats& stats) const {
  return quickFilterT<true>(c, pose, stats);
}

template <bool InBox>
bool MobmDocker::quickFilterT(const Conf& c, const Pose& pose, Dock9bStats& stats) const {
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
    if constexpr (!InBox)
      for (int d = 0; d < 3; ++d)
        if (w[d] < boxLo[d] || w[d] > boxHi[d]) {
          ++stats.outsideBox;
          return false;
        }
    if (!heavy_[k]) continue;
    const float dd = field_.distance<float>(w, nullptr);
    if (dd < cd && ++overlap >= p_.maxOverlapAtomNum && InBox) {
      ++stats.overlap;
      return false;
    }
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
  const V gLo = castVec<T>(field_.geom().origin), gHi = castVec<T>(field_.geom().max());
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
      const T q = static_cast<T>(charge_[k]);
      V gw{}, dg{};
      pe.inter += field_.energy<T>(map_[k], q, w, g ? &gw : nullptr);
      if (heavy_[k]) {
        const T d = field_.distance<T>(w, g ? &dg : nullptr);
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

std::vector<DockProbe> MobmDocker::probes() const {
  std::vector<DockProbe> v(probe_.size());
  for (std::size_t k = 0; k < probe_.size(); ++k) v[k] = {probe_[k], map_[k], charge_[k], heavy_[k] ? 1 : 0};
  return v;
}

DockView MobmDocker::dockView(const std::vector<DockProbe>& probes, const FieldView& field) const {
  return {static_cast<int>(probes.size()), probes.data(), p_.clashDistance, p_.clashBase, p_.clashK, field, nullptr};
}

bool MobmDocker::minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats) const {
  return minimize(c, start, out, stats, [&](const Pose& pose, std::span<const double> a, std::vector<double>* g) {
    return evaluate(c, pose, a, g);
  });
}

bool MobmDocker::minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats,
                          const PoseEval& ev) const {
  const std::size_t nt = c.anglesDeg.size();
  // adadelta: AutoDock-GPU units (translation in grid spacings, angles in degrees).
  const bool ada = p_.minimizer.method == MinMethod::Adadelta;
  const double uT = ada ? field_.geom().spacing : 1.0, uA = ada ? kDeg : 1.0;
  const double sT = p_.scaleTrans * uT, sR = p_.scaleRot * uA, sA = p_.scaleTors * uA;
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
  Point trial, best;
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
    trial.pe = ev(trial.pose, trial.angles, &g);
    if (trial.pe.outside > 0) return false;
    if (!started) {
      best = trial;
      started = true;
    }
    e = trial.pe.objective();
    grad.resize(g.size());
    for (std::size_t i = 0; i < g.size(); ++i) grad[i] = g[i] * (i < 3 ? sT : i < 6 ? sR : sA);
    return true;
  };
  VecMinParams lp = p_.minimizer;
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
  // An accepted point is the last one evaluated (a line search stops at it).
  lp.onAccept = [&](std::vector<double>& v, bool isBest) {
    if (isBest) best = trial;
    q0 = trial.pose.q;
    v[3] = v[4] = v[5] = 0;
  };
  const VecMinResult r = minimizeVec(f, z, lp);
  stats.steps += r.steps;
  stats.evaluations += r.evaluations;
  if (r.aborted) {
    ++stats.leftGrid;
    return false;
  }
  const PoseEnergy& pe = best.pe;
  out.pose = best.pose;
  out.angles.assign(best.angles.begin(), best.angles.end());
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
  return dockPair(mobAnglesDeg, stream, std::numeric_limits<std::size_t>::max(), [&](std::vector<DockSample>& chunk) {
    out.insert(out.end(), std::make_move_iterator(chunk.begin()), std::make_move_iterator(chunk.end()));
  });
}

MobmDocker::Prepared MobmDocker::prepare(std::span<const double> mobAnglesDeg, std::uint64_t stream,
                                         bool posm1Only) const {
  Prepared pr;
  pr.c = conformation(mobAnglesDeg);
  const Conf& c = pr.c;
  const long long n = l1Count(c);
  const long long first = posm1Only && p_.posm == Posm::Posm4 ? c.sepaPoses() : 0;
  Dock9bStats& total = pr.stats;
  std::vector<long long>& l2 = pr.l2;
  total.l1 = n - first;
  if (first == 0) {
    total.sepaPairs = static_cast<long long>(c.pairs.size());
    total.sepaSingles = static_cast<long long>(c.singles.size());
    total.sepaTriples = static_cast<long long>(c.triples.size());
  }
  if (!p_.evaluator || !p_.evaluator->quickFilter(*this, c, stream, first, n, l2, total)) {
    // posm3 part: per rotation, the box check of its translations from the rotated extent (spec 15.1.21).
    const long long n3 = p_.posm == Posm::Posm3 || p_.posm == Posm::Posm4 ? std::min(n, c.sepaPoses()) : first;
    const long long nt = std::max(posm3Translations(), 1);
#pragma omp parallel
    {
      Dock9bStats st;
      std::vector<long long> mine;
#pragma omp for schedule(static) nowait
      for (long long r = first / nt; r < (n3 + nt - 1) / nt; ++r) {
        const Posm3Rot rot = posm3Rot(c, r);
        const RotBox box = rotBox(c, rot.q);
        for (long long k = std::max(first, r * nt); k < std::min(n3, (r + 1) * nt); ++k) {
          const Pose pose = posm3Pose(rot, static_cast<int>(k - r * nt));
          if (!fitsBox(box, pose)) ++st.outsideBox;
          else if (quickFilterInBox(c, pose, st)) mine.push_back(k);
        }
      }
#pragma omp for schedule(static) nowait
      for (long long k = std::max(first, n3); k < n; ++k) {
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
  }
  if (p_.evaluator) p_.evaluator->prepared();
  std::sort(l2.begin(), l2.end());
  total.l2 = static_cast<long long>(l2.size());
  return pr;
}

Dock9bStats MobmDocker::dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream, std::size_t chunkSize,
                                 const std::function<void(std::vector<DockSample>&)>& sink, bool posm1Only) const {
  return dockPrepared(prepare(mobAnglesDeg, stream, posm1Only), stream, chunkSize, sink);
}

Dock9bStats MobmDocker::dockPrepared(Prepared pr, std::uint64_t stream, std::size_t chunkSize,
                                     const std::function<void(std::vector<DockSample>&)>& sink) const {
  const Conf& c = pr.c;
  const std::vector<long long>& l2 = pr.l2;
  Dock9bStats& total = pr.stats;
  chunkSize = std::max<std::size_t>(chunkSize, 1);
  for (std::size_t begin = 0; begin < l2.size(); begin += chunkSize) {
    const std::size_t m = std::min(chunkSize, l2.size() - begin);
    std::vector<DockSample> res(m);
    std::vector<char> ok(m, 0);
    if (p_.evaluator) {
      std::vector<BatchItem> items(m);
      for (std::size_t i = 0; i < m; ++i) items[i] = {this, 0, &c, stream, l2[begin + i]};
      std::vector<Dock9bStats> st(m);
      const MobmDocker* self = this;
      minimizeBatched(*p_.evaluator, std::span<const MobmDocker* const>(&self, 1), items, res, ok, st);
      for (const auto& s : st) {
        total.leftGrid += s.leftGrid;
        total.steps += s.steps;
        total.evaluations += s.evaluations;
      }
    } else {
#pragma omp parallel
      {
        Dock9bStats st;
#pragma omp for schedule(dynamic, 1) nowait
        for (std::size_t i = 0; i < m; ++i) ok[i] = minimize(c, *l1Pose(c, stream, l2[begin + i]), res[i], st);
#pragma omp critical
        {
          total.leftGrid += st.leftGrid;
          total.steps += st.steps;
          total.evaluations += st.evaluations;
        }
      }
    }
    std::vector<DockSample> chunk;
    for (std::size_t i = 0; i < m; ++i)
      if (ok[i]) chunk.push_back(std::move(res[i]));
    sink(chunk);
  }
  return total;
}

}  // namespace od
