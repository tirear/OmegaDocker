// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "csgm1.hpp"

#include <algorithm>
#include <cmath>

#include "log.hpp"

namespace od {

RamaSampler::RamaSampler(const RamaGrid& g, double floor) {
  cdf_.resize(g.value.size());
  double s = 0;
  for (std::size_t c = 0; c < g.value.size(); ++c) {
    s += g.value[c] + floor;
    cdf_[c] = s;
  }
}

void RamaSampler::sample(Rng& rng, double& phi, double& psi) const {
  const double u = rng.uniform() * cdf_.back();
  const auto c = static_cast<int>(std::min<std::size_t>(std::upper_bound(cdf_.begin(), cdf_.end(), u) - cdf_.begin(),
                                                        cdf_.size() - 1));
  phi = wrapDeg(-180.0 + 2.0 * (c / kRamaBins) + rng.uniform(0.0, 2.0));
  psi = wrapDeg(-180.0 + 2.0 * (c % kRamaBins) + rng.uniform(0.0, 2.0));
}

double RamaSampler::probability(double phi, double psi) const {
  const int c = RamaGrid::bin(phi) * kRamaBins + RamaGrid::bin(psi);
  return (cdf_[c] - (c > 0 ? cdf_[c - 1] : 0.0)) / cdf_.back();
}

Csgm1Sampler::Csgm1Sampler(const Molecule& m, const MolSettings& s, const MolConfig& config, const Csgm1Params& params)
    : ssTypes_(config.ssTypes),
      halfWidth_(config.ssHalfWidth),
      originalProb_(config.ssOriginalProb),
      params_(params),
      rama_(std::make_shared<std::vector<std::unique_ptr<RamaSampler>>>(kRamaClassCount)) {
  if (!(params.angleStep > 0 && params.angleStep <= 180)) fail("csgm1AngleStep must be within (0, 180]");
  stepMin_ = static_cast<int>(std::floor(-180.0 / params.angleStep)) + 1;
  stepCount_ = static_cast<int>(std::floor(180.0 / params.angleStep)) - stepMin_ + 1;

  const Topology& t = s.topology;
  const auto& tree = s.tree.torsions;
  std::vector<int> treeOf(t.torsions.size(), -1);
  for (std::size_t i = 0; i < tree.size(); ++i) {
    treeOf[tree[i].torsion] = static_cast<int>(i);
    input_.push_back(tree[i].inputAngle);
    kind_.push_back(s.mode[tree[i].torsion] == TorsionMode::CisTrans ? Kind::CisTrans : Kind::Uniform);
  }
  if (t.backboneType == BackboneType::Protein && s.csgm1 == Csgm1::B) initPairs(m, s, treeOf);
  std::vector<int> pairOf(input_.size(), -1);
  for (std::size_t k = 0; k < pairs_.size(); ++k) {
    const Pair& p = pairs_[k];
    std::vector<int> seg{p.phi, p.psi};
    if (p.rule == Rule::RamaAuto && p.omegaBefore >= 0) seg.push_back(p.omegaBefore);
    std::sort(seg.begin(), seg.end());
    for (int i : seg) pairOf[i] = static_cast<int>(k);
    segments_.push_back(std::move(seg));
    segmentPair_.push_back(static_cast<int>(k));
  }
  for (std::size_t i = 0; i < input_.size(); ++i)
    if (pairOf[i] < 0) {
      segments_.push_back({static_cast<int>(i)});
      segmentPair_.push_back(-1);
    }
  std::vector<std::size_t> order(segments_.size());
  for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
  std::sort(order.begin(), order.end(), [&](auto x, auto y) { return segments_[x].front() < segments_[y].front(); });
  std::vector<std::vector<int>> seg;
  std::vector<int> pair;
  for (auto k : order) {
    seg.push_back(std::move(segments_[k]));
    pair.push_back(segmentPair_[k]);
  }
  segments_ = std::move(seg);
  segmentPair_ = std::move(pair);
}

void Csgm1Sampler::initPairs(const Molecule& m, const MolSettings& s, const std::vector<int>& treeOf) {
  const Topology& t = s.topology;
  std::vector<int> phiOf(t.units.size(), -1), psiOf(t.units.size(), -1), omegaOf(t.units.size(), -1);
  for (std::size_t i = 0; i < t.torsions.size(); ++i) {
    const auto& tor = t.torsions[i];
    if (tor.backbone == BackboneAngle::Phi) phiOf[tor.unit] = static_cast<int>(i);
    if (tor.backbone == BackboneAngle::Psi) psiOf[tor.unit] = static_cast<int>(i);
    if (tor.backbone == BackboneAngle::Omega) omegaOf[tor.unit] = static_cast<int>(i);
  }
  auto resName = [&](int u) { return u >= 0 && t.units[u].residue >= 0 ? m.residues[t.units[u].residue].info.name : ""; };
  auto need = [&](RamaClass c) {
    auto& r = (*rama_)[static_cast<int>(c)];
    if (!r) r = std::make_unique<RamaSampler>(ramaGrid(c), params_.ramaFloor);
  };
  for (std::size_t u = 0; u < t.units.size(); ++u) {
    if (phiOf[u] < 0 || psiOf[u] < 0 || treeOf[phiOf[u]] < 0 || treeOf[psiOf[u]] < 0) continue;
    const std::string& name = s.ssType[u];
    if (name == "R") continue;
    Pair p;
    p.phi = treeOf[phiOf[u]];
    p.psi = treeOf[psiOf[u]];
    if (name == "O") {
      p.rule = Rule::O;
    } else if (name == "RAMA_AUTO") {
      p.rule = Rule::RamaAuto;
      p.resName = resName(static_cast<int>(u));
      p.nextResName = resName(t.units[u].next);
      if (const int prev = t.units[u].prev; prev >= 0 && omegaOf[prev] >= 0) {
        const int o = omegaOf[prev];
        if (s.mode[o] == TorsionMode::CisTrans) p.omegaBefore = treeOf[o];
        else p.cisBefore = std::abs(t.torsions[o].angle) < 90.0;
      }
      need(ramaClassOf(p.resName, p.nextResName, p.cisBefore));
      if (p.omegaBefore >= 0) need(ramaClassOf(p.resName, p.nextResName, !p.cisBefore));
    } else if (auto c = ramaClassFromTypeName(name)) {
      p.rule = Rule::Rama;
      p.rama = *c;
      need(*c);
    } else {
      p.rule = Rule::Custom;
      for (std::size_t k = 0; k < ssTypes_.size(); ++k)
        if (ssTypes_[k].name == name) p.custom = static_cast<int>(k);
      if (p.custom < 0) fail("ssType {} is not defined", name);
    }
    kind_[p.phi] = Kind::PairPhi;
    kind_[p.psi] = Kind::PairPsi;
    pairs_.push_back(std::move(p));
  }
}

double Csgm1Sampler::uniformAngle(Rng& rng, int i) const {
  const int k = stepMin_ + static_cast<int>(rng.below(static_cast<std::uint64_t>(stepCount_)));
  return wrapDeg(input_[i] + k * params_.angleStep);
}

const RamaSampler& Csgm1Sampler::rama(RamaClass c) const { return *(*rama_)[static_cast<int>(c)]; }

void Csgm1Sampler::samplePair(Rng& rng, const Pair& p, std::span<double> a) const {
  const double hw = halfWidth_;
  auto box = [&](double phi0, double psi0) {
    a[p.phi] = wrapDeg(phi0 + rng.uniform(-hw, hw));
    a[p.psi] = wrapDeg(psi0 + rng.uniform(-hw, hw));
  };
  auto uniform = [&] {
    a[p.phi] = uniformAngle(rng, p.phi);
    a[p.psi] = uniformAngle(rng, p.psi);
  };
  switch (p.rule) {
    case Rule::R: uniform(); return;
    case Rule::O:
      if (rng.uniform() < originalProb_) box(input_[p.phi], input_[p.psi]);
      else uniform();
      return;
    case Rule::Rama: rama(p.rama).sample(rng, a[p.phi], a[p.psi]); return;
    case Rule::RamaAuto: {
      const bool cis = p.omegaBefore >= 0 ? std::abs(a[p.omegaBefore]) < 90.0 : p.cisBefore;
      rama(ramaClassOf(p.resName, p.nextResName, cis)).sample(rng, a[p.phi], a[p.psi]);
      return;
    }
    case Rule::Custom: {
      const auto& pts = ssTypes_[p.custom].points;
      double u = rng.uniform(), total = 0;
      for (const auto& pt : pts) total += pt.prob;
      for (const auto& pt : pts) {
        if (u < pt.prob) return box(pt.phi, pt.psi);
        u -= pt.prob;
      }
      // Outside all point ranges, uniform; if those cover (almost) everything, a point by its relative probability.
      for (int tries = 0; tries < 100000; ++tries) {
        const double phi = rng.uniform(-180.0, 180.0), psi = rng.uniform(-180.0, 180.0);
        bool inside = false;
        for (const auto& pt : pts)
          inside = inside || (std::abs(wrapDeg(phi - pt.phi)) <= hw && std::abs(wrapDeg(psi - pt.psi)) <= hw);
        if (!inside) {
          a[p.phi] = wrapDeg(phi);
          a[p.psi] = wrapDeg(psi);
          return;
        }
      }
      u = rng.uniform() * total;
      for (const auto& pt : pts) {
        if (u < pt.prob) return box(pt.phi, pt.psi);
        u -= pt.prob;
      }
      return box(pts.back().phi, pts.back().psi);
    }
  }
}

void Csgm1Sampler::sampleSegment(Rng& rng, int segment, std::span<double> a) const {
  for (int i : segments_[segment]) {
    if (kind_[i] == Kind::Uniform) a[i] = uniformAngle(rng, i);
    else if (kind_[i] == Kind::CisTrans) a[i] = rng.uniform() < 0.5 ? 0.0 : 180.0;
  }
  if (const int k = segmentPair_[segment]; k >= 0) samplePair(rng, pairs_[k], a);
}

void Csgm1Sampler::sample(Rng& rng, std::span<double> a) const {
  for (std::size_t i = 0; i < kind_.size(); ++i) {
    if (kind_[i] == Kind::Uniform) a[i] = uniformAngle(rng, static_cast<int>(i));
    else if (kind_[i] == Kind::CisTrans) a[i] = rng.uniform() < 0.5 ? 0.0 : 180.0;
  }
  for (const Pair& p : pairs_) samplePair(rng, p, a);
}

int Csgm1GroupSampler::size() const {
  int n = 0;
  for (const auto& p : parts_) n += p.size();
  return n;
}

std::vector<std::vector<int>> Csgm1GroupSampler::segments() const {
  std::vector<std::vector<int>> all;
  int o = 0;
  for (const auto& p : parts_) {
    for (auto seg : p.segments()) {
      for (int& i : seg) i += o;
      all.push_back(std::move(seg));
    }
    o += p.size();
  }
  return all;
}

void Csgm1GroupSampler::sampleSegment(Rng& rng, int segment, std::span<double> a) const {
  int o = 0;
  for (const auto& p : parts_) {
    const int n = static_cast<int>(p.segments().size());
    if (segment < n) return p.sampleSegment(rng, segment, a.subspan(o, p.size()));
    segment -= n;
    o += p.size();
  }
}

void Csgm1GroupSampler::sample(Rng& rng, std::span<double> a) const {
  int o = 0;
  for (const auto& p : parts_) {
    p.sample(rng, a.subspan(o, p.size()));
    o += p.size();
  }
}

}  // namespace od
