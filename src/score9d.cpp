// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "score9d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "ad4_pair.hpp"

namespace od {

std::vector<Vec3f> dockWorld(std::span<const Vec3f> local, const DockSample& s) {
  const QuatF q = toFloat(s.pose.q);
  const Vec3f c{static_cast<float>(s.center[0]), static_cast<float>(s.center[1]), static_cast<float>(s.center[2])};
  const Vec3f t{static_cast<float>(s.pose.t[0]), static_cast<float>(s.pose.t[1]), static_cast<float>(s.pose.t[2])};
  std::vector<Vec3f> w(local.size());
  for (std::size_t i = 0; i < local.size(); ++i) {
    const Vec3f r = rotate(q, {local[i][0] - c[0], local[i][1] - c[1], local[i][2] - c[2]});
    w[i] = {r[0] + t[0], r[1] + t[1], r[2] + t[2]};
  }
  return w;
}

InterTerms exactInter(const ShakeAtoms& stam, const ShakeAtoms& mob, const EnergyParams& ep) {
  const Ad4Coeffs& c = ep.coeffs;
  const float half = static_cast<float>(0.5 * ep.smooth);
  const float rMin = static_cast<float>(ep.elecMinDistance);
  const float cd = static_cast<float>(ep.clashDistance);
  InterTerms out;
  double vdw = 0, hbond = 0, elec = 0, desolv = 0, clash = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : vdw, hbond, elec, desolv, clash)
  for (std::size_t i = 0; i < mob.size(); ++i) {
    const Ad4TypeParams& pi = ad4Params(mob.type[i]);
    const float Si = static_cast<float>(pi.solpar + c.qasp * std::abs(mob.charge[i]));
    float nearest = std::numeric_limits<float>::infinity();
    double sv = 0, sh = 0, se = 0, sd = 0;
    for (std::size_t j = 0; j < stam.size(); ++j) {
      const float dx = mob.xyz[i][0] - stam.xyz[j][0], dy = mob.xyz[i][1] - stam.xyz[j][1],
                  dz = mob.xyz[i][2] - stam.xyz[j][2];
      const float r2 = dx * dx + dy * dy + dz * dz;
      if (r2 >= static_cast<float>(kElecCutoff * kElecCutoff)) continue;
      const float r = std::sqrt(r2);
      if (stam.heavy[j] && ad4Clashes(stam.type[j])) nearest = std::min(nearest, r);
      const Ad4TypeParams& pj = ad4Params(stam.type[j]);
      const float Sj = static_cast<float>(pj.solpar + c.qasp * std::abs(stam.charge[j]));
      float ee, de;
      if (r < static_cast<float>(kVdwCutoff)) {
        const Ad4PairParams q = ad4PairParams(mob.type[i], stam.type[j], c);
        ad4VdwHb<float>(r, static_cast<float>(q.AC), static_cast<float>(q.BD), q.m, q.n, static_cast<float>(q.reqm),
                        half, ee, de);
        (q.n == 10 ? sh : sv) += ee;
      }
      ad4Elec<float>(r, static_cast<float>(c.elec) * mob.charge[i] * stam.charge[j], rMin, ee, de);
      se += ee;
      ad4Desolv<float>(r, static_cast<float>(c.desolv) * (Si * static_cast<float>(pj.volume) + Sj * static_cast<float>(pi.volume)),
                       rMin, ee, de);
      sd += ee;
    }
    vdw += sv;
    hbond += sh;
    elec += se;
    desolv += sd;
    if (mob.heavy[i] && ad4Clashes(mob.type[i]) && nearest < cd) clash += ep.clashBase + ep.clashK * (cd - nearest) * (cd - nearest);
  }
  out.vdw = vdw;
  out.hbond = hbond;
  out.elec = elec;
  out.desolv = desolv;
  out.clash = clash;
  return out;
}

ComplexLnF complexSumLnF(const ShakeAtoms& stam, std::span<const double> stamAlone, const ShakeAtoms& mob,
                         const EnergyParams& ep, const ShakeParams& p, const ShakeFn& fn) {
  const ShakeAtoms all = joinShakeAtoms(stam, mob);
  const double reach = kElecCutoff + p.radius;
  std::vector<char> shake(all.size(), 1);
  ComplexLnF out;
  for (std::size_t j = 0; j < stam.size(); ++j) {
    if (!stam.heavy[j]) continue;
    bool nearMob = false;
    for (std::size_t i = 0; i < mob.size() && !nearMob; ++i) {
      const double dx = mob.xyz[i][0] - stam.xyz[j][0], dy = mob.xyz[i][1] - stam.xyz[j][1],
                   dz = mob.xyz[i][2] - stam.xyz[j][2];
      nearMob = dx * dx + dy * dy + dz * dz < reach * reach;
    }
    if (!nearMob) {
      shake[j] = 0;
      out.stam += stamAlone[j];
    }
  }
  const auto l = fn ? fn(all, shake) : shakeLnF(all, ep, p, shake);
  out.stam += std::accumulate(l.begin(), l.begin() + static_cast<long>(stam.size()), 0.0);
  out.mobm += std::accumulate(l.begin() + static_cast<long>(stam.size()), l.end(), 0.0);
  return out;
}

}  // namespace od
