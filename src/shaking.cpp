// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "shaking.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "ad4_pair.hpp"

namespace od {

ShakeAtoms shakeAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts, std::span<const Vec3f> coords) {
  ShakeAtoms s;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const Molecule& m = *parts[k].mol;
    const Ad4Typing& t = sys.typing(static_cast<int>(k));
    const int off = sys.coordOffset(static_cast<int>(k));
    const int n = static_cast<int>(m.atoms.size());
    std::vector<int> index(n, -1);
    for (int a = 0; a < n; ++a) {
      if (t.type[a] == Ad4Type::Merged) continue;
      index[a] = static_cast<int>(s.size());
      s.xyz.push_back(coords[off + a]);
      s.type.push_back(t.type[a]);
      s.charge.push_back(static_cast<float>(t.charge[a]));
      s.heavy.push_back(m.atoms[a].atomicNum > 1);
    }
    // Within three bonds over all bonds (as the 9a pair list, spec 4.5).
    std::vector<std::vector<int>> adj(n);
    for (const Bond& b : m.bonds) {
      adj[b.a].push_back(b.b);
      adj[b.b].push_back(b.a);
    }
    s.near.resize(s.size());
    std::vector<int> depth(n, -1);
    for (int a = 0; a < n; ++a) {
      if (index[a] < 0) continue;
      std::vector<int> q{a};
      depth[a] = 0;
      for (std::size_t i = 0; i < q.size(); ++i)
        if (depth[q[i]] < 3)
          for (int o : adj[q[i]])
            if (depth[o] < 0) {
              depth[o] = depth[q[i]] + 1;
              q.push_back(o);
            }
      auto& nr = s.near[index[a]];
      for (int o : q) {
        if (index[o] >= 0) nr.push_back(index[o]);
        depth[o] = -1;
      }
      std::sort(nr.begin(), nr.end());
    }
  }
  return s;
}

ShakeAtoms joinShakeAtoms(const ShakeAtoms& a, const ShakeAtoms& b) {
  ShakeAtoms s = a;
  const int off = static_cast<int>(a.size());
  s.xyz.insert(s.xyz.end(), b.xyz.begin(), b.xyz.end());
  s.type.insert(s.type.end(), b.type.begin(), b.type.end());
  s.charge.insert(s.charge.end(), b.charge.begin(), b.charge.end());
  s.heavy.insert(s.heavy.end(), b.heavy.begin(), b.heavy.end());
  for (const auto& nr : b.near) {
    s.near.push_back(nr);
    for (int& x : s.near.back()) x += off;
  }
  return s;
}

std::vector<Vec3f> shakePoints(double radius, double step) {
  std::vector<Vec3f> pts{{0, 0, 0}};
  const int n = static_cast<int>(std::floor(radius / step + 1e-9));
  const double r2 = radius * radius * (1 + 1e-9);
  for (int x = -n; x <= n; ++x)
    for (int y = -n; y <= n; ++y)
      for (int z = -n; z <= n; ++z) {
        if (x == 0 && y == 0 && z == 0) continue;
        const double dx = x * step, dy = y * step, dz = z * step;
        if (dx * dx + dy * dy + dz * dz <= r2)
          pts.push_back({static_cast<float>(dx), static_cast<float>(dy), static_cast<float>(dz)});
      }
  return pts;
}

std::vector<double> shakeLnF(const ShakeAtoms& atoms, const EnergyParams& ep, const ShakeParams& p,
                             std::span<const char> shake) {
  const std::size_t na = atoms.size();
  std::vector<double> lnF(na, 0.0);
  const std::vector<Vec3f> pts = shakePoints(p.radius, p.step);
  const Ad4Coeffs& c = ep.coeffs;
  const float half = static_cast<float>(0.5 * ep.smooth);
  const float rMin = static_cast<float>(ep.elecMinDistance);
  const double rrt = p.rrt();
  const double reach = kElecCutoff + p.radius;
  std::vector<float> S(na), V(na);
  for (std::size_t a = 0; a < na; ++a) {
    S[a] = static_cast<float>(ad4Params(atoms.type[a]).solpar + c.qasp * std::abs(atoms.charge[a]));
    V[a] = static_cast<float>(ad4Params(atoms.type[a]).volume);
  }
  struct Partner {
    Vec3f x;
    float AC, BD, reqm, q, sv;
    int m, n;
  };
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t i = 0; i < na; ++i) {
    if (!atoms.heavy[i] || (!shake.empty() && !shake[i])) continue;
    const Vec3f xi = atoms.xyz[i];
    std::vector<Partner> nb;
    const auto& ex = atoms.near[i];
    for (std::size_t j = 0; j < na; ++j) {
      if (j == i || std::binary_search(ex.begin(), ex.end(), static_cast<int>(j))) continue;
      const double dx = atoms.xyz[j][0] - xi[0], dy = atoms.xyz[j][1] - xi[1], dz = atoms.xyz[j][2] - xi[2];
      if (dx * dx + dy * dy + dz * dz >= reach * reach) continue;
      const Ad4PairParams q = ad4PairParams(atoms.type[i], atoms.type[j], c);
      nb.push_back({atoms.xyz[j], static_cast<float>(q.AC), static_cast<float>(q.BD), static_cast<float>(q.reqm),
                    static_cast<float>(c.elec) * atoms.charge[i] * atoms.charge[j],
                    static_cast<float>(c.desolv) * (S[i] * V[j] + S[j] * V[i]), q.m, q.n});
    }
    std::vector<double> e(pts.size());
    for (std::size_t k = 0; k < pts.size(); ++k) {
      const Vec3f x{xi[0] + pts[k][0], xi[1] + pts[k][1], xi[2] + pts[k][2]};
      double sum = 0;
      for (const Partner& b : nb) {
        const float dx = x[0] - b.x[0], dy = x[1] - b.x[1], dz = x[2] - b.x[2];
        const float r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= static_cast<float>(kElecCutoff * kElecCutoff)) continue;
        const float r = std::sqrt(r2);
        float t = 0, ee, de;
        if (r < static_cast<float>(kVdwCutoff)) {
          ad4VdwHb<float>(r, b.AC, b.BD, b.m, b.n, b.reqm, half, ee, de);
          t += ee;
        }
        ad4Elec<float>(r, b.q, rMin, ee, de);
        t += ee;
        ad4Desolv<float>(r, b.sv, rMin, ee, de);
        sum += t + ee;
      }
      e[k] = sum;
    }
    const double emin = *std::min_element(e.begin(), e.end());
    double w = 0;
    for (double v : e) w += std::exp(-(v - emin) * rrt);
    lnF[i] = std::log(w / static_cast<double>(pts.size()));
  }
  return lnF;
}

}  // namespace od
