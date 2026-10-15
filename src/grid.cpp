// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "grid.hpp"

#include <algorithm>
#include <cmath>

#include "ad4_pair.hpp"
#include "log.hpp"

namespace od {

namespace {

// Uniform cell list over atom positions.
class CellList {
 public:
  CellList(std::span<const Vec3> xyz, double cell) : cell_(cell) {
    if (xyz.empty()) return;
    lo_ = xyz[0];
    Vec3 hi = xyz[0];
    for (const Vec3& p : xyz)
      for (int c = 0; c < 3; ++c) {
        lo_[c] = std::min(lo_[c], p[c]);
        hi[c] = std::max(hi[c], p[c]);
      }
    for (int c = 0; c < 3; ++c) n_[c] = static_cast<int>((hi[c] - lo_[c]) / cell_) + 1;
    start_.assign(static_cast<std::size_t>(n_[0]) * n_[1] * n_[2] + 1, 0);
    std::vector<std::size_t> of(xyz.size());
    for (std::size_t a = 0; a < xyz.size(); ++a) {
      of[a] = flat(cellOf(xyz[a]));
      ++start_[of[a] + 1];
    }
    for (std::size_t c = 1; c < start_.size(); ++c) start_[c] += start_[c - 1];
    atoms_.resize(xyz.size());
    std::vector<int> fill(start_.begin(), start_.end() - 1);
    for (std::size_t a = 0; a < xyz.size(); ++a) atoms_[fill[of[a]]++] = static_cast<int>(a);
  }

  // Calls f(atom) for every atom in the cells overlapping the cube p ± r.
  template <class F>
  void forNear(const Vec3& p, double r, F&& f) const {
    if (atoms_.empty()) return;
    std::array<int, 3> a{}, b{};
    for (int c = 0; c < 3; ++c) {
      a[c] = std::max(0, static_cast<int>(std::floor((p[c] - r - lo_[c]) / cell_)));
      b[c] = std::min(n_[c] - 1, static_cast<int>(std::floor((p[c] + r - lo_[c]) / cell_)));
      if (a[c] > b[c]) return;
    }
    for (int z = a[2]; z <= b[2]; ++z)
      for (int y = a[1]; y <= b[1]; ++y) {
        const std::size_t row = flat({0, y, z});
        for (int k = start_[row + a[0]]; k < start_[row + b[0] + 1]; ++k) f(atoms_[k]);
      }
  }

 private:
  std::array<int, 3> cellOf(const Vec3& p) const {
    std::array<int, 3> c{};
    for (int k = 0; k < 3; ++k) c[k] = std::clamp(static_cast<int>((p[k] - lo_[k]) / cell_), 0, n_[k] - 1);
    return c;
  }
  std::size_t flat(const std::array<int, 3>& c) const {
    return (static_cast<std::size_t>(c[2]) * n_[1] + c[1]) * n_[0] + c[0];
  }
  double cell_;
  Vec3 lo_{};
  std::array<int, 3> n_{};
  std::vector<int> start_;
  std::vector<int> atoms_;
};

double dist2(const Vec3& a, const Vec3& b) {
  const double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return x * x + y * y + z * z;
}

}  // namespace

bool GridGeometry::inside(const Vec3& p) const {
  const Vec3 hi = max();
  for (int c = 0; c < 3; ++c)
    if (!(p[c] >= origin[c] && p[c] <= hi[c])) return false;
  return true;
}

GridGeometry boxGrid(const Vec3& center, const Vec3& size, double spacing) {
  GridGeometry g;
  g.spacing = spacing;
  for (int c = 0; c < 3; ++c) {
    int intervals = static_cast<int>(std::ceil(size[c] / spacing - 1e-9));
    intervals += intervals % 2;
    g.n[c] = intervals + 1;
    g.origin[c] = center[c] - 0.5 * intervals * spacing;
  }
  return g;
}

GridAtoms gridAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts, std::span<const Vec3> coords) {
  GridAtoms g;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const Molecule& m = *parts[k].mol;
    const Ad4Typing& t = sys.typing(static_cast<int>(k));
    const int off = sys.coordOffset(static_cast<int>(k));
    for (std::size_t a = 0; a < m.atoms.size(); ++a) {
      if (t.type[a] == Ad4Type::Merged) continue;
      g.xyz.push_back(coords[off + a]);
      g.type.push_back(t.type[a]);
      g.charge.push_back(t.charge[a]);
      g.heavy.push_back(m.atoms[a].atomicNum > 1);
    }
  }
  return g;
}

int EnergyGrids::mapOf(Ad4Type t) const {
  const auto it = std::find(types.begin(), types.end(), t);
  return it == types.end() ? -1 : static_cast<int>(it - types.begin());
}

EnergyGrids buildEnergyGrids(const GridGeometry& geom, const GridAtoms& stam, std::span<const Ad4Type> probeTypes,
                             const EnergyParams& params) {
  EnergyGrids g;
  g.geom = geom;
  for (Ad4Type t : probeTypes)
    if (t != Ad4Type::Merged && g.mapOf(t) < 0) g.types.push_back(t);
  const int nt = static_cast<int>(g.types.size());
  const int nTypes = static_cast<int>(Ad4Type::Count);
  const Ad4Coeffs& c = params.coeffs;
  // Per (probe type, StaM type): pair parameters.
  std::vector<Ad4PairParams> pp(static_cast<std::size_t>(nt) * nTypes);
  for (int i = 0; i < nt; ++i)
    for (int j = 0; j < nTypes; ++j) pp[i * nTypes + j] = ad4PairParams(g.types[i], static_cast<Ad4Type>(j), c);
  // Per StaM atom: S_j = solpar + qasp |q|, V_j.
  const std::size_t na = stam.size();
  std::vector<double> S(na), V(na);
  for (std::size_t a = 0; a < na; ++a) {
    S[a] = ad4Params(stam.type[a]).solpar + c.qasp * std::abs(stam.charge[a]);
    V[a] = ad4Params(stam.type[a]).volume;
  }
  std::vector<double> solT(nt), volT(nt);
  for (int i = 0; i < nt; ++i) {
    solT[i] = ad4Params(g.types[i]).solpar;
    volT[i] = ad4Params(g.types[i]).volume;
  }

  const std::size_t np = geom.size();
  g.typeMaps.assign(nt, std::vector<float>(np));
  g.elecMap.assign(np, 0.0f);
  g.desolvMap.assign(np, 0.0f);
  const CellList cells(stam.xyz, 4.0);
  const double half = 0.5 * params.smooth;
  const double rMin = params.elecMinDistance;
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
  for (int k = 0; k < geom.n[2]; ++k)
    for (int j = 0; j < geom.n[1]; ++j) {
      std::vector<double> e(nt);
      for (int i = 0; i < geom.n[0]; ++i) {
        const Vec3 p = geom.point(i, j, k);
        std::fill(e.begin(), e.end(), 0.0);
        double el = 0, ds = 0;
        cells.forNear(p, kElecCutoff, [&](int a) {
          const double r2 = dist2(p, stam.xyz[a]);
          if (r2 >= kElecCutoff * kElecCutoff) return;
          const double r = std::sqrt(r2);
          double ee, de;
          ad4Elec(r, c.elec * stam.charge[a], rMin, ee, de);
          el += ee;
          double gauss;
          ad4Desolv(r, 1.0, rMin, gauss, de);
          ds += c.desolv * c.qasp * V[a] * gauss;
          const int ta = static_cast<int>(stam.type[a]);
          for (int t = 0; t < nt; ++t) {
            e[t] += c.desolv * (solT[t] * V[a] + S[a] * volT[t]) * gauss;
            if (r < kVdwCutoff) {
              const Ad4PairParams& q = pp[t * nTypes + ta];
              ad4VdwHb(r, q.AC, q.BD, q.m, q.n, q.reqm, half, ee, de);
              e[t] += ee;
            }
          }
        });
        const std::size_t idx = geom.index(i, j, k);
        for (int t = 0; t < nt; ++t) g.typeMaps[t][idx] = static_cast<float>(std::min(e[t], kGridClamp));
        g.elecMap[idx] = static_cast<float>(el);
        g.desolvMap[idx] = static_cast<float>(ds);
      }
    }
  return g;
}

DistanceGrid buildDistanceGrid(const GridGeometry& geom, const GridAtoms& stam, double cap) {
  DistanceGrid g;
  g.geom = geom;
  g.cap = cap;
  g.d.assign(geom.size(), static_cast<float>(cap));
  std::vector<Vec3> heavy;
  for (std::size_t a = 0; a < stam.size(); ++a)
    if (stam.heavy[a]) heavy.push_back(stam.xyz[a]);
  const CellList cells(heavy, std::max(cap, 1.0));
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
  for (int k = 0; k < geom.n[2]; ++k)
    for (int j = 0; j < geom.n[1]; ++j)
      for (int i = 0; i < geom.n[0]; ++i) {
        const Vec3 p = geom.point(i, j, k);
        double best = cap * cap;
        cells.forNear(p, cap, [&](int a) { best = std::min(best, dist2(p, heavy[a])); });
        g.d[geom.index(i, j, k)] = static_cast<float>(std::sqrt(best));
      }
  return g;
}

namespace {

template <class T>
T interpolateT(const GridGeometry& geom, std::span<const float> map, const Vec3T<T>& p, Vec3T<T>* grad) {
  const T spacing = static_cast<T>(geom.spacing);
  std::array<int, 3> i0{};
  std::array<T, 3> f{};
  for (int c = 0; c < 3; ++c) {
    const T u = (p[c] - static_cast<T>(geom.origin[c])) / spacing;
    i0[c] = std::clamp(static_cast<int>(std::floor(u)), 0, geom.n[c] - 2);
    f[c] = u - i0[c];
  }
  T v[2][2][2];
  for (int z = 0; z < 2; ++z)
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 2; ++x) v[z][y][x] = map[geom.index(i0[0] + x, i0[1] + y, i0[2] + z)];
  const T fx = f[0], fy = f[1], fz = f[2];
  auto lerp = [](T a, T b, T t) { return a + (b - a) * t; };
  const T c00 = lerp(v[0][0][0], v[0][0][1], fx), c10 = lerp(v[0][1][0], v[0][1][1], fx);
  const T c01 = lerp(v[1][0][0], v[1][0][1], fx), c11 = lerp(v[1][1][0], v[1][1][1], fx);
  const T c0 = lerp(c00, c10, fy), c1 = lerp(c01, c11, fy);
  if (grad) {
    const T dx0 = lerp(v[0][0][1] - v[0][0][0], v[0][1][1] - v[0][1][0], fy);
    const T dx1 = lerp(v[1][0][1] - v[1][0][0], v[1][1][1] - v[1][1][0], fy);
    (*grad)[0] = lerp(dx0, dx1, fz) / spacing;
    (*grad)[1] = lerp(c10 - c00, c11 - c01, fz) / spacing;
    (*grad)[2] = (c1 - c0) / spacing;
  }
  return lerp(c0, c1, fz);
}

}  // namespace

double interpolate(const GridGeometry& geom, std::span<const float> map, const Vec3& p, Vec3* grad) {
  return interpolateT(geom, map, p, grad);
}

float interpolate(const GridGeometry& geom, std::span<const float> map, const Vec3f& p, Vec3f* grad) {
  return interpolateT(geom, map, p, grad);
}

InterEnergy gridEnergy(const EnergyGrids& g, const GridAtoms& probe, std::vector<Vec3>* grad) {
  InterEnergy r;
  if (grad) grad->assign(probe.size(), Vec3{});
  for (std::size_t a = 0; a < probe.size(); ++a) {
    const Vec3& p = probe.xyz[a];
    if (!g.geom.inside(p)) {
      ++r.outside;
      continue;
    }
    const int m = g.mapOf(probe.type[a]);
    if (m < 0) fail("no grid map for AD4 type {}", ad4TypeName(probe.type[a]));
    const double q = probe.charge[a];
    Vec3 gt{}, ge{}, gd{};
    r.energy += interpolate(g.geom, g.typeMaps[m], p, grad ? &gt : nullptr) +
                q * interpolate(g.geom, g.elecMap, p, grad ? &ge : nullptr) +
                std::abs(q) * interpolate(g.geom, g.desolvMap, p, grad ? &gd : nullptr);
    if (grad)
      for (int c = 0; c < 3; ++c) (*grad)[a][c] = gt[c] + q * ge[c] + std::abs(q) * gd[c];
  }
  return r;
}

double pairInterEnergy(const GridAtoms& stam, const GridAtoms& probe, const EnergyParams& params) {
  const Ad4Coeffs& c = params.coeffs;
  const double half = 0.5 * params.smooth;
  double sum = 0;
  for (std::size_t i = 0; i < probe.size(); ++i) {
    const Ad4TypeParams& pi = ad4Params(probe.type[i]);
    const double Si = pi.solpar + c.qasp * std::abs(probe.charge[i]);
    for (std::size_t j = 0; j < stam.size(); ++j) {
      const double r2 = dist2(probe.xyz[i], stam.xyz[j]);
      if (r2 >= kElecCutoff * kElecCutoff) continue;
      const double r = std::sqrt(r2);
      const Ad4TypeParams& pj = ad4Params(stam.type[j]);
      const double Sj = pj.solpar + c.qasp * std::abs(stam.charge[j]);
      double e = 0, ee, de;
      if (r < kVdwCutoff) {
        const Ad4PairParams q = ad4PairParams(probe.type[i], stam.type[j], c);
        ad4VdwHb(r, q.AC, q.BD, q.m, q.n, q.reqm, half, ee, de);
        e += ee;
      }
      ad4Elec(r, c.elec * probe.charge[i] * stam.charge[j], params.elecMinDistance, ee, de);
      e += ee;
      ad4Desolv(r, c.desolv * (Si * pj.volume + Sj * pi.volume), params.elecMinDistance, ee, de);
      sum += e + ee;
    }
  }
  return sum;
}

}  // namespace od
