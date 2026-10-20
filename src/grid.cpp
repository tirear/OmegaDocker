// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "grid.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>

#include "ad4_pair.hpp"
#include "log.hpp"

namespace od {

namespace {

double dist2(const Vec3& a, const Vec3& b) {
  const double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return x * x + y * y + z * z;
}

}  // namespace

namespace detail {

CellList::CellList(std::span<const Vec3> xyz, double cell) : cell_(cell) {
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
  auto flatOf = [&](const Vec3& p) {
    std::array<int, 3> c{};
    for (int k = 0; k < 3; ++k) c[k] = std::clamp(static_cast<int>((p[k] - lo_[k]) / cell_), 0, n_[k] - 1);
    return (static_cast<std::size_t>(c[2]) * n_[1] + c[1]) * n_[0] + c[0];
  };
  std::vector<std::size_t> of(xyz.size());
  for (std::size_t a = 0; a < xyz.size(); ++a) {
    of[a] = flatOf(xyz[a]);
    ++start_[of[a] + 1];
  }
  for (std::size_t c = 1; c < start_.size(); ++c) start_[c] += start_[c - 1];
  atoms_.resize(xyz.size());
  std::vector<int> fill(start_.begin(), start_.end() - 1);
  for (std::size_t a = 0; a < xyz.size(); ++a) atoms_[fill[of[a]]++] = static_cast<int>(a);
}

}  // namespace detail

using detail::CellList;

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

std::vector<Vec3> changedPositions(const GridAtoms& ref, const GridAtoms& cur) {
  if (ref.size() != cur.size()) fail("changedPositions: different atoms");
  std::vector<Vec3> out;
  for (std::size_t a = 0; a < cur.size(); ++a)
    if (cur.xyz[a] != ref.xyz[a]) {
      out.push_back(ref.xyz[a]);
      out.push_back(cur.xyz[a]);
    }
  return out;
}

int EnergyGrids::mapOf(Ad4Type t) const {
  const auto it = std::find(types.begin(), types.end(), t);
  return it == types.end() ? -1 : static_cast<int>(it - types.begin());
}

GridBuildData::GridBuildData(const GridGeometry& g, const GridAtoms& stam, std::span<const Ad4Type> probeTypes,
                             const EnergyParams& params)
    : geom(g) {
  for (Ad4Type t : probeTypes)
    if (t != Ad4Type::Merged && std::find(types.begin(), types.end(), t) == types.end()) types.push_back(t);
  const int nt = static_cast<int>(types.size());
  if (nt > kGridMaxTypeMaps) fail("energy grids: more than {} probe types", kGridMaxTypeMaps);
  const int nTypes = static_cast<int>(Ad4Type::Count);
  const Ad4Coeffs& c = params.coeffs;
  // Per (probe type, StaM type): pair parameters.
  for (int i = 0; i < nt; ++i)
    for (int j = 0; j < nTypes; ++j) {
      const Ad4PairParams q = ad4PairParams(types[i], static_cast<Ad4Type>(j), c);
      AC.push_back(q.AC);
      BD.push_back(q.BD);
      reqm.push_back(q.reqm);
      pm.push_back(q.m);
      pn.push_back(q.n);
    }
  // Per StaM atom: S_j = solpar + qasp |q|, V_j.
  for (std::size_t a = 0; a < stam.size(); ++a) {
    for (int k = 0; k < 3; ++k) xyz.push_back(stam.xyz[a][k]);
    type.push_back(static_cast<int>(stam.type[a]));
    elecQ.push_back(c.elec * stam.charge[a]);
    S.push_back(ad4Params(stam.type[a]).solpar + c.qasp * std::abs(stam.charge[a]));
    V.push_back(ad4Params(stam.type[a]).volume);
  }
  for (int i = 0; i < nt; ++i) {
    solT.push_back(ad4Params(types[i]).solpar);
    volT.push_back(ad4Params(types[i]).volume);
  }
  for (std::size_t a = 0; a < stam.size(); ++a) {
    dsW.push_back(c.desolv * c.qasp * V[a]);
    for (int i = 0; i < nt; ++i) desolvW.push_back(c.desolv * (solT[i] * V[a] + S[a] * volT[i]));
  }
  cells = CellList(stam.xyz, 4.0);
  GridBuildView& v = view;
  for (int k = 0; k < 3; ++k) {
    v.origin[k] = geom.origin[k];
    v.n[k] = geom.n[k];
  }
  v.spacing = geom.spacing;
  v.typeMapCount = nt;
  v.nTypes = nTypes;
  v.xyz = xyz.data();
  v.type = type.data();
  v.elecQ = elecQ.data();
  v.S = S.data();
  v.V = V.data();
  v.solT = solT.data();
  v.volT = volT.data();
  v.desolvW = desolvW.data();
  v.dsW = dsW.data();
  v.AC = AC.data();
  v.BD = BD.data();
  v.reqm = reqm.data();
  v.pm = pm.data();
  v.pn = pn.data();
  v.cDesolv = c.desolv;
  v.cDesolvQasp = c.desolv * c.qasp;
  v.halfSmooth = 0.5 * params.smooth;
  v.rMin = params.elecMinDistance;
  v.cells = cells.flat();
}

EnergyGrids buildEnergyGrids(const GridGeometry& geom, const GridAtoms& stam, std::span<const Ad4Type> probeTypes,
                             const EnergyParams& params) {
  const GridBuildData b(geom, stam, probeTypes, params);
  EnergyGrids g;
  g.geom = geom;
  g.types = b.types;
  const int nt = static_cast<int>(g.types.size());
  const std::size_t np = geom.size();
  g.typeMaps.assign(nt, std::vector<float>(np));
  g.elecMap.assign(np, 0.0f);
  g.desolvMap.assign(np, 0.0f);
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
  for (int k = 0; k < geom.n[2]; ++k)
    for (int j = 0; j < geom.n[1]; ++j) {
      double e[kGridMaxTypeMaps];
      for (int i = 0; i < geom.n[0]; ++i) {
        double el, ds;
        gridPoint(b.view, i, j, k, e, el, ds);
        const std::size_t idx = geom.index(i, j, k);
        for (int t = 0; t < nt; ++t) g.typeMaps[t][idx] = static_cast<float>(std::min(e[t], kGridClamp));
        g.elecMap[idx] = static_cast<float>(el);
        g.desolvMap[idx] = static_cast<float>(ds);
      }
    }
  return g;
}

DistBuildData::DistBuildData(const GridGeometry& geom, const GridAtoms& stam, double cap) {
  std::vector<Vec3> heavy;
  for (std::size_t a = 0; a < stam.size(); ++a)
    if (stam.heavy[a]) {
      heavy.push_back(stam.xyz[a]);
      for (int k = 0; k < 3; ++k) xyz.push_back(stam.xyz[a][k]);
    }
  cells = CellList(heavy, std::max(cap, 1.0));
  for (int k = 0; k < 3; ++k) {
    view.origin[k] = geom.origin[k];
    view.n[k] = geom.n[k];
  }
  view.spacing = geom.spacing;
  view.cap = cap;
  view.xyz = xyz.data();
  view.cells = cells.flat();
}

DistanceGrid buildDistanceGrid(const GridGeometry& geom, const GridAtoms& stam, double cap) {
  DistanceGrid g;
  g.geom = geom;
  g.cap = cap;
  g.d.assign(geom.size(), static_cast<float>(cap));
  const DistBuildData b(geom, stam, cap);
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
  for (int k = 0; k < geom.n[2]; ++k)
    for (int j = 0; j < geom.n[1]; ++j)
      for (int i = 0; i < geom.n[0]; ++i) g.d[geom.index(i, j, k)] = distPoint(b.view, i, j, k);
  return g;
}

namespace {

// at(i, j, k): map value at a grid point.
template <class T, class At>
T interpolateAt(const GridGeometry& geom, const At& at, const Vec3T<T>& p, Vec3T<T>* grad) {
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
      for (int x = 0; x < 2; ++x) v[z][y][x] = at(i0[0] + x, i0[1] + y, i0[2] + z);
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

template <class T>
T interpolateT(const GridGeometry& geom, std::span<const float> map, const Vec3T<T>& p, Vec3T<T>* grad) {
  return interpolateAt(geom, [&](int i, int j, int k) { return map[geom.index(i, j, k)]; }, p, grad);
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

GridAtoms subset(const GridAtoms& atoms, std::span<const char> mask, bool keep) {
  GridAtoms g;
  for (std::size_t a = 0; a < atoms.size(); ++a) {
    if ((mask[a] != 0) != keep) continue;
    g.xyz.push_back(atoms.xyz[a]);
    g.type.push_back(atoms.type[a]);
    g.charge.push_back(atoms.charge[a]);
    g.heavy.push_back(atoms.heavy[a]);
  }
  return g;
}

std::vector<char> movedGridAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts) {
  std::vector<char> mask;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const Molecule& m = *parts[k].mol;
    const TorsionTree& tree = *parts[k].tree;
    std::vector<char> moved(m.atoms.size(), 0);
    for (const TreeTorsion& tt : tree.torsions)
      for (int c : tt.moving)
        if (c < static_cast<int>(m.atoms.size())) moved[c] = 1;
    const Ad4Typing& t = sys.typing(static_cast<int>(k));
    for (std::size_t a = 0; a < m.atoms.size(); ++a)
      if (t.type[a] != Ad4Type::Merged) mask.push_back(moved[a]);
  }
  return mask;
}

const char* gridModeName(GridMode m) {
  switch (m) {
    case GridMode::Full: return "full";
    case GridMode::Flexres: return "flexres";
    case GridMode::Delta: return "delta";
  }
  return "";
}

FieldCells detail::CellList::flat() const {
  FieldCells c{};
  c.cell = cell_;
  for (int k = 0; k < 3; ++k) {
    c.lo[k] = lo_[k];
    c.n[k] = n_[k];
  }
  c.start = start_.data();
  c.atoms = atoms_.data();
  c.atomCount = static_cast<int>(atoms_.size());
  return c;
}

StamField StamField::full(std::shared_ptr<const EnergyGrids> grids, std::shared_ptr<const DistanceGrid> dist) {
  static std::atomic<std::uint64_t> counter{0};
  StamField f;
  f.mode_ = GridMode::Full;
  f.id_ = ++counter;
  f.grids_ = std::move(grids);
  f.dist_ = std::move(dist);
  f.mapCount_ = static_cast<int>(f.grids_->types.size()) + 3;
  return f;
}

StamField StamField::flexres(std::shared_ptr<const EnergyGrids> fixed, std::shared_ptr<const DistanceGrid> fixedDist,
                             GridAtoms moved, const EnergyParams& params) {
  StamField f = full(std::move(fixed), std::move(fixedDist));
  f.mode_ = GridMode::Flexres;
  const Ad4Coeffs& c = params.coeffs;
  f.cDesolv_ = c.desolv;
  f.qasp_ = c.qasp;
  f.halfSmooth_ = 0.5 * params.smooth;
  f.rMin_ = params.elecMinDistance;
  const auto& types = f.grids_->types;
  const int nTypes = static_cast<int>(Ad4Type::Count);
  for (Ad4Type t : types) {
    f.solT_.push_back(ad4Params(t).solpar);
    f.volT_.push_back(ad4Params(t).volume);
    for (int j = 0; j < nTypes; ++j) {
      const Ad4PairParams q = ad4PairParams(t, static_cast<Ad4Type>(j), c);
      f.pairAC_.push_back(q.AC);
      f.pairBD_.push_back(q.BD);
      f.pairReqm_.push_back(q.reqm);
      f.pairM_.push_back(q.m);
      f.pairN_.push_back(q.n);
    }
  }
  std::vector<Vec3> xyz, heavyXyz;
  for (std::size_t a = 0; a < moved.size(); ++a) {
    Moved m;
    m.xyz = moved.xyz[a];
    m.type = static_cast<int>(moved.type[a]);
    m.elecQ = c.elec * moved.charge[a];
    m.S = ad4Params(moved.type[a]).solpar + c.qasp * std::abs(moved.charge[a]);
    m.V = ad4Params(moved.type[a]).volume;
    m.heavy = moved.heavy[a] != 0;
    if (m.heavy) {
      f.movedHeavy_.push_back(static_cast<int>(a));
      heavyXyz.push_back(m.xyz);
    }
    xyz.push_back(m.xyz);
    f.moved_.push_back(m);
  }
  f.cells_ = CellList(xyz, 4.0);
  f.heavyCells_ = CellList(heavyXyz, std::max(f.dist_->cap, 1.0));
  return f;
}

StamField StamField::delta(std::shared_ptr<const EnergyGrids> ref, std::shared_ptr<const DistanceGrid> refDist,
                           const EnergyGrids& grids, const DistanceGrid& dist, double threshold,
                           const std::vector<Vec3>* changed) {
  StamField f = full(std::move(ref), std::move(refDist));
  f.mode_ = GridMode::Delta;
  const GridGeometry& g = f.grids_->geom;
  const int nt = static_cast<int>(f.grids_->types.size());
  if (grids.types != f.grids_->types || grids.geom.n != g.n || dist.geom.n != g.n)
    fail("delta grids: the sample's grids differ in layout from the reference");
  auto sampleMap = [&](int m) -> std::span<const float> {
    if (m < nt) return grids.typeMaps[m];
    if (m == nt) return grids.elecMap;
    if (m == nt + 1) return grids.desolvMap;
    return dist.d;
  };
  auto refMap = [&](int m) -> std::span<const float> {
    if (m < nt) return f.grids_->typeMaps[m];
    if (m == nt) return f.grids_->elecMap;
    if (m == nt + 1) return f.grids_->desolvMap;
    return f.dist_->d;
  };
  constexpr int B = kDeltaBlock;
  for (int c = 0; c < 3; ++c) f.nb_[c] = (g.n[c] + B - 1) / B;
  f.block_.assign(static_cast<std::size_t>(f.nb_[0]) * f.nb_[1] * f.nb_[2], -1);
  const std::size_t bs = static_cast<std::size_t>(B) * B * B;
  const double reach2 = deltaReach2(g.spacing, std::max(kElecCutoff, dist.cap));
  const double origin[3] = {g.origin[0], g.origin[1], g.origin[2]};
  std::vector<float> data;
  for (int bz = 0; bz < f.nb_[2]; ++bz)
    for (int by = 0; by < f.nb_[1]; ++by)
      for (int bx = 0; bx < f.nb_[0]; ++bx) {
        if (changed && std::none_of(changed->begin(), changed->end(), [&](const Vec3& p) {
              return deltaBlockReached(origin, g.spacing, bx, by, bz, p.data(), reach2);
            }))
          continue;
        bool differs = false;
        for (int m = 0; m < f.mapCount_ && !differs; ++m) {
          const auto s = sampleMap(m), r = refMap(m);
          for (int z = bz * B; z < std::min(g.n[2], bz * B + B) && !differs; ++z)
            for (int y = by * B; y < std::min(g.n[1], by * B + B) && !differs; ++y)
              for (int x = bx * B; x < std::min(g.n[0], bx * B + B); ++x) {
                const std::size_t i = g.index(x, y, z);
                if (std::abs(static_cast<double>(s[i]) - r[i]) > threshold) {
                  differs = true;
                  break;
                }
              }
        }
        if (!differs) continue;
        const std::size_t b = (static_cast<std::size_t>(bz) * f.nb_[1] + by) * f.nb_[0] + bx;
        f.block_[b] = static_cast<int>(data.size() / (bs * f.mapCount_));
        const std::size_t base = data.size();
        data.resize(base + bs * f.mapCount_, 0.0f);
        for (int m = 0; m < f.mapCount_; ++m) {
          const auto s = sampleMap(m);
          for (int z = bz * B; z < std::min(g.n[2], bz * B + B); ++z)
            for (int y = by * B; y < std::min(g.n[1], by * B + B); ++y)
              for (int x = bx * B; x < std::min(g.n[0], bx * B + B); ++x)
                data[base + m * bs + ((z - bz * B) * B + (y - by * B)) * B + (x - bx * B)] = s[g.index(x, y, z)];
        }
      }
  f.data_ = std::make_shared<const std::vector<float>>(std::move(data));
  return f;
}

StamField StamField::delta(std::shared_ptr<const EnergyGrids> ref, std::shared_ptr<const DistanceGrid> refDist,
                           std::vector<int> block, std::vector<float> data,
                           std::shared_ptr<const DeviceBlocks> device) {
  StamField f = full(std::move(ref), std::move(refDist));
  f.mode_ = GridMode::Delta;
  const GridGeometry& g = f.grids_->geom;
  constexpr int B = kDeltaBlock;
  for (int c = 0; c < 3; ++c) f.nb_[c] = (g.n[c] + B - 1) / B;
  const std::size_t bs = static_cast<std::size_t>(B) * B * B;
  const std::size_t stored = static_cast<std::size_t>(std::count_if(block.begin(), block.end(), [](int b) { return b >= 0; }));
  const std::size_t floats = stored * bs * f.mapCount_;
  if (block.size() != static_cast<std::size_t>(f.nb_[0]) * f.nb_[1] * f.nb_[2] ||
      (device ? device->floats != floats || (!data.empty() && data.size() != floats) : data.size() != floats))
    fail("delta grids: block layout differs from the reference");
  f.block_ = std::move(block);
  if (!device || !data.empty()) f.data_ = std::make_shared<const std::vector<float>>(std::move(data));
  f.device_ = std::move(device);
  return f;
}

FieldModel StamField::model() const {
  FieldModel fm;
  fm.id = id_;
  fm.base = grids_;
  fm.baseDist = dist_;
  const GridGeometry& g = grids_->geom;
  fm.points = g.size();
  FieldView& v = fm.view;
  v.mode = mode_ == GridMode::Full ? 0 : mode_ == GridMode::Flexres ? 1 : 2;
  for (int c = 0; c < 3; ++c) {
    v.origin[c] = g.origin[c];
    v.n[c] = g.n[c];
    v.nb[c] = nb_[c];
  }
  v.spacing = g.spacing;
  v.typeMapCount = static_cast<int>(grids_->types.size());
  if (v.typeMapCount + 3 > kFieldMaxMaps) fail("StamField::model: too many maps");
  for (int m = 0; m < v.typeMapCount; ++m) v.maps[m] = grids_->typeMaps[m].data();
  v.maps[v.typeMapCount] = grids_->elecMap.data();
  v.maps[v.typeMapCount + 1] = grids_->desolvMap.data();
  v.maps[v.typeMapCount + 2] = dist_->d.data();
  v.mapCount = mapCount_;
  fm.block = block_;
  fm.data = data_;
  fm.deviceData = device_;
  for (const Moved& m : moved_) {
    FieldMoved x{};
    for (int c = 0; c < 3; ++c) x.xyz[c] = m.xyz[c];
    x.type = m.type;
    x.elecQ = m.elecQ;
    x.S = m.S;
    x.V = m.V;
    x.heavy = m.heavy;
    fm.moved.push_back(x);
  }
  fm.movedHeavy = movedHeavy_;
  v.movedCount = static_cast<int>(moved_.size());
  v.movedHeavyCount = static_cast<int>(movedHeavy_.size());
  v.nTypes = static_cast<int>(Ad4Type::Count);
  v.cells = cells_.flat();
  v.heavyCells = heavyCells_.flat();
  fm.cellStart = cells_.start();
  fm.cellAtoms = cells_.atoms();
  fm.heavyStart = heavyCells_.start();
  fm.heavyAtoms = heavyCells_.atoms();
  fm.solT = solT_;
  fm.volT = volT_;
  fm.pairAC = pairAC_;
  fm.pairBD = pairBD_;
  fm.pairReqm = pairReqm_;
  fm.pairM = pairM_;
  fm.pairN = pairN_;
  v.cDesolv = cDesolv_;
  v.qasp = qasp_;
  v.halfSmooth = halfSmooth_;
  v.rMin = rMin_;
  v.cap = dist_->cap;
  fm.bind();
  return fm;
}

std::size_t StamField::deltaBlockCount() const {
  return static_cast<std::size_t>(std::count_if(block_.begin(), block_.end(), [](int b) { return b >= 0; }));
}

template <class T>
T StamField::map(int m, const Vec3T<T>& p, Vec3T<T>* grad) const {
  const int nt = static_cast<int>(grids_->types.size());
  const std::span<const float> r = m < nt        ? std::span<const float>(grids_->typeMaps[m])
                                   : m == nt     ? std::span<const float>(grids_->elecMap)
                                   : m == nt + 1 ? std::span<const float>(grids_->desolvMap)
                                                 : std::span<const float>(dist_->d);
  const GridGeometry& g = grids_->geom;
  if (mode_ != GridMode::Delta) return interpolateT(g, r, p, grad);
  constexpr int B = kDeltaBlock;
  const std::size_t bs = static_cast<std::size_t>(B) * B * B;
  const float* data = data_ ? data_->data() : nullptr;
  return interpolateAt(
      g,
      [&](int i, int j, int k) {
        const int b = block_[(static_cast<std::size_t>(k / B) * nb_[1] + j / B) * nb_[0] + i / B];
        if (b < 0) return r[g.index(i, j, k)];
        if (!data) fail("delta grids: the blocks are only in device memory");
        return data[(static_cast<std::size_t>(b) * mapCount_ + m) * bs + ((k % B) * B + j % B) * B + i % B];
      },
      p, grad);
}

template <class T>
T StamField::energy(int m, T q, const Vec3T<T>& p, Vec3T<T>* grad) const {
  using V = Vec3T<T>;
  const int nt = static_cast<int>(grids_->types.size());
  V a{}, b{}, s{};
  const T aq = std::abs(q);
  const T grid = map(m, p, grad ? &a : nullptr) + q * map(nt, p, grad ? &b : nullptr) + aq * map(nt + 1, p, grad ? &s : nullptr);
  if (grad)
    for (int i = 0; i < 3; ++i) (*grad)[i] = a[i] + q * b[i] + aq * s[i];
  if (mode_ != GridMode::Flexres || moved_.empty()) return grid;
  return static_cast<T>(grid + movedSum(m, q, p, grad));
}

template <class T>
T StamField::movedSum(int m, T q, const Vec3T<T>& p, Vec3T<T>* grad) const {
  using V = Vec3T<T>;
  const T aq = std::abs(q);
  // Pair sum over the moved StaM atoms: the grid terms of spec 8.1 per atom (no clamp).
  const int nTypes = static_cast<int>(Ad4Type::Count);
  const T rMin = static_cast<T>(rMin_), half = static_cast<T>(halfSmooth_);
  const T cut2 = static_cast<T>(kElecCutoff * kElecCutoff), vdwCut = static_cast<T>(kVdwCutoff);
  const T solT = static_cast<T>(solT_[m]), volT = static_cast<T>(volT_[m]);
  const T cd = static_cast<T>(cDesolv_), cdq = static_cast<T>(cDesolv_ * qasp_);
  double sum = 0;
  Vec3 gs{};
  cells_.forNear({static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2])}, kElecCutoff,
                 [&](int j) {
                   const Moved& at = moved_[j];
                   const V d{p[0] - static_cast<T>(at.xyz[0]), p[1] - static_cast<T>(at.xyz[1]),
                             p[2] - static_cast<T>(at.xyz[2])};
                   const T r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                   if (r2 >= cut2) return;
                   const T r = std::sqrt(r2);
                   T gauss, dgauss, el, del;
                   ad4Desolv(r, T(1), rMin, gauss, dgauss);
                   ad4Elec(r, static_cast<T>(at.elecQ), rMin, el, del);
                   const T V_ = static_cast<T>(at.V), S_ = static_cast<T>(at.S);
                   const T w = cd * (solT * V_ + S_ * volT) + aq * cdq * V_;
                   T e = w * gauss + q * el, dedr = w * dgauss + q * del;
                   if (r < vdwCut) {
                     const std::size_t k = static_cast<std::size_t>(m) * nTypes + at.type;
                     T ev, dv;
                     ad4VdwHb(r, static_cast<T>(pairAC_[k]), static_cast<T>(pairBD_[k]), pairM_[k], pairN_[k],
                              static_cast<T>(pairReqm_[k]), half, ev, dv);
                     e += ev;
                     dedr += dv;
                   }
                   sum += e;
                   if (grad && r > 0)
                     for (int i = 0; i < 3; ++i) gs[i] += dedr / r * d[i];
                 });
  if (grad)
    for (int i = 0; i < 3; ++i) (*grad)[i] = static_cast<T>((*grad)[i] + gs[i]);
  return static_cast<T>(sum);
}

std::vector<float> StamField::pointMap(int m) const {
  const GridGeometry& g = grids_->geom;
  const int nt = typeMapCount();
  if (m < 0 || m > nt + 1) fail("pointMap: map {} out of range", m);
  const std::span<const float> r = m < nt ? std::span<const float>(grids_->typeMaps[m])
                                   : m == nt ? std::span<const float>(grids_->elecMap)
                                             : std::span<const float>(grids_->desolvMap);
  std::vector<float> out(r.begin(), r.end());
  if (mode_ == GridMode::Delta) {
    if (!data_ && deltaBlockCount() > 0) fail("delta grids: the blocks are only in device memory");
    constexpr int B = kDeltaBlock;
    const std::size_t bs = static_cast<std::size_t>(B) * B * B;
    for (int k = 0; k < g.n[2]; ++k)
      for (int j = 0; j < g.n[1]; ++j)
        for (int i = 0; i < g.n[0]; ++i) {
          const int b = block_[(static_cast<std::size_t>(k / B) * nb_[1] + j / B) * nb_[0] + i / B];
          if (b >= 0)
            out[g.index(i, j, k)] = (*data_)[(static_cast<std::size_t>(b) * mapCount_ + m) * bs + ((k % B) * B + j % B) * B + i % B];
        }
  }
  if (mode_ == GridMode::Flexres && !moved_.empty()) {
#pragma omp parallel for schedule(static)
    for (int k = 0; k < g.n[2]; ++k)
      for (int j = 0; j < g.n[1]; ++j)
        for (int i = 0; i < g.n[0]; ++i) {
          const Vec3 p = g.point(i, j, k);
          double v;
          if (m < nt) {
            v = movedSum<double>(m, 0.0, p, nullptr);
          } else {
            const double ep = movedSum<double>(0, 1.0, p, nullptr), en = movedSum<double>(0, -1.0, p, nullptr);
            v = m == nt ? 0.5 * (ep - en) : 0.5 * (ep + en) - movedSum<double>(0, 0.0, p, nullptr);
          }
          out[g.index(i, j, k)] += static_cast<float>(v);
        }
  }
  return out;
}

template <class T>
T StamField::distance(const Vec3T<T>& p, Vec3T<T>* grad) const {
  const int nt = static_cast<int>(grids_->types.size());
  T d = map(nt + 2, p, grad);
  if (mode_ != GridMode::Flexres || movedHeavy_.empty()) return d;
  const double cap = dist_->cap;
  int best = -1;
  T best2 = d * d;
  heavyCells_.forNear({static_cast<double>(p[0]), static_cast<double>(p[1]), static_cast<double>(p[2])}, cap, [&](int h) {
    const Vec3& x = moved_[movedHeavy_[h]].xyz;
    const T dx = p[0] - static_cast<T>(x[0]), dy = p[1] - static_cast<T>(x[1]), dz = p[2] - static_cast<T>(x[2]);
    const T r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < best2) {
      best2 = r2;
      best = movedHeavy_[h];
    }
  });
  if (best < 0) return d;
  const T r = std::sqrt(best2);
  if (grad) {
    const Vec3& x = moved_[best].xyz;
    for (int i = 0; i < 3; ++i) (*grad)[i] = r > 0 ? (p[i] - static_cast<T>(x[i])) / r : T(0);
  }
  return r;
}

template float StamField::energy<float>(int, float, const Vec3f&, Vec3f*) const;
template double StamField::energy<double>(int, double, const Vec3&, Vec3*) const;
template float StamField::distance<float>(const Vec3f&, Vec3f*) const;
template double StamField::distance<double>(const Vec3&, Vec3*) const;

}  // namespace od
