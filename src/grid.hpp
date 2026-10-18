// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

#include "ad4.hpp"
#include "energy.hpp"
#include "field_kernel.hpp"
#include "grid_kernel.hpp"
#include "molecule.hpp"

namespace od {

// StaM–MobM intermolecular energy grids (spec 8.1) and the distance grid
// (spec 8.2). CPU reference implementation.

constexpr double kGridClamp = 1.0e5;  // kcal/mol, upper bound of a type-map value (AutoGrid EINTCLAMP)

struct GridGeometry {
  Vec3 origin{};                // coordinates of point (0, 0, 0)
  double spacing = 0.375;       // Å
  std::array<int, 3> n{};       // points per axis
  std::size_t size() const { return static_cast<std::size_t>(n[0]) * n[1] * n[2]; }
  std::size_t index(int i, int j, int k) const {
    return (static_cast<std::size_t>(k) * n[1] + j) * n[0] + i;
  }
  Vec3 point(int i, int j, int k) const {
    return {origin[0] + i * spacing, origin[1] + j * spacing, origin[2] + k * spacing};
  }
  Vec3 max() const { return point(n[0] - 1, n[1] - 1, n[2] - 1); }
  bool inside(const Vec3& p) const;
};

// Grid covering the docking box: an even number of intervals per axis
// (at least size / spacing), centred on the box centre.
GridGeometry boxGrid(const Vec3& center, const Vec3& size, double spacing);

// Atoms that generate or probe a grid: real atoms of the parts (merged
// hydrogens and virtual closure atoms excluded).
struct GridAtoms {
  std::vector<Vec3> xyz;
  std::vector<Ad4Type> type;
  std::vector<double> charge;
  std::vector<char> heavy;
  std::size_t size() const { return xyz.size(); }
};
// coords: coordinates of sys.evaluate for the same parts.
GridAtoms gridAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts, std::span<const Vec3> coords);

// For probe atom i at p: E = typeMap[type](p) + q * elecMap(p) + |q| * desolvMap(p).
struct EnergyGrids {
  GridGeometry geom;
  std::vector<Ad4Type> types;  // one type map per probe type
  std::vector<std::vector<float>> typeMaps;
  std::vector<float> elecMap;
  std::vector<float> desolvMap;
  int mapOf(Ad4Type t) const;  // -1 when the type has no map
};

// Probe types: the union of the MobM types (spec 11). Same pair functions and
// cutoffs as the 9a energy; no H-bond directionality.
EnergyGrids buildEnergyGrids(const GridGeometry& geom, const GridAtoms& stam, std::span<const Ad4Type> probeTypes,
                             const EnergyParams& params);

// Each point holds the distance to the nearest StaM heavy atom, at most cap.
struct DistanceGrid {
  GridGeometry geom;
  std::vector<float> d;
  double cap = 0;
};
DistanceGrid buildDistanceGrid(const GridGeometry& geom, const GridAtoms& stam, double cap);

// Trilinear interpolation of a map at p (p inside geom); grad (optional): dvalue/dp.
double interpolate(const GridGeometry& geom, std::span<const float> map, const Vec3& p, Vec3* grad = nullptr);
// Single precision (9b, spec 8.6): all arithmetic in float.
float interpolate(const GridGeometry& geom, std::span<const float> map, const Vec3f& p, Vec3f* grad = nullptr);

struct InterEnergy {
  double energy = 0;
  int outside = 0;  // probe atoms outside the grid (not counted in energy)
};
// Grid energy of the probe atoms; grad (optional, per probe atom, kcal/mol/Å).
InterEnergy gridEnergy(const EnergyGrids& g, const GridAtoms& probe, std::vector<Vec3>* grad = nullptr);

// Reference: direct sum over StaM–probe atom pairs with the grid's terms (no clamp).
double pairInterEnergy(const GridAtoms& stam, const GridAtoms& probe, const EnergyParams& params);

// Atoms of a GridAtoms list with mask[a] == keep.
GridAtoms subset(const GridAtoms& atoms, std::span<const char> mask, bool keep);
// Per atom of gridAtoms(sys, parts, ...): 1 when a torsion of the parts' trees moves it.
std::vector<char> movedGridAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts);

namespace detail {

// Uniform cell list over atom positions.
class CellList {
 public:
  CellList() = default;
  CellList(std::span<const Vec3> xyz, double cell);
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
        const std::size_t row = (static_cast<std::size_t>(z) * n_[1] + y) * n_[0];
        for (int k = start_[row + a[0]]; k < start_[row + b[0] + 1]; ++k) f(atoms_[k]);
      }
  }
  // Flat form (field::forNear); the arrays point into this list.
  FieldCells flat() const;
  const std::vector<int>& start() const { return start_; }
  const std::vector<int>& atoms() const { return atoms_; }

 private:
  double cell_ = 1;
  Vec3 lo_{};
  std::array<int, 3> n_{};
  std::vector<int> start_;
  std::vector<int> atoms_;
};

}  // namespace detail

// Inputs of gridPoint for buildEnergyGrids on the CPU or the GPU (spec 15.1.10); view points
// into the members.
struct GridBuildData {
  GridBuildData(const GridGeometry& geom, const GridAtoms& stam, std::span<const Ad4Type> probeTypes,
                const EnergyParams& params);
  GridBuildData(const GridBuildData&) = delete;
  GridBuildData& operator=(const GridBuildData&) = delete;

  GridGeometry geom;
  std::vector<Ad4Type> types;  // EnergyGrids::types
  GridBuildView view{};
  std::vector<double> xyz, elecQ, S, V, solT, volT, AC, BD, reqm;
  std::vector<int> type, pm, pn;
  detail::CellList cells;
};

// Inputs of distPoint for buildDistanceGrid (spec 15.1.10).
struct DistBuildData {
  DistBuildData(const GridGeometry& geom, const GridAtoms& stam, double cap);
  DistBuildData(const DistBuildData&) = delete;
  DistBuildData& operator=(const DistBuildData&) = delete;

  DistBuildView view{};
  std::vector<double> xyz;
  detail::CellList cells;
};

// gridMode (spec 8.1).
enum class GridMode : std::uint8_t { Full, Flexres, Delta };
const char* gridModeName(GridMode m);  // "full" | "flexres" | "delta"

// The StaM side of the 9b energy for one StaM conformation sample (spec 8.1, 8.2):
// full: energy and distance grids of the sample; flexres: grids of the fixed atoms
// plus a pair sum over the moved atoms; delta: the grids of a reference conformation
// with the blocks that differ replaced by the sample's values.
class StamField {
 public:
  StamField() = default;
  static StamField full(std::shared_ptr<const EnergyGrids> grids, std::shared_ptr<const DistanceGrid> dist);
  // fixed: grids of the fixed atoms; moved: the other StaM grid atoms of the sample.
  static StamField flexres(std::shared_ptr<const EnergyGrids> fixed, std::shared_ptr<const DistanceGrid> fixedDist,
                           GridAtoms moved, const EnergyParams& params);
  // Blocks with a point differing from the reference by more than threshold in any map
  // (type maps, elec, desolv, distance) store the sample's values.
  static StamField delta(std::shared_ptr<const EnergyGrids> ref, std::shared_ptr<const DistanceGrid> refDist,
                         const EnergyGrids& grids, const DistanceGrid& dist, double threshold);
  // delta from its blocks (e.g. computed on the GPU, spec 15.1.10): block per block (bz, by, bx
  // order) -1 or the index of its values in data; data: per stored block, the maps (type maps,
  // elec, desolv, distance) × kDeltaBlock^3 values (z, y, x order; 0 outside the grid).
  static StamField delta(std::shared_ptr<const EnergyGrids> ref, std::shared_ptr<const DistanceGrid> refDist,
                         std::vector<int> block, std::vector<float> data);

  GridMode mode() const { return mode_; }
  // Unique per StamField created by full, flexres or delta (copies share it).
  std::uint64_t id() const { return id_; }
  // Flat form for field::energy / field::distance and the GPU (spec 15.1.4).
  FieldModel model() const;
  const GridGeometry& geom() const { return grids_->geom; }
  int mapOf(Ad4Type t) const { return grids_->mapOf(t); }
  double distanceCap() const { return dist_->cap; }
  std::size_t movedAtomCount() const { return moved_.size(); }  // flexres
  std::size_t blockCount() const { return block_.size(); }      // delta: all blocks
  std::size_t deltaBlockCount() const;                          // delta: blocks with their own values

  // Energy of a probe atom of map `map` and charge q at p (inside geom()); grad: dE/dp.
  template <class T>
  T energy(int map, T q, const Vec3T<T>& p, Vec3T<T>* grad) const;
  // Distance to the nearest StaM heavy atom (at most distanceCap()); grad: dd/dp.
  template <class T>
  T distance(const Vec3T<T>& p, Vec3T<T>* grad) const;

  int typeMapCount() const { return static_cast<int>(grids_->types.size()); }
  // Values of map m at the grid points (m: type maps, then elec, desolv), the moved atoms
  // of flexres included: E = typeMap + q * elec + |q| * desolv at every point (posm2_fft).
  std::vector<float> pointMap(int m) const;

 private:
  // flexres: pair sum over the moved StaM atoms; adds dE/dp to grad.
  template <class T>
  T movedSum(int m, T q, const Vec3T<T>& p, Vec3T<T>* grad) const;
  template <class T>
  T map(int m, const Vec3T<T>& p, Vec3T<T>* grad) const;  // m: type maps, then elec, desolv, distance

  GridMode mode_ = GridMode::Full;
  std::uint64_t id_ = 0;
  std::shared_ptr<const EnergyGrids> grids_;
  std::shared_ptr<const DistanceGrid> dist_;
  // flexres
  struct Moved {
    Vec3 xyz{};
    int type = 0;           // Ad4Type
    double elecQ = 0;       // c_elec q_j
    double S = 0, V = 0;    // S_j, V_j
    bool heavy = false;
  };
  std::vector<Moved> moved_;
  std::vector<int> movedHeavy_;          // indices into moved_
  detail::CellList cells_, heavyCells_;  // over moved_, over movedHeavy_
  std::vector<double> solT_, volT_;      // per map
  std::vector<double> pairAC_, pairBD_, pairReqm_;  // per (map, Ad4Type)
  std::vector<int> pairM_, pairN_;
  double cDesolv_ = 0, qasp_ = 0, halfSmooth_ = 0, rMin_ = 0;
  // delta
  std::array<int, 3> nb_{};
  std::vector<int> block_;   // per block: -1 (reference) or index of its values in data_
  std::vector<float> data_;  // per stored block: maps x kDeltaBlock^3 values
  int mapCount_ = 0;         // type maps + elec + desolv + distance
};

}  // namespace od
