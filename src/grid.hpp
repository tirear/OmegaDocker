// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "ad4.hpp"
#include "energy.hpp"
#include "molecule.hpp"

namespace od {

// StaM–MobM intermolecular energy grids (spec 8.1, gridMode = full) and the
// distance grid (spec 8.2). CPU reference implementation.

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

struct InterEnergy {
  double energy = 0;
  int outside = 0;  // probe atoms outside the grid (not counted in energy)
};
// Grid energy of the probe atoms; grad (optional, per probe atom, kcal/mol/Å).
InterEnergy gridEnergy(const EnergyGrids& g, const GridAtoms& probe, std::vector<Vec3>* grad = nullptr);

// Reference: direct sum over StaM–probe atom pairs with the grid's terms (no clamp).
double pairInterEnergy(const GridAtoms& stam, const GridAtoms& probe, const EnergyParams& params);

}  // namespace od
