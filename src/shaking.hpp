// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>
#include <span>
#include <vector>

#include "energy.hpp"

namespace od {

// Atom shaking entropy of od_0100 (spec 10.2): each heavy atom is moved over
// the points of a sphere around its position while all other atoms stay; the
// pair energies at the points give its accessible fraction
//   f_i = sum_k exp(-(E_k - E_min) * rrtScaling) / (number of points),
// and -TS = -RT * sum_i ln f_i (>= 0; the coefficient is applied in 9d).

// Non-merged atoms of one or more molecules. Pairs of atoms of the same
// molecule within three bonds are excluded; all other pairs interact.
struct ShakeAtoms {
  std::vector<Vec3f> xyz;
  std::vector<Ad4Type> type;
  std::vector<float> charge;
  std::vector<char> heavy;
  std::vector<std::vector<int>> near;  // sorted indices of the excluded partners
  std::size_t size() const { return xyz.size(); }
};
// coords: coordinates of sys.evaluate for the same parts (any frame).
ShakeAtoms shakeAtoms(const EnergySystem& sys, const std::vector<EnergyPart>& parts, std::span<const Vec3f> coords);
// Atoms of a followed by atoms of b; no exclusions between them.
ShakeAtoms joinShakeAtoms(const ShakeAtoms& a, const ShakeAtoms& b);

struct ShakeParams {
  double radius = 1.4;         // Å, --shakingRadius
  double step = 0.2;           // Å, --shakingStep
  double temperature = 298.15;  // K, --temperature (RT of -TS)
  double rrtScaling = 4.0;      // mol/kcal, --rrtScaling
};

// Points (offsets) of the sphere: a cubic lattice with spacing step centred on
// the atom, inside radius (the centre is point 0).
std::vector<Vec3f> shakePoints(double radius, double step);

// ln f_i of every atom (0 for atoms that are not shaken: hydrogens and atoms
// with shake[i] == 0 when shake is not empty).
std::vector<double> shakeLnF(const ShakeAtoms& atoms, const EnergyParams& ep, const ShakeParams& p,
                             std::span<const char> shake = {});

// shakeLnF(atoms, ep, p, shake) of a back end (e.g. GpuShaker); empty: shakeLnF on the CPU.
using ShakeFn = std::function<std::vector<double>(const ShakeAtoms& atoms, std::span<const char> shake)>;

constexpr double kGasConstant = 1.987204259e-3;  // kcal/(mol K)
inline double minusTS(double sumLnF, const ShakeParams& p) { return -kGasConstant * p.temperature * sumLnF; }

}  // namespace od
