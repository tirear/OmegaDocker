// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "molecule.hpp"

namespace od {

class EnergySystem;
class Csgm1GroupSampler;
class LoopCloser;
struct TorsionTree;

// Docking box (spec 3.3): center and size in box coordinates, whose axes are the rows of axes in
// input coordinates (x_box = axes x_input). Without box.autoRotate the axes are the input axes.
struct BoxFrame {
  Vec3 center{}, size{};
  std::array<Vec3, 3> axes{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  bool rotated = false;
  Vec3 toBox(const Vec3& x) const;
  Vec3 toInput(const Vec3& x) const;
};

// Coordinates of the heavy atoms.
std::vector<Vec3> heavyCoords(const Molecule& m);
// Largest distance of a point from their unweighted mean (0 when empty).
double radiusOf(std::span<const Vec3> x);

// box.mobmRadius = max (spec 3.3): upper bound, over all conformations of the tree, of the distance of
// a heavy atom from the heavy-atom centre: max_i (1/N) sum_j d_ub(i, j), with d_ub the sum of the
// fixed distances along the path through the rigid parts of the tree. With ring operations (spec 4.9),
// a distance in a rigid part is the largest over the combinations of the ring conformers that change
// distances in it (all of them up to maxCombos per part; above, each ring operation alone, and
// *overflow, when given, counts those parts: the result is then not guaranteed to be a bound).
struct RingOps;
double mobmRadiusMax(const Molecule& m, const TorsionTree& tree, const RingOps* rings = nullptr,
                     int* overflow = nullptr, long long maxCombos = 100000);

// box.mobmRadius = sampled: the largest radius of the input conformation and of the csgm1 samples
// (sample k with Rng(seed, stream, k), k < sampleNum) without clashing heavy-atom pairs. The ring
// slots of sys (spec 4.9) are uniform, from Rng(seed, 2^61 + stream, k) as in 9a. With an
// active closer (loop, crosslink, metal or macrocycle breaks), each sample is closed first (spec 7.11);
// samples that fail to close, or whose closure deviation then exceeds maxDeviation (> 0), are left out.
struct SampledRadius {
  double radius = 0;
  int accepted = 0;     // samples kept
  int notClosed = 0;    // samples left out by the closure
};
SampledRadius mobmRadiusSampled(const Molecule& m, const EnergySystem& sys, const Csgm1GroupSampler& sampler,
                                int sampleNum, std::uint64_t seed, std::uint64_t stream,
                                const LoopCloser* closer = nullptr, double maxDeviation = 0);

// Principal axes of the points (rows, by decreasing variance; right-handed, the largest component
// of each of the first two axes positive). The input axes for fewer than two points.
std::array<Vec3, 3> principalAxes(std::span<const Vec3> x);

struct AutoBoxParams {
  double stamMargin = 4.0;  // Å, per axis in total
  double mobmRadius = 0.0;  // r_max (Å)
  bool autoRotate = true;
  std::optional<Vec3> center;  // input coordinates
  std::optional<Vec3> size;    // along the box axes
};

// Box of box.mode = sidechainAuto / fixAuto (spec 3.3) around the StaM heavy atoms. Axes
// (autoRotate): the principal axes, then the rotation nearby with the smallest volume of the
// extents plus stamMargin + 2 r_max per axis. Center: middle of the extents, or the given center;
// size: the extents (with a given center: twice the largest extent from it) plus stamMargin + 2 r_max,
// or the given size.
BoxFrame autoBox(std::span<const Vec3> stamHeavy, const AutoBoxParams& p);

}  // namespace od
