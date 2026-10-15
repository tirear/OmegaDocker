// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <vector>

#include "molecule.hpp"
#include "rng.hpp"

namespace od {

// Orientations as unit quaternions (spec 8.3).
struct Quat {
  double w = 1, x = 0, y = 0, z = 0;
};

Quat operator*(const Quat& a, const Quat& b);
Quat normalized(const Quat& q);
Vec3 rotate(const Quat& q, const Vec3& v);
// Rotation by |v| radians about v / |v|.
Quat fromRotationVector(const Vec3& v);
// Angle (radians, in [0, pi]) of the rotation taking a to b.
double rotationAngle(const Quat& a, const Quat& b);

// Uniformly distributed random orientation (Shoemake 1992).
Quat randomQuat(Rng& rng);

// Hopf-fibration grid of SO(3) (Yershova et al. 2010): HEALPix on S2 with
// Nside = 2^level and 6 * 2^level points on S1; 72 * 8^level orientations.
std::vector<Quat> hopfGrid(int level);
double hopfSpacingDeg(int level);  // 60 / 2^level
// Level whose spacing is closest to stepDeg (levels 0..maxLevel).
int hopfLevelForStep(double stepDeg, int maxLevel = 6);

}  // namespace od
