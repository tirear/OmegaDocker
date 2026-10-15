// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "orient.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace od {

namespace {

constexpr double kPi = std::numbers::pi;

// Centre (theta, phi) of HEALPix pixel p in the RING scheme.
void healpixCenter(long long nside, long long p, double& theta, double& phi) {
  const long long npix = 12 * nside * nside;
  const long long ncap = 2 * nside * (nside - 1);
  const double n2 = static_cast<double>(nside * nside);
  double z;
  if (p < ncap) {
    const long long iring = static_cast<long long>((1 + std::sqrt(1.0 + 2.0 * static_cast<double>(p))) / 2);
    const long long iphi = p + 1 - 2 * iring * (iring - 1);
    z = 1 - static_cast<double>(iring * iring) / (3 * n2);
    phi = (static_cast<double>(iphi) - 0.5) * kPi / (2.0 * static_cast<double>(iring));
  } else if (p < npix - ncap) {
    const long long ip = p - ncap;
    const long long iring = ip / (4 * nside) + nside;
    const long long iphi = ip % (4 * nside) + 1;
    const double fodd = ((iring + nside) & 1) ? 1.0 : 0.5;
    z = static_cast<double>(2 * nside - iring) * 2.0 / (3.0 * static_cast<double>(nside));
    phi = (static_cast<double>(iphi) - fodd) * kPi / (2.0 * static_cast<double>(nside));
  } else {
    const long long ip = npix - p;
    const long long iring = static_cast<long long>((1 + std::sqrt(2.0 * static_cast<double>(ip) - 1)) / 2);
    const long long iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
    z = -1 + static_cast<double>(iring * iring) / (3 * n2);
    phi = (static_cast<double>(iphi) - 0.5) * kPi / (2.0 * static_cast<double>(iring));
  }
  theta = std::acos(std::clamp(z, -1.0, 1.0));
}

}  // namespace

Quat operator*(const Quat& a, const Quat& b) {
  return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z, a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
          a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x, a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

Quat normalized(const Quat& q) {
  const double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
  return {q.w / n, q.x / n, q.y / n, q.z / n};
}

Vec3 rotate(const Quat& q, const Vec3& v) {
  // v + 2 u x (u x v + w v), u = (x, y, z)
  const double tx = q.y * v[2] - q.z * v[1] + q.w * v[0];
  const double ty = q.z * v[0] - q.x * v[2] + q.w * v[1];
  const double tz = q.x * v[1] - q.y * v[0] + q.w * v[2];
  return {v[0] + 2 * (q.y * tz - q.z * ty), v[1] + 2 * (q.z * tx - q.x * tz), v[2] + 2 * (q.x * ty - q.y * tx)};
}

Quat fromRotationVector(const Vec3& v) {
  const double a = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (a < 1e-12) return normalized({1, 0.5 * v[0], 0.5 * v[1], 0.5 * v[2]});
  const double s = std::sin(0.5 * a) / a;
  return {std::cos(0.5 * a), s * v[0], s * v[1], s * v[2]};
}

double rotationAngle(const Quat& a, const Quat& b) {
  const double d = std::abs(a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z);
  return 2 * std::acos(std::min(1.0, d));
}

Quat randomQuat(Rng& rng) {
  const double u1 = rng.uniform(), u2 = rng.uniform(), u3 = rng.uniform();
  const double a = std::sqrt(1 - u1), b = std::sqrt(u1);
  return {b * std::cos(2 * kPi * u3), a * std::sin(2 * kPi * u2), a * std::cos(2 * kPi * u2), b * std::sin(2 * kPi * u3)};
}

std::vector<Quat> hopfGrid(int level) {
  const long long nside = 1LL << level;
  const long long npix = 12 * nside * nside;
  const long long npsi = 6 * nside;
  std::vector<Quat> q;
  q.reserve(static_cast<std::size_t>(npix * npsi));
  for (long long p = 0; p < npix; ++p) {
    double theta, phi;
    healpixCenter(nside, p, theta, phi);
    for (long long k = 0; k < npsi; ++k) {
      const double psi = (static_cast<double>(k) + 0.5) * 2 * kPi / static_cast<double>(npsi);
      const double c = std::cos(0.5 * theta), s = std::sin(0.5 * theta);
      q.push_back({c * std::cos(0.5 * psi), c * std::sin(0.5 * psi), s * std::cos(phi + 0.5 * psi),
                   s * std::sin(phi + 0.5 * psi)});
    }
  }
  return q;
}

double hopfSpacingDeg(int level) { return 60.0 / static_cast<double>(1 << level); }

int hopfLevelForStep(double stepDeg, int maxLevel) {
  int best = 0;
  for (int l = 1; l <= maxLevel; ++l)
    if (std::abs(hopfSpacingDeg(l) - stepDeg) < std::abs(hopfSpacingDeg(best) - stepDeg)) best = l;
  return best;
}

}  // namespace od
