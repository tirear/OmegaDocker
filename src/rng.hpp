// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>

namespace od {

// Counter-based random streams: the generator of a sample depends only on
// (seed, stream, index), so results do not depend on threads or batching.
class Rng {
 public:
  Rng(std::uint64_t seed, std::uint64_t stream, std::uint64_t index)
      : s_(mix(mix(seed ^ 0x6a09e667f3bcc909ULL) ^ mix(stream + 0xbb67ae8584caa73bULL) ^ index)) {}

  std::uint64_t next() {
    s_ += 0x9e3779b97f4a7c15ULL;
    return mix(s_);
  }
  // [0, 1)
  double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
  // [lo, hi)
  double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }
  // [0, n)
  std::uint64_t below(std::uint64_t n) { return static_cast<std::uint64_t>(uniform() * static_cast<double>(n)); }
  // Internal state, e.g. to continue the stream on a GPU (spec 15.1.18).
  std::uint64_t state() const { return s_; }
  void setState(std::uint64_t s) { s_ = s; }

 private:
  static std::uint64_t mix(std::uint64_t z) {  // SplitMix64 finalizer
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  std::uint64_t s_;
};

// Angle wrapped to (-180, 180].
inline double wrapDeg(double a) {
  a = a - 360.0 * static_cast<double>(static_cast<long long>(a / 360.0));
  if (a > 180.0) a -= 360.0;
  if (a <= -180.0) a += 360.0;
  return a;
}

}  // namespace od
