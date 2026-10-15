// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#include "energy.hpp"
#include "minimize.hpp"
#include "rng.hpp"

namespace od {

// 9a: conformation samples of a StaM group or one MobM (spec 7.3, 7.4, 7.9).

struct ConfSample {
  std::vector<float> angles;  // degrees, one per tree torsion of the system
  float energy = 0;           // EnergyTerms::filterEnergy after 9a minimization
};

using ConfGenerator = std::function<void(Rng& rng, std::span<double> anglesDeg)>;

struct Fill9aParams {
  int minL2 = 0;                       // minL2confSampleNum
  long long maxL1Batch = 65536;        // maxL1BatchNum
  long long maxL1Total = 1000000000;   // maxL1TotalNum
  double Emax = 1000.0;
  double maxBreakDist = 0.5;
  MinParams min;
  std::uint64_t seed = 0;
  std::uint64_t stream = 0;            // one per StaM group / MobM
};

struct Fill9aStats {
  long long l1 = 0;         // samples generated
  long long added = 0;      // samples added to the pool
  long long overEmax = 0;   // rejected: filterEnergy >= Emax
  long long notClosed = 0;  // rejected: closure deviation > maxBreakDist (energy below Emax)
  int batches = 0;
  long long minSteps = 0;   // minimization steps over all samples
  long long nextL1 = 0;     // L1 count of a following refill (spec 7.4)
};

// Generates batches of L1 samples, minimizes them and adds those passing the
// filters to pool until pool.size() >= minL2 (spec 7.3, 7.4). l1 is the size of
// the first batch; counter numbers the random streams of the samples and is
// advanced by the samples generated. Throws od::Error when maxL1Total samples
// are generated without reaching minL2.
Fill9aStats fillL2Pool(const EnergySystem& sys, const ConfGenerator& gen, long long l1, const Fill9aParams& p,
                       std::uint64_t& counter, std::vector<ConfSample>& pool);

// One sample: minimizes anglesDeg in place; returns the energy terms at the result.
EnergyTerms minimizeConf(const EnergySystem& sys, std::span<double> anglesDeg, const MinParams& p, int* steps = nullptr);

}  // namespace od
