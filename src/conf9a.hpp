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

// Index of a generated sample: global numbers the random streams (counter of
// fillL2Pool); local counts from 0 within the current L1 batch (all of its parts).
struct GenIndex {
  std::uint64_t global = 0;
  long long local = 0;
};

// Outcome of generating one sample; only Built samples are minimized.
enum class GenResult : std::uint8_t { Built, GrowFailed, CloseFailed, SideGrowFailed, PreNotClosed };

// Writes one sample to anglesDeg; GrowFailed / CloseFailed / SideGrowFailed:
// discarded before minimization (csgmBuild = grow, loop closure of spec 7.11,
// building after the closure with csgmBuild = growClose).
// rng: Rng(seed, stream, index.global).
using ConfGenerator = std::function<GenResult(Rng& rng, const GenIndex& index, std::span<double> anglesDeg)>;

// Batched 9a energy evaluation of an EnergySystem, e.g. on a GPU (spec 15.1.3).
class ConfEvaluator {
 public:
  virtual ~ConfEvaluator() = default;
  virtual int maxBatch() const = 0;  // most samples per call
  // anglesDeg: n × torsionCount; terms: n; gradRad (optional): n × torsionCount, as EnergySystem::evaluate.
  virtual void evaluate(const double* anglesDeg, int n, EnergyTerms* terms, double* gradRad) = 0;
  // Called when a fillL2Pool call ends (e.g. to free work memory).
  virtual void idle() {}
  // Asynchronous form (spec 15.1.11), as PoseEvaluator: with lanes() > 1, up to lanes()
  // evaluations with gradients are in flight; complete stores the results of the lane's last submit.
  virtual int lanes() const { return 1; }
  virtual void submit(int /*lane*/, const double* /*anglesDeg*/, int /*n*/) {}
  virtual void complete(int /*lane*/, EnergyTerms* /*terms*/, double* /*gradRad*/) {}
};

struct Fill9aParams {
  int minL2 = 0;                       // minL2confSampleNum
  long long maxL1Batch = 65536;        // --stamMaxL1BatchNum / --mobmMaxL1BatchNum
  long long maxL1Total = 1000000000;   // --stamMaxL1TotalNum / --mobmMaxL1TotalNum
  double Emax = 1000.0;
  double maxBreakDist = 0.5;
  // > 0: samples whose closure deviation exceeds preCloseFactor × maxBreakDist before the
  // minimization are discarded without minimizing them (conf9a.preCloseFactor).
  double preCloseFactor = 1.0;
  MinParams min;
  std::uint64_t seed = 0;
  std::uint64_t stream = 0;            // one per StaM group / MobM
  long long fixedL1 = 0;               // > 0: every L1 batch has this size (csgm3) instead of spec 7.4
  ConfEvaluator* evaluator = nullptr;  // minimization energies; null: EnergySystem::evaluate on the CPU
};

struct Fill9aStats {
  long long l1 = 0;         // samples generated
  long long added = 0;      // samples added to the pool
  long long overEmax = 0;   // rejected: filterEnergy >= Emax
  long long notClosed = 0;  // rejected: closure deviation > maxBreakDist (energy below Emax)
  long long growFailed = 0; // discarded by the generator (csgmBuild = grow)
  long long closeFailed = 0; // discarded by the generator (loop closure, spec 7.11)
  long long sideGrowFailed = 0; // discarded by the generator (csgmBuild = growClose, after the closure)
  long long preNotClosed = 0;   // discarded before the minimization (preCloseFactor)
  int batches = 0;
  long long minSteps = 0;   // minimization steps over all samples
  long long nextL1 = 0;     // L1 count of a following refill (spec 7.4)
};

// Generates batches of L1 samples, minimizes them and adds those passing the
// filters to pool until pool.size() >= minL2 (spec 7.3, 7.4). l1 is the size of
// the first batch (fixedL1 when set); counter numbers the random streams of the samples and is
// advanced by the samples generated. Throws od::Error when maxL1Total samples
// are generated without reaching minL2.
Fill9aStats fillL2Pool(const EnergySystem& sys, const ConfGenerator& gen, long long l1, const Fill9aParams& p,
                       std::uint64_t& counter, std::vector<ConfSample>& pool);

// Samples anglesDeg (count × torsionCount) minimized in place with ev; terms and steps: count each.
void minimizeConfBatch(ConfEvaluator& ev, int torsionCount, std::span<double> anglesDeg, const MinParams& p,
                       EnergyTerms* terms, int* steps);

// One sample: minimizes anglesDeg in place; returns the energy terms at the result.
EnergyTerms minimizeConf(const EnergySystem& sys, std::span<double> anglesDeg, const MinParams& p, int* steps = nullptr);

}  // namespace od
