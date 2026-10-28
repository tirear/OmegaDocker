// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <span>
#include <vector>

#include "rng.hpp"

namespace od {

// csgm2 (spec 7.6): genetic operators of AutoDock-GPU 1.6 (kernel4) on the
// torsion genes only. Parents are the dcss samples; fitness is the lowest
// docking energy of the round (lower is better).

struct Csgm2Params {
  double tournamentRate = 0.6;  // probability that the fitter of two candidates wins
  double crossoverRate = 0.8;
  double mutationRate = 0.02;   // per gene
  double absMaxDang = 90.0;     // mutation: uniform in [-absMaxDang, absMaxDang) degrees
};

struct Csgm2Parent {
  const std::vector<float>* angles = nullptr;  // degrees, one per torsion
  float fitness = 0;
};

// One offspring written to anglesDeg (size of the parents' angles). cisTrans:
// genes of cisTrans omegas, which mutate by switching between 0 and 180 degrees.
// ringK: conformer counts of the ring slots after the torsions (spec 4.9); a mutated ring gene
// takes a uniform conformer index.
void csgm2Offspring(Rng& rng, std::span<const Csgm2Parent> parents, const std::vector<bool>& cisTrans,
                    const Csgm2Params& p, std::span<double> anglesDeg, std::span<const int> ringK = {});

}  // namespace od
