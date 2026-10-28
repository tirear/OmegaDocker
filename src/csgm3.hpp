// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <span>
#include <vector>

#include "rng.hpp"

namespace od {

// csgm3, Path Exploring (spec 7.7): samples on lines through the dcss samples
// in torsion space. Every ordered pair (i, j), i != j, of the k end points gives
// one line with peSampleNumPerLine samples.

struct Csgm3Params {
  int maxOffset = 30;         // peMaxOffset (degrees, integer)
  double coverAngle = 360.0;  // peCoverAngle
  int perLine = 4;            // peSampleNumPerLine
};

// Samples of one pass over all lines: k (k - 1) perLine.
long long csgm3SampleNum(std::size_t k, const Csgm3Params& p);

// Sample m (0 .. perLine - 1) of line `line` (0 .. k (k - 1) - 1) written to
// anglesDeg. lineRng: the same generator state for every sample of the line.
// cisTrans: genes of cisTrans omegas, set to the nearer of 0 and 180 degrees. Angles after the
// torsions (cisTrans.size()) are ring slots (spec 4.9): the conformer of end point i for
// m < perLine / 2, otherwise of end point j.
void csgm3Sample(Rng lineRng, std::span<const std::vector<float>* const> ends, long long line, int m,
                 const std::vector<bool>& cisTrans, const Csgm3Params& p, std::span<double> anglesDeg);

}  // namespace od
