// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// Command-line parameters (spec section 14). Options are added as the
// corresponding stages are implemented.
struct Options {
  std::string sta;
  std::vector<std::string> mob;
  std::string molConfig;
  std::vector<std::string> ccd;
  std::optional<Vec3> boxCenter;
  std::optional<Vec3> boxSize;
  std::string outPrefix;
  bool outPdb = true;
  bool outCif = false;
  bool outSdf = false;
  bool outMol2 = false;
  double omegaTolerance = 30.0;
  bool acceptNonPlanarOmega = false;
  double staFlexMargin = 4.0;
  bool staFixedResMovable = false;
  double maxBreakDist = 0.5;
  double closureG = 200.0;
  double clashDistance = 1.8;
  double Emax = 1000.0;
  double EclashBase = 10000.0;
  double Kclash = 1000.0;
  std::optional<int> dockingConfSampleNum;
  std::optional<int> minL2confSampleNum;
  std::optional<int> initL1confSampleNum;
  std::optional<int> eliteConfSampleNum;
  long long maxL1TotalNum = 1000000000;
  long long maxL1BatchNum = 65536;
  double csgm1AngleStep = 1.0;
  double ramaFloor = 1e-4;
  int confMinMaxStep = 100;
  double confMinDeltaE = 10.0;
  double confMaxDihedralStep = 10.0;
  std::optional<std::uint64_t> seed;  // random when not given

  bool showHelp = false;
  bool showVersion = false;
};

// Parses argv; throws od::Error on invalid input.
Options parseCommandLine(int argc, const char* const* argv);

// Sample numbers of spec 7.2 with their defaults applied.
struct ConfSampleNums {
  int docking = 0, minL2 = 0, initL1 = 0, elite = 0;
};
ConfSampleNums confSampleNums(const Options& o);

// Checks required parameters and value ranges; throws od::Error.
void validateOptions(const Options& o);

std::string usageText();

}  // namespace od
