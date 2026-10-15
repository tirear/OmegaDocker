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
// Sample numbers of spec 7.2 as given on the command line.
struct ConfSampleOpts {
  std::optional<int> docking, minL2, initL1, elite;
};

// 9a batch limits (spec 7.4) and csgm parameters (spec 7.5, 7.8), one set for
// the StaM group and one shared by all MobM.
struct CsgmOpts {
  long long maxL1TotalNum = 1000000000;
  long long maxL1BatchNum = 65536;
  bool grow = false;  // --<set>CsgmBuild grow
  int growMaxRetry = 16;
  double csgm1AngleStep = 1.0;
  double ramaFloor = 1e-4;
};

struct Options {
  std::string stamFile;                // --stam
  std::vector<std::string> mobmFiles;  // --mobm
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
  double stamFlexMargin = 4.0;
  bool stamFixedResMovable = false;
  double maxBreakDist = 0.5;
  double closureG = 200.0;
  double reachProtein = 3.8;  // --reachPerRes <protein>,<nucleic>
  double reachNucleic = 7.0;
  double clashDistance = 1.8;
  double Emax = 1000.0;
  double EclashBase = 10000.0;
  double Kclash = 1000.0;
  ConfSampleOpts stam;  // --stam*ConfSampleNum: the StaM group
  ConfSampleOpts mobm;  // --mobm*ConfSampleNum: shared by all MobM
  CsgmOpts stamCsgm;  // --stamCsgmBuild, --stamGrowMaxRetry, ...
  CsgmOpts mobmCsgm;  // --mobmCsgmBuild, --mobmGrowMaxRetry, ...
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
enum class SampleSet : std::uint8_t { StaM, MobM };
ConfSampleNums confSampleNums(const Options& o, SampleSet set);

// Checks that the set is given and valid (spec 7.2); used when the StaM group
// or any MobM has rotatable torsions. Throws od::Error.
void requireConfSampleNums(const Options& o, SampleSet set);

// Checks required parameters and value ranges; throws od::Error.
void validateOptions(const Options& o);

std::string usageText();

}  // namespace od
