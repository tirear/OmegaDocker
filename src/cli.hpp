// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "feparams.hpp"
#include "molecule.hpp"

namespace od {

// Run parameters (spec section 14): built-in defaults, then --paramFile (TOML),
// then the command line, where --<table>_<key>=<value> overrides <table>.<key>.
// Sample numbers of spec 7.2 as given.
struct ConfSampleOpts {
  std::optional<int> docking, minL2, initL1, elite;
};

// 9a batch limits (spec 7.4) and csgm parameters (spec 7.5, 7.8), one set for
// the StaM group and one shared by all MobM.
struct CsgmOpts {
  long long maxL1TotalNum = 1000000000;
  long long maxL1BatchNum = 65536;
  std::string build = "whole";  // <set>.csgmBuild: whole | grow | growClose
  int growMaxRetry = 16;
  int growMaxBacktrack = 4;
  double csgm1AngleStep = 1.0;
  double ramaFloor = 1e-4;
  // csgm2 (spec 7.6); rates in percent as AutoDock-GPU.
  double tournamentRate = 60.0;
  double crossoverRate = 80.0;
  double mutationRate = 2.0;
  double absMaxDang = 90.0;
  bool pathExploring = false;  // <set>.csgmGA = "csgm3"
  // csgm3 (spec 7.7).
  int peMaxOffset = 30;
  double peCoverAngle = 360.0;
  int peSampleNumPerLine = 4;
  // Loop closure (spec 7.11).
  std::string closure = "rcd";  // none | rcd | kic
  double rcdTolerance = 0.25;
  int rcdMaxIter = 1000;
  int rcdMaxRestart = 10;
  int kicMaxTrial = 100;
};

struct Options {
  // Command line only.
  std::string paramFile;               // --paramFile
  std::string stamFile;                // --stamFile
  std::vector<std::string> mobmFiles;  // --mobmFiles
  std::string molConfig;
  std::string outPrefix;
  std::string resume;                  // --resume: output directory of the interrupted run (spec 16.4)
  // --paramFile and command line.
  std::vector<std::string> ccd;
  std::optional<Vec3> boxCenter;
  std::optional<Vec3> boxSize;
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
  double reachProtein = 3.8;  // topology.reachPerRes = [protein, nucleic]
  double reachNucleic = 7.0;
  double clashDistance = 1.8;
  double Emax = 1000.0;
  double EclashBase = 10000.0;
  double Kclash = 1000.0;
  ConfSampleOpts stam;  // stam.*ConfSampleNum: the StaM group
  ConfSampleOpts mobm;  // mobm.*ConfSampleNum: shared by all MobM
  double gridSpacing = 0.375;
  double maxNearDistance = 6.0;
  double rotStep = 7.5;
  std::string posm;  // posm1 | posm2
  double posStep = 1.0;
  long long posm1SampleNum = 1000;
  int maxOverlapAtomNum = 1;
  int dockMinMaxStep = 1000;
  int dockMinConvSteps = 10;
  double dockMinDeltaE = 1e-3;
  int lbfgsMemory = 7;
  double scaleTrans = 1, scaleRot = 1, scaleTors = 1;
  double maxTransStep = 0.5;
  double maxAngleStep = 10.0;
  int maxRound = 10;
  std::optional<int> topSampleNum;  // default: eliteStaM x eliteMobM (spec 9)
  double diversityRmsd = 2.0;
  int sampleBufferNum = 65536;
  std::optional<int> midTopSampleFactor;  // default: 2 with midKeepSkipped, otherwise 1 (spec 9)
  bool midKeepSkipped = false;
  int symMaxPerm = 256;  // raised to 4 x topSampleNum when smaller (spec 9)
  bool feAdgpu16 = false;  // dock9d.freeEnergy (spec 10.2)
  bool feOd0100 = true;
  int outMinPerRank = 40;
  FeParams fe;  // [od_0100]
  CsgmOpts stamCsgm;  // stam.csgmBuild, stam.growMaxRetry, ...
  CsgmOpts mobmCsgm;  // mobm.csgmBuild, mobm.growMaxRetry, ...
  int confMinMaxStep = 100;
  double confMinDeltaE = 10.0;
  double confMaxDihedralStep = 10.0;
  std::optional<std::uint64_t> seed;  // random when not given

  bool showHelp = false;
  bool showVersion = false;
};

// Parses argv and reads --paramFile; throws od::Error on invalid input (all
// value problems of the file and the command line are listed together).
Options parseCommandLine(int argc, const char* const* argv);

// The parameter file with the values of o (spec 14); unset optional values
// are written as comments.
std::string paramsText(const Options& o);

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

// dock9c.midTopSampleFactor, or its default.
inline int midTopSampleFactor(const Options& o) { return o.midTopSampleFactor.value_or(o.midKeepSkipped ? 2 : 1); }

std::string usageText();

}  // namespace od
