// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "cli.hpp"

#include <algorithm>
#include <charconv>
#include <format>
#include <functional>
#include <map>

#include "log.hpp"
#include "version.hpp"

namespace od {
namespace {

enum class Kind { Flag, Value };

struct OptionDef {
  std::string name;
  Kind kind;
  std::string valueHint;
  std::string defaultText;
  std::string help;
  std::function<void(Options&, const std::string&)> apply;
};

double parseDouble(const std::string& opt, const std::string& s) {
  double v = 0;
  auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc() || p != s.data() + s.size()) fail("--{}: '{}' is not a number", opt, s);
  return v;
}

Vec3 parseVec3(const std::string& opt, const std::string& s) {
  Vec3 v{};
  std::size_t start = 0;
  for (int k = 0; k < 3; ++k) {
    std::size_t comma = s.find(',', start);
    if ((k < 2) != (comma != std::string::npos)) fail("--{}: expected x,y,z but got '{}'", opt, s);
    v[k] = parseDouble(opt, s.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
    start = comma + 1;
  }
  return v;
}

long long parseInt(const std::string& opt, const std::string& s) {
  long long v = 0;
  auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc() || p != s.data() + s.size()) fail("--{}: '{}' is not an integer", opt, s);
  return v;
}

int parseInt32(const std::string& opt, const std::string& s) {
  const long long v = parseInt(opt, s);
  if (v < -2147483647LL || v > 2147483647LL) fail("--{}: '{}' is out of range", opt, s);
  return static_cast<int>(v);
}

bool parseBool(const std::string& opt, const std::string& s) {
  if (s == "true") return true;
  if (s == "false") return false;
  fail("--{}: expected true or false but got '{}'", opt, s);
}

// --<set><Name>ConfSampleNum, e.g. --stamDockingConfSampleNum, --mobmMinL2confSampleNum.
OptionDef sampleOption(const std::string& set, const std::string& name, std::optional<int> ConfSampleOpts::* field,
                       const std::string& defaultText, const std::string& help) {
  const std::string opt = set + name + (name == "MinL2" || name == "InitL1" ? "confSampleNum" : "ConfSampleNum");
  return {opt, Kind::Value, "<n>", defaultText, help, [opt, set, field](Options& o, const std::string& v) {
            (set == "stam" ? o.stam : o.mobm).*field = parseInt32(opt, v);
          }};
}

// --stam<Name> and --mobm<Name> for one CsgmOpts field.
template <typename T>
void csgmOption(std::vector<OptionDef>& table, const std::string& name, const std::string& metavar,
                const std::string& defaultText, const std::string& help,
                std::function<T(const std::string& opt, const std::string& v)> parse, T CsgmOpts::* field) {
  for (const std::string set : {"stam", "mobm"}) {
    const std::string opt = set + name;
    table.push_back({opt, Kind::Value, metavar, defaultText,
                     (set == "stam" ? "StaM group, " : "Every MobM, ") + help,
                     [opt, set, parse, field](Options& o, const std::string& v) {
                       (set == "stam" ? o.stamCsgm : o.mobmCsgm).*field = parse(opt, v);
                     }});
  }
}

std::vector<OptionDef> csgmOptions() {
  std::vector<OptionDef> t;
  csgmOption<long long>(t, "MaxL1TotalNum", "<n>", "1000000000", "9a: largest number of samples generated in one refill",
                        parseInt, &CsgmOpts::maxL1TotalNum);
  csgmOption<long long>(t, "MaxL1BatchNum", "<n>", "65536", "9a: largest number of samples processed at once", parseInt,
                        &CsgmOpts::maxL1BatchNum);
  csgmOption<bool>(t, "CsgmBuild", "whole|grow", "whole", "9a: apply all torsions at once, or build segment by segment",
                   [](const std::string& opt, const std::string& v) {
                     if (v != "whole" && v != "grow") fail("--{}: expected whole or grow but got '{}'", opt, v);
                     return v == "grow";
                   },
                   &CsgmOpts::grow);
  csgmOption<int>(t, "GrowMaxRetry", "<n>", "16", "CsgmBuild = grow: resamplings of a failing segment before discarding",
                  parseInt32, &CsgmOpts::growMaxRetry);
  csgmOption<double>(t, "Csgm1AngleStep", "<deg>", "1", "csgm1a / R: sampled angles are multiples of this step",
                     parseDouble, &CsgmOpts::csgm1AngleStep);
  csgmOption<double>(t, "RamaFloor", "<w>", "1e-4", "csgm1b: weight added to every Ramachandran cell", parseDouble,
                     &CsgmOpts::ramaFloor);
  return t;
}

// Inserts the 9a batch-limit and csgm options after the sample-number options.
std::vector<OptionDef> withCsgmOptions(std::vector<OptionDef> t) {
  const auto at =
      std::find_if(t.begin(), t.end(), [](const OptionDef& d) { return d.name == "mobmEliteConfSampleNum"; });
  const auto c = csgmOptions();
  t.insert(at + 1, c.begin(), c.end());
  return t;
}

const std::vector<OptionDef>& optionTable() {
  static const std::vector<OptionDef> table = withCsgmOptions({
      {"stamFile", Kind::Value, "<file>", "required", "Stationary molecules (StaM); may be given once",
       [](Options& o, const std::string& v) {
         if (!o.stamFile.empty()) fail("--stamFile may be given only once");
         o.stamFile = v;
       }},
      {"mobmFiles", Kind::Value, "<file>", "required", "Mobile molecule (MobM), one molecule per file; repeatable",
       [](Options& o, const std::string& v) { o.mobmFiles.push_back(v); }},
      {"molConfig", Kind::Value, "<file>", "none", "Molecule configuration file (TOML)",
       [](Options& o, const std::string& v) { o.molConfig = v; }},
      {"ccd", Kind::Value, "<file>", "none", "Additional CCD definitions (mmCIF); repeatable",
       [](Options& o, const std::string& v) { o.ccd.push_back(v); }},
      {"boxCenter", Kind::Value, "x,y,z", "required", "Docking box center (Å)",
       [](Options& o, const std::string& v) { o.boxCenter = parseVec3("boxCenter", v); }},
      {"boxSize", Kind::Value, "x,y,z", "required", "Docking box size (Å)",
       [](Options& o, const std::string& v) { o.boxSize = parseVec3("boxSize", v); }},
      {"outPrefix", Kind::Value, "<path>", "required", "Output prefix",
       [](Options& o, const std::string& v) { o.outPrefix = v; }},
      {"outPdb", Kind::Value, "true|false", "true", "Write PDB output",
       [](Options& o, const std::string& v) { o.outPdb = parseBool("outPdb", v); }},
      {"outCif", Kind::Value, "true|false", "false", "Write mmCIF output",
       [](Options& o, const std::string& v) { o.outCif = parseBool("outCif", v); }},
      {"outSdf", Kind::Value, "true|false", "false", "Write SDF output",
       [](Options& o, const std::string& v) { o.outSdf = parseBool("outSdf", v); }},
      {"outMol2", Kind::Value, "true|false", "false", "Write mol2 output",
       [](Options& o, const std::string& v) { o.outMol2 = parseBool("outMol2", v); }},
      {"omegaTolerance", Kind::Value, "<deg>", "30", "Allowed deviation of a fixed omega from 0 or 180 degrees",
       [](Options& o, const std::string& v) { o.omegaTolerance = parseDouble("omegaTolerance", v); }},
      {"acceptNonPlanarOmega", Kind::Flag, "", "", "Continue when a fixed omega exceeds --omegaTolerance",
       [](Options& o, const std::string&) { o.acceptNonPlanarOmega = true; }},
      {"stamFlexMargin", Kind::Value, "<dist>", "4", "StaM torsions are rotatable by default within the box extended by this margin (Å)",
       [](Options& o, const std::string& v) { o.stamFlexMargin = parseDouble("stamFlexMargin", v); }},
      {"stamFixedResMovable", Kind::Value, "true|false", "false",
       "StaM: false keeps the backbone of fixed residues in place; true lets the whole chain move by torsions",
       [](Options& o, const std::string& v) { o.stamFixedResMovable = parseBool("stamFixedResMovable", v); }},
      {"maxBreakDist", Kind::Value, "<dist>", "0.5", "Largest closure deviation of a broken bond kept after 9a (Å)",
       [](Options& o, const std::string& v) { o.maxBreakDist = parseDouble("maxBreakDist", v); }},
      {"closureG", Kind::Value, "<k>", "200", "Closure constraint stiffness G (kcal/mol/Å²)",
       [](Options& o, const std::string& v) { o.closureG = parseDouble("closureG", v); }},
      {"reachPerRes", Kind::Value, "<p>,<n>", "3.8,7.0",
       "CsgmBuild = grow: reach per residue of protein / nucleic-acid chains toward a closure (Å)",
       [](Options& o, const std::string& v) {
         const auto comma = v.find(',');
         if (comma == std::string::npos || v.find(',', comma + 1) != std::string::npos)
           fail("--reachPerRes: expected <protein>,<nucleic> but got '{}'", v);
         o.reachProtein = parseDouble("reachPerRes", v.substr(0, comma));
         o.reachNucleic = parseDouble("reachPerRes", v.substr(comma + 1));
       }},
      {"gridSpacing", Kind::Value, "<dist>", "0.375", "9b: spacing of the energy and distance grids (Å)",
       [](Options& o, const std::string& v) { o.gridSpacing = parseDouble("gridSpacing", v); }},
      {"posm", Kind::Value, "posm1|posm2", "required", "9b: position and orientation sampling method",
       [](Options& o, const std::string& v) {
         if (v == "posm2_fft" || v == "posm3" || v == "posm4") fail("--posm {} is not implemented yet", v);
         if (v != "posm1" && v != "posm2") fail("--posm: expected posm1 or posm2 but got '{}'", v);
         o.posm = v;
       }},
      {"posStep", Kind::Value, "<dist>", "1", "9b: spacing of the sampled positions of the MobM heavy-atom centre (Å)",
       [](Options& o, const std::string& v) { o.posStep = parseDouble("posStep", v); }},
      {"posm1SampleNum", Kind::Value, "<n>", "1000", "9b, posm1: L1poss per conformation pair",
       [](Options& o, const std::string& v) { o.posm1SampleNum = parseInt("posm1SampleNum", v); }},
      {"maxOverlapAtomNum", Kind::Value, "<n>", "1", "9b: poses with this many MobM heavy atoms closer than --clashDistance to the StaM are removed",
       [](Options& o, const std::string& v) { o.maxOverlapAtomNum = parseInt32("maxOverlapAtomNum", v); }},
      {"dockMinMaxStep", Kind::Value, "<n>", "1000", "9b minimization: largest number of steps",
       [](Options& o, const std::string& v) { o.dockMinMaxStep = parseInt32("dockMinMaxStep", v); }},
      {"dockMinConvSteps", Kind::Value, "<n>", "10", "9b minimization: stop after this many consecutive steps with |dE| < --dockMinDeltaE",
       [](Options& o, const std::string& v) { o.dockMinConvSteps = parseInt32("dockMinConvSteps", v); }},
      {"dockMinDeltaE", Kind::Value, "<E>", "1e-3", "9b minimization: energy change regarded as converged (kcal/mol)",
       [](Options& o, const std::string& v) { o.dockMinDeltaE = parseDouble("dockMinDeltaE", v); }},
      {"lbfgsMemory", Kind::Value, "<n>", "7", "L-BFGS: number of stored steps",
       [](Options& o, const std::string& v) { o.lbfgsMemory = parseInt32("lbfgsMemory", v); }},
      {"scaleTrans", Kind::Value, "<s>", "1", "9b minimization: scale of the translation variables (Å)",
       [](Options& o, const std::string& v) { o.scaleTrans = parseDouble("scaleTrans", v); }},
      {"scaleRot", Kind::Value, "<s>", "1", "9b minimization: scale of the rotation variables (rad)",
       [](Options& o, const std::string& v) { o.scaleRot = parseDouble("scaleRot", v); }},
      {"scaleTors", Kind::Value, "<s>", "1", "9b minimization: scale of the torsion variables (rad)",
       [](Options& o, const std::string& v) { o.scaleTors = parseDouble("scaleTors", v); }},
      {"maxTransStep", Kind::Value, "<dist>", "0.5", "9b minimization: largest translation per step (Å)",
       [](Options& o, const std::string& v) { o.maxTransStep = parseDouble("maxTransStep", v); }},
      {"maxAngleStep", Kind::Value, "<deg>", "10", "9b minimization: largest rotation and torsion change per step",
       [](Options& o, const std::string& v) { o.maxAngleStep = parseDouble("maxAngleStep", v); }},
      {"maxRound", Kind::Value, "<n>", "10", "Number of 9a-9c rounds",
       [](Options& o, const std::string& v) { o.maxRound = parseInt32("maxRound", v); }},
      {"topSampleNum", Kind::Value, "<n>", "eliteStaM x eliteMobM", "9c: size of the topSampleSet of each MobM",
       [](Options& o, const std::string& v) { o.topSampleNum = parseInt32("topSampleNum", v); }},
      {"diversityRmsd", Kind::Value, "<dist>", "2.0", "9c: samples with the same StaM conformation closer than this RMSD to a kept sample are skipped (Å; 0 disables)",
       [](Options& o, const std::string& v) { o.diversityRmsd = parseDouble("diversityRmsd", v); }},
      {"sampleBufferNum", Kind::Value, "<n>", "65536", "9b/9c: run 9c whenever this many docking samples are waiting (at least (midTopSampleFactor + 2) x topSampleNum)",
       [](Options& o, const std::string& v) { o.sampleBufferNum = parseInt32("sampleBufferNum", v); }},
      {"midTopSampleFactor", Kind::Value, "<n>", "2", "9c: the 9c runs before the last one of a round keep this many x topSampleNum samples (1: keep topSampleNum)",
       [](Options& o, const std::string& v) { o.midTopSampleFactor = parseInt32("midTopSampleFactor", v); }},
      {"symMaxPerm", Kind::Value, "<n>", "256", "9d: symmetry-corrected RMSD over the graph automorphisms when there are at most this many; otherwise one-to-one pairing within each AD4 type",
       [](Options& o, const std::string& v) { o.symMaxPerm = parseInt32("symMaxPerm", v); }},
      {"rotStep", Kind::Value, "<deg>", "7.5", "9b: spacing of the Hopf orientation grid (the closest level is used)",
       [](Options& o, const std::string& v) { o.rotStep = parseDouble("rotStep", v); }},
      {"maxNearDistance", Kind::Value, "<dist>", "6", "9b: poses whose heavy atoms are all farther from the StaM are removed (Å)",
       [](Options& o, const std::string& v) { o.maxNearDistance = parseDouble("maxNearDistance", v); }},
      {"clashDistance", Kind::Value, "<dist>", "1.8", "Heavy-atom pairs closer than this clash (Å)",
       [](Options& o, const std::string& v) { o.clashDistance = parseDouble("clashDistance", v); }},
      {"Emax", Kind::Value, "<E>", "1000", "9a: conformations with energy below this enter the L2 pool (kcal/mol)",
       [](Options& o, const std::string& v) { o.Emax = parseDouble("Emax", v); }},
      {"EclashBase", Kind::Value, "<E>", "10000", "Constant part of the clash penalty; must exceed --Emax (kcal/mol)",
       [](Options& o, const std::string& v) { o.EclashBase = parseDouble("EclashBase", v); }},
      {"Kclash", Kind::Value, "<k>", "1000", "Quadratic part of the clash penalty (kcal/mol/Å²)",
       [](Options& o, const std::string& v) { o.Kclash = parseDouble("Kclash", v); }},
      sampleOption("stam", "Docking", &ConfSampleOpts::docking, "required if rotatable",
                   "StaM group: conformation samples docked per round (dcss size, >= 4)"),
      sampleOption("stam", "MinL2", &ConfSampleOpts::minL2, "2 x docking",
                   "StaM group, 9a: refill the L2 pool to at least this many samples"),
      sampleOption("stam", "InitL1", &ConfSampleOpts::initL1, "minL2", "StaM group, 9a: samples generated in the first batch"),
      sampleOption("stam", "Elite", &ConfSampleOpts::elite, "max(2, ceil(docking/4))",
                   "StaM group: csgm1 is used while the dcss has fewer samples"),
      sampleOption("mobm", "Docking", &ConfSampleOpts::docking, "required if rotatable",
                   "Every MobM: conformation samples docked per round (dcss size, >= 4)"),
      sampleOption("mobm", "MinL2", &ConfSampleOpts::minL2, "2 x docking",
                   "Every MobM, 9a: refill the L2 pool to at least this many samples"),
      sampleOption("mobm", "InitL1", &ConfSampleOpts::initL1, "minL2", "Every MobM, 9a: samples generated in the first batch"),
      sampleOption("mobm", "Elite", &ConfSampleOpts::elite, "max(2, ceil(docking/4))",
                   "Every MobM: csgm1 is used while the dcss has fewer samples"),
      {"confMinMaxStep", Kind::Value, "<n>", "100", "9a minimization: largest number of steps",
       [](Options& o, const std::string& v) { o.confMinMaxStep = parseInt32("confMinMaxStep", v); }},
      {"confMinDeltaE", Kind::Value, "<E>", "10", "9a minimization: stop when one step changes the energy by less (kcal/mol)",
       [](Options& o, const std::string& v) { o.confMinDeltaE = parseDouble("confMinDeltaE", v); }},
      {"confMaxDihedralStep", Kind::Value, "<deg>", "10", "9a minimization: largest change of one torsion per step",
       [](Options& o, const std::string& v) { o.confMaxDihedralStep = parseDouble("confMaxDihedralStep", v); }},
      {"seed", Kind::Value, "<n>", "random", "Random seed (written to the log)",
       [](Options& o, const std::string& v) {
         const long long s = parseInt("seed", v);
         if (s < 0) fail("--seed must be >= 0");
         o.seed = static_cast<std::uint64_t>(s);
       }},
      {"help", Kind::Flag, "", "", "Show this help", [](Options& o, const std::string&) { o.showHelp = true; }},
      {"version", Kind::Flag, "", "", "Show version", [](Options& o, const std::string&) { o.showVersion = true; }},
  });;
  return table;
}

}  // namespace

Options parseCommandLine(int argc, const char* const* argv) {
  std::map<std::string, const OptionDef*> byName;
  for (const auto& d : optionTable()) byName[d.name] = &d;

  Options o;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) fail("Unexpected argument '{}'", arg);
    std::string name = arg.substr(2), value;
    bool hasInlineValue = false;
    if (auto eq = name.find('='); eq != std::string::npos) {
      value = name.substr(eq + 1);
      name = name.substr(0, eq);
      hasInlineValue = true;
    }
    auto it = byName.find(name);
    if (it == byName.end()) fail("Unknown option --{}", name);
    const OptionDef& d = *it->second;
    if (d.kind == Kind::Flag) {
      if (hasInlineValue) fail("--{} does not take a value", name);
    } else if (!hasInlineValue) {
      if (i + 1 >= argc) fail("--{} requires a value {}", name, d.valueHint);
      value = argv[++i];
    }
    d.apply(o, value);
  }
  return o;
}

ConfSampleNums confSampleNums(const Options& o, SampleSet set) {
  const ConfSampleOpts& s = set == SampleSet::StaM ? o.stam : o.mobm;
  ConfSampleNums n;
  n.docking = s.docking.value_or(0);
  n.minL2 = s.minL2.value_or(2 * n.docking);
  n.initL1 = s.initL1.value_or(n.minL2);
  n.elite = s.elite.value_or(std::max(2, (n.docking + 3) / 4));
  return n;
}

namespace {

void checkConfSampleNums(const Options& o, SampleSet set) {
  const std::string p = set == SampleSet::StaM ? "stam" : "mobm";
  const ConfSampleNums n = confSampleNums(o, set);
  if (n.docking < 4) fail("--{}DockingConfSampleNum must be >= 4", p);
  if (n.minL2 < 2 * n.docking)
    fail("--{0}MinL2confSampleNum must be >= 2 x --{0}DockingConfSampleNum ({1})", p, 2 * n.docking);
  if (n.initL1 < n.minL2) fail("--{0}InitL1confSampleNum must be >= --{0}MinL2confSampleNum ({1})", p, n.minL2);
  if (n.elite < 2 || n.elite > n.docking / 2)
    fail("--{0}EliteConfSampleNum must be within [2, --{0}DockingConfSampleNum / 2] = [2, {1}]", p, n.docking / 2);
}

}  // namespace

void requireConfSampleNums(const Options& o, SampleSet set) {
  const bool stam = set == SampleSet::StaM;
  if (!(stam ? o.stam : o.mobm).docking)
    fail("{}", stam ? "--stamDockingConfSampleNum is required: the StaM group has rotatable torsions"
                    : "--mobmDockingConfSampleNum is required: a MobM has rotatable torsions");
  checkConfSampleNums(o, set);
}

void validateOptions(const Options& o) {
  if (o.stamFile.empty()) fail("--stamFile is required");
  if (o.mobmFiles.empty()) fail("At least one --mobmFiles is required");
  if (!o.boxCenter) fail("--boxCenter is required");
  if (!o.boxSize) fail("--boxSize is required");
  for (double v : *o.boxSize)
    if (!(v > 0)) fail("--boxSize: all components must be > 0");
  if (o.outPrefix.empty()) fail("--outPrefix is required");
  if (!(o.omegaTolerance >= 0 && o.omegaTolerance <= 90)) fail("--omegaTolerance must be within [0, 90]");
  if (!(o.stamFlexMargin >= 0)) fail("--stamFlexMargin must be >= 0");
  if (!(o.maxBreakDist > 0)) fail("--maxBreakDist must be > 0");
  if (!(o.closureG > 0)) fail("--closureG must be > 0");
  if (!(o.reachProtein > 0 && o.reachNucleic > 0)) fail("--reachPerRes values must be > 0");
  if (!(o.clashDistance > 0)) fail("--clashDistance must be > 0");
  if (!(o.gridSpacing > 0)) fail("--gridSpacing must be > 0");
  if (!(o.maxNearDistance > 0)) fail("--maxNearDistance must be > 0");
  if (!(o.rotStep > 0)) fail("--rotStep must be > 0");
  if (o.posm.empty()) fail("--posm is required");
  if (!(o.posStep > 0)) fail("--posStep must be > 0");
  if (o.posm1SampleNum < 1) fail("--posm1SampleNum must be >= 1");
  if (o.maxOverlapAtomNum < 1) fail("--maxOverlapAtomNum must be >= 1");
  if (o.dockMinMaxStep < 0) fail("--dockMinMaxStep must be >= 0");
  if (o.dockMinConvSteps < 1) fail("--dockMinConvSteps must be >= 1");
  if (!(o.dockMinDeltaE >= 0)) fail("--dockMinDeltaE must be >= 0");
  if (o.lbfgsMemory < 1) fail("--lbfgsMemory must be >= 1");
  if (!(o.scaleTrans > 0 && o.scaleRot > 0 && o.scaleTors > 0)) fail("--scaleTrans, --scaleRot, --scaleTors must be > 0");
  if (!(o.maxTransStep > 0)) fail("--maxTransStep must be > 0");
  if (!(o.maxAngleStep > 0)) fail("--maxAngleStep must be > 0");
  if (o.maxRound < 1) fail("--maxRound must be >= 1");
  if (o.topSampleNum && *o.topSampleNum < 1) fail("--topSampleNum must be >= 1");
  if (!(o.diversityRmsd >= 0)) fail("--diversityRmsd must be >= 0");
  if (o.sampleBufferNum < 1) fail("--sampleBufferNum must be >= 1");
  if (o.midTopSampleFactor < 1) fail("--midTopSampleFactor must be >= 1");
  if (o.symMaxPerm < 0) fail("--symMaxPerm must be >= 0");
  if (!(o.Kclash >= 0)) fail("--Kclash must be >= 0");
  if (!(o.EclashBase > o.Emax)) fail("--EclashBase must be greater than --Emax");
  for (SampleSet set : {SampleSet::StaM, SampleSet::MobM}) {
    const ConfSampleOpts& s = set == SampleSet::StaM ? o.stam : o.mobm;
    if (s.docking) checkConfSampleNums(o, set);
    else if (s.minL2 || s.initL1 || s.elite)
      fail("--{}DockingConfSampleNum is required when other --{}*ConfSampleNum options are given",
           set == SampleSet::StaM ? "stam" : "mobm", set == SampleSet::StaM ? "stam" : "mobm");
  }
  for (const auto& [set, c] : {std::pair{"stam", &o.stamCsgm}, std::pair{"mobm", &o.mobmCsgm}}) {
    if (c->maxL1TotalNum < 1) fail("--{}MaxL1TotalNum must be >= 1", set);
    if (c->maxL1BatchNum < 1) fail("--{}MaxL1BatchNum must be >= 1", set);
    if (c->growMaxRetry < 0) fail("--{}GrowMaxRetry must be >= 0", set);
    if (!(c->csgm1AngleStep > 0 && c->csgm1AngleStep <= 180)) fail("--{}Csgm1AngleStep must be within (0, 180]", set);
    if (!(c->ramaFloor >= 0)) fail("--{}RamaFloor must be >= 0", set);
  }
  if (o.confMinMaxStep < 0) fail("--confMinMaxStep must be >= 0");
  if (!(o.confMinDeltaE >= 0)) fail("--confMinDeltaE must be >= 0");
  if (!(o.confMaxDihedralStep > 0)) fail("--confMaxDihedralStep must be > 0");
  if (!o.outPdb && !o.outCif && !o.outSdf && !o.outMol2)
    fail("All output formats are disabled (--outPdb, --outCif, --outSdf, --outMol2)");
}

std::string usageText() {
  std::string s = std::format("OmegaDocker {}\n\nUsage: omegadocker --stamFile <file> --mobmFiles <file> [--mobmFiles <file> ...] "
                              "--boxCenter x,y,z --boxSize x,y,z --outPrefix <path> [options]\n\nOptions:\n",
                              kVersion);
  for (const auto& d : optionTable()) {
    std::string left = "--" + d.name + (d.valueHint.empty() ? "" : " " + d.valueHint);
    s += std::format("  {:<26} {}", left, d.help);
    if (!d.defaultText.empty()) s += std::format(" [{}]", d.defaultText);
    s += '\n';
  }
  return s;
}

}  // namespace od
