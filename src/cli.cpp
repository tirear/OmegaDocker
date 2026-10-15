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

const std::vector<OptionDef>& optionTable() {
  static const std::vector<OptionDef> table = {
      {"sta", Kind::Value, "<file>", "required", "Stationary molecules (StaM); may be given once",
       [](Options& o, const std::string& v) {
         if (!o.sta.empty()) fail("--sta may be given only once");
         o.sta = v;
       }},
      {"mob", Kind::Value, "<file>", "required", "Mobile molecule (MobM), one molecule per file; repeatable",
       [](Options& o, const std::string& v) { o.mob.push_back(v); }},
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
      {"staFlexMargin", Kind::Value, "<dist>", "4", "StaM torsions are rotatable by default within the box extended by this margin (Å)",
       [](Options& o, const std::string& v) { o.staFlexMargin = parseDouble("staFlexMargin", v); }},
      {"staFixedResMovable", Kind::Value, "true|false", "false",
       "StaM: false keeps the backbone of fixed residues in place; true lets the whole chain move by torsions",
       [](Options& o, const std::string& v) { o.staFixedResMovable = parseBool("staFixedResMovable", v); }},
      {"maxBreakDist", Kind::Value, "<dist>", "0.5", "Largest closure deviation of a broken bond kept after 9a (Å)",
       [](Options& o, const std::string& v) { o.maxBreakDist = parseDouble("maxBreakDist", v); }},
      {"closureG", Kind::Value, "<k>", "200", "Closure constraint stiffness G (kcal/mol/Å²)",
       [](Options& o, const std::string& v) { o.closureG = parseDouble("closureG", v); }},
      {"clashDistance", Kind::Value, "<dist>", "1.8", "Heavy-atom pairs closer than this clash (Å)",
       [](Options& o, const std::string& v) { o.clashDistance = parseDouble("clashDistance", v); }},
      {"Emax", Kind::Value, "<E>", "1000", "9a: conformations with energy below this enter the L2 pool (kcal/mol)",
       [](Options& o, const std::string& v) { o.Emax = parseDouble("Emax", v); }},
      {"EclashBase", Kind::Value, "<E>", "10000", "Constant part of the clash penalty; must exceed --Emax (kcal/mol)",
       [](Options& o, const std::string& v) { o.EclashBase = parseDouble("EclashBase", v); }},
      {"Kclash", Kind::Value, "<k>", "1000", "Quadratic part of the clash penalty (kcal/mol/Å²)",
       [](Options& o, const std::string& v) { o.Kclash = parseDouble("Kclash", v); }},
      {"dockingConfSampleNum", Kind::Value, "<n>", "required", "Conformation samples docked per round (dcss size, >= 4)",
       [](Options& o, const std::string& v) { o.dockingConfSampleNum = parseInt32("dockingConfSampleNum", v); }},
      {"minL2confSampleNum", Kind::Value, "<n>", "2 x docking", "9a: refill the L2 pool to at least this many samples",
       [](Options& o, const std::string& v) { o.minL2confSampleNum = parseInt32("minL2confSampleNum", v); }},
      {"initL1confSampleNum", Kind::Value, "<n>", "minL2", "9a: samples generated in the first batch",
       [](Options& o, const std::string& v) { o.initL1confSampleNum = parseInt32("initL1confSampleNum", v); }},
      {"eliteConfSampleNum", Kind::Value, "<n>", "max(2, ceil(docking/4))", "csgm1 is used while the dcss has fewer samples",
       [](Options& o, const std::string& v) { o.eliteConfSampleNum = parseInt32("eliteConfSampleNum", v); }},
      {"maxL1TotalNum", Kind::Value, "<n>", "1000000000", "9a: largest number of samples generated in one refill",
       [](Options& o, const std::string& v) { o.maxL1TotalNum = parseInt("maxL1TotalNum", v); }},
      {"maxL1BatchNum", Kind::Value, "<n>", "65536", "9a: largest number of samples processed at once",
       [](Options& o, const std::string& v) { o.maxL1BatchNum = parseInt("maxL1BatchNum", v); }},
      {"csgm1AngleStep", Kind::Value, "<deg>", "1", "csgm1a: angle changes are multiples of this",
       [](Options& o, const std::string& v) { o.csgm1AngleStep = parseDouble("csgm1AngleStep", v); }},
      {"ramaFloor", Kind::Value, "<w>", "1e-4", "csgm1b: weight added to every Ramachandran cell",
       [](Options& o, const std::string& v) { o.ramaFloor = parseDouble("ramaFloor", v); }},
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
  };
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

ConfSampleNums confSampleNums(const Options& o) {
  ConfSampleNums n;
  n.docking = o.dockingConfSampleNum.value_or(0);
  n.minL2 = o.minL2confSampleNum.value_or(2 * n.docking);
  n.initL1 = o.initL1confSampleNum.value_or(n.minL2);
  n.elite = o.eliteConfSampleNum.value_or(std::max(2, (n.docking + 3) / 4));
  return n;
}

void validateOptions(const Options& o) {
  if (o.sta.empty()) fail("--sta is required");
  if (o.mob.empty()) fail("At least one --mob is required");
  if (!o.boxCenter) fail("--boxCenter is required");
  if (!o.boxSize) fail("--boxSize is required");
  for (double v : *o.boxSize)
    if (!(v > 0)) fail("--boxSize: all components must be > 0");
  if (o.outPrefix.empty()) fail("--outPrefix is required");
  if (!(o.omegaTolerance >= 0 && o.omegaTolerance <= 90)) fail("--omegaTolerance must be within [0, 90]");
  if (!(o.staFlexMargin >= 0)) fail("--staFlexMargin must be >= 0");
  if (!(o.maxBreakDist > 0)) fail("--maxBreakDist must be > 0");
  if (!(o.closureG > 0)) fail("--closureG must be > 0");
  if (!(o.clashDistance > 0)) fail("--clashDistance must be > 0");
  if (!(o.Kclash >= 0)) fail("--Kclash must be >= 0");
  if (!(o.EclashBase > o.Emax)) fail("--EclashBase must be greater than --Emax");
  if (!o.dockingConfSampleNum) fail("--dockingConfSampleNum is required");
  const ConfSampleNums n = confSampleNums(o);
  if (n.docking < 4) fail("--dockingConfSampleNum must be >= 4");
  if (n.minL2 < 2 * n.docking) fail("--minL2confSampleNum must be >= 2 x --dockingConfSampleNum ({})", 2 * n.docking);
  if (n.initL1 < n.minL2) fail("--initL1confSampleNum must be >= --minL2confSampleNum ({})", n.minL2);
  if (n.elite < 2 || n.elite > n.docking / 2)
    fail("--eliteConfSampleNum must be within [2, --dockingConfSampleNum / 2] = [2, {}]", n.docking / 2);
  if (o.maxL1TotalNum < 1) fail("--maxL1TotalNum must be >= 1");
  if (o.maxL1BatchNum < 1) fail("--maxL1BatchNum must be >= 1");
  if (!(o.csgm1AngleStep > 0 && o.csgm1AngleStep <= 180)) fail("--csgm1AngleStep must be within (0, 180]");
  if (!(o.ramaFloor >= 0)) fail("--ramaFloor must be >= 0");
  if (o.confMinMaxStep < 0) fail("--confMinMaxStep must be >= 0");
  if (!(o.confMinDeltaE >= 0)) fail("--confMinDeltaE must be >= 0");
  if (!(o.confMaxDihedralStep > 0)) fail("--confMaxDihedralStep must be > 0");
  if (!o.outPdb && !o.outCif && !o.outSdf && !o.outMol2)
    fail("All output formats are disabled (--outPdb, --outCif, --outSdf, --outMol2)");
}

std::string usageText() {
  std::string s = std::format("OmegaDocker {}\n\nUsage: omegadocker --sta <file> --mob <file> [--mob <file> ...] "
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
