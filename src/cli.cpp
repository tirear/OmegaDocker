// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "cli.hpp"

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
