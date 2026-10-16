// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "cli.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <functional>
#include <map>
#include <optional>

#include <toml.hpp>

#include "log.hpp"
#include "version.hpp"

namespace od {
namespace {

constexpr int kParamFormatVersion = 1;

using Tokens = std::vector<std::string>;

// Value type of a parameter: in the file a TOML value of this type; on the
// command line --<option>=<text>, arrays comma separated.
enum class Type { Real, Int, Bool, Str, RealArr, StrArr };

struct Param {
  std::string key;  // "<table>.<name>", or "<name>" at the top level
  Type type;
  std::string hint;
  std::string defaultText;
  std::string help;
  std::function<void(Options&, const Tokens&)> set;               // throws od::Error (message without the name)
  std::function<std::optional<Tokens>(const Options&)> get;        // std::nullopt when not set
};

std::string optionName(const std::string& key) {
  std::string s = key;
  std::replace(s.begin(), s.end(), '.', '_');
  return s;
}

double toReal(const std::string& s) {
  double v = 0;
  auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc() || p != s.data() + s.size()) fail("'{}' is not a number", s);
  if (!std::isfinite(v)) fail("'{}' is not finite", s);
  return v;
}

long long toInt(const std::string& s) {
  long long v = 0;
  auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc() || p != s.data() + s.size()) fail("'{}' is not an integer", s);
  return v;
}

int toInt32(const std::string& s) {
  const long long v = toInt(s);
  if (v < -2147483647LL || v > 2147483647LL) fail("'{}' is out of range", s);
  return static_cast<int>(v);
}

bool toBool(const std::string& s) {
  if (s == "true" || s == "1" || s == "yes") return true;
  if (s == "false" || s == "0" || s == "no") return false;
  fail("expected true/false, 1/0 or yes/no but got '{}'", s);
}

const std::string& one(const Tokens& t) {
  if (t.size() != 1) fail("expected one value");
  return t.front();
}

std::string fmtReal(double v) { return std::format("{}", v); }

// Accessor of an Options member, usable with Options& and const Options&.
#define OD_ACC(member) [](auto& o) -> auto& { return o.member; }

template <typename Acc>
Param real(std::string key, std::string hint, std::string def, std::string help, Acc acc) {
  return {std::move(key), Type::Real, std::move(hint), std::move(def), std::move(help),
          [acc](Options& o, const Tokens& t) { acc(o) = toReal(one(t)); },
          [acc](const Options& o) -> std::optional<Tokens> { return Tokens{fmtReal(acc(o))}; }};
}

template <typename Acc>
Param integer(std::string key, std::string hint, std::string def, std::string help, Acc acc) {
  return {std::move(key), Type::Int, std::move(hint), std::move(def), std::move(help),
          [acc](Options& o, const Tokens& t) {
            using T = std::remove_reference_t<decltype(acc(o))>;
            if constexpr (std::is_same_v<T, int>) acc(o) = toInt32(one(t));
            else if constexpr (std::is_same_v<T, std::optional<int>>) acc(o) = toInt32(one(t));
            else acc(o) = toInt(one(t));
          },
          [acc](const Options& o) -> std::optional<Tokens> {
            const auto& v = acc(o);
            using T = std::remove_cvref_t<decltype(v)>;
            if constexpr (std::is_same_v<T, std::optional<int>>) {
              if (!v) return std::nullopt;
              return Tokens{std::to_string(*v)};
            } else {
              return Tokens{std::to_string(v)};
            }
          }};
}

template <typename Acc>
Param boolean(std::string key, std::string def, std::string help, Acc acc) {
  return {std::move(key), Type::Bool, "true|false", std::move(def), std::move(help),
          [acc](Options& o, const Tokens& t) { acc(o) = toBool(one(t)); },
          [acc](const Options& o) -> std::optional<Tokens> { return Tokens{acc(o) ? "true" : "false"}; }};
}

template <typename Acc>
Param vec3(std::string key, std::string help, Acc acc) {
  return {std::move(key), Type::RealArr, "x,y,z", "required", std::move(help),
          [acc](Options& o, const Tokens& t) {
            if (t.size() != 3) fail("expected 3 numbers x,y,z");
            acc(o) = Vec3{toReal(t[0]), toReal(t[1]), toReal(t[2])};
          },
          [acc](const Options& o) -> std::optional<Tokens> {
            const auto& v = acc(o);
            if (!v) return std::nullopt;
            return Tokens{fmtReal((*v)[0]), fmtReal((*v)[1]), fmtReal((*v)[2])};
          }};
}

// 9a sample numbers and csgm parameters of one set (spec 7.0, 7.2, 7.8).
void addSetParams(std::vector<Param>& t, const std::string& set) {
  const bool stam = set == "stam";
  const std::string who = stam ? "StaM group: " : "Every MobM: ";
  auto acc = [stam](auto& o) -> auto& { return stam ? o.stam : o.mobm; };
  auto csgm = [stam](auto& o) -> auto& { return stam ? o.stamCsgm : o.mobmCsgm; };
  if (stam) {
    t.push_back(real("stam.flexMargin", "<dist>", "4",
                     "StaM torsions are rotatable by default within the box extended by this margin (Å)",
                     OD_ACC(stamFlexMargin)));
    t.push_back(boolean("stam.fixedResMovable", "false",
                        "false keeps the backbone of fixed residues in place; true lets the whole chain move by torsions",
                        OD_ACC(stamFixedResMovable)));
  }
  t.push_back(integer(set + ".dockingConfSampleNum", "<n>", "required if rotatable",
                      who + "conformation samples docked per round (dcss size, >= 4)",
                      [acc](auto& o) -> auto& { return acc(o).docking; }));
  t.push_back(integer(set + ".minL2confSampleNum", "<n>", "2 x dockingConfSampleNum",
                      who + "9a: refill the L2 pool to at least this many samples",
                      [acc](auto& o) -> auto& { return acc(o).minL2; }));
  t.push_back(integer(set + ".initL1confSampleNum", "<n>", "minL2confSampleNum",
                      who + "9a: samples generated in the first batch", [acc](auto& o) -> auto& { return acc(o).initL1; }));
  t.push_back(integer(set + ".eliteConfSampleNum", "<n>", "max(2, ceil(dockingConfSampleNum / 4))",
                      who + "csgm1 is used while the dcss has fewer samples",
                      [acc](auto& o) -> auto& { return acc(o).elite; }));
  t.push_back(integer(set + ".maxL1TotalNum", "<n>", "1000000000", who + "9a: largest number of samples generated in one refill",
                      [csgm](auto& o) -> auto& { return csgm(o).maxL1TotalNum; }));
  t.push_back(integer(set + ".maxL1BatchNum", "<n>", "65536", who + "9a: largest number of samples processed at once",
                      [csgm](auto& o) -> auto& { return csgm(o).maxL1BatchNum; }));
  t.push_back({set + ".csgmBuild", Type::Str, "whole|grow", "whole",
               who + "9a: apply all torsions at once, or build segment by segment",
               [csgm](Options& o, const Tokens& v) {
                 const std::string& s = one(v);
                 if (s != "whole" && s != "grow") fail("expected whole or grow but got '{}'", s);
                 csgm(o).grow = s == "grow";
               },
               [csgm](const Options& o) -> std::optional<Tokens> { return Tokens{csgm(o).grow ? "grow" : "whole"}; }});
  t.push_back({set + ".csgmGA", Type::Str, "csgm2", "csgm2",
               who + "9a: method once the dcss has eliteConfSampleNum samples (csgm2 or csgm3)",
               [csgm](Options& o, const Tokens& v) {
                 const std::string& s = one(v);
                 if (s != "csgm2" && s != "csgm3") fail("expected csgm2 or csgm3 but got '{}'", s);
                 csgm(o).pathExploring = s == "csgm3";
               },
               [csgm](const Options& o) -> std::optional<Tokens> {
                 return Tokens{csgm(o).pathExploring ? "csgm3" : "csgm2"};
               }});
  t.push_back(real(set + ".tournamentRate", "<%>", "60", who + "csgm2: probability (%) that the fitter of two candidates is chosen",
                   [csgm](auto& o) -> auto& { return csgm(o).tournamentRate; }));
  t.push_back(real(set + ".crossoverRate", "<%>", "80", who + "csgm2: probability (%) of two-point crossover",
                   [csgm](auto& o) -> auto& { return csgm(o).crossoverRate; }));
  t.push_back(real(set + ".mutationRate", "<%>", "2", who + "csgm2: mutation probability (%) per torsion",
                   [csgm](auto& o) -> auto& { return csgm(o).mutationRate; }));
  t.push_back(real(set + ".absMaxDang", "<deg>", "90", who + "csgm2: largest torsion change of a mutation",
                   [csgm](auto& o) -> auto& { return csgm(o).absMaxDang; }));
  t.push_back(integer(set + ".peMaxOffset", "<deg>", "30", who + "csgm3: largest offset of a line end point (integer degrees)",
                      [csgm](auto& o) -> auto& { return csgm(o).peMaxOffset; }));
  t.push_back(real(set + ".peCoverAngle", "<deg>", "360", who + "csgm3: span of a line along its largest component",
                   [csgm](auto& o) -> auto& { return csgm(o).peCoverAngle; }));
  t.push_back(integer(set + ".peSampleNumPerLine", "<n>", "4", who + "csgm3: samples per line",
                      [csgm](auto& o) -> auto& { return csgm(o).peSampleNumPerLine; }));
  t.push_back(integer(set + ".growMaxRetry", "<n>", "16", who + "csgmBuild = grow, csgm1: resamplings of a failing segment before discarding",
                      [csgm](auto& o) -> auto& { return csgm(o).growMaxRetry; }));
  t.push_back(real(set + ".csgm1AngleStep", "<deg>", "1", who + "csgm1a / R: sampled angles are multiples of this step",
                   [csgm](auto& o) -> auto& { return csgm(o).csgm1AngleStep; }));
  t.push_back(real(set + ".ramaFloor", "<w>", "1e-4", who + "csgm1b: weight added to every Ramachandran cell",
                   [csgm](auto& o) -> auto& { return csgm(o).ramaFloor; }));
}

const std::vector<Param>& params() {
  static const std::vector<Param> table = [] {
    std::vector<Param> t;
    t.push_back({"seed", Type::Int, "<n>", "random", "Random seed (written to the log and to params.used.toml)",
                 [](Options& o, const Tokens& v) {
                   const long long s = toInt(one(v));
                   if (s < 0) fail("must be >= 0");
                   o.seed = static_cast<std::uint64_t>(s);
                 },
                 [](const Options& o) -> std::optional<Tokens> {
                   if (!o.seed) return std::nullopt;
                   return Tokens{std::to_string(*o.seed)};
                 }});
    t.push_back({"ccd", Type::StrArr, "<file>,...", "none", "Additional CCD definitions (mmCIF)",
                 [](Options& o, const Tokens& v) { o.ccd = v; },
                 [](const Options& o) -> std::optional<Tokens> { return o.ccd; }});

    t.push_back(vec3("box.center", "Docking box center (Å)", OD_ACC(boxCenter)));
    t.push_back(vec3("box.size", "Docking box size (Å)", OD_ACC(boxSize)));

    t.push_back(real("grid.spacing", "<dist>", "0.375", "9b: spacing of the energy and distance grids (Å)", OD_ACC(gridSpacing)));

    t.push_back(boolean("output.pdb", "true", "Write PDB output", OD_ACC(outPdb)));
    t.push_back(boolean("output.cif", "false", "Write mmCIF output", OD_ACC(outCif)));
    t.push_back(boolean("output.sdf", "false", "Write SDF output", OD_ACC(outSdf)));
    t.push_back(boolean("output.mol2", "false", "Write mol2 output", OD_ACC(outMol2)));
    t.push_back(integer("output.minPerRank", "<n>", "40", "9d: output at least the first n samples of every ranking",
                        OD_ACC(outMinPerRank)));

    t.push_back(real("topology.omegaTolerance", "<deg>", "30", "Allowed deviation of a fixed omega from 0 or 180 degrees",
                     OD_ACC(omegaTolerance)));
    t.push_back(boolean("topology.acceptNonPlanarOmega", "false",
                        "Continue when a fixed omega exceeds topology.omegaTolerance", OD_ACC(acceptNonPlanarOmega)));
    t.push_back(real("topology.maxBreakDist", "<dist>", "0.5", "Largest closure deviation of a broken bond kept after 9a (Å)",
                     OD_ACC(maxBreakDist)));
    t.push_back({"topology.reachPerRes", Type::RealArr, "<p>,<n>", "3.8,7.0",
                 "csgmBuild = grow: reach per residue of protein / nucleic-acid chains toward a closure (Å)",
                 [](Options& o, const Tokens& v) {
                   if (v.size() != 2) fail("expected 2 numbers <protein>,<nucleic>");
                   o.reachProtein = toReal(v[0]);
                   o.reachNucleic = toReal(v[1]);
                 },
                 [](const Options& o) -> std::optional<Tokens> {
                   return Tokens{fmtReal(o.reachProtein), fmtReal(o.reachNucleic)};
                 }});

    t.push_back(real("energy.closureG", "<k>", "200", "Closure constraint stiffness G (kcal/mol/Å²)", OD_ACC(closureG)));
    t.push_back(real("energy.clashDistance", "<dist>", "1.8", "Heavy-atom pairs closer than this clash (Å)", OD_ACC(clashDistance)));
    t.push_back(real("energy.Emax", "<E>", "1000", "9a: conformations with energy below this enter the L2 pool (kcal/mol)",
                     OD_ACC(Emax)));
    t.push_back(real("energy.EclashBase", "<E>", "10000", "Constant part of the clash penalty; must exceed energy.Emax (kcal/mol)",
                     OD_ACC(EclashBase)));
    t.push_back(real("energy.Kclash", "<k>", "1000", "Quadratic part of the clash penalty (kcal/mol/Å²)", OD_ACC(Kclash)));

    addSetParams(t, "stam");
    addSetParams(t, "mobm");

    t.push_back(integer("conf9a.minMaxStep", "<n>", "100", "9a minimization: largest number of steps", OD_ACC(confMinMaxStep)));
    t.push_back(real("conf9a.minDeltaE", "<E>", "10", "9a minimization: stop when one step changes the energy by less (kcal/mol)",
                     OD_ACC(confMinDeltaE)));
    t.push_back(real("conf9a.maxDihedralStep", "<deg>", "10", "9a minimization: largest change of one torsion per step",
                     OD_ACC(confMaxDihedralStep)));

    t.push_back({"dock9b.posm", Type::Str, "posm1|posm2", "required", "9b: position and orientation sampling method",
                 [](Options& o, const Tokens& v) {
                   const std::string& s = one(v);
                   if (s == "posm2_fft" || s == "posm3" || s == "posm4") fail("{} is not implemented yet", s);
                   if (s != "posm1" && s != "posm2") fail("expected posm1 or posm2 but got '{}'", s);
                   o.posm = s;
                 },
                 [](const Options& o) -> std::optional<Tokens> {
                   if (o.posm.empty()) return std::nullopt;
                   return Tokens{o.posm};
                 }});
    t.push_back(real("dock9b.posStep", "<dist>", "1", "9b: spacing of the sampled positions of the MobM heavy-atom centre (Å)",
                     OD_ACC(posStep)));
    t.push_back(integer("dock9b.posm1SampleNum", "<n>", "1000", "9b, posm1: L1poss per conformation pair", OD_ACC(posm1SampleNum)));
    t.push_back(real("dock9b.rotStep", "<deg>", "7.5", "9b: spacing of the Hopf orientation grid (the closest level is used)",
                     OD_ACC(rotStep)));
    t.push_back(real("dock9b.maxNearDistance", "<dist>", "6", "9b: poses whose heavy atoms are all farther from the StaM are removed (Å)",
                     OD_ACC(maxNearDistance)));
    t.push_back(integer("dock9b.maxOverlapAtomNum", "<n>", "1",
                        "9b: poses with this many MobM heavy atoms closer than energy.clashDistance to the StaM are removed",
                        OD_ACC(maxOverlapAtomNum)));
    t.push_back(integer("dock9b.minMaxStep", "<n>", "1000", "9b minimization: largest number of steps", OD_ACC(dockMinMaxStep)));
    t.push_back(integer("dock9b.minConvSteps", "<n>", "10",
                        "9b minimization: stop after this many consecutive steps with |dE| < dock9b.minDeltaE",
                        OD_ACC(dockMinConvSteps)));
    t.push_back(real("dock9b.minDeltaE", "<E>", "1e-3", "9b minimization: energy change regarded as converged (kcal/mol)",
                     OD_ACC(dockMinDeltaE)));
    t.push_back(integer("dock9b.lbfgsMemory", "<n>", "7", "L-BFGS: number of stored steps", OD_ACC(lbfgsMemory)));
    t.push_back(real("dock9b.scaleTrans", "<s>", "1", "9b minimization: scale of the translation variables (Å)", OD_ACC(scaleTrans)));
    t.push_back(real("dock9b.scaleRot", "<s>", "1", "9b minimization: scale of the rotation variables (rad)", OD_ACC(scaleRot)));
    t.push_back(real("dock9b.scaleTors", "<s>", "1", "9b minimization: scale of the torsion variables (rad)", OD_ACC(scaleTors)));
    t.push_back(real("dock9b.maxTransStep", "<dist>", "0.5", "9b minimization: largest translation per step (Å)",
                     OD_ACC(maxTransStep)));
    t.push_back(real("dock9b.maxAngleStep", "<deg>", "10", "9b minimization: largest rotation and torsion change per step",
                     OD_ACC(maxAngleStep)));

    t.push_back(integer("dock9c.maxRound", "<n>", "10", "Number of 9a-9c rounds", OD_ACC(maxRound)));
    t.push_back(integer("dock9c.topSampleNum", "<n>", "eliteStaM x eliteMobM", "9c: size of the topSampleSet of each MobM",
                        OD_ACC(topSampleNum)));
    t.push_back(real("dock9c.diversityRmsd", "<dist>", "2.0",
                     "9c: samples with the same StaM conformation closer than this RMSD to a kept sample are skipped (Å; 0 disables)",
                     OD_ACC(diversityRmsd)));
    t.push_back(integer("dock9c.sampleBufferNum", "<n>", "65536",
                        "9b/9c: run 9c whenever this many docking samples are waiting (at least max(4, midTopSampleFactor + 2) x topSampleNum)",
                        OD_ACC(sampleBufferNum)));
    t.push_back(integer("dock9c.midTopSampleFactor", "<n>", "1 (2 with midKeepSkipped)",
                        "9c: the 9c runs before the last one of a round keep this many x topSampleNum samples",
                        OD_ACC(midTopSampleFactor)));
    t.push_back(boolean("dock9c.midKeepSkipped", "false",
                        "9c: the 9c runs before the last one of a round also keep the skipped samples up to the last kept one "
                        "(at most sampleBufferNum samples in total; the kept samples are always kept)",
                        OD_ACC(midKeepSkipped)));

    t.push_back(integer("dock9d.symMaxPerm", "<n>", "256",
                        "9d: symmetry-corrected RMSD over the graph automorphisms when there are at most this many; otherwise "
                        "one-to-one pairing within each AD4 type",
                        OD_ACC(symMaxPerm)));
    t.push_back({"dock9d.freeEnergy", Type::StrArr, "<list>", "od_0100",
                 "9d: binding free energies to compute and rank by: adgpu_16, od_0100 (none or an empty list: neither)",
                 [](Options& o, const Tokens& v) {
                   o.feAdgpu16 = o.feOd0100 = false;
                   if (v.size() == 1 && v.front() == "none") return;
                   for (const auto& m : v) {
                     if (m == "adgpu_16") o.feAdgpu16 = true;
                     else if (m == "od_0100") o.feOd0100 = true;
                     else fail("unknown method '{}' (adgpu_16, od_0100, none)", m);
                   }
                 },
                 [](const Options& o) -> std::optional<Tokens> {
                   Tokens v;
                   if (o.feAdgpu16) v.push_back("adgpu_16");
                   if (o.feOd0100) v.push_back("od_0100");
                   return v;
                 }});

    t.push_back(real("od_0100.c0", "<c>", "0", "constant term (kcal/mol)", OD_ACC(fe.c0)));
    t.push_back(real("od_0100.tors", "<c>", "0", "coefficient of the MobM rotatable torsion number", OD_ACC(fe.tors)));
    t.push_back(real("od_0100.temperature", "<K>", "298.15", "temperature T of -TS = -RT sum ln f_i (K)", OD_ACC(fe.temperature)));
    t.push_back(integer("od_0100.unboundCandidateNum", "<n>", "20",
                        "conformation samples with the lowest 9a energy evaluated for the unbound free energy",
                        OD_ACC(fe.unboundCandidateNum)));
    for (const char* g : {"inter", "intra"})
      for (int k = 0; k < FeParams::kTerms; ++k) {
        const bool inter = g == std::string("inter");
        t.push_back(real(std::format("od_0100.{}.{}", g, kFeTermNames[k]), "<c>", "0.7564",
                         std::format("coefficient of the {} {} energy", inter ? "intermolecular" : "intramolecular deformation",
                                     kFeTermNames[k]),
                         [inter, k](auto& o) -> auto& { return (inter ? o.fe.inter : o.fe.intra)[k]; }));
      }
    t.push_back(real("od_0100.entropy.stam", "<c>", "0.9", "coefficient of -TS of the StaM group", OD_ACC(fe.entropyStam)));
    t.push_back(real("od_0100.entropy.mobm", "<c>", "0.9", "coefficient of -TS of the MobM", OD_ACC(fe.entropyMobm)));
    t.push_back(real("od_0100.shaking.radius", "<dist>", "1.4", "radius of the sphere each heavy atom is moved over (Å)",
                     OD_ACC(fe.shakingRadius)));
    t.push_back(real("od_0100.shaking.step", "<dist>", "0.2", "lattice spacing of the points in the sphere (Å)",
                     OD_ACC(fe.shakingStep)));
    t.push_back(real("od_0100.shaking.rrtScaling", "<c>", "4",
                     "reverse RT scaling of the shaking weights w_k = exp(-(E_k - E_min) x rrtScaling) (mol/kcal; not calibrated)",
                     OD_ACC(fe.rrtScaling)));
    return t;
  }();
  return table;
}

#undef OD_ACC

const Param* findParam(const std::string& key) {
  for (const Param& p : params())
    if (p.key == key) return &p;
  return nullptr;
}

// Options that exist only on the command line.
struct CliOption {
  std::string name, hint, help;
  bool flag;
  std::function<void(Options&, const std::string&)> apply;
};

const std::vector<CliOption>& cliOptions() {
  static const std::vector<CliOption> t{
      {"paramFile", "<file>", "Parameter file (TOML); the options below override it; may be given once", false,
       [](Options& o, const std::string& v) {
         if (!o.paramFile.empty()) fail("--paramFile may be given only once");
         o.paramFile = v;
       }},
      {"stamFile", "<file>", "Stationary molecules (StaM); required; may be given once", false,
       [](Options& o, const std::string& v) {
         if (!o.stamFile.empty()) fail("--stamFile may be given only once");
         o.stamFile = v;
       }},
      {"mobmFiles", "<file>", "Mobile molecule (MobM), one molecule per file; required; repeatable", false,
       [](Options& o, const std::string& v) { o.mobmFiles.push_back(v); }},
      {"molConfig", "<file>", "Molecule configuration file (TOML)", false,
       [](Options& o, const std::string& v) { o.molConfig = v; }},
      {"outPrefix", "<path>", "Output directory; required", false, [](Options& o, const std::string& v) { o.outPrefix = v; }},
      {"help", "", "Show this help", true, [](Options& o, const std::string&) { o.showHelp = true; }},
      {"version", "", "Show version", true, [](Options& o, const std::string&) { o.showVersion = true; }},
  };
  return t;
}

Tokens splitComma(const std::string& s) {
  Tokens t;
  if (s.empty()) return t;
  std::size_t start = 0;
  while (true) {
    const std::size_t c = s.find(',', start);
    t.push_back(s.substr(start, c == std::string::npos ? std::string::npos : c - start));
    if (c == std::string::npos) break;
    start = c + 1;
  }
  return t;
}

// TOML value to tokens, checking its TOML type.
Tokens tomlTokens(const toml::node& v, Type type) {
  auto scalar = [](const toml::node& n, Type ty) -> std::string {
    switch (ty) {
      case Type::Real:
        if (n.is_integer()) return std::to_string(*n.value<std::int64_t>());
        if (n.is_floating_point()) return fmtReal(*n.value<double>());
        fail("must be a number");
      case Type::Int:
        if (n.is_integer()) return std::to_string(*n.value<std::int64_t>());
        fail("must be an integer");
      case Type::Bool:
        if (n.is_boolean()) return *n.value<bool>() ? "true" : "false";
        fail("must be true or false");
      default:
        if (n.is_string()) return *n.value<std::string>();
        fail("must be a string");
    }
  };
  if (type == Type::RealArr || type == Type::StrArr) {
    if (!v.is_array()) fail("must be an array");
    Tokens t;
    for (const auto& e : *v.as_array()) t.push_back(scalar(e, type == Type::RealArr ? Type::Real : Type::Str));
    return t;
  }
  return {scalar(v, type)};
}

int lineOf(const toml::node& n) { return static_cast<int>(n.source().begin.line); }

void walk(const toml::table& t, const std::string& where, Options& o, std::multimap<int, std::string>& problems) {
  for (const auto& [kk, v] : t) {
    const std::string key = where.empty() ? std::string(kk.str()) : where + "." + std::string(kk.str());
    auto add = [&](const std::string& msg) {
      problems.emplace(lineOf(v), std::format("line {}: {}: {}", lineOf(v), key, msg));
    };
    if (key == "formatVersion") {
      if (!v.is_integer() || *v.value<std::int64_t>() != kParamFormatVersion)
        add(std::format("must be {}", kParamFormatVersion));
      continue;
    }
    bool table = false;
    for (const Param& p : params()) table |= p.key.starts_with(key + ".");
    if (table) {
      if (v.is_table()) walk(*v.as_table(), key, o, problems);
      else add("must be a table");
      continue;
    }
    const Param* p = findParam(key);
    if (!p) {
      add("unknown key");
      continue;
    }
    try {
      p->set(o, tomlTokens(v, p->type));
    } catch (const Error& e) {
      add(e.what());
    }
  }
}

std::string tomlString(const std::string& s) {
  std::string r = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') r += '\\';
    r += c;
  }
  return r + "\"";
}

std::string tomlValue(const Param& p, const Tokens& t) {
  switch (p.type) {
    case Type::Str: return tomlString(t.front());
    case Type::RealArr:
    case Type::StrArr: {
      std::string s = "[";
      for (std::size_t i = 0; i < t.size(); ++i)
        s += (i ? ", " : "") + (p.type == Type::StrArr ? tomlString(t[i]) : t[i]);
      return s + "]";
    }
    default: return t.front();
  }
}

// Name of a parameter in messages.
std::string nm(const std::string& key) { return std::format("{} (--{})", key, optionName(key)); }

}  // namespace

Options parseCommandLine(int argc, const char* const* argv) {
  std::map<std::string, const CliOption*> cli;
  for (const auto& c : cliOptions()) cli[c.name] = &c;
  std::map<std::string, const Param*> byOption;
  for (const auto& p : params()) byOption[optionName(p.key)] = &p;

  Options o;
  std::vector<std::tuple<const Param*, std::string, std::string>> overrides;  // param, argument, value
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) fail("Unexpected argument '{}' (options are --<name>=<value>)", arg);
    const auto eq = arg.find('=');
    const std::string name = arg.substr(2, eq == std::string::npos ? std::string::npos : eq - 2);
    const auto c = cli.find(name);
    const auto p = byOption.find(name);
    if (c == cli.end() && p == byOption.end()) fail("Unknown option --{}", name);
    if (c != cli.end() && c->second->flag) {
      if (eq != std::string::npos) fail("--{} does not take a value", name);
      c->second->apply(o, "");
      continue;
    }
    const std::string hint = c != cli.end() ? c->second->hint : p->second->hint;
    if (eq == std::string::npos) fail("--{} requires a value: --{}={}", name, name, hint);
    const std::string value = arg.substr(eq + 1);
    if (c != cli.end()) c->second->apply(o, value);
    else overrides.emplace_back(p->second, arg, value);
  }
  if (o.showHelp || o.showVersion) return o;

  std::vector<std::string> problems;
  if (!o.paramFile.empty()) {
    if (!std::filesystem::exists(o.paramFile)) fail("{}: file not found", o.paramFile);
    toml::table root;
    try {
      root = toml::parse_file(o.paramFile);
    } catch (const toml::parse_error& ex) {
      fail("{}: line {}: {}", o.paramFile, ex.source().begin.line, ex.description());
    }
    std::multimap<int, std::string> p;
    walk(root, "", o, p);
    for (auto& [line, s] : p) problems.push_back(o.paramFile + ": " + s);
  }
  for (const auto& [p, arg, value] : overrides) {
    try {
      const bool array = p->type == Type::RealArr || p->type == Type::StrArr;
      p->set(o, array ? splitComma(value) : Tokens{value});
    } catch (const Error& e) {
      problems.push_back(std::format("command line: {}: {}", arg, e.what()));
    }
  }
  if (!problems.empty()) {
    std::string msg = std::format("parameters: {} problem(s) found:", problems.size());
    for (const auto& s : problems) msg += "\n  " + s;
    fail("{}", msg);
  }
  return o;
}

std::string paramsText(const Options& o) {
  std::string s = std::format(
      "# OmegaDocker {} parameter file (spec 14). Keys not given keep the defaults;\n"
      "# --<table>_<key>=<value> on the command line overrides <table>.<key>.\n"
      "formatVersion = {}\n",
      kVersion, kParamFormatVersion);
  std::string table;
  for (const Param& p : params()) {
    const auto dot = p.key.rfind('.');
    const std::string t = dot == std::string::npos ? "" : p.key.substr(0, dot);
    if (t != table) {
      if (t == "od_0100")
        s += "\n# od_0100 (spec 10.2; not calibrated yet):\n"
             "#   dG = c0 + sum_t inter.t * E_inter,t + sum_t intra.t * dE_intra,t\n"
             "#        + entropy.stam * d(-TS)_StaM + entropy.mobm * d(-TS)_MobM + tors * N_tors\n"
             "#   t = vdw, hbond, elec, desolv (AD4 terms, kcal/mol, AD4 weights included)\n"
             "#   d(-TS)_X = -R * temperature * sum ln f_i over the heavy atoms of X, complex - unbound\n"
             "#   f_i = sum_k exp(-(E_k - E_min) * shaking.rrtScaling) / points within shaking.radius\n";
      s += std::format("\n[{}]\n", t);
      table = t;
    }
    const std::string name = dot == std::string::npos ? p.key : p.key.substr(dot + 1);
    s += std::format("# {} [{}]\n", p.help, p.defaultText);
    if (const auto v = p.get(o)) s += std::format("{} = {}\n", name, tomlValue(p, *v));
    else s += std::format("# {} =\n", name);
  }
  return s;
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
  if (n.docking < 4) fail("{} must be >= 4", nm(p + ".dockingConfSampleNum"));
  if (n.minL2 < 2 * n.docking)
    fail("{} must be >= 2 x {}.dockingConfSampleNum ({})", nm(p + ".minL2confSampleNum"), p, 2 * n.docking);
  if (n.initL1 < n.minL2) fail("{} must be >= {}.minL2confSampleNum ({})", nm(p + ".initL1confSampleNum"), p, n.minL2);
  if (n.elite < 2 || n.elite > n.docking / 2)
    fail("{} must be within [2, {}.dockingConfSampleNum / 2] = [2, {}]", nm(p + ".eliteConfSampleNum"), p, n.docking / 2);
}

}  // namespace

void requireConfSampleNums(const Options& o, SampleSet set) {
  const bool stam = set == SampleSet::StaM;
  if (!(stam ? o.stam : o.mobm).docking)
    fail("{} is required: {} rotatable torsions", nm(stam ? "stam.dockingConfSampleNum" : "mobm.dockingConfSampleNum"),
         stam ? "the StaM group has" : "a MobM has");
  checkConfSampleNums(o, set);
}

void validateOptions(const Options& o) {
  if (o.stamFile.empty()) fail("--stamFile is required");
  if (o.mobmFiles.empty()) fail("At least one --mobmFiles is required");
  if (o.outPrefix.empty()) fail("--outPrefix is required");
  if (!o.boxCenter) fail("{} is required", nm("box.center"));
  if (!o.boxSize) fail("{} is required", nm("box.size"));
  for (double v : *o.boxSize)
    if (!(v > 0)) fail("{}: all components must be > 0", nm("box.size"));
  if (o.posm.empty()) fail("{} is required", nm("dock9b.posm"));
  auto check = [](bool ok, const std::string& key, const std::string& rule) {
    if (!ok) fail("{} must be {}", nm(key), rule);
  };
  check(o.omegaTolerance >= 0 && o.omegaTolerance <= 90, "topology.omegaTolerance", "within [0, 90]");
  check(o.stamFlexMargin >= 0, "stam.flexMargin", ">= 0");
  check(o.maxBreakDist > 0, "topology.maxBreakDist", "> 0");
  check(o.reachProtein > 0 && o.reachNucleic > 0, "topology.reachPerRes", "> 0 (both values)");
  check(o.closureG > 0, "energy.closureG", "> 0");
  check(o.clashDistance > 0, "energy.clashDistance", "> 0");
  check(o.Kclash >= 0, "energy.Kclash", ">= 0");
  check(o.EclashBase > o.Emax, "energy.EclashBase", "greater than energy.Emax");
  check(o.gridSpacing > 0, "grid.spacing", "> 0");
  check(o.maxNearDistance > 0, "dock9b.maxNearDistance", "> 0");
  check(o.rotStep > 0, "dock9b.rotStep", "> 0");
  check(o.posStep > 0, "dock9b.posStep", "> 0");
  check(o.posm1SampleNum >= 1, "dock9b.posm1SampleNum", ">= 1");
  check(o.maxOverlapAtomNum >= 1, "dock9b.maxOverlapAtomNum", ">= 1");
  check(o.dockMinMaxStep >= 0, "dock9b.minMaxStep", ">= 0");
  check(o.dockMinConvSteps >= 1, "dock9b.minConvSteps", ">= 1");
  check(o.dockMinDeltaE >= 0, "dock9b.minDeltaE", ">= 0");
  check(o.lbfgsMemory >= 1, "dock9b.lbfgsMemory", ">= 1");
  check(o.scaleTrans > 0, "dock9b.scaleTrans", "> 0");
  check(o.scaleRot > 0, "dock9b.scaleRot", "> 0");
  check(o.scaleTors > 0, "dock9b.scaleTors", "> 0");
  check(o.maxTransStep > 0, "dock9b.maxTransStep", "> 0");
  check(o.maxAngleStep > 0, "dock9b.maxAngleStep", "> 0");
  check(o.maxRound >= 1, "dock9c.maxRound", ">= 1");
  check(!o.topSampleNum || *o.topSampleNum >= 1, "dock9c.topSampleNum", ">= 1");
  check(o.diversityRmsd >= 0, "dock9c.diversityRmsd", ">= 0");
  check(o.sampleBufferNum >= 1, "dock9c.sampleBufferNum", ">= 1");
  check(!o.midTopSampleFactor || *o.midTopSampleFactor >= 1, "dock9c.midTopSampleFactor", ">= 1");
  check(o.symMaxPerm >= 0, "dock9d.symMaxPerm", ">= 0");
  check(o.outMinPerRank >= 1, "output.minPerRank", ">= 1");
  check(o.outPdb || o.outCif || o.outSdf || o.outMol2, "output.pdb, output.cif, output.sdf or output.mol2", "true (all output formats are disabled)");
  for (SampleSet set : {SampleSet::StaM, SampleSet::MobM}) {
    const ConfSampleOpts& s = set == SampleSet::StaM ? o.stam : o.mobm;
    const std::string p = set == SampleSet::StaM ? "stam" : "mobm";
    if (s.docking) checkConfSampleNums(o, set);
    else if (s.minL2 || s.initL1 || s.elite)
      fail("{} is required when other {}.*ConfSampleNum values are given", nm(p + ".dockingConfSampleNum"), p);
  }
  for (const auto& [set, c] : {std::pair{std::string("stam"), &o.stamCsgm}, std::pair{std::string("mobm"), &o.mobmCsgm}}) {
    check(c->maxL1TotalNum >= 1, set + ".maxL1TotalNum", ">= 1");
    check(c->maxL1BatchNum >= 1, set + ".maxL1BatchNum", ">= 1");
    check(c->growMaxRetry >= 0, set + ".growMaxRetry", ">= 0");
    check(c->csgm1AngleStep > 0 && c->csgm1AngleStep <= 180, set + ".csgm1AngleStep", "within (0, 180]");
    check(c->ramaFloor >= 0, set + ".ramaFloor", ">= 0");
    check(c->tournamentRate >= 0 && c->tournamentRate <= 100, set + ".tournamentRate", "within [0, 100]");
    check(c->crossoverRate >= 0 && c->crossoverRate <= 100, set + ".crossoverRate", "within [0, 100]");
    check(c->mutationRate >= 0 && c->mutationRate <= 100, set + ".mutationRate", "within [0, 100]");
    check(c->absMaxDang > 0 && c->absMaxDang <= 180, set + ".absMaxDang", "within (0, 180]");
    check(c->peMaxOffset >= 0 && c->peMaxOffset <= 180, set + ".peMaxOffset", "within [0, 180]");
    check(c->peCoverAngle > 0 && c->peCoverAngle <= 360, set + ".peCoverAngle", "within (0, 360]");
    check(c->peSampleNumPerLine >= 1, set + ".peSampleNumPerLine", ">= 1");
  }
  check(o.confMinMaxStep >= 0, "conf9a.minMaxStep", ">= 0");
  check(o.confMinDeltaE >= 0, "conf9a.minDeltaE", ">= 0");
  check(o.confMaxDihedralStep > 0, "conf9a.maxDihedralStep", "> 0");
  check(o.fe.temperature > 0, "od_0100.temperature", "> 0");
  check(o.fe.unboundCandidateNum >= 1, "od_0100.unboundCandidateNum", ">= 1");
  check(o.fe.shakingRadius > 0, "od_0100.shaking.radius", "> 0");
  check(o.fe.shakingStep > 0, "od_0100.shaking.step", "> 0");
  check(o.fe.rrtScaling >= 0, "od_0100.shaking.rrtScaling", ">= 0");
}

std::string usageText() {
  std::string s = std::format(
      "OmegaDocker {}\n\n"
      "Usage: omegadocker --stamFile=<file> --mobmFiles=<file> [--mobmFiles=<file> ...] --outPrefix=<path>\n"
      "                   [--paramFile=<file>] [--<table>_<key>=<value> ...]\n\n"
      "Every option takes its value after '='. A parameter <table>.<key> of the parameter file is\n"
      "overridden by --<table>_<key>=<value>; arrays are comma separated; booleans are true/false, 1/0 or yes/no.\n"
      "Priority: built-in default < --paramFile < command line.\n\nCommand line only:\n",
      kVersion);
  auto line = [&](const std::string& left, const std::string& help, const std::string& def) {
    s += std::format("  {:<40} {}", left, help);
    if (!def.empty()) s += std::format(" [{}]", def);
    s += '\n';
  };
  for (const auto& c : cliOptions()) line("--" + c.name + (c.flag ? "" : "=" + c.hint), c.help, "");
  std::string table = "-";
  for (const Param& p : params()) {
    const auto dot = p.key.rfind('.');
    const std::string t = dot == std::string::npos ? "" : p.key.substr(0, dot);
    if (t != table) {
      s += t.empty() ? "\nTop level:\n" : std::format("\n[{}]:\n", t);
      table = t;
    }
    line("--" + optionName(p.key) + "=" + p.hint, p.help, p.defaultText);
  }
  return s;
}

}  // namespace od
