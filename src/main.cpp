// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <format>
#include <numeric>
#include <optional>
#include <random>
#include <span>

#include "ccd.hpp"
#include "cli.hpp"
#include "log.hpp"
#include "read_input.hpp"
#include "conf9a.hpp"
#include "csgm1.hpp"
#include "dssp.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "grow.hpp"
#include "orient.hpp"
#include "molconfig.hpp"
#include "topology.hpp"
#include "version.hpp"

namespace {

void printMolecule(const od::Molecule& m) {
  int netCharge = 0;
  double sumPartial = 0.0;
  for (const auto& a : m.atoms) {
    netCharge += a.formalCharge;
    sumPartial += a.partialCharge;
  }
  od::logInfo(std::format("  {:<40} atoms {:>7}  bonds {:>7}  residues {:>6}  formal charge {:>+4}  "
                          "sum partial charge {:>+8.3f}{}",
                          m.id, m.atoms.size(), m.bonds.size(), m.residues.size(), netCharge, sumPartial,
                          m.chargesFromFile ? " (file)" : " (Gasteiger)"));
}

od::TopologyParams topologyParams(const od::Options& opt) {
  od::TopologyParams p;
  p.box = od::Box{*opt.boxCenter, *opt.boxSize};
  p.stamFlexMargin = opt.stamFlexMargin;
  p.stamFixedResMovable = opt.stamFixedResMovable;
  return p;
}

od::Topology analyze(const od::Molecule& m, const od::Options& opt) {
  od::Topology t = od::analyzeTopology(m, topologyParams(opt));
  int rotatable = 0;
  for (const auto& tor : t.torsions) rotatable += tor.rotatable;
  od::logInfo(std::format("    backbone {} ({} units), torsions {}, rotatable by default {}",
                          od::backboneTypeName(t.backboneType), t.units.size(), t.torsions.size(), rotatable));
  return t;
}

// Spec 4.3, for omega angles set to fixed.
void checkOmegas(const od::Molecule& m, const od::MolSettings& s, const od::Options& opt) {
  const od::Topology& t = s.topology;
  std::vector<int> bad;
  for (int i : od::nonPlanarOmegas(t, opt.omegaTolerance))
    if (s.mode[i] == od::TorsionMode::Fixed) bad.push_back(i);
  if (bad.empty()) return;
  std::string list;
  for (int i : bad) {
    const od::Torsion& tor = t.torsions[i];
    const od::Atom& c = m.atoms[tor.atoms[1]];
    std::string where = c.residue >= 0 ? std::format("{}:{}{}:{}", m.residues[c.residue].info.chain,
                                                     m.residues[c.residue].info.seq,
                                                     m.residues[c.residue].info.iCode, m.residues[c.residue].info.name)
                                       : std::format("#{}", c.fileSerial);
    list += std::format("\n  {} omega {:.1f}", where, tor.angle);
  }
  if (!opt.acceptNonPlanarOmega)
    od::fail("{}: {} fixed omega angle(s) differ from 0 and 180 degrees by more than {} degrees "
             "(--acceptNonPlanarOmega continues):{}",
             m.id, bad.size(), opt.omegaTolerance, list);
  od::logWarn(std::format("{}: {} non-planar fixed omega angle(s):{}", m.id, bad.size(), list));
}

// Torsion tree, loops and added loop breaks (spec 4.7, 4.8).
void printTree(const od::Molecule& m, const od::MolSettings& s) {
  const od::Topology& t = s.topology;
  const od::TorsionTree& tree = s.tree;
  auto atomRef = [&](int a) {
    const od::Atom& x = m.atoms[a];
    if (x.residue < 0 || x.name.empty()) return std::format("#{}", x.fileSerial);
    const auto& r = m.residues[x.residue].info;
    return std::format("{}:{}{}:{}", r.chain, r.seq, r.iCode, x.name);
  };
  auto unitRef = [&](int u) {
    const int r = t.units[u].residue;
    return r >= 0 ? std::format("{}:{}{}", m.residues[r].info.chain, m.residues[r].info.seq, m.residues[r].info.iCode)
                  : std::format("unit {}", u + 1);
  };
  int depth = 0;
  std::vector<int> level(tree.torsions.size(), 0);
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    if (tree.torsions[i].parent >= 0) level[i] = level[tree.torsions[i].parent] + 1;
    depth = std::max(depth, level[i] + 1);
  }
  od::logInfo(std::format("    torsion tree: {} torsions, depth {}, anchor {}, closure pairs {}, distance restraints {}",
                          tree.torsions.size(), depth,
                          tree.anchors.size() == 1 ? atomRef(tree.anchors.front())
                                                   : std::format("{} fixed backbone atoms", tree.anchors.size()),
                          tree.virtuals.size(), tree.restraints.size()));
  for (const od::Loop& l : tree.loops) {
    std::string where = l.before >= 0 && l.after >= 0 ? "between fixed residues"
                        : l.before >= 0               ? "C-terminal, free end"
                        : l.after >= 0                ? "N-terminal, free end"
                                                      : "no fixed residue in the chain";
    std::string brk;
    if (l.breakBond >= 0)
      brk = std::format(", break [{}, {}]", atomRef(m.bonds[l.breakBond].a), atomRef(m.bonds[l.breakBond].b));
    od::logInfo(std::format("    loop {}-{} ({} residues, {}{})", unitRef(l.units.front()), unitRef(l.units.back()),
                            l.units.size(), where, brk));
  }
  for (const auto& c : s.chainBreaks)
    if (c.line == 0)
      od::logInfo(std::format("    added chainBreak [{}, {}] source = \"auto-loop\"", atomRef(m.bonds[c.bond].a),
                              atomRef(m.bonds[c.bond].b)));
}

od::EnergyParams energyParams(const od::Options& opt) {
  od::EnergyParams p;
  p.clashDistance = opt.clashDistance;
  p.clashBase = opt.EclashBase;
  p.clashK = opt.Kclash;
  p.closureG = opt.closureG;
  return p;
}

// 9a energy of the input conformation (spec 7.9).
void printInputEnergy(const std::string& name, const std::vector<od::EnergyPart>& parts, const od::Options& opt) {
  const od::EnergySystem sys(parts, energyParams(opt));
  int merged = 0, hd = 0;
  for (std::size_t k = 0; k < parts.size(); ++k)
    for (od::Ad4Type t : sys.typing(static_cast<int>(k)).type) {
      merged += t == od::Ad4Type::Merged;
      hd += t == od::Ad4Type::HD;
    }
  const auto e = sys.evaluate(sys.inputAngles());
  od::logInfo(std::format("  {:<40} torsions {:>5}  pairs {:>9}  merged H {:>5}  HD {:>5}", name, sys.torsionCount(),
                          sys.pairCount(), merged, hd));
  od::logInfo(std::format("    input: E {:.3f} (vdW {:.3f}  hbond {:.3f}  elec {:.3f}  desolv {:.3f}  clash {:.1f} "
                          "in {} pairs)  closure {:.3g}",
                          e.filterEnergy(), e.vdw, e.hbond, e.elec, e.desolv, e.clash, e.clashPairs, e.closure));
}

// First 9a refill of the L2 pool with csgm1 (spec 7.3).
void runFirst9a(const std::string& name, const std::vector<od::EnergyPart>& parts,
                const std::vector<const od::MolSettings*>& settings, const od::MolConfig& cfg, const od::Options& opt,
                od::SampleSet set, std::uint64_t seed, std::uint64_t stream) {
  const od::EnergySystem sys(parts, energyParams(opt));
  if (sys.torsionCount() == 0) {
    od::logInfo(std::format("  {}: no rotatable torsions; the input conformation is the only sample", name));
    return;
  }
  const od::CsgmOpts& csgm = set == od::SampleSet::StaM ? opt.stamCsgm : opt.mobmCsgm;
  od::Csgm1Params cp;
  cp.angleStep = csgm.csgm1AngleStep;
  cp.ramaFloor = csgm.ramaFloor;
  od::Csgm1GroupSampler sampler;
  for (std::size_t k = 0; k < parts.size(); ++k) sampler.add(od::Csgm1Sampler(*parts[k].mol, *settings[k], cfg, cp));
  const od::ConfSampleNums n = od::confSampleNums(opt, set);
  od::Fill9aParams p;
  p.minL2 = n.minL2;
  p.maxL1Batch = csgm.maxL1BatchNum;
  p.maxL1Total = csgm.maxL1TotalNum;
  p.Emax = opt.Emax;
  p.maxBreakDist = opt.maxBreakDist;
  p.min = {opt.confMinMaxStep, opt.confMinDeltaE, opt.confMaxDihedralStep};
  p.seed = seed;
  p.stream = stream;
  std::optional<od::GrowBuilder> grow;
  if (csgm.grow) {
    std::vector<const od::Topology*> topos;
    for (const auto* s : settings) topos.push_back(&s->topology);
    grow.emplace(sys, parts, topos, sampler.segments(),
                 od::GrowParams{csgm.growMaxRetry, opt.clashDistance, opt.reachProtein, opt.reachNucleic});
    od::logInfo(std::format("  {:<40} grow: {} segments, {} reach checks", name, grow->stepCount(),
                            grow->reachCheckCount()));
  }
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { sampler.sampleSegment(r, seg, a); };
  const od::ConfGenerator gen = [&](od::Rng& r, std::span<double> a) {
    sampler.sample(r, a);
    return !grow || grow->build(r, a, resample).ok;
  };
  std::uint64_t counter = 0;
  std::vector<od::ConfSample> pool;
  const auto t0 = std::chrono::steady_clock::now();
  const auto st = od::fillL2Pool(sys, gen, n.initL1, p, counter, pool);
  const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::vector<float> e;
  for (const auto& s : pool) e.push_back(s.energy);
  std::sort(e.begin(), e.end());
  const long long minimized = st.l1 - st.growFailed;
  od::logInfo(std::format("  {:<40} L1 {:>9} in {} batch(es)  -> L2 {:>7}  (rejected: grow {}, E >= Emax {}, "
                          "not closed {})  steps/sample {:.1f}  {:.2f} s",
                          name, st.l1, st.batches, st.added, st.growFailed, st.overEmax, st.notClosed,
                          minimized > 0 ? static_cast<double>(st.minSteps) / static_cast<double>(minimized) : 0.0, sec));
  od::logInfo(std::format("    L2 energies: min {:.3f}  median {:.3f}  max {:.3f}", e.front(), e[e.size() / 2], e.back()));
}

// DSSP residue counts, e.g. "H 12  E 10  . 5".
void printDssp(const od::Topology& t) {
  if (t.dsspSS.empty()) return;
  std::string line;
  for (char c : std::string("HGIPEBTS.")) {
    const auto n = std::count(t.dsspSS.begin(), t.dsspSS.end(), c);
    if (n > 0) line += std::format("  {} {}", c, n);
  }
  od::logInfo("    DSSP" + line);
}

int run(int argc, char** argv) {
  od::Options opt = od::parseCommandLine(argc, argv);
  if (opt.showHelp) {
    std::fputs(od::usageText().c_str(), stdout);
    return 0;
  }
  if (opt.showVersion) {
    od::logInfo(std::format("OmegaDocker {}", od::kVersion));
    return 0;
  }
  od::validateOptions(opt);

  od::CcdLibrary ccd;
  for (const auto& f : opt.ccd) ccd.addFile(f);

  auto sta = od::loadMolecules(opt.stamFile, od::Role::StaM, ccd);
  od::logInfo(std::format("StaM group from {}: {} molecule(s)", opt.stamFile, sta.size()));
  std::vector<od::Topology> staTopo;
  for (const auto& m : sta) {
    printMolecule(m);
    staTopo.push_back(analyze(m, opt));
  }
  {
    std::vector<const od::Molecule*> mols;
    std::vector<od::Topology*> topos;
    for (std::size_t i = 0; i < sta.size(); ++i) {
      mols.push_back(&sta[i]);
      topos.push_back(&staTopo[i]);
    }
    od::assignDssp(mols, topos);
    for (const auto& t : staTopo) printDssp(t);
  }
  std::vector<od::Molecule> mobs;
  std::vector<od::Topology> mobTopo;
  for (const auto& f : opt.mobmFiles) {
    auto mob = od::loadMolecules(f, od::Role::MobM, ccd);
    if (mob.size() != 1) od::fail("{}: --mobm file must contain exactly one molecule, found {}", f, mob.size());
    od::logInfo(std::format("MobM from {}:", f));
    printMolecule(mob.front());
    od::Topology t = analyze(mob.front(), opt);
    od::assignDssp({&mob.front()}, {&t});
    printDssp(t);
    mobs.push_back(std::move(mob.front()));
    mobTopo.push_back(std::move(t));
  }

  // Molecule configuration (spec 5.1).
  std::vector<od::MolInput> inputs;
  for (std::size_t i = 0; i < sta.size(); ++i) inputs.push_back({&sta[i], &staTopo[i]});
  for (std::size_t i = 0; i < mobs.size(); ++i) inputs.push_back({&mobs[i], &mobTopo[i]});
  if (opt.molConfig.empty()) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs);
    od::logInfo(std::format("Default molecule configuration written to {}\n"
                            "Review it (a copy may be edited) and run again with --molConfig <file>.",
                            path));
    return 0;
  }
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(opt.molConfig, inputs, topologyParams(opt), errors);
  if (!cfg) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs);
    std::string list;
    for (const auto& e : errors) list += "\n  " + e;
    od::fail("{}: {} problem(s) found:{}\nDefault molecule configuration written to {}", opt.molConfig, errors.size(),
             list, path);
  }
  od::logInfo(std::format("Molecule configuration {}:", opt.molConfig));
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    const od::MolSettings& s = cfg->molecules[i];
    checkOmegas(*inputs[i].mol, s, opt);
    const auto free = std::count(s.mode.begin(), s.mode.end(), od::TorsionMode::Free);
    const auto cisTrans = std::count(s.mode.begin(), s.mode.end(), od::TorsionMode::CisTrans);
    od::logInfo(std::format("  {:<40} rotatable {:>5}  cisTrans omega {:>4}  ringBreak {:>3}  chainBreak {:>3}",
                            inputs[i].mol->id, free, cisTrans, s.ringBreaks.size(), s.chainBreaks.size()));
    printTree(*inputs[i].mol, s);
  }

  {
    std::size_t staTorsions = 0, mobTorsions = 0;
    for (std::size_t i = 0; i < inputs.size(); ++i)
      (i < sta.size() ? staTorsions : mobTorsions) += cfg->molecules[i].tree.torsions.size();
    if (staTorsions > 0) od::requireConfSampleNums(opt, od::SampleSet::StaM);
    if (mobTorsions > 0) od::requireConfSampleNums(opt, od::SampleSet::MobM);
  }

  od::logInfo("9a energy:");
  {
    std::vector<od::EnergyPart> parts;
    for (std::size_t i = 0; i < sta.size(); ++i) parts.push_back({&sta[i], &cfg->molecules[i].tree});
    printInputEnergy("StaM group", parts, opt);
  }
  for (std::size_t i = 0; i < mobs.size(); ++i)
    printInputEnergy(mobs[i].id, {{&mobs[i], &cfg->molecules[sta.size() + i].tree}}, opt);

  const std::uint64_t seed = opt.seed ? *opt.seed : (static_cast<std::uint64_t>(std::random_device{}()) << 32) ^ std::random_device{}();
  od::logInfo(std::format("9a first refill (csgm1), seed {}:", seed));
  {
    std::vector<od::EnergyPart> parts;
    std::vector<const od::MolSettings*> settings;
    for (std::size_t i = 0; i < sta.size(); ++i) {
      parts.push_back({&sta[i], &cfg->molecules[i].tree});
      settings.push_back(&cfg->molecules[i]);
    }
    runFirst9a("StaM group", parts, settings, *cfg, opt, od::SampleSet::StaM, seed, 0);
  }
  for (std::size_t i = 0; i < mobs.size(); ++i) {
    const auto& s = cfg->molecules[sta.size() + i];
    runFirst9a(mobs[i].id, {{&mobs[i], &s.tree}}, {&s}, *cfg, opt, od::SampleSet::MobM, seed, 1 + i);
  }

  {
    const int level = od::hopfLevelForStep(opt.rotStep);
    od::logInfo(std::format("9b orientations: Hopf level {}, {} orientations, spacing {} deg (--rotStep {})", level,
                            72LL << (3 * level), od::hopfSpacingDeg(level), opt.rotStep));
  }
  od::logInfo("9b grids of the input StaM conformation:");
  {
    std::vector<od::EnergyPart> parts;
    for (std::size_t i = 0; i < sta.size(); ++i) parts.push_back({&sta[i], &cfg->molecules[i].tree});
    const od::EnergySystem sys(parts, energyParams(opt));
    std::vector<od::Vec3> x;
    sys.evaluate(sys.inputAngles(), nullptr, &x);
    const od::GridAtoms stam = od::gridAtoms(sys, parts, x);
    std::vector<od::Ad4Type> probe;
    for (std::size_t i = 0; i < mobs.size(); ++i) {
      const od::EnergySystem ms({{&mobs[i], &cfg->molecules[sta.size() + i].tree}}, energyParams(opt));
      for (od::Ad4Type t : ms.typing(0).type)
        if (t != od::Ad4Type::Merged && std::find(probe.begin(), probe.end(), t) == probe.end()) probe.push_back(t);
    }
    const od::GridGeometry geom = od::boxGrid(*opt.boxCenter, *opt.boxSize, opt.gridSpacing);
    const auto t0 = std::chrono::steady_clock::now();
    const od::EnergyGrids grids = od::buildEnergyGrids(geom, stam, probe, energyParams(opt));
    const auto t1 = std::chrono::steady_clock::now();
    const od::DistanceGrid dist = od::buildDistanceGrid(geom, stam, std::max(opt.maxNearDistance, opt.clashDistance) + 2 * opt.gridSpacing);
    const auto t2 = std::chrono::steady_clock::now();
    std::string types;
    for (od::Ad4Type t : grids.types) types += std::string(" ") + od::ad4TypeName(t);
    od::logInfo(std::format("  {} x {} x {} points, spacing {} Å, StaM atoms {}, type maps{} + elec + desolv: {:.2f} s; "
                            "distance grid: {:.2f} s",
                            geom.n[0], geom.n[1], geom.n[2], geom.spacing, stam.size(), types,
                            std::chrono::duration<double>(t1 - t0).count(),
                            std::chrono::duration<double>(t2 - t1).count()));
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const od::Error& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (const std::exception& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "ERROR: unexpected: %s\n", e.what());
    return 1;
  }
}
