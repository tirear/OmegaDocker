// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <cstdio>
#include <algorithm>
#include <format>
#include <numeric>

#include "ccd.hpp"
#include "cli.hpp"
#include "log.hpp"
#include "read_input.hpp"
#include "dssp.hpp"
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
  p.staFlexMargin = opt.staFlexMargin;
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

  auto sta = od::loadMolecules(opt.sta, od::Role::StaM, ccd);
  od::logInfo(std::format("StaM group from {}: {} molecule(s)", opt.sta, sta.size()));
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
  for (const auto& f : opt.mob) {
    auto mob = od::loadMolecules(f, od::Role::MobM, ccd);
    if (mob.size() != 1) od::fail("{}: --mob file must contain exactly one molecule, found {}", f, mob.size());
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
