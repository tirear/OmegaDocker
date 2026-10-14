// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <cstdio>
#include <format>
#include <numeric>

#include "ccd.hpp"
#include "cli.hpp"
#include "log.hpp"
#include "read_input.hpp"
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

od::Topology analyze(const od::Molecule& m, const od::Options& opt) {
  od::TopologyParams p;
  p.box = od::Box{*opt.boxCenter, *opt.boxSize};
  p.staFlexMargin = opt.staFlexMargin;
  od::Topology t = od::analyzeTopology(m, p);
  int rotatable = 0;
  for (const auto& tor : t.torsions) rotatable += tor.rotatable;
  od::logInfo(std::format("    backbone {} ({} units), torsions {}, rotatable by default {}",
                          od::backboneTypeName(t.backboneType), t.units.size(), t.torsions.size(), rotatable));

  const auto bad = od::nonPlanarOmegas(t, opt.omegaTolerance);
  if (!bad.empty()) {
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
  return t;
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
  for (const auto& m : sta) {
    printMolecule(m);
    analyze(m, opt);
  }
  for (const auto& f : opt.mob) {
    auto mob = od::loadMolecules(f, od::Role::MobM, ccd);
    if (mob.size() != 1) od::fail("{}: --mob file must contain exactly one molecule, found {}", f, mob.size());
    od::logInfo(std::format("MobM from {}:", f));
    printMolecule(mob.front());
    analyze(mob.front(), opt);
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
