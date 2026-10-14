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
  for (const auto& m : sta) printMolecule(m);
  for (const auto& f : opt.mob) {
    auto mob = od::loadMolecules(f, od::Role::MobM, ccd);
    if (mob.size() != 1) od::fail("{}: --mob file must contain exactly one molecule, found {}", f, mob.size());
    od::logInfo(std::format("MobM from {}:", f));
    printMolecule(mob.front());
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
