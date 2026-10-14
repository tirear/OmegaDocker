// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "read_input.hpp"

#include <algorithm>
#include <filesystem>
#include <format>
#include <map>
#include <set>

#include "chem.hpp"
#include "log.hpp"
#include "read_mmcif.hpp"
#include "read_rdkit.hpp"

namespace od {
namespace {

std::string lowerExtension(const std::string& path) {
  std::string e = std::filesystem::path(path).extension().string();
  std::transform(e.begin(), e.end(), e.begin(), [](unsigned char c) { return std::tolower(c); });
  return e;
}

std::string chainKey(const Molecule& m) {
  std::vector<std::string> chains;
  for (const Residue& r : m.residues)
    if (std::find(chains.begin(), chains.end(), r.info.chain) == chains.end()) chains.push_back(r.info.chain);
  std::string k;
  for (const auto& c : chains) k += (k.empty() ? "" : "+") + c;
  return k;
}

}  // namespace

std::vector<RawStructure> readStructureFile(const std::string& path, const CcdLibrary& ccd) {
  if (!std::filesystem::exists(path)) fail("File not found: {}", path);
  const std::string e = lowerExtension(path);
  if (e == ".cif" || e == ".mmcif") return {readMmcif(path, ccd)};
  if (e == ".sdf" || e == ".sd" || e == ".mol") return readSdf(path);
  if (e == ".mol2") return readMol2(path);
  fail("{}: unsupported input format '{}' (supported: mmCIF, SDF/MOL, mol2)", path, e);
}

std::vector<Molecule> loadMolecules(const std::string& path, Role role, const CcdLibrary& ccd) {
  const std::string fileName = std::filesystem::path(path).filename().string();
  std::vector<Molecule> mols;
  for (const RawStructure& s : readStructureFile(path, ccd))
    for (Molecule& m : splitMolecules(s)) mols.push_back(std::move(m));

  std::map<std::string, int> keyCount;
  for (const Molecule& m : mols)
    if (m.hasResidueInfo) ++keyCount[chainKey(m)];
  std::set<std::string> used;
  for (std::size_t k = 0; k < mols.size(); ++k) {
    Molecule& m = mols[k];
    m.role = role;
    std::string tail = m.hasResidueInfo ? chainKey(m) : std::string();
    if (!tail.empty()) {
      if (keyCount[tail] > 1) {
        const ResidueInfo& r = m.residues.front().info;
        tail += std::format(":{}{}{}", r.name, r.seq, r.iCode);
      }
    } else {
      tail = std::to_string(k + 1);
    }
    m.id = std::format("{}:{}:{}", roleTag(role), fileName, tail);
    if (!used.insert(m.id).second) fail("{}: molecule id {} is not unique", path, m.id);
    perceiveChemistry(m);
  }
  return mols;
}

}  // namespace od
