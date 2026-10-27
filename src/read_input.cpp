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

namespace {

void nameAndPerceive(std::vector<Molecule>& mols, const std::string& path, Role role) {
  const std::string fileName = std::filesystem::path(path).filename().string();

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
}

std::set<std::pair<int, int>> bondSet(const Molecule& m) {
  std::set<std::pair<int, int>> s;
  for (const Bond& b : m.bonds) s.insert(std::minmax(b.a, b.b));
  return s;
}

}  // namespace

std::vector<Molecule> loadMolecules(const std::string& path, Role role, const CcdLibrary& ccd) {
  std::vector<Molecule> mols;
  for (const RawStructure& s : readStructureFile(path, ccd))
    for (Molecule& m : splitMolecules(s)) mols.push_back(std::move(m));
  nameAndPerceive(mols, path, role);
  return mols;
}

PoseRecords loadPoseRecords(const std::string& path, Role role, const CcdLibrary& ccd) {
  if (!std::filesystem::exists(path)) fail("File not found: {}", path);
  const std::string e = lowerExtension(path);
  const std::vector<RawStructure> recs =
      e == ".cif" || e == ".mmcif" ? readMmcifModels(path, ccd) : readStructureFile(path, ccd);
  PoseRecords out;
  std::set<std::pair<int, int>> bonds;
  for (std::size_t k = 0; k < recs.size(); ++k) {
    std::vector<Molecule> mols = splitMolecules(recs[k]);
    if (mols.size() != 1) fail("{}: record {} contains {} molecules; each record must be one pose of the same molecule", path, k + 1, mols.size());
    Molecule& m = mols.front();
    std::vector<Vec3> x;
    for (const Atom& a : m.atoms) x.push_back(a.xyz);
    if (k == 0) {
      bonds = bondSet(m);
      nameAndPerceive(mols, path, role);
      out.mol = std::move(mols.front());
    } else {
      const Molecule& f = out.mol;
      if (m.atoms.size() != f.atoms.size())
        fail("{}: record {} has {} atoms, record 1 has {}", path, k + 1, m.atoms.size(), f.atoms.size());
      for (std::size_t i = 0; i < m.atoms.size(); ++i)
        if (m.atoms[i].atomicNum != f.atoms[i].atomicNum || m.atoms[i].name != f.atoms[i].name)
          fail("{}: atom {} of record {} ({} {}) differs from record 1 ({} {})", path, i + 1, k + 1, m.atoms[i].name,
               elementSymbol(m.atoms[i].atomicNum), f.atoms[i].name, elementSymbol(f.atoms[i].atomicNum));
      if (bondSet(m) != bonds) fail("{}: the bonds of record {} differ from record 1", path, k + 1);
    }
    out.xyz.push_back(std::move(x));
  }
  return out;
}

Molecule loadReference(const std::string& path, const CcdLibrary& ccd) {
  if (!std::filesystem::exists(path)) fail("File not found: {}", path);
  const std::string e = lowerExtension(path);
  const std::vector<RawStructure> recs =
      e == ".cif" || e == ".mmcif" ? std::vector<RawStructure>{readMmcif(path, ccd, true)} : readStructureFile(path, ccd);
  std::vector<Molecule> mols = splitMolecules(recs.front());
  if (mols.size() != 1) fail("{}: a reference structure must be one molecule, found {}", path, mols.size());
  if (recs.size() > 1) logWarn(std::format("{}: {} records; the first is the reference", path, recs.size()));
  mols.front().id = std::format("ref:{}", std::filesystem::path(path).filename().string());
  return std::move(mols.front());
}

}  // namespace od
