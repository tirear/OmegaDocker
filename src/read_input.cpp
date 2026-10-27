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

namespace {

// SDF / mol2 record: hydrogens replaced (spec 3.4); partial charges are then computed.
void protonateRecord(RawStructure& s, const std::string& path, std::size_t record, const ProtonationOptions& prot,
                     ProtonationReport& report, ProtonationPlan* planOut, const ProtonationPlan* planIn) {
  std::vector<std::string> errors;
  protonate(s, nullptr, prot, report, nullptr, nullptr, errors, planOut, planIn);
  if (!errors.empty()) {
    std::string msg = std::format("{} record {}: protonation failed:", path, record + 1);
    for (const auto& e : errors) msg += "\n  " + e;
    throw Error(msg);
  }
  if (s.hasPartialCharges)
    logWarn(std::format("{} record {}: partial charges of the file are not used after protonation; Gasteiger charges will be computed",
                        path, record + 1));
  s.hasPartialCharges = false;
}

}  // namespace

std::vector<RawStructure> readStructureFile(const std::string& path, const CcdLibrary& ccd, const ProtonationOptions* prot,
                                            ProtonationReport* report) {
  if (!std::filesystem::exists(path)) fail("File not found: {}", path);
  const std::string e = lowerExtension(path);
  if (e == ".cif" || e == ".mmcif") return {readMmcif(path, ccd, false, prot, report)};
  std::vector<RawStructure> v;
  if (e == ".sdf" || e == ".sd" || e == ".mol") v = readSdf(path);
  else if (e == ".mol2") v = readMol2(path);
  else fail("{}: unsupported input format '{}' (supported: mmCIF, SDF/MOL, mol2)", path, e);
  if (prot) {
    ProtonationReport scratch;
    for (std::size_t k = 0; k < v.size(); ++k) protonateRecord(v[k], path, k, *prot, report ? *report : scratch, nullptr, nullptr);
  }
  return v;
}

namespace {

void nameAndPerceive(std::vector<Molecule>& mols, const std::string& path, Role role) {
  const std::string fileName = std::filesystem::path(path).filename().string();

  std::map<std::string, int> keyCount;
  for (const Molecule& m : mols)
    if (m.hasResidueInfo) ++keyCount[chainKey(m)];
  std::vector<std::string> tails(mols.size());
  for (std::size_t k = 0; k < mols.size(); ++k) {
    const Molecule& m = mols[k];
    std::string& tail = tails[k];
    tail = m.hasResidueInfo ? chainKey(m) : std::string();
    if (!tail.empty()) {
      if (keyCount[tail] > 1) {
        const ResidueInfo& r = m.residues.front().info;
        tail += std::format(":{}{}{}", r.name, r.seq, r.iCode);
      }
    } else {
      tail = std::to_string(k + 1);
    }
  }
  // Fragments of one residue (e.g. missing atoms): the name of the first atom, then the order in the file.
  std::map<std::string, int> tailCount;
  for (const auto& t : tails) ++tailCount[t];
  for (std::size_t k = 0; k < mols.size(); ++k) {
    if (tailCount[tails[k]] < 2) continue;
    const Molecule& m = mols[k];
    const int a = m.residues.empty() || m.residues.front().atoms.empty() ? 0 : m.residues.front().atoms.front();
    tails[k] += ":" + m.atoms[a].name;
  }
  std::map<std::string, int> nameCount, ordinal;
  for (const auto& t : tails) ++nameCount[t];
  for (std::size_t k = 0; k < mols.size(); ++k)
    if (nameCount[tails[k]] > 1) tails[k] += std::format("#{}", ++ordinal[tails[k]]);
  std::set<std::string> used;
  for (std::size_t k = 0; k < mols.size(); ++k) {
    Molecule& m = mols[k];
    m.role = role;
    m.id = std::format("{}:{}:{}", roleTag(role), fileName, tails[k]);
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

std::vector<Molecule> loadMolecules(const std::string& path, Role role, const CcdLibrary& ccd, const ProtonationOptions* prot,
                                    ProtonationReport* report) {
  std::vector<Molecule> mols;
  for (const RawStructure& s : readStructureFile(path, ccd, prot, report))
    for (Molecule& m : splitMolecules(s)) mols.push_back(std::move(m));
  nameAndPerceive(mols, path, role);
  return mols;
}

PoseRecords loadPoseRecords(const std::string& path, Role role, const CcdLibrary& ccd, const ProtonationOptions* prot,
                            ProtonationReport* report) {
  if (!std::filesystem::exists(path)) fail("File not found: {}", path);
  const std::string e = lowerExtension(path);
  std::vector<RawStructure> recs;
  if (e == ".cif" || e == ".mmcif") {
    recs = readMmcifModels(path, ccd, prot, report);
  } else {
    recs = readStructureFile(path, ccd);
    if (prot) {
      ProtonationPlan plan;
      ProtonationReport scratch;
      for (std::size_t k = 0; k < recs.size(); ++k)
        protonateRecord(recs[k], path, k, *prot, k == 0 && report ? *report : scratch, k == 0 ? &plan : nullptr,
                        k == 0 ? nullptr : &plan);
    }
  }
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
