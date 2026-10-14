// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "read_rdkit.hpp"

#include <cctype>
#include <format>
#include <fstream>
#include <map>
#include <sstream>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>

#include "log.hpp"

namespace od {
namespace {

RawStructure fromRdkit(const RDKit::ROMol& mol, const std::string& path, const std::string& recordName) {
  RawStructure s;
  s.sourceFile = path;
  RDKit::RWMol kek(mol);
  try {
    RDKit::MolOps::Kekulize(kek, false);
  } catch (const std::exception& e) {
    fail("{}{}: cannot kekulize: {}", path, recordName, e.what());
  }
  const bool hasConf = mol.getNumConformers() > 0;
  if (!hasConf) fail("{}{}: no coordinates", path, recordName);
  const RDKit::Conformer& conf = mol.getConformer();
  int implicitH = 0;
  for (const RDKit::Atom* at : mol.atoms()) {
    Atom a;
    a.atomicNum = at->getAtomicNum();
    a.formalCharge = at->getFormalCharge();
    a.aromatic = at->getIsAromatic();
    const auto& p = conf.getAtomPos(at->getIdx());
    a.xyz = {p.x, p.y, p.z};
    a.fileSerial = static_cast<int>(at->getIdx()) + 1;
    if (a.atomicNum == 0) fail("{}{}: atom #{} has no element", path, recordName, a.fileSerial);
    implicitH += static_cast<int>(at->getNumImplicitHs());
    s.atoms.push_back(std::move(a));
  }
  if (implicitH > 0)
    logWarn(std::format("{}{}: {} hydrogen(s) implied by valence are missing from the file and are not added",
                        path, recordName, implicitH));
  for (const RDKit::Bond* b : kek.bonds()) {
    Bond nb;
    nb.a = static_cast<int>(b->getBeginAtomIdx());
    nb.b = static_cast<int>(b->getEndAtomIdx());
    nb.aromatic = mol.getBondWithIdx(b->getIdx())->getIsAromatic();
    switch (b->getBondType()) {
      case RDKit::Bond::SINGLE: nb.order = BondOrder::Single; break;
      case RDKit::Bond::DOUBLE: nb.order = BondOrder::Double; break;
      case RDKit::Bond::TRIPLE: nb.order = BondOrder::Triple; break;
      case RDKit::Bond::DATIVE: nb.order = BondOrder::Metal; break;
      default:
        fail("{}{}: bond #{}-#{} has an unsupported bond type", path, recordName, nb.a + 1, nb.b + 1);
    }
    s.bonds.push_back(nb);
  }
  return s;
}

std::string trim(const std::string& x) {
  auto b = x.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return {};
  return x.substr(b, x.find_last_not_of(" \t\r\n") - b + 1);
}

// Fields RDKit's mol2 parser does not keep: atom names, substructures and charges.
struct Mol2Extra {
  struct AtomRec {
    std::string name;
    int substId = 0;
    std::string substName;
    bool hasCharge = false;
    double charge = 0.0;
  };
  std::vector<AtomRec> atoms;
  std::map<int, std::string> substChain;
  bool noCharges = false;
};

Mol2Extra parseMol2Extra(const std::string& block) {
  Mol2Extra x;
  std::istringstream in(block);
  std::string line, section;
  int moleculeLine = 0;
  while (std::getline(in, line)) {
    std::string t = trim(line);
    if (t.empty() || t[0] == '#') continue;
    if (t.rfind("@<TRIPOS>", 0) == 0) {
      section = t.substr(9);
      moleculeLine = 0;
      continue;
    }
    std::istringstream f(t);
    if (section == "MOLECULE") {
      if (++moleculeLine == 4 && t == "NO_CHARGES") x.noCharges = true;
    } else if (section == "ATOM") {
      Mol2Extra::AtomRec r;
      std::string id, xs, ys, zs, type, charge;
      f >> id >> r.name >> xs >> ys >> zs >> type;
      if (f >> r.substId >> r.substName) {
        if (f >> charge) {
          r.hasCharge = true;
          r.charge = std::stod(charge);
        }
      }
      x.atoms.push_back(std::move(r));
    } else if (section == "SUBSTRUCTURE") {
      int id;
      std::string name, root, type, dictType, chain;
      if (f >> id >> name >> root) {
        if (f >> type >> dictType >> chain && chain != "****") x.substChain[id] = chain;
      }
    }
  }
  return x;
}

// "ALA12" -> ("ALA", 12); falls back to the substructure id for the number.
void splitSubstName(const std::string& name, int substId, std::string& resName, int& seq) {
  std::size_t i = name.size();
  while (i > 0 && std::isdigit(static_cast<unsigned char>(name[i - 1]))) --i;
  if (i > 0 && name[i - 1] == '-') --i;
  if (i == name.size() || i == 0) {
    resName = name;
    seq = substId;
  } else {
    resName = name.substr(0, i);
    seq = std::stoi(name.substr(i));
  }
}

}  // namespace

std::vector<RawStructure> readSdf(const std::string& path) {
  auto stream = std::make_unique<std::ifstream>(path);
  if (!*stream) fail("Cannot open {}", path);
  RDKit::v2::FileParsers::MolFileParserParams params;
  params.sanitize = true;
  params.removeHs = false;
  RDKit::v2::FileParsers::ForwardSDMolSupplier sup(stream.release(), true, params);
  std::vector<RawStructure> out;
  for (int rec = 1; !sup.atEnd(); ++rec) {
    std::unique_ptr<RDKit::RWMol> mol;
    try {
      mol = sup.next();
    } catch (const std::exception& e) {
      fail("{} record {}: {}", path, rec, e.what());
    }
    if (!mol) {
      if (sup.atEnd()) break;
      fail("{} record {}: cannot be parsed or sanitized", path, rec);
    }
    out.push_back(fromRdkit(*mol, path, std::format(" record {}", rec)));
  }
  if (out.empty()) fail("{}: no molecules", path);
  return out;
}

std::vector<RawStructure> readMol2(const std::string& path) {
  std::ifstream in(path);
  if (!in) fail("Cannot open {}", path);
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string text = ss.str();
  std::vector<std::string> blocks;
  const std::string tag = "@<TRIPOS>MOLECULE";
  for (std::size_t p = text.find(tag); p != std::string::npos;) {
    std::size_t q = text.find(tag, p + tag.size());
    blocks.push_back(text.substr(p, q == std::string::npos ? std::string::npos : q - p));
    p = q;
  }
  if (blocks.empty()) fail("{}: no @<TRIPOS>MOLECULE record", path);

  std::vector<RawStructure> out;
  for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
    std::string recName = std::format(" record {}", bi + 1);
    RDKit::v2::FileParsers::Mol2ParserParams params;
    params.sanitize = true;
    params.removeHs = false;
    std::unique_ptr<RDKit::RWMol> mol;
    try {
      mol = RDKit::v2::FileParsers::MolFromMol2Block(blocks[bi], params);
    } catch (const std::exception& e) {
      fail("{}{}: {}", path, recName, e.what());
    }
    if (!mol) fail("{}{}: cannot be parsed or sanitized", path, recName);
    RawStructure s = fromRdkit(*mol, path, recName);
    Mol2Extra x = parseMol2Extra(blocks[bi]);
    if (x.atoms.size() != s.atoms.size())
      fail("{}{}: {} ATOM lines but {} atoms parsed", path, recName, x.atoms.size(), s.atoms.size());

    bool allCharges = !x.noCharges;
    bool allSubst = true;
    for (const auto& r : x.atoms) {
      allCharges = allCharges && r.hasCharge;
      allSubst = allSubst && !r.substName.empty();
    }
    s.hasPartialCharges = allCharges;
    if (!allCharges)
      logWarn(std::format("{}{}: partial charges missing; Gasteiger charges will be computed", path, recName));
    std::map<int, int> residueOfSubst;
    for (std::size_t i = 0; i < x.atoms.size(); ++i) {
      const auto& r = x.atoms[i];
      s.atoms[i].name = r.name;
      if (allCharges) s.atoms[i].partialCharge = r.charge;
      if (!allSubst) continue;
      auto [it, inserted] = residueOfSubst.try_emplace(r.substId, static_cast<int>(s.residues.size()));
      if (inserted) {
        ResidueInfo info;
        splitSubstName(r.substName, r.substId, info.name, info.seq);
        if (auto c = x.substChain.find(r.substId); c != x.substChain.end()) info.chain = c->second;
        s.residues.push_back(Residue{info, {}});
      }
      s.residues[it->second].atoms.push_back(static_cast<int>(i));
      s.atoms[i].residue = it->second;
    }
    s.hasResidueInfo = allSubst;
    out.push_back(std::move(s));
  }
  return out;
}

}  // namespace od
