// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "read_mmcif.hpp"

#include <cmath>
#include <format>
#include <map>
#include <set>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#include "bond_orders.hpp"
#include "log.hpp"

namespace od {
namespace {

constexpr double kPeptideBondMax = 2.0;  // Å, C(i)–N(i+1)
constexpr double kPhosphoBondMax = 2.2;  // Å, O3'(i)–P(i+1)
constexpr std::size_t kMaxListedErrors = 50;

double distance(const Vec3& a, const Vec3& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}

std::string icodeOf(char c) { return (c == ' ' || c == '\0' || c == '?' || c == '.') ? "" : std::string(1, c); }

// Terminal atoms that CCD definitions of linking residues do not list:
// peptide N-terminal H/H1/H2/H3 and nucleotide 5'-terminal HO5'.
std::string terminalPartner(const CcdComp& c, const std::string& atom) {
  bool peptide = c.find("N") >= 0 && c.find("CA") >= 0 && c.find("C") >= 0;
  bool nucleotide = c.find("O5'") >= 0 && c.find("C5'") >= 0;
  if (peptide && (atom == "H" || atom == "H1" || atom == "H2" || atom == "H3")) return "N";
  if (nucleotide && atom == "HO5'") return "O5'";
  return {};
}

struct ResidueRec {
  const CcdComp* comp = nullptr;
  std::map<std::string, int> atomByName;  // atom name -> structure atom index
};

}  // namespace

RawStructure readMmcif(const std::string& path, const CcdLibrary& ccd) {
  gemmi::Structure st;
  try {
    gemmi::cif::Document doc = gemmi::cif::read_file(path);
    if (doc.blocks.empty()) fail("{}: no data block", path);
    st = gemmi::make_structure_from_block(doc.blocks.at(0));
  } catch (const Error&) {
    throw;
  } catch (const std::exception& e) {
    fail("Cannot read mmCIF file {}: {}", path, e.what());
  }
  if (st.models.empty()) fail("{}: no atoms", path);
  if (st.models.size() > 1) logWarn(std::format("{}: {} models found; only the first is used", path, st.models.size()));

  RawStructure s;
  s.sourceFile = path;
  s.hasResidueInfo = true;
  std::vector<std::string> errors;
  std::vector<bool> aromaticCandidate;
  std::vector<int> ccdCharge;
  std::set<std::pair<int, int>> bonded;
  std::map<std::string, int> atomByAddress;  // "chain|seq|icode|atom"
  std::vector<ResidueRec> recs;
  int altlocSkipped = 0, duplicateResidues = 0;

  auto addBond = [&](int a, int b, BondOrder o, bool aromatic) {
    auto key = std::minmax(a, b);
    if (!bonded.insert(key).second) return;
    s.bonds.push_back(Bond{a, b, o, false});
    aromaticCandidate.push_back(aromatic);
  };

  const gemmi::Model& model = st.models[0];
  for (const gemmi::Chain& chain : model.chains) {
    const gemmi::Residue* prev = nullptr;
    for (const gemmi::Residue& res : chain.residues) {
      if (prev && prev->seqid == res.seqid) {  // microheterogeneity: keep the first residue
        ++duplicateResidues;
        continue;
      }
      prev = &res;
      ResidueInfo info{chain.name, res.seqid.num.has_value() ? res.seqid.num.value : 0,
                       icodeOf(res.seqid.icode), res.name};
      const CcdComp* comp = ccd.find(res.name);
      if (!comp) {
        errors.push_back(std::format("{}:{}{}:{}: residue not found in CCD (use --ccd)", info.chain,
                                     info.seq, info.iCode, info.name));
        continue;
      }
      int ri = static_cast<int>(s.residues.size());
      s.residues.push_back(Residue{info, {}});
      ResidueRec rec;
      rec.comp = comp;
      char keptAltloc = '\0';
      for (const gemmi::Atom& at : res.atoms) {
        if (at.altloc != '\0') {
          if (keptAltloc == '\0') keptAltloc = at.altloc;
          if (at.altloc != keptAltloc) {
            ++altlocSkipped;
            continue;
          }
        }
        if (rec.atomByName.contains(at.name)) {
          errors.push_back(std::format("{}:{}{}:{}:{}: duplicate atom name", info.chain, info.seq,
                                       info.iCode, info.name, at.name));
          continue;
        }
        Atom a;
        a.name = at.name;
        a.atomicNum = atomicNumber(at.element.name());
        a.xyz = {at.pos.x, at.pos.y, at.pos.z};
        a.residue = ri;
        a.fileSerial = static_cast<int>(s.atoms.size()) + 1;
        int idx = static_cast<int>(s.atoms.size());
        int ci = comp->find(at.name);
        int charge = 0;
        if (ci >= 0) {
          if (comp->atoms[ci].atomicNum != a.atomicNum)
            errors.push_back(std::format("{}:{}{}:{}:{}: element {} differs from CCD ({})", info.chain,
                                         info.seq, info.iCode, info.name, at.name,
                                         elementSymbol(a.atomicNum), elementSymbol(comp->atoms[ci].atomicNum)));
          charge = comp->atoms[ci].charge;
        } else if (terminalPartner(*comp, at.name).empty()) {
          errors.push_back(std::format("{}:{}{}:{}:{}: atom not defined in CCD component {}", info.chain,
                                       info.seq, info.iCode, info.name, at.name, comp->id));
        }
        s.atoms.push_back(std::move(a));
        ccdCharge.push_back(charge);
        s.residues[ri].atoms.push_back(idx);
        rec.atomByName[at.name] = idx;
        atomByAddress[std::format("{}|{}|{}|{}", info.chain, info.seq, info.iCode, at.name)] = idx;
      }
      // Intra-residue bonds from the CCD definition.
      for (const CcdComp::Bond& cb : comp->bonds) {
        auto ia = rec.atomByName.find(comp->atoms[cb.a].name);
        auto ib = rec.atomByName.find(comp->atoms[cb.b].name);
        if (ia != rec.atomByName.end() && ib != rec.atomByName.end())
          addBond(ia->second, ib->second, static_cast<BondOrder>(cb.order), cb.aromatic);
      }
      for (const auto& [name, idx] : rec.atomByName) {
        std::string partner = terminalPartner(*comp, name);
        if (partner.empty() || comp->find(name) >= 0) continue;
        if (auto ip = rec.atomByName.find(partner); ip != rec.atomByName.end())
          addBond(idx, ip->second, BondOrder::Single, false);
      }
      recs.push_back(std::move(rec));
    }
  }

  // Polymer links between consecutive residues of a chain.
  for (std::size_t r = 0; r + 1 < s.residues.size(); ++r) {
    if (s.residues[r].info.chain != s.residues[r + 1].info.chain) continue;
    const auto& a = recs[r].atomByName;
    const auto& b = recs[r + 1].atomByName;
    auto link = [&](const char* an, const char* bn, double maxDist) {
      auto ia = a.find(an), ib = b.find(bn);
      if (ia == a.end() || ib == b.end()) return;
      if (distance(s.atoms[ia->second].xyz, s.atoms[ib->second].xyz) <= maxDist)
        addBond(ia->second, ib->second, BondOrder::Single, false);
    };
    link("C", "N", kPeptideBondMax);
    link("O3'", "P", kPhosphoBondMax);
  }

  // _struct_conn: covalent, disulfide and metal coordination bonds.
  for (const gemmi::Connection& c : st.connections) {
    if (c.type != gemmi::Connection::Covale && c.type != gemmi::Connection::Disulf &&
        c.type != gemmi::Connection::MetalC)
      continue;
    auto find = [&](const gemmi::AtomAddress& ad) {
      std::string key = std::format("{}|{}|{}|{}", ad.chain_name,
                                    ad.res_id.seqid.num.has_value() ? ad.res_id.seqid.num.value : 0,
                                    icodeOf(ad.res_id.seqid.icode), ad.atom_name);
      auto it = atomByAddress.find(key);
      return it == atomByAddress.end() ? -1 : it->second;
    };
    int ia = find(c.partner1), ib = find(c.partner2);
    if (ia < 0 || ib < 0) {
      logWarn(std::format("{}: struct_conn {} refers to an atom not read; ignored", path, c.name));
      continue;
    }
    if (ia == ib) continue;
    bool metal = c.type == gemmi::Connection::MetalC || isMetal(s.atoms[ia].atomicNum) ||
                 isMetal(s.atoms[ib].atomicNum);
    addBond(ia, ib, metal ? BondOrder::Metal : BondOrder::Single, false);
  }

  if (errors.empty()) assignChargesAndKekule(s, aromaticCandidate, ccdCharge, errors);
  if (!errors.empty()) {
    std::string msg = std::format("{}: {} problem(s) found:", path, errors.size());
    for (std::size_t i = 0; i < errors.size() && i < kMaxListedErrors; ++i) msg += "\n  " + errors[i];
    if (errors.size() > kMaxListedErrors) msg += std::format("\n  ... ({} more)", errors.size() - kMaxListedErrors);
    throw Error(msg);
  }

  // A thiolate close to another sulfur usually means a disulfide missing from _struct_conn.
  for (std::size_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].atomicNum != 16 || s.atoms[i].formalCharge != -1) continue;
    for (std::size_t j = 0; j < s.atoms.size(); ++j)
      if (j != i && s.atoms[j].atomicNum == 16 && distance(s.atoms[i].xyz, s.atoms[j].xyz) < 2.5) {
        logWarn(std::format("{}: {} is a thiolate 2.5 Å or closer to {}; disulfide missing from _struct_conn?",
                            path, describeAtom(s, static_cast<int>(i)), describeAtom(s, static_cast<int>(j))));
        break;
      }
  }
  if (altlocSkipped > 0)
    logWarn(std::format("{}: {} alternate-location atom(s) skipped (first altloc of each residue kept)", path,
                        altlocSkipped));
  if (duplicateResidues > 0)
    logWarn(std::format("{}: {} residue(s) with a repeated sequence number skipped", path, duplicateResidues));
  return s;
}

}  // namespace od
