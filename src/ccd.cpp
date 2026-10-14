// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "ccd.hpp"

#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>

#include "log.hpp"
#include "molecule.hpp"

namespace od {
namespace embedded {
extern const unsigned char kCcdData[];
extern const std::size_t kCcdSize;
}  // namespace embedded

namespace {

int bondOrderOf(gemmi::BondType t) {
  switch (t) {
    case gemmi::BondType::Single: return 1;
    case gemmi::BondType::Double: return 2;
    case gemmi::BondType::Triple: return 3;
    default: return 0;
  }
}

}  // namespace

CcdLibrary::CcdLibrary() {
  auto doc = gemmi::cif::read_memory(reinterpret_cast<const char*>(embedded::kCcdData),
                                     embedded::kCcdSize, "embedded CCD");
  addDocument(&doc, "embedded CCD", false);
}

void CcdLibrary::addFile(const std::string& path) {
  gemmi::cif::Document doc;
  try {
    doc = gemmi::cif::read_file(path);
  } catch (const std::exception& e) {
    fail("Cannot read CCD file {}: {}", path, e.what());
  }
  addDocument(&doc, path, true);
}

void CcdLibrary::addDocument(const void* docPtr, const std::string& origin, bool fromUser) {
  const auto& doc = *static_cast<const gemmi::cif::Document*>(docPtr);
  for (const gemmi::cif::Block& block : doc.blocks) {
    if (!const_cast<gemmi::cif::Block&>(block).find_values("_chem_comp_atom.atom_id")) continue;
    gemmi::ChemComp cc = gemmi::make_chemcomp_from_block(block);
    CcdComp c;
    c.id = cc.name;
    c.type = cc.type_or_group;
    for (const auto& a : cc.atoms) {
      CcdComp::Atom ca;
      ca.name = a.id;
      ca.atomicNum = atomicNumber(a.el.name());
      ca.charge = static_cast<int>(a.charge);
      if (ca.atomicNum == 0) fail("{}: component {}: unknown element for atom {}", origin, c.id, a.id);
      c.atomIndex.emplace(ca.name, static_cast<int>(c.atoms.size()));
      c.atoms.push_back(std::move(ca));
    }
    for (const auto& b : cc.rt.bonds) {
      CcdComp::Bond cb;
      cb.a = c.find(b.id1.atom);
      cb.b = c.find(b.id2.atom);
      cb.order = bondOrderOf(b.type);
      cb.aromatic = b.aromatic;
      if (cb.a < 0 || cb.b < 0)
        fail("{}: component {}: bond {} refers to an undefined atom", origin, c.id, b.str());
      if (cb.order == 0)
        fail("{}: component {}: bond {} has an unsupported value_order", origin, c.id, b.str());
      c.bonds.push_back(cb);
    }
    if (fromUser && comps_.contains(c.id))
      logWarn(std::format("{}: component {} replaces the embedded definition", origin, c.id));
    comps_[c.id] = std::move(c);
  }
}

const CcdComp* CcdLibrary::find(const std::string& id) const {
  auto it = comps_.find(id);
  return it == comps_.end() ? nullptr : &it->second;
}

}  // namespace od
