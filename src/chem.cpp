// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "chem.hpp"

#include <cmath>
#include <format>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>

#include "log.hpp"

namespace od {

std::unique_ptr<RDKit::RWMol> toRdkit(const Molecule& m) {
  auto mol = std::make_unique<RDKit::RWMol>();
  auto* conf = new RDKit::Conformer(static_cast<unsigned>(m.atoms.size()));
  conf->set3D(true);
  for (std::size_t i = 0; i < m.atoms.size(); ++i) {
    const Atom& a = m.atoms[i];
    RDKit::Atom at(a.atomicNum);
    at.setFormalCharge(a.formalCharge);
    at.setNoImplicit(true);
    mol->addAtom(&at, true, false);
    conf->setAtomPos(static_cast<unsigned>(i), RDGeom::Point3D(a.xyz[0], a.xyz[1], a.xyz[2]));
  }
  mol->addConformer(conf, true);
  for (const Bond& b : m.bonds) {
    RDKit::Bond::BondType t = RDKit::Bond::SINGLE;
    switch (b.order) {
      case BondOrder::Single: t = RDKit::Bond::SINGLE; break;
      case BondOrder::Double: t = RDKit::Bond::DOUBLE; break;
      case BondOrder::Triple: t = RDKit::Bond::TRIPLE; break;
      case BondOrder::Metal:  t = RDKit::Bond::ZERO; break;
    }
    mol->addBond(static_cast<unsigned>(b.a), static_cast<unsigned>(b.b), t);
  }
  return mol;
}

void perceiveChemistry(Molecule& m) {
  auto mol = toRdkit(m);
  try {
    RDKit::MolOps::sanitizeMol(*mol);
  } catch (const std::exception& e) {
    fail("{}: chemistry check failed: {}", m.id, e.what());
  }
  for (std::size_t i = 0; i < m.atoms.size(); ++i)
    m.atoms[i].aromatic = mol->getAtomWithIdx(static_cast<unsigned>(i))->getIsAromatic();
  for (std::size_t i = 0; i < m.bonds.size(); ++i)
    m.bonds[i].aromatic = mol->getBondWithIdx(static_cast<unsigned>(i))->getIsAromatic();

  if (m.chargesFromFile) return;
  // Metal bonds are removed so that atoms without parameters stay isolated;
  // metals take their formal charge.
  for (int i = static_cast<int>(m.bonds.size()) - 1; i >= 0; --i)
    if (isMetalBond(m, m.bonds[i])) mol->removeBond(static_cast<unsigned>(m.bonds[i].a), static_cast<unsigned>(m.bonds[i].b));
  RDKit::computeGasteigerCharges(*mol, 12, false);
  int unparameterized = 0;
  for (std::size_t i = 0; i < m.atoms.size(); ++i) {
    if (isMetal(m.atoms[i].atomicNum)) {
      m.atoms[i].partialCharge = m.atoms[i].formalCharge;
      continue;
    }
    double q = 0;
    mol->getAtomWithIdx(static_cast<unsigned>(i))->getPropIfPresent(RDKit::common_properties::_GasteigerCharge, q);
    if (!std::isfinite(q)) {
      q = m.atoms[i].formalCharge;
      ++unparameterized;
    }
    m.atoms[i].partialCharge = q;
  }
  if (unparameterized > 0)
    logWarn(std::format("{}: {} atom(s) without Gasteiger parameters use their formal charge", m.id,
                        unparameterized));
}

}  // namespace od
