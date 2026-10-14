// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "ccd.hpp"
#include "log.hpp"
#include "read_input.hpp"

namespace {

const std::string kData = OD_TEST_DATA;

const od::CcdLibrary& ccd() {
  static od::CcdLibrary lib;
  return lib;
}

int netFormalCharge(const od::Molecule& m) {
  int q = 0;
  for (const auto& a : m.atoms) q += a.formalCharge;
  return q;
}

double sumPartialCharge(const od::Molecule& m) {
  double q = 0;
  for (const auto& a : m.atoms) q += a.partialCharge;
  return q;
}

int aromaticBondCount(const od::Molecule& m) {
  int n = 0;
  for (const auto& b : m.bonds) n += b.aromatic;
  return n;
}

const od::Atom& atomOf(const od::Molecule& m, const std::string& chain, int seq, const std::string& name) {
  for (const auto& r : m.residues)
    if (r.info.chain == chain && r.info.seq == seq)
      for (int i : r.atoms)
        if (m.atoms[i].name == name) return m.atoms[i];
  FAIL("atom not found: " << chain << ":" << seq << ":" << name);
  return m.atoms.front();
}

std::string errorOf(const std::string& path, od::Role role = od::Role::StaM) {
  try {
    od::loadMolecules(path, role, ccd());
  } catch (const od::Error& e) {
    return e.what();
  }
  return {};
}

}  // namespace

TEST_CASE("ccd: embedded subset") {
  CHECK(ccd().size() == 29);
  const od::CcdComp* his = ccd().find("HIS");
  REQUIRE(his);
  CHECK(his->find("NE2") >= 0);
  CHECK(ccd().find("PCA") == nullptr);
}

TEST_CASE("mmCIF: Trp-cage 1L2Y") {
  auto mols = od::loadMolecules(kData + "/1L2Y.cif", od::Role::StaM, ccd());
  REQUIRE(mols.size() == 1);
  const auto& m = mols[0];
  CHECK(m.id == "sta:1L2Y.cif:A");
  CHECK(m.atoms.size() == 304);
  CHECK(m.bonds.size() == 310);  // 303 tree bonds + 7 rings (TYR 1, TRP 2, PRO 4)
  CHECK(m.residues.size() == 20);
  CHECK(netFormalCharge(m) == 1);  // N-term +1, LYS +1, ARG +1, ASP -1, C-term -1
  CHECK(sumPartialCharge(m) == doctest::Approx(1.0).epsilon(1e-6));
  CHECK(aromaticBondCount(m) == 16);  // TYR 6 + TRP 10
  CHECK(atomOf(m, "A", 1, "N").formalCharge == 1);
  CHECK(atomOf(m, "A", 20, "OXT").formalCharge == -1);
}

TEST_CASE("mmCIF: HIS tautomers and charged HIS") {
  // 1D3Z HIS68 carries HD1; 2KJJ B:HIS5 carries HE2.
  auto ub = od::loadMolecules(kData + "/1D3Z.cif", od::Role::StaM, ccd());
  REQUIRE(ub.size() == 1);
  CHECK(netFormalCharge(ub[0]) == 0);
  CHECK(atomOf(ub[0], "A", 68, "ND1").formalCharge == 0);
  CHECK(atomOf(ub[0], "A", 68, "NE2").formalCharge == 0);
  CHECK(atomOf(ub[0], "A", 68, "CG").aromatic);

  auto ins = od::loadMolecules(kData + "/2KJJ.cif", od::Role::StaM, ccd());
  REQUIRE(ins.size() == 1);  // chains A and B joined by disulfides
  CHECK(ins[0].id == "sta:2KJJ.cif:A+B");
  CHECK(atomOf(ins[0], "B", 5, "NE2").formalCharge == 0);
  CHECK(atomOf(ins[0], "B", 5, "ND1").formalCharge == 0);

  // Add HE2 to HIS68 of 1D3Z: the ring becomes imidazolium with net charge +1.
  std::ifstream in(kData + "/1D3Z.cif");
  std::stringstream ss;
  ss << in.rdbuf();
  std::string text = ss.str();
  const std::string he1 = "ATOM 1095  H HE1  . HIS A 1 68 ? 43.796 -74.846 2.342   1.00 0.00  ? 68 HIS A HE1  1  \n";
  auto pos = text.find(he1);
  REQUIRE(pos != std::string::npos);
  text.insert(pos + he1.size(),
              "ATOM 9999  H HE2  . HIS A 1 68 ? 42.100 -75.200 0.000   1.00 0.00  ? 68 HIS A HE2  1  \n");
  auto tmp = std::filesystem::temp_directory_path() / "od_test_1D3Z_HIP.cif";
  std::ofstream(tmp) << text;
  auto hip = od::loadMolecules(tmp.string(), od::Role::StaM, ccd());
  std::filesystem::remove(tmp);
  REQUIRE(hip.size() == 1);
  CHECK(netFormalCharge(hip[0]) == 1);
  CHECK(atomOf(hip[0], "A", 68, "ND1").formalCharge + atomOf(hip[0], "A", 68, "NE2").formalCharge == 1);
  CHECK(sumPartialCharge(hip[0]) == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("mmCIF: disulfide bonds from struct_conn") {
  auto mols = od::loadMolecules(kData + "/2KJJ.cif", od::Role::StaM, ccd());
  REQUIRE(mols.size() == 1);
  int ss = 0;
  for (const auto& b : mols[0].bonds)
    ss += mols[0].atoms[b.a].atomicNum == 16 && mols[0].atoms[b.b].atomicNum == 16;
  CHECK(ss == 3);
}

TEST_CASE("mmCIF: reported problems") {
  CHECK(errorOf(kData + "/1PIT.cif").find("A:58:ALA:C: element C with bond-order sum 3") != std::string::npos);
  CHECK(errorOf(kData + "/2KOC.cif").find("A:1:G:O5'") != std::string::npos);
  CHECK(errorOf(kData + "/2LZX.cif").find("A:1:PCA: residue not found in CCD") != std::string::npos);
  CHECK(errorOf(kData + "/missing.cif").find("File not found") != std::string::npos);
  auto pdb = std::filesystem::temp_directory_path() / "od_test.pdb";
  std::ofstream(pdb) << "END\n";
  CHECK(errorOf(pdb.string()).find("unsupported input format") != std::string::npos);
  std::filesystem::remove(pdb);
}

TEST_CASE("SDF: biotin") {
  auto mols = od::loadMolecules(kData + "/BTN_ideal.sdf", od::Role::MobM, ccd());
  REQUIRE(mols.size() == 1);
  CHECK(mols[0].id == "mob:BTN_ideal.sdf:1");
  CHECK(mols[0].atoms.size() == 32);
  CHECK(mols[0].bonds.size() == 33);
  CHECK_FALSE(mols[0].hasResidueInfo);
  CHECK(netFormalCharge(mols[0]) == 0);
  CHECK(std::abs(sumPartialCharge(mols[0])) < 1e-6);
}

TEST_CASE("mol2: residues, aromaticity, multiple records") {
  auto bz = od::loadMolecules(kData + "/benzene.mol2", od::Role::MobM, ccd());
  REQUIRE(bz.size() == 1);
  CHECK(bz[0].id == "mob:benzene.mol2:1");  // no chain information
  CHECK(bz[0].hasResidueInfo);
  CHECK(bz[0].residues.size() == 1);
  CHECK(bz[0].residues[0].info.name == "BENZENE");
  CHECK(bz[0].atoms[0].name == "C1");
  CHECK(aromaticBondCount(bz[0]) == 6);
  CHECK_FALSE(bz[0].chargesFromFile);  // NO_CHARGES

  auto two = od::loadMolecules(kData + "/pyrazole_pyridine.mol2", od::Role::StaM, ccd());
  CHECK(two.size() == 2);
}
