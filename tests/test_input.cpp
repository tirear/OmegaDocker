// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <limits>
#include <map>
#include <numeric>
#include <fstream>
#include <numbers>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <unordered_map>

#include <omp.h>
#include <unistd.h>

#include "ad4.hpp"
#include "ad4_pair.hpp"
#include "ccd.hpp"
#include "conf9a.hpp"
#include "csgm1.hpp"
#include "dock9b.hpp"
#include "dock9c.hpp"
#include "log.hpp"
#include "minimize.hpp"
#include "dssp.hpp"
#include "symrmsd.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "grow.hpp"
#include "orient.hpp"
#include "feparams.hpp"
#include "score9d.hpp"
#include "shaking.hpp"
#include "write_struct.hpp"
#include "molconfig.hpp"
#include "rama.hpp"
#include "read_input.hpp"
#include "topology.hpp"
#include "torsion_tree.hpp"

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
  CHECK(m.id == "stam:1L2Y.cif:A");
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
  CHECK(ins[0].id == "stam:2KJJ.cif:A+B");
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
  CHECK(mols[0].id == "mobm:BTN_ideal.sdf:1");
  CHECK(mols[0].atoms.size() == 32);
  CHECK(mols[0].bonds.size() == 33);
  CHECK_FALSE(mols[0].hasResidueInfo);
  CHECK(netFormalCharge(mols[0]) == 0);
  CHECK(std::abs(sumPartialCharge(mols[0])) < 1e-6);
}

TEST_CASE("mol2: residues, aromaticity, multiple records") {
  auto bz = od::loadMolecules(kData + "/benzene.mol2", od::Role::MobM, ccd());
  REQUIRE(bz.size() == 1);
  CHECK(bz[0].id == "mobm:benzene.mol2:1");  // no chain information
  CHECK(bz[0].hasResidueInfo);
  CHECK(bz[0].residues.size() == 1);
  CHECK(bz[0].residues[0].info.name == "BENZENE");
  CHECK(bz[0].atoms[0].name == "C1");
  CHECK(aromaticBondCount(bz[0]) == 6);
  CHECK_FALSE(bz[0].chargesFromFile);  // NO_CHARGES

  auto two = od::loadMolecules(kData + "/pyrazole_pyridine.mol2", od::Role::StaM, ccd());
  CHECK(two.size() == 2);
}


namespace {

int countAngle(const od::Topology& t, od::BackboneAngle a) {
  int n = 0;
  for (const auto& tor : t.torsions) n += tor.backbone == a;
  return n;
}

int countRotatable(const od::Topology& t) {
  int n = 0;
  for (const auto& tor : t.torsions) n += tor.rotatable;
  return n;
}

const od::Torsion& torsionAt(const od::Molecule& m, const od::Topology& t, int seq, const std::string& b,
                             const std::string& c) {
  for (const auto& tor : t.torsions) {
    const auto& x = m.atoms[tor.atoms[1]];
    const auto& y = m.atoms[tor.atoms[2]];
    if (x.residue >= 0 && m.residues[x.residue].info.seq == seq && x.name == b && y.name == c) return tor;
  }
  FAIL("torsion not found: " << seq << ":" << b << "-" << c);
  return t.torsions.front();
}

}  // namespace

TEST_CASE("topology: dihedral sign convention") {
  CHECK(od::dihedralDeg({1, 0, 0}, {0, 0, 0}, {0, 0, 1}, {0, 1, 1}) == doctest::Approx(90));
  CHECK(od::dihedralDeg({1, 0, 0}, {0, 0, 0}, {0, 0, 1}, {-1, 0, 1}) == doctest::Approx(180));
}

TEST_CASE("topology: protein backbone and side chains (1L2Y)") {
  auto mob = od::loadMolecules(kData + "/1L2Y.cif", od::Role::MobM, ccd());
  const auto& m = mob[0];
  auto t = od::analyzeTopology(m, {});
  CHECK(t.backboneType == od::BackboneType::Protein);
  CHECK(t.units.size() == 20);
  CHECK(countAngle(t, od::BackboneAngle::Phi) == 19);
  CHECK(countAngle(t, od::BackboneAngle::Psi) == 19);
  CHECK(countAngle(t, od::BackboneAngle::Omega) == 19);

  const auto& phi3 = torsionAt(m, t, 3, "N", "CA");
  CHECK(phi3.backbone == od::BackboneAngle::Phi);
  CHECK(phi3.angle == doctest::Approx(-66.5).epsilon(0.002));
  CHECK(phi3.rotatable);
  const auto& phiPro = torsionAt(m, t, 12, "N", "CA");  // proline phi lies in the ring
  CHECK(phiPro.category == od::TorsionCategory::Ring);
  CHECK_FALSE(phiPro.rotatable);
  const auto& omega = torsionAt(m, t, 5, "C", "N");
  CHECK(omega.backbone == od::BackboneAngle::Omega);
  CHECK(omega.category == od::TorsionCategory::Amide);
  CHECK_FALSE(omega.rotatable);
  CHECK(torsionAt(m, t, 16, "NE", "CZ").category == od::TorsionCategory::Resonance);
  CHECK(torsionAt(m, t, 9, "CB", "CG").category == od::TorsionCategory::Symmetric);  // carboxylate
  CHECK(torsionAt(m, t, 3, "CG", "CD1").category == od::TorsionCategory::Aromatic);
  CHECK(torsionAt(m, t, 8, "CD", "CE").rotatable);

  CHECK(od::nonPlanarOmegas(t, 30).empty());
  CHECK_FALSE(od::nonPlanarOmegas(t, 5).empty());

  // StaM: backbone fixed; side chains rotatable only near the box.
  auto sta = od::loadMolecules(kData + "/1L2Y.cif", od::Role::StaM, ccd());
  auto all = od::analyzeTopology(sta[0], {});
  CHECK_FALSE(torsionAt(sta[0], all, 3, "N", "CA").rotatable);
  CHECK(torsionAt(sta[0], all, 8, "CD", "CE").rotatable);
  od::TopologyParams p;
  p.box = od::Box{atomOf(sta[0], "A", 6, "CA").xyz, {1, 1, 1}};
  p.stamFlexMargin = 0;
  auto near = od::analyzeTopology(sta[0], p);
  CHECK(countRotatable(near) == 2);  // TRP6 chi1, chi2
  CHECK(torsionAt(sta[0], near, 6, "CA", "CB").rotatable);
  CHECK_FALSE(torsionAt(sta[0], near, 8, "CD", "CE").rotatable);
}

TEST_CASE("topology: nucleic-acid backbone (1D18 DNA, 1RNG RNA)") {
  auto dna = od::loadMolecules(kData + "/1D18.cif", od::Role::MobM, ccd());
  REQUIRE(dna.size() == 2);
  auto t = od::analyzeTopology(dna[0], {});
  CHECK(t.backboneType == od::BackboneType::Nucleic);
  CHECK(t.units.size() == 8);
  CHECK(countAngle(t, od::BackboneAngle::Alpha) == 7);  // the 5' residue has no phosphate
  CHECK(countAngle(t, od::BackboneAngle::Beta) == 7);
  CHECK(countAngle(t, od::BackboneAngle::Gamma) == 8);
  CHECK(countAngle(t, od::BackboneAngle::Delta) == 8);
  CHECK(countAngle(t, od::BackboneAngle::Epsilon) == 7);
  CHECK(countAngle(t, od::BackboneAngle::Zeta) == 7);
  const auto& delta = torsionAt(dna[0], t, 2, "C4'", "C3'");
  CHECK(delta.category == od::TorsionCategory::Ring);
  CHECK_FALSE(delta.rotatable);
  CHECK(torsionAt(dna[0], t, 2, "P", "O5'").angle == doctest::Approx(-64.8).epsilon(0.002));
  CHECK(countRotatable(t) == 7 + 7 + 8 + 7 + 7 + 8);  // alpha..zeta except delta, plus chi

  auto rna = od::loadMolecules(kData + "/1RNG.cif", od::Role::MobM, ccd());
  REQUIRE(rna.size() == 1);
  auto r = od::analyzeTopology(rna[0], {});
  CHECK(r.backboneType == od::BackboneType::Nucleic);
  CHECK(r.units.size() == 12);
}

TEST_CASE("topology: small molecule (biotin)") {
  auto mols = od::loadMolecules(kData + "/BTN_ideal.sdf", od::Role::MobM, ccd());
  auto t = od::analyzeTopology(mols[0], {});
  CHECK(t.backboneType == od::BackboneType::None);
  CHECK(t.torsions.size() == 14);
  CHECK(countRotatable(t) == 5);
}


namespace {

std::string dsspString(const std::string& file) {
  auto mols = od::loadMolecules(kData + "/" + file, od::Role::StaM, ccd());
  std::vector<od::Topology> topo;
  for (const auto& m : mols) topo.push_back(od::analyzeTopology(m, {}));
  std::vector<const od::Molecule*> mp;
  std::vector<od::Topology*> tp;
  for (std::size_t i = 0; i < mols.size(); ++i) {
    mp.push_back(&mols[i]);
    tp.push_back(&topo[i]);
  }
  od::assignDssp(mp, tp);
  std::string s;
  for (const auto& t : topo)
    for (char c : t.dsspSS)
      if (c != ' ') s += c;
  return s;
}

}  // namespace

// Reference strings from mkdssp 4.6.1 on the same files.
TEST_CASE("dssp: agreement with DSSP 4") {
  CHECK(dsspString("1L2Y.cif") == ".HHHHHHHTTGGGGT.PPP.");
  CHECK(dsspString("1D3Z.cif") ==
        ".EEEEE.TTS.EEEEE..TT.BHHHHHHHHHHHH...GGGEEEEETTEEPPTTSBTGGGT..TT.EEEEEE..S..");
  CHECK(dsspString("2KJJ.cif") == ".HIIIIISS.B.HHHHHTTB....B...SHHHHHHHHHHTGGG.B......");
}

TEST_CASE("dssp: molecules without a protein backbone") {
  auto mols = od::loadMolecules(kData + "/BTN_ideal.sdf", od::Role::MobM, ccd());
  auto t = od::analyzeTopology(mols[0], {});
  od::assignDssp({&mols[0]}, {&t});
  CHECK(t.dsspSS.empty());
  auto dna = od::loadMolecules(kData + "/1D18.cif", od::Role::MobM, ccd());
  auto d = od::analyzeTopology(dna[0], {});
  od::assignDssp({&dna[0]}, {&d});
  CHECK(d.dsspSS.empty());
}

TEST_CASE("topology: crosslinks are broken without disconnecting the molecule (2KJJ)") {
  auto mols = od::loadMolecules(kData + "/2KJJ.cif", od::Role::MobM, ccd());
  const auto& m = mols[0];
  auto t = od::analyzeTopology(m, {});
  // Three disulfides; one joins chains A and B and is kept.
  REQUIRE(t.crosslinkBreaks.size() == 2);
  CHECK(t.brokenBonds == t.crosslinkBreaks);
  for (int bi : t.crosslinkBreaks) {
    CHECK(m.atoms[m.bonds[bi].a].atomicNum == 16);
    CHECK(m.atoms[m.bonds[bi].b].atomicNum == 16);
  }
  int brokenTorsions = 0;
  for (const auto& tor : t.torsions) {
    brokenTorsions += tor.broken;
    if (tor.broken) CHECK_FALSE(tor.rotatable);
  }
  CHECK(brokenTorsions == 2);
  auto ringPhi = [](const od::Topology& x) {
    int n = 0;
    for (const auto& tor : x.torsions)
      n += tor.backbone == od::BackboneAngle::Phi && tor.category == od::TorsionCategory::Ring;
    return n;
  };
  od::TopologyParams p;
  p.brokenBonds = std::vector<int>{};
  auto whole = od::analyzeTopology(m, p);
  CHECK(whole.crosslinkBreaks == t.crosslinkBreaks);
  CHECK(ringPhi(t) == 1);  // the proline
  CHECK(ringPhi(whole) > ringPhi(t));
}

namespace {

struct ConfigFixture {
  std::vector<od::Molecule> mols;
  std::vector<od::Topology> topos;
  std::vector<od::MolInput> inputs;
  std::string path;

  ConfigFixture() {
    for (const char* f : {"1L2Y.cif", "BTN_ideal.sdf", "1RNG.cif"}) {
      mols.push_back(od::loadMolecules(kData + "/" + f, od::Role::MobM, ccd()).front());
      topos.push_back(od::analyzeTopology(mols.back(), {}));
      od::assignDssp({&mols.back()}, {&topos.back()});
    }
    for (std::size_t i = 0; i < mols.size(); ++i) inputs.push_back({&mols[i], &topos[i]});
    path = (std::filesystem::temp_directory_path() / std::format("od_test_molconfig_{}.toml", ::getpid())).string();
  }
  ~ConfigFixture() { std::filesystem::remove(path); }

  std::vector<std::string> errorsOf(const std::string& text) {
    std::ofstream(path) << text;
    std::vector<std::string> errors;
    auto cfg = od::readMolConfig(path, inputs, {}, errors);
    CHECK(cfg.has_value() == errors.empty());
    return errors;
  }
};

// Replaces the first occurrence of from after the first occurrence of anchor.
std::string replaceAfter(std::string s, const std::string& anchor, const std::string& from, const std::string& to) {
  const auto a = s.find(anchor);
  REQUIRE(a != std::string::npos);
  const auto i = s.find(from, a);
  REQUIRE(i != std::string::npos);
  return s.replace(i, from.size(), to);
}

bool hasError(const std::vector<std::string>& errors, const std::string& part) {
  for (const auto& e : errors)
    if (e.find(part) != std::string::npos) return true;
  return false;
}

}  // namespace

TEST_CASE("molconfig: default file is accepted and keeps the defaults") {
  ConfigFixture f;
  const std::string text = od::defaultMolConfigText(f.inputs);
  std::ofstream(f.path) << text;
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  REQUIRE(cfg);
  CHECK(errors.empty());
  REQUIRE(cfg->molecules.size() == 3);
  for (std::size_t i = 0; i < 3; ++i) {
    const auto& s = cfg->molecules[i];
    REQUIRE(s.mode.size() == f.topos[i].torsions.size());
    for (std::size_t k = 0; k < s.mode.size(); ++k)
      CHECK((s.mode[k] == od::TorsionMode::Free) == f.topos[i].torsions[k].rotatable);
  }
  CHECK(cfg->molecules[0].csgm1 == od::Csgm1::A);
  // 1L2Y: helix residues default to O, others to R.
  const auto& t0 = cfg->molecules[0].topology;
  for (std::size_t u = 0; u < t0.units.size(); ++u) {
    const char c = t0.dsspSS[t0.units[u].residue];
    CHECK(cfg->molecules[0].ssType[u] == (c == 'H' || c == 'G' || c == 'I' ? "O" : "R"));
  }
  CHECK(text.find("dsspSS") != std::string::npos);
  CHECK(text.find("alpha   = true") != std::string::npos);
  CHECK(text.find("delta   = false") != std::string::npos);
}

TEST_CASE("molconfig: validation errors") {
  ConfigFixture f;
  const std::string text = od::defaultMolConfigText(f.inputs);

  CHECK(hasError(f.errorsOf("formatVersion = 1\n[[molecule\n"), "line 2:"));
  CHECK(hasError(f.errorsOf(replaceAfter(text, "formatVersion", "1", "2")), "formatVersion: must be 1"));

  // Proline phi is in the pyrrolidine ring.
  auto e = f.errorsOf(replaceAfter(text, "resName = \"PRO\"", "phi     = false", "phi     = true"));
  CHECK(hasError(e, "residue A:12.phi: category 'ring' is not rotatable"));
  // Nucleic delta is in the ribose ring.
  e = f.errorsOf(replaceAfter(text, "mobm:1RNG", "delta   = false", "delta   = true"));
  CHECK(hasError(e, "residue A:1.delta: category 'ring'"));

  e = f.errorsOf(replaceAfter(text, "resSeq  = 5", "ssType  = \"O\"", "ssType  = \"Q\""));
  CHECK(hasError(e, "ssType 'Q' is not defined"));
  e = f.errorsOf(replaceAfter(text, "resSeq  = 6", "omega   = \"fixed\"", "omega   = \"trans\""));
  CHECK(hasError(e, "omega: must be fixed or cisTrans"));
  CHECK(f.errorsOf(replaceAfter(text, "resSeq  = 6", "omega   = \"fixed\"", "omega   = \"cisTrans\"")).empty());
  e = f.errorsOf(replaceAfter(text, "resSeq  = 7", "psi     = true", "psi     = 1"));
  CHECK(hasError(e, "psi: must be true or false"));
  e = f.errorsOf(replaceAfter(text, "resSeq  = 3", "iCode", "colour = 1\n  iCode"));
  CHECK(hasError(e, "residue A:3.colour: unknown key"));
  e = f.errorsOf(replaceAfter(text, "mobm:BTN", "backboneType = \"none\"", "backboneType = \"none\"\ncsgm1 = \"1a\""));
  CHECK(hasError(e, "csgm1: only allowed for protein molecules"));
  e = f.errorsOf(replaceAfter(text, "mobm:1L2Y", "csgm1        = \"1a\"", "csgm1        = \"1c\""));
  CHECK(hasError(e, "csgm1: must be 1a or 1b"));
  e = f.errorsOf(replaceAfter(text, "[[molecule]]\nid           = \"mobm:BTN", "mobm:BTN_ideal.sdf:1", "mobm:BTN_ideal.sdf:2"));
  CHECK(hasError(e, "id 'mobm:BTN_ideal.sdf:2' does not match any input molecule"));
  CHECK(hasError(e, "molecule: id 'mobm:BTN_ideal.sdf:1' is missing"));

  // Residue and torsion lists must match the detected ones.
  e = f.errorsOf(replaceAfter(text, "resSeq  = 2", "resSeq  = 2", "resSeq  = 99"));
  CHECK(hasError(e, "residue A:99 is not a backbone residue"));
  CHECK(hasError(e, "residue A:2 is missing"));
  e = f.errorsOf(replaceAfter(text, "mobm:BTN", "[\"#2\", \"#1\", \"#4\", \"#5\"]", "[\"#2\", \"#1\", \"#4\", \"#99\"]"));
  CHECK(hasError(e, "atom '#99' not found"));
  CHECK(hasError(e, "torsion [#2, #1, #4, #5] is missing"));

  // Custom ssType rules.
  const std::string custom = "[[ssType]]\nname = \"A\"\n[[ssType.point]]\nphi = -58\npsi = -48\nprob = 0.9\n"
                             "[[ssType.point]]\nphi = -58\npsi = -48\nprob = 0.2\n\n[[ssType]]\nname    = \"R\"";
  e = f.errorsOf(replaceAfter(text, "[[ssType]]", "[[ssType]]\nname    = \"R\"", custom));
  CHECK(hasError(e, "sum of prob is 1.1"));
  CHECK(hasError(e, "(phi, psi) = (-58, -48) is repeated"));
  e = f.errorsOf(replaceAfter(text, "name    = \"O\"", "name    = \"O\"\nbuiltin = true", "name    = \"O\"\nbuiltin = false"));
  CHECK(hasError(e, "built-in type 'O' is required"));
}

TEST_CASE("molconfig: ringBreak makes ring bonds rotatable (biotin)") {
  ConfigFixture f;
  std::string text = od::defaultMolConfigText(f.inputs);
  auto setTrue = [&](const std::string& atoms) {
    text = replaceAfter(text, "atoms     = " + atoms, "rotatable = false", "rotatable = true");
  };
  // Atom #9 is the thioether S of the tetrahydrothiophene ring.
  setTrue("[\"#7\", \"#8\", \"#9\", \"#10\"]");
  auto e = f.errorsOf(text);
  CHECK(hasError(e, "category 'ring' is not rotatable (a ringBreak in the ring changes this)"));
  text = replaceAfter(text, "mobm:BTN", "backboneType = \"none\"",
                      "backboneType = \"none\"\n[[molecule.ringBreak]]\nbond = [\"#9\", \"#10\"]");
  CHECK(f.errorsOf(text).empty());
  setTrue("[\"#8\", \"#9\", \"#10\", \"#11\"]");  // the broken bond itself
  setTrue("[\"#10\", \"#11\", \"#16\", \"#8\"]");  // still in the imidazolidinone ring
  e = f.errorsOf(text);
  CHECK(e.size() == 2);
  CHECK(hasError(e, "the bond is broken"));
  CHECK(hasError(e, "category 'ring'"));
  text = replaceAfter(text, "[[molecule.ringBreak]]", "bond = [\"#9\", \"#10\"]",
                      "bond = [\"#9\", \"#10\"]\n[[molecule.ringBreak]]\nbond = [\"#1\", \"#2\"]");
  CHECK(hasError(f.errorsOf(text), "ringBreak.bond: the bond is not in a ring"));
}

TEST_CASE("molconfig: default file name does not overwrite") {
  ConfigFixture f;
  const std::string prefix = (std::filesystem::temp_directory_path() / std::format("od_test_{}", ::getpid())).string();
  const std::string p0 = od::writeDefaultMolConfig(prefix, f.inputs);
  const std::string p1 = od::writeDefaultMolConfig(prefix, f.inputs);
  CHECK(p0 == prefix + ".molconfig.default.toml");
  CHECK(p1 == prefix + ".molconfig.default.1.toml");
  std::filesystem::remove(p0);
  std::filesystem::remove(p1);
}

namespace {

od::Molecule zincFinger() {
  static od::CcdLibrary lib;
  static bool loaded = false;
  if (!loaded) {
    lib.addFile(kData + "/1ZNF_ccd.cif");
    loaded = true;
  }
  return od::loadMolecules(kData + "/1ZNF.cif", od::Role::MobM, lib).front();
}

int ringBackbone(const od::Topology& t) {
  int n = 0;
  for (const auto& tor : t.torsions)
    n += (tor.backbone == od::BackboneAngle::Phi || tor.backbone == od::BackboneAngle::Psi) &&
         tor.category == od::TorsionCategory::Ring;
  return n;
}

// Replaces the value of the n-th occurrence of key (from the anchor) up to the closing quote.
std::string setQuoted(std::string s, const std::string& anchor, const std::string& key, const std::string& value) {
  auto i = s.find(anchor);
  REQUIRE(i != std::string::npos);
  i = s.find(key + "\"", i);
  REQUIRE(i != std::string::npos);
  i += key.size() + 1;
  const auto j = s.find('"', i);
  return s.replace(i, j - i, value);
}

}  // namespace

TEST_CASE("chem: Gasteiger charges with metal bonds (1ZNF)") {
  const auto m = zincFinger();
  for (const auto& a : m.atoms) {
    CHECK(std::isfinite(a.partialCharge));
    if (a.atomicNum == 30) CHECK(a.partialCharge == 2.0);
  }
  CHECK(sumPartialCharge(m) == doctest::Approx(netFormalCharge(m)).epsilon(1e-6));
}

TEST_CASE("topology: metal bonds in rings (1ZNF)") {
  const auto m = zincFinger();
  auto t = od::analyzeTopology(m, {});
  CHECK(t.crosslinkBreaks.empty());
  // Zn bonded to Cys3 SG, Cys6 SG, His19 NE2, His23 NE2: three bonds close rings.
  REQUIRE(t.metalBreaks.size() == 3);
  CHECK(t.brokenBonds.empty());
  CHECK(od::hasMetalRing(m));
  // Residue 3 psi through residue 23 phi lie in rings through Zn.
  CHECK(ringBackbone(t) == 40);
  od::TopologyParams brk;
  brk.brokenBonds = t.metalBreaks;
  CHECK(ringBackbone(od::analyzeTopology(m, brk)) == 0);
}

namespace {

// Sets phi of 1ZNF residue 10 to true.
std::string phi10(std::string s) {
  const auto r = s.find("resSeq  = 10\n", s.find("mobm:1ZNF.cif:A"));
  const auto p = s.find("phi     = false", r);
  REQUIRE(p < s.find("[[molecule.residue]]", r));
  return s.replace(p, 15, "phi     = true");
}

std::string uncommentMetal(std::string s) {
  for (const char* x : {"  # [[molecule.chainBreak]]", "  # bond   = ", "  # source = \"auto-metal\""})
    for (std::size_t p; (p = s.find(x)) != std::string::npos;) s.erase(p + 2, 2);
  return s;
}

}  // namespace

TEST_CASE("molconfig: metalBond") {
  ConfigFixture f;
  f.mols.push_back(zincFinger());
  f.topos.push_back(od::analyzeTopology(f.mols.back(), {}));
  od::assignDssp({&f.mols.back()}, {&f.topos.back()});
  f.inputs.clear();
  for (std::size_t i = 0; i < f.mols.size(); ++i) f.inputs.push_back({&f.mols[i], &f.topos[i]});
  const std::string text = od::defaultMolConfigText(f.inputs);
  const std::string zn = "id           = \"mobm:1ZNF.cif:A\"";
  // Default keep: the auto-metal entries are commented and the ring torsions are not rotatable.
  CHECK(text.find("metalBond    = \"keep\"") != std::string::npos);
  CHECK(text.find("  # bond   = [\"A:6:SG\", \"A:27:ZN\"]") != std::string::npos);
  CHECK(f.errorsOf(text).empty());
  CHECK(hasError(f.errorsOf(phi10(text)),
                 "residue A:10.phi: category 'ring' is not rotatable"));

  // keep with the entries uncommented: conflict.
  const std::string entries = uncommentMetal(text);
  auto e = f.errorsOf(entries);
  CHECK(std::count_if(e.begin(), e.end(), [](const std::string& x) {
          return x.find("metalBond = \"keep\" does not break metal bonds") != std::string::npos;
        }) == 3);

  // break: valid with the entries; the backbone may then be set rotatable.
  const std::string brk = setQuoted(entries, zn, "metalBond    = ", "break");
  CHECK(f.errorsOf(brk).empty());
  CHECK(f.errorsOf(phi10(brk)).empty());
  // break with the entries still commented or one missing: Zn is still in a ring.
  CHECK(hasError(f.errorsOf(setQuoted(text, zn, "metalBond    = ", "break")), "is still in a ring"));
  std::string missing = brk;
  const auto i = missing.find("  [[molecule.chainBreak]]\n  bond   = [\"A:23:NE2\"");
  missing.erase(i, missing.find("auto-metal\"\n", i) + 12 - i);
  CHECK(hasError(f.errorsOf(missing), "metal bond [A:3:SG, A:27:ZN] is still in a ring"));

  CHECK(hasError(f.errorsOf(setQuoted(text, zn, "metalBond    = ", "remove")), "metalBond: must be keep or break"));
  CHECK(hasError(f.errorsOf(setQuoted(brk, "bond   = [\"A:6:SG\"", "source = ", "user")),
                 "source: must be auto-metal for a metal bond"));
  e = f.errorsOf(replaceAfter(text, zn, "metalBond    = \"keep\"", ""));
  CHECK(hasError(e, "metalBond: missing"));
  e = f.errorsOf(replaceAfter(text, "mobm:BTN", "backboneType = \"none\"",
                              "backboneType = \"none\"\nmetalBond = \"keep\"\n"
                              "[[molecule.chainBreak]]\nbond = [\"#9\", \"#10\"]\nsource = \"auto-metal\"\n"
                              "[[molecule.chainBreak]]\nbond = [\"#1\", \"#2\"]\nsource = \"user\""));
  CHECK(hasError(e, "metalBond: only allowed when a metal bond lies in a ring"));
  CHECK(hasError(e, "source: auto-metal is only for metal bonds"));
  CHECK(hasError(e, "split the molecule into 2 parts"));
}

TEST_CASE("rama: embedded Top8000 grids") {
  const int listed[] = {27151, 29072, 9276, 6510, 17410, 21272};
  for (int i = 0; i < od::kRamaClassCount; ++i) {
    const auto& g = od::ramaGrid(static_cast<od::RamaClass>(i));
    CHECK(g.listed == listed[i]);
    CHECK(g.value.size() == od::kRamaBins * od::kRamaBins);
  }
  const auto& gen = od::ramaGrid(od::RamaClass::General);
  CHECK(gen.at(-179, -173) == doctest::Approx(0.0027118446100793205));
  CHECK(gen.at(-180, -174) == gen.at(-179, -173));
  CHECK(gen.at(180, 186) == gen.at(-179, -173));
  CHECK(gen.at(-63, -41) > 0.5f);  // alpha-helix region
  CHECK(od::RamaGrid::bin(-180) == 0);
  CHECK(od::RamaGrid::bin(180) == 0);
  CHECK(od::RamaGrid::bin(179.9) == 179);

  CHECK(od::ramaClassOf("GLY", "PRO", false) == od::RamaClass::Gly);
  CHECK(od::ramaClassOf("PRO", "PRO", true) == od::RamaClass::CisPro);
  CHECK(od::ramaClassOf("ILE", "PRO", false) == od::RamaClass::PrePro);
  CHECK(od::ramaClassOf("VAL", "", false) == od::RamaClass::IleVal);
  CHECK(od::ramaClassOf("ALA", "", false) == od::RamaClass::General);
  CHECK(od::ramaClassFromTypeName("RAMA_PREPRO") == od::RamaClass::PrePro);
  CHECK_FALSE(od::ramaClassFromTypeName("RAMA_AUTO"));
}

TEST_CASE("molconfig: RAMA types must match the residue (1L2Y)") {
  ConfigFixture f;
  const std::string text = od::defaultMolConfigText(f.inputs);
  auto ss = [&](const std::string& base, int seq, const std::string& type) {
    return setQuoted(base, std::format("resSeq  = {}\n", seq), "ssType  = ", type);
  };
  // NLYIQWLKDGGPSSGRPPPS
  std::string ok = ss(text, 2, "RAMA_GENERAL");
  ok = ss(ok, 4, "RAMA_ILEVAL");
  ok = ss(ok, 11, "RAMA_GLY");
  ok = ss(ok, 12, "RAMA_TRANSPRO");
  ok = ss(ok, 16, "RAMA_PREPRO");
  ok = ss(ok, 17, "RAMA_TRANSPRO");
  ok = ss(ok, 13, "RAMA_AUTO");
  CHECK(f.errorsOf(ok).empty());

  auto e = f.errorsOf(ss(text, 2, "RAMA_GLY"));
  CHECK(hasError(e, "residue A:2.ssType: RAMA_GLY is for Gly residues; residue A:2 LEU is general"));
  CHECK(hasError(f.errorsOf(ss(text, 4, "RAMA_GENERAL")), "use RAMA_ILEVAL or RAMA_AUTO"));
  CHECK(hasError(f.errorsOf(ss(text, 12, "RAMA_CISPRO")), "residue A:12 PRO is trans-Pro"));
  CHECK(hasError(f.errorsOf(ss(text, 16, "RAMA_GENERAL")), "use RAMA_PREPRO or RAMA_AUTO"));
  const std::string cisTrans =
      replaceAfter(ss(text, 12, "RAMA_TRANSPRO"), "resSeq  = 11\n", "omega   = \"fixed\"", "omega   = \"cisTrans\"");
  CHECK(hasError(f.errorsOf(cisTrans), "the omega before this Pro is cisTrans; use RAMA_AUTO"));
  CHECK(f.errorsOf(ss(cisTrans, 12, "RAMA_AUTO")).empty());

  // Chain ends default to R.
  CHECK(setQuoted(text, "resSeq  = 1\n", "ssType  = ", "X") == ss(text, 1, "X"));
  const auto r1 = text.find("ssType  = ", text.find("resSeq  = 1\n"));
  const auto r20 = text.find("ssType  = ", text.find("resSeq  = 20\n"));
  CHECK(text.substr(r1, 13) == "ssType  = \"R\"");
  CHECK(text.substr(r20, 13) == "ssType  = \"R\"");
}

TEST_CASE("molconfig: ssHalfWidth and ssOriginalProb") {
  ConfigFixture f;
  const std::string text = od::defaultMolConfigText(f.inputs);
  CHECK(text.find("ssHalfWidth    = 15.0\n") != std::string::npos);
  CHECK(text.find("ssOriginalProb = 0.9\n") != std::string::npos);
  auto read = [&](const std::string& s) {
    std::ofstream(f.path) << s;
    std::vector<std::string> errors;
    auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
    REQUIRE(cfg);
    return *cfg;
  };
  auto cfg = read(text);
  CHECK(cfg.ssHalfWidth == 15.0);
  CHECK(cfg.ssOriginalProb == 0.9);
  cfg = read(replaceAfter(replaceAfter(text, "ssHalfWidth", "15.0", "10"), "ssOriginalProb", "0.9", "1"));
  CHECK(cfg.ssHalfWidth == 10.0);
  CHECK(cfg.ssOriginalProb == 1.0);
  std::string none = replaceAfter(text, "ssHalfWidth    = 15.0\n", "ssOriginalProb = 0.9\n", "");
  none = replaceAfter(none, "# ssOriginalProb", "ssHalfWidth    = 15.0\n", "");
  cfg = read(none);
  CHECK(cfg.ssHalfWidth == 15.0);
  CHECK(cfg.ssOriginalProb == 0.9);
  CHECK(hasError(f.errorsOf(replaceAfter(text, "ssHalfWidth", "15.0", "0")), "ssHalfWidth: must be > 0 and <= 90"));
  CHECK(hasError(f.errorsOf(replaceAfter(text, "ssHalfWidth", "15.0", "91")), "ssHalfWidth: must be > 0 and <= 90"));
  CHECK(hasError(f.errorsOf(replaceAfter(text, "ssOriginalProb", "0.9", "1.5")), "ssOriginalProb: must be > 0 and <= 1"));
  CHECK(hasError(f.errorsOf(replaceAfter(text, "ssOriginalProb", "0.9", "\"x\"")), "line 17: ssOriginalProb: must be a number"));
}

namespace {

double wrapDeg(double x) {
  x = std::fmod(x + 180.0, 360.0);
  if (x < 0) x += 360.0;
  return x - 180.0;
}

double distance(const od::Vec3& a, const od::Vec3& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

std::vector<double> inputAngles(const od::TorsionTree& tree) {
  std::vector<double> v;
  for (const auto& tt : tree.torsions) v.push_back(tt.inputAngle);
  return v;
}

// One StaM (2KJJ, chains A+B) and one MobM (1L2Y).
struct StaFixture {
  std::vector<od::Molecule> mols;
  std::vector<od::Topology> topos;
  std::vector<od::MolInput> inputs;
  od::TopologyParams params;
  std::string path;
  std::string text;

  explicit StaFixture(bool stamFixedResMovable = false) {
    params.stamFixedResMovable = stamFixedResMovable;
    mols.push_back(od::loadMolecules(kData + "/2KJJ.cif", od::Role::StaM, ccd()).front());
    mols.push_back(od::loadMolecules(kData + "/1L2Y.cif", od::Role::MobM, ccd()).front());
    for (auto& m : mols) {
      topos.push_back(od::analyzeTopology(m, params));
      od::assignDssp({&m}, {&topos.back()});
    }
    for (std::size_t i = 0; i < mols.size(); ++i) inputs.push_back({&mols[i], &topos[i]});
    path = (std::filesystem::temp_directory_path() / std::format("od_test_tree_{}.toml", ::getpid())).string();
    text = od::defaultMolConfigText(inputs);
  }
  ~StaFixture() { std::filesystem::remove(path); }

  std::optional<od::MolConfig> read(const std::string& s, std::vector<std::string>& errors) {
    std::ofstream(path) << s;
    errors.clear();
    return od::readMolConfig(path, inputs, params, errors);
  }
  od::MolConfig ok(const std::string& s) {
    std::vector<std::string> errors;
    auto cfg = read(s, errors);
    for (const auto& x : errors) MESSAGE(x);
    REQUIRE(cfg);
    return *cfg;
  }
  std::vector<std::string> errorsOf(const std::string& s) {
    std::vector<std::string> errors;
    CHECK_FALSE(read(s, errors));
    return errors;
  }
};

// Sets phi and psi of StaM residues chain:first..last to rotatable.
std::string freeResidues(std::string s, const std::string& chain, int first, int last) {
  const auto start = s.find("id           = \"stam:");
  for (int seq = first; seq <= last; ++seq) {
    auto i = s.find(std::format("chain   = \"{}\"\n  resName", chain), start);
    i = s.find(std::format("resSeq  = {}\n", seq), i);
    REQUIRE(i != std::string::npos);
    const auto j = s.find("[[molecule.residue]]", i);
    for (const char* a : {"phi     = false", "psi     = false"})
      if (auto p = s.find(a, i); p < j) s.replace(p, 15, std::string(a).replace(10, 5, "true"));
  }
  return s;
}

std::string addBreak(std::string s, const std::string& a, const std::string& b, const std::string& source) {
  const auto i = s.find("csgm1        = \"1a\"\n", s.find("id           = \"stam:"));
  REQUIRE(i != std::string::npos);
  return s.insert(i + 20, std::format("\n  [[molecule.chainBreak]]\n  bond   = [\"{}\", \"{}\"]\n  source = \"{}\"\n", a, b,
                                      source));
}

}  // namespace

TEST_CASE("torsion tree: MobM coordinates (1L2Y)") {
  ConfigFixture f;
  const std::string text = od::defaultMolConfigText(f.inputs);
  std::ofstream(f.path) << text;
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  REQUIRE(cfg);
  const od::Molecule& m = f.mols[0];
  const od::MolSettings& s = cfg->molecules[0];
  const od::TorsionTree& tree = s.tree;
  CHECK(tree.torsions.size() == static_cast<std::size_t>(std::count(s.mode.begin(), s.mode.end(), od::TorsionMode::Free)));
  REQUIRE(tree.anchors.size() == 1);
  CHECK(tree.loops.empty());
  CHECK(tree.virtuals.empty());
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    const auto& tt = tree.torsions[i];
    CHECK(!std::binary_search(tt.moving.begin(), tt.moving.end(), tree.anchors.front()));
    if (tt.parent >= 0) {
      CHECK(static_cast<std::size_t>(tt.parent) < i);
      const auto& pm = tree.torsions[tt.parent].moving;
      CHECK(std::includes(pm.begin(), pm.end(), tt.moving.begin(), tt.moving.end()));
    }
  }

  // Input angles reproduce the input coordinates.
  auto xyz = od::buildCoordinates(m, tree, inputAngles(tree));
  double maxShift = 0;
  for (std::size_t a = 0; a < m.atoms.size(); ++a) maxShift = std::max(maxShift, distance(xyz[a], m.atoms[a].xyz));
  CHECK(maxShift < 1e-9);

  // Arbitrary angles: every tree torsion takes its value; bonds and the anchor are unchanged.
  std::vector<double> angles;
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) angles.push_back(wrapDeg(-170.0 + 37.0 * i));
  xyz = od::buildCoordinates(m, tree, angles);
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    const auto& x = s.topology.torsions[tree.torsions[i].torsion].atoms;
    CHECK(wrapDeg(od::dihedralDeg(xyz[x[0]], xyz[x[1]], xyz[x[2]], xyz[x[3]]) - angles[i]) == doctest::Approx(0).epsilon(1e-6));
  }
  double bondErr = 0;
  for (const auto& b : m.bonds)
    bondErr = std::max(bondErr, std::abs(distance(xyz[b.a], xyz[b.b]) - distance(m.atoms[b.a].xyz, m.atoms[b.b].xyz)));
  CHECK(bondErr < 1e-9);
  CHECK(distance(xyz[tree.anchors.front()], m.atoms[tree.anchors.front()].xyz) < 1e-12);
  // Fixed torsions keep their input values.
  for (std::size_t ti = 0; ti < s.topology.torsions.size(); ++ti) {
    if (s.mode[ti] != od::TorsionMode::Fixed) continue;
    const auto& tor = s.topology.torsions[ti];
    const auto& x = tor.atoms;
    CHECK(wrapDeg(od::dihedralDeg(xyz[x[0]], xyz[x[1]], xyz[x[2]], xyz[x[3]]) - tor.angle) == doctest::Approx(0).epsilon(1e-6));
  }
}

TEST_CASE("torsion tree: single-precision coordinates of a long chain (1D3Z)") {
  std::vector<od::Molecule> mols{od::loadMolecules(kData + "/1D3Z.cif", od::Role::MobM, ccd()).front()};
  std::vector<od::Topology> topos{od::analyzeTopology(mols[0], {})};
  od::assignDssp({&mols[0]}, {&topos[0]});
  const std::vector<od::MolInput> inputs{{&mols[0], &topos[0]}};
  const auto path = (std::filesystem::temp_directory_path() / std::format("od_test_1d3z_{}.toml", ::getpid())).string();
  std::ofstream(path) << od::defaultMolConfigText(inputs);
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(path, inputs, {}, errors);
  std::filesystem::remove(path);
  REQUIRE(cfg);
  const od::TorsionTree& tree = cfg->molecules[0].tree;
  int depth = 0;
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    int d = 0;
    for (int p = static_cast<int>(i); p >= 0; p = tree.torsions[p].parent) ++d;
    depth = std::max(depth, d);
  }
  CHECK(depth > 50);
  double worst = 0, worstBond = 0;
  for (std::uint64_t s = 0; s < 5; ++s) {
    od::Rng rng(11, 0, s);
    std::vector<double> ad(tree.torsions.size());
    for (double& v : ad) v = static_cast<float>(rng.uniform() * 360 - 180);
    const std::vector<float> af(ad.begin(), ad.end());
    const auto xd = od::buildCoordinates(mols[0], tree, ad);
    const auto xf = od::buildCoordinates(mols[0], tree, std::span<const float>(af));
    for (std::size_t i = 0; i < xd.size(); ++i)
      for (int k = 0; k < 3; ++k) worst = std::max(worst, std::abs(xd[i][k] - xf[i][k]));
    for (const auto& b : mols[0].bonds) {
      auto len = [&](const auto& x) {
        const double dx = x[b.a][0] - x[b.b][0], dy = x[b.a][1] - x[b.b][1], dz = x[b.a][2] - x[b.b][2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
      };
      worstBond = std::max(worstBond, std::abs(len(xd) - len(xf)));
    }
  }
  MESSAGE("1D3Z: ", tree.torsions.size(), " torsions, depth ", depth, ", max |float - double| ", worst,
          " Å, max bond-length change ", worstBond, " Å");
  CHECK(worst < 1e-2);
  CHECK(worstBond < 1e-3);
}

TEST_CASE("torsion tree: StaM defaults (2KJJ)") {
  StaFixture f;
  const auto cfg = f.ok(f.text);
  const od::MolSettings& s = cfg.molecules[0];
  // Every backbone residue is fixed: anchors are all backbone atoms, no loops.
  CHECK(s.tree.loops.empty());
  CHECK(s.tree.anchors.size() == 4 * s.topology.units.size());
  // Two disulfides are broken (A6-A11 and one of the two inter-chain ones): 3 pairs each.
  CHECK(s.topology.crosslinkBreaks.size() == 2);
  CHECK(s.tree.virtuals.size() == 6);
  auto xyz = od::buildCoordinates(f.mols[0], s.tree, inputAngles(s.tree));
  const auto c = od::evaluateClosure(s.tree, xyz, 200);
  CHECK(c.energy < 1e-12);
  CHECK(c.maxDeviation < 1e-9);
  // The kept inter-chain disulfide: its torsions have backbone atoms on both sides and are not rotatable by default.
  int restricted = -1;
  for (std::size_t i = 0; i < s.topology.torsions.size(); ++i) {
    const auto& tor = s.topology.torsions[i];
    if (tor.backbone == od::BackboneAngle::None && tor.category == od::TorsionCategory::Single && !tor.rotatable &&
        f.mols[0].atoms[tor.atoms[1]].atomicNum == 16 && f.mols[0].atoms[tor.atoms[2]].atomicNum == 16)
      restricted = static_cast<int>(i);
  }
  REQUIRE(restricted >= 0);
  // A7 chi2 (kept disulfide A7 SG - B7 SG).
  const std::string chi2 = "atoms     = [\"CA\", \"CB\", \"SG\", \"B:7:SG\"]\n    rotatable = false";
  const auto p = f.text.find(chi2);
  REQUIRE(p != std::string::npos);
  std::string edited = f.text;
  edited.replace(p + chi2.size() - 5, 5, "true");
  CHECK(hasError(f.errorsOf(edited), "cannot rotate: fixed atoms on both sides"));
}

TEST_CASE("torsion tree: StaM loops (2KJJ)") {
  StaFixture f;
  const od::Molecule& m = f.mols[0];
  // A10-A14 between fixed residues: the break is added at A14 C - A15 N.
  const std::string loop = freeResidues(f.text, "A", 10, 14);
  auto cfg = f.ok(loop);
  const od::MolSettings& s = cfg.molecules[0];
  REQUIRE(s.tree.loops.size() == 1);
  const od::Loop& l = s.tree.loops.front();
  CHECK(l.units.size() == 5);
  CHECK(l.before >= 0);
  CHECK(l.after >= 0);
  REQUIRE(l.breakBond >= 0);
  CHECK(m.atoms[m.bonds[l.breakBond].a].name + m.atoms[m.bonds[l.breakBond].b].name == "CN");
  CHECK(m.residues[m.atoms[m.bonds[l.breakBond].a].residue].info.seq == 14);
  REQUIRE(s.chainBreaks.size() == 3);
  CHECK(s.chainBreaks.back().source == "auto-loop");
  CHECK(s.chainBreaks.back().line == 0);
  CHECK(s.tree.virtuals.size() == 9);

  // Rotating loop torsions moves loop atoms only; the fixed backbone stays; the closure opens.
  std::vector<double> angles = inputAngles(s.tree);
  for (std::size_t i = 0; i < angles.size(); ++i)
    if (s.topology.torsions[s.tree.torsions[i].torsion].backbone != od::BackboneAngle::None) angles[i] += 20;
  auto xyz = od::buildCoordinates(m, s.tree, angles);
  for (int a : s.tree.anchors) CHECK(distance(xyz[a], m.atoms[a].xyz) < 1e-12);
  const auto c = od::evaluateClosure(s.tree, xyz, 200);
  CHECK(c.maxDeviation > 0.5);
  CHECK(c.energy >= 200 * c.maxDeviation * c.maxDeviation);

  // A user break elsewhere in the loop replaces the default.
  cfg = f.ok(addBreak(loop, "A:11:C", "A:12:N", "user"));
  CHECK(cfg.molecules[0].chainBreaks.size() == 3);
  CHECK(cfg.molecules[0].tree.loops.front().breakBond == cfg.molecules[0].chainBreaks.front().bond);
  CHECK(cfg.molecules[0].chainBreaks.front().source == "user");
  // The added break written into the file is accepted.
  cfg = f.ok(addBreak(loop, "A:14:C", "A:15:N", "auto-loop"));
  CHECK(cfg.molecules[0].chainBreaks.size() == 3);
  CHECK(cfg.molecules[0].chainBreaks.front().line > 0);
  // Two breaks in one loop leave a free part.
  CHECK(hasError(f.errorsOf(addBreak(addBreak(loop, "A:11:C", "A:12:N", "user"), "A:12:C", "A:13:N", "user")),
                 "leave a part without fixed atoms"));
  // auto-loop outside a loop.
  CHECK(hasError(f.errorsOf(addBreak(loop, "A:3:C", "A:4:N", "auto-loop")), "auto-loop is only for a backbone bond"));
  // A cisTrans omega at the default break.
  CHECK(hasError(f.errorsOf(replaceAfter(loop, "resSeq  = 14\n", "omega   = \"fixed\"", "omega   = \"cisTrans\"")),
                 "the default break [A:14:C, A:15:N] has a rotatable omega"));

  // N-terminal loop: free end, no break.
  cfg = f.ok(freeResidues(f.text, "A", 1, 3));
  REQUIRE(cfg.molecules[0].tree.loops.size() == 1);
  CHECK(cfg.molecules[0].tree.loops.front().before < 0);
  CHECK(cfg.molecules[0].tree.loops.front().breakBond < 0);
  CHECK(cfg.molecules[0].chainBreaks.size() == 2);
}

TEST_CASE("torsion tree: stamFixedResMovable (2KJJ)") {
  StaFixture f(true);
  const auto cfg = f.ok(freeResidues(f.text, "A", 10, 14));
  const od::MolSettings& s = cfg.molecules[0];
  CHECK(s.tree.loops.empty());
  REQUIRE(s.tree.anchors.size() == 1);
  CHECK(f.mols[0].atoms[s.tree.anchors.front()].name == "CA");
  CHECK(s.chainBreaks.size() == 2);
}

TEST_CASE("torsion tree: metal restraints (1ZNF)") {
  ConfigFixture f;
  f.mols.push_back(zincFinger());
  f.topos.push_back(od::analyzeTopology(f.mols.back(), {}));
  od::assignDssp({&f.mols.back()}, {&f.topos.back()});
  f.inputs.clear();
  for (std::size_t i = 0; i < f.mols.size(); ++i) f.inputs.push_back({&f.mols[i], &f.topos[i]});
  const std::string text = od::defaultMolConfigText(f.inputs);
  const std::string brk = setQuoted(uncommentMetal(text), "id           = \"mobm:1ZNF.cif:A\"", "metalBond    = ", "break");
  std::ofstream(f.path) << phi10(brk);
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  REQUIRE(cfg);
  const auto& s = cfg->molecules.back();
  REQUIRE(s.tree.restraints.size() == 3);
  CHECK(s.tree.virtuals.empty());
  auto angles = inputAngles(s.tree);
  auto xyz = od::buildCoordinates(f.mols.back(), s.tree, angles);
  CHECK(od::evaluateClosure(s.tree, xyz, 200).maxDeviation < 1e-9);
  for (double& a : angles) a += 15;
  xyz = od::buildCoordinates(f.mols.back(), s.tree, angles);
  CHECK(od::evaluateClosure(s.tree, xyz, 200).maxDeviation > 0.5);
}

namespace {

int atomOf(const od::Molecule& m, int seq, const std::string& name) {
  for (std::size_t a = 0; a < m.atoms.size(); ++a)
    if (m.atoms[a].residue >= 0 && m.residues[m.atoms[a].residue].info.seq == seq && m.atoms[a].name == name)
      return static_cast<int>(a);
  return -1;
}

// Number of torsions whose analytic gradient differs from the numeric one.
int gradientMismatches(const od::EnergySystem& sys, std::vector<double> angles) {
  std::vector<double> g;
  sys.evaluate(angles, &g);
  const double h = 1e-4;  // degrees
  int bad = 0;
  for (std::size_t i = 0; i < angles.size(); ++i) {
    const double a = angles[i];
    auto at = [&](double d) {
      angles[i] = a + d;
      return sys.evaluate(angles).total();
    };
    // Fourth-order central difference: clashing pairs give large higher derivatives.
    const double fd = (8 * (at(h) - at(-h)) - (at(2 * h) - at(-2 * h))) / (12 * h) * 180.0 / std::numbers::pi;
    angles[i] = a;
    if (std::abs(fd - g[i]) > 1e-5 * std::max(1.0, std::abs(g[i])) + 1e-6) {
      ++bad;
      MESSAGE(std::format("torsion {}: analytic {} numeric {}", i, g[i], fd));
    }
  }
  return bad;
}

std::vector<double> perturbed(std::vector<double> v, double width, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-width, width);
  for (double& x : v) x += u(rng);
  return v;
}

}  // namespace

TEST_CASE("energy: AD4 pair functions") {
  const od::Ad4Coeffs c;
  double e, d;
  auto vdw = od::ad4PairParams(od::Ad4Type::C, od::Ad4Type::OA, c);
  CHECK(vdw.n == 6);
  CHECK(vdw.reqm == doctest::Approx(3.6));
  od::ad4VdwHb(vdw.reqm, vdw.AC, vdw.BD, vdw.m, vdw.n, vdw.reqm, 0.25, e, d);
  CHECK(e == doctest::Approx(-0.1662 * std::sqrt(0.15 * 0.2)));
  CHECK(d == 0);
  auto hb = od::ad4PairParams(od::Ad4Type::HD, od::Ad4Type::OA, c);
  CHECK(hb.n == 10);
  CHECK(hb.reqm == doctest::Approx(1.9));
  od::ad4VdwHb(2.0, hb.AC, hb.BD, hb.m, hb.n, hb.reqm, 0.25, e, d);
  CHECK(e == doctest::Approx(-0.1209 * 5.0));
  CHECK(od::ad4PairParams(od::Ad4Type::HD, od::Ad4Type::SA, c).reqm == doctest::Approx(2.5));
  CHECK(od::ad4PairParams(od::Ad4Type::HD, od::Ad4Type::N, c).n == 6);
  // Desolvation: fitted Gaussian is 1 at r = 0.
  od::ad4Desolv(0.0, 2.0, 0.0, e, d);
  CHECK(e == doctest::Approx(2.0));
  // Derivatives by central differences.
  auto fd = [](auto f, double r) {
    double ep, em, dd;
    f(r + 1e-6, ep, dd);
    f(r - 1e-6, em, dd);
    return (ep - em) / 2e-6;
  };
  for (double r : {0.9, 1.5, 3.0, 4.2, 7.5, 12.0, 19.0}) {
    auto fv = [&](double x, double& ee, double& dd) { od::ad4VdwHb(x, vdw.AC, vdw.BD, 12, 6, vdw.reqm, 0.25, ee, dd); };
    auto fe = [&](double x, double& ee, double& dd) { od::ad4Elec(x, -10.0, 0.01, ee, dd); };
    auto fs = [&](double x, double& ee, double& dd) { od::ad4Desolv(x, 0.5, 0.01, ee, dd); };
    fv(r, e, d);
    CHECK(d == doctest::Approx(fd(fv, r)).epsilon(1e-6));
    fe(r, e, d);
    CHECK(d == doctest::Approx(fd(fe, r)).epsilon(1e-6));
    fs(r, e, d);
    CHECK(d == doctest::Approx(fd(fs, r)).epsilon(1e-6));
  }
  od::clashPenalty(1.0, 1.8, 10000.0, 1000.0, e, d);
  CHECK(e == doctest::Approx(10000 + 1000 * 0.64));
  CHECK(d == doctest::Approx(-1600));
}

TEST_CASE("energy: AD4 types (1L2Y, biotin)") {
  const auto m = od::loadMolecules(kData + "/1L2Y.cif", od::Role::MobM, ccd()).front();
  const auto t = od::assignAd4Types(m);
  CHECK(t.unsupported.empty());
  auto type = [&](int seq, const char* name) {
    const int a = atomOf(m, seq, name);
    REQUIRE(a >= 0);
    return std::string(od::ad4TypeName(t.type[a]));
  };
  CHECK(type(1, "N") == "N");     // N-terminal NH3+
  CHECK(type(2, "N") == "N");     // amide
  CHECK(type(2, "O") == "OA");
  CHECK(type(1, "ND2") == "N");   // Asn amide
  CHECK(type(1, "OD1") == "OA");
  CHECK(type(5, "NE2") == "N");   // Gln amide
  CHECK(type(3, "OH") == "OA");   // Tyr
  CHECK(type(3, "CZ") == "A");
  CHECK(type(3, "CB") == "C");
  CHECK(type(6, "NE1") == "N");   // Trp pyrrole N-H
  CHECK(type(6, "CD1") == "A");
  CHECK(type(8, "NZ") == "N");    // Lys NH3+
  CHECK(type(9, "OD1") == "OA");
  CHECK(type(16, "NH1") == "N");  // Arg
  CHECK(type(16, "NE") == "N");
  CHECK(type(6, "HE1") == "HD");
  CHECK(type(13, "HG") == "HD");  // Ser OG-H
  CHECK(type(6, "HA") == "-");
  // Merged hydrogens move their charge to the heavy atom.
  double sumIn = 0, sumOut = 0;
  for (std::size_t a = 0; a < m.atoms.size(); ++a) {
    sumIn += m.atoms[a].partialCharge;
    sumOut += t.charge[a];
    if (t.type[a] == od::Ad4Type::Merged) CHECK(t.charge[a] == 0);
  }
  CHECK(sumOut == doctest::Approx(sumIn));
  const int ca = atomOf(m, 6, "CA"), ha = atomOf(m, 6, "HA");
  CHECK(t.charge[ca] == doctest::Approx(m.atoms[ca].partialCharge + m.atoms[ha].partialCharge));

  const auto btn = od::loadMolecules(kData + "/BTN_ideal.sdf", od::Role::MobM, ccd()).front();
  const auto tb = od::assignAd4Types(btn);
  for (std::size_t a = 0; a < btn.atoms.size(); ++a) {
    const int z = btn.atoms[a].atomicNum;
    if (z == 16) CHECK(tb.type[a] == od::Ad4Type::SA);  // thioether
    if (z == 7) CHECK(tb.type[a] == od::Ad4Type::N);    // ureido
    if (z == 8) CHECK(tb.type[a] == od::Ad4Type::OA);
  }
}

TEST_CASE("energy: pairs, gradient and clash (1L2Y)") {
  ConfigFixture f;
  std::ofstream(f.path) << od::defaultMolConfigText(f.inputs);
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  REQUIRE(cfg);
  const od::Molecule& m = f.mols[0];
  const od::TorsionTree& tree = cfg->molecules[0].tree;
  const od::EnergySystem sys({{&m, &tree}}, {});
  CHECK(sys.torsionCount() == static_cast<int>(tree.torsions.size()));

  // Pairs: no atoms within three bonds; every excluded pair keeps its distance.
  const int n = static_cast<int>(m.atoms.size());
  std::vector<std::vector<int>> adj(n);
  for (const auto& b : m.bonds) {
    adj[b.a].push_back(b.b);
    adj[b.b].push_back(b.a);
  }
  std::vector<std::vector<int>> hops(n, std::vector<int>(n, 99));
  for (int s = 0; s < n; ++s) {
    std::vector<int> q{s};
    hops[s][s] = 0;
    for (std::size_t i = 0; i < q.size(); ++i)
      for (int o : adj[q[i]])
        if (hops[s][o] == 99) {
          hops[s][o] = hops[s][q[i]] + 1;
          q.push_back(o);
        }
  }
  std::vector<std::vector<bool>> listed(n, std::vector<bool>(n, false));
  for (const auto& p : sys.pairs()) {
    CHECK(hops[p.i][p.j] > 3);
    listed[p.i][p.j] = listed[p.j][p.i] = true;
  }
  std::vector<od::Vec3> x0, x1;
  sys.evaluate(perturbed(sys.inputAngles(), 180, 1), nullptr, &x0);
  sys.evaluate(perturbed(sys.inputAngles(), 180, 2), nullptr, &x1);
  const auto& type = sys.typing(0).type;
  int excluded = 0, notConstant = 0;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) {
      if (type[i] == od::Ad4Type::Merged || type[j] == od::Ad4Type::Merged || hops[i][j] <= 3 || listed[i][j]) continue;
      ++excluded;
      if (std::abs(distance(x0[i], x0[j]) - distance(x1[i], x1[j])) > 1e-9) ++notConstant;
    }
  CHECK(excluded > 0);
  CHECK(notConstant == 0);

  // Input conformation: no clash; analytic gradient matches numeric.
  const auto e0 = sys.evaluate(sys.inputAngles());
  CHECK(e0.clashPairs == 0);
  CHECK(e0.closure == 0);
  CHECK(e0.intra() < 0);
  CHECK(gradientMismatches(sys, sys.inputAngles()) == 0);
  CHECK(gradientMismatches(sys, perturbed(sys.inputAngles(), 20, 3)) == 0);

  // A random conformation with mild clashes (severe overlaps make the numeric gradient inaccurate).
  bool found = false;
  for (unsigned seed = 10; seed < 2000 && !found; ++seed) {
    const auto a = perturbed(sys.inputAngles(), 40, seed);
    const auto e = sys.evaluate(a);
    if (e.clashPairs == 0 || e.total() > 1e5) continue;
    found = true;
    CHECK(e.clash >= 10000.0 * e.clashPairs);
    CHECK(e.filterEnergy() > 1000);
    CHECK(gradientMismatches(sys, a) == 0);
  }
  CHECK(found);
}

TEST_CASE("energy: StaM loop closure (2KJJ)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 10, 14));
  const od::EnergySystem sys({{&f.mols[0], &cfg.molecules[0].tree}}, {});
  const auto e0 = sys.evaluate(sys.inputAngles());
  CHECK(e0.closure < 1e-12);
  CHECK(e0.closureMaxDeviation < 1e-9);
  // Pairs between two fixed backbone atoms are constant and not listed.
  std::vector<bool> anchor(f.mols[0].atoms.size(), false);
  for (int a : cfg.molecules[0].tree.anchors) anchor[a] = true;
  for (const auto& p : sys.pairs()) CHECK_FALSE((anchor[p.i] && anchor[p.j]));
  const auto a = perturbed(sys.inputAngles(), 10, 5);
  const auto e = sys.evaluate(a);
  CHECK(e.closure > 0);
  CHECK(e.total() == doctest::Approx(e.filterEnergy() + e.closure));
  CHECK(gradientMismatches(sys, a) == 0);
}

TEST_CASE("energy: metal restraints (1ZNF)") {
  ConfigFixture f;
  f.mols.push_back(zincFinger());
  f.topos.push_back(od::analyzeTopology(f.mols.back(), {}));
  od::assignDssp({&f.mols.back()}, {&f.topos.back()});
  f.inputs.clear();
  for (std::size_t i = 0; i < f.mols.size(); ++i) f.inputs.push_back({&f.mols[i], &f.topos[i]});
  const std::string text = od::defaultMolConfigText(f.inputs);
  std::ofstream(f.path) << phi10(setQuoted(uncommentMetal(text), "id           = \"mobm:1ZNF.cif:A\"", "metalBond    = ", "break"));
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  REQUIRE(cfg);
  const od::EnergySystem sys({{&f.mols.back(), &cfg->molecules.back().tree}}, {});
  const auto& t = sys.typing(0);
  for (std::size_t a = 0; a < f.mols.back().atoms.size(); ++a)
    if (f.mols.back().atoms[a].atomicNum == 30) CHECK(t.type[a] == od::Ad4Type::Zn);
  const auto a = perturbed(sys.inputAngles(), 10, 7);
  CHECK(sys.evaluate(a).closure > 0);
  CHECK(gradientMismatches(sys, a) == 0);
}

TEST_CASE("energy: clash between fixed atoms in the input (2KJJ)") {
  StaFixture f;
  const auto cfg = f.ok(f.text);
  od::Molecule m = f.mols[0];
  const auto& tree = cfg.molecules[0].tree;
  CHECK_NOTHROW(od::EnergySystem({{&m, &tree}}, {}));
  // Move one fixed backbone atom onto a distant fixed backbone atom.
  const int a = tree.anchors.front();
  int b = -1;
  for (int c : tree.anchors)
    if (distance(m.atoms[a].xyz, m.atoms[c].xyz) > 10) b = c;
  REQUIRE(b >= 0);
  m.atoms[a].xyz = {m.atoms[b].xyz[0] + 1.0, m.atoms[b].xyz[1], m.atoms[b].xyz[2]};
  try {
    od::EnergySystem sys({{&m, &tree}}, {});
    FAIL("no error");
  } catch (const od::Error& e) {
    CHECK(std::string(e.what()).find("distance cannot change") != std::string::npos);
    CHECK(std::string(e.what()).find("1.000 Å") != std::string::npos);
  }
}

namespace {

// 1L2Y as MobM with the omega of Gly 11 (before Pro 12) set to cisTrans, and
// the Pro 12 ring broken at CD-N so that its phi is rotatable.
struct Csgm1Fixture {
  ConfigFixture f;
  std::optional<od::MolConfig> cfg;
  Csgm1Fixture() {
    std::string text = od::defaultMolConfigText(f.inputs);
    const auto mol = text.find("id           = \"mobm:1L2Y.cif:A\"");
    auto i = text.find("resSeq  = 11\n", mol);
    REQUIRE(i != std::string::npos);
    i = text.find("omega   = \"fixed\"", i);
    REQUIRE(i != std::string::npos);
    text.replace(i, 17, "omega   = \"cisTrans\"");
    i = text.find("resSeq  = 12\n", mol);
    REQUIRE(i != std::string::npos);
    i = text.find("phi     = false", i);
    REQUIRE(i != std::string::npos);
    text.replace(i, 15, "phi     = true");
    i = text.find("csgm1        = \"1a\"\n", mol);
    REQUIRE(i != std::string::npos);
    text.insert(i + 20, "\n  [[molecule.ringBreak]]\n  bond = [\"A:12:CD\", \"A:12:N\"]\n");
    std::ofstream(f.path) << text;
    std::vector<std::string> errors;
    cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
    for (const auto& e : errors) MESSAGE(e);
    REQUIRE(cfg);
  }
  od::MolSettings& s() { return cfg->molecules[0]; }
  const od::Molecule& m() { return f.mols[0]; }
  int unitOf(int seq) {
    const auto& t = s().topology;
    for (std::size_t u = 0; u < t.units.size(); ++u)
      if (m().residues[t.units[u].residue].info.seq == seq) return static_cast<int>(u);
    return -1;
  }
  int treeIndex(int seq, od::BackboneAngle a) {
    const auto& t = s().topology;
    for (std::size_t k = 0; k < s().tree.torsions.size(); ++k) {
      const auto& tor = t.torsions[s().tree.torsions[k].torsion];
      if (tor.backbone == a && m().residues[t.units[tor.unit].residue].info.seq == seq) return static_cast<int>(k);
    }
    return -1;
  }
};

}  // namespace

TEST_CASE("csgm1: csgm1a angles and cisTrans omega") {
  Csgm1Fixture x;
  const int omega = x.treeIndex(11, od::BackboneAngle::Omega);
  REQUIRE(omega >= 0);
  for (double step : {1.0, 7.5, 100.0}) {
    const od::Csgm1Sampler smp(x.m(), x.s(), *x.cfg, {step, 1e-4});
    REQUIRE(smp.size() == static_cast<int>(x.s().tree.torsions.size()));
    std::vector<double> a(smp.size());
    int cis = 0, offGrid = 0;
    double sumCos = 0;
    const int n = 4000;
    std::set<long> deltas;
    for (int k = 0; k < n; ++k) {
      od::Rng rng(7, 1, k);
      smp.sample(rng, a);
      for (int i = 0; i < smp.size(); ++i) {
        CHECK((a[i] > -180.0 && a[i] <= 180.0));
        if (i == omega) continue;
        const double d = od::wrapDeg(a[i] - x.s().tree.torsions[i].inputAngle);
        const double q = d / step;
        if (std::abs(q - std::round(q)) > 1e-9 && std::abs(std::abs(d) - 180.0) > 1e-9) ++offGrid;
        if (i == 0) {
          deltas.insert(std::lround(q));
          sumCos += std::cos(d * std::numbers::pi / 180.0);
        }
      }
      CHECK((a[omega] == 0.0 || a[omega] == 180.0));
      cis += a[omega] == 0.0;
    }
    CHECK(offGrid == 0);
    CHECK(cis == doctest::Approx(n / 2).epsilon(0.08));
    if (step == 100.0) CHECK(deltas.size() == 3);  // -100, 0, 100
    else CHECK(std::abs(sumCos / n) < 0.05);
  }
}

TEST_CASE("csgm1: csgm1b types") {
  Csgm1Fixture x;
  auto& s = x.s();
  s.csgm1 = od::Csgm1::B;
  od::SsType custom{"A", false, {{-58, -48, 0.6}, {-120, 130, 0.3}}};
  x.cfg->ssTypes.push_back(custom);
  x.cfg->ssHalfWidth = 10;
  x.cfg->ssOriginalProb = 0.8;
  s.ssType[x.unitOf(3)] = "O";
  s.ssType[x.unitOf(4)] = "A";
  s.ssType[x.unitOf(5)] = "RAMA_GENERAL";
  s.ssType[x.unitOf(10)] = "RAMA_GLY";
  s.ssType[x.unitOf(12)] = "RAMA_AUTO";  // Pro after the cisTrans omega
  s.ssType[x.unitOf(16)] = "RAMA_AUTO";  // Arg before Pro 17: pre-Pro
  s.ssType[x.unitOf(6)] = "R";
  const od::Csgm1Sampler smp(x.m(), s, *x.cfg, {});
  auto idx = [&](int seq, od::BackboneAngle a) {
    const int i = x.treeIndex(seq, a);
    REQUIRE(i >= 0);
    return i;
  };
  const int phi3 = idx(3, od::BackboneAngle::Phi), psi3 = idx(3, od::BackboneAngle::Psi);
  const int phi4 = idx(4, od::BackboneAngle::Phi), psi4 = idx(4, od::BackboneAngle::Psi);
  const int phi5 = idx(5, od::BackboneAngle::Phi), psi5 = idx(5, od::BackboneAngle::Psi);
  const int phi10 = idx(10, od::BackboneAngle::Phi), psi10 = idx(10, od::BackboneAngle::Psi);
  const int phi12 = idx(12, od::BackboneAngle::Phi), psi12 = idx(12, od::BackboneAngle::Psi);
  const int phi16 = idx(16, od::BackboneAngle::Phi), psi16 = idx(16, od::BackboneAngle::Psi);
  const int omega = idx(11, od::BackboneAngle::Omega), phi6 = idx(6, od::BackboneAngle::Phi);
  const auto& tt = s.tree.torsions;
  auto in = [](double a, double c, double hw) { return std::abs(od::wrapDeg(a - c)) <= hw + 1e-9; };
  auto rs = [](od::RamaClass c) { return od::RamaSampler(od::ramaGrid(c), 1e-4); };
  const auto general = rs(od::RamaClass::General), gly = rs(od::RamaClass::Gly), prePro = rs(od::RamaClass::PrePro);
  const auto transPro = rs(od::RamaClass::TransPro), cisPro = rs(od::RamaClass::CisPro);

  const int n = 20000;
  int oIn = 0, p1 = 0, p2 = 0, outside = 0, phi6Grid = 0;
  double ll[6] = {};  // general/Gly at 5, general/Gly at 10, general/pre-Pro at 16
  double llPro[2][2] = {};  // [cis][trans-Pro, cis-Pro]
  int count[2] = {0, 0};
  std::vector<double> a(smp.size());
  for (int k = 0; k < n; ++k) {
    od::Rng rng(11, 2, k);
    smp.sample(rng, a);
    oIn += in(a[phi3], tt[phi3].inputAngle, 10) && in(a[psi3], tt[psi3].inputAngle, 10);
    const bool i1 = in(a[phi4], -58, 10) && in(a[psi4], -48, 10);
    const bool i2 = in(a[phi4], -120, 10) && in(a[psi4], 130, 10);
    p1 += i1;
    p2 += i2;
    outside += !i1 && !i2;
    ll[0] += std::log(general.probability(a[phi5], a[psi5]));
    ll[1] += std::log(gly.probability(a[phi5], a[psi5]));
    ll[2] += std::log(general.probability(a[phi10], a[psi10]));
    ll[3] += std::log(gly.probability(a[phi10], a[psi10]));
    ll[4] += std::log(general.probability(a[phi16], a[psi16]));
    ll[5] += std::log(prePro.probability(a[phi16], a[psi16]));
    const int c = a[omega] == 0.0;
    ++count[c];
    llPro[c][0] += std::log(transPro.probability(a[phi12], a[psi12]));
    llPro[c][1] += std::log(cisPro.probability(a[phi12], a[psi12]));
    const double d6 = od::wrapDeg(a[phi6] - tt[phi6].inputAngle);
    phi6Grid += std::abs(d6 - std::round(d6)) < 1e-9;
  }
  // O: 0.8 within the box, plus the uniform part falling into it.
  CHECK(oIn / double(n) == doctest::Approx(0.8 + 0.2 * (21.0 / 360) * (21.0 / 360)).epsilon(0.03));
  CHECK(p1 / double(n) == doctest::Approx(0.6).epsilon(0.03));
  CHECK(p2 / double(n) == doctest::Approx(0.3).epsilon(0.05));
  CHECK(outside / double(n) == doctest::Approx(0.1).epsilon(0.1));
  CHECK(ll[0] > ll[1]);
  CHECK(ll[3] > ll[2]);
  CHECK(ll[5] > ll[4]);
  CHECK(count[0] > n / 3);
  CHECK(count[1] > n / 3);
  CHECK(llPro[0][0] > llPro[0][1]);  // trans omega: trans-Pro
  CHECK(llPro[1][1] > llPro[1][0]);  // cis omega: cis-Pro
  CHECK(phi6Grid == n);              // R: csgm1a
}

TEST_CASE("minimize: FIRE") {
  // E = sum k_i (x_i - c_i)^2 with x in radians.
  const std::vector<double> c{30, -150, 170}, k{1, 10, 100};
  constexpr double rad = std::numbers::pi / 180.0;
  std::vector<std::vector<double>> visited;
  const od::TorsionObjective f = [&](std::span<const double> x, std::vector<double>& g) {
    visited.emplace_back(x.begin(), x.end());
    g.assign(x.size(), 0.0);
    double e = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      const double d = od::wrapDeg(x[i] - c[i]) * rad;
      e += k[i] * d * d;
      g[i] = 2 * k[i] * d;
    }
    return e;
  };
  std::vector<double> x{-100, 100, -170};
  od::MinParams p{1000, 1e-12, 10.0};
  const auto r = od::minimizeFire(f, x, p);
  for (std::size_t i = 0; i < x.size(); ++i) CHECK(std::abs(od::wrapDeg(x[i] - c[i])) < 0.01);
  CHECK(r.energy < 1e-6);
  CHECK(r.evaluations == r.steps + 1);
  double big = 0;
  for (std::size_t s = 1; s < visited.size(); ++s)
    for (std::size_t i = 0; i < x.size(); ++i) big = std::max(big, std::abs(visited[s][i] - visited[s - 1][i]));
  CHECK(big <= 10.0 + 1e-9);
  CHECK(big > 9.0);

  // Step limit and energy-change stop.
  x = {-100, 100, -170};
  CHECK(od::minimizeFire(f, x, {7, 0.0, 10.0}).steps == 7);
  x = {-100, 100, -170};
  const auto r2 = od::minimizeFire(f, x, {1000, 1.0, 10.0});
  CHECK(r2.steps < 1000);
  x = {-100, 100, -170};
  const auto r0 = od::minimizeFire(f, x, {0, 1.0, 10.0});
  CHECK(r0.steps == 0);
  CHECK(x[0] == -100);
}

TEST_CASE("9a: L2 pool refill (1L2Y)") {
  Csgm1Fixture x;
  const od::Molecule& m = x.m();
  const od::EnergySystem sys({{&m, &x.s().tree}}, {});
  const od::Csgm1Sampler smp(m, x.s(), *x.cfg, {});
  const od::ConfGenerator gen = [&](od::Rng& r, std::span<double> a) {
    smp.sample(r, a);
    return true;
  };
  od::Fill9aParams p;
  p.minL2 = 24;
  p.maxL1Batch = 10;
  p.seed = 5;
  p.stream = 1;
  auto run = [&](int threads, std::vector<od::ConfSample>& pool) {
    const int before = omp_get_max_threads();
    omp_set_num_threads(threads);
    std::uint64_t counter = 0;
    const auto st = od::fillL2Pool(sys, gen, 6, p, counter, pool);
    omp_set_num_threads(before);
    CHECK(counter == static_cast<std::uint64_t>(st.l1));
    return st;
  };
  std::vector<od::ConfSample> a, b;
  const auto st = run(1, a);
  run(4, b);
  CHECK(static_cast<int>(a.size()) >= p.minL2);
  CHECK(st.added == static_cast<long long>(a.size()));
  CHECK(st.l1 == st.added + st.overEmax + st.notClosed);
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    CHECK(a[i].angles == b[i].angles);
    CHECK(a[i].energy == b[i].energy);
  }
  for (const auto& s : a) {
    const std::vector<double> ang(s.angles.begin(), s.angles.end());
    const auto e = sys.evaluate(ang);
    CHECK(e.filterEnergy() < p.Emax);
    CHECK(e.filterEnergy() == doctest::Approx(s.energy).epsilon(1e-4));
    CHECK(e.closureMaxDeviation <= p.maxBreakDist + 1e-4);
  }
  // Minimization lowers the energy of a random sample.
  std::vector<double> ang(sys.torsionCount());
  od::Rng rng(5, 9, 0);
  smp.sample(rng, ang);
  const double e0 = sys.evaluate(ang).total();
  od::minimizeConf(sys, ang, {});
  CHECK(sys.evaluate(ang).total() <= e0);

  // No sample can pass: stops at maxL1Total.
  p.Emax = -1e9;
  p.maxL1Total = 30;
  std::vector<od::ConfSample> c;
  std::uint64_t counter = 0;
  CHECK_THROWS_AS(od::fillL2Pool(sys, gen, 6, p, counter, c), od::Error);
  CHECK(counter == 30);
}

namespace {

double maxCoordDiff(const std::vector<od::Vec3>& a, const std::vector<od::Vec3>& b) {
  REQUIRE(a.size() == b.size());
  double d = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    for (int k = 0; k < 3; ++k) d = std::max(d, std::abs(a[i][k] - b[i][k]));
  return d;
}

}  // namespace

TEST_CASE("csgm1: segments") {
  Csgm1Fixture x;
  auto& s = x.s();
  s.csgm1 = od::Csgm1::B;
  s.ssType[x.unitOf(5)] = "RAMA_GENERAL";
  s.ssType[x.unitOf(12)] = "RAMA_AUTO";
  const od::Csgm1Sampler smp(x.m(), s, *x.cfg, {});
  const auto& segs = smp.segments();
  std::vector<int> seen(smp.size(), 0);
  for (std::size_t k = 0; k < segs.size(); ++k) {
    CHECK(std::is_sorted(segs[k].begin(), segs[k].end()));
    if (k > 0) CHECK(segs[k - 1].front() < segs[k].front());
    for (int i : segs[k]) ++seen[i];
  }
  CHECK(std::count(seen.begin(), seen.end(), 1) == smp.size());
  auto segOf = [&](int i) {
    for (std::size_t k = 0; k < segs.size(); ++k)
      if (std::find(segs[k].begin(), segs[k].end(), i) != segs[k].end()) return static_cast<int>(k);
    return -1;
  };
  const int phi5 = x.treeIndex(5, od::BackboneAngle::Phi), psi5 = x.treeIndex(5, od::BackboneAngle::Psi);
  const int phi12 = x.treeIndex(12, od::BackboneAngle::Phi), psi12 = x.treeIndex(12, od::BackboneAngle::Psi);
  const int omega = x.treeIndex(11, od::BackboneAngle::Omega);
  CHECK(segs[segOf(phi5)].size() == 2);
  CHECK(segOf(phi5) == segOf(psi5));
  CHECK(segs[segOf(phi12)].size() == 3);
  CHECK(segOf(phi12) == segOf(omega));
  CHECK(segOf(psi12) == segOf(omega));

  // sampleSegment changes only the angles of the segment.
  std::vector<double> a(smp.size());
  od::Rng rng(3, 4, 5);
  smp.sample(rng, a);
  for (int seg : {segOf(phi12), segOf(phi5), 0}) {
    const auto before = a;
    for (int k = 0; k < 20; ++k) smp.sampleSegment(rng, seg, a);
    for (int i = 0; i < smp.size(); ++i) {
      const bool member = std::find(segs[seg].begin(), segs[seg].end(), i) != segs[seg].end();
      if (!member) CHECK(a[i] == before[i]);
      else CHECK(a[i] != before[i]);
    }
  }
  od::Csgm1GroupSampler g;
  g.add(smp);
  g.add(smp);
  const auto gs = g.segments();
  CHECK(gs.size() == 2 * segs.size());
  CHECK(gs[segs.size()].front() == smp.size());
}

TEST_CASE("grow: coordinates, clashes and retries (1L2Y)") {
  Csgm1Fixture x;
  const od::Molecule& m = x.m();
  const std::vector<od::EnergyPart> parts{{&m, &x.s().tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(m, x.s(), *x.cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };

  // Without checks every sample is built, with the coordinates of buildCoordinates.
  const od::GrowBuilder open(sys, parts, {&x.s().topology}, smp.segments(), {0, 0.0, 1e9, 1e9});
  CHECK(open.stepCount() == static_cast<int>(smp.segments().size()));
  std::vector<double> a(smp.size());
  for (int k = 0; k < 50; ++k) {
    od::Rng rng(21, 0, k);
    smp.sample(rng, a);
    std::vector<od::Vec3> grown, ref;
    REQUIRE(open.build(rng, a, resample, &grown).ok);
    sys.evaluate(a, nullptr, &ref);
    CHECK(maxCoordDiff(grown, ref) < 1e-9);
  }

  // Default checks: built samples have no clash; retries resample whole segments.
  const od::GrowBuilder grow(sys, parts, {&x.s().topology}, smp.segments(), {});
  int ok = 0, retried = 0;
  for (int k = 0; k < 300; ++k) {
    od::Rng rng(22, 0, k);
    smp.sample(rng, a);
    const auto r = grow.build(rng, a, resample);
    CHECK(r.reachFails == 0);
    retried += r.retries > 0;
    if (!r.ok) {
      CHECK(r.clashFails == r.retries + 1);
      CHECK(r.retries >= 16);  // the last segment failed its first try and 16 retries
      continue;
    }
    ++ok;
    CHECK(r.clashFails == r.retries);
    CHECK(sys.evaluate(a).clashPairs == 0);
  }
  CHECK(ok > 0);
  CHECK(retried > 0);
  // No retries: a sample is either built at once or discarded.
  const od::GrowBuilder noRetry(sys, parts, {&x.s().topology}, smp.segments(), {0});
  for (int k = 0; k < 100; ++k) {
    od::Rng rng(22, 0, k);
    smp.sample(rng, a);
    const auto r = noRetry.build(rng, a, resample);
    CHECK(r.retries == 0);
    CHECK(r.clashFails == (r.ok ? 0 : 1));
  }
}

TEST_CASE("grow: reach toward a loop closure (2KJJ)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 10, 14));
  const od::MolSettings& s = cfg.molecules[0];
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &s.tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  const od::GrowBuilder reachOnly(sys, parts, {&s.topology}, smp.segments(), {0, 0.0});
  CHECK(reachOnly.reachCheckCount() > 0);

  // The input conformation and closed minimized conformations pass the reach checks.
  std::vector<double> a = sys.inputAngles();
  od::Rng rng(1, 0, 0);
  CHECK(reachOnly.build(rng, a, resample).ok);
  int closed = 0;
  for (int k = 0; k < 10; ++k) {
    auto x = perturbed(sys.inputAngles(), 30, 100 + k);
    const auto e = od::minimizeConf(sys, x, {});
    if (e.closureMaxDeviation > 0.5) continue;
    ++closed;
    CHECK(reachOnly.build(rng, x, resample).ok);
  }
  CHECK(closed > 0);

  // Random loops: some fail the reach checks; a tiny reach rejects every sample.
  int reachFails = 0;
  for (int k = 0; k < 200; ++k) {
    od::Rng r(3, 0, k);
    smp.sample(r, a);
    reachFails += reachOnly.build(r, a, resample).reachFails > 0;
  }
  CHECK(reachFails > 0);
  const od::GrowBuilder tiny(sys, parts, {&s.topology}, smp.segments(), {0, 0.0, 0.01, 0.01});
  for (int k = 0; k < 20; ++k) {
    od::Rng r(4, 0, k);
    smp.sample(r, a);
    const auto res = tiny.build(r, a, resample);
    CHECK_FALSE(res.ok);
    CHECK(res.reachFails == 1);
  }
}

namespace {

// StaM 2KJJ in its input conformation and MobM 1L2Y as grid atoms.
struct GridFixture {
  StaFixture f;
  od::MolConfig cfg;
  std::vector<od::EnergyPart> staParts, mobParts;
  std::unique_ptr<od::EnergySystem> staSys, mobSys;
  od::GridAtoms stam, mobm;

  GridFixture() : cfg(f.ok(f.text)) {
    staParts = {{&f.mols[0], &cfg.molecules[0].tree}};
    mobParts = {{&f.mols[1], &cfg.molecules[1].tree}};
    staSys = std::make_unique<od::EnergySystem>(staParts, od::EnergyParams{});
    mobSys = std::make_unique<od::EnergySystem>(mobParts, od::EnergyParams{});
    std::vector<od::Vec3> x;
    staSys->evaluate(staSys->inputAngles(), nullptr, &x);
    stam = od::gridAtoms(*staSys, staParts, x);
    mobSys->evaluate(mobSys->inputAngles(), nullptr, &x);
    mobm = od::gridAtoms(*mobSys, mobParts, x);
  }
  od::Vec3 stamCenter() const {
    od::Vec3 c{};
    for (const auto& p : stam.xyz)
      for (int k = 0; k < 3; ++k) c[k] += p[k] / static_cast<double>(stam.size());
    return c;
  }
};

}  // namespace

TEST_CASE("grid: geometry") {
  const auto g = od::boxGrid({1, 2, 3}, {10, 7.4, 0.375}, 0.375);
  CHECK(g.n[0] == 29);  // 26.67 -> 27 -> 28 intervals
  CHECK(g.n[1] == 21);  // 19.73 -> 20
  CHECK(g.n[2] == 3);   // 1 -> 2
  for (int c = 0; c < 3; ++c) CHECK(g.origin[c] + 0.5 * (g.n[c] - 1) * g.spacing == doctest::Approx(c + 1.0));
  CHECK(g.inside(g.origin));
  CHECK(g.inside(g.max()));
  CHECK_FALSE(g.inside({g.origin[0] - 1e-6, 2, 3}));
  // Trilinear interpolation reproduces a linear function and its gradient.
  std::vector<float> map(g.size());
  for (int k = 0; k < g.n[2]; ++k)
    for (int j = 0; j < g.n[1]; ++j)
      for (int i = 0; i < g.n[0]; ++i) {
        const auto p = g.point(i, j, k);
        map[g.index(i, j, k)] = static_cast<float>(2 * p[0] - 3 * p[1] + 0.5 * p[2]);
      }
  od::Vec3 grad;
  const od::Vec3 p{2.1, 1.3, 3.2};
  CHECK(od::interpolate(g, map, p, &grad) == doctest::Approx(2 * 2.1 - 3 * 1.3 + 0.5 * 3.2).epsilon(1e-5));
  CHECK(grad[0] == doctest::Approx(2).epsilon(1e-4));
  CHECK(grad[1] == doctest::Approx(-3).epsilon(1e-4));
  CHECK(grad[2] == doctest::Approx(0.5).epsilon(1e-4));
  CHECK(od::interpolate(g, map, g.max()) == doctest::Approx(map.back()));
}

TEST_CASE("grid: energy and distance grids (2KJJ, 1L2Y)") {
  GridFixture x;
  const od::Vec3 c = x.stamCenter();
  const auto geom = od::boxGrid(c, {20, 20, 20}, 0.375);
  const auto grids = od::buildEnergyGrids(geom, x.stam, x.mobm.type, {});
  CHECK(grids.types.size() < x.mobm.size());
  for (od::Ad4Type t : x.mobm.type) CHECK(grids.mapOf(t) >= 0);

  // At grid points the grid energy equals the pair sum (points away from clashes).
  std::mt19937 rng(5);
  std::uniform_int_distribution<int> ui(0, geom.n[0] - 1);
  od::GridAtoms probe;
  const auto dist = od::buildDistanceGrid(geom, x.stam, 6.0);
  int checked = 0;
  for (int k = 0; k < 4000 && checked < 300; ++k) {
    const int i = ui(rng), j = ui(rng), l = ui(rng);
    const std::size_t idx = geom.index(i, j, l);
    if (dist.d[idx] < 2.5) continue;
    od::GridAtoms one;
    const std::size_t a = static_cast<std::size_t>(k) % x.mobm.size();
    one.xyz = {geom.point(i, j, l)};
    one.type = {x.mobm.type[a]};
    one.charge = {x.mobm.charge[a]};
    one.heavy = {x.mobm.heavy[a]};
    const double eg = od::gridEnergy(grids, one).energy;
    const double ep = od::pairInterEnergy(x.stam, one, {});
    CHECK(eg == doctest::Approx(ep).epsilon(1e-5).scale(1e-3));
    ++checked;
  }
  CHECK(checked >= 100);

  // Distance grid: brute force at grid points, capped.
  for (int k = 0; k < 500; ++k) {
    const int i = ui(rng), j = ui(rng), l = ui(rng);
    const auto p = geom.point(i, j, l);
    double best = 6.0;
    for (std::size_t a = 0; a < x.stam.size(); ++a)
      if (x.stam.heavy[a]) {
        const double dx = p[0] - x.stam.xyz[a][0], dy = p[1] - x.stam.xyz[a][1], dz = p[2] - x.stam.xyz[a][2];
        best = std::min(best, std::sqrt(dx * dx + dy * dy + dz * dz));
      }
    CHECK(dist.d[geom.index(i, j, l)] == doctest::Approx(best).epsilon(1e-5));
  }

  // Probe atoms at random points 3-6 Å from the StaM: interpolated energy
  // close to the pair sum, analytic gradient equal to finite differences.
  std::uniform_real_distribution<double> u(0, 1);
  od::GridAtoms m;
  while (m.size() < 40) {
    od::Vec3 p;
    for (int k = 0; k < 3; ++k) p[k] = geom.origin[k] + u(rng) * (geom.n[k] - 1) * geom.spacing;
    const double d = od::interpolate(geom, dist.d, p);
    if (d < 3.0 || d >= 6.0) continue;
    const std::size_t a = m.size() % x.mobm.size();
    m.xyz.push_back(p);
    m.type.push_back(x.mobm.type[a]);
    m.charge.push_back(x.mobm.charge[a]);
    m.heavy.push_back(x.mobm.heavy[a]);
  }
  std::vector<od::Vec3> g;
  const auto em = od::gridEnergy(grids, m, &g);
  CHECK(em.outside == 0);
  CHECK(em.energy == doctest::Approx(od::pairInterEnergy(x.stam, m, {})).epsilon(0.15));  // interpolation error
  for (std::size_t a = 0; a < m.size(); ++a) {
    od::GridAtoms one;
    one.xyz = {m.xyz[a]};
    one.type = {m.type[a]};
    one.charge = {m.charge[a]};
    one.heavy = {m.heavy[a]};
    CHECK(od::gridEnergy(grids, one).energy == doctest::Approx(od::pairInterEnergy(x.stam, one, {})).epsilon(0.15).scale(1.0));
    const double h = 1e-5;
    for (int k = 0; k < 3; ++k) {
      od::GridAtoms mp = m, mm = m;
      mp.xyz[a][k] += h;
      mm.xyz[a][k] -= h;
      const double fd = (od::gridEnergy(grids, mp).energy - od::gridEnergy(grids, mm).energy) / (2 * h);
      CHECK(g[a][k] == doctest::Approx(fd).epsilon(1e-3).scale(1e-2));
    }
  }
  // Atoms outside the grid are counted and skipped.
  od::GridAtoms far = x.mobm;
  for (auto& p : far.xyz) p[0] += 1000;
  const auto e = od::gridEnergy(grids, far);
  CHECK(e.outside == static_cast<int>(far.size()));
  CHECK(e.energy == 0);
}

TEST_CASE("orient: quaternions, Shoemake and Hopf grid") {
  // Rotation vector, composition and rotation angle.
  const od::Quat a = od::fromRotationVector({0, 0, std::numbers::pi / 2});
  const od::Vec3 v = od::rotate(a, {1, 0, 0});
  CHECK(v[0] == doctest::Approx(0).scale(1));
  CHECK(v[1] == doctest::Approx(1));
  const od::Quat b = od::fromRotationVector({0.3, -0.2, 0.5});
  const od::Vec3 w = od::rotate(a * b, {0.2, 1.1, -0.7});
  const od::Vec3 w2 = od::rotate(a, od::rotate(b, {0.2, 1.1, -0.7}));
  for (int k = 0; k < 3; ++k) CHECK(w[k] == doctest::Approx(w2[k]));
  CHECK(od::rotationAngle(od::Quat{}, b) == doctest::Approx(std::sqrt(0.09 + 0.04 + 0.25)));
  CHECK(od::rotationAngle(b, od::Quat{-b.w, -b.x, -b.y, -b.z}) == doctest::Approx(0).scale(1));

  CHECK(od::hopfLevelForStep(7.5) == 3);
  CHECK(od::hopfLevelForStep(10) == 3);
  CHECK(od::hopfLevelForStep(50) == 0);
  CHECK(od::hopfLevelForStep(0.1) == 6);

  // Random orientations: unit and uniform (mean of |w| on S3 is 4 / (3 pi)).
  double meanW = 0;
  const int nr = 20000;
  std::vector<od::Quat> rnd;
  for (int i = 0; i < nr; ++i) {
    od::Rng rng(7, 0, i);
    const od::Quat q = od::randomQuat(rng);
    CHECK(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z == doctest::Approx(1));
    meanW += std::abs(q.w) / nr;
    rnd.push_back(q);
  }
  CHECK(meanW == doctest::Approx(4 / (3 * std::numbers::pi)).epsilon(0.01));

  for (int level = 0; level <= 2; ++level) {
    const auto g = od::hopfGrid(level);
    REQUIRE(g.size() == static_cast<std::size_t>(72) << (3 * level));
    double minPair = 10, cover = 0;
    for (std::size_t i = 0; i < g.size(); ++i) {
      CHECK(g[i].w * g[i].w + g[i].x * g[i].x + g[i].y * g[i].y + g[i].z * g[i].z == doctest::Approx(1));
      for (std::size_t j = 0; j < i; ++j) minPair = std::min(minPair, od::rotationAngle(g[i], g[j]));
    }
    for (std::size_t r = 0; r < rnd.size(); r += 4) {
      double best = 10;
      for (const auto& q : g) best = std::min(best, od::rotationAngle(rnd[r], q));
      cover = std::max(cover, best);
    }
    const double s = od::hopfSpacingDeg(level) * std::numbers::pi / 180;
    MESSAGE("level ", level, ": min pair ", minPair * 180 / std::numbers::pi, " deg, covering ",
            cover * 180 / std::numbers::pi, " deg");
    CHECK(minPair > 0.3 * s);
    CHECK(cover < 1.0 * s);
  }
}

TEST_CASE("lbfgs: quadratic and Rosenbrock") {
  const od::VecObjective quad = [](std::span<const double> x, double& e, std::vector<double>& g) {
    const double a[3] = {1, 10, 100};
    e = 0;
    g.assign(3, 0);
    for (int i = 0; i < 3; ++i) {
      e += 0.5 * a[i] * (x[i] - i) * (x[i] - i);
      g[i] = a[i] * (x[i] - i);
    }
    return true;
  };
  od::LbfgsParams p;
  p.deltaE = 1e-12;
  p.convSteps = 2;
  std::vector<double> x{5, -3, 7};
  auto r = od::minimizeLbfgs(quad, x, p);
  CHECK_FALSE(r.aborted);
  for (int i = 0; i < 3; ++i) CHECK(x[i] == doctest::Approx(i).scale(1).epsilon(1e-5));
  const od::VecObjective rosen = [](std::span<const double> x, double& e, std::vector<double>& g) {
    e = 100 * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1 - x[0], 2);
    g = {-400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]), 200 * (x[1] - x[0] * x[0])};
    return true;
  };
  x = {-1.2, 1};
  r = od::minimizeLbfgs(rosen, x, p);
  CHECK(x[0] == doctest::Approx(1).epsilon(1e-4));
  CHECK(x[1] == doctest::Approx(1).epsilon(1e-4));
  // Step limit, accept hook and abort.
  int accepted = 0;
  p.stepLimit = [](std::span<const double> d) { return 0.1 / std::max(0.1, std::abs(d[0])); };
  p.onAccept = [&](std::vector<double>&) { ++accepted; };
  x = {5, -3, 7};
  double last = 5;
  const od::VecObjective quad2 = [&](std::span<const double> v, double& e, std::vector<double>& g) {
    CHECK(std::abs(v[0] - last) <= 0.1 + 1e-12);
    return quad(v, e, g);
  };
  p.onAccept = [&](std::vector<double>& v) { ++accepted; last = v[0]; };
  r = od::minimizeLbfgs(quad2, x, p);
  CHECK(accepted == r.steps);
  CHECK(r.steps >= 50);
  const od::VecObjective stop = [&](std::span<const double> v, double& e, std::vector<double>& g) {
    quad(v, e, g);
    return v[0] > 4;
  };
  x = {5, -3, 7};
  p.stepLimit = {};
  p.onAccept = {};
  // An infeasible trial point aborts.
  CHECK(od::minimizeLbfgs(stop, x, p).aborted);
  x = {3, -3, 7};
  CHECK(od::minimizeLbfgs(stop, x, p).aborted);
}

namespace {

struct DockFixture {
  GridFixture x;
  od::Dock9bParams p;
  od::EnergyGrids grids;
  od::DistanceGrid dist;
  DockFixture(od::Posm posm = od::Posm::Posm1) {
    p.posm = posm;
    p.boxCenter = x.stamCenter();
    p.boxSize = {48, 48, 48};
    p.hopfLevel = 0;
    p.posm1SampleNum = 400;
    p.seed = 3;
    const auto geom = od::boxGrid(p.boxCenter, p.boxSize, 0.375);
    grids = od::buildEnergyGrids(geom, x.stam, x.mobm.type, {});
    dist = od::buildDistanceGrid(geom, x.stam, std::max(p.maxNearDistance, p.clashDistance) + 0.75);
  }
};

}  // namespace

TEST_CASE("dock9b: poses, quick filter and gradient (2KJJ, 1L2Y)") {
  DockFixture f;
  const od::MobmDocker dk(*f.x.mobSys, f.grids, f.dist, f.p);
  const auto angles = f.x.mobSys->inputAngles();
  const auto c = dk.conformation(angles);
  // Heavy-atom centre.
  od::Vec3 hc{};
  int nh = 0;
  for (std::size_t a = 0; a < f.x.mobm.size(); ++a)
    if (f.x.mobm.heavy[a]) {
      for (int k = 0; k < 3; ++k) hc[k] += f.x.mobm.xyz[a][k];
      ++nh;
    }
  for (int k = 0; k < 3; ++k) CHECK(c.center[k] == doctest::Approx(hc[k] / nh));

  // Probe atoms (all real atoms) relative to the heavy-atom centre.
  std::vector<od::Vec3> rel;
  for (const auto& y : f.x.mobm.xyz) rel.push_back({y[0] - c.center[0], y[1] - c.center[1], y[2] - c.center[2]});
  // shrink: distance kept from the box faces.
  const auto allInside = [&](const od::Pose& pose, const od::Vec3& bc, double shrink = 0) {
    for (const auto& y : rel) {
      const auto r = od::rotate(pose.q, y);
      for (int d = 0; d < 3; ++d)
        if (std::abs(r[d] + pose.t[d] - bc[d]) > 24 - shrink + 1e-9) return false;
    }
    return true;
  };
  const double inner = od::MobmDocker::kPosMargin - 1e-4;  // float rounding of the range
  const auto onLattice = [&](const od::Pose& pose, const od::Vec3& bc) {
    for (int d = 0; d < 3; ++d) {
      const double u = pose.t[d] - (bc[d] - 24);
      if (std::abs(u - std::round(u)) > 1e-9 || u < -1e-9 || u > 48 + 1e-9) return false;
    }
    return true;
  };

  // posm1: lattice positions with every probe atom in the box, reproducible.
  CHECK(dk.l1Count(c) == 400);
  for (long long k = 0; k < 400; ++k) {
    const auto a = dk.l1Pose(c, 9, k), b = dk.l1Pose(c, 9, k);
    REQUIRE(a.has_value());  // 1L2Y fits in a 48 Å box in every orientation
    CHECK(a->q.w == b->q.w);
    CHECK(a->t == b->t);
    CHECK(onLattice(*a, f.p.boxCenter));
    CHECK(allInside(*a, f.p.boxCenter, inner));
  }
  CHECK(dk.l1Pose(c, 9, 0)->t != dk.l1Pose(c, 10, 0)->t);
  // The range is exact (up to kPosMargin): one more lattice step comes closer than kPosMargin to a box face.
  for (long long k = 0; k < 50; ++k) {
    const auto q = dk.l1Pose(c, 9, k)->q;
    const auto pr = dk.posRange(c, q);
    for (int d = 0; d < 3; ++d) {
      od::Pose lo{f.p.boxCenter, q}, hi{f.p.boxCenter, q};
      lo.t[d] = f.p.boxCenter[d] - 24 + pr.lo[d] - 1;
      hi.t[d] = f.p.boxCenter[d] - 24 + pr.hi[d] + 1;
      od::Pose in = lo;
      in.t[d] += 1;
      CHECK(allInside(in, f.p.boxCenter, inner) == (pr.lo[d] <= pr.hi[d]));
      if (pr.lo[d] > 0) CHECK_FALSE(allInside(lo, f.p.boxCenter, od::MobmDocker::kPosMargin + 1e-4));
      if (pr.hi[d] < 48) CHECK_FALSE(allInside(hi, f.p.boxCenter, od::MobmDocker::kPosMargin + 1e-4));
    }
  }
  // posm2: per Hopf orientation, exactly the lattice positions keeping the MobM in the box.
  DockFixture f2(od::Posm::Posm2);
  const od::MobmDocker dk2(*f2.x.mobSys, f2.grids, f2.dist, f2.p);
  const auto c2 = dk2.conformation(angles);
  const auto hopf = od::hopfGrid(0);
  REQUIRE(c2.hopfRange.size() == 72);
  long long total = 0;
  for (std::size_t j = 0; j < 72; ++j) {
    const auto& pr = c2.hopfRange[j];
    CHECK(pr.count() == (pr.hi[0] - pr.lo[0] + 1) * (pr.hi[1] - pr.lo[1] + 1) * (pr.hi[2] - pr.lo[2] + 1));
    CHECK(c2.hopfOffset[j] == total);
    total += pr.count();
  }
  CHECK(dk2.l1Count(c2) == total);
  CHECK(total < 49LL * 49 * 49 * 72);
  for (long long k = 0; k < total; k += 997) {
    const auto pose = dk2.l1Pose(c2, 0, k);
    REQUIRE(pose.has_value());
    CHECK(onLattice(*pose, f2.p.boxCenter));
    CHECK(allInside(*pose, f2.p.boxCenter, inner));
  }
  const auto first1 = dk2.l1Pose(c2, 0, c2.hopfOffset[1]);
  CHECK(first1->q.x == hopf[1].x);
  for (int d = 0; d < 3; ++d) CHECK(first1->t[d] == doctest::Approx(f2.p.boxCenter[d] - 24 + c2.hopfRange[1].lo[d]));
  const auto last0 = dk2.l1Pose(c2, 0, c2.hopfOffset[1] - 1);
  CHECK(last0->q.x == hopf[0].x);
  for (int d = 0; d < 3; ++d) CHECK(last0->t[d] == doctest::Approx(f2.p.boxCenter[d] - 24 + c2.hopfRange[0].hi[d]));

  // Quick filter against brute force.
  od::Dock9bStats st;
  int pass = 0;
  for (long long k = 0; k < 400; ++k) {
    const auto pose = *dk.l1Pose(c, 1, k);
    int overlap = 0;
    double dmin = 1e9;
    for (std::size_t a = 0; a < f.x.mobm.size(); ++a) {
      if (!f.x.mobm.heavy[a]) continue;
      const auto r = od::rotate(pose.q, rel[a]);
      const od::Vec3 w{r[0] + pose.t[0], r[1] + pose.t[1], r[2] + pose.t[2]};
      double da = 1e9;
      for (std::size_t b = 0; b < f.x.stam.size(); ++b)
        if (f.x.stam.heavy[b]) {
          const double dx = w[0] - f.x.stam.xyz[b][0], dy = w[1] - f.x.stam.xyz[b][1], dz = w[2] - f.x.stam.xyz[b][2];
          da = std::min(da, std::sqrt(dx * dx + dy * dy + dz * dz));
        }
      dmin = std::min(dmin, da);
      if (da < 1.7) ++overlap;  // margin for interpolation of the distance grid
    }
    const bool ok = dk.quickFilter(c, pose, st);
    pass += ok;
    if (ok) {
      CHECK(overlap == 0);
      CHECK(dmin < 6.3);
    }
  }
  CHECK(st.outsideBox == 0);
  CHECK(st.overlap + st.far + pass == 400);
  MESSAGE("quick filter: pass ", pass, ", overlap ", st.overlap, ", far ", st.far);
  CHECK(pass > 0);

  // Gradient against finite differences at passing poses (double reference),
  // and the single-precision evaluation against the double reference.
  int checked = 0;
  for (long long k = 0; k < 400 && checked < 3; ++k) {
    const auto pose = *dk.l1Pose(c, 1, k);
    od::Dock9bStats s;
    if (!dk.quickFilter(c, pose, s)) continue;
    ++checked;
    std::vector<double> g, gf;
    const auto pe = dk.evaluateDouble(c, pose, angles, &g);
    CHECK(pe.outside == 0);
    const auto pf = dk.evaluate(c, pose, angles, &gf);
    CHECK(pf.outside == 0);
    CHECK(pf.objective() == doctest::Approx(pe.objective()).epsilon(1e-4).scale(1));
    REQUIRE(gf.size() == g.size());
    // Relative to the largest component: float cancellation in sums of large pair terms.
    double gmax = 1;
    for (double v : g) gmax = std::max(gmax, std::abs(v));
    for (std::size_t i = 0; i < g.size(); ++i) CHECK(std::abs(gf[i] - g[i]) <= 1e-4 * gmax);
    const double h = 1e-6;
    for (int i = 0; i < 3; ++i) {
      od::Pose a = pose, b = pose;
      a.t[i] += h;
      b.t[i] -= h;
      const double fd = (dk.evaluateDouble(c, a, angles).objective() - dk.evaluateDouble(c, b, angles).objective()) / (2 * h);
      CHECK(g[i] == doctest::Approx(fd).epsilon(1e-3).scale(1e-2));
      od::Vec3 v{};
      v[i] = h;
      a = pose;
      b = pose;
      a.q = od::fromRotationVector(v) * pose.q;
      v[i] = -h;
      b.q = od::fromRotationVector(v) * pose.q;
      const double fr = (dk.evaluateDouble(c, a, angles).objective() - dk.evaluateDouble(c, b, angles).objective()) / (2 * h);
      CHECK(g[3 + i] == doctest::Approx(fr).epsilon(1e-3).scale(1e-2));
    }
    for (std::size_t t = 0; t < angles.size(); t += 5) {
      auto ap = angles, am = angles;
      const double hd = h * 180 / std::numbers::pi;
      ap[t] += hd;
      am[t] -= hd;
      const double fd = (dk.evaluateDouble(c, pose, ap).objective() - dk.evaluateDouble(c, pose, am).objective()) / (2 * h);
      CHECK(g[6 + t] == doctest::Approx(fd).epsilon(1e-3).scale(1e-2));
    }
  }
  CHECK(checked == 3);
}

TEST_CASE("dock9b: docking one conformation pair (2KJJ, 1L2Y)") {
  DockFixture f;
  f.p.posm1SampleNum = 1500;
  f.p.lbfgs.maxSteps = 300;
  const od::MobmDocker dk(*f.x.mobSys, f.grids, f.dist, f.p);
  const auto angles = f.x.mobSys->inputAngles();
  const auto c = dk.conformation(angles);
  std::vector<od::DockSample> out;
  const auto st = dk.dockPair(angles, 5, out);
  CHECK(st.l1 == 1500);
  CHECK(st.l1 == st.l2 + st.overlap + st.far + st.outsideBox);
  CHECK(static_cast<long long>(out.size()) == st.l2 - st.leftGrid);
  REQUIRE(!out.empty());
  float best = 1e30f;
  for (const auto& s : out) {
    best = std::min(best, s.energy);
    CHECK(s.energy == doctest::Approx(s.inter + s.interClash + s.intra).epsilon(1e-4));
    CHECK(s.steps <= 300);
    // Re-evaluating the stored float pose and angles reproduces the stored energy exactly.
    for (double v : s.pose.t) CHECK(v == static_cast<double>(static_cast<float>(v)));
    std::vector<double> a(s.angles.begin(), s.angles.end());
    const auto pe = dk.evaluate(c, s.pose, a);
    CHECK(pe.outside == 0);
    CHECK(static_cast<float>(pe.ranking()) == s.energy);
  }
  // Minimization lowers the energy of every kept pose.
  std::vector<long long> l2;
  for (long long k = 0; k < 1500; ++k) {
    od::Dock9bStats s;
    const auto pose = dk.l1Pose(c, 5, k);
    if (pose && dk.quickFilter(c, *pose, s)) l2.push_back(k);
  }
  std::size_t j = 0;
  for (long long k : l2) {
    od::DockSample s;
    od::Dock9bStats s2;
    const auto start = *dk.l1Pose(c, 5, k);
    if (!dk.minimize(c, start, s, s2)) continue;
    REQUIRE(j < out.size());
    CHECK(s.energy == out[j].energy);  // same result as the parallel run, same order
    CHECK(s.energy <= dk.evaluate(c, start, angles).ranking() + 1e-3);
    ++j;
  }
  CHECK(j == out.size());
  MESSAGE("9b: L1 ", st.l1, " L2 ", st.l2, " left grid ", st.leftGrid, " steps/pose ",
          static_cast<double>(st.steps) / std::max<long long>(1, st.l2), " best ", best);
}

TEST_CASE("dock9c: RMSD, topSampleSet and dcss filter (2KJJ, 1L2Y)") {
  const std::vector<od::Vec3f> a{{0, 0, 0}, {1, 0, 0}, {0, 2, 0}};
  CHECK(od::fixedRmsd(a, a) == 0);
  const std::vector<od::Vec3f> sw{{1, 0, 0}, {0, 0, 0}, {0, 2, 0}};  // atoms 0 and 1 swapped
  CHECK(od::fixedRmsd(a, sw) == doctest::Approx(std::sqrt(2.0 / 3)));
  std::vector<od::Vec3f> sh = a;
  for (auto& p : sh) p[2] += 3;
  CHECK(od::fixedRmsd(a, sh) == doctest::Approx(3));

  // Docking samples of one conformation pair.
  DockFixture f;
  f.p.posm1SampleNum = 600;
  f.p.lbfgs.maxSteps = 100;
  const od::MobmDocker dk(*f.x.mobSys, f.grids, f.dist, f.p);
  const auto angles = f.x.mobSys->inputAngles();
  std::vector<od::DockSample> out;
  dk.dockPair(angles, 9, out);
  REQUIRE(out.size() >= 10);
  const od::MobmHeavyAtoms heavy(*f.x.mobSys);
  CHECK(heavy.size() == static_cast<std::size_t>(std::count_if(f.x.mobm.heavy.begin(), f.x.mobm.heavy.end(), [](char h) { return h != 0; })));
  // World coordinates: inside the docking box, and the energy of the sample re-evaluates from the same pose.
  for (const auto& p : heavy.world(out[0]))
    for (int d = 0; d < 3; ++d) CHECK(std::abs(p[d] - f.p.boxCenter[d]) <= 0.5 * f.p.boxSize[d] + 1e-3);
  CHECK(od::fixedRmsd(heavy.world(out[0]), heavy.world(out[0])) == 0);

  std::vector<od::RoundSample> merged;
  for (std::size_t k = 0; k < out.size(); ++k) merged.push_back({static_cast<int>(k % 2), 0, out[k]});
  std::vector<float> e;
  for (const auto& m : merged) e.push_back(m.s.energy);
  std::sort(e.begin(), e.end());

  // Without the diversity filter: the topNum lowest energies.
  auto top = od::selectTopSamples(merged, heavy, 5, 0);
  REQUIRE(top.size() == 5);
  for (int k = 0; k < 5; ++k) CHECK(top[k].s.energy == e[k]);
  CHECK(od::topMaxEnergy(top) == e[4]);
  {
    // Ranking energy includes the StaM term (spec 9): StaM 1 shifted below StaM 0.
    auto shifted = merged;
    for (auto& m : shifted) m.stamRel = m.stam == 1 ? -1000.0f : 0.0f;
    const auto t2 = od::selectTopSamples(shifted, heavy, static_cast<int>(shifted.size()), 0);
    for (std::size_t k = 1; k < t2.size(); ++k) CHECK(t2[k - 1].rank() <= t2[k].rank());
    CHECK(t2.front().stam == 1);
    CHECK(t2.back().stam == 0);
    CHECK(od::topMaxEnergy(t2) == t2.back().rank());
  }

  // With the filter: sorted, and samples of the same StaM conformation at least diversityRmsd apart.
  // Threshold: median fixed RMSD between samples with the same StaM id.
  std::vector<double> pr;
  for (std::size_t k = 0; k < merged.size(); ++k)
    for (std::size_t j = 0; j < k; ++j)
      if (merged[j].stam == merged[k].stam) pr.push_back(od::fixedRmsd(heavy.world(merged[k].s), heavy.world(merged[j].s)));
  std::nth_element(pr.begin(), pr.begin() + pr.size() / 2, pr.end());
  const double rmsd = pr[pr.size() / 2];
  top = od::selectTopSamples(merged, heavy, 1000, rmsd);
  REQUIRE(!top.empty());
  CHECK(top.size() < merged.size());
  CHECK(top[0].s.energy == e[0]);
  for (std::size_t k = 1; k < top.size(); ++k) CHECK(top[k - 1].s.energy <= top[k].s.energy);
  for (std::size_t k = 0; k < top.size(); ++k)
    for (std::size_t j = 0; j < k; ++j)
      if (top[j].stam == top[k].stam) CHECK(od::fixedRmsd(heavy.world(top[k].s), heavy.world(top[j].s)) >= rmsd);
  // Every skipped sample is close to a kept sample of the same StaM conformation with a lower or equal energy.
  for (const auto& m : merged) {
    bool kept = false, covered = false;
    for (const auto& t : top) {
      kept |= t.s.energy == m.s.energy && t.stam == m.stam && t.s.pose.t == m.s.pose.t;
      covered |= t.stam == m.stam && t.s.energy <= m.s.energy &&
                 od::fixedRmsd(heavy.world(m.s), heavy.world(t.s)) < rmsd;
    }
    CHECK((kept || covered));
  }
  // Ties keep the merged order (previous topSampleSet first).
  std::vector<od::RoundSample> tie{{0, 1, out[0]}, {0, 2, out[0]}};
  top = od::selectTopSamples(tie, heavy, 2, 0);
  CHECK(top[0].mobm == 1);
  CHECK(top[1].mobm == 2);
  top = od::selectTopSamples(tie, heavy, 2, rmsd);
  REQUIRE(top.size() == 1);
  CHECK(top[0].mobm == 1);

  // Chunked dockPair: the same docking samples in the same order.
  std::vector<od::DockSample> chunked;
  std::vector<std::size_t> sizes;
  const auto stc = dk.dockPair(angles, 9, 7, [&](std::vector<od::DockSample>& c) {
    sizes.push_back(c.size());
    chunked.insert(chunked.end(), c.begin(), c.end());
  });
  REQUIRE(chunked.size() == out.size());
  CHECK(sizes.size() == static_cast<std::size_t>((stc.l2 + 6) / 7));
  for (std::size_t k = 0; k < out.size(); ++k) {
    CHECK(chunked[k].energy == out[k].energy);
    CHECK(chunked[k].pose.t == out[k].pose.t);
  }
  // 9c in batches without the diversity filter: the same topSampleSet as at once.
  std::vector<od::RoundSample> inc;
  for (std::size_t b = 0; b < merged.size(); b += 20) {
    std::vector<od::RoundSample> m2 = inc;
    m2.insert(m2.end(), merged.begin() + b, merged.begin() + std::min(merged.size(), b + 20));
    inc = od::selectTopSamples(m2, heavy, 5, 0);
  }
  const auto once = od::selectTopSamples(merged, heavy, 5, 0);
  REQUIRE(inc.size() == once.size());
  for (std::size_t k = 0; k < once.size(); ++k) CHECK(inc[k].s.energy == once[k].s.energy);

  // 9d filter: never keeps more than the 9c filter with the same threshold, and keeps the lowest energy.
  top = od::selectTopSamples(merged, heavy, 1000, rmsd);
  const od::SymmetryRmsd sym(f.x.f.mols[1], heavy.atoms(), heavy.types(), 256);
  const auto fin = od::symmetryFilter(top, heavy, sym, rmsd);
  CHECK(!fin.empty());
  CHECK(fin.size() <= top.size());
  CHECK(fin[0].s.energy == top[0].s.energy);
  for (std::size_t k = 0; k < fin.size(); ++k)
    for (std::size_t j = 0; j < k; ++j)
      if (fin[j].stam == fin[k].stam) CHECK(sym.rmsd(heavy.world(fin[k].s), heavy.world(fin[j].s)) >= rmsd);

  // midKeepSkipped: batches through selectPool give the single-9c result when no
  // sample below the final largest energy was removed.
  for (std::size_t topN : {3, 5, 8}) {
    const auto once = od::selectTopSamples(merged, heavy, static_cast<int>(topN), rmsd);
    std::vector<od::RoundSample> pool;
    float dropped = std::numeric_limits<float>::infinity();
    const std::size_t bs = merged.size() / 4 + 1;
    for (std::size_t b = 0; b < merged.size(); b += bs) {
      std::vector<od::RoundSample> m2 = std::move(pool);
      m2.insert(m2.end(), merged.begin() + b, merged.begin() + std::min(b + bs, merged.size()));
      if (b + bs >= merged.size()) {
        pool = od::selectTopSamples(std::move(m2), heavy, static_cast<int>(topN), rmsd);
        break;
      }
      auto r = od::selectPool(std::move(m2), heavy, 2 * topN, rmsd, 100000);
      for (std::size_t k = 1; k < r.pool.size(); ++k) CHECK(r.pool[k - 1].s.energy <= r.pool[k].s.energy);
      if (!r.pool.empty()) CHECK(r.pool.back().s.energy <= r.droppedMin);
      dropped = std::min(dropped, r.droppedMin);
      pool = std::move(r.pool);
    }
    const bool same = std::isinf(dropped) || (pool.size() == topN && od::topMaxEnergy(pool) < dropped);
    CHECK(same);
    if (same) {
      REQUIRE(pool.size() == once.size());
      for (std::size_t k = 0; k < once.size(); ++k) {
        CHECK(pool[k].s.energy == once[k].s.energy);
        CHECK(pool[k].stam == once[k].stam);
        CHECK(pool[k].mobm == once[k].mobm);
      }
    }
  }
  {
    // Over the cap: the kept samples stay, skipped samples fill the rest.
    const auto kept = od::selectTopSamples(merged, heavy, 4, rmsd);
    REQUIRE(kept.size() == 4);
    auto r = od::selectPool(merged, heavy, 4, rmsd, 5);
    auto full = od::selectPool(merged, heavy, 4, rmsd, 100000);
    REQUIRE(full.pool.size() > 5);
    CHECK(r.pool.size() == 5);
    for (const auto& k : kept)
      CHECK(std::any_of(r.pool.begin(), r.pool.end(), [&](const auto& p) { return p.s.energy == k.s.energy && p.mobm == k.mobm && p.stam == k.stam; }));
    for (std::size_t k = 1; k < r.pool.size(); ++k) CHECK(r.pool[k - 1].s.energy <= r.pool[k].s.energy);
    CHECK(r.droppedMin <= full.droppedMin);
    CHECK(od::selectTopSamples(r.pool, heavy, 4, rmsd).size() == 4);
  }

  // dcss filter.
  std::vector<int> dcss{1, 2, 3, 4};
  const std::unordered_map<int, float> lowest{{1, -5.f}, {2, -2.f}, {4, -3.f}};
  CHECK(od::filterDcss(dcss, lowest, -3.f) == 2);  // 2 above E_top_max, 3 without docking samples
  CHECK(dcss == std::vector<int>{1, 4});
  CHECK(od::filterDcss(dcss, lowest, -std::numeric_limits<float>::infinity()) == 2);
  CHECK(dcss.empty());
}

TEST_CASE("symrmsd: automorphisms, Hungarian pairing, symmetry-corrected RMSD") {
  // Labelled graphs.
  const std::vector<std::vector<int>> ring{{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}};
  const std::vector<int> same(6, 0);
  CHECK(od::graphAutomorphisms(ring, same, 1000).size() == 12);
  const std::vector<int> oneLabel{1, 0, 0, 0, 0, 0};
  CHECK(od::graphAutomorphisms(ring, oneLabel, 1000).size() == 2);
  CHECK(od::graphAutomorphisms(ring, same, 5).size() == 6);  // stops after limit + 1
  // Star with 4 equal leaves: 4! automorphisms.
  const std::vector<std::vector<int>> star{{1, 2, 3, 4}, {0}, {0}, {0}, {0}};
  CHECK(od::graphAutomorphisms(star, std::vector<int>(5, 0), 1000).size() == 24);
  // k independent 2-fold swaps: 2^k.
  std::vector<std::vector<int>> chain;
  std::vector<int> lab;
  const int k = 6;
  for (int i = 0; i < k; ++i) {
    const int c = static_cast<int>(chain.size());
    chain.push_back({});
    lab.push_back(10 + i);  // distinct backbone labels
    chain.push_back({c});
    chain.push_back({c});
    lab.push_back(0);
    lab.push_back(0);
    chain[c].push_back(c + 1);
    chain[c].push_back(c + 2);
    if (i > 0) {
      chain[c].push_back(c - 3);
      chain[c - 3].push_back(c);
    }
  }
  CHECK(od::graphAutomorphisms(chain, lab, 1000).size() == 64);

  // Hungarian pairing against brute force.
  std::mt19937 gen(3);
  std::uniform_real_distribution<float> u(-3, 3);
  for (int trial = 0; trial < 20; ++trial) {
    const int n = 6;
    std::vector<od::Vec3f> p(n), q(n);
    for (int i = 0; i < n; ++i) {
      p[i] = {u(gen), u(gen), u(gen)};
      q[i] = {u(gen), u(gen), u(gen)};
    }
    std::vector<int> perm{0, 1, 2, 3, 4, 5};
    double best = 1e30;
    do {
      double s = 0;
      for (int i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d) s += (p[i][d] - q[perm[i]][d]) * (p[i][d] - q[perm[i]][d]);
      best = std::min(best, s);
    } while (std::next_permutation(perm.begin(), perm.end()));
    CHECK(od::hungarianSquaredSum(p, q, {{0, 1, 2, 3, 4, 5}}) == doctest::Approx(best).epsilon(1e-6));
  }

  // Benzene (united atom: 6 aromatic carbons): rotating the labels by one atom gives RMSD 0.
  const auto mols = od::loadMolecules(kData + "/benzene.mol2", od::Role::MobM, ccd());
  REQUIRE(mols.size() == 1);
  const auto typing = od::assignAd4Types(mols[0]);
  std::vector<int> atoms, types;
  std::vector<od::Vec3f> x;
  for (std::size_t i = 0; i < mols[0].atoms.size(); ++i)
    if (typing.type[i] == od::Ad4Type::A) {
      atoms.push_back(static_cast<int>(i));
      types.push_back(static_cast<int>(typing.type[i]));
      const auto& c = mols[0].atoms[i].xyz;
      x.push_back({static_cast<float>(c[0]), static_cast<float>(c[1]), static_cast<float>(c[2])});
    }
  REQUIRE(atoms.size() == 6);
  const od::SymmetryRmsd c1(mols[0], atoms, types, 256), c2(mols[0], atoms, types, 4);
  CHECK(c1.usesPermutations());
  CHECK(c1.permutationCount() == 12);
  CHECK_FALSE(c2.usesPermutations());
  // Neighbours along the ring of atom 0, then rotate the labels along the ring.
  std::vector<int> ringOrder{0};
  while (ringOrder.size() < 6)
    for (const auto& b : mols[0].bonds) {
      const auto ia = std::find(atoms.begin(), atoms.end(), b.a), ib = std::find(atoms.begin(), atoms.end(), b.b);
      if (ia == atoms.end() || ib == atoms.end()) continue;
      const int pa = static_cast<int>(ia - atoms.begin()), pb = static_cast<int>(ib - atoms.begin());
      const int last = ringOrder.back();
      const int next = pa == last ? pb : pb == last ? pa : -1;
      if (next >= 0 && std::find(ringOrder.begin(), ringOrder.end(), next) == ringOrder.end()) {
        ringOrder.push_back(next);
        break;
      }
    }
  std::vector<od::Vec3f> rot(6);
  for (int i = 0; i < 6; ++i) rot[ringOrder[i]] = x[ringOrder[(i + 1) % 6]];
  CHECK(od::fixedRmsd(x, rot) > 1.0);
  CHECK(c1.rmsd(x, rot) == doctest::Approx(0).epsilon(1e-5));
  CHECK(c2.rmsd(x, rot) == doctest::Approx(0).epsilon(1e-5));
  // One-to-one: RMSD >= centroid distance.
  std::vector<od::Vec3f> moved = rot;
  for (auto& p : moved) p[0] += 2.5f;
  CHECK(c1.rmsd(x, moved) >= 2.5 - 1e-4);
  CHECK(c2.rmsd(x, moved) >= 2.5 - 1e-4);
}

TEST_CASE("score9d: pair recalculation and atom shaking entropy (2KJJ, 1L2Y)") {
  GridFixture x;
  const od::EnergyParams ep;
  std::vector<od::Vec3f> sx, mx;
  {
    std::vector<od::Vec3> d;
    x.staSys->evaluate(x.staSys->inputAngles(), nullptr, &d);
    for (const auto& v : d) sx.push_back({float(v[0]), float(v[1]), float(v[2])});
    x.mobSys->evaluate(x.mobSys->inputAngles(), nullptr, &d);
    for (const auto& v : d) mx.push_back({float(v[0]), float(v[1]), float(v[2])});
  }
  const od::ShakeAtoms stam = od::shakeAtoms(*x.staSys, x.staParts, sx);
  REQUIRE(stam.size() == x.stam.size());

  // dockWorld: identity rotation and t = centre keeps the coordinates.
  od::DockSample s0;
  s0.center = {1.5, -2, 3};
  s0.pose.t = s0.center;
  const auto w0 = od::dockWorld(mx, s0);
  for (std::size_t k = 0; k < mx.size(); ++k)
    for (int d = 0; d < 3; ++d) CHECK(w0[k][d] == doctest::Approx(mx[k][d]).epsilon(1e-5));

  for (float shift : {12.0f, 25.0f}) {
    std::vector<od::Vec3f> moved = mx;
    for (auto& v : moved) v[0] += shift;
    const od::ShakeAtoms mob = od::shakeAtoms(*x.mobSys, x.mobParts, moved);
    od::GridAtoms md = x.mobm;
    for (auto& v : md.xyz) v[0] += shift;
    const od::InterTerms it = od::exactInter(stam, mob, ep);
    const double ref = od::pairInterEnergy(x.stam, md, ep);
    CHECK(it.inter() == doctest::Approx(ref).epsilon(1e-4).scale(1.0));
    CHECK(it.inter() == it.vdw + it.hbond + it.elec + it.desolv);
    if (shift == 12.0f) {
      CHECK(it.vdw != 0);
      CHECK(it.elec != 0);
      CHECK(it.desolv != 0);
    }
    // Clash: nearest StaM heavy atom of every MobM heavy atom.
    double clash = 0;
    for (std::size_t i = 0; i < md.size(); ++i) {
      if (!md.heavy[i]) continue;
      double nearest = 1e30;
      for (std::size_t j = 0; j < x.stam.size(); ++j)
        if (x.stam.heavy[j]) nearest = std::min(nearest, std::hypot(md.xyz[i][0] - x.stam.xyz[j][0], md.xyz[i][1] - x.stam.xyz[j][1], md.xyz[i][2] - x.stam.xyz[j][2]));
      if (nearest < ep.clashDistance) clash += ep.clashBase + ep.clashK * (ep.clashDistance - nearest) * (ep.clashDistance - nearest);
    }
    CHECK(it.clash == doctest::Approx(clash).epsilon(1e-3));
  }

  // Sphere points.
  const auto pts = od::shakePoints(1.4, 0.2);
  CHECK(pts.size() == 1419);
  CHECK(pts[0] == od::Vec3f{0, 0, 0});
  for (const auto& p : pts) CHECK(std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) <= 1.4 + 1e-5);
  CHECK(od::shakePoints(1.0, 0.4).size() == 81);

  od::ShakeParams sp;
  sp.step = 0.4;
  // A single atom, or two atoms within three bonds, move freely: ln f = 0.
  {
    od::ShakeAtoms one;
    one.xyz = {{0, 0, 0}, {1.5f, 0, 0}};
    one.type = {od::Ad4Type::C, od::Ad4Type::C};
    one.charge = {0.1f, -0.1f};
    one.heavy = {1, 1};
    one.near = {{0, 1}, {0, 1}};
    for (double v : od::shakeLnF(one, ep, sp)) CHECK(v == doctest::Approx(0).epsilon(1e-12));
    one.near = {{0}, {1}};
    const auto l = od::shakeLnF(one, ep, sp);
    for (double v : l) {
      CHECK(v < 0);
      CHECK(v >= std::log(1.0 / 179) - 1e-9);
    }
    // Boltzmann weights exp(-dE * rrtScaling): independent of temperature; 0 makes every point equal (ln f = 0).
    CHECK(sp.rrtScaling == 4.0);
    od::ShakeParams s2 = sp;
    s2.temperature = 500;
    const auto l2 = od::shakeLnF(one, ep, s2);
    for (std::size_t i = 0; i < l.size(); ++i) CHECK(l2[i] == l[i]);
    s2.rrtScaling = 0.0;
    for (double v : od::shakeLnF(one, ep, s2)) CHECK(v == doctest::Approx(0).epsilon(1e-12));
    s2.rrtScaling = 2 * sp.rrtScaling;
    const auto l3 = od::shakeLnF(one, ep, s2);
    for (std::size_t i = 0; i < l.size(); ++i) CHECK(l3[i] < l[i]);
  }
  // Complex: StaM atoms out of reach of the MobM keep their values alone.
  std::vector<od::Vec3f> far = mx;
  for (auto& v : far) v[0] += 30;
  const od::ShakeAtoms mob = od::shakeAtoms(*x.mobSys, x.mobParts, far);
  const auto alone = od::shakeLnF(stam, ep, sp);
  for (std::size_t a = 0; a < stam.size(); ++a) {
    if (stam.heavy[a]) CHECK(alone[a] <= 1e-12);
    else CHECK(alone[a] == 0);
  }
  const auto full = od::shakeLnF(od::joinShakeAtoms(stam, mob), ep, sp);
  const double sumFull = std::accumulate(full.begin(), full.end(), 0.0);
  const od::ComplexLnF cl = od::complexSumLnF(stam, alone, mob, ep, sp);
  CHECK(cl.stam + cl.mobm == doctest::Approx(sumFull).epsilon(1e-9));
  CHECK(cl.stam == doctest::Approx(std::accumulate(full.begin(), full.begin() + static_cast<long>(stam.size()), 0.0)).epsilon(1e-9));
  CHECK(cl.mobm == doctest::Approx(std::accumulate(full.begin() + static_cast<long>(stam.size()), full.end(), 0.0)).epsilon(1e-9));
  CHECK(od::minusTS(sumFull, sp) > 0);
  CHECK(od::minusTS(-1.0, sp) == doctest::Approx(od::kGasConstant * 298.15));
}

TEST_CASE("write_struct: PDB, mmCIF, SDF and mol2") {
  const auto dir = std::filesystem::temp_directory_path() / std::format("od_test_out_{}", ::getpid());
  std::filesystem::create_directories(dir);
  for (const std::string& file : {std::string("1L2Y.cif"), std::string("2KJJ.cif"), std::string("BTN_ideal.sdf")}) {
    const auto mols = od::loadMolecules(kData + "/" + file, od::Role::MobM, ccd());
    REQUIRE(mols.size() == 1);
    const od::Molecule& m = mols[0];
    od::StructPart part{&m, {}};
    for (const auto& a : m.atoms) part.xyz.push_back({float(a.xyz[0] + 1), float(a.xyz[1]), float(a.xyz[2])});
    const std::vector<od::StructPart> parts{part};
    CHECK(od::fitsPdb(parts));
    bool fb = true;
    const auto files = od::writeStructure((dir / "x").string(), parts, {true, true, true, true}, "test", fb);
    CHECK_FALSE(fb);
    CHECK(files == std::vector<std::string>{"x.pdb", "x.cif", "x.sdf", "x.mol2"});
    for (const char* ext : {".cif", ".sdf", ".mol2"}) {
      if (std::string(ext) == ".cif" && !m.hasResidueInfo) {
        // Without residue information (UNL, no CCD entry): bonds are given in _chem_comp_bond.
        std::ifstream c(dir / "x.cif");
        const std::string text((std::istreambuf_iterator<char>(c)), std::istreambuf_iterator<char>());
        CHECK(text.find("_chem_comp_bond.comp_id") != std::string::npos);
        continue;
      }
      const auto back = od::loadMolecules((dir / (std::string("x") + ext)).string(), od::Role::MobM, ccd());
      REQUIRE(back.size() == 1);
      REQUIRE(back[0].atoms.size() == m.atoms.size());
      CHECK(back[0].bonds.size() == m.bonds.size());
      for (std::size_t a = 0; a < m.atoms.size(); ++a) {
        CHECK(back[0].atoms[a].atomicNum == m.atoms[a].atomicNum);
        if (std::string(ext) != ".mol2") CHECK(back[0].atoms[a].formalCharge == m.atoms[a].formalCharge);
        for (int d = 0; d < 3; ++d) CHECK(back[0].atoms[a].xyz[d] == doctest::Approx(part.xyz[a][d]).epsilon(1e-3).scale(1.0));
      }
      if (std::string(ext) == ".mol2") {
        // mol2 has no formal charges: the reader places them on an equivalent resonance atom.
        int q0 = 0, q1 = 0;
        for (std::size_t a = 0; a < m.atoms.size(); ++a) {
          q0 += m.atoms[a].formalCharge;
          q1 += back[0].atoms[a].formalCharge;
        }
        CHECK(q0 == q1);
      }
      if (std::string(ext) == ".mol2")
        for (std::size_t a = 0; a < m.atoms.size(); ++a)
          CHECK(back[0].atoms[a].partialCharge == doctest::Approx(m.atoms[a].partialCharge).epsilon(1e-3).scale(1.0));
    }
    // PDB: fixed columns.
    std::ifstream in(dir / "x.pdb");
    std::string line;
    std::size_t n = 0;
    while (std::getline(in, line)) {
      if (line.rfind("ATOM", 0) != 0 && line.rfind("HETATM", 0) != 0) continue;
      REQUIRE(line.size() == 80);
      CHECK(std::stod(line.substr(30, 8)) == doctest::Approx(part.xyz[n][0]).epsilon(1e-3).scale(1.0));
      CHECK(std::stod(line.substr(46, 8)) == doctest::Approx(part.xyz[n][2]).epsilon(1e-3).scale(1.0));
      ++n;
    }
    CHECK(n == m.atoms.size());
    if (file == "BTN_ideal.sdf") {
      std::map<std::string, int> types;
      for (std::size_t a = 0; a < m.atoms.size(); ++a) ++types[od::sybylType(m, static_cast<int>(a))];
      CHECK(types["N.am"] == 2);
      CHECK(types["S.3"] == 1);
      CHECK(types["O.2"] == 2);
      CHECK(types["O.3"] == 1);
    }
    // Over the PDB limits (two-character chain): mmCIF instead.
    od::Molecule big = m;
    big.hasResidueInfo = true;
    if (big.residues.empty()) {
      big.residues.push_back({{"AB", 1, "", "LIG"}, {}});
      for (auto& a : big.atoms) a.residue = 0;
    } else {
      big.residues[0].info.chain = "AB";
    }
    const std::vector<od::StructPart> bp{{&big, part.xyz}};
    CHECK_FALSE(od::fitsPdb(bp));
    const auto f2 = od::writeStructure((dir / "y").string(), bp, {true, false, false, false}, "test", fb);
    CHECK(fb);
    CHECK(f2 == std::vector<std::string>{"y.cif"});
  }
  std::filesystem::remove_all(dir);
}

TEST_CASE("feparams: defaults, file, command-line overrides and errors") {
  const od::FeParams d;
  CHECK(d.c0 == 0);
  CHECK(d.tors == 0);
  for (int k = 0; k < od::FeParams::kTerms; ++k) {
    CHECK(d.inter[k] == 0.7564);
    CHECK(d.intra[k] == 0.7564);
  }
  CHECK(d.entropyStam == 0.9);
  CHECK(d.entropyMobm == 0.9);
  CHECK(d.temperature == 298.15);
  CHECK(d.unboundCandidateNum == 20);
  CHECK(d.shakingRadius == 1.4);
  CHECK(d.shakingStep == 0.2);
  CHECK(d.rrtScaling == 4.0);
  CHECK(od::loadFeParams("", {}).inter == d.inter);
  const auto path = (std::filesystem::temp_directory_path() / std::format("od_test_fe_{}.toml", ::getpid())).string();
  auto write = [&](const std::string& s) { std::ofstream(path) << s; };
  // The written file reads back to the same values.
  od::FeParams x;
  x.c0 = -1.25;
  x.tors = 0.2983;
  x.inter = {0.1, 0.2, 0.3, 0.4};
  x.intra = {0.5, 0.6, 0.7, 0.8};
  x.entropyStam = 1.1;
  x.entropyMobm = 0.7;
  x.temperature = 310;
  x.unboundCandidateNum = 5;
  x.shakingRadius = 1.2;
  x.shakingStep = 0.25;
  x.rrtScaling = 1.6878;
  write(od::feParamsText(x));
  const od::FeParams y = od::loadFeParams(path, {});
  CHECK(y.c0 == x.c0);
  CHECK(y.tors == x.tors);
  CHECK(y.inter == x.inter);
  CHECK(y.intra == x.intra);
  CHECK(y.entropyStam == x.entropyStam);
  CHECK(y.entropyMobm == x.entropyMobm);
  CHECK(y.temperature == x.temperature);
  CHECK(y.unboundCandidateNum == x.unboundCandidateNum);
  CHECK(y.shakingRadius == x.shakingRadius);
  CHECK(y.shakingStep == x.shakingStep);
  CHECK(y.rrtScaling == x.rrtScaling);
  write(od::feParamsText(d));
  CHECK(od::loadFeParams(path, {}).rrtScaling == 4.0);
  // Partial file: other keys keep the defaults; integers are accepted for numbers.
  write("[od_0100]\ntors = 1\n[od_0100.inter]\nelec = 0.5\n[od_0100.shaking]\nrrtScaling = 2\n");
  const od::FeParams z = od::loadFeParams(path, {});
  CHECK(z.tors == 1);
  CHECK(z.inter == std::array<double, 4>{0.7564, 0.7564, 0.5, 0.7564});
  CHECK(z.intra == d.intra);
  CHECK(z.rrtScaling == 2);
  CHECK(z.shakingRadius == 1.4);
  // Command line overrides the file, in order; the last one wins.
  const od::FeParams o = od::loadFeParams(path, {{"od_0100.inter.elec", "0.25"}, {"od_0100.shaking.rrtScaling", "3"},
                                                 {"od_0100.shaking.rrtScaling", "5"}, {"od_0100.unboundCandidateNum", "7"}});
  CHECK(o.inter[2] == 0.25);
  CHECK(o.rrtScaling == 5);
  CHECK(o.unboundCandidateNum == 7);
  CHECK(o.tors == 1);
  // Errors: all listed, file problems with lines, command-line problems with their key.
  write("[od_0100]\nc0 = \"x\"\nfoo = 1\nunboundCandidateNum = 2.5\n[od_0100.inter]\nvdW = 1\n"
        "[od_0100.entropy]\nstam = nan\n[od_0100.shaking]\nstep = 0\nrrtScaling = -1\n[other]\n");
  try {
    od::loadFeParams(path, {{"od_0100.temperature", "abc"}, {"od_0100.shaking.radius", "-1"}});
    FAIL("no error");
  } catch (const od::Error& e) {
    const std::string m = e.what();
    CHECK(m.find("10 problem(s)") != std::string::npos);
    CHECK(m.find("line 2: od_0100.c0: must be a number") != std::string::npos);
    CHECK(m.find("line 3: od_0100.foo: unknown key") != std::string::npos);
    CHECK(m.find("line 4: od_0100.unboundCandidateNum: must be an integer") != std::string::npos);
    CHECK(m.find("line 6: od_0100.inter.vdW: unknown key") != std::string::npos);
    CHECK(m.find("line 8: od_0100.entropy.stam: must be finite") != std::string::npos);
    CHECK(m.find("line 10: od_0100.shaking.step: must be > 0") != std::string::npos);
    CHECK(m.find("line 11: od_0100.shaking.rrtScaling: must be >= 0") != std::string::npos);
    CHECK(m.find("other: unknown key") != std::string::npos);
    CHECK(m.find("command line: od_0100.temperature (abc): must be a number") != std::string::npos);
    CHECK(m.find("command line: od_0100.shaking.radius (-1): must be > 0") != std::string::npos);
  }
  write("[od_0100\n");
  CHECK_THROWS_AS(od::loadFeParams(path, {}), od::Error);
  std::filesystem::remove(path);
  CHECK_THROWS_AS(od::loadFeParams(path, {}), od::Error);
}
