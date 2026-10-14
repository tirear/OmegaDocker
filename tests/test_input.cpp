// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <numbers>
#include <random>
#include <fstream>
#include <sstream>
#include <string>

#include <unistd.h>

#include "ad4.hpp"
#include "ad4_pair.hpp"
#include "ccd.hpp"
#include "log.hpp"
#include "dssp.hpp"
#include "energy.hpp"
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
  p.staFlexMargin = 0;
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
  e = f.errorsOf(replaceAfter(text, "mob:1RNG", "delta   = false", "delta   = true"));
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
  e = f.errorsOf(replaceAfter(text, "mob:BTN", "backboneType = \"none\"", "backboneType = \"none\"\ncsgm1 = \"1a\""));
  CHECK(hasError(e, "csgm1: only allowed for protein molecules"));
  e = f.errorsOf(replaceAfter(text, "mob:1L2Y", "csgm1        = \"1a\"", "csgm1        = \"1c\""));
  CHECK(hasError(e, "csgm1: must be 1a or 1b"));
  e = f.errorsOf(replaceAfter(text, "[[molecule]]\nid           = \"mob:BTN", "mob:BTN_ideal.sdf:1", "mob:BTN_ideal.sdf:2"));
  CHECK(hasError(e, "id 'mob:BTN_ideal.sdf:2' does not match any input molecule"));
  CHECK(hasError(e, "molecule: id 'mob:BTN_ideal.sdf:1' is missing"));

  // Residue and torsion lists must match the detected ones.
  e = f.errorsOf(replaceAfter(text, "resSeq  = 2", "resSeq  = 2", "resSeq  = 99"));
  CHECK(hasError(e, "residue A:99 is not a backbone residue"));
  CHECK(hasError(e, "residue A:2 is missing"));
  e = f.errorsOf(replaceAfter(text, "mob:BTN", "[\"#2\", \"#1\", \"#4\", \"#5\"]", "[\"#2\", \"#1\", \"#4\", \"#99\"]"));
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
  text = replaceAfter(text, "mob:BTN", "backboneType = \"none\"",
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
  const auto r = s.find("resSeq  = 10\n", s.find("mob:1ZNF.cif:A"));
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
  const std::string zn = "id           = \"mob:1ZNF.cif:A\"";
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
  e = f.errorsOf(replaceAfter(text, "mob:BTN", "backboneType = \"none\"",
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

  explicit StaFixture(bool staFixedResMovable = false) {
    params.staFixedResMovable = staFixedResMovable;
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
  const auto start = s.find("id           = \"sta:");
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
  const auto i = s.find("csgm1        = \"1a\"\n", s.find("id           = \"sta:"));
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

TEST_CASE("torsion tree: staFixedResMovable (2KJJ)") {
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
  const std::string brk = setQuoted(uncommentMetal(text), "id           = \"mob:1ZNF.cif:A\"", "metalBond    = ", "break");
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
  std::ofstream(f.path) << phi10(setQuoted(uncommentMetal(text), "id           = \"mob:1ZNF.cif:A\"", "metalBond    = ", "break"));
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
