// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

#include <unistd.h>

#include "ccd.hpp"
#include "log.hpp"
#include "dssp.hpp"
#include "molconfig.hpp"
#include "rama.hpp"
#include "read_input.hpp"
#include "topology.hpp"

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
