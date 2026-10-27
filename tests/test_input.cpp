// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <algorithm>
#include <atomic>
#include <thread>
#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <fstream>
#include <numbers>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <tuple>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <omp.h>
#include <unistd.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "ad4.hpp"
#include "box.hpp"
#include "ad4_pair.hpp"
#include "ccd.hpp"
#include "conf9a.hpp"
#include "csgm1.hpp"
#include "csgm2.hpp"
#include "csgm3.hpp"
#include "checkpoint.hpp"
#include "dock9b.hpp"
#include "dock9c.hpp"
#include "fiber_batch.hpp"
#include "log.hpp"
#include "minimize.hpp"
#include "dssp.hpp"
#include "energy_kernel.hpp"
#ifdef OD_HAVE_CUDA
#include "conf9a_gpu.hpp"
#include "dock9b_gpu.hpp"
#include "rcd_gpu.hpp"
#include "energy_gpu.hpp"
#include "grid_gpu.hpp"
#include "shake_gpu.hpp"
#endif
#include "symrmsd.hpp"
#include "energy.hpp"
#include "grid.hpp"
#include "grow.hpp"
#include "closure.hpp"
#include "orient.hpp"
#include "propka.hpp"
#include "dimorphite.hpp"
#include "protonate.hpp"
#include "read_mmcif.hpp"
#include "score9d.hpp"
#include "shaking.hpp"
#include "sysmem.hpp"
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

TEST_CASE("csgm2: tournament, crossover and mutation") {
  const std::vector<float> g0(8, 10.0f), g1(8, -20.0f), g2(8, 50.0f);
  const std::vector<bool> noCis(8, false);
  std::vector<double> a(8);
  // Without crossover and mutation the offspring copies the first parent.
  od::Csgm2Params copy{1.0, 0.0, 0.0, 90.0};
  std::vector<od::Csgm2Parent> pp{{&g0, -5.0f}, {&g1, -1.0f}};
  for (std::uint64_t k = 0; k < 200; ++k) {
    od::Rng r(1, 2, k);
    od::csgm2Offspring(r, pp, noCis, copy, a);
    // tournamentRate 1: the fitter candidate always wins, so g1 only when both candidates are g1.
    CHECK((a[0] == 10.0 || a[0] == -20.0));
    for (double x : a) CHECK(x == a[0]);
  }
  int g0Count = 0;
  for (std::uint64_t k = 0; k < 4000; ++k) {
    od::Rng r(1, 2, k);
    od::csgm2Offspring(r, pp, noCis, copy, a);
    g0Count += a[0] == 10.0;
  }
  CHECK(g0Count == doctest::Approx(3000).epsilon(0.05));  // P(fitter in the pair) = 3/4
  // tournamentRate 0: the less fit candidate wins.
  od::Csgm2Params worse{0.0, 0.0, 0.0, 90.0};
  int g1Count = 0;
  for (std::uint64_t k = 0; k < 4000; ++k) {
    od::Rng r(1, 3, k);
    od::csgm2Offspring(r, pp, noCis, worse, a);
    g1Count += a[0] == -20.0;
  }
  CHECK(g1Count == doctest::Approx(3000).epsilon(0.05));
  // Crossover: two-point (first, second, first parent) or single-point when the points coincide (first, second).
  od::Csgm2Params cross{0.5, 1.0, 0.0, 90.0};
  std::vector<od::Csgm2Parent> p3{{&g0, 0.0f}, {&g1, 0.0f}, {&g2, 0.0f}};
  int mixed = 0;
  for (std::uint64_t k = 0; k < 500; ++k) {
    od::Rng r(4, 5, k);
    od::csgm2Offspring(r, p3, noCis, cross, a);
    int changes = 0;
    for (int i = 1; i < 8; ++i) changes += a[i] != a[i - 1];
    CHECK(changes <= 2);
    mixed += changes > 0;
    if (changes == 2) CHECK(a[7] == a[0]);  // two-point: both ends from the first parent
  }
  CHECK(mixed > 100);
  // Mutation: every gene within absMaxDang of the parent value, wrapped; cisTrans genes switch 0 <-> 180.
  const std::vector<float> h{170.0f, 0.0f, 180.0f, -100.0f};
  std::vector<bool> cis{false, true, true, false};
  std::vector<od::Csgm2Parent> ph{{&h, 0.0f}};
  od::Csgm2Params mut{0.6, 0.0, 1.0, 30.0};
  std::vector<double> b(4);
  for (std::uint64_t k = 0; k < 500; ++k) {
    od::Rng r(6, 7, k);
    od::csgm2Offspring(r, ph, cis, mut, b);
    CHECK(b[0] > -180.0);
    CHECK(b[0] <= 180.0);
    CHECK(std::abs(od::wrapDeg(b[0] - 170.0)) <= 30.0);
    CHECK(b[1] == 180.0);
    CHECK(b[2] == 0.0);
    CHECK(std::abs(b[3] + 100.0) <= 30.0);
  }
  // Same (seed, stream, index): same offspring.
  std::vector<double> c1(8), c2(8);
  od::Rng r1(9, 9, 9), r2(9, 9, 9);
  od::csgm2Offspring(r1, p3, noCis, od::Csgm2Params{}, c1);
  od::csgm2Offspring(r2, p3, noCis, od::Csgm2Params{}, c2);
  CHECK(c1 == c2);
  od::Rng r3(1, 1, 1);
  CHECK_THROWS_AS(od::csgm2Offspring(r3, std::vector<od::Csgm2Parent>{}, noCis, od::Csgm2Params{}, c1), od::Error);
}

TEST_CASE("csgm3: lines through the end points") {
  using V = std::vector<float>;
  od::Csgm3Params p0{0, 360.0, 4};
  CHECK(od::csgm3SampleNum(3, p0) == 24);
  CHECK(od::csgm3SampleNum(2000, p0) == 2000LL * 1999 * 4);
  const std::vector<bool> noCis(2, false);
  auto line = [&](std::span<const V* const> ends, long long l, const od::Csgm3Params& p, std::uint64_t seed,
                  const std::vector<bool>& cis) {
    std::vector<std::vector<double>> s;
    for (int m = 0; m < p.perLine; ++m) {
      std::vector<double> a(ends[0]->size());
      od::csgm3Sample(od::Rng(seed, 1, 100), ends, l, m, cis, p, a);
      s.push_back(a);
    }
    return s;
  };
  // No offset: the samples cover 360 degrees along the line in steps of 90, from either end point.
  const V e0{0.0f, 0.0f}, e1{90.0f, 0.0f};
  const std::vector<const V*> two{&e0, &e1};
  for (std::uint64_t seed = 0; seed < 20; ++seed)
    for (long long l = 0; l < 2; ++l) {
      auto s = line(two, l, p0, seed, noCis);
      std::vector<double> x;
      for (const auto& a : s) {
        CHECK(a[1] == doctest::Approx(0.0).epsilon(1e-12));
        x.push_back(a[0]);
      }
      std::sort(x.begin(), x.end());
      CHECK(x == std::vector<double>{-90.0, 0.0, 90.0, 180.0});
      CHECK((s[0][0] == 0.0 || s[0][0] == 90.0));
      // Line 0 goes from end 0 to end 1 (+90 per step), line 1 the other way.
      CHECK(od::wrapDeg(s[1][0] - s[0][0]) == doctest::Approx(l == 0 ? 90.0 : -90.0));
    }
  // The shortest periodic difference sets the direction: 170 -> -170 goes through 180.
  const V f0{170.0f, 0.0f}, f1{-170.0f, 0.0f};
  const std::vector<const V*> wrapEnds{&f0, &f1};
  od::Csgm3Params p8{0, 360.0, 8};
  auto s = line(wrapEnds, 0, p8, 3, noCis);
  for (int m = 1; m < 8; ++m) CHECK(od::wrapDeg(s[m][0] - s[m - 1][0]) == doctest::Approx(45.0));
  for (const auto& a : s) {
    CHECK(a[0] > -180.0);
    CHECK(a[0] <= 180.0);
  }
  // Offsets: perpendicular to the end-point direction, integer length <= peMaxOffset.
  const V g0{10.0f, 20.0f, -30.0f, 40.0f}, g1{50.0f, -60.0f, 70.0f, 0.0f}, g2{-120.0f, 160.0f, 5.0f, -90.0f};
  const std::vector<const V*> three{&g0, &g1, &g2};
  od::Csgm3Params p30{30, 360.0, 4};
  const std::vector<bool> noCis4(4, false);
  int nonZero = 0;
  for (std::uint64_t seed = 0; seed < 200; ++seed)
    for (long long l = 0; l < 6; ++l) {
      std::vector<double> a(4);
      od::csgm3Sample(od::Rng(seed, 2, 0), three, l, 0, noCis4, p30, a);
      const long long i = l / 2;
      long long j = l % 2;
      if (j >= i) ++j;
      std::vector<double> u(4);
      double un = 0;
      for (int t = 0; t < 4; ++t) {
        u[t] = od::wrapDeg(static_cast<double>((*three[j])[t]) - (*three[i])[t]);
        un += u[t] * u[t];
      }
      bool found = false;
      for (const V* e : {three[i], three[j]}) {
        double dot = 0, len = 0;
        for (int t = 0; t < 4; ++t) {
          const double v = od::wrapDeg(a[t] - (*e)[t]);
          dot += v * u[t];
          len += v * v;
        }
        len = std::sqrt(len);
        if (std::abs(dot) / std::sqrt(un) < 1e-6 && len <= 30.0 + 1e-6 && std::abs(len - std::round(len)) < 1e-6) {
          found = true;
          nonZero += len > 0.5;
        }
      }
      CHECK(found);
    }
  CHECK(nonZero > 1000);
  // Same generator state: same samples; cisTrans genes are 0 or 180.
  std::vector<double> c1(4), c2(4);
  od::csgm3Sample(od::Rng(7, 7, 7), three, 3, 2, noCis4, p30, c1);
  od::csgm3Sample(od::Rng(7, 7, 7), three, 3, 2, noCis4, p30, c2);
  CHECK(c1 == c2);
  const std::vector<bool> cis{false, true, false, true};
  for (std::uint64_t seed = 0; seed < 50; ++seed)
    for (int m = 0; m < 4; ++m) {
      std::vector<double> a(4);
      od::csgm3Sample(od::Rng(seed, 3, 0), three, seed % 6, m, cis, p30, a);
      CHECK((a[1] == 0.0 || a[1] == 180.0));
      CHECK((a[3] == 0.0 || a[3] == 180.0));
    }
  // One torsion: no perpendicular direction, so no offset.
  const V h0{-60.0f}, h1{60.0f};
  const std::vector<const V*> one{&h0, &h1};
  std::vector<double> b(1);
  od::csgm3Sample(od::Rng(1, 1, 1), one, 0, 0, std::vector<bool>{false}, p30, b);
  CHECK((b[0] == -60.0 || b[0] == 60.0));
  CHECK_THROWS_AS(od::csgm3Sample(od::Rng(1, 1, 1), std::vector<const V*>{&h0}, 0, 0, std::vector<bool>{false}, p30, b),
                  od::Error);
  CHECK_THROWS_AS(od::csgm3Sample(od::Rng(1, 1, 1), one, 2, 0, std::vector<bool>{false}, p30, b), od::Error);
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
  od::MinParams p{1000, 1e-12, 10.0, od::MinMethod::Fire};
  const auto r = od::minimizeTorsions(f, x, p);
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
  CHECK(od::minimizeTorsions(f, x, {7, 0.0, 10.0}).steps == 7);
  x = {-100, 100, -170};
  const auto r2 = od::minimizeTorsions(f, x, {1000, 1.0, 10.0});
  CHECK(r2.steps < 1000);
  x = {-100, 100, -170};
  const auto r0 = od::minimizeTorsions(f, x, {0, 1.0, 10.0});
  CHECK(r0.steps == 0);
  CHECK(x[0] == -100);
  // Every method: converges, respects the angle step limit, returns the lowest point.
  for (const char* name : {"sd", "cg", "fire", "lbfgs", "adadelta"}) {
    const std::string nameS = name;
    CAPTURE(nameS);
    visited.clear();
    // ADADELTA (rho 0.8, epsilon 1e-2, variables in degrees) takes small steps near the minimum
    // and ends close to it only.
    const bool ada = nameS == "adadelta";
    x = {-100, 100, -170};
    const auto rm = od::minimizeTorsions(f, x, {ada ? 50000 : 5000, 1e-14, 10.0, od::parseMinMethod(name)});
    CHECK(rm.energy < (ada ? 1e-2 : 1e-4));
    for (std::size_t i = 0; i < x.size(); ++i) CHECK(std::abs(od::wrapDeg(x[i] - c[i])) < (ada ? 3.0 : 0.5));
    std::vector<double> gg;
    CHECK(f(x, gg) == doctest::Approx(rm.energy).epsilon(1e-9).scale(1e-9));
    double bigM = 0;
    for (std::size_t s = 1; s < visited.size(); ++s)
      for (std::size_t i = 0; i < x.size(); ++i) bigM = std::max(bigM, std::abs(od::wrapDeg(visited[s][i] - visited[s - 1][i])));
    CHECK(bigM <= 10.0 + 1e-9);
    CHECK(rm.evaluations >= rm.steps + 1);
  }
}

namespace {

std::array<double, 7> termArray(const od::EnergyTerms& e) {
  return {e.vdw, e.hbond, e.elec, e.desolv, e.clash, e.closure, e.closureMaxDeviation};
}
std::array<double, 7> termArray(const od::FlatTerms& e) {
  return {e.vdw, e.hbond, e.elec, e.desolv, e.clash, e.closure, e.closureMaxDeviation};
}

// ConfEvaluator on the CPU (EnergySystem::evaluate), recording the batch sizes.
// With nLanes > 1 the asynchronous form (submit evaluates, complete returns the stored results).
struct CpuConfEvaluator : od::ConfEvaluator {
  const od::EnergySystem& sys;
  int batch, nLanes;
  int calls = 0, largest = 0, idles = 0;
  struct Done {
    bool pending = false;
    std::vector<od::EnergyTerms> terms;
    std::vector<double> grad;
  };
  std::vector<Done> done;
  CpuConfEvaluator(const od::EnergySystem& s, int b, int l = 1) : sys(s), batch(b), nLanes(l), done(l) {}
  int maxBatch() const override { return batch; }
  int lanes() const override { return nLanes; }
  void submit(int lane, const double* a, int n) override {
    Done& d = done[lane];
    REQUIRE_FALSE(d.pending);
    d.pending = true;
    d.terms.resize(n);
    d.grad.resize(static_cast<std::size_t>(n) * sys.torsionCount());
    evaluate(a, n, d.terms.data(), d.grad.data());
  }
  void complete(int lane, od::EnergyTerms* terms, double* grad) override {
    Done& d = done[lane];
    REQUIRE(d.pending);
    d.pending = false;
    std::copy(d.terms.begin(), d.terms.end(), terms);
    std::copy(d.grad.begin(), d.grad.end(), grad);
  }
  void evaluate(const double* a, int n, od::EnergyTerms* terms, double* grad) override {
    ++calls;
    largest = std::max(largest, n);
    const int d = sys.torsionCount();
    std::vector<double> g;
    for (int i = 0; i < n; ++i) {
      terms[i] = sys.evaluate(std::span<const double>(a + static_cast<std::size_t>(i) * d, d), grad ? &g : nullptr);
      if (grad) std::copy(g.begin(), g.end(), grad + static_cast<std::size_t>(i) * d);
    }
  }
  void idle() override { ++idles; }
};

}  // namespace

TEST_CASE("9a: batched minimization equals per-sample minimization (1L2Y)") {
  Csgm1Fixture x;
  const od::Molecule& m = x.m();
  const od::EnergySystem sys({{&m, &x.s().tree}}, {});
  const od::Csgm1Sampler smp(m, x.s(), *x.cfg, {});
  const int d = sys.torsionCount(), count = 20;
  std::vector<double> start(static_cast<std::size_t>(count) * d);
  for (int i = 0; i < count; ++i) {
    od::Rng r(3, 0, i);
    smp.sample(r, std::span<double>(start.data() + static_cast<std::size_t>(i) * d, d));
  }
  for (const char* method : {"sd", "cg", "fire", "lbfgs", "adadelta"}) {
    od::MinParams mp;
    mp.method = od::parseMinMethod(method);
    std::vector<double> ref = start;
    std::vector<od::EnergyTerms> refTerms(count);
    std::vector<int> refSteps(count);
    for (int i = 0; i < count; ++i)
      refTerms[i] = od::minimizeConf(sys, std::span<double>(ref.data() + static_cast<std::size_t>(i) * d, d), mp, &refSteps[i]);
    for (int batch : {1, 7, 64})
      for (int lanes : {1, 2, 3}) {
        CpuConfEvaluator ev(sys, batch, lanes);
        std::vector<double> a = start;
        std::vector<od::EnergyTerms> terms(count);
        std::vector<int> steps(count);
        od::minimizeConfBatch(ev, d, a, mp, terms.data(), steps.data());
        INFO(method, " batch ", batch, " lanes ", lanes);
        CHECK(a == ref);
        CHECK(steps == refSteps);
        int diff = 0;
        for (int i = 0; i < count; ++i) diff += termArray(terms[i]) != termArray(refTerms[i]) || terms[i].clashPairs != refTerms[i].clashPairs;
        CHECK(diff == 0);
        // The final evaluation (no gradients) uses evaluate with up to batch samples.
        CHECK(ev.largest == std::min(batch, count));
      }
  }
}

TEST_CASE("9a: L2 pool refill with a ConfEvaluator equals the CPU path (1L2Y)") {
  Csgm1Fixture x;
  const od::Molecule& m = x.m();
  const od::EnergySystem sys({{&m, &x.s().tree}}, {});
  const od::Csgm1Sampler smp(m, x.s(), *x.cfg, {});
  int calls = 0;
  const od::ConfGenerator gen = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> a) {
    smp.sample(r, a);
#pragma omp atomic
    ++calls;
    return idx.global % 5 == 3 ? od::GenResult::GrowFailed : od::GenResult::Built;
  };
  od::Fill9aParams p;
  p.minL2 = 24;
  p.maxL1Batch = 10;
  p.seed = 5;
  p.stream = 1;
  std::vector<od::ConfSample> a, b;
  std::uint64_t ca = 0, cb = 0;
  const auto sa = od::fillL2Pool(sys, gen, 6, p, ca, a);
  CpuConfEvaluator ev(sys, 4);
  p.evaluator = &ev;
  const auto sb = od::fillL2Pool(sys, gen, 6, p, cb, b);
  CHECK(ca == cb);
  CHECK(sa.l1 == sb.l1);
  CHECK(sa.growFailed == sb.growFailed);
  CHECK(sa.growFailed > 0);
  CHECK(sa.minSteps == sb.minSteps);
  CHECK(ev.idles == 1);
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    CHECK(a[i].angles == b[i].angles);
    CHECK(a[i].energy == b[i].energy);
  }
  // GenBatch: the generators run as coroutines that wait for a flush; the same samples. Batches
  // smaller than genBatchMin run on OpenMP threads.
  struct Batch : od::GenBatch {
    int flushes = 0, waits = 0;
    int slots() const override { return 3; }
    void flush(std::span<const int> s) override {
      ++flushes;
      waits += static_cast<int>(s.size());
    }
  } gb;
  const od::ConfGenerator waiting = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> a) {
    if (idx.slot >= 0) {
      CHECK(idx.slot < gb.slots());
      od::batchWait();
    }
    return gen(r, idx, a);
  };
  for (long long minBatch : {1LL, 7LL}) {
    std::vector<od::ConfSample> c;
    std::uint64_t cc = 0;
    p.evaluator = nullptr;
    p.genBatch = &gb;
    p.genBatchMin = minBatch;
    gb.flushes = gb.waits = 0;
    const auto sc = od::fillL2Pool(sys, waiting, 6, p, cc, c);
    CHECK(cc == ca);
    CHECK(sc.l1 == sa.l1);
    CHECK(sc.growFailed == sa.growFailed);
    REQUIRE(c.size() == a.size());
    for (std::size_t i = 0; i < a.size(); ++i) CHECK(c[i].angles == a[i].angles);
    // L1 batches: 6, then batches of at most 10 (maxL1Batch).
    if (minBatch == 1) CHECK(gb.waits == sa.l1);
    else CHECK(gb.waits == sa.l1 - 6);
    CHECK(gb.flushes > 0);
  }
}

TEST_CASE("9a: L2 pool refill (1L2Y)") {
  Csgm1Fixture x;
  const od::Molecule& m = x.m();
  const od::EnergySystem sys({{&m, &x.s().tree}}, {});
  const od::Csgm1Sampler smp(m, x.s(), *x.cfg, {});
  const od::ConfGenerator gen = [&](od::Rng& r, const od::GenIndex&, std::span<double> a) {
    smp.sample(r, a);
    return od::GenResult::Built;
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
  CHECK(st.l1 == st.added + st.overEmax + st.notClosed + st.preNotClosed);
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

  // fixedL1: every batch has that size; the last one is cut at maxL1Total. GenIndex.local restarts per batch.
  p.fixedL1 = 7;
  p.maxL1Total = 18;
  std::vector<std::pair<std::uint64_t, long long>> seen;
  const od::ConfGenerator rec = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> a) {
#pragma omp critical
    seen.emplace_back(idx.global, idx.local);
    smp.sample(r, a);
    return od::GenResult::Built;
  };
  counter = 5;
  CHECK_THROWS_AS(od::fillL2Pool(sys, rec, 100, p, counter, c), od::Error);
  CHECK(counter == 5 + 18);
  std::sort(seen.begin(), seen.end());
  REQUIRE(seen.size() == 18);
  for (std::size_t i = 0; i < seen.size(); ++i) {
    CHECK(seen[i].first == 5 + i);
    CHECK(seen[i].second == static_cast<long long>(i % 7));
  }
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

  // Without backtracking: built samples have no clash; retries resample whole segments.
  od::GrowParams noBack;
  noBack.maxBacktrack = 0;
  const od::GrowBuilder flat(sys, parts, {&x.s().topology}, smp.segments(), noBack);
  int okFlat = 0, retried = 0;
  for (int k = 0; k < 300; ++k) {
    od::Rng rng(22, 0, k);
    smp.sample(rng, a);
    const auto r = flat.build(rng, a, resample);
    CHECK(r.reachFails == 0);
    CHECK(r.backtracks == 0);
    retried += r.retries > 0;
    if (!r.ok) {
      CHECK(r.clashFails == r.retries + 1);
      CHECK(r.retries >= 16);  // the last segment failed its first try and 16 retries
      continue;
    }
    ++okFlat;
    CHECK(r.clashFails == r.retries);
    CHECK(sys.evaluate(a).clashPairs == 0);
  }
  CHECK(okFlat > 0);
  CHECK(retried > 0);
  // Default (growMaxBacktrack 4): every failed check is followed by a retry, a backtrack or the discard.
  const od::GrowBuilder grow(sys, parts, {&x.s().topology}, smp.segments(), {});
  int ok = 0, backtracked = 0;
  for (int k = 0; k < 300; ++k) {
    od::Rng rng(22, 0, k);
    smp.sample(rng, a);
    const auto r = grow.build(rng, a, resample);
    backtracked += r.backtracks > 0;
    CHECK(r.clashFails == r.retries + r.backtracks + (r.ok ? 0 : 1));
    if (!r.ok) continue;
    ++ok;
    CHECK(sys.evaluate(a).clashPairs == 0);
  }
  MESSAGE(std::format("grow 1L2Y: built {} of 300 without backtracking, {} with ({} backtracked)", okFlat, ok, backtracked));
  CHECK(backtracked > 0);
  CHECK(ok >= okFlat);
  // No retries: a sample is either built at once or discarded.
  od::GrowParams once;
  once.maxRetry = 0;
  once.maxBacktrack = 0;
  const od::GrowBuilder noRetry(sys, parts, {&x.s().topology}, smp.segments(), once);
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
  // With loop closure the reach checks are off (spec 7.8).
  od::GrowParams noReach{0, 0.0};
  noReach.reach = false;
  CHECK(od::GrowBuilder(sys, parts, {&s.topology}, smp.segments(), noReach).reachCheckCount() == 0);
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

double dist(const od::Vec3& a, const od::Vec3& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

// Largest distance of the closure atom pairs of one link (single-part system).
double linkDeviation(const od::TorsionTree& tree, const od::ClosureLink& l, const std::vector<od::Vec3>& x) {
  if (l.restraint >= 0) {
    const auto& r = tree.restraints[l.restraint];
    return std::abs(dist(x[r.a], x[r.b]) - r.r0);
  }
  double d = 0;
  for (int v : {l.copyB, l.copyY, l.copyA})
    if (v >= 0) d = std::max(d, dist(x[tree.atomCount + v], x[tree.virtuals[v].copyOf]));
  return d;
}

}  // namespace

TEST_CASE("9a: closure deviation checked before the minimization (2KJJ)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 10, 14));
  const od::MolSettings& s = cfg.molecules[0];
  const od::EnergySystem sys({{&f.mols[0], &s.tree}}, {});
  REQUIRE(sys.hasClosures());
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  std::vector<double> a(sys.torsionCount());
  for (int k = 0; k < 50; ++k) {
    od::Rng r(4, 0, k);
    smp.sample(r, a);
    CHECK(sys.closureDeviation(a) == doctest::Approx(sys.evaluate(a).closureMaxDeviation).epsilon(1e-12));
  }
  CHECK(sys.closureDeviation(sys.inputAngles()) < 1e-6);
  // Random loops (open) and, for every third sample, the input conformation with small changes (closed).
  const std::vector<double> in = sys.inputAngles();
  const od::ConfGenerator gen = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> x) {
    smp.sample(r, x);
    if (idx.global % 3 == 0)
      for (std::size_t q = 0; q < x.size(); ++q) x[q] = in[q] + 2.0 * (r.uniform() - 0.5);
    return od::GenResult::Built;
  };
  od::Fill9aParams p;
  p.minL2 = 3;
  p.maxL1Batch = 40;
  p.seed = 2;
  p.min.maxSteps = 20;
  const auto fill = [&](double factor, std::vector<od::ConfSample>& pool) {
    p.preCloseFactor = factor;
    std::uint64_t counter = 0;
    return od::fillL2Pool(sys, gen, 40, p, counter, pool);
  };
  p.maxL1Total = 400;
  std::vector<od::ConfSample> on, off;
  const auto a1 = fill(1.0, on);
  const auto a0 = fill(0.0, off);
  CHECK(a1.preNotClosed > 0);
  CHECK(a0.preNotClosed == 0);
  CHECK(a1.l1 == a1.added + a1.overEmax + a1.notClosed + a1.preNotClosed);
  // Only the samples discarded before the minimization differ: the others give the same results.
  CHECK(a0.notClosed >= a1.notClosed);
  CHECK(on.size() == off.size());
  for (std::size_t i = 0; i < on.size() && i < off.size(); ++i) CHECK(on[i].angles == off[i].angles);
}

TEST_CASE("closure: RCD and KIC close a loop (2KJJ A:12-17)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 12, 17));
  const od::MolSettings& s = cfg.molecules[0];
  const od::TorsionTree& tree = s.tree;
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  // Link 0: the loop break (12 torsions); then two disulfides with 4 torsions each (best effort), not coupled to it.
  REQUIRE(tree.links.size() == 3);
  const auto& loop = tree.links[0];
  CHECK(f.mols[0].atoms[loop.a].name == "C");
  CHECK(f.mols[0].atoms[loop.b].name == "N");
  // The random side chains of csgm1a often clash with the closed loop; without the clash check every sample closes.
  for (const auto& [method, clash] : {std::pair{od::ClosureMethod::Rcd, 1.8}, std::pair{od::ClosureMethod::Kic, 1.8},
                                      std::pair{od::ClosureMethod::Rcd, 0.0}, std::pair{od::ClosureMethod::Kic, 0.0}}) {
    CAPTURE(static_cast<int>(method));
    CAPTURE(clash);
    od::ClosureParams cp;
    cp.method = method;
    cp.clashDistance = clash;
    const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
    REQUIRE(closer.active());
    const auto lines = closer.describe();
    for (const auto& l : lines) MESSAGE(l);
    REQUIRE(lines.size() == 3);
    CHECK(lines[0].find(method == od::ClosureMethod::Kic ? "12 torsions, 6 KIC pivots, kic" : "12 torsions, 6 KIC pivots, rcd") !=
          std::string::npos);
    CHECK(lines[1].find("best effort") != std::string::npos);
    CHECK(lines[2].find("best effort") != std::string::npos);
    int closed = 0;
    const int total = 40;
    for (int k = 0; k < total; ++k) {
      od::Rng r(7, 0, k);
      std::vector<double> a = sys.inputAngles();
      smp.sample(r, a);
      const auto x0 = sys.coordinates(std::span<const double>(a));
      CHECK(linkDeviation(tree, loop, x0) > 1.0);
      const std::vector<double> before = a;
      if (!closer.close(r, a, resample)) {
        CHECK(a == before);
        continue;
      }
      ++closed;
      const auto x = sys.coordinates(std::span<const double>(a));
      CHECK(linkDeviation(tree, loop, x) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
      if (method == od::ClosureMethod::Kic) CHECK(linkDeviation(tree, loop, x) < 1e-4);
      // Best effort does not increase the deviation of the disulfides.
      for (int li : {1, 2}) CHECK(linkDeviation(tree, tree.links[li], x) <= linkDeviation(tree, tree.links[li], x0) + 1e-6);
      // Torsions off the closure paths are unchanged; cisTrans omegas are never adjusted.
      int changed = 0;
      for (int g = 0; g < sys.torsionCount(); ++g)
        if (a[g] != before[g]) {
          ++changed;
          CHECK_FALSE(smp.cisTrans(g));
        }
      CHECK(changed <= 20);
    }
    MESSAGE(std::format("closed {} of {}", closed, total));
    CHECK(closed >= (clash > 0 ? total / 4 : total * 9 / 10));
  }
  // rcdStallTolerance: a best-effort disulfide stops at a stall; the deviation reached equals
  // that of rcdMaxIter updates, in a fraction of the time. Compared for the first disulfide when
  // the loop ended alike (fewer updates draw fewer random numbers, so later links differ).
  {
    od::ClosureParams cs, c0;
    cs.clashDistance = c0.clashDistance = 0.0;
    c0.rcdStallTolerance = 0.0;
    CHECK(cs.rcdStallTolerance == 1e-4);
    const od::LoopCloser stall(sys, parts, freeMask, smp.segments(), cs), full(sys, parts, freeMask, smp.segments(), c0);
    double tStall = 0, tFull = 0, worst = 0;
    int compared = 0;
    for (int k = 0; k < 20; ++k) {
      od::Rng r(11, 0, k);
      std::vector<double> a = sys.inputAngles();
      smp.sample(r, a);
      std::vector<double> as = a, af = a;
      od::Rng rs = r, rf = r;
      auto t0 = std::chrono::steady_clock::now();
      const bool okS = stall.close(rs, as, resample);
      auto t1 = std::chrono::steady_clock::now();
      const bool okF = full.close(rf, af, resample);
      auto t2 = std::chrono::steady_clock::now();
      tStall += std::chrono::duration<double>(t1 - t0).count();
      tFull += std::chrono::duration<double>(t2 - t1).count();
      if (!okS || !okF) continue;
      const auto xs = sys.coordinates(std::span<const double>(as)), xf = sys.coordinates(std::span<const double>(af));
      if (std::abs(linkDeviation(tree, loop, xs) - linkDeviation(tree, loop, xf)) > 1e-9) continue;
      ++compared;
      worst = std::max(worst, std::abs(linkDeviation(tree, tree.links[1], xs) - linkDeviation(tree, tree.links[1], xf)));
    }
    MESSAGE(std::format("rcdStallTolerance 1e-4: {:.3f} s, 0: {:.3f} s; {} compared, largest disulfide deviation difference {:.2e} Å",
                        tStall, tFull, compared, worst));
    CHECK(compared >= 5);
    CHECK(worst < 1e-3);
    CHECK(tStall < tFull);
  }
  // Same random stream: same result.
  od::ClosureParams cp;
  cp.method = od::ClosureMethod::Kic;
  const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
  std::vector<double> a1 = sys.inputAngles(), a2;
  od::Rng r1(9, 0, 1);
  smp.sample(r1, a1);
  a2 = a1;
  od::Rng r2 = r1;
  CHECK(closer.close(r1, a1, resample) == closer.close(r2, a2, resample));
  CHECK(a1 == a2);
  // The input conformation is already closed: nothing changes.
  std::vector<double> in = sys.inputAngles();
  od::Rng r3(1, 0, 0);
  CHECK(closer.close(r3, in, resample));
  CHECK(in == sys.inputAngles());
  // none: no closure.
  cp.method = od::ClosureMethod::None;
  CHECK_FALSE(od::LoopCloser(sys, parts, freeMask, smp.segments(), cp).active());
}

TEST_CASE("closure: coupled links closed together (2KJJ A:8-15 with the disulfide A:6-A:11)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 8, 15));
  const od::MolSettings& s = cfg.molecules[0];
  const od::TorsionTree& tree = s.tree;
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  // The loop break (16 torsions) and the disulfide A:6-A:11 (Cys 11 in the loop) share torsions: one group of
  // 20 torsions, closed by RCD (KIC is not applicable to a group); the disulfide A:20-B:19 alone is best effort.
  REQUIRE(tree.links.size() == 3);
  const auto& loop = tree.links[0];
  CHECK(f.mols[0].atoms[loop.a].name == "C");
  CHECK(f.mols[0].atoms[loop.b].name == "N");
  // The random side chains of csgm1a often clash with the closed loop; without the clash check every sample closes.
  // growClose: the clash check leaves out the side chains built after the closure.
  struct Run {
    od::ClosureMethod method;
    double clash;
    bool growClose;
  };
  for (const auto& [method, clash, growClose] : {Run{od::ClosureMethod::Rcd, 1.8, false}, Run{od::ClosureMethod::Kic, 1.8, false},
                                                 Run{od::ClosureMethod::Rcd, 1.8, true}, Run{od::ClosureMethod::Rcd, 0.0, false},
                                                 Run{od::ClosureMethod::Kic, 0.0, false}}) {
    CAPTURE(static_cast<int>(method));
    CAPTURE(clash);
    CAPTURE(growClose);
    od::ClosureParams cp;
    cp.method = method;
    cp.clashDistance = clash;
    cp.growClose = growClose;
    const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
    REQUIRE(closer.active());
    const auto lines = closer.describe();
    for (const auto& l : lines) MESSAGE(l);
    REQUIRE(lines.size() == 2);
    CHECK(lines[0].find(method == od::ClosureMethod::Kic
                            ? "A:15:C - A:16:N + A:6:SG - A:11:SG: coupled, 2 links, 20 torsions, rcd (kic not applicable"
                            : "A:15:C - A:16:N + A:6:SG - A:11:SG: coupled, 2 links, 20 torsions, rcd") != std::string::npos);
    CHECK(lines[1].find("A:20:SG - B:19:SG") != std::string::npos);
    CHECK(lines[1].find("best effort") != std::string::npos);
    int closed = 0;
    const int total = 40;
    for (int k = 0; k < total; ++k) {
      od::Rng r(7, 0, k);
      std::vector<double> a = sys.inputAngles();
      smp.sample(r, a);
      const auto x0 = sys.coordinates(std::span<const double>(a));
      CHECK(linkDeviation(tree, loop, x0) > 1.0);
      const std::vector<double> before = a;
      if (!closer.close(r, a, resample)) {
        CHECK(a == before);
        continue;
      }
      ++closed;
      const auto x = sys.coordinates(std::span<const double>(a));
      // Both links of the group are closed; best effort does not increase the deviation of A:20-B:19.
      for (int li : {0, 1}) CHECK(linkDeviation(tree, tree.links[li], x) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
      CHECK(linkDeviation(tree, tree.links[2], x) <= linkDeviation(tree, tree.links[2], x0) + 1e-6);
      // Torsions off the closure paths are unchanged; cisTrans omegas are never adjusted.
      int changed = 0;
      for (int g = 0; g < sys.torsionCount(); ++g)
        if (a[g] != before[g]) {
          ++changed;
          CHECK_FALSE(smp.cisTrans(g));
        }
      CHECK(changed <= 24);
    }
    MESSAGE(std::format("closed {} of {}", closed, total));
    // With the clash check most closed groups clash with residues 1-6 (1.4-1.7 Å), packing rather than geometry.
    CHECK(closed >= (clash == 0 ? total * 3 / 4 : 0));
  }
}

TEST_CASE("box: automatic box (spec 3.3)") {
  // An elongated cloud along a turned axis.
  std::vector<od::Vec3> x;
  od::Rng r(5, 0, 0);
  const od::Vec3 u{0.6, 0.64, 0.48}, v{-0.8, 0.48, 0.36}, w{0.0, -0.6, 0.8};
  for (int k = 0; k < 400; ++k) {
    const double a = r.uniform(-15, 15), b = r.uniform(-4, 4), c = r.uniform(-2, 2);
    x.push_back({10 + a * u[0] + b * v[0] + c * w[0], -3 + a * u[1] + b * v[1] + c * w[1], 7 + a * u[2] + b * v[2] + c * w[2]});
  }
  od::AutoBoxParams p;
  p.stamMargin = 4;
  p.mobmRadius = 3;
  const double add = 4 + 2 * 3;
  auto volume = [](const od::BoxFrame& b) { return b.size[0] * b.size[1] * b.size[2]; };
  auto contains = [&](const od::BoxFrame& b, double margin) {
    for (const auto& q : x) {
      const od::Vec3 y = b.toBox(q);
      for (int k = 0; k < 3; ++k)
        if (std::abs(y[k] - b.center[k]) > b.size[k] / 2 - margin + 1e-9) return false;
    }
    return true;
  };
  p.autoRotate = false;
  const od::BoxFrame plain = od::autoBox(x, p);
  CHECK_FALSE(plain.rotated);
  CHECK(contains(plain, add / 2));
  p.autoRotate = true;
  const od::BoxFrame rot = od::autoBox(x, p);
  CHECK(rot.rotated);
  CHECK(contains(rot, add / 2));
  CHECK(volume(rot) < 0.5 * volume(plain));
  // About the generating axes: extents 30 x 8 x 4 plus the margins.
  CHECK(volume(rot) <= (30 + add) * (8 + add) * (4 + add) * 1.01);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double d = rot.axes[i][0] * rot.axes[j][0] + rot.axes[i][1] * rot.axes[j][1] + rot.axes[i][2] * rot.axes[j][2];
      CHECK(d == doctest::Approx(i == j ? 1.0 : 0.0).epsilon(1e-12));
    }
  const od::Vec3 back = rot.toInput(rot.toBox(x[7]));
  for (int k = 0; k < 3; ++k) CHECK(back[k] == doctest::Approx(x[7][k]).epsilon(1e-12));
  // A given center: the box is symmetric about it and still holds the points with the margins.
  p.center = od::Vec3{12, -2, 8};
  const od::BoxFrame c = od::autoBox(x, p);
  const od::Vec3 cIn = c.toInput(c.center);
  for (int k = 0; k < 3; ++k) CHECK(cIn[k] == doctest::Approx((*p.center)[k]).epsilon(1e-12));
  CHECK(contains(c, add / 2));
  // A given size replaces the computed one.
  p.size = od::Vec3{20, 21, 22};
  const od::BoxFrame s = od::autoBox(x, p);
  CHECK(s.size == *p.size);
}

TEST_CASE("box: MobM radius (box.mobmRadius) and fixAuto") {
  auto mols = od::loadMolecules(kData + "/BTN_ideal.sdf", od::Role::MobM, ccd());
  const od::Molecule& m = mols[0];
  od::Topology topo = od::analyzeTopology(m, {});
  const std::vector<od::MolInput> in{{&m, &topo}};
  std::vector<std::string> errors;
  auto cfg = od::parseMolConfig(od::defaultMolConfigText(in), in, {}, errors);
  REQUIRE(cfg);
  const od::MolSettings& s = cfg->molecules[0];
  REQUIRE(s.tree.torsions.size() >= 4);
  const double rInput = od::radiusOf(od::heavyCoords(m));
  const double rMax = od::mobmRadiusMax(m, s.tree);
  const od::EnergySystem sys({{&m, &s.tree}}, {});
  od::Csgm1GroupSampler smp;
  smp.add(od::Csgm1Sampler(m, s, *cfg, {}));
  const od::SampledRadius rs = od::mobmRadiusSampled(m, sys, smp, 500, 3, 9);
  MESSAGE(std::format("BTN radius: input {:.3f}, sampled {:.3f} ({} of 500 without clashes), max {:.3f} Å", rInput,
                      rs.radius, rs.accepted, rMax));
  CHECK(rs.accepted > 0);
  CHECK(rs.radius >= rInput);
  CHECK(rMax >= rs.radius);
  // Also for every sample with clashes: the bound holds for any conformation.
  for (int k = 0; k < 200; ++k) {
    od::Rng r(11, 0, k);
    std::vector<double> a = sys.inputAngles();
    smp.sample(r, a);
    const auto c = sys.coordinates(std::span<const double>(a));
    std::vector<od::Vec3> h;
    for (std::size_t i = 0; i < m.atoms.size(); ++i)
      if (m.atoms[i].atomicNum > 1) h.push_back(c[sys.coordOffset(0) + i]);
    CHECK(od::radiusOf(h) <= rMax + 1e-9);
  }
  // Without torsions: max_i (1/N) sum_j |x_i - x_j|.
  const od::TorsionTree rigid{.atomCount = static_cast<int>(m.atoms.size())};
  const auto h = od::heavyCoords(m);
  double ref = 0;
  for (const auto& p : h) {
    double sum = 0;
    for (const auto& q : h) sum += std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
    ref = std::max(ref, sum / static_cast<double>(h.size()));
  }
  CHECK(od::mobmRadiusMax(m, rigid) == doctest::Approx(ref).epsilon(1e-12));
  CHECK(od::mobmRadiusMax(m, rigid) >= rInput);

  // fixAuto: no StaM torsion is rotatable by default; MobM unchanged.
  const od::Molecule sta = od::loadMolecules(kData + "/2KJJ.cif", od::Role::StaM, ccd()).front();
  od::TopologyParams tp;
  tp.stamAllFixed = true;
  const od::Topology ts = od::analyzeTopology(sta, tp);
  CHECK(std::none_of(ts.torsions.begin(), ts.torsions.end(), [](const od::Torsion& x) { return x.rotatable; }));
  const od::Topology tsDefault = od::analyzeTopology(sta, {});
  CHECK(std::any_of(tsDefault.torsions.begin(), tsDefault.torsions.end(), [](const od::Torsion& x) { return x.rotatable; }));
  const od::Topology tm = od::analyzeTopology(m, tp);
  CHECK(std::count_if(tm.torsions.begin(), tm.torsions.end(), [](const od::Torsion& x) { return x.rotatable; }) ==
        std::count_if(topo.torsions.begin(), topo.torsions.end(), [](const od::Torsion& x) { return x.rotatable; }));
}

#ifdef OD_HAVE_CUDA
TEST_CASE("closure: RCD on the GPU (2KJJ A:8-15 with the disulfide A:6-A:11)") {
  if (od::gpuDeviceCount() == 0) return;
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 8, 15));
  const od::MolSettings& s = cfg.molecules[0];
  const od::TorsionTree& tree = s.tree;
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  od::ClosureParams cp;
  cp.clashDistance = 0.0;
  const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
  const std::vector<int> dev0{0};
  od::GpuRcd gpu(closer, dev0, 64);
  REQUIRE(gpu.active());
  // Both links can be closed on the GPU (no clash term, no metal distance).
  const od::RcdFlat fl = closer.flat();
  REQUIRE(fl.eligible.size() == 2);
  CHECK(fl.eligible[0] == 1);
  CHECK(fl.eligible[1] == 1);
  const int total = 200;
  std::vector<std::vector<double>> ac(total), ag(total);
  std::vector<char> okC(total), okG(total);
  for (int k = 0; k < total; ++k) {
    od::Rng r(7, 0, k);
    ac[k] = sys.inputAngles();
    smp.sample(r, ac[k]);
    ag[k] = ac[k];
    okC[k] = closer.close(r, ac[k], resample);
  }
  od::runBatched(total, gpu.slots(), [&](int k, int slot) {
    od::Rng r(7, 0, k);
    std::vector<double> a = sys.inputAngles();
    smp.sample(r, a);
    okG[k] = closer.close(r, a, resample, {&gpu, slot});
    ag[k] = a;
  }, [&](std::span<const int> sl) { gpu.flush(sl); });
  int closedC = 0, closedG = 0;
  double devC = 0, devG = 0;
  for (int k = 0; k < total; ++k) {
    const auto x0 = sys.coordinates(std::span<const double>(ag[k]));
    if (okG[k]) {
      ++closedG;
      for (int li : {0, 1}) CHECK(linkDeviation(tree, tree.links[li], x0) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
      devG += linkDeviation(tree, tree.links[2], x0);
    }
    if (okC[k]) {
      ++closedC;
      devC += linkDeviation(tree, tree.links[2], sys.coordinates(std::span<const double>(ac[k])));
    }
  }
  MESSAGE(std::format("closed CPU {} GPU {} of {}; mean best-effort deviation of A:20-B:19 CPU {:.4f} GPU {:.4f} Å; "
                      "GPU {} attempts, {} updates in {} launches",
                      closedC, closedG, total, devC / std::max(1, closedC), devG / std::max(1, closedG), gpu.attempts(),
                      gpu.updates(), gpu.kernels()));
  CHECK(std::abs(closedC - closedG) <= total / 10);
  CHECK(devG / std::max(1, closedG) <= devC / std::max(1, closedC) * 1.1 + 0.01);
  // The slots divided among two devices (the same device twice with one GPU): the same results.
  const std::vector<int> two{0, od::gpuDeviceCount() > 1 ? 1 : 0};
  od::GpuRcd gpu2(closer, two, 64);
  REQUIRE(gpu2.active());
  std::vector<std::vector<double>> a2(total);
  std::vector<char> ok2(total);
  od::runBatched(total, gpu2.slots(), [&](int k, int slot) {
    od::Rng r(7, 0, k);
    std::vector<double> a = sys.inputAngles();
    smp.sample(r, a);
    ok2[k] = closer.close(r, a, resample, {&gpu2, slot});
    a2[k] = a;
  }, [&](std::span<const int> sl) { gpu2.flush(sl); });
  CHECK(ok2 == okG);
  CHECK(a2 == ag);
  CHECK(gpu2.attempts() == gpu.attempts());
  CHECK(gpu2.updates() == gpu.updates());
  CHECK(gpu2.kernels() > gpu.kernels());
}
#endif

TEST_CASE("closure: coupling is transitive (2KJJ A:3-6 + A:8-15)") {
  // Loops A:3-6 (break A:6:C-A:7:N) and A:8-15 (break A:15:C-A:16:N) share no torsion, but each shares
  // torsions with the disulfide A:6-A:11: one group of three links.
  StaFixture f;
  const auto cfg = f.ok(freeResidues(freeResidues(f.text, "A", 3, 6), "A", 8, 15));
  const od::MolSettings& s = cfg.molecules[0];
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &s.tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  REQUIRE(s.tree.links.size() == 4);
  const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), {});
  const auto lines = closer.describe();
  for (const auto& l : lines) MESSAGE(l);
  REQUIRE(lines.size() == 2);
  CHECK(lines[0].find("A:6:C - A:7:N + A:15:C - A:16:N + A:6:SG - A:11:SG: coupled, 3 links") != std::string::npos);
  CHECK(lines[1].find("A:20:SG - B:19:SG") != std::string::npos);
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  for (const double clash : {0.0, 1.8}) {
    od::ClosureParams cp;
    cp.clashDistance = clash;
    const od::LoopCloser c(sys, parts, freeMask, smp.segments(), cp);
    int closed = 0;
    for (int k = 0; k < 40; ++k) {
      od::Rng r(7, 0, k);
      std::vector<double> a = sys.inputAngles();
      smp.sample(r, a);
      if (!c.close(r, a, resample)) continue;
      ++closed;
      const auto x = sys.coordinates(std::span<const double>(a));
      for (int li : {0, 1, 2}) CHECK(linkDeviation(s.tree, s.tree.links[li], x) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
    }
    MESSAGE(std::format("clashDistance {}: closed {} of 40", clash, closed));
    if (clash == 0) CHECK(closed > 0);
  }
}

TEST_CASE("closure: clash term in the RCD objective (2KJJ A:8-15 with the disulfide A:6-A:11)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 8, 15));
  const od::MolSettings& s = cfg.molecules[0];
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &s.tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  int closed[3] = {0, 0, 0};
  for (const od::RcdClash mode : {od::RcdClash::None, od::RcdClash::Repair, od::RcdClash::Full}) {
    od::ClosureParams cp;
    cp.growClose = true;
    cp.rcdMaxIter = 1000;
    cp.rcdClash = mode;
    cp.rcdClashWeight = 10;
    const od::LoopCloser c(sys, parts, freeMask, smp.segments(), cp);
    // full is slow (spec 7.11): fewer samples.
    for (int k = 0; k < (mode == od::RcdClash::Full ? 4 : 20); ++k) {
      od::Rng r(7, 0, k);
      std::vector<double> a = sys.inputAngles();
      smp.sample(r, a);
      if (!c.close(r, a, resample)) continue;
      ++closed[static_cast<int>(mode)];
      const auto x = sys.coordinates(std::span<const double>(a));
      for (int li : {0, 1}) CHECK(linkDeviation(s.tree, s.tree.links[li], x) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
    }
  }
  MESSAGE(std::format("closed: none {} of 20, repair {} of 20, full {} of 4", closed[0], closed[1], closed[2]));
  CHECK(closed[1] >= closed[0]);
}

TEST_CASE("closure: growClose builds the rest after the closure (2KJJ A:12-17)") {
  StaFixture f;
  const auto cfg = f.ok(freeResidues(f.text, "A", 12, 17));
  const od::MolSettings& s = cfg.molecules[0];
  const od::TorsionTree& tree = s.tree;
  const std::vector<od::EnergyPart> parts{{&f.mols[0], &tree}};
  const od::EnergySystem sys(parts, {});
  const od::Csgm1Sampler smp(f.mols[0], s, cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  od::ClosureParams cp;
  cp.growClose = true;
  const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
  const std::vector<bool>& path = closer.pathTorsions();
  REQUIRE(static_cast<int>(path.size()) == sys.torsionCount());
  CHECK(std::count(path.begin(), path.end(), true) >= 12);
  od::GrowParams gp;
  gp.reach = false;
  gp.fixed = path;
  const od::GrowBuilder grow(sys, parts, {&s.topology}, smp.segments(), gp);
  CHECK(grow.fixedStepCount() > 0);
  CHECK(grow.fixedStepCount() < grow.stepCount());
  CHECK(grow.reachCheckCount() == 0);
  int closed = 0, built = 0;
  for (int k = 0; k < 40; ++k) {
    od::Rng r(11, 0, k);
    std::vector<double> a = sys.inputAngles();
    smp.sample(r, a);
    if (!closer.close(r, a, resample)) continue;
    ++closed;
    const std::vector<double> afterClose = a;
    if (!grow.build(r, a, resample).ok) continue;
    ++built;
    // The closure paths keep their closed values; the loop stays closed.
    for (int g = 0; g < sys.torsionCount(); ++g)
      if (path[g]) CHECK(a[g] == afterClose[g]);
    const auto x = sys.coordinates(std::span<const double>(a));
    CHECK(linkDeviation(tree, tree.links[0], x) <= cp.rcdTolerance * std::sqrt(3.0) + 1e-9);
  }
  MESSAGE(std::format("growClose 2KJJ A:12-17: closed {} of 40, built {}", closed, built));
  CHECK(closed > 0);
  CHECK(built > 0);
}

TEST_CASE("closure: metal distance restraints (1ZNF, metalBond = break)") {
  ConfigFixture f;
  f.mols.push_back(zincFinger());
  f.topos.push_back(od::analyzeTopology(f.mols.back(), {}));
  od::assignDssp({&f.mols.back()}, {&f.topos.back()});
  f.inputs.clear();
  for (std::size_t i = 0; i < f.mols.size(); ++i) f.inputs.push_back({&f.mols[i], &f.topos[i]});
  const std::string zn = "id           = \"mobm:1ZNF.cif:A\"";
  std::string text = setQuoted(uncommentMetal(od::defaultMolConfigText(f.inputs)), zn, "metalBond    = ", "break");
  // The torsions of the rings through Zn become rotatable after the break.
  for (std::size_t p; (p = text.find("rotatable = false\n    category  = \"ring\"", text.find(zn))) != std::string::npos;)
    text.replace(p, 17, "rotatable = true ");
  text = phi10(text);
  std::ofstream(f.path) << text;
  std::vector<std::string> errors;
  const auto cfg = od::readMolConfig(f.path, f.inputs, {}, errors);
  for (const auto& e : errors) MESSAGE(e);
  REQUIRE(cfg);
  const od::MolSettings& s = cfg->molecules[3];
  const std::vector<od::EnergyPart> parts{{&f.mols[3], &s.tree}};
  const od::EnergySystem sys(parts, {});
  REQUIRE(s.tree.restraints.size() == 3);
  const od::Csgm1Sampler smp(f.mols[3], s, *cfg, {});
  const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { smp.sampleSegment(r, seg, a); };
  std::vector<bool> freeMask(sys.torsionCount());
  for (int g = 0; g < sys.torsionCount(); ++g) freeMask[g] = !smp.cisTrans(g);
  od::ClosureParams cp;
  cp.method = od::ClosureMethod::Kic;  // not applicable to distance restraints: RCD
  cp.clashDistance = 0;
  const od::LoopCloser closer(sys, parts, freeMask, smp.segments(), cp);
  const auto lines = closer.describe();
  // The three distances share torsions: one best-effort group.
  REQUIRE(lines.size() == 1);
  CHECK(lines[0].find("A:6:SG - A:27:ZN + A:19:NE2 - A:27:ZN + A:23:NE2 - A:27:ZN: coupled, 3 links") != std::string::npos);
  CHECK(lines[0].find("torsions, rcd, best effort (metal distance)") != std::string::npos);
  // Best effort: closure never fails; the deviation does not grow and often reaches the tolerance.
  int reached = 0;
  for (int k = 0; k < 20; ++k) {
    od::Rng r(5, 0, k);
    std::vector<double> a = sys.inputAngles();
    smp.sample(r, a);
    const auto x0 = sys.coordinates(std::span<const double>(a));
    CHECK(closer.close(r, a, resample));
    const auto x = sys.coordinates(std::span<const double>(a));
    for (const auto& l : s.tree.links) {
      CHECK(linkDeviation(s.tree, l, x) <= linkDeviation(s.tree, l, x0) + 1e-9);
      reached += linkDeviation(s.tree, l, x) <= cp.rcdTolerance + 1e-9;
    }
  }
  CHECK(reached >= 20);
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
  od::VecMinParams p;
  p.deltaE = 1e-12;
  p.convSteps = 2;
  std::vector<double> x{5, -3, 7};
  auto r = od::minimizeVec(quad, x, p);
  CHECK_FALSE(r.aborted);
  for (int i = 0; i < 3; ++i) CHECK(x[i] == doctest::Approx(i).scale(1).epsilon(1e-5));
  const od::VecObjective rosen = [](std::span<const double> x, double& e, std::vector<double>& g) {
    e = 100 * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1 - x[0], 2);
    g = {-400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]), 200 * (x[1] - x[0] * x[0])};
    return true;
  };
  x = {-1.2, 1};
  r = od::minimizeVec(rosen, x, p);
  CHECK(x[0] == doctest::Approx(1).epsilon(1e-4));
  CHECK(x[1] == doctest::Approx(1).epsilon(1e-4));
  // Step limit, accept hook and abort.
  int accepted = 0;
  p.stepLimit = [](std::span<const double> d) { return 0.1 / std::max(0.1, std::abs(d[0])); };
  p.onAccept = [&](std::vector<double>&, bool) { ++accepted; };
  x = {5, -3, 7};
  double last = 5;
  const od::VecObjective quad2 = [&](std::span<const double> v, double& e, std::vector<double>& g) {
    CHECK(std::abs(v[0] - last) <= 0.1 + 1e-12);
    return quad(v, e, g);
  };
  p.onAccept = [&](std::vector<double>& v, bool) { ++accepted; last = v[0]; };
  r = od::minimizeVec(quad2, x, p);
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
  CHECK(od::minimizeVec(stop, x, p).aborted);
  x = {3, -3, 7};
  CHECK(od::minimizeVec(stop, x, p).aborted);
}

TEST_CASE("minimize: line search options (lineSearchInit, lineSearchReduce, maxBacktrack)") {
  CHECK(od::parseLineSearchInit("unit") == od::LineSearchInit::Unit);
  CHECK(od::parseLineSearchInit("prev") == od::LineSearchInit::Prev);
  CHECK_THROWS_AS(od::parseLineSearchInit("one"), od::Error);
  CHECK(od::parseLineSearchReduce("halve") == od::LineSearchReduce::Halve);
  CHECK(od::parseLineSearchReduce("interp") == od::LineSearchReduce::Interp);
  CHECK_THROWS_AS(od::parseLineSearchReduce("cubic"), od::Error);
  int evals = 0;
  std::vector<double> trials;  // first coordinate of each evaluated point
  const od::VecObjective rosen = [&](std::span<const double> x, double& e, std::vector<double>& g) {
    ++evals;
    trials.push_back(x[0]);
    e = 100 * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1 - x[0], 2);
    g = {-400 * x[0] * (x[1] - x[0] * x[0]) - 2 * (1 - x[0]), 200 * (x[1] - x[0] * x[0])};
    return true;
  };
  for (od::MinMethod m : {od::MinMethod::Cg, od::MinMethod::Lbfgs})
    for (auto init : {od::LineSearchInit::Unit, od::LineSearchInit::Prev})
      for (auto red : {od::LineSearchReduce::Halve, od::LineSearchReduce::Interp}) {
        INFO(std::string(od::minMethodName(m)), " ", static_cast<int>(init), " ", static_cast<int>(red));
        od::VecMinParams p;
        p.method = m;
        p.maxSteps = 20000;
        p.deltaE = 1e-14;
        p.convSteps = 3;
        p.lineSearchInit = init;
        p.lineSearchReduce = red;
        std::vector<double> x{-1.2, 1};
        const auto r = od::minimizeVec(rosen, x, p);
        CHECK_FALSE(r.aborted);
        if (m == od::MinMethod::Lbfgs) {
          CHECK(x[0] == doctest::Approx(1).epsilon(1e-3));
          CHECK(x[1] == doctest::Approx(1).epsilon(1e-3));
        } else {  // cg stops on the valley floor before the minimum
          CHECK(r.energy < 0.05);
        }
      }
  // Along a quadratic valley: halving gives trial factors 1, 1/2, 1/4, ...; interp jumps to the
  // minimum of the parabola (within [0.1, 0.5] of the factor).
  const od::VecObjective steep = [&](std::span<const double> x, double& e, std::vector<double>& g) {
    ++evals;
    trials.push_back(x[0]);
    e = 50 * x[0] * x[0];
    g = {100 * x[0]};
    return true;
  };
  for (auto red : {od::LineSearchReduce::Halve, od::LineSearchReduce::Interp}) {
    od::VecMinParams p;
    p.method = od::MinMethod::Lbfgs;
    p.maxSteps = 1;
    p.lineSearchReduce = red;
    std::vector<double> x{1};
    trials.clear();
    od::minimizeVec(steep, x, p);
    // Start, then trial steps -100, then -50 (halve) or -10 (interp: the exact minimum is at -1).
    REQUIRE(trials.size() >= 3);
    CHECK(trials[1] == doctest::Approx(-99));
    CHECK(trials[2] == doctest::Approx(red == od::LineSearchReduce::Halve ? -49.0 : -9.0));
  }
  // maxBacktrack: the trial points of a step that is never accepted.
  const od::VecObjective rising = [&](std::span<const double> x, double& e, std::vector<double>& g) {
    ++evals;
    e = x[0] == 0 ? 0 : 1;
    g = {1};
    return true;
  };
  for (int mb : {1, 5, 20}) {
    od::VecMinParams p;
    p.method = od::MinMethod::Lbfgs;
    p.maxBacktrack = mb;
    std::vector<double> x{0};
    const auto r = od::minimizeVec(rising, x, p);
    CHECK(r.evaluations == 1 + mb);
    CHECK(r.steps == 0);
  }
}

TEST_CASE("minimize: all methods (quadratic, step limit, best point, abort)") {
  const double a[3] = {1, 10, 100};
  int evals = 0;
  const od::VecObjective quad = [&](std::span<const double> x, double& e, std::vector<double>& g) {
    ++evals;
    e = 0;
    g.assign(3, 0);
    for (int i = 0; i < 3; ++i) {
      e += 0.5 * a[i] * (x[i] - i) * (x[i] - i);
      g[i] = a[i] * (x[i] - i);
    }
    return true;
  };
  for (const char* name : {"sd", "cg", "fire", "lbfgs", "adadelta"}) {
    const std::string nameS = name;
    CAPTURE(nameS);
    const od::MinMethod m = od::parseMinMethod(name);
    CHECK(std::string(od::minMethodName(m)) == name);
    const bool monotone = m == od::MinMethod::Sd || m == od::MinMethod::Cg || m == od::MinMethod::Lbfgs;
    od::VecMinParams p;
    p.method = m;
    p.maxSteps = 20000;
    p.deltaE = 1e-14;
    p.convSteps = 5;
    // Step limit 0.1 on every variable; accepted energies, best flags.
    std::vector<double> x{5, -3, 7}, prev = x, accE;
    double bestSeen = 1e300;
    int bestFlags = 0, stepErr = 0;
    p.stepLimit = [](std::span<const double> d) {
      double big = 0;
      for (double v : d) big = std::max(big, std::abs(v));
      return big > 0.1 ? 0.1 / big : 1.0;
    };
    p.onAccept = [&](std::vector<double>& v, bool best) {
      for (int i = 0; i < 3; ++i)
        if (std::abs(v[i] - prev[i]) > 0.1 + 1e-12) ++stepErr;
      prev = v;
      double e;
      std::vector<double> g;
      quad(v, e, g);
      accE.push_back(e);
      if (best != (e < bestSeen)) ++bestFlags;
      bestSeen = std::min(bestSeen, e);
    };
    evals = 0;
    auto r = od::minimizeVec(quad, x, p);
    CHECK_FALSE(r.aborted);
    CHECK(stepErr == 0);
    CHECK(bestFlags == 0);
    CHECK(r.steps == static_cast<int>(accE.size()));
    CHECK(r.energy < 1e-6);
    CHECK(r.energy == doctest::Approx(std::min(bestSeen, 0.5 * (25 + 10 * 25 + 100 * 25))));
    double e0;
    std::vector<double> g0;
    quad(x, e0, g0);
    CHECK(e0 == r.energy);
    if (monotone)
      for (std::size_t k = 1; k < accE.size(); ++k) CHECK(accE[k] <= accE[k - 1]);
    // An infeasible point aborts.
    const od::VecObjective stop = [&](std::span<const double> v, double& e, std::vector<double>& g) {
      quad(v, e, g);
      return v[0] > 4;
    };
    p.onAccept = {};
    x = {5, -3, 7};
    CHECK(od::minimizeVec(stop, x, p).aborted);
    // maxSteps.
    p.maxSteps = 3;
    p.stepLimit = {};
    x = {5, -3, 7};
    r = od::minimizeVec(quad, x, p);
    CHECK(r.steps <= 3);
    p.maxSteps = 0;
    x = {5, -3, 7};
    r = od::minimizeVec(quad, x, p);
    CHECK(r.steps == 0);
    CHECK(r.evaluations == 1);
    CHECK(x[0] == 5);
  }
  CHECK_THROWS_AS(od::parseMinMethod("bfgs"), od::Error);
}

TEST_CASE("grid: flexres and delta fields (2KJJ, 1L2Y)") {
  GridFixture x;
  const od::Vec3 c = x.stamCenter();
  const auto geom = od::boxGrid(c, {16, 16, 16}, 0.375);
  const double cap = 6.75;
  CHECK(od::movedGridAtoms(*x.staSys, x.staParts).size() == x.stam.size());
  std::vector<char> mask(x.stam.size());
  for (std::size_t a = 0; a < mask.size(); ++a) mask[a] = a % 4 == 0;
  auto shared = [&](const od::GridAtoms& a) {
    return std::pair{std::make_shared<const od::EnergyGrids>(od::buildEnergyGrids(geom, a, x.mobm.type, {})),
                     std::make_shared<const od::DistanceGrid>(od::buildDistanceGrid(geom, a, cap))};
  };
  const auto [g, d] = shared(x.stam);
  const auto full = od::StamField::full(g, d);
  const od::GridAtoms fixedAtoms = od::subset(x.stam, mask, false), moved = od::subset(x.stam, mask, true);
  CHECK(fixedAtoms.size() + moved.size() == x.stam.size());
  const auto [gf, df] = shared(fixedAtoms);
  const auto flex = od::StamField::flexres(gf, df, moved, {});
  CHECK(flex.mode() == od::GridMode::Flexres);
  CHECK(flex.movedAtomCount() == moved.size());

  auto probeAt = [&](const od::Vec3& p, int m, double q) {
    od::GridAtoms one;
    one.xyz = {p};
    one.type = {g->types[m]};
    one.charge = {q};
    one.heavy = {1};
    return one;
  };
  std::mt19937 rng(11);
  std::uniform_int_distribution<int> ui(2, geom.n[0] - 3);
  std::uniform_real_distribution<double> ur(-0.5, 0.5);
  // At grid points away from clashes both equal the pair sum.
  int checked = 0;
  for (int k = 0; k < 3000 && checked < 100; ++k) {
    const od::Vec3 p = geom.point(ui(rng), ui(rng), ui(rng));
    if (full.distance<double>(p, nullptr) < 2.5) continue;
    ++checked;
    CHECK(flex.distance<double>(p, nullptr) == doctest::Approx(full.distance<double>(p, nullptr)).epsilon(1e-5));
    for (int m = 0; m < static_cast<int>(g->types.size()); ++m)
      for (double q : {-0.4, 0.3})
        CHECK(flex.energy<double>(m, q, p, nullptr) ==
              doctest::Approx(full.energy<double>(m, q, p, nullptr)).epsilon(1e-4).scale(1.0));
  }
  CHECK(checked == 100);
  // Between grid points the pair sum of the moved atoms is exact: closer to the pair reference.
  double errFull = 0, errFlex = 0;
  checked = 0;
  for (int k = 0; k < 3000 && checked < 200; ++k) {
    const auto q0 = geom.point(ui(rng), ui(rng), ui(rng));
    const od::Vec3 p{q0[0] + 0.375 * ur(rng), q0[1] + 0.375 * ur(rng), q0[2] + 0.375 * ur(rng)};
    if (full.distance<double>(p, nullptr) < 2.5) continue;
    ++checked;
    const int m = checked % static_cast<int>(g->types.size());
    const double ref = od::pairInterEnergy(x.stam, probeAt(p, m, 0.3), {});
    errFull += std::abs(full.energy<double>(m, 0.3, p, nullptr) - ref);
    errFlex += std::abs(flex.energy<double>(m, 0.3, p, nullptr) - ref);
    // Gradients (finite differences) of the flexres energy and distance.
    od::Vec3 ge, gd;
    flex.energy<double>(m, 0.3, p, &ge);
    flex.distance<double>(p, &gd);
    const double h = 1e-5;
    for (int i = 0; i < 3; ++i) {
      od::Vec3 a = p, b = p;
      a[i] += h;
      b[i] -= h;
      const double fe = (flex.energy<double>(m, 0.3, a, nullptr) - flex.energy<double>(m, 0.3, b, nullptr)) / (2 * h);
      CHECK(ge[i] == doctest::Approx(fe).epsilon(1e-3).scale(1e-2));
      const double fd = (flex.distance<double>(a, nullptr) - flex.distance<double>(b, nullptr)) / (2 * h);
      CHECK(gd[i] == doctest::Approx(fd).epsilon(1e-3).scale(1e-3));
    }
  }
  MESSAGE(std::format("mean |E - pair sum| between grid points: full {:.4f}, flexres {:.4f}", errFull / checked, errFlex / checked));
  CHECK(errFlex < errFull);
  // Float evaluation follows the double one.
  const od::Vec3f pf{static_cast<float>(c[0] + 0.1), static_cast<float>(c[1] + 0.2), static_cast<float>(c[2] + 7.3)};
  const od::Vec3 pd{pf[0], pf[1], pf[2]};
  CHECK(flex.energy<float>(0, 0.3f, pf, nullptr) == doctest::Approx(flex.energy<double>(0, 0.3, pd, nullptr)).epsilon(1e-4).scale(1e-3));

  // delta: moved atoms on one side shifted by 0.7 Å.
  od::GridAtoms sample = x.stam;
  for (std::size_t a = 0; a < sample.size(); ++a)
    if (mask[a] && sample.xyz[a][0] > c[0] + 4) sample.xyz[a][0] += 0.7;
  const auto [gs, ds] = shared(sample);
  const auto fullS = od::StamField::full(gs, ds);
  const auto d0 = od::StamField::delta(g, d, *gs, *ds, 0.0);
  const auto d1 = od::StamField::delta(g, d, *gs, *ds, 1e-2);
  CHECK(d0.blockCount() == 6 * 6 * 6);  // 45 points per axis, blocks of 8
  CHECK(d1.deltaBlockCount() > 0);
  CHECK(d1.deltaBlockCount() < d0.deltaBlockCount());
  // Electrostatics reach 20.48 Å: in a 16 Å box nearly every block differs by more than 1e-4 (spec 8.1).
  MESSAGE(std::format("delta blocks stored: threshold 0: {} of {}; 1e-2: {}", d0.deltaBlockCount(), d0.blockCount(),
                      d1.deltaBlockCount()));
  for (int k = 0; k < 300; ++k) {
    const auto q0 = geom.point(ui(rng), ui(rng), ui(rng));
    const od::Vec3f p{static_cast<float>(q0[0] + 0.375 * ur(rng)), static_cast<float>(q0[1] + 0.375 * ur(rng)),
                      static_cast<float>(q0[2] + 0.375 * ur(rng))};
    const int m = k % static_cast<int>(g->types.size());
    od::Vec3f a, b;
    CHECK(d0.energy<float>(m, -0.35f, p, &a) == fullS.energy<float>(m, -0.35f, p, &b));
    CHECK(a == b);
    CHECK(d0.distance<float>(p, nullptr) == fullS.distance<float>(p, nullptr));
    CHECK(std::abs(d1.energy<float>(m, -0.35f, p, nullptr) - fullS.energy<float>(m, -0.35f, p, nullptr)) <= 1e-2f * 1.7f + 1e-4f);
  }
}

TEST_CASE("dock9b: flexres pose energy and gradient (2KJJ, 1L2Y)") {
  GridFixture x;
  od::Dock9bParams p;
  p.boxCenter = x.stamCenter();
  p.boxSize = {48, 48, 48};
  const auto geom = od::boxGrid(p.boxCenter, p.boxSize, 0.375);
  const double cap = std::max(p.maxNearDistance, p.clashDistance) + 0.75;
  std::vector<char> mask(x.stam.size());
  for (std::size_t a = 0; a < mask.size(); ++a) mask[a] = a % 3 == 0;
  const od::GridAtoms fixedAtoms = od::subset(x.stam, mask, false);
  const auto flex = od::StamField::flexres(std::make_shared<const od::EnergyGrids>(od::buildEnergyGrids(geom, fixedAtoms, x.mobm.type, {})),
                                           std::make_shared<const od::DistanceGrid>(od::buildDistanceGrid(geom, fixedAtoms, cap)),
                                           od::subset(x.stam, mask, true), {});
  const od::MobmDocker dk(*x.mobSys, flex, p);
  const auto angles = x.mobSys->inputAngles();
  const auto c = dk.conformation(angles);
  od::Pose pose{{p.boxCenter[0] + 9, p.boxCenter[1] + 2, p.boxCenter[2] - 3}, od::normalized(od::Quat{0.9, 0.1, -0.3, 0.2})};
  std::vector<double> g;
  const auto pe = dk.evaluateDouble(c, pose, angles, &g);
  CHECK(pe.outside == 0);
  CHECK(dk.evaluate(c, pose, angles).objective() == doctest::Approx(pe.objective()).epsilon(1e-4).scale(1.0));
  const double h = 1e-5;
  for (int i = 0; i < 3; ++i) {
    od::Pose a = pose, b = pose;
    a.t[i] += h;
    b.t[i] -= h;
    const double fd = (dk.evaluateDouble(c, a, angles).objective() - dk.evaluateDouble(c, b, angles).objective()) / (2 * h);
    CHECK(g[i] == doctest::Approx(fd).epsilon(1e-3).scale(1e-2));
  }
}

namespace {

struct DockFixture {
  GridFixture x;
  od::Dock9bParams p;
  od::EnergyGrids grids;
  od::DistanceGrid dist;
  od::StamField field;
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
    field = od::StamField::full(std::make_shared<const od::EnergyGrids>(grids), std::make_shared<const od::DistanceGrid>(dist));
  }
};

}  // namespace

namespace {

// full, flexres and delta fields of 2KJJ in a 48 Å box, with 1L2Y as MobM (spec 15.1.4 tests).
struct FieldSet {
  GridFixture x;
  od::Dock9bParams p;
  std::vector<std::pair<std::string, od::StamField>> fields;
  FieldSet() {
    p.boxCenter = x.stamCenter();
    p.boxSize = {48, 48, 48};
    p.hopfLevel = 0;
    p.posm1SampleNum = 2000;
    p.seed = 3;
    const auto geom = od::boxGrid(p.boxCenter, p.boxSize, 0.375);
    const double cap = std::max(p.maxNearDistance, p.clashDistance) + 0.75;
    auto shared = [&](const od::GridAtoms& a) {
      return std::pair{std::make_shared<const od::EnergyGrids>(od::buildEnergyGrids(geom, a, x.mobm.type, {})),
                       std::make_shared<const od::DistanceGrid>(od::buildDistanceGrid(geom, a, cap))};
    };
    std::vector<char> mask(x.stam.size());
    for (std::size_t a = 0; a < mask.size(); ++a) mask[a] = a % 3 == 0;
    const auto [g, d] = shared(x.stam);
    fields.emplace_back("full", od::StamField::full(g, d));
    const auto [gf, df] = shared(od::subset(x.stam, mask, false));
    fields.emplace_back("flexres", od::StamField::flexres(gf, df, od::subset(x.stam, mask, true), {}));
    od::GridAtoms sample = x.stam;
    const od::Vec3 c = x.stamCenter();
    for (std::size_t a = 0; a < sample.size(); ++a)
      if (mask[a] && sample.xyz[a][0] > c[0] + 4) sample.xyz[a][0] += 0.7;
    const auto [gs, ds] = shared(sample);
    fields.emplace_back("delta", od::StamField::delta(g, d, *gs, *ds, 1e-2));
  }
};

// PoseEvaluator on the CPU (MobmDocker::evaluate), recording the batch sizes; with nLanes > 1
// the asynchronous form (submit evaluates, complete returns the stored results).
struct CpuPoseEvaluator : od::PoseEvaluator {
  const od::MobmDocker* dk = nullptr;
  std::vector<const od::MobmDocker*> dks;  // bindAll; PoseRequest::field
  const od::MobmDocker::Conf* conf = nullptr;
  int batch, nLanes, largest = 0, binds = 0, nFields = 1;
  struct Done {
    bool pending = false;
    std::vector<od::MobmDocker::PoseEnergy> out;
    std::vector<double> grad;
  };
  std::vector<Done> done;
  explicit CpuPoseEvaluator(int b, int l = 1) : batch(b), nLanes(l), done(l) {}
  int maxBatch() const override { return batch; }
  int lanes() const override { return nLanes; }
  void submit(int lane, int n, const od::PoseRequest* req) override {
    Done& d = done[lane];
    REQUIRE_FALSE(d.pending);
    d.pending = true;
    d.out.resize(n);
    d.grad.resize(n * (6 + dk->sys().torsionCount()));
    evaluate(n, req, d.out.data(), d.grad.data());
  }
  void complete(int lane, od::MobmDocker::PoseEnergy* out, double* grad) override {
    Done& d = done[lane];
    REQUIRE(d.pending);
    d.pending = false;
    std::copy(d.out.begin(), d.out.end(), out);
    std::copy(d.grad.begin(), d.grad.end(), grad);
  }
  void bind(const od::MobmDocker& d) override {
    dk = &d;
    dks.assign(1, &d);
    ++binds;
  }
  int fields() const override { return nFields; }
  void bindAll(std::span<const od::MobmDocker* const> d) override {
    REQUIRE(static_cast<int>(d.size()) <= nFields);
    dk = d[0];
    dks.assign(d.begin(), d.end());
    ++binds;
  }
  void evaluate(int n, const od::PoseRequest* req, od::MobmDocker::PoseEnergy* out, double* grad) override {
    largest = std::max(largest, n);
    const std::size_t ng = 6 + dk->sys().torsionCount();
    std::vector<double> g;
    for (int i = 0; i < n; ++i) {
      od::MobmDocker::Conf c;
      c.center = req[i].center;
      const std::span<const double> a(req[i].anglesDeg, dk->sys().torsionCount());
      out[i] = dks[static_cast<std::size_t>(req[i].field)]->evaluate(c, req[i].pose, a, &g);
      std::copy(g.begin(), g.end(), grad + i * ng);
    }
  }
};

// Poses near the StaM surface: quick-filtered L1poss of the input conformation.
std::vector<od::Pose> filteredPoses(const od::MobmDocker& dk, const od::MobmDocker::Conf& c, std::size_t want) {
  std::vector<od::Pose> v;
  od::Dock9bStats st;
  for (long long k = 0; k < dk.l1Count(c) && v.size() < want; ++k) {
    const auto pose = dk.l1Pose(c, 7, k);
    if (pose && dk.quickFilter(c, *pose, st)) v.push_back(*pose);
  }
  return v;
}

}  // namespace

TEST_CASE("dock9b: flat field and pose energy equal StamField and MobmDocker::evaluate") {
  FieldSet fs;
  std::mt19937 rng(17);
  const auto geom = od::boxGrid(fs.p.boxCenter, fs.p.boxSize, 0.375);
  std::uniform_real_distribution<double> ux(geom.origin[0] + 0.1, geom.max()[0] - 0.1), uy(geom.origin[1] + 0.1, geom.max()[1] - 0.1),
      uz(geom.origin[2] + 0.1, geom.max()[2] - 0.1);
  const od::EnergyModel mob = fs.x.mobSys->model();
  for (const auto& [name, field] : fs.fields) {
    INFO(name);
    const od::FieldModel fm = field.model();
    const od::FieldView& v = fm.view;
    int diff = 0;
    for (int k = 0; k < 400; ++k) {
      const od::Vec3 pd{ux(rng), uy(rng), uz(rng)};
      const od::Vec3f pf{static_cast<float>(pd[0]), static_cast<float>(pd[1]), static_cast<float>(pd[2])};
      const int m = k % field.typeMapCount();
      const float q = k % 3 == 0 ? -0.35f : 0.42f;
      od::Vec3f ga, gb;
      od::Vec3 da, db;
      diff += field.energy<float>(m, q, pf, &ga) != od::field::energy<float>(v, m, q, pf.data(), gb.data()) || ga != gb;
      diff += field.distance<float>(pf, &ga) != od::field::distance<float>(v, pf.data(), gb.data()) || ga != gb;
      diff += field.energy<double>(m, q, pd, &da) != od::field::energy<double>(v, m, q, pd.data(), db.data()) || da != db;
      diff += field.distance<double>(pd, &da) != od::field::distance<double>(v, pd.data(), db.data()) || da != db;
    }
    CHECK(diff == 0);
    const od::MobmDocker dk(*fs.x.mobSys, field, fs.p);
    const auto probes = dk.probes();
    const od::DockView dv = dk.dockView(probes, v);
    const auto start = perturbed(fs.x.mobSys->inputAngles(), 20, 3);
    std::vector<double> angles(start.size());
    for (std::size_t i = 0; i < start.size(); ++i) angles[i] = static_cast<float>(start[i]);
    const auto c = dk.conformation(fs.x.mobSys->inputAngles());
    const auto poses = filteredPoses(dk, c, 40);
    REQUIRE(poses.size() == 40);
    std::vector<float> af(angles.begin(), angles.end()), x(3 * mob.coordCount);
    std::vector<double> g(3 * mob.coordCount), grad(6 + angles.size()), ref;
    int bad = 0, outside = 0;
    for (const auto& pose : poses) {
      const auto pe = dk.evaluate(c, pose, angles, &ref);
      od::DockPose ps{};
      for (int i = 0; i < 3; ++i) {
        ps.t[i] = static_cast<float>(pose.t[i]);
        ps.center[i] = static_cast<float>(c.center[i]);
      }
      ps.q[0] = static_cast<float>(pose.q.w);
      ps.q[1] = static_cast<float>(pose.q.x);
      ps.q[2] = static_cast<float>(pose.q.y);
      ps.q[3] = static_cast<float>(pose.q.z);
      od::DockTerms dt;
      const od::FlatTerms ft = od::dock::evaluate(mob.view(), dv, af.data(), ps, x.data(), g.data(), grad.data(), dt);
      outside += pe.outside;
      bad += dt.inter != pe.inter || dt.interClash != pe.interClash || dt.outside != pe.outside || ft.vdw != pe.intra.vdw ||
             ft.elec != pe.intra.elec || ft.desolv != pe.intra.desolv || ft.clash != pe.intra.clash || grad != ref;
    }
    CHECK(bad == 0);
    CHECK(outside < 40 * static_cast<int>(probes.size()));  // probe atoms outside the grid are compared too
  }
}

TEST_CASE("dock9b: batched minimization with a PoseEvaluator equals the CPU path") {
  FieldSet fs;
  fs.p.posm1SampleNum = 200;
  for (const auto& [name, field] : fs.fields) {
    INFO(name);
    std::vector<od::DockSample> a, b;
    const od::MobmDocker dk(*fs.x.mobSys, field, fs.p);
    const auto sa = dk.dockPair(fs.x.mobSys->inputAngles(), 5, 16, [&](std::vector<od::DockSample>& ch) {
      a.insert(a.end(), ch.begin(), ch.end());
    });
    for (int lanes : {1, 2, 3}) {
      INFO("lanes ", lanes);
      b.clear();
      CpuPoseEvaluator ev(7, lanes);
      od::Dock9bParams pe = fs.p;
      pe.evaluator = &ev;
      const od::MobmDocker dke(*fs.x.mobSys, field, pe);
      const auto sb = dke.dockPair(fs.x.mobSys->inputAngles(), 5, 16, [&](std::vector<od::DockSample>& ch) {
        b.insert(b.end(), ch.begin(), ch.end());
      });
      CHECK(sa.l2 > 0);
      CHECK(sa.l2 == sb.l2);
      CHECK(sa.steps == sb.steps);
      CHECK(sa.evaluations == sb.evaluations);
      CHECK(sa.leftGrid == sb.leftGrid);
      CHECK(ev.largest == (7 + lanes - 1) / lanes);
      CHECK(ev.binds == static_cast<int>((sa.l2 + 15) / 16));
      REQUIRE(a.size() == b.size());
      int diff = 0;
      for (std::size_t i = 0; i < a.size(); ++i)
        diff += a[i].energy != b[i].energy || a[i].angles != b[i].angles || a[i].pose.t != b[i].pose.t || a[i].steps != b[i].steps;
      CHECK(diff == 0);
    }
  }
}

TEST_CASE("dock9b: dockPairs delivers the chunks of sequential dockPair calls") {
  FieldSet fs;
  fs.p.posm1SampleNum = 120;
  const od::StamField& field = fs.fields[0].second;
  std::vector<std::vector<double>> angles;
  for (unsigned j = 0; j < 4; ++j) angles.push_back(perturbed(fs.x.mobSys->inputAngles(), j == 2 ? 0 : 10, 60 + j));
  // (job, chunk) sequence: sizes and energies.
  using Seq = std::vector<std::pair<int, std::vector<float>>>;
  auto record = [](Seq& s, int j) {
    return [&s, j](std::vector<od::DockSample>& ch) {
      std::vector<float> e;
      for (const auto& d : ch) e.push_back(d.energy);
      s.emplace_back(j, std::move(e));
    };
  };
  Seq ref, got;
  std::vector<od::Dock9bStats> refStats;
  const od::MobmDocker dk(*fs.x.mobSys, field, fs.p);
  for (int j = 0; j < 4; ++j) refStats.push_back(dk.dockPair(angles[j], 11 + j, 5, record(ref, j)));
  for (int lanes : {1, 2}) {
    INFO("lanes ", lanes);
    got.clear();
    CpuPoseEvaluator ev(9, lanes);
    od::Dock9bParams pe = fs.p;
    pe.evaluator = &ev;
    const od::MobmDocker dke(*fs.x.mobSys, field, pe);
    std::vector<od::MobmDocker::PairJob> jobs(4);
    for (int j = 0; j < 4; ++j) {
      jobs[j].mobAnglesDeg = angles[j];
      jobs[j].stream = 11 + j;
      jobs[j].sink = record(got, j);
    }
    dke.dockPairs(jobs, 5);
    CHECK(ref.size() > 4);
    CHECK(got == ref);
    for (int j = 0; j < 4; ++j) {
      CHECK(jobs[j].stats.l2 == refStats[j].l2);
      CHECK(jobs[j].stats.steps == refStats[j].steps);
      CHECK(jobs[j].stats.leftGrid == refStats[j].leftGrid);
      CHECK(jobs[j].stats.overlap == refStats[j].overlap);
    }
    CHECK(ev.largest == (9 + lanes - 1) / lanes);
  }
}

TEST_CASE("dock9b: dockPairs over several StaM fields delivers the chunks of sequential calls") {
  FieldSet fs;
  fs.p.posm1SampleNum = 120;
  std::vector<std::vector<double>> angles;
  for (unsigned j = 0; j < 3; ++j) angles.push_back(perturbed(fs.x.mobSys->inputAngles(), j == 1 ? 0 : 10, 70 + j));
  // (field, job, chunk) sequence: sizes and energies.
  using Seq = std::vector<std::tuple<int, int, std::vector<float>>>;
  auto record = [](Seq& s, int f, int j) {
    return [&s, f, j](std::vector<od::DockSample>& ch) {
      std::vector<float> e;
      for (const auto& d : ch) e.push_back(d.energy);
      s.emplace_back(f, j, std::move(e));
    };
  };
  const int nf = static_cast<int>(fs.fields.size());
  Seq ref, got;
  std::vector<od::Dock9bStats> refStats;
  for (int f = 0; f < nf; ++f) {
    const od::MobmDocker dk(*fs.x.mobSys, fs.fields[f].second, fs.p);
    for (int j = 0; j < 3; ++j) refStats.push_back(dk.dockPair(angles[j], 21 + j + 3 * f, 4, record(ref, f, j)));
  }
  // lanes 0: without an evaluator; prepared: the jobs prepared ahead with prepareJobs.
  for (int lanes : {0, 1, 2})
  for (bool prepared : {false, true}) {
    INFO("lanes ", lanes, " prepared ", prepared);
    got.clear();
    CpuPoseEvaluator ev(11, std::max(lanes, 1));
    ev.nFields = nf;
    od::Dock9bParams pe = fs.p;
    pe.evaluator = lanes > 0 ? &ev : nullptr;
    std::vector<std::unique_ptr<od::MobmDocker>> dks;
    std::vector<std::vector<od::MobmDocker::PairJob>> jobs(nf, std::vector<od::MobmDocker::PairJob>(3));
    std::vector<od::MobmDocker::PairGroup> groups;
    for (int f = 0; f < nf; ++f) {
      dks.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, fs.fields[f].second, pe));
      for (int j = 0; j < 3; ++j) {
        jobs[f][j].mobAnglesDeg = angles[j];
        jobs[f][j].stream = 21 + j + 3 * f;
        jobs[f][j].sink = record(got, f, j);
      }
      groups.push_back({dks.back().get(), &jobs[f]});
      if (prepared) {
        dks.back()->prepareJobs(jobs[f]);
        for (int j = 0; j < 3; ++j) CHECK(jobs[f][j].stats.l2 == static_cast<long long>(jobs[f][j].prepared->l2.size()));
      }
    }
    od::MobmDocker::dockPairs(groups, 4);
    CHECK(ref.size() > 9);
    CHECK(got == ref);
    for (int f = 0; f < nf; ++f)
      for (int j = 0; j < 3; ++j) CHECK_FALSE(jobs[f][j].prepared);
    for (int f = 0; f < nf; ++f)
      for (int j = 0; j < 3; ++j) {
        CHECK(jobs[f][j].stats.l2 == refStats[3 * f + j].l2);
        CHECK(jobs[f][j].stats.steps == refStats[3 * f + j].steps);
        CHECK(jobs[f][j].stats.leftGrid == refStats[3 * f + j].leftGrid);
      }
    if (lanes == 0) continue;
    CHECK(ev.largest == (11 + lanes - 1) / lanes);
    // More groups than the evaluator binds together.
    ev.nFields = nf - 1;
    CHECK_THROWS_AS(od::MobmDocker::dockPairs(groups, 4), od::Error);
  }
}

TEST_CASE("sysmem: cgroup memory limit, resident memory and MemAvailable") {
  namespace fs = std::filesystem;
  const fs::path root = fs::temp_directory_path() / std::format("od_test_cgroup_{}", ::getpid());
  fs::remove_all(root);
  const auto put = [](const fs::path& f, const std::string& s) {
    fs::create_directories(f.parent_path());
    std::ofstream(f) << s;
  };
  const std::string proc = (root / "proc").string(), cg = (root / "cg").string();
  std::string src;
  // cgroup v2: the smallest limit of the cgroup and its ancestors.
  put(root / "proc/self/cgroup", "0::/job_1/step_0/task_0\n");
  put(root / "cg/job_1/step_0/task_0/memory.max", "max\n");
  put(root / "cg/job_1/step_0/memory.max", "3221225472\n");
  put(root / "cg/job_1/memory.max", "4294967296\n");
  auto v = od::cgroupMemLimit(&src, proc, cg);
  REQUIRE(v);
  CHECK(*v == 3221225472ULL);
  CHECK(src == cg + "/job_1/step_0/memory.max");
  // No limit anywhere.
  put(root / "cg/job_1/step_0/memory.max", "max\n");
  put(root / "cg/job_1/memory.max", "max\n");
  CHECK_FALSE(od::cgroupMemLimit(nullptr, proc, cg));
  // cgroup v1: the memory controller; the huge value means unlimited.
  put(root / "proc/self/cgroup", "5:cpuacct,cpu:/x\n4:memory:/slurm/uid_1/job_2\n0::/\n");
  put(root / "cg/memory/slurm/uid_1/job_2/memory.limit_in_bytes", "9223372036854771712\n");
  put(root / "cg/memory/slurm/uid_1/memory.limit_in_bytes", "2147483648\n");
  v = od::cgroupMemLimit(&src, proc, cg);
  REQUIRE(v);
  CHECK(*v == 2147483648ULL);
  CHECK(src == cg + "/memory/slurm/uid_1/memory.limit_in_bytes");
  // No cgroup file.
  CHECK_FALSE(od::cgroupMemLimit(nullptr, (root / "none").string(), cg));
  put(root / "proc/meminfo", "MemTotal:       197426540 kB\nMemFree:        1000 kB\nMemAvailable:   185812345 kB\n");
  CHECK(od::memAvailableBytes(proc) == 185812345ULL * 1024);
  CHECK(od::memAvailableBytes((root / "none").string()) == 0);
  // This process.
  CHECK(od::residentBytes() > 0);
  CHECK(od::hostMemLimit(2048).bytes == 2048ULL << 20);
  CHECK(od::hostMemLimit(2048).source == "memory.hostLimit");
  CHECK(od::hostMemLimit(0).bytes > od::residentBytes());
  fs::remove_all(root);
}

TEST_CASE("sysmem: cgroup CPU quota, CPU order by core, CPU list text, cpu.threadNum limit") {
  namespace fs = std::filesystem;
  const fs::path root = fs::temp_directory_path() / std::format("od_test_cgcpu_{}", ::getpid());
  fs::remove_all(root);
  const auto put = [](const fs::path& f, const std::string& s) {
    fs::create_directories(f.parent_path());
    std::ofstream(f) << s;
  };
  const std::string proc = (root / "proc").string(), cg = (root / "cg").string();
  std::string src;
  // cgroup v2: the smallest quota of the cgroup and its ancestors, rounded down, at least 1.
  put(root / "proc/self/cgroup", "0::/job_1/step_0\n");
  put(root / "cg/job_1/step_0/cpu.max", "max 100000\n");
  put(root / "cg/job_1/cpu.max", "650000 100000\n");
  put(root / "cg/cpu.max", "1600000 100000\n");
  auto v = od::cgroupCpuLimit(&src, proc, cg);
  REQUIRE(v);
  CHECK(*v == 6);
  CHECK(src == cg + "/job_1/cpu.max");
  put(root / "cg/job_1/cpu.max", "50000 100000\n");
  CHECK(od::cgroupCpuLimit(nullptr, proc, cg) == 1);
  put(root / "cg/job_1/cpu.max", "max 100000\n");
  put(root / "cg/cpu.max", "max 100000\n");
  CHECK_FALSE(od::cgroupCpuLimit(nullptr, proc, cg));
  // cgroup v1: the cpu controller; quota -1 means unlimited.
  put(root / "proc/self/cgroup", "5:cpuacct,cpu:/slurm/job_2\n4:memory:/slurm/job_2\n0::/\n");
  put(root / "cg/cpu/slurm/job_2/cpu.cfs_quota_us", "-1\n");
  put(root / "cg/cpu/slurm/job_2/cpu.cfs_period_us", "100000\n");
  put(root / "cg/cpu/slurm/cpu.cfs_quota_us", "400000\n");
  put(root / "cg/cpu/slurm/cpu.cfs_period_us", "100000\n");
  v = od::cgroupCpuLimit(&src, proc, cg);
  REQUIRE(v);
  CHECK(*v == 4);
  CHECK(src == cg + "/cpu/slurm/cpu.cfs_quota_us");
  CHECK_FALSE(od::cgroupCpuLimit(nullptr, (root / "none").string(), cg));
  // 4 cores with 2 threads each (cpu c and c + 4): one thread of every core first.
  for (int c = 0; c < 8; ++c)
    put(root / std::format("sys/cpu{}/topology/thread_siblings_list", c), std::format("{},{}\n", c % 4, c % 4 + 4));
  const std::string sys = (root / "sys").string();
  CHECK(od::orderCpusByCore({0, 1, 2, 3, 4, 5, 6, 7}, sys) == std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7});
  CHECK(od::orderCpusByCore({1, 2, 5, 6, 7}, sys) == std::vector<int>{1, 2, 7, 5, 6});
  // Adjacent threads of a core (cpu 2k and 2k + 1), as ranges.
  for (int c = 0; c < 8; ++c)
    put(root / std::format("sys/cpu{}/topology/thread_siblings_list", c), std::format("{}-{}\n", c / 2 * 2, c / 2 * 2 + 1));
  CHECK(od::orderCpusByCore({0, 1, 2, 3, 4, 5, 6, 7}, sys) == std::vector<int>{0, 2, 4, 6, 1, 3, 5, 7});
  CHECK(od::cpuListText({0, 1, 2, 3, 8, 10, 11}) == "0-3,8,10-11");
  CHECK(od::cpuListText({}) == "unknown");
  fs::remove_all(root);
  // cpu.threadNum above the CPUs this process may use: the error names the limit and its sources.
  od::CpuLimit l;
  l.cpus = {0, 1, 2, 3, 4, 5, 6, 7};
  l.count = 8;
  CHECK(od::cpuThreadsError(8, l).empty());
  CHECK(od::cpuThreadsError(9, l) ==
        "cpu.threadNum 9 is more than the 8 CPU(s) this process may use (affinity mask: 8 CPU(s) 0-7; no cgroup CPU quota)");
  l.quota = 3;
  l.quotaSource = "/sys/fs/cgroup/cpu.max";
  l.count = 3;
  CHECK(od::cpuThreadsError(3, l).empty());
  CHECK(od::cpuThreadsError(4, l) == "cpu.threadNum 4 is more than the 3 CPU(s) this process may use (affinity mask: 8 CPU(s) "
                                     "0-7; cgroup CPU quota 3 (/sys/fs/cgroup/cpu.max))");
  CHECK(od::cpuThreadsError(100, od::CpuLimit{}).empty());
  const od::CpuLimit real = od::cpuLimit();
  CHECK(real.count > 0);
  CHECK_THROWS_WITH_AS(od::applyCpuThreads(real.count + 1), doctest::Contains("CPU(s) this process may use"), od::Error);
}

TEST_CASE("runBatchedAsync: every task completes; submit and complete alternate per group") {
  for (int groups : {1, 2, 3}) {
    INFO("groups ", groups);
    constexpr int kCount = 50, kActive = 8;
    std::vector<int> req(kActive), res(kActive), out(kCount, -1);
    std::vector<int> pending(groups, 0), submits(groups, 0);
    std::vector<std::vector<int>> inFlight(groups);
    bool ok = true;
    od::runBatchedAsync(
        kCount, kActive, groups,
        [&](int i, int slot) {
          int v = i;
          for (int k = 0; k < i % 5; ++k) {
            req[slot] = v;
            od::batchWait();
            v = res[slot];
          }
          out[i] = v;
        },
        [&](int g, std::span<const int> slots) {
          ok = ok && !pending[g];
          pending[g] = 1;
          ++submits[g];
          inFlight[g].clear();
          for (int s : slots) {
            ok = ok && kActive * g / groups <= s && s < kActive * (g + 1) / groups;  // contiguous groups
            inFlight[g].push_back(req[s] + 1000);
          }
          inFlight[g].insert(inFlight[g].begin(), slots.begin(), slots.end());
        },
        [&](int g) {
          ok = ok && pending[g];
          pending[g] = 0;
          const std::size_t n = inFlight[g].size() / 2;
          for (std::size_t k = 0; k < n; ++k) res[inFlight[g][k]] = inFlight[g][n + k];
        });
    CHECK(ok);
    for (int g = 0; g < groups; ++g) CHECK(pending[g] == 0);
    for (int i = 0; i < kCount; ++i) CHECK(out[i] == i + 1000 * (i % 5));
    if (groups > 1) CHECK(submits[1] > 0);
  }
  // A failing task: pending evaluations are completed before the exception is rethrown.
  int pendingCount = 0;
  CHECK_THROWS_AS(od::runBatchedAsync(
                      20, 6, 2,
                      [&](int i, int) {
                        od::batchWait();
                        if (i == 7) throw od::Error("task failed");
                        od::batchWait();
                      },
                      [&](int, std::span<const int>) { ++pendingCount; }, [&](int) { --pendingCount; }),
                  od::Error);
  CHECK(pendingCount == 0);
}

TEST_CASE("dock9b: posm3 regions and poses (2KJJ, 1L2Y)") {
  DockFixture f(od::Posm::Posm3);
  CHECK(f.p.sepaPosStep == 2.0);
  CHECK(f.p.sepaPosRange == 2.0);
  CHECK(f.p.sepaTriple);
  CHECK(f.p.sepaPair);
  CHECK_FALSE(f.p.sepaSingle);
  od::Dock9bParams pDef = f.p;
  f.p.sepaPosStep = 1.0;
  f.p.sepaPosRange = 4.0;
  f.p.sepaSingle = true;
  f.p.sepaTriple = false;
  od::StamSepaInfo info;
  const auto sepas = od::stamSepas(f.field, f.p, &info);
  const auto sepas26 = od::stamSepas(f.field, f.p, nullptr, 26);
  MESSAGE(std::format("posm3: shell points {}, selected {}, cut {:.3f}; regions {} (6 neighbours), {} (26 neighbours)",
                      info.shellPoints, info.selected, info.potCut, sepas.size(), sepas26.size()));
  REQUIRE(!sepas.empty());
  CHECK(sepas26.size() <= sepas.size());
  CHECK(info.selected == doctest::Approx(0.05 * static_cast<double>(info.shellPoints)).epsilon(0.02));
  long long pts = 0;
  for (const auto& s : sepas) {
    pts += s.points;
    CHECK(std::abs(s.sign) == 1);
    const double phi = 0.5 * (f.field.energy<double>(0, 1.0, s.xyz, nullptr) - f.field.energy<double>(0, -1.0, s.xyz, nullptr));
    if (s.points >= 8) CHECK(phi * s.sign > 0);  // a compact region keeps its sign at the centroid
  }
  CHECK(pts == info.selected);
  // Threshold replaces the percentile.
  od::Dock9bParams pt = f.p;
  pt.sepaPotThreshold = info.potCut;
  od::StamSepaInfo info2;
  CHECK(od::stamSepas(f.field, pt, &info2).size() == sepas.size());
  CHECK(info2.selected == info.selected);

  f.p.hopfLevel = 0;
  const od::MobmDocker dk(*f.x.mobSys, f.field, f.p, &sepas);
  const auto c = dk.conformation(f.x.mobSys->inputAngles());
  MESSAGE(std::format("posm3: matched pairs {}, single placements {}, L1poss {}", c.pairs.size(), c.singles.size(), dk.l1Count(c)));
  REQUIRE(!c.pairs.empty());
  REQUIRE(!c.singles.empty());
  CHECK(dk.l1Count(c) == static_cast<long long>(c.pairs.size()) * 48 * 729 + static_cast<long long>(c.singles.size()) * 72 * 729);
  const auto& mt = f.x.mobSys->typing(0);
  std::vector<int> probeOf;  // probe index -> coordinate index
  for (std::size_t a = 0; a < mt.type.size(); ++a)
    if (mt.type[a] != od::Ad4Type::Merged) probeOf.push_back(f.x.mobSys->coordOffset(0) + static_cast<int>(a));
  auto world = [&](const od::Pose& pose, int probe) {
    const od::Vec3f& y = c.coords[probeOf[probe]];
    const od::Vec3 r = od::rotate(pose.q, {y[0] - c.center[0], y[1] - c.center[1], y[2] - c.center[2]});
    return od::Vec3{r[0] + pose.t[0], r[1] + pose.t[1], r[2] + pose.t[2]};
  };
  int checkedPairs = 0, checkedSingles = 0;
  std::vector<long long> ks;
  for (long long i = 0; i < 2000; ++i) {
    ks.push_back(i * (c.pairPoses / 2000) + i % 729);
    ks.push_back(c.pairPoses + i * (c.singlePoses / 2000) + (i * 7) % 729);
  }
  for (long long k : ks) {
    const auto pose = dk.l1Pose(c, 0, k);
    if (!pose) continue;
    const long long tr = k % 729;
    const od::Vec3 delta{static_cast<double>(tr % 9 - 4), static_cast<double>(tr / 9 % 9 - 4), static_cast<double>(tr / 81 - 4)};
    if (k < c.pairPoses) {
      const auto& pm = c.pairs[static_cast<std::size_t>(k / 729 / 48)];
      const auto w1 = world(*pose, pm.m1), w2 = world(*pose, pm.m2);
      const auto& s1 = sepas[pm.s1].xyz;
      const auto& s2 = sepas[pm.s2].xyz;
      od::Vec3 us{s2[0] - s1[0], s2[1] - s1[1], s2[2] - s1[2]}, um{w2[0] - w1[0], w2[1] - w1[1], w2[2] - w1[2]};
      const double ls = std::sqrt(us[0] * us[0] + us[1] * us[1] + us[2] * us[2]), lm = std::sqrt(um[0] * um[0] + um[1] * um[1] + um[2] * um[2]);
      CHECK((us[0] * um[0] + us[1] * um[1] + us[2] * um[2]) / (ls * lm) == doctest::Approx(1.0).epsilon(1e-9));
      CHECK(std::abs(ls - lm) <= f.p.sepaDistTolerance + 1e-6);
      for (int d = 0; d < 3; ++d) CHECK(0.5 * (w1[d] + w2[d]) == doctest::Approx(0.5 * (s1[d] + s2[d]) + delta[d]).epsilon(1e-9));
      ++checkedPairs;
    } else {
      const auto& sp = c.singles[static_cast<std::size_t>((k - c.pairPoses) / 729 / 72)];
      const auto w = world(*pose, sp[1]);
      for (int d = 0; d < 3; ++d) CHECK(w[d] == doctest::Approx(sepas[sp[0]].xyz[d] + delta[d]).epsilon(1e-9));
      ++checkedSingles;
    }
  }
  CHECK(checkedPairs > 0);
  CHECK(checkedSingles > 0);
  // Pair alignments only.
  od::Dock9bParams pp = f.p;
  pp.sepaSingle = false;
  const od::MobmDocker dkp(*f.x.mobSys, f.field, pp, &sepas);
  const auto cp = dkp.conformation(f.x.mobSys->inputAngles());
  CHECK(cp.singles.empty());
  CHECK(cp.pairs.size() == c.pairs.size());
  CHECK(dkp.l1Count(cp) == c.pairPoses);
  for (long long k = 0; k < c.pairPoses; k += c.pairPoses / 50 + 1) {
    const auto a = dkp.l1Pose(cp, 0, k), b = dk.l1Pose(c, 0, k);
    REQUIRE(a.has_value() == b.has_value());
    if (a) CHECK(a->t == b->t);
  }
  // Triple alignments only.
  od::Dock9bParams p3t = f.p;
  p3t.sepaPair = false;
  p3t.sepaSingle = false;
  p3t.sepaTriple = true;
  p3t.sepaPosRange = 1.0;
  const od::MobmDocker dkt(*f.x.mobSys, f.field, p3t, &sepas);
  const auto ct = dkt.conformation(f.x.mobSys->inputAngles());
  MESSAGE(std::format("posm3 triple: triples {}, L1poss {}", ct.triples.size(), dkt.l1Count(ct)));
  REQUIRE(!ct.triples.empty());
  CHECK(ct.pairs.empty());
  CHECK(ct.singles.empty());
  CHECK(dkt.l1Count(ct) == static_cast<long long>(ct.triples.size()) * 27);
  int tripleChecked = 0;
  double maxRmsd = 0;
  for (std::size_t r = 0; r < ct.triples.size(); r += ct.triples.size() / 300 + 1) {
    const auto& tm = ct.triples[r];
    CHECK(tm.s[0] < tm.s[1]);
    CHECK(tm.s[1] < tm.s[2]);
    const long long k = static_cast<long long>(r) * 27 + 13;  // zero translation
    const auto pose = dkt.l1Pose(ct, 0, k);
    if (!pose) continue;
    od::Vec3 cs{}, cm{};
    double sq = 0;
    for (int v = 0; v < 3; ++v) {
      const auto w = world(*pose, tm.m[v]);
      const auto& s = sepas[tm.s[v]].xyz;
      for (int d = 0; d < 3; ++d) {
        cs[d] += s[d] / 3;
        cm[d] += w[d] / 3;
        sq += (w[d] - s[d]) * (w[d] - s[d]);
      }
      CHECK(sepas[tm.s[v]].sign * (f.x.mobSys->typing(0).charge[probeOf[tm.m[v]] - f.x.mobSys->coordOffset(0)] > 0 ? 1 : -1) == -1);
    }
    for (int d = 0; d < 3; ++d) CHECK(cm[d] == doctest::Approx(cs[d]).epsilon(1e-9));
    maxRmsd = std::max(maxRmsd, std::sqrt(sq / 3));
    ++tripleChecked;
  }
  MESSAGE(std::format("posm3 triple: checked {}, max RMSD of the three anchor atoms {:.3f} A", tripleChecked, maxRmsd));
  CHECK(tripleChecked > 0);
  CHECK(maxRmsd <= f.p.sepaDistTolerance);
  // Defaults: triples + pairs, step 2 Å, range 2 Å (27 translations).
  const od::MobmDocker dkr(*f.x.mobSys, f.field, pDef, &sepas);
  const auto cr = dkr.conformation(f.x.mobSys->inputAngles());
  CHECK(cr.pairs.size() == c.pairs.size());
  CHECK(cr.triples.size() == ct.triples.size());
  CHECK(cr.singles.empty());
  CHECK(cr.pairPoses == static_cast<long long>(cr.pairs.size()) * 48 * 27);
  CHECK(dkr.l1Count(cr) == cr.pairPoses + static_cast<long long>(cr.triples.size()) * 27);
  for (long long r = 0; r < static_cast<long long>(cr.triples.size()); r += 997) {  // zero translations
    const auto a = dkr.l1Pose(cr, 0, cr.pairPoses + r * 27 + 13), b = dkt.l1Pose(ct, 0, r * 27 + 13);
    REQUIRE(a.has_value() == b.has_value());
    if (a) {
      CHECK(a->q.w == doctest::Approx(b->q.w));
      CHECK(a->t == b->t);
    }
  }
  // posm4 = posm3 followed by posm1.
  od::Dock9bParams p4 = f.p;
  p4.posm = od::Posm::Posm4;
  const od::MobmDocker dk4(*f.x.mobSys, f.field, p4, &sepas);
  const auto c4 = dk4.conformation(f.x.mobSys->inputAngles());
  CHECK(dk4.l1Count(c4) == dk.l1Count(c) + 400);
  DockFixture f1;
  const od::MobmDocker dk1(*f1.x.mobSys, f1.field, f1.p);
  const auto c1 = dk1.conformation(f1.x.mobSys->inputAngles());
  for (long long k = 0; k < 20; ++k) {
    const auto a = dk4.l1Pose(c4, 9, dk.l1Count(c) + k), b = dk1.l1Pose(c1, 9, k);
    REQUIRE(a.has_value() == b.has_value());
    if (a) CHECK(a->t == b->t);
  }
  // L1poss per posm3 rotation (spec 15.1.21): the L2poss and counters of l1Pose + quickFilter.
  {
    od::Dock9bParams pr4 = pDef;
    pr4.posm = od::Posm::Posm4;
    pr4.boxSize = {34, 34, 34};
    pr4.sepaSingle = true;  // all three kinds of rotation
    pr4.sepaPotPercentile = 1.0;
    pr4.sepaChargeThreshold = 0.5;
    pr4.maxOverlapAtomNum = 6;
    const auto sepasB = od::stamSepas(f.field, pr4);
    const od::MobmDocker dkb(*f.x.mobSys, f.field, pr4, &sepasB);
    const auto cb = dkb.conformation(f.x.mobSys->inputAngles());
    const long long n = dkb.l1Count(cb), nt = dkb.posm3Translations();
    CHECK(nt == 27);
    od::Dock9bStats ns;
    std::vector<long long> naive;
    long long inBoxSame = 0, sepaFit = 0;
    for (long long k = 0; k < n; ++k) {
      const auto pose = dkb.l1Pose(cb, 9, k);
      if (!pose) {
        ++ns.outsideBox;
        continue;
      }
      if (k < cb.sepaPoses()) {
        const od::Pose q = dkb.posm3Pose(dkb.posm3Rot(cb, k / nt), static_cast<int>(k % nt));
        CHECK(q.t == pose->t);
        od::Dock9bStats s1, s2;
        const bool r1 = dkb.quickFilter(cb, *pose, s1), r2 = dkb.quickFilterInBox(cb, *pose, s2);
        ++sepaFit;
        inBoxSame += r1 == r2 && s1.overlap == s2.overlap && s1.far == s2.far && s1.outsideBox == 0 && s2.outsideBox == 0;
      }
      if (dkb.quickFilter(cb, *pose, ns)) naive.push_back(k);
    }
    MESSAGE(std::format("posm4 34 A box: pairs {}, singles {}, triples {}; L1poss {}, L2poss {}, outside {}, overlap {}, far {}",
                        cb.pairs.size(), cb.singles.size(), cb.triples.size(), n, naive.size(), ns.outsideBox,
                        ns.overlap, ns.far));
    CHECK(!cb.pairs.empty());
    CHECK(!cb.singles.empty());
    CHECK(!cb.triples.empty());
    CHECK(ns.outsideBox > 0);
    CHECK(!naive.empty());
    CHECK(sepaFit > 0);
    CHECK(inBoxSame == sepaFit);
    std::vector<od::MobmDocker::PairJob> jobs(1);
    const auto ang = f.x.mobSys->inputAngles();
    jobs[0].mobAnglesDeg = ang;
    jobs[0].stream = 9;
    dkb.prepareJobs(jobs);
    REQUIRE(jobs[0].prepared);
    CHECK(jobs[0].prepared->l2 == naive);
    const od::Dock9bStats& ps = jobs[0].prepared->stats;
    CHECK(ps.outsideBox == ns.outsideBox);
    CHECK(ps.overlap == ns.overlap);
    CHECK(ps.far == ns.far);
  }
  // posm4, pair docked before: only the posm1 part.
  std::vector<od::DockSample> s4;
  const auto st4 = dk4.dockPair(f.x.mobSys->inputAngles(), 9, 1 << 20, [&](std::vector<od::DockSample>& ch) { s4.insert(s4.end(), ch.begin(), ch.end()); }, true);
  std::vector<od::DockSample> s1;
  const auto st1 = dk1.dockPair(f1.x.mobSys->inputAngles(), 9, s1);
  CHECK(st4.l1 == 400);
  CHECK(st4.sepaPairs == 0);
  CHECK(st4.l2 == st1.l2);
  REQUIRE(s4.size() == s1.size());
  for (std::size_t i = 0; i < s4.size(); ++i) CHECK(s4[i].energy == s1[i].energy);
}

TEST_CASE("dock9b: posm2_fft against several StaM conformations equals each alone (spec 15.1.22)") {
  FieldSet fs;
  od::Dock9bParams p = fs.p;
  p.posm = od::Posm::Posm2Fft;
  const auto angles = fs.x.mobSys->inputAngles();
  for (int top : {1, 3}) {
    INFO(top);
    p.fftTopPerRot = top;
    std::vector<od::FftStam> ss;
    for (const auto& [name, field] : fs.fields) ss.push_back(od::buildFftStam(field, 10.0));
    std::vector<std::unique_ptr<od::MobmDocker>> own;
    std::vector<const od::MobmDocker*> dks;
    for (std::size_t s = 0; s < ss.size(); ++s) {
      own.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, fs.fields[s].second, p, nullptr, &ss[s]));
      dks.push_back(own.back().get());
    }
    // Reversed and repeated order as well.
    for (const auto& order : {std::vector<std::size_t>{0, 1, 2}, std::vector<std::size_t>{2, 0, 1, 2}}) {
      std::vector<const od::MobmDocker*> sel;
      for (std::size_t s : order) sel.push_back(dks[s]);
      const auto cs = od::MobmDocker::conformations(sel, angles);
      REQUIRE(cs.size() == order.size());
      for (std::size_t k = 0; k < order.size(); ++k) {
        INFO(fs.fields[order[k]].first);
        const auto one = dks[order[k]]->conformation(angles);
        REQUIRE(!one.fftPoses.empty());
        REQUIRE(cs[k].fftPoses.size() == one.fftPoses.size());
        CHECK(cs[k].fftEnergy == one.fftEnergy);
        CHECK(cs[k].coords == one.coords);
        for (std::size_t i = 0; i < one.fftPoses.size(); ++i) {
          CHECK(cs[k].fftPoses[i].t == one.fftPoses[i].t);
          CHECK(cs[k].fftPoses[i].q.w == one.fftPoses[i].q.w);
        }
      }
    }
    // The StaM conformations give different translations.
    CHECK(dks[0]->conformation(angles).fftEnergy != dks[2]->conformation(angles).fftEnergy);
  }
  // Dockers of different boxes are rejected.
  od::Dock9bParams p2 = p;
  p2.boxSize = {40, 40, 40};
  const od::FftStam s0 = od::buildFftStam(fs.fields[0].second, 10.0);
  const od::MobmDocker a(*fs.x.mobSys, fs.fields[0].second, p, nullptr, &s0), b(*fs.x.mobSys, fs.fields[0].second, p2, nullptr, &s0);
  const std::vector<const od::MobmDocker*> mixed{&a, &b};
  CHECK_THROWS_AS(od::MobmDocker::conformations(mixed, angles), od::Error);
  // Prepared from a given conformation: as from the angles.
  std::vector<od::MobmDocker::PairJob> j1(1), j2(1);
  j1[0].mobAnglesDeg = j2[0].mobAnglesDeg = angles;
  j2[0].conf = a.conformation(angles);
  a.prepareJobs(j1);
  a.prepareJobs(j2);
  CHECK(!j2[0].conf);
  CHECK(j1[0].prepared->l2 == j2[0].prepared->l2);
  CHECK(j1[0].stats.overlap == j2[0].stats.overlap);
}

TEST_CASE("dock9b: posm2_fft (2KJJ, 1L2Y)") {
  CHECK(od::fftSize(1) == 1);
  CHECK(od::fftSize(11) == 12);
  CHECK(od::fftSize(97) == 98);
  CHECK(od::fftSize(129) == 135);
  DockFixture f(od::Posm::Posm2Fft);
  const od::FftStam s = od::buildFftStam(f.field, 10.0);
  const od::GridGeometry& g = s.geom;
  for (int d = 0; d < 3; ++d) CHECK(s.pad[d] == od::fftSize(g.n[d]));
  const int nt = f.field.typeMapCount();
  // Clamped maps (spec 8.4.3) and the direct trilinear sum.
  std::vector<std::vector<float>> cm(static_cast<std::size_t>(nt) + 2);
  for (int m = 0; m < nt + 2; ++m) {
    cm[m] = f.field.pointMap(m);
    CHECK(cm[m] == (m < nt ? f.grids.typeMaps[m] : m == nt ? f.grids.elecMap : f.grids.desolvMap));
    for (float& v : cm[m]) v = m == nt ? std::clamp(v, -10.0f, 10.0f) : std::min(v, 10.0f);
  }
  const od::MobmDocker dk(*f.x.mobSys, f.field, f.p, nullptr, &s);
  const auto c = dk.conformation(f.x.mobSys->inputAngles());
  const auto& mt = f.x.mobSys->typing(0);
  std::vector<int> probeOf, mapOf;
  std::vector<double> qOf;
  for (std::size_t a = 0; a < mt.type.size(); ++a)
    if (mt.type[a] != od::Ad4Type::Merged) {
      probeOf.push_back(f.x.mobSys->coordOffset(0) + static_cast<int>(a));
      mapOf.push_back(f.field.mapOf(mt.type[a]));
      qOf.push_back(mt.charge[a]);
    }
  auto direct = [&](const od::Pose& pose, bool* inBox = nullptr) {
    double e = 0;
    const od::QuatF qf = od::toFloat(pose.q);
    for (std::size_t a = 0; a < probeOf.size(); ++a) {
      const od::Vec3f& y = c.coords[probeOf[a]];
      const od::Vec3f r = od::rotate(qf, {y[0] - static_cast<float>(c.center[0]), y[1] - static_cast<float>(c.center[1]),
                                          y[2] - static_cast<float>(c.center[2])});
      const od::Vec3 w{r[0] + pose.t[0], r[1] + pose.t[1], r[2] + pose.t[2]};
      for (int d = 0; d < 3; ++d)
        if (inBox && (w[d] < f.p.boxCenter[d] - 0.5 * f.p.boxSize[d] + 0.01 || w[d] > f.p.boxCenter[d] + 0.5 * f.p.boxSize[d] - 0.01))
          *inBox = false;
      if (!g.inside(w)) return 1e30;
      e += od::interpolate(g, cm[mapOf[a]], w) + qOf[a] * od::interpolate(g, cm[nt], w) +
           std::abs(qOf[a]) * od::interpolate(g, cm[nt + 1], w);
    }
    return e;
  };
  REQUIRE(!c.fftPoses.empty());
  CHECK(c.fftPoses.size() <= 72);
  CHECK(dk.l1Count(c) == static_cast<long long>(c.fftPoses.size()));
  double maxErr = 0;
  for (std::size_t k = 0; k < c.fftPoses.size(); ++k) {
    const od::Pose& pose = c.fftPoses[k];
    CHECK(dk.l1Pose(c, 0, static_cast<long long>(k))->t == pose.t);
    bool in = true;
    const double e = direct(pose, &in);
    CHECK(in);
    maxErr = std::max(maxErr, std::abs(e - c.fftEnergy[k]));
    CHECK(c.fftEnergy[k] == doctest::Approx(e).epsilon(1e-4).scale(50.0));
    // On a grid point.
    for (int d = 0; d < 3; ++d) {
      const double u = (pose.t[d] - g.origin[d]) / g.spacing;
      CHECK(std::abs(u - std::round(u)) < 1e-9);
    }
    // Lowest among the in-box neighbouring translations.
    if (k % 6 == 0)
      for (int dz = -2; dz <= 2; ++dz)
        for (int dy = -2; dy <= 2; ++dy)
          for (int dx = -2; dx <= 2; ++dx) {
            od::Pose p2 = pose;
            p2.t = {pose.t[0] + dx * g.spacing, pose.t[1] + dy * g.spacing, pose.t[2] + dz * g.spacing};
            bool in2 = true;
            const double e2 = direct(p2, &in2);
            if (in2) CHECK(e2 >= e - 1e-2);
          }
  }
  MESSAGE(std::format("posm2_fft: pad {}x{}x{}, kept {}, max |FFT - direct| {:.2e}, lowest {:.3f}", s.pad[0], s.pad[1],
                      s.pad[2], c.fftPoses.size(), maxErr, *std::min_element(c.fftEnergy.begin(), c.fftEnergy.end())));
  // fftStride 2 (spec 15.1.23): every second grid point; energies of the trilinear sum on those points.
  {
    const od::FftStam s2 = od::buildFftStam(f.field, 10.0, 2);
    const od::GridGeometry& g2 = s2.geom;
    CHECK(s2.stride == 2);
    CHECK(g2.spacing == 2 * g.spacing);
    CHECK(g2.origin == g.origin);
    for (int d = 0; d < 3; ++d) {
      CHECK(g2.n[d] == (g.n[d] - 1) / 2 + 1);
      CHECK(s2.pad[d] == od::fftSize(g2.n[d]));
    }
    std::vector<std::vector<float>> cm2(cm.size(), std::vector<float>(g2.size()));
    for (std::size_t m = 0; m < cm.size(); ++m)
      for (int k = 0; k < g2.n[2]; ++k)
        for (int j = 0; j < g2.n[1]; ++j)
          for (int i = 0; i < g2.n[0]; ++i) cm2[m][g2.index(i, j, k)] = cm[m][g.index(2 * i, 2 * j, 2 * k)];
    const od::MobmDocker dk2(*f.x.mobSys, f.field, f.p, nullptr, &s2);
    const auto c2 = dk2.conformation(f.x.mobSys->inputAngles());
    REQUIRE(!c2.fftPoses.empty());
    CHECK(c2.fftPoses.size() <= c.fftPoses.size());
    double err2 = 0;
    for (std::size_t k = 0; k < c2.fftPoses.size(); ++k) {
      const od::Pose& pose = c2.fftPoses[k];
      const od::QuatF qf = od::toFloat(pose.q);
      double e = 0;
      for (std::size_t a = 0; a < probeOf.size(); ++a) {
        const od::Vec3f& y = c2.coords[probeOf[a]];
        const od::Vec3f r = od::rotate(qf, {y[0] - static_cast<float>(c2.center[0]), y[1] - static_cast<float>(c2.center[1]),
                                            y[2] - static_cast<float>(c2.center[2])});
        const od::Vec3 w{r[0] + pose.t[0], r[1] + pose.t[1], r[2] + pose.t[2]};
        REQUIRE(g2.inside(w));
        e += od::interpolate(g2, cm2[mapOf[a]], w) + qOf[a] * od::interpolate(g2, cm2[nt], w) +
             std::abs(qOf[a]) * od::interpolate(g2, cm2[nt + 1], w);
      }
      err2 = std::max(err2, std::abs(e - c2.fftEnergy[k]));
      CHECK(c2.fftEnergy[k] == doctest::Approx(e).epsilon(1e-4).scale(50.0));
      // On a point of the FFT grid, i.e. on every second point of the field grid.
      for (int d = 0; d < 3; ++d) {
        const double u = (pose.t[d] - g.origin[d]) / g2.spacing;
        CHECK(std::abs(u - std::round(u)) < 1e-9);
      }
    }
    MESSAGE(std::format("posm2_fft fftStride 2: pad {}x{}x{}, kept {}, max |FFT - direct| {:.2e}", s2.pad[0], s2.pad[1], s2.pad[2],
                        c2.fftPoses.size(), err2));
    CHECK_THROWS_AS(od::buildFftStam(f.field, 10.0, 0), od::Error);
  }
  // fftTopPerRot 3: ascending per orientation, the first one as with 1.
  od::Dock9bParams p3 = f.p;
  p3.fftTopPerRot = 3;
  const od::MobmDocker dk3(*f.x.mobSys, f.field, p3, nullptr, &s);
  const auto c3 = dk3.conformation(f.x.mobSys->inputAngles());
  REQUIRE(c3.fftPoses.size() == 3 * c.fftPoses.size());
  for (std::size_t k = 0; k < c.fftPoses.size(); ++k) {
    CHECK(c3.fftPoses[3 * k].t == c.fftPoses[k].t);
    CHECK(c3.fftEnergy[3 * k] <= c3.fftEnergy[3 * k + 1]);
    CHECK(c3.fftEnergy[3 * k + 1] <= c3.fftEnergy[3 * k + 2]);
  }
  // 9b of the pair.
  std::vector<od::DockSample> out;
  const auto st = dk.dockPair(f.x.mobSys->inputAngles(), 0, out);
  CHECK(st.l1 == dk.l1Count(c));
  CHECK(st.l1 == st.overlap + st.far + st.l2);
  CHECK(static_cast<long long>(out.size()) == st.l2 - st.leftGrid);
  // flexres: pointMap includes the moved atoms.
  const double cap = std::max(f.p.maxNearDistance, f.p.clashDistance) + 0.75;
  std::vector<char> mask(f.x.stam.size());
  for (std::size_t a = 0; a < mask.size(); ++a) mask[a] = a % 3 == 0;
  const od::GridAtoms fixedAtoms = od::subset(f.x.stam, mask, false);
  const auto flex = od::StamField::flexres(std::make_shared<const od::EnergyGrids>(od::buildEnergyGrids(g, fixedAtoms, f.x.mobm.type, {})),
                                           std::make_shared<const od::DistanceGrid>(od::buildDistanceGrid(g, fixedAtoms, cap)),
                                           od::subset(f.x.stam, mask, true), {});
  std::vector<std::vector<float>> fm;
  for (int m = 0; m < nt + 2; ++m) fm.push_back(flex.pointMap(m));
  std::mt19937 rng(5);
  for (int r = 0; r < 300; ++r) {
    const int i = 1 + static_cast<int>(rng() % (g.n[0] - 2)), j = 1 + static_cast<int>(rng() % (g.n[1] - 2)),
              k = 1 + static_cast<int>(rng() % (g.n[2] - 2));
    const od::Vec3 p = g.point(i, j, k);
    const int m = r % nt;
    const std::size_t idx = g.index(i, j, k);
    for (double q : {-0.4, 0.0, 0.7}) {
      const double want = flex.energy<double>(m, q, p, nullptr);
      const double got = fm[m][idx] + q * fm[nt][idx] + std::abs(q) * fm[nt + 1][idx];
      CHECK(got == doctest::Approx(want).epsilon(1e-4).scale(1.0));
    }
  }
}

TEST_CASE("dock9b: poses, quick filter and gradient (2KJJ, 1L2Y)") {
  DockFixture f;
  const od::MobmDocker dk(*f.x.mobSys, f.field, f.p);
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
  const od::MobmDocker dk2(*f2.x.mobSys, f2.field, f2.p);
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
  f.p.minimizer.maxSteps = 300;
  const od::MobmDocker dk(*f.x.mobSys, f.field, f.p);
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
  f.p.minimizer.maxSteps = 100;
  const od::MobmDocker dk(*f.x.mobSys, f.field, f.p);
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
#ifdef OD_HAVE_CUDA
  // GPU shaking (spec 15.1.17): the CPU values within the summation differences.
  if (od::gpuDeviceCount() > 0) {
    for (double step : {0.4, 0.2}) {
      od::ShakeParams g = sp;
      g.step = step;
      od::GpuShaker sh(ep, g, 0);
      std::vector<od::Vec3f> near = mx;
      for (auto& v : near) v[0] += 6;
      const od::ShakeAtoms all = od::joinShakeAtoms(stam, od::shakeAtoms(*x.mobSys, x.mobParts, near));
      std::vector<char> shake(all.size(), 1);
      for (std::size_t a = 0; a < all.size(); a += 3) shake[a] = 0;
      for (const auto& s : {std::vector<char>{}, shake}) {
        const auto cpu = od::shakeLnF(all, ep, g, s);
        const auto gpu = sh.lnF(all, s);
        REQUIRE(gpu.size() == cpu.size());
        double worst = 0, sumC = 0, sumG = 0;
        for (std::size_t a = 0; a < cpu.size(); ++a) {
          worst = std::max(worst, std::abs(cpu[a] - gpu[a]));
          sumC += cpu[a];
          sumG += gpu[a];
          if (cpu[a] == 0) CHECK(gpu[a] == 0);
        }
        INFO("step ", step, " shake mask ", !s.empty());
        MESSAGE(std::format("GPU shaking step {} mask {}: largest ln f difference {:.2e}, sum {:.6f} vs {:.6f}", step, !s.empty(), worst, sumG, sumC));
        CHECK(worst < 1e-4);
        CHECK(sumG == doctest::Approx(sumC).epsilon(1e-6));
        const od::ComplexLnF lc = od::complexSumLnF(stam, alone, od::shakeAtoms(*x.mobSys, x.mobParts, near), ep, g);
        const od::ComplexLnF lg = od::complexSumLnF(stam, alone, od::shakeAtoms(*x.mobSys, x.mobParts, near), ep, g,
                                                    [&](const od::ShakeAtoms& a, std::span<const char> m) { return sh.lnF(a, m); });
        CHECK(lg.stam == doctest::Approx(lc.stam).epsilon(1e-6));
        CHECK(lg.mobm == doctest::Approx(lc.mobm).epsilon(1e-6));
      }
    }
  }
#endif
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


TEST_CASE("checkpoint: SHA-256, binary state and directories") {
  CHECK(od::sha256Hex(std::string()) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  CHECK(od::sha256Hex(std::string("abc")) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  CHECK(od::sha256Hex(std::string("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")) ==
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
  CHECK(od::sha256Hex(std::string(1000000, 'a')) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
  const auto dir = std::filesystem::temp_directory_path() / std::format("od_ckpt_{}", ::getpid());
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  {
    std::ofstream f(dir / "in.txt", std::ios::binary);
    f << "abc";
  }
  CHECK(od::sha256File((dir / "in.txt").string()) == od::sha256Hex(std::string("abc")));

  od::BinWriter w;
  od::putStateHeader(w);
  w.put<int>(-7);
  w.putVec(std::vector<float>{1.5f, -2.25f});
  w.putVec(std::vector<int>{});
  w.put<double>(0.1);
  od::BinReader r(w.data());
  od::checkStateHeader(r);
  CHECK(r.get<int>() == -7);
  CHECK(r.getVec<float>() == std::vector<float>{1.5f, -2.25f});
  CHECK(r.getVec<int>().empty());
  CHECK(r.get<double>() == 0.1);
  CHECK(r.atEnd());
  CHECK_THROWS_AS(r.get<char>(), od::Error);
  od::BinReader cut(w.data().substr(0, w.data().size() - 3));
  od::checkStateHeader(cut);
  cut.get<int>();
  cut.getVec<float>();
  cut.getVec<int>();
  CHECK_THROWS_AS(cut.get<double>(), od::Error);
  od::BinReader junk(std::string(16, 'x'));
  CHECK_THROWS_AS(od::checkStateHeader(junk), od::Error);

  const od::Manifest m{"0.40.0", 3, {{"stamFile", "/a b/s.cif", "11"}, {"mobmFiles", "/m'1.sdf", "22"}}, "ff"};
  const std::string out = dir.string();
  CHECK(od::latestCheckpoint(out) == 0);
  for (int round : {1, 2, 3}) {
    od::Manifest mr = m;
    mr.round = round;
    const auto d = od::commitCheckpoint(out, round, {{"state.bin", "S"}, {"manifest.toml", od::manifestText(mr)}},
                                        [](const std::string& tmp) { std::ofstream(tmp + "/extra.txt") << "e"; });
    CHECK(d == od::checkpointDir(out, round));
    CHECK(std::filesystem::exists(d + "/extra.txt"));
  }
  CHECK(std::filesystem::read_symlink(dir / "checkpoint" / "last") == "round_3");
  CHECK(od::latestCheckpoint(out) == 3);
  // An incomplete directory (no manifest) and a temporary one are ignored.
  std::filesystem::create_directories(dir / "checkpoint" / "round_9");
  std::filesystem::create_directories(dir / "checkpoint" / "round_10.tmp");
  CHECK(od::latestCheckpoint(out) == 3);
  const od::Manifest back = od::readManifest(od::checkpointDir(out, 3) + "/manifest.toml");
  CHECK(back.version == "0.40.0");
  CHECK(back.round == 3);
  CHECK(back.stateSha256 == "ff");
  REQUIRE(back.inputs.size() == 2);
  CHECK(back.inputs[0].path == "/a b/s.cif");
  CHECK(back.inputs[1].role == "mobmFiles");
  CHECK(back.inputs[1].path == "/m'1.sdf");
  CHECK(back.inputs[1].sha256 == "22");
  CHECK(od::readFile(od::checkpointDir(out, 2) + "/state.bin") == "S");
  // Rewriting a round replaces it.
  od::commitCheckpoint(out, 2, {{"state.bin", "T"}, {"manifest.toml", od::manifestText(m)}}, {});
  CHECK(od::readFile(od::checkpointDir(out, 2) + "/state.bin") == "T");
  CHECK_FALSE(std::filesystem::exists(od::checkpointDir(out, 2) + "/extra.txt"));
  CHECK(std::filesystem::read_symlink(dir / "checkpoint" / "last") == "round_2");
  std::filesystem::remove_all(dir);
}

namespace {

// Energy systems of the flat-form tests: 1L2Y (MobM), 2KJJ with a free loop
// (closure virtual atoms), 1ZNF with metalBond = "break" (distance restraints)
// and 1D3Z (287 torsions).
struct FlatCase {
  std::string name;
  std::unique_ptr<od::EnergySystem> sys;
};

struct FlatFixture {
  ConfigFixture cf;
  StaFixture sf;
  std::vector<od::Molecule> big;
  std::vector<od::Topology> bigTopo;
  std::optional<od::MolConfig> cfgA, cfgZ, cfgBig;
  od::MolConfig cfgS;
  std::vector<FlatCase> cases;

  FlatFixture() {
    std::vector<std::string> errors;
    std::ofstream(cf.path) << od::defaultMolConfigText(cf.inputs);
    cfgA = od::readMolConfig(cf.path, cf.inputs, {}, errors);
    REQUIRE(cfgA);
    cases.push_back({"1L2Y", std::make_unique<od::EnergySystem>(
                                 std::vector<od::EnergyPart>{{&cf.mols[0], &cfgA->molecules[0].tree}}, od::EnergyParams{})});

    cfgS = sf.ok(freeResidues(sf.text, "A", 10, 14));
    cases.push_back({"2KJJ loop", std::make_unique<od::EnergySystem>(
                                      std::vector<od::EnergyPart>{{&sf.mols[0], &cfgS.molecules[0].tree}}, od::EnergyParams{})});

    cf.mols.push_back(zincFinger());
    cf.topos.push_back(od::analyzeTopology(cf.mols.back(), {}));
    od::assignDssp({&cf.mols.back()}, {&cf.topos.back()});
    cf.inputs.clear();
    for (std::size_t i = 0; i < cf.mols.size(); ++i) cf.inputs.push_back({&cf.mols[i], &cf.topos[i]});
    std::ofstream(cf.path) << phi10(setQuoted(uncommentMetal(od::defaultMolConfigText(cf.inputs)),
                                              "id           = \"mobm:1ZNF.cif:A\"", "metalBond    = ", "break"));
    cfgZ = od::readMolConfig(cf.path, cf.inputs, {}, errors);
    REQUIRE(cfgZ);
    cases.push_back({"1ZNF metal break", std::make_unique<od::EnergySystem>(
                                             std::vector<od::EnergyPart>{{&cf.mols.back(), &cfgZ->molecules.back().tree}},
                                             od::EnergyParams{})});

    big.push_back(od::loadMolecules(kData + "/1D3Z.cif", od::Role::MobM, ccd()).front());
    bigTopo.push_back(od::analyzeTopology(big[0], {}));
    od::assignDssp({&big[0]}, {&bigTopo[0]});
    const std::vector<od::MolInput> in{{&big[0], &bigTopo[0]}};
    std::ofstream(cf.path) << od::defaultMolConfigText(in);
    cfgBig = od::readMolConfig(cf.path, in, {}, errors);
    REQUIRE(cfgBig);
    cases.push_back({"1D3Z", std::make_unique<od::EnergySystem>(
                                 std::vector<od::EnergyPart>{{&big[0], &cfgBig->molecules[0].tree}}, od::EnergyParams{})});
  }
};

// Angle sets: input angles, then perturbations of growing width.
template <class T>
std::vector<std::vector<T>> flatAngleSets(const od::EnergySystem& sys, int n) {
  std::vector<std::vector<T>> out;
  for (int s = 0; s < n; ++s) {
    const auto a = s == 0 ? sys.inputAngles() : perturbed(sys.inputAngles(), 5.0 * s, 100 + s);
    out.emplace_back(a.begin(), a.end());
  }
  return out;
}

}  // namespace

TEST_CASE("energy kernel: flat form equals EnergySystem::evaluate") {
  FlatFixture f;
  for (const auto& c : f.cases) {
    const od::EnergyModel model = c.sys->model();
    const od::EnergyModelView v = model.view();
    CHECK(v.torsionCount == c.sys->torsionCount());
    CHECK(v.pairCount == static_cast<int>(c.sys->pairCount()));
    std::vector<double> xd(3 * v.coordCount), gd(3 * v.coordCount), gradD(v.torsionCount), grad;
    std::vector<float> xf(3 * v.coordCount);
    int mismatches = 0;
    for (const auto& a : flatAngleSets<double>(*c.sys, 8)) {
      std::vector<od::Vec3> coords;
      const auto ref = c.sys->evaluate(a, &grad, &coords);
      const auto e = od::flat::evaluate(v, a.data(), xd.data(), gd.data(), gradD.data());
      mismatches += termArray(ref) != termArray(e) || ref.clashPairs != e.clashPairs || grad != gradD;
      for (int i = 0; i < v.coordCount; ++i)
        for (int k = 0; k < 3; ++k) mismatches += coords[i][k] != xd[3 * i + k];

      const std::vector<float> af(a.begin(), a.end());
      std::vector<od::Vec3f> coordsF;
      const auto refF = c.sys->evaluate(std::span<const float>(af), &grad, &coordsF, nullptr);
      const auto eF = od::flat::evaluate(v, af.data(), xf.data(), gd.data(), gradD.data());
      mismatches += termArray(refF) != termArray(eF) || refF.clashPairs != eF.clashPairs || grad != gradD;
      for (int i = 0; i < v.coordCount; ++i)
        for (int k = 0; k < 3; ++k) mismatches += coordsF[i][k] != xf[3 * i + k];
    }
    INFO(c.name);
    CHECK(mismatches == 0);
  }
}

#ifdef OD_HAVE_CUDA
namespace {

// Largest differences of one evaluation from a reference: energy terms
// relative to the sum of their magnitudes, gradient relative to its largest
// component, coordinates in Å.
struct EvalDiff {
  double energy = 0, gradient = 0, coord = 0;
  void add(const od::EnergyTerms& ref, const std::vector<double>& gRef, const std::vector<od::Vec3>& xRef,
           const std::array<double, 7>& e, const double* g, auto xAt) {
    const auto r = termArray(ref);
    const double scale = std::max(1.0, std::abs(ref.vdw) + std::abs(ref.hbond) + std::abs(ref.elec) + std::abs(ref.desolv) +
                                           std::abs(ref.clash) + std::abs(ref.closure));
    for (int k = 0; k < 6; ++k) energy = std::max(energy, std::abs(r[k] - e[k]) / scale);
    double gMax = 1;
    for (double v : gRef) gMax = std::max(gMax, std::abs(v));
    for (std::size_t t = 0; t < gRef.size(); ++t) gradient = std::max(gradient, std::abs(gRef[t] - g[t]) / gMax);
    for (std::size_t i = 0; i < xRef.size(); ++i)
      for (int k = 0; k < 3; ++k) coord = std::max(coord, std::abs(xRef[i][k] - xAt(i, k)));
  }
};

}  // namespace

// Double: GPU against the CPU double reference. Float: GPU float and CPU
// float, both against the CPU double reference (spec 12.1).
TEST_CASE("gpu: energy and gradient match the CPU reference") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  const auto info = od::gpuDeviceInfo(0);
  MESSAGE(std::format("{} (sm_{}{}), {} MiB", info.name, info.major, info.minor, info.totalMem >> 20));
  FlatFixture f;
  constexpr int kSets = 64;
  for (const auto& c : f.cases) {
    const od::EnergyModel model = c.sys->model();
    od::GpuEnergy gpu(model);
    const int nt = c.sys->torsionCount(), nc = 3 * model.coordCount;
    const auto sets = flatAngleSets<double>(*c.sys, kSets);
    std::vector<od::EnergyTerms> ref(kSets);
    std::vector<std::vector<double>> gRef(kSets);
    std::vector<std::vector<od::Vec3>> xRef(kSets);
    for (int s = 0; s < kSets; ++s) ref[s] = c.sys->evaluate(sets[s], &gRef[s], &xRef[s]);

    // The sets repeated kRep times: more samples than resident blocks (blocks evaluate several).
    constexpr int kRep = 40, kAll = kSets * kRep;
    auto onGpu = [&]<class T>(T, EvalDiff& d, int& clashDiff) {
      std::vector<T> angles;
      for (int r = 0; r < kRep; ++r)
        for (const auto& a : sets) angles.insert(angles.end(), a.begin(), a.end());
      std::vector<od::FlatTerms> terms(kAll);
      std::vector<double> grad(static_cast<std::size_t>(kAll) * nt);
      std::vector<T> coords(static_cast<std::size_t>(kAll) * nc);
      gpu.evaluate(angles.data(), kAll, terms.data(), grad.data(), coords.data());
      for (int s = 0; s < kAll; ++s) {
        const T* x = coords.data() + static_cast<std::size_t>(s) * nc;
        d.add(ref[s % kSets], gRef[s % kSets], xRef[s % kSets], termArray(terms[s]),
              grad.data() + static_cast<std::size_t>(s) * nt, [&](std::size_t i, int k) { return static_cast<double>(x[3 * i + k]); });
        clashDiff += ref[s % kSets].clashPairs != terms[s].clashPairs;
      }
    };
    EvalDiff gd, gf, cf;
    int clashD = 0, clashF = 0, clashCpuF = 0;
    onGpu(0.0, gd, clashD);
    onGpu(0.0f, gf, clashF);
    for (int s = 0; s < kSets; ++s) {
      const std::vector<float> a(sets[s].begin(), sets[s].end());
      std::vector<double> g;
      std::vector<od::Vec3f> x;
      const auto e = c.sys->evaluate(std::span<const float>(a), &g, &x, nullptr);
      cf.add(ref[s], gRef[s], xRef[s], termArray(e), g.data(), [&](std::size_t i, int k) { return static_cast<double>(x[i][k]); });
      clashCpuF += ref[s].clashPairs != e.clashPairs;
    }
    MESSAGE(std::format("{}: {} torsions, {} pairs. GPU double: energy {:.1e}, gradient {:.1e}, coord {:.1e} Å. "
                        "Float vs double, GPU / CPU: energy {:.1e} / {:.1e}, gradient {:.1e} / {:.1e}, coord {:.1e} / {:.1e} Å, "
                        "clash count {} / {}",
                        c.name, nt, model.pairs.size(), gd.energy, gd.gradient, gd.coord, gf.energy, cf.energy, gf.gradient,
                        cf.gradient, gf.coord, cf.coord, clashF, clashCpuF));
    INFO(c.name);
    CHECK(clashD == 0);
    CHECK(gd.energy < 1e-10);
    CHECK(gd.gradient < 1e-10);
    CHECK(gd.coord < 1e-10);
    CHECK(gf.energy <= 2 * cf.energy + 1e-5);
    CHECK(gf.gradient <= 2 * cf.gradient + 1e-5);
    CHECK(gf.coord <= 2 * cf.coord + 1e-5);
  }
}
#endif

#ifdef OD_HAVE_CUDA
// 9a minimization with GPU energies against the CPU path (spec 15.1.3). Minimization paths of
// large clashing systems are sensitive to rounding: the control is the CPU path from start
// angles changed by 1e-10 degrees, and the GPU must agree with the CPU at least as often.
TEST_CASE("gpu: 9a minimization matches the CPU path") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FlatFixture f;
  for (const auto& c : f.cases) {
    const od::EnergySystem& sys = *c.sys;
    const int d = sys.torsionCount(), count = c.name == "1D3Z" ? 64 : 256;
    std::vector<double> start;
    for (int i = 0; i < count; ++i) {
      const auto a = perturbed(sys.inputAngles(), 30, 500 + i);
      start.insert(start.end(), a.begin(), a.end());
    }
    od::MinParams mp;
    std::vector<double> ref = start, ctl = start;
    for (double& v : ctl) v += 1e-10;
    std::vector<od::EnergyTerms> refTerms(count), ctlTerms(count);
    std::vector<int> refSteps(count), ctlSteps(count);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < count; ++i) {
      refTerms[i] = od::minimizeConf(sys, std::span<double>(ref.data() + static_cast<std::size_t>(i) * d, d), mp, &refSteps[i]);
      ctlTerms[i] = od::minimizeConf(sys, std::span<double>(ctl.data() + static_cast<std::size_t>(i) * d, d), mp, &ctlSteps[i]);
    }
    od::GpuConfEvaluator ev(sys.model(), 0, 8192);
    std::vector<double> a = start;
    std::vector<od::EnergyTerms> terms(count);
    std::vector<int> steps(count);
    od::minimizeConfBatch(ev, d, a, mp, terms.data(), steps.data());
    int sameSteps = 0, close = 0, ctlClose = 0;
    double worst = 0, maxAng = 0;
    auto relDiff = [](const od::EnergyTerms& x, const od::EnergyTerms& r) {
      return std::abs(x.filterEnergy() - r.filterEnergy()) / std::max(1.0, std::abs(r.filterEnergy()));
    };
    for (int i = 0; i < count; ++i) {
      sameSteps += steps[i] == refSteps[i];
      const double rel = relDiff(terms[i], refTerms[i]);
      worst = std::max(worst, rel);
      close += rel < 1e-6;
      ctlClose += relDiff(ctlTerms[i], refTerms[i]) < 1e-6;
      for (int k = 0; k < d; ++k) {
        double dd = std::abs(a[static_cast<std::size_t>(i) * d + k] - ref[static_cast<std::size_t>(i) * d + k]);
        maxAng = std::max(maxAng, std::min(dd, 360 - dd));
      }
    }
    MESSAGE(std::format("{}: {} samples; same step count {}, energy within 1e-6 (relative): GPU {}, CPU control {}; "
                        "worst {:.2e}, max angle diff {:.2e} deg",
                        c.name, count, sameSteps, close, ctlClose, worst, maxAng));
    INFO(c.name);
    CHECK(close >= ctlClose - count / 20);
    // gpu.conf9aPrecision = float (spec 15.1.12): final energies in double precision, against the
    // CPU (ref) and the CPU control (a perturbation of 1e-10 deg) as the scale of trajectory changes.
    od::GpuConfEvaluator evF(sys.model(), 0, 8192, true);
    std::vector<double> af = start;
    std::vector<od::EnergyTerms> termsF(count);
    std::vector<int> stepsF(count);
    od::minimizeConfBatch(evF, d, af, mp, termsF.data(), stepsF.data());
    // Samples below energy.Emax (1000) and, for those below it in both, dE = E - E_ref.
    auto stats = [&](const std::vector<od::EnergyTerms>& x, const std::vector<int>& st) {
      constexpr double kEmax = 1000;
      std::vector<double> dE, absE;
      int pass = 0, stepSum = 0;
      for (int i = 0; i < count; ++i) {
        stepSum += st[i];
        const double e = x[i].filterEnergy(), r = refTerms[i].filterEnergy();
        pass += e < kEmax;
        if (e < kEmax && r < kEmax) {
          dE.push_back(e - r);
          absE.push_back(std::abs(e - r));
        }
      }
      std::sort(absE.begin(), absE.end());
      auto pct = [&](double q) { return absE.empty() ? 0.0 : absE[static_cast<std::size_t>(q * (absE.size() - 1))]; };
      const int lower = static_cast<int>(std::count_if(dE.begin(), dE.end(), [](double v) { return v < -0.01; }));
      const int higher = static_cast<int>(std::count_if(dE.begin(), dE.end(), [](double v) { return v > 0.01; }));
      return std::format("< Emax {}; both {}: mean dE {:+.3f}, |dE| median {:.4f} p90 {:.3f} max {:.2f}, lower/higher by > 0.01 {}/{}; steps/sample {:.1f}",
                         pass, dE.size(), dE.empty() ? 0.0 : std::accumulate(dE.begin(), dE.end(), 0.0) / dE.size(),
                         pct(0.5), pct(0.9), pct(1.0), lower, higher, static_cast<double>(stepSum) / count);
    };
    const int refPass = static_cast<int>(std::count_if(refTerms.begin(), refTerms.end(), [](const od::EnergyTerms& e) { return e.filterEnergy() < 1000; }));
    MESSAGE(std::format("{}: CPU double: < Emax {}, steps/sample {:.1f}\n  GPU float:   {}\n  CPU control: {}", c.name,
                        refPass, std::accumulate(refSteps.begin(), refSteps.end(), 0.0) / count, stats(termsF, stepsF), stats(ctlTerms, ctlSteps)));
    for (int i = 0; i < count; ++i) REQUIRE(std::isfinite(termsF[i].filterEnergy()));
  }
}
#endif

#ifdef OD_HAVE_CUDA
// 9b pose energies on the GPU against the CPU float evaluation (spec 15.1.4).
TEST_CASE("gpu: 9b pose energies and minimization match the CPU path") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  auto gfield = std::make_shared<od::GpuField>(0);
  od::GpuPoseEvaluator ev(fs.x.mobSys->model(), gfield, 0, 4096);
  const int nt = fs.x.mobSys->torsionCount();
  for (const auto& [name, field] : fs.fields) {
    INFO(name);
    const od::MobmDocker dk(*fs.x.mobSys, field, fs.p);
    ev.bind(dk);
    const auto c = dk.conformation(fs.x.mobSys->inputAngles());
    const auto poses = filteredPoses(dk, c, 128);
    REQUIRE(poses.size() == 128);
    // Single evaluations.
    std::vector<std::vector<double>> angles;
    std::vector<od::PoseRequest> req;
    for (std::size_t i = 0; i < poses.size(); ++i) {
      auto a = perturbed(fs.x.mobSys->inputAngles(), 20, 40 + static_cast<unsigned>(i));
      for (double& v : a) v = static_cast<float>(v);
      angles.push_back(std::move(a));
    }
    for (std::size_t i = 0; i < poses.size(); ++i) req.push_back({poses[i], c.center, angles[i].data()});
    std::vector<od::MobmDocker::PoseEnergy> out(req.size());
    std::vector<double> grad(req.size() * (6 + nt)), ref;
    ev.evaluate(static_cast<int>(req.size()), req.data(), out.data(), grad.data());
    double worstE = 0, worstG = 0;
    for (std::size_t i = 0; i < req.size(); ++i) {
      const auto pe = dk.evaluate(c, poses[i], angles[i], &ref);
      CHECK(out[i].outside == pe.outside);
      worstE = std::max(worstE, std::abs(out[i].objective() - pe.objective()) / std::max(1.0, std::abs(pe.objective())));
      double gmax = 1;
      for (double v : ref) gmax = std::max(gmax, std::abs(v));
      for (std::size_t k = 0; k < ref.size(); ++k) worstG = std::max(worstG, std::abs(grad[i * (6 + nt) + k] - ref[k]) / gmax);
    }
    MESSAGE(std::format("{}: worst relative difference: energy {:.2e}, gradient {:.2e}", name, worstE, worstG));
    // Both are single precision with different rounding (torsion transforms on the GPU, spec 15.1.15).
    CHECK(worstE < 2e-4);
    CHECK(worstG < 1e-3);
    // Minimization: GPU energies against the CPU, with the CPU from starts shifted by 1e-5 Å as control.
    const int n = static_cast<int>(poses.size());
    std::vector<od::DockSample> gpuOut(n), cpuOut(n), ctlOut(n);
    std::vector<char> gpuOk(n), cpuOk(n), ctlOk(n);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n; ++i) {
      od::Dock9bStats s;
      cpuOk[i] = dk.minimize(c, poses[i], cpuOut[i], s);
      od::Pose shifted = poses[i];
      shifted.t[0] += 1e-5;
      ctlOk[i] = dk.minimize(c, shifted, ctlOut[i], s);
    }
    const std::size_t ng = 6 + nt;
    std::vector<od::PoseRequest> slotReq(n);
    std::vector<od::MobmDocker::PoseEnergy> slotPe(n);
    std::vector<double> slotGrad(n * ng);
    od::runBatched(
        n, n,
        [&](int i, int slot) {
          od::Dock9bStats st;
          gpuOk[i] = dk.minimize(c, poses[i], gpuOut[i], st, [&, slot](const od::Pose& p, std::span<const double> a, std::vector<double>* g) {
            slotReq[slot] = {p, c.center, a.data()};
            od::batchWait();
            g->assign(slotGrad.begin() + slot * ng, slotGrad.begin() + (slot + 1) * ng);
            return slotPe[slot];
          });
        },
        [&](std::span<const int> slots) {
          std::vector<od::PoseRequest> r;
          for (int s : slots) r.push_back(slotReq[s]);
          std::vector<od::MobmDocker::PoseEnergy> o(r.size());
          std::vector<double> gg(r.size() * ng);
          ev.evaluate(static_cast<int>(r.size()), r.data(), o.data(), gg.data());
          for (std::size_t k = 0; k < slots.size(); ++k) {
            slotPe[slots[k]] = o[k];
            std::copy_n(gg.begin() + k * ng, ng, slotGrad.begin() + slots[k] * ng);
          }
        });
    int agree = 0, ctlAgree = 0;
    double gpuMin = 1e300, cpuMin = 1e300;
    auto same = [](bool okA, const od::DockSample& x, bool okB, const od::DockSample& y) {
      if (okA != okB) return false;
      return !okA || std::abs(x.energy - y.energy) <= 1e-4 * std::max(1.0f, std::abs(y.energy));
    };
    for (int i = 0; i < n; ++i) {
      agree += same(gpuOk[i], gpuOut[i], cpuOk[i], cpuOut[i]);
      ctlAgree += same(ctlOk[i], ctlOut[i], cpuOk[i], cpuOut[i]);
      if (gpuOk[i]) gpuMin = std::min(gpuMin, static_cast<double>(gpuOut[i].energy));
      if (cpuOk[i]) cpuMin = std::min(cpuMin, static_cast<double>(cpuOut[i].energy));
    }
    MESSAGE(std::format("{}: {} poses minimized; same result as the CPU (1e-4 relative): GPU {}, CPU control {}; "
                        "lowest energy GPU {:.3f}, CPU {:.3f}",
                        name, n, agree, ctlAgree, gpuMin, cpuMin));
    CHECK(agree >= ctlAgree - n / 10);
    // Whole minimizations on the device (L-BFGS, spec 15.1.26).
    REQUIRE(fs.p.minimizer.method == od::MinMethod::Lbfgs);
    const od::DockMinParams mp = dk.deviceMinParams();
    std::vector<od::PoseRequest> startReq;
    for (int i = 0; i < n; ++i) startReq.push_back({poses[i], c.center, c.anglesDeg.data()});
    std::vector<od::DockMinResult> dres(n);
    std::vector<float> dang(static_cast<std::size_t>(n) * nt);
    REQUIRE(ev.minimize(n, startReq.data(), mp, dres.data(), dang.data()));
    int devAgree = 0, cpuSteps = 0, devSteps = 0;
    double devMin = 1e300;
    for (int i = 0; i < n; ++i) {
      od::DockSample ds;
      ds.energy = static_cast<float>(dres[i].dterms.inter + dres[i].dterms.interClash + dres[i].terms.vdw +
                                     dres[i].terms.hbond + dres[i].terms.elec + dres[i].terms.desolv + dres[i].terms.clash);
      devAgree += same(dres[i].ok != 0, ds, cpuOk[i], cpuOut[i]);
      if (dres[i].ok) {
        devMin = std::min(devMin, static_cast<double>(ds.energy));
        devSteps += dres[i].steps;
        // The returned angles and pose give the returned energy.
        od::Pose rp;
        for (int k = 0; k < 3; ++k) rp.t[k] = dres[i].t[k];
        rp.q = {dres[i].q[0], dres[i].q[1], dres[i].q[2], dres[i].q[3]};
        const std::vector<double> ra(dang.begin() + static_cast<std::ptrdiff_t>(i) * nt,
                                     dang.begin() + static_cast<std::ptrdiff_t>(i + 1) * nt);
        CHECK(std::abs(dk.evaluate(c, rp, ra).ranking() - ds.energy) <= 1e-3 * std::max(1.0f, std::abs(ds.energy)));
      }
      if (cpuOk[i]) cpuSteps += cpuOut[i].steps;
    }
    MESSAGE(std::format("{}: device minimizer: same result as the CPU {}, lowest energy {:.3f}; steps device {}, CPU {}",
                        name, devAgree, devMin, devSteps, cpuSteps));
    CHECK(devAgree >= ctlAgree - n / 10);
    // Memory above the device limit: minimized on the CPU (the evaluator declines).
    od::DockMinParams big = mp;
    big.memory = od::kDeviceLbfgsMaxMemory + 1;
    CHECK_FALSE(ev.minimize(n, startReq.data(), big, dres.data(), dang.data()));
  }
  ev.idle();
}

TEST_CASE("gpu: 9b device minimizer with each method matches the CPU") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  auto gfield = std::make_shared<od::GpuField>(0);
  od::GpuPoseEvaluator ev(fs.x.mobSys->model(), gfield, 0, 4096);
  const int nt = fs.x.mobSys->torsionCount();
  const auto& field = fs.fields[0].second;
  struct Case {
    const char *method, *init, *reduce;
    int maxBacktrack;
  };
  for (const Case cs : {Case{"sd", "unit", "halve", 20}, Case{"cg", "unit", "halve", 20}, Case{"fire", "unit", "halve", 20},
                        Case{"lbfgs", "unit", "halve", 20}, Case{"adadelta", "unit", "halve", 20},
                        Case{"lbfgs", "prev", "interp", 10}, Case{"lbfgs", "prev", "halve", 15},
                        Case{"cg", "unit", "interp", 10}, Case{"sd", "unit", "halve", 5}}) {
    const std::string method = std::format("{} {} {} {}", cs.method, cs.init, cs.reduce, cs.maxBacktrack);
    INFO(method);
    od::Dock9bParams p = fs.p;
    p.minimizer.method = od::parseMinMethod(cs.method);
    p.minimizer.lineSearchInit = od::parseLineSearchInit(cs.init);
    p.minimizer.lineSearchReduce = od::parseLineSearchReduce(cs.reduce);
    p.minimizer.maxBacktrack = cs.maxBacktrack;
    p.minimizer.maxSteps = 200;
    const od::MobmDocker dk(*fs.x.mobSys, field, p);
    ev.bind(dk);
    const auto c = dk.conformation(fs.x.mobSys->inputAngles());
    const auto poses = filteredPoses(dk, c, 128);
    REQUIRE(poses.size() == 128);
    const int n = static_cast<int>(poses.size());
    // The CPU, and the CPU from starts shifted by 1e-5 Å as control (spec 15.1.26).
    std::vector<od::DockSample> cpuOut(n), ctlOut(n);
    std::vector<char> cpuOk(n), ctlOk(n);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n; ++i) {
      od::Dock9bStats s;
      cpuOk[i] = dk.minimize(c, poses[i], cpuOut[i], s);
      od::Pose shifted = poses[i];
      shifted.t[0] += 1e-5;
      ctlOk[i] = dk.minimize(c, shifted, ctlOut[i], s);
    }
    std::vector<od::PoseRequest> req;
    for (int i = 0; i < n; ++i) req.push_back({poses[i], c.center, c.anglesDeg.data()});
    std::vector<od::DockMinResult> dres(n);
    std::vector<float> dang(static_cast<std::size_t>(n) * nt);
    REQUIRE(ev.minimize(n, req.data(), dk.deviceMinParams(), dres.data(), dang.data()));
    const auto same = [](bool okA, const od::DockSample& a, bool okB, const od::DockSample& b) {
      if (okA != okB) return 0;
      if (!okA) return 1;
      return std::abs(a.energy - b.energy) <= 1e-4 * std::max(1.0f, std::abs(b.energy)) ? 1 : 0;
    };
    int devAgree = 0, ctlAgree = 0, devSteps = 0, cpuSteps = 0, nBoth = 0;
    double devMin = 1e300, cpuMin = 1e300, ctlMin = 1e300, devSum = 0, cpuSum = 0, ctlSum = 0;
    for (int i = 0; i < n; ++i) {
      od::DockSample ds;
      ds.energy = static_cast<float>(dres[i].dterms.inter + dres[i].dterms.interClash + dres[i].terms.vdw +
                                     dres[i].terms.hbond + dres[i].terms.elec + dres[i].terms.desolv + dres[i].terms.clash);
      devAgree += same(dres[i].ok != 0, ds, cpuOk[i], cpuOut[i]);
      ctlAgree += same(ctlOk[i], ctlOut[i], cpuOk[i], cpuOut[i]);
      if (dres[i].ok) {
        devMin = std::min(devMin, static_cast<double>(ds.energy));
        devSteps += dres[i].steps;
        CHECK(dres[i].steps <= p.minimizer.maxSteps);
        // The returned angles and pose give the returned energy.
        od::Pose rp;
        for (int k = 0; k < 3; ++k) rp.t[k] = dres[i].t[k];
        rp.q = {dres[i].q[0], dres[i].q[1], dres[i].q[2], dres[i].q[3]};
        const std::vector<double> ra(dang.begin() + static_cast<std::ptrdiff_t>(i) * nt,
                                     dang.begin() + static_cast<std::ptrdiff_t>(i + 1) * nt);
        CHECK(std::abs(dk.evaluate(c, rp, ra).ranking() - ds.energy) <= 1e-3 * std::max(1.0f, std::abs(ds.energy)));
      }
      if (cpuOk[i]) {
        cpuMin = std::min(cpuMin, static_cast<double>(cpuOut[i].energy));
        cpuSteps += cpuOut[i].steps;
      }
      if (ctlOk[i]) ctlMin = std::min(ctlMin, static_cast<double>(ctlOut[i].energy));
      // Mean energy over the starts minimized by all three (each capped at 0: clashing starts).
      if (dres[i].ok && cpuOk[i] && ctlOk[i]) {
        ++nBoth;
        devSum += std::min(0.0, static_cast<double>(ds.energy));
        cpuSum += std::min(0.0, static_cast<double>(cpuOut[i].energy));
        ctlSum += std::min(0.0, static_cast<double>(ctlOut[i].energy));
      }
    }
    REQUIRE(nBoth > 0);
    const double devMean = devSum / nBoth, cpuMean = cpuSum / nBoth, ctlMean = ctlSum / nBoth;
    MESSAGE(std::format("{}: same result as the CPU: device {}, CPU control {}; lowest energy device {:.3f}, CPU {:.3f}, "
                        "control {:.3f}; mean energy device {:.3f}, CPU {:.3f}, control {:.3f}; steps device {}, CPU {}",
                        method, devAgree, ctlAgree, devMin, cpuMin, ctlMin, devMean, cpuMean, ctlMean, devSteps, cpuSteps));
    CHECK(devAgree >= ctlAgree - n / 10);
    // The lowest energy of 128 trajectories varies by more than 1 between the CPU and its control;
    // the mean is compared (as far from the CPU as the control, with a margin of 0.5).
    CHECK(std::abs(devMean - cpuMean) <= std::abs(ctlMean - cpuMean) + 0.5);
  }
  ev.idle();
}

TEST_CASE("gpu: 9b poses of several resident fields in one call") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  const int nf = static_cast<int>(fs.fields.size());
  auto gfield = std::make_shared<od::GpuField>(0);
  gfield->setSlots(nf);
  od::GpuPoseEvaluator ev(fs.x.mobSys->model(), gfield, 0, 4096);
  CHECK(ev.fields() == nf);
  const int nt = fs.x.mobSys->torsionCount();
  const std::size_t ng = 6 + static_cast<std::size_t>(nt);
  std::vector<std::unique_ptr<od::MobmDocker>> dks;
  std::vector<const od::MobmDocker*> ptr;
  for (const auto& [name, field] : fs.fields) {
    dks.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, field, fs.p));
    ptr.push_back(dks.back().get());
  }
  const auto c = dks[0]->conformation(fs.x.mobSys->inputAngles());
  const auto poses = filteredPoses(*dks[0], c, 96);
  REQUIRE(poses.size() == 96);
  std::vector<std::vector<double>> angles;
  for (std::size_t i = 0; i < poses.size(); ++i) {
    auto a = perturbed(fs.x.mobSys->inputAngles(), 20, 140 + static_cast<unsigned>(i));
    for (double& v : a) v = static_cast<float>(v);
    angles.push_back(std::move(a));
  }
  // Each field alone (one-field kernel).
  const std::size_t n = poses.size();
  std::vector<std::vector<od::MobmDocker::PoseEnergy>> single(nf, std::vector<od::MobmDocker::PoseEnergy>(n));
  std::vector<std::vector<double>> singleGrad(nf, std::vector<double>(n * ng));
  for (int f = 0; f < nf; ++f) {
    std::vector<od::PoseRequest> req;
    for (std::size_t i = 0; i < n; ++i) req.push_back({poses[i], c.center, angles[i].data()});
    ev.bind(*dks[f]);
    ev.evaluate(static_cast<int>(n), req.data(), single[f].data(), singleGrad[f].data());
  }
  CHECK(gfield->slotOf(fs.fields[0].second.id()) >= 0);  // all still resident
  // All fields in one call, pose i against field i % nf.
  ev.bindAll(ptr);
  std::vector<od::PoseRequest> req;
  for (std::size_t i = 0; i < n; ++i) req.push_back({poses[i], c.center, angles[i].data(), static_cast<int>(i % nf)});
  std::vector<od::MobmDocker::PoseEnergy> out(n);
  std::vector<double> grad(n * ng);
  ev.evaluate(static_cast<int>(n), req.data(), out.data(), grad.data());
  double worstE = 0, worstG = 0;
  int outsideDiff = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& s = single[i % nf][i];
    outsideDiff += out[i].outside != s.outside;
    worstE = std::max(worstE, std::abs(out[i].objective() - s.objective()) / std::max(1.0, std::abs(s.objective())));
    double gmax = 1;
    for (std::size_t k = 0; k < ng; ++k) gmax = std::max(gmax, std::abs(singleGrad[i % nf][i * ng + k]));
    for (std::size_t k = 0; k < ng; ++k)
      worstG = std::max(worstG, std::abs(grad[i * ng + k] - singleGrad[i % nf][i * ng + k]) / gmax);
  }
  MESSAGE(std::format("{} fields: worst relative difference to one field per call: energy {:.2e}, gradient {:.2e}", nf,
                      worstE, worstG));
  CHECK(outsideDiff == 0);
  // Same computation; only the order of the gradient atomics differs.
  CHECK(worstE < 1e-9);
  CHECK(worstG < 1e-6);
  // The fields differ: the poses against another field give other energies.
  int differ = 0;
  for (std::size_t i = 0; i < n; ++i) differ += single[(i + 1) % nf][i].objective() != out[i].objective();
  CHECK(differ > 0);
  ev.idle();
}

// gpu.packedMaps (spec 15.1.32), full grids only in device memory (spec 15.1.34) and minimizations
// divided among evaluators (spec 15.1.33; with one device, two evaluators on it): the same results.
TEST_CASE("gpu: packed reference maps, device maps and several evaluators give the same 9b minimizations") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  const int nt = fs.x.mobSys->torsionCount();
  const int dev1 = od::gpuDeviceCount() > 1 ? 1 : 0;
  for (const auto& [name, field] : fs.fields) {
    INFO(name);
    od::Dock9bParams p = fs.p;
    p.minimizer.maxSteps = 100;
    const od::MobmDocker dk(*fs.x.mobSys, field, p);
    const auto c = dk.conformation(fs.x.mobSys->inputAngles());
    const auto poses = filteredPoses(dk, c, 301);
    const int n = static_cast<int>(poses.size());
    std::vector<od::PoseRequest> req;
    for (int i = 0; i < n; ++i) req.push_back({poses[i], c.center, c.anglesDeg.data()});
    const auto run = [&](od::PoseEvaluator& ev, std::vector<od::DockMinResult>& r, std::vector<float>& a) {
      r.assign(n, {});
      a.assign(static_cast<std::size_t>(n) * nt, 0.0f);
      ev.bind(dk);
      REQUIRE(ev.minimize(n, req.data(), dk.deviceMinParams(), r.data(), a.data()));
      ev.idle();
    };
    std::vector<od::DockMinResult> r0, r1, r2, r3, r4;
    std::vector<float> a0, a1, a2, a3, a4;
    {
      od::GpuPoseEvaluator ev(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, 4096);
      run(ev, r0, a0);
    }
    {
      auto f = std::make_shared<od::GpuField>(0);
      f->setPackedMaps(od::GpuField::Packed::On);
      od::GpuPoseEvaluator ev(fs.x.mobSys->model(), f, 0, 4096);
      run(ev, r1, a1);
    }
    {
      std::vector<std::unique_ptr<od::GpuPoseEvaluator>> evs;
      evs.push_back(std::make_unique<od::GpuPoseEvaluator>(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, 4096));
      evs.push_back(std::make_unique<od::GpuPoseEvaluator>(fs.x.mobSys->model(), std::make_shared<od::GpuField>(dev1), dev1, 4096));
      od::MultiGpuPoseEvaluator ev(std::move(evs));
      run(ev, r2, a2);
    }
    // full: the maps only in device memory (spec 15.1.34), used in place on their device and
    // copied (cudaMemcpyPeer) to another one.
    if (name == "full") {
      od::GpuGridBuilder gb(0);
      const od::GridBuildData b(field.geom(), fs.x.stam, fs.x.mobm.type, {});
      auto g = std::make_shared<od::EnergyGrids>();
      g->geom = field.geom();
      g->types = b.types;
      auto d = std::make_shared<od::DistanceGrid>();
      d->geom = field.geom();
      d->cap = field.distanceCap();
      const od::StamField devField = od::StamField::full(
          g, d, gb.buildDevice(b.view, od::DistBuildData(field.geom(), fs.x.stam, field.distanceCap()).view, od::kGridClamp));
      CHECK_FALSE(devField.hasHostData());
      CHECK(devField.model().view.maps[0] == nullptr);
      const od::MobmDocker dkDev(*fs.x.mobSys, devField, p);
      const auto runDev = [&](int dev, std::vector<od::DockMinResult>& r, std::vector<float>& a) {
        od::GpuPoseEvaluator ev(fs.x.mobSys->model(), std::make_shared<od::GpuField>(dev), dev, 4096);
        r.assign(n, {});
        a.assign(static_cast<std::size_t>(n) * nt, 0.0f);
        ev.bind(dkDev);
        REQUIRE(ev.minimize(n, req.data(), dkDev.deviceMinParams(), r.data(), a.data()));
        ev.idle();
      };
      runDev(0, r3, a3);
      runDev(dev1, r4, a4);
    }
    // Packed maps: the same values in another layout (identical results up to the order of the
    // gradient atomics); the divided call: each minimization as alone.
    int same1 = 0, same2 = 0, same3 = 0, same4 = 0;
    for (int i = 0; i < n; ++i) {
      const auto e = [](const od::DockMinResult& r) { return r.dterms.inter + r.terms.vdw + r.terms.hbond + r.terms.elec + r.terms.desolv; };
      CHECK(r1[i].ok == r0[i].ok);
      CHECK(r2[i].ok == r0[i].ok);
      if (!r0[i].ok) continue;
      same1 += std::abs(e(r1[i]) - e(r0[i])) <= 1e-3 * std::max(1.0, std::abs(e(r0[i])));
      same2 += std::abs(e(r2[i]) - e(r0[i])) <= 1e-3 * std::max(1.0, std::abs(e(r0[i])));
      if (name != "full") continue;
      CHECK(r3[i].ok == r0[i].ok);
      CHECK(r4[i].ok == r0[i].ok);
      same3 += std::abs(e(r3[i]) - e(r0[i])) <= 1e-3 * std::max(1.0, std::abs(e(r0[i])));
      same4 += std::abs(e(r4[i]) - e(r0[i])) <= 1e-3 * std::max(1.0, std::abs(e(r0[i])));
    }
    int ok = 0;
    for (const auto& r : r0) ok += r.ok;
    MESSAGE(std::format("{}: {} minimized; same energy: packed {}, divided {}, device maps {}, device maps on device {} {}",
                        name, ok, same1, same2, same3, dev1, same4));
    CHECK(same1 >= ok * 9 / 10);
    CHECK(same2 >= ok * 9 / 10);
    if (name == "full") {
      CHECK(same3 >= ok * 9 / 10);
      CHECK(same4 >= ok * 9 / 10);
    }
  }
}

TEST_CASE("gpu: 9b quick filter and posm2_fft match the CPU path") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  auto gfield = std::make_shared<od::GpuField>(0);
  auto prep = std::make_shared<od::GpuDockPrep>(0);
  od::GpuPoseEvaluator ev(fs.x.mobSys->model(), gfield, 0, 1024, prep);
  const auto angles = fs.x.mobSys->inputAngles();
  // Quick filter (spec 8.5): the same L2poss and counters as the CPU.
  for (od::Posm posm : {od::Posm::Posm1, od::Posm::Posm2}) {
    for (const auto& [name, field] : fs.fields) {
      INFO(name);
      od::Dock9bParams p = fs.p;
      p.posm = posm;
      p.posStep = 2.0;
      const od::MobmDocker dk(*fs.x.mobSys, field, p);
      const auto c = dk.conformation(angles);
      const long long n = dk.l1Count(c);
      od::Dock9bStats cs, gs;
      std::vector<long long> cpu, gpu;
      for (long long k = 0; k < n; ++k) {
        const auto pose = dk.l1Pose(c, 17, k);
        if (!pose) ++cs.outsideBox;
        else if (dk.quickFilter(c, *pose, cs)) cpu.push_back(k);
      }
      REQUIRE(ev.quickFilter(dk, c, 17, 0, n, gpu, gs));
      std::sort(gpu.begin(), gpu.end());
      MESSAGE(std::format("{} {}: L1poss {}, L2poss CPU {} GPU {}; outside {}/{}, overlap {}/{}, far {}/{}",
                          posm == od::Posm::Posm1 ? "posm1" : "posm2", name, n, cpu.size(), gpu.size(), cs.outsideBox,
                          gs.outsideBox, cs.overlap, gs.overlap, cs.far, gs.far));
      CHECK(!cpu.empty());
      CHECK(gpu == cpu);
      CHECK(gs.outsideBox == cs.outsideBox);
      CHECK(gs.overlap == cs.overlap);
      CHECK(gs.far == cs.far);
    }
  }
  // posm3 and posm4: the poses generated on the device from the rotations (spec 15.1.21).
  for (od::Posm posm : {od::Posm::Posm3, od::Posm::Posm4}) {
    for (const auto& [name, field] : fs.fields) {
      INFO(name);
      od::Dock9bParams p = fs.p;
      p.posm = posm;
      p.boxSize = {34, 34, 34};
      p.posm1SampleNum = 300;
      p.sepaSingle = true;
      p.sepaPotPercentile = 1.0;
      p.sepaChargeThreshold = 0.5;
      p.maxOverlapAtomNum = 6;
      const auto sepas = od::stamSepas(field, p, nullptr);
      const od::MobmDocker dk(*fs.x.mobSys, field, p, &sepas);
      const auto c = dk.conformation(angles);
      const long long n = dk.l1Count(c);
      od::Dock9bStats cs, gs;
      std::vector<long long> cpu, gpu;
      for (long long k = 0; k < n; ++k) {
        const auto pose = dk.l1Pose(c, 17, k);
        if (!pose) ++cs.outsideBox;
        else if (dk.quickFilter(c, *pose, cs)) cpu.push_back(k);
      }
      // Part of the range: not aligned to the rotations.
      const long long k0 = 5 + c.sepaPoses() / 3;
      REQUIRE(ev.quickFilter(dk, c, 17, 0, k0, gpu, gs));
      REQUIRE(ev.quickFilter(dk, c, 17, k0, n, gpu, gs));
      std::sort(gpu.begin(), gpu.end());
      MESSAGE(std::format("{} {}: L1poss {}, L2poss CPU {} GPU {}; outside {}/{}, overlap {}/{}, far {}/{}",
                          posm == od::Posm::Posm3 ? "posm3" : "posm4", name, n, cpu.size(), gpu.size(), cs.outsideBox,
                          gs.outsideBox, cs.overlap, gs.overlap, cs.far, gs.far));
      CHECK(!cpu.empty());
      CHECK(cs.outsideBox > 0);
      CHECK(gpu == cpu);
      CHECK(gs.outsideBox == cs.outsideBox);
      CHECK(gs.overlap == cs.overlap);
      CHECK(gs.far == cs.far);
    }
  }
  // posm2_fft (spec 8.4.3): energies of the kept translations as the CPU (cuFFT vs pocketfft rounding),
  // with the direct and the FFT probe spectra.
  for (auto mode : {od::GpuDockPrep::Spectra::Direct, od::GpuDockPrep::Spectra::Fft})
  for (int top : {1, 3, 20}) {
    prep->setSpectra(mode);
    const char* mname = mode == od::GpuDockPrep::Spectra::Direct ? "direct" : "fft";
    for (const auto& [name, field] : fs.fields) {
      INFO(name);
      INFO(top);
      INFO(mname);
      const od::FftStam s = od::buildFftStam(field, 10.0);
      od::Dock9bParams p = fs.p;
      p.posm = od::Posm::Posm2Fft;
      p.fftTopPerRot = top;
      const od::MobmDocker dkc(*fs.x.mobSys, field, p, nullptr, &s);
      p.evaluator = &ev;
      const od::MobmDocker dkg(*fs.x.mobSys, field, p, nullptr, &s);
      const auto cc = dkc.conformation(angles);
      const auto cg = dkg.conformation(angles);
      REQUIRE(cg.fftPoses.size() == cc.fftPoses.size());
      REQUIRE(!cc.fftPoses.empty());
      std::size_t same = 0;
      double worst = 0;
      for (std::size_t k = 0; k < cc.fftPoses.size(); ++k) {
        same += cg.fftPoses[k].t == cc.fftPoses[k].t;
        CHECK(cg.fftPoses[k].q.w == cc.fftPoses[k].q.w);
        worst = std::max(worst, std::abs(static_cast<double>(cg.fftEnergy[k]) - cc.fftEnergy[k]));
        CHECK(cg.fftEnergy[k] == doctest::Approx(cc.fftEnergy[k]).epsilon(1e-4).scale(50.0));
      }
      MESSAGE(std::format("posm2_fft {} spectra {} top {}: kept {}, same translation {}, max |GPU - CPU| energy {:.2e}", name,
                          mname, top, cc.fftPoses.size(), same, worst));
      CHECK(same >= cc.fftPoses.size() * 9 / 10);
    }
  }
  // Several StaM conformations at once (spec 15.1.22): as each alone on the CPU; fftStride 2 (15.1.23).
  for (const auto& [top, stride] : {std::pair{1, 1}, std::pair{3, 1}, std::pair{20, 1}, std::pair{1, 2}}) {
    INFO(top);
    INFO(stride);
    std::vector<od::FftStam> ss;
    for (const auto& [name, field] : fs.fields) ss.push_back(od::buildFftStam(field, 10.0, stride));
    od::Dock9bParams p = fs.p;
    p.posm = od::Posm::Posm2Fft;
    p.fftTopPerRot = top;
    od::Dock9bParams pg = p;
    pg.evaluator = &ev;
    std::vector<std::unique_ptr<od::MobmDocker>> own;
    std::vector<const od::MobmDocker*> dkc, dkg;
    for (std::size_t s = 0; s < ss.size(); ++s) {
      own.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, fs.fields[s].second, p, nullptr, &ss[s]));
      dkc.push_back(own.back().get());
      own.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, fs.fields[s].second, pg, nullptr, &ss[s]));
      dkg.push_back(own.back().get());
    }
    const auto cc = od::MobmDocker::conformations(dkc, angles);
    const auto cg = od::MobmDocker::conformations(dkg, angles);
    for (std::size_t s = 0; s < ss.size(); ++s) {
      INFO(fs.fields[s].first);
      REQUIRE(cg[s].fftPoses.size() == cc[s].fftPoses.size());
      REQUIRE(!cc[s].fftPoses.empty());
      std::size_t same = 0;
      double worst = 0;
      for (std::size_t k = 0; k < cc[s].fftPoses.size(); ++k) {
        same += cg[s].fftPoses[k].t == cc[s].fftPoses[k].t;
        CHECK(cg[s].fftPoses[k].q.w == cc[s].fftPoses[k].q.w);
        worst = std::max(worst, std::abs(static_cast<double>(cg[s].fftEnergy[k]) - cc[s].fftEnergy[k]));
        CHECK(cg[s].fftEnergy[k] == doctest::Approx(cc[s].fftEnergy[k]).epsilon(1e-4).scale(50.0));
      }
      MESSAGE(std::format("posm2_fft {} StaM sets, {} top {} stride {}: kept {}, same translation {}, max |GPU - CPU| energy {:.2e}",
                          ss.size(), fs.fields[s].first, top, stride, cc[s].fftPoses.size(), same, worst));
      CHECK(same >= cc[s].fftPoses.size() * 9 / 10);
    }
    // The orientations divided among two evaluators (spec 15.1.35): the same as one.
    {
      const int dev1 = od::gpuDeviceCount() > 1 ? 1 : 0;
      // The probe spectra as ev (FFT, the last mode set above).
      std::vector<std::unique_ptr<od::GpuPoseEvaluator>> evs;
      for (int d : {0, dev1}) {
        auto pr = std::make_shared<od::GpuDockPrep>(d);
        pr->setSpectra(od::GpuDockPrep::Spectra::Fft);
        evs.push_back(std::make_unique<od::GpuPoseEvaluator>(fs.x.mobSys->model(), std::make_shared<od::GpuField>(d), d, 1024, pr));
      }
      od::MultiGpuPoseEvaluator multi(std::move(evs));
      od::Dock9bParams pm = p;
      pm.evaluator = &multi;
      std::vector<const od::MobmDocker*> dkm;
      for (std::size_t s = 0; s < ss.size(); ++s) {
        own.push_back(std::make_unique<od::MobmDocker>(*fs.x.mobSys, fs.fields[s].second, pm, nullptr, &ss[s]));
        dkm.push_back(own.back().get());
      }
      const auto cm = od::MobmDocker::conformations(dkm, angles);
      for (std::size_t s = 0; s < ss.size(); ++s) {
        REQUIRE(cm[s].fftPoses.size() == cg[s].fftPoses.size());
        std::size_t same = 0;
        for (std::size_t k = 0; k < cg[s].fftPoses.size(); ++k)
          same += cm[s].fftPoses[k].t == cg[s].fftPoses[k].t && cm[s].fftPoses[k].q.w == cg[s].fftPoses[k].q.w &&
                  cm[s].fftEnergy[k] == cg[s].fftEnergy[k];
        CHECK(same == cg[s].fftPoses.size());
      }
      multi.idle();
    }
  }
  ev.idle();
}

// Minimizations taken over by another device (spec 15.1.35): the poses given away by a running
// kernel (GpuEnergy::stealDock) and minimized elsewhere give the results of one call; with two
// devices, DockStealHub between evaluators.
TEST_CASE("gpu: 9b minimizations taken over by another device give the same results") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  FieldSet fs;
  const int nt = fs.x.mobSys->torsionCount();
  const auto& field = fs.fields[0].second;
  od::Dock9bParams p = fs.p;
  p.minimizer.maxSteps = 300;
  const od::MobmDocker dk(*fs.x.mobSys, field, p);
  const auto c = dk.conformation(fs.x.mobSys->inputAngles());
  const auto base = filteredPoses(dk, c, 400);
  // Enough poses for the kernel to run while the requests arrive.
  std::vector<od::PoseRequest> req;
  for (int rep = 0; rep < 40; ++rep)
    for (const auto& ps : base) req.push_back({ps, c.center, c.anglesDeg.data()});
  const int n = static_cast<int>(req.size());
  const auto e = [](const od::DockMinResult& r) {
    return static_cast<double>(r.terms.vdw + r.terms.hbond + r.terms.elec + r.terms.desolv + r.dterms.inter);
  };
  std::vector<od::DockMinResult> r0(n);
  std::vector<float> a0(static_cast<std::size_t>(n) * nt);
  {
    od::GpuPoseEvaluator ev(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, n);
    ev.bind(dk);
    REQUIRE(ev.minimize(n, req.data(), dk.deviceMinParams(), r0.data(), a0.data()));
  }
  const auto compare = [&](const std::vector<od::DockMinResult>& r, const char* what) {
    int ok = 0, same = 0;
    for (int i = 0; i < n; ++i) {
      CHECK(r[i].ok == r0[i].ok);
      ok += r0[i].ok != 0;
      same += r0[i].ok && std::abs(e(r[i]) - e(r0[i])) <= 1e-3 * std::max(1.0, std::abs(e(r0[i])));
    }
    MESSAGE(std::format("{}: {} minimized, same energy {}", what, ok, same));
    CHECK(same >= ok * 9 / 10);
  };
  // GpuEnergy: a request while the kernel runs gives away a tail; [0, end) from the call.
  {
    od::GpuPoseEvaluator owner(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, n);
    od::GpuPoseEvaluator helper(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, n);
    owner.bind(dk);
    helper.bind(dk);
    // Nothing runs: no poses.
    CHECK(owner.energy().stealDock(0.5, 1) == std::pair{0, 0});
    std::vector<od::DockMinResult> r(n);
    std::vector<float> a(static_cast<std::size_t>(n) * nt);
    auto hub = std::make_shared<od::DockStealHub>();
    owner.setSteal(hub, 0);
    std::pair<int, int> got{0, 0};
    std::thread th([&] {
      int taken = 0, end = 0;
      double sec = 0;
      while (!owner.energy().dockProgress(taken, end, sec) || taken == 0) std::this_thread::sleep_for(std::chrono::microseconds(50));
      CHECK(end == n);
      got = owner.energy().stealDock(0.5, 1);
      if (got.first < got.second)
        REQUIRE(helper.minimize(got.second - got.first, req.data() + got.first, dk.deviceMinParams(), r.data() + got.first,
                                a.data() + static_cast<std::size_t>(got.first) * nt));
    });
    REQUIRE(owner.minimize(n, req.data(), dk.deviceMinParams(), r.data(), a.data()));
    th.join();
    MESSAGE(std::format("stealDock: poses [{}, {}) of {}", got.first, got.second, n));
    CHECK(got.first > 0);
    CHECK(got.second == n);
    CHECK(got.first < n);
    compare(r, "stealDock");
  }
  if (od::gpuDeviceCount() < 2) {
    MESSAGE("one CUDA device; DockStealHub between devices skipped");
    return;
  }
  // DockStealHub: device 1 takes over poses of the call on device 0.
  {
    auto hub = std::make_shared<od::DockStealHub>();
    od::GpuPoseEvaluator owner(fs.x.mobSys->model(), std::make_shared<od::GpuField>(0), 0, n);
    od::GpuPoseEvaluator helper(fs.x.mobSys->model(), std::make_shared<od::GpuField>(1), 1, n);
    owner.setSteal(hub, 0);
    helper.setSteal(hub, 0);
    owner.bind(dk);
    std::vector<od::DockMinResult> r(n);
    std::vector<float> a(static_cast<std::size_t>(n) * nt);
    std::atomic<bool> done{false};
    std::thread th([&] {
      while (!done.load())
        if (!hub->helpOnce(1, [&](int) { return &helper; })) hub->wait(1);
    });
    REQUIRE(owner.minimize(n, req.data(), dk.deviceMinParams(), r.data(), a.data()));
    done = true;
    th.join();
    const auto st = hub->stats();
    MESSAGE(std::format("DockStealHub: {} poses taken over in {} part(s)", st.poses, st.steals));
    CHECK(st.steals >= 1);
    CHECK(st.poses > 0);
    CHECK(!hub->anyCall());
    compare(r, "DockStealHub");
  }
}

TEST_CASE("gpu: energy grids equal the CPU grids") {
  if (od::gpuDeviceCount() == 0) {
    MESSAGE("no CUDA device; skipped");
    return;
  }
  GridFixture x;
  const od::Vec3 c = x.stamCenter();
  const od::EnergyParams params;
  for (double size : {12.0, 26.0}) {
    INFO("box ", size);
    const auto geom = od::boxGrid(c, {size, size, size * 0.75}, 0.375);
    const auto ref = od::buildEnergyGrids(geom, x.stam, x.mobm.type, params);
    // The default slab and slabs that split the grid unevenly.
    for (std::size_t slab : {std::size_t(0), std::size_t(7919)}) {
      INFO("slab ", slab);
      od::GpuGridBuilder gb(0, slab);
      const od::GridBuildData b(geom, x.stam, x.mobm.type, params);
      REQUIRE(b.types == ref.types);
      std::vector<std::vector<float>> maps(b.types.size(), std::vector<float>(geom.size(), -1.0f));
      std::vector<float*> ptr;
      for (auto& m : maps) ptr.push_back(m.data());
      std::vector<float> elec(geom.size(), -1.0f), desolv(geom.size(), -1.0f);
      gb.build(b.view, od::kGridClamp, ptr.data(), elec.data(), desolv.data());
      std::size_t diff = 0;
      for (std::size_t t = 0; t < maps.size(); ++t) diff += maps[t] != ref.typeMaps[t];
      CHECK(diff == 0);
      CHECK(elec == ref.elecMap);
      CHECK(desolv == ref.desolvMap);
      // Distance grid.
      const double cap = 4.25;
      const auto refDist = od::buildDistanceGrid(geom, x.stam, cap);
      std::vector<float> dist(geom.size(), -1.0f);
      gb.buildDistance(od::DistBuildData(geom, x.stam, cap).view, dist.data());
      CHECK(dist == refDist.d);
      // Full grids kept on the device (spec 15.1.34): the same values in the FieldView::maps layout.
      {
        const auto dm = gb.buildDevice(b.view, od::DistBuildData(geom, x.stam, cap).view, od::kGridClamp);
        REQUIRE(dm->floats == (maps.size() + 3) * geom.size());
        CHECK(dm->device == 0);
        std::vector<float> all(dm->floats), want;
        od::gpuDownload(dm->ptr, all.data(), all.size());
        for (const auto& m : ref.typeMaps) want.insert(want.end(), m.begin(), m.end());
        want.insert(want.end(), ref.elecMap.begin(), ref.elecMap.end());
        want.insert(want.end(), ref.desolvMap.begin(), ref.desolvMap.end());
        want.insert(want.end(), refDist.d.begin(), refDist.d.end());
        CHECK(all == want);
      }
      // delta blocks of a sample with moved atoms against the CPU StamField::delta.
      od::GridAtoms moved = x.stam;
      for (std::size_t a = 0; a < moved.size(); a += 7) moved.xyz[a][0] += 0.8;
      const auto sg = od::buildEnergyGrids(geom, moved, x.mobm.type, params);
      const auto sd = od::buildDistanceGrid(geom, moved, cap);
      auto rg = std::make_shared<const od::EnergyGrids>(ref);
      auto rd = std::make_shared<const od::DistanceGrid>(refDist);
      std::vector<const float*> rmaps;
      for (const auto& m : ref.typeMaps) rmaps.push_back(m.data());
      rmaps.push_back(ref.elecMap.data());
      rmaps.push_back(ref.desolvMap.data());
      rmaps.push_back(refDist.d.data());
      gb.setReference(42, rmaps.data(), static_cast<int>(rmaps.size()), geom.size());
      CHECK(gb.referenceId() == 42);
      for (double thr : {1e-4, 1e-2, 0.5}) {
        INFO("threshold ", thr);
        const auto cpu = od::StamField::delta(rg, rd, sg, sd, thr);
        std::vector<int> block;
        std::vector<float> data;
        gb.buildDelta(od::GridBuildData(geom, moved, x.mobm.type, params).view, od::DistBuildData(geom, moved, cap).view,
                      od::kGridClamp, thr, block, data);
        const auto gpu = od::StamField::delta(rg, rd, std::move(block), std::move(data));
        CHECK(cpu.deltaBlockCount() > 0);
        if (thr == 0.5) CHECK(cpu.deltaBlockCount() < cpu.blockCount());
        const od::FieldModel a = cpu.model(), b2 = gpu.model();
        CHECK(a.block == b2.block);
        REQUIRE(a.data);
        REQUIRE(b2.data);
        CHECK(*a.data == *b2.data);
        // Only the blocks reached by the changed positions (spec 15.1.19); device-only blocks.
        const std::vector<od::Vec3> changed = od::changedPositions(x.stam, moved);
        CHECK(changed.size() > 0);
        CHECK(changed.size() < 2 * moved.size());
        const auto cpuR = od::StamField::delta(rg, rd, sg, sd, thr, &changed);
        std::vector<int> blockR;
        std::vector<float> dataR;
        auto dev = gb.buildDelta(od::GridBuildData(geom, moved, x.mobm.type, params).view,
                                 od::DistBuildData(geom, moved, cap).view, od::kGridClamp, thr, changed[0].data(),
                                 static_cast<int>(changed.size()), blockR, &dataR);
        std::vector<int> blockD;
        auto devOnly = gb.buildDelta(od::GridBuildData(geom, moved, x.mobm.type, params).view,
                                     od::DistBuildData(geom, moved, cap).view, od::kGridClamp, thr, changed[0].data(),
                                     static_cast<int>(changed.size()), blockD, nullptr);
        CHECK(blockD == blockR);
        // submitDelta / finishDelta (spec 15.1.20): the inputs may be released after submitDelta.
        {
          std::vector<int> blockS;
          std::vector<float> dataS;
          {
            const od::GridBuildData gbd(geom, moved, x.mobm.type, params);
            const od::DistBuildData dbd(geom, moved, cap);
            const std::vector<od::Vec3> ch = changed;
            gb.submitDelta(gbd.view, dbd.view, od::kGridClamp, thr, ch[0].data(), static_cast<int>(ch.size()));
          }
          CHECK_THROWS(gb.submitDelta(od::GridBuildData(geom, moved, x.mobm.type, params).view,
                                      od::DistBuildData(geom, moved, cap).view, od::kGridClamp, thr, nullptr, -1));
          gb.finishDelta(blockS, &dataS);
          CHECK(blockS == blockR);
          CHECK(dataS == dataR);
          CHECK_THROWS(gb.finishDelta(blockS, &dataS));
        }
        REQUIRE(dev);
        REQUIRE(devOnly);
        CHECK(dev->floats == dataR.size());
        std::vector<float> back(devOnly->floats);
        if (!back.empty()) od::gpuDownload(devOnly->ptr, back.data(), back.size());
        CHECK(back == dataR);
        const auto gpuR = od::StamField::delta(rg, rd, std::move(blockR), std::move(dataR), dev);
        const od::FieldModel ar = cpuR.model(), br = gpuR.model();
        CHECK(ar.block == br.block);
        CHECK(*ar.data == *br.data);
        CHECK(br.deviceData == dev);
        // Reached blocks include all blocks that differ beyond the reach of rounding.
        CHECK(cpuR.deltaBlockCount() <= cpu.deltaBlockCount());
        CHECK(cpuR.deltaBlockCount() >= cpu.deltaBlockCount() * 9 / 10);
        const auto gpuD = od::StamField::delta(rg, rd, std::move(blockD), {}, devOnly);
        CHECK_FALSE(gpuD.hasHostData());
        CHECK(gpuD.deltaBlockCount() == cpuR.deltaBlockCount());
        CHECK_THROWS_AS(gpuD.pointMap(0), od::Error);
      }
    }
  }
}
#endif

#ifdef OD_HAVE_CUDA
template <class T>
static double zbenchRun(od::GpuEnergy& ev, const std::vector<T>& a, int n, std::vector<od::FlatTerms>& t,
                        std::vector<double>& g, int reps) {
  ev.evaluate(a.data(), n, t.data(), g.data(), static_cast<T*>(nullptr));
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) ev.evaluate(a.data(), n, t.data(), g.data(), static_cast<T*>(nullptr));
  return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / reps;
}
template <class T>
static void zbenchAll(const char* label) {
  FlatFixture f;
  const char* only = std::getenv("ZB_CASE");
  const char* ns = std::getenv("ZB_N");
  for (const auto& c : f.cases) {
    if (only && c.name != only) continue;
    const od::EnergySystem& sys = *c.sys;
    const od::EnergyModel model = sys.model();
    const int d = sys.torsionCount();
    std::string line = std::format("{} {} (tors {}, pairs {}, coords {}):", label, c.name, d, model.pairs.size(), model.coordCount);
    std::vector<int> nl{1, 32, 256, 4096};
    if (ns) nl = {std::atoi(ns)};
    for (int n : nl) {
      std::vector<T> a;
      for (int i = 0; i < n; ++i) {
        const auto v = perturbed(sys.inputAngles(), 30, 900 + i);
        a.insert(a.end(), v.begin(), v.end());
      }
      od::GpuEnergy ev(model, 0);
      std::vector<od::FlatTerms> t(n);
      std::vector<double> g(static_cast<std::size_t>(n) * d);
      line += std::format("  n={} {:.0f} us", n, zbenchRun(ev, a, n, t, g, ns ? 3 : (n >= 4096 ? 10 : 200)));
    }
    MESSAGE(line);
  }
}
// Time per energy evaluation of the 9b device minimizer (each field of FieldSet, lbfgs and
// adadelta, 8192 poses, 200 steps); ZB_REPS repetitions (default 5).
TEST_CASE("zbench: 9b device minimizer" * doctest::skip()) {
  if (od::gpuDeviceCount() == 0) return;
  FieldSet fs;
  auto gfield = std::make_shared<od::GpuField>(0);
  od::GpuPoseEvaluator ev(fs.x.mobSys->model(), gfield, 0, 8192);
  const int nt = fs.x.mobSys->torsionCount(), reps = std::getenv("ZB_REPS") ? std::atoi(std::getenv("ZB_REPS")) : 5;
  for (const auto& [name, field] : fs.fields) {
    for (const char* method : {"lbfgs", "adadelta"}) {
      od::Dock9bParams p = fs.p;
      p.minimizer.method = od::parseMinMethod(method);
      p.minimizer.maxSteps = 200;
      const od::MobmDocker dk(*fs.x.mobSys, field, p);
      ev.bind(dk);
      const auto c = dk.conformation(fs.x.mobSys->inputAngles());
      auto poses = filteredPoses(dk, c, 2048);
      while (poses.size() < 8192) poses.insert(poses.end(), poses.begin(), poses.begin() + std::min<std::size_t>(poses.size(), 8192 - poses.size()));
      const int n = static_cast<int>(poses.size());
      std::vector<od::PoseRequest> req;
      for (int i = 0; i < n; ++i) req.push_back({poses[i], c.center, c.anglesDeg.data()});
      std::vector<od::DockMinResult> dres(n);
      std::vector<float> dang(static_cast<std::size_t>(n) * nt);
      REQUIRE(ev.minimize(n, req.data(), dk.deviceMinParams(), dres.data(), dang.data()));
      double best = 1e300;
      long long evals = 0;
      for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        ev.minimize(n, req.data(), dk.deviceMinParams(), dres.data(), dang.data());
        best = std::min(best, std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count());
      }
      for (const auto& d : dres) evals += d.evaluations;
      MESSAGE(std::format("{} {}: {:.0f} us, {} evaluations, {:.4f} us per evaluation", name, method, best, evals, best / evals));
    }
  }
  ev.idle();
}

TEST_CASE("zbench: latency" * doctest::skip()) {
  zbenchAll<float>("float");
  if (!std::getenv("ZB_FLOAT")) zbenchAll<double>("double");
}
#endif

TEST_CASE("score: pose records, reference atom map, Kabsch superposition") {
  const std::string sdf = kData + "/BTN_ideal.sdf";
  std::ifstream in(sdf);
  const std::string rec((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  const auto tmp = std::filesystem::temp_directory_path();
  // Three records: the input, the same atoms moved by +1 Å along x, and a record with one atom less.
  std::string moved;
  {
    std::istringstream ss(rec);
    std::string line;
    int n = 0, k = 0;
    while (std::getline(ss, line)) {
      ++k;
      if (k == 4) n = std::stoi(line.substr(0, 3));
      if (k > 4 && k <= 4 + n) {
        const double x = std::stod(line.substr(0, 10)) + 1.0;
        line = std::format("{:10.4f}", x) + line.substr(10);
      }
      moved += line + "\n";
    }
  }
  const std::string two = (tmp / std::format("od_test_poses_{}.sdf", ::getpid())).string();
  std::ofstream(two) << rec << moved;
  const od::PoseRecords p = od::loadPoseRecords(two, od::Role::MobM, ccd());
  REQUIRE(p.xyz.size() == 2);
  REQUIRE(p.xyz[1].size() == p.mol.atoms.size());
  for (std::size_t i = 0; i < p.mol.atoms.size(); ++i) {
    CHECK(p.xyz[0][i][0] == doctest::Approx(p.mol.atoms[i].xyz[0]));
    CHECK(p.xyz[1][i][0] == doctest::Approx(p.mol.atoms[i].xyz[0] + 1.0));
  }
  // A second record with another element for atom 1.
  const std::string bad = (tmp / std::format("od_test_poses_bad_{}.sdf", ::getpid())).string();
  {
    std::istringstream ss(rec);
    std::string line, out;
    int k = 0;
    while (std::getline(ss, line)) {
      ++k;
      if (k == 5) line.replace(31, 3, line.substr(31, 1) == "O" ? "C  " : "O  ");
      out += line + "\n";
    }
    std::ofstream(bad) << rec << out;
  }
  CHECK_THROWS_AS(od::loadPoseRecords(bad, od::Role::MobM, ccd()), od::Error);

  // Reference in another atom order: the map recovers the correspondence.
  const od::Molecule mol = od::loadMolecules(sdf, od::Role::MobM, ccd()).front();
  std::vector<int> heavy;
  for (std::size_t i = 0; i < mol.atoms.size(); ++i)
    if (mol.atoms[i].atomicNum != 1) heavy.push_back(static_cast<int>(i));
  od::Molecule ref = od::loadReference(sdf, ccd());
  {
    const auto same = od::referenceAtomMap(mol, heavy, ref);
    CHECK(same == heavy);
  }
  const int n = static_cast<int>(ref.atoms.size());
  std::vector<int> perm(n);  // new index of old atom i
  for (int i = 0; i < n; ++i) perm[i] = (i * 7 + 3) % n;
  REQUIRE(std::gcd(7, n) == 1);
  od::Molecule shuffled = ref;
  for (int i = 0; i < n; ++i) shuffled.atoms[perm[i]] = ref.atoms[i];
  for (auto& b : shuffled.bonds) {
    b.a = perm[b.a];
    b.b = perm[b.b];
  }
  const auto map = od::referenceAtomMap(mol, heavy, shuffled);
  REQUIRE(map.size() == heavy.size());
  for (std::size_t k = 0; k < heavy.size(); ++k) CHECK(shuffled.atoms[map[k]].atomicNum == mol.atoms[heavy[k]].atomicNum);
  const auto rmsdOf = [&](const std::vector<int>& m) {
    std::vector<od::Vec3f> a, b;
    for (std::size_t k = 0; k < heavy.size(); ++k) {
      const auto& p1 = mol.atoms[heavy[k]].xyz;
      const auto& p2 = shuffled.atoms[m[k]].xyz;
      a.push_back({float(p1[0]), float(p1[1]), float(p1[2])});
      b.push_back({float(p2[0]), float(p2[1]), float(p2[2])});
    }
    std::vector<int> types;
    for (int h : heavy) types.push_back(mol.atoms[h].atomicNum);
    return od::SymmetryRmsd(mol, heavy, types, 256).rmsd(a, b);
  };
  CHECK(rmsdOf(map) < 1e-6);
  // Different molecule: no map.
  od::Molecule fewer = ref;
  fewer.atoms[heavy[0]].atomicNum = 9;
  CHECK_THROWS_AS(od::referenceAtomMap(mol, heavy, fewer), od::Error);

  // Kabsch: a rotated and translated copy is recovered exactly.
  std::vector<od::Vec3> from, to;
  const od::Quat q = od::normalized({0.8, 0.1, -0.5, 0.3});
  for (const auto& a : mol.atoms) {
    from.push_back(a.xyz);
    const od::Vec3 r = od::rotate(q, a.xyz);
    to.push_back({r[0] + 3, r[1] - 2, r[2] + 1});
  }
  double fit = 1;
  const od::Quat k = od::kabschQuat(from, to, fit);
  CHECK(fit < 1e-9);
  CHECK(od::rotationAngle(k, q) < 1e-6);
  std::filesystem::remove(two);
  std::filesystem::remove(bad);
}


TEST_CASE("protonation: PROPKA 3 port, Dimorphite-DL port") {
  // PROPKA 3 (Python, commit 58ceb7b) pKa values of 1L2Y.
  const od::RawStructure s = od::readMmcif(kData + "/1L2Y.cif", ccd(), true);
  std::vector<od::propka::InputAtom> in;
  for (const auto& a : s.atoms) {
    if (a.atomicNum == 1) continue;
    const auto& r = s.residues[a.residue].info;
    in.push_back({a.name, r.name, r.chain, r.seq, r.iCode, od::elementSymbol(a.atomicNum), a.hetatm, a.xyz});
  }
  const auto g = od::propka::calculate(in);
  const std::vector<std::pair<std::string, double>> ref{{"N+    1 A", 7.734941}, {"TYR   3 A", 10.404878}, {"LYS   8 A", 10.272701},
                                                        {"ASP   9 A", 4.100234}, {"ARG  16 A", 12.266749}, {"C-   20 A", 3.348626}};
  REQUIRE(g.size() == ref.size());
  for (std::size_t k = 0; k < g.size(); ++k) {
    CHECK(g[k].label == ref[k].first);
    CHECK(g[k].pka == doctest::Approx(ref[k].second).epsilon(1e-6));
  }

  // Dimorphite-DL 2.1.0 (Python) at one pH.
  auto dimo = [](const std::string& smi, double pH) {
    std::unique_ptr<RDKit::RWMol> m(RDKit::SmilesToMol(smi));
    std::vector<bool> none(m->getNumAtoms(), false);
    auto r = od::dimorphite::protonate(*m, pH, none, none);
    RDKit::MolOps::sanitizeMol(*r.mol);
    return RDKit::MolToSmiles(*r.mol);
  };
  auto canon = [](const std::string& smi) {
    std::unique_ptr<RDKit::RWMol> m(RDKit::SmilesToMol(smi));
    return RDKit::MolToSmiles(*m);
  };
  const std::vector<std::tuple<std::string, std::string, std::string>> cases{
      {"CC(=O)O", "CC(=O)[O-]", "CC(=O)[O-]"},
      {"CN", "C[NH3+]", "C[NH3+]"},
      {"c1c[nH]cn1", "c1c[nH]cn1", "c1c[nH+]c[nH]1"},
      {"Oc1ccccc1", "Oc1ccccc1", "Oc1ccccc1"},
      {"NC(=N)N", "NC(N)=[NH2+]", "NC(N)=[NH2+]"},
      {"OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1O", "Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])[O-])C(O)C1O",
       "Nc1ncnc2c1ncn2C1OC(COP(=O)([O-])O)C(O)C1O"},
      {"c1nn[nH]n1", "c1nn[n-]n1", "c1nn[nH]n1"},
      {"CC(=O)NC", "CNC(C)=O", "CNC(C)=O"}};
  for (const auto& [in7, at74, at4] : cases) {
    CHECK(dimo(in7, 7.4) == canon(at74));
    CHECK(dimo(in7, 4.0) == canon(at4));
  }
}

TEST_CASE("protonation: hydrogens of 1L2Y, overrides, pose records") {
  // The NMR hydrogens are removed and placed again: same count and names, the same positions
  // for hydrogens that cannot rotate.
  od::ProtonationOptions opt;
  od::ProtonationReport rep;
  const od::RawStructure ref = od::readMmcif(kData + "/1L2Y.cif", ccd());
  const od::RawStructure s = od::readMmcif(kData + "/1L2Y.cif", ccd(), false, &opt, &rep);
  CHECK(rep.removedH == 150);
  CHECK(rep.addedH == 150);
  REQUIRE(s.atoms.size() == ref.atoms.size());
  std::map<std::string, od::Vec3> refH;
  for (const auto& a : ref.atoms)
    if (a.atomicNum == 1) refH[std::format("{}:{}", ref.residues[a.residue].info.seq, a.name)] = a.xyz;
  int close = 0, total = 0;
  std::set<std::string> nterm;
  for (const auto& a : s.atoms) {
    if (a.atomicNum != 1) continue;
    const auto& r = s.residues[a.residue].info;
    if (r.seq == 1 && (a.name == "H1" || a.name == "H2" || a.name == "H3")) nterm.insert(a.name);
    const auto it = refH.find(std::format("{}:{}", r.seq, a.name));
    REQUIRE(it != refH.end());
    ++total;
    const double d = std::sqrt(std::pow(a.xyz[0] - it->second[0], 2) + std::pow(a.xyz[1] - it->second[1], 2) +
                               std::pow(a.xyz[2] - it->second[2], 2));
    if (d < 0.2) ++close;
  }
  CHECK(nterm.size() == 3);
  CHECK(close >= 110);  // methyl / NH3+ / OH hydrogens may rotate, Asn / Arg NH2 names follow the CCD
  std::map<std::string, std::string> state;
  for (const auto& x : rep.sites) state[x.atom] = x.state;
  CHECK(state["A:9:ASP:OD2"] == "ASP");
  CHECK(state["A:8:LYS:NZ"] == "LYS");
  CHECK(state["A:20:SER:OXT"] == "COO-");
  CHECK(state["A:1:ASN:N"] == "NH3+");
  int q = 0;
  for (const auto& a : s.atoms) q += a.formalCharge;
  CHECK(q == 1);

  // Overrides.
  od::ProtonationOptions o2;
  o2.overrides = {{"x", "A:9", "ASH", std::nullopt, "t:1"}, {"x", "A:1:N", "", 0, "t:2"}};
  od::ProtonationReport rep2;
  const od::RawStructure s2 = od::readMmcif(kData + "/1L2Y.cif", ccd(), false, &o2, &rep2);
  int q2 = 0;
  for (const auto& a : s2.atoms) q2 += a.formalCharge;
  CHECK(q2 == 1 + 1 - 1);  // ASH (+1) and a neutral N terminus (-1)
  o2.overrides = {{"x", "A:9", "HIP", std::nullopt, "t:1"}};
  CHECK_THROWS_AS(od::readMmcif(kData + "/1L2Y.cif", ccd(), false, &o2, &rep2), od::Error);
  o2.overrides = {{"x", "A:99:N", "", 1, "t:1"}};
  CHECK_THROWS_AS(od::readMmcif(kData + "/1L2Y.cif", ccd(), false, &o2, &rep2), od::Error);

  // Zinc finger: metal-bonded Cys and His.
  od::CcdLibrary lib;
  lib.addFile(kData + "/1ZNF_ccd.cif");
  od::ProtonationReport rep3;
  od::readMmcif(kData + "/1ZNF.cif", lib, false, &opt, &rep3);
  std::map<std::string, std::string> st3;
  for (const auto& x : rep3.sites) st3[x.atom] = x.method + " " + x.state;
  CHECK(st3["A:3:CYS:SG"] == "metal CYM");
  CHECK(st3["A:19:HIS:CG"] == "metal HID");

  // SDF: carboxylic acid of biotin deprotonated; pose records get the same hydrogens.
  const std::string sdf = kData + "/BTN_ideal.sdf";
  std::ifstream f(sdf);
  const std::string rec((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  const std::string two = (std::filesystem::temp_directory_path() / std::format("od_test_prot_{}.sdf", ::getpid())).string();
  std::ofstream(two) << rec << rec;
  od::ProtonationReport rep4;
  const od::PoseRecords p = od::loadPoseRecords(two, od::Role::MobM, ccd(), &opt, &rep4);
  REQUIRE(p.xyz.size() == 2);
  CHECK(p.xyz[0].size() == p.mol.atoms.size());
  CHECK(p.xyz[1].size() == p.mol.atoms.size());
  int qb = 0;
  for (const auto& a : p.mol.atoms) qb += a.formalCharge;
  CHECK(qb == -1);
  CHECK(rep4.removedH == 16);
  CHECK(rep4.addedH == 15);
  std::filesystem::remove(two);
}
