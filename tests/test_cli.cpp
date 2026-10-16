// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <filesystem>
#include <format>
#include <fstream>
#include <vector>

#include <unistd.h>

#include "cli.hpp"
#include "log.hpp"

namespace {

od::Options parse(std::vector<const char*> args) {
  args.insert(args.begin(), "omegadocker");
  return od::parseCommandLine(static_cast<int>(args.size()), args.data());
}

}  // namespace

TEST_CASE("cli: required options and values") {
  auto o = parse({"--stamFile=r.cif", "--dock9b_posm=posm1", "--mobmFiles=a.sdf", "--mobmFiles=b.mol2", "--box_center=1,-2.5,3", "--box_size=20,20,30", "--outPrefix=out", "--output_cif=true", "--output_pdb=false", "--mobm_dockingConfSampleNum=8"});
  CHECK(o.stamFile == "r.cif");
  REQUIRE(o.mobmFiles.size() == 2);
  CHECK(o.mobmFiles[1] == "b.mol2");
  CHECK((*o.boxCenter)[1] == doctest::Approx(-2.5));
  CHECK((*o.boxSize)[2] == doctest::Approx(30));
  CHECK(o.outCif);
  CHECK_FALSE(o.outPdb);
  CHECK(o.omegaTolerance == doctest::Approx(30));
  CHECK_FALSE(o.acceptNonPlanarOmega);
  CHECK(o.stamFlexMargin == doctest::Approx(4));
  CHECK(o.posm == "posm1");
  CHECK(o.posStep == doctest::Approx(1));
  CHECK(o.posm1SampleNum == 1000);
  CHECK(o.dockMinMaxStep == 1000);
  CHECK(o.dockMinConvSteps == 10);
  CHECK(o.dockMinDeltaE == doctest::Approx(1e-3));
  CHECK(o.lbfgsMemory == 7);
  CHECK(o.maxTransStep == doctest::Approx(0.5));
  CHECK(o.maxAngleStep == doctest::Approx(10));
  CHECK(o.maxRound == 10);
  CHECK_FALSE(o.topSampleNum);
  CHECK(o.diversityRmsd == doctest::Approx(2.0));
  CHECK(o.sampleBufferNum == 65536);
  CHECK(od::midTopSampleFactor(o) == 1);
  CHECK_FALSE(o.midKeepSkipped);
  CHECK(o.symMaxPerm == 256);
  CHECK_NOTHROW(od::validateOptions(o));

  auto t = parse({"--topology_omegaTolerance=20", "--topology_acceptNonPlanarOmega=true", "--stam_flexMargin=6.5"});
  CHECK(t.omegaTolerance == doctest::Approx(20));
  CHECK(t.acceptNonPlanarOmega);
  CHECK(t.stamFlexMargin == doctest::Approx(6.5));
}

TEST_CASE("cli: errors") {
  CHECK_THROWS_AS(parse({"--stamFile=a", "--stamFile=b"}), od::Error);
  CHECK_THROWS_AS(parse({"--unknown"}), od::Error);
  CHECK_THROWS_AS(parse({"--box_center=1,2"}), od::Error);
  CHECK_THROWS_AS(parse({"--box_center=1,2,x"}), od::Error);
  CHECK_THROWS_AS(parse({"--output_pdb=maybe"}), od::Error);
  CHECK_THROWS_AS(parse({"--dock9c_maxRound", "3"}), od::Error);  // only --x=value
  CHECK_THROWS_AS(parse({"--dock9c_midKeepSkipped"}), od::Error);  // booleans need a value
  CHECK_THROWS_AS(parse({"--stamFile", "a"}), od::Error);
  CHECK_THROWS_AS(parse({"--stamCsgmBuild=grow"}), od::Error);  // old names are gone
  CHECK_THROWS_AS(parse({"--feParamFile=x"}), od::Error);
  CHECK_THROWS_AS(parse({"--outPrefix"}), od::Error);
  CHECK_THROWS_AS(parse({"--help=1"}), od::Error);
  CHECK_THROWS_AS(parse({"--dock9b_posm=posm9"}), od::Error);
  CHECK_THROWS_AS(parse({"--dock9b_posm=posm2_fft"}), od::Error);
  CHECK_THROWS_AS(od::validateOptions(parse({"--stamFile=r.cif", "--mobmFiles=a.sdf", "--box_center=0,0,0", "--box_size=10,10,10", "--outPrefix=o"})),
                  od::Error);  // --posm missing

  auto o = parse({"--stamFile=r.cif", "--dock9b_posm=posm1", "--mobmFiles=a.sdf", "--box_center=0,0,0", "--box_size=10,10,10", "--outPrefix=o", "--output_pdb=false", "--mobm_dockingConfSampleNum=8"});
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);  // all output formats disabled
  o.outMol2 = true;
  CHECK_NOTHROW(od::validateOptions(o));
  o.omegaTolerance = 95;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.omegaTolerance = 30;
  o.stamFlexMargin = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stamFlexMargin = 4;
  o.boxSize = od::Vec3{10, 0, 10};
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.boxSize = od::Vec3{10, 10, 10};
  auto r = parse({"--dock9c_maxRound=3", "--dock9c_topSampleNum=20", "--dock9c_diversityRmsd=0"});
  CHECK(r.maxRound == 3);
  CHECK(r.topSampleNum == 20);
  CHECK(r.diversityRmsd == 0);
  o.maxRound = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.maxRound = 1;
  o.topSampleNum = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.topSampleNum.reset();
  o.diversityRmsd = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.diversityRmsd = 2;
  CHECK(parse({"--dock9c_sampleBufferNum=100"}).sampleBufferNum == 100);
  CHECK(od::midTopSampleFactor(parse({"--dock9c_midTopSampleFactor=2"})) == 2);
  CHECK(od::midTopSampleFactor(parse({"--dock9c_midKeepSkipped=true"})) == 2);
  CHECK(od::midTopSampleFactor(parse({"--dock9c_midKeepSkipped=true", "--dock9c_midTopSampleFactor=1"})) == 1);
  o.midTopSampleFactor = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.midTopSampleFactor.reset();
  o.sampleBufferNum = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.sampleBufferNum = 100;
  CHECK(parse({"--dock9d_symMaxPerm=0"}).symMaxPerm == 0);
  o.symMaxPerm = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
}

TEST_CASE("cli: 9a sample numbers and minimization") {
  auto o = parse({"--stamFile=r.cif", "--dock9b_posm=posm1", "--mobmFiles=a.sdf", "--box_center=0,0,0", "--box_size=10,10,10", "--outPrefix=o"});
  CHECK_NOTHROW(od::validateOptions(o));  // required only with rotatable torsions
  CHECK_THROWS_AS(od::requireConfSampleNums(o, od::SampleSet::StaM), od::Error);
  CHECK_THROWS_AS(od::requireConfSampleNums(o, od::SampleSet::MobM), od::Error);
  o.mobm.docking = 10;
  auto n = od::confSampleNums(o, od::SampleSet::MobM);
  CHECK(n.minL2 == 20);
  CHECK(n.initL1 == 20);
  CHECK(n.elite == 3);
  CHECK_NOTHROW(od::validateOptions(o));
  CHECK_NOTHROW(od::requireConfSampleNums(o, od::SampleSet::MobM));
  CHECK_THROWS_AS(od::requireConfSampleNums(o, od::SampleSet::StaM), od::Error);
  o.stam.minL2 = 50;  // without --stamDockingConfSampleNum
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stam.docking = 20;
  CHECK_NOTHROW(od::validateOptions(o));
  CHECK(od::confSampleNums(o, od::SampleSet::StaM).initL1 == 50);
  CHECK(od::confSampleNums(o, od::SampleSet::StaM).elite == 5);
  CHECK(od::confSampleNums(o, od::SampleSet::MobM).minL2 == 20);
  o.stam.minL2 = 39;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stam.minL2.reset();
  o.mobm.docking = 3;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobm.docking = 10;
  o.mobm.minL2 = 19;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobm.minL2 = 30;
  o.mobm.initL1 = 29;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobm.initL1.reset();
  CHECK(od::confSampleNums(o, od::SampleSet::MobM).initL1 == 30);
  o.mobm.elite = 6;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobm.elite = 5;
  CHECK_NOTHROW(od::validateOptions(o));
  o.mobmCsgm.mutationRate = 101;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobmCsgm.mutationRate = 2;
  o.stamCsgm.tournamentRate = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stamCsgm.tournamentRate = 60;
  o.stamCsgm.absMaxDang = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stamCsgm.absMaxDang = 180;
  CHECK_NOTHROW(od::validateOptions(o));

  auto s = parse({"--stam_dockingConfSampleNum=12", "--stam_minL2confSampleNum=30", "--stam_initL1confSampleNum=40", "--stam_eliteConfSampleNum=4", "--mobm_dockingConfSampleNum=8", "--mobm_minL2confSampleNum=16", "--mobm_initL1confSampleNum=17", "--mobm_eliteConfSampleNum=2"});
  CHECK(*s.stam.docking == 12);
  CHECK(*s.stam.minL2 == 30);
  CHECK(*s.stam.initL1 == 40);
  CHECK(*s.stam.elite == 4);
  CHECK(*s.mobm.docking == 8);
  CHECK(*s.mobm.minL2 == 16);
  CHECK(*s.mobm.initL1 == 17);
  CHECK(*s.mobm.elite == 2);
  CHECK_THROWS_AS(parse({"--dockingConfSampleNum", "8"}), od::Error);
  o.stamCsgm.csgm1AngleStep = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.stamCsgm.csgm1AngleStep = 1;
  o.mobmCsgm.growMaxRetry = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.mobmCsgm.growMaxRetry = 0;
  o.confMaxDihedralStep = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);

  auto t = parse({"--seed=42", "--conf9a_minMaxStep=50", "--mobm_maxL1BatchNum=1000", "--mobm_ramaFloor=0", "--stam_csgmBuild=grow", "--mobm_csgm1AngleStep=5", "--stam_growMaxRetry=3"});
  CHECK(*t.seed == 42);
  CHECK(t.confMinMaxStep == 50);
  CHECK(t.mobmCsgm.maxL1BatchNum == 1000);
  CHECK(t.stamCsgm.maxL1BatchNum == 65536);
  CHECK_THROWS_AS(parse({"--maxL1BatchNum", "10"}), od::Error);
  CHECK(t.mobmCsgm.ramaFloor == 0);
  CHECK(t.stamCsgm.ramaFloor == doctest::Approx(1e-4));
  CHECK(t.stamCsgm.grow);
  CHECK_FALSE(t.mobmCsgm.grow);
  CHECK(t.mobmCsgm.csgm1AngleStep == 5);
  CHECK(t.stamCsgm.csgm1AngleStep == 1);
  CHECK(t.stamCsgm.growMaxRetry == 3);
  CHECK(t.mobmCsgm.growMaxRetry == 16);
  CHECK_THROWS_AS(parse({"--csgmBuild", "grow"}), od::Error);
  CHECK_THROWS_AS(parse({"--mobm_csgmBuild=partial"}), od::Error);
  auto ga = parse({"--stam_tournamentRate=70", "--mobm_crossoverRate=50", "--mobm_mutationRate=5", "--mobm_absMaxDang=45",
                   "--stam_csgmGA=csgm2"});
  CHECK(ga.stamCsgm.tournamentRate == 70);
  CHECK(ga.mobmCsgm.tournamentRate == 60);
  CHECK(ga.mobmCsgm.crossoverRate == 50);
  CHECK(ga.mobmCsgm.mutationRate == 5);
  CHECK(ga.mobmCsgm.absMaxDang == 45);
  CHECK_THROWS_AS(parse({"--mobm_csgmGA=csgm3"}), od::Error);  // not implemented yet
  CHECK_THROWS_AS(parse({"--mobm_csgmGA=csgm1"}), od::Error);
  CHECK_THROWS_AS(parse({"--seed=-1"}), od::Error);
  CHECK_THROWS_AS(parse({"--conf9a_minMaxStep=1.5"}), od::Error);
}

TEST_CASE("cli: booleans, arrays and od_0100 parameters") {
  for (const char* v : {"true", "1", "yes"}) CHECK(parse({std::format("--dock9c_midKeepSkipped={}", v).c_str()}).midKeepSkipped);
  for (const char* v : {"false", "0", "no"}) CHECK_FALSE(parse({std::format("--output_pdb={}", v).c_str()}).outPdb);
  CHECK_THROWS_AS(parse({"--output_pdb=True"}), od::Error);
  auto a = parse({"--ccd=a.cif,b.cif", "--topology_reachPerRes=4,8", "--dock9d_freeEnergy=adgpu_16,od_0100"});
  CHECK(a.ccd == std::vector<std::string>{"a.cif", "b.cif"});
  CHECK(a.reachProtein == 4);
  CHECK(a.reachNucleic == 8);
  CHECK(a.feAdgpu16);
  CHECK(a.feOd0100);
  auto n = parse({"--dock9d_freeEnergy=none"});
  CHECK_FALSE(n.feAdgpu16);
  CHECK_FALSE(n.feOd0100);
  CHECK_FALSE(parse({"--dock9d_freeEnergy="}).feOd0100);
  CHECK_THROWS_AS(parse({"--topology_reachPerRes=4"}), od::Error);
  CHECK_THROWS_AS(parse({"--dock9b_posStep=inf"}), od::Error);
  const od::FeParams d;
  CHECK(d.rrtScaling == 4.0);
  CHECK(d.inter[0] == 0.7564);
  auto f = parse({"--od_0100_c0=-1.5", "--od_0100_inter_elec=0.25", "--od_0100_shaking_rrtScaling=3",
                  "--od_0100_shaking_rrtScaling=5", "--od_0100_unboundCandidateNum=7", "--od_0100_entropy_mobm=0.7"});
  CHECK(f.fe.c0 == -1.5);
  CHECK(f.fe.inter == std::array<double, 4>{0.7564, 0.7564, 0.25, 0.7564});
  CHECK(f.fe.rrtScaling == 5);  // the last one wins
  CHECK(f.fe.unboundCandidateNum == 7);
  CHECK(f.fe.entropyMobm == 0.7);
  auto o = parse({"--stamFile=r.cif", "--dock9b_posm=posm1", "--mobmFiles=a.sdf", "--box_center=0,0,0", "--box_size=10,10,10",
                  "--outPrefix=o"});
  CHECK_NOTHROW(od::validateOptions(o));
  o.fe.rrtScaling = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.fe.rrtScaling = 0;
  o.fe.shakingStep = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
}

TEST_CASE("cli: --paramFile, command-line overrides and params.used.toml") {
  const auto path = (std::filesystem::temp_directory_path() / std::format("od_test_param_{}.toml", ::getpid())).string();
  auto write = [&](const std::string& s) { std::ofstream(path) << s; };
  const std::string file = "--paramFile=" + path;
  write("formatVersion = 1\nseed = 7\nccd = [\"x.cif\"]\n[box]\ncenter = [1, 2, 3.5]\nsize = [20, 20, 20]\n"
        "[stam]\ncsgmBuild = \"grow\"\nfixedResMovable = true\n[dock9b]\nposm = \"posm2\"\nposStep = 2\n"
        "[dock9c]\nmidKeepSkipped = true\n[od_0100.shaking]\nrrtScaling = 2\n[od_0100.inter]\nvdw = 0.5\n");
  auto o = parse({file.c_str(), "--stamFile=r.cif", "--mobmFiles=a.sdf", "--outPrefix=o"});
  CHECK(*o.seed == 7);
  CHECK(o.ccd == std::vector<std::string>{"x.cif"});
  CHECK((*o.boxCenter)[2] == 3.5);
  CHECK(o.stamCsgm.grow);
  CHECK(o.stamFixedResMovable);
  CHECK(o.posm == "posm2");
  CHECK(o.posStep == 2);
  CHECK(o.midKeepSkipped);
  CHECK(o.fe.rrtScaling == 2);
  CHECK(o.fe.inter[0] == 0.5);
  CHECK(o.fe.inter[1] == 0.7564);
  CHECK_NOTHROW(od::validateOptions(o));
  // The command line overrides the file, wherever --paramFile is given.
  auto c = parse({"--dock9c_midKeepSkipped=false", "--dock9b_posStep=1.5", file.c_str(), "--od_0100_shaking_rrtScaling=6"});
  CHECK_FALSE(c.midKeepSkipped);
  CHECK(c.posStep == 1.5);
  CHECK(c.fe.rrtScaling == 6);
  CHECK(c.stamCsgm.grow);
  // params.used.toml reads back to the same values.
  o.fe.tors = 0.2983;
  o.mobm.docking = 8;
  o.feAdgpu16 = true;
  const std::string used = od::paramsText(o);
  CHECK(used.find("# topSampleNum =") != std::string::npos);  // unset optional values are comments
  write(used);
  auto r = parse({file.c_str()});
  CHECK(*r.seed == 7);
  CHECK(r.ccd == o.ccd);
  CHECK(*r.boxCenter == *o.boxCenter);
  CHECK(r.stamCsgm.grow);
  CHECK(r.posStep == 2);
  CHECK(r.fe.tors == 0.2983);
  CHECK(r.fe.inter == o.fe.inter);
  CHECK(*r.mobm.docking == 8);
  CHECK(r.feAdgpu16);
  CHECK(r.feOd0100);
  CHECK_FALSE(r.topSampleNum);
  CHECK(od::paramsText(r) == used);
  // All problems of the file and the command line are listed together.
  write("formatVersion = 2\n[dock9b]\nposStep = \"x\"\nfoo = 1\n[dock9c]\nmaxRound = 2.5\nmidKeepSkipped = \"yes\"\n"
        "[box]\ncenter = [1, 2]\n[od_0100.entropy]\nstam = nan\n[other]\n");
  try {
    parse({file.c_str(), "--dock9b_rotStep=abc", "--output_pdb=maybe"});
    FAIL("no error");
  } catch (const od::Error& e) {
    const std::string m = e.what();
    CHECK(m.find("10 problem(s)") != std::string::npos);
    CHECK(m.find("line 1: formatVersion: must be 1") != std::string::npos);
    CHECK(m.find("line 3: dock9b.posStep: must be a number") != std::string::npos);
    CHECK(m.find("line 4: dock9b.foo: unknown key") != std::string::npos);
    CHECK(m.find("line 6: dock9c.maxRound: must be an integer") != std::string::npos);
    CHECK(m.find("line 7: dock9c.midKeepSkipped: must be true or false") != std::string::npos);
    CHECK(m.find("line 9: box.center: expected 3 numbers") != std::string::npos);
    CHECK(m.find("line 11: od_0100.entropy.stam: 'nan' is not finite") != std::string::npos);
    CHECK(m.find("other: unknown key") != std::string::npos);
    CHECK(m.find("command line: --dock9b_rotStep=abc: 'abc' is not a number") != std::string::npos);
    CHECK(m.find("command line: --output_pdb=maybe") != std::string::npos);
  }
  write("[dock9b\n");
  CHECK_THROWS_AS(parse({file.c_str()}), od::Error);
  std::filesystem::remove(path);
  CHECK_THROWS_AS(parse({file.c_str()}), od::Error);
  CHECK_THROWS_AS(parse({file.c_str(), file.c_str()}), od::Error);
}
