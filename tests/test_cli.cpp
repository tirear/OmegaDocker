// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <doctest.h>

#include <vector>

#include "cli.hpp"
#include "log.hpp"

namespace {

od::Options parse(std::vector<const char*> args) {
  args.insert(args.begin(), "omegadocker");
  return od::parseCommandLine(static_cast<int>(args.size()), args.data());
}

}  // namespace

TEST_CASE("cli: required options and values") {
  auto o = parse({"--stamFile", "r.cif", "--posm", "posm1", "--mobmFiles", "a.sdf", "--mobmFiles=b.mol2", "--boxCenter", "1,-2.5,3",
                  "--boxSize", "20,20,30", "--outPrefix", "out", "--outCif", "true", "--outPdb=false",
                  "--mobmDockingConfSampleNum", "8"});
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

  auto t = parse({"--omegaTolerance", "20", "--acceptNonPlanarOmega", "--stamFlexMargin=6.5"});
  CHECK(t.omegaTolerance == doctest::Approx(20));
  CHECK(t.acceptNonPlanarOmega);
  CHECK(t.stamFlexMargin == doctest::Approx(6.5));
}

TEST_CASE("cli: errors") {
  CHECK_THROWS_AS(parse({"--stamFile", "a", "--stamFile", "b"}), od::Error);
  CHECK_THROWS_AS(parse({"--unknown"}), od::Error);
  CHECK_THROWS_AS(parse({"--boxCenter", "1,2"}), od::Error);
  CHECK_THROWS_AS(parse({"--boxCenter", "1,2,x"}), od::Error);
  CHECK_THROWS_AS(parse({"--outPdb", "yes"}), od::Error);
  CHECK_THROWS_AS(parse({"--outPrefix"}), od::Error);
  CHECK_THROWS_AS(parse({"--help=1"}), od::Error);
  CHECK_THROWS_AS(parse({"--posm", "posm9"}), od::Error);
  CHECK_THROWS_AS(parse({"--posm", "posm2_fft"}), od::Error);
  CHECK_THROWS_AS(od::validateOptions(parse({"--stamFile", "r.cif", "--mobmFiles", "a.sdf", "--boxCenter", "0,0,0",
                                             "--boxSize", "10,10,10", "--outPrefix", "o"})),
                  od::Error);  // --posm missing

  auto o = parse({"--stamFile", "r.cif", "--posm", "posm1", "--mobmFiles", "a.sdf", "--boxCenter", "0,0,0", "--boxSize", "10,10,10",
                  "--outPrefix", "o", "--outPdb", "false", "--mobmDockingConfSampleNum", "8"});
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
  auto r = parse({"--maxRound", "3", "--topSampleNum", "20", "--diversityRmsd", "0"});
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
  CHECK(parse({"--sampleBufferNum", "100"}).sampleBufferNum == 100);
  CHECK(od::midTopSampleFactor(parse({"--midTopSampleFactor", "2"})) == 2);
  CHECK(od::midTopSampleFactor(parse({"--midKeepSkipped"})) == 2);
  CHECK(od::midTopSampleFactor(parse({"--midKeepSkipped", "--midTopSampleFactor", "1"})) == 1);
  o.midTopSampleFactor = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.midTopSampleFactor.reset();
  o.sampleBufferNum = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.sampleBufferNum = 100;
  CHECK(parse({"--symMaxPerm", "0"}).symMaxPerm == 0);
  o.symMaxPerm = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
}

TEST_CASE("cli: 9a sample numbers and minimization") {
  auto o = parse({"--stamFile", "r.cif", "--posm", "posm1", "--mobmFiles", "a.sdf", "--boxCenter", "0,0,0", "--boxSize", "10,10,10", "--outPrefix", "o"});
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

  auto s = parse({"--stamDockingConfSampleNum", "12", "--stamMinL2confSampleNum=30", "--stamInitL1confSampleNum", "40",
                  "--stamEliteConfSampleNum", "4", "--mobmDockingConfSampleNum", "8", "--mobmMinL2confSampleNum", "16",
                  "--mobmInitL1confSampleNum", "17", "--mobmEliteConfSampleNum", "2"});
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

  auto t = parse({"--seed", "42", "--confMinMaxStep", "50", "--mobmMaxL1BatchNum=1000", "--mobmRamaFloor", "0",
                  "--stamCsgmBuild", "grow", "--mobmCsgm1AngleStep=5", "--stamGrowMaxRetry", "3"});
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
  CHECK_THROWS_AS(parse({"--mobmCsgmBuild", "partial"}), od::Error);
  CHECK_THROWS_AS(parse({"--seed", "-1"}), od::Error);
  CHECK_THROWS_AS(parse({"--confMinMaxStep", "1.5"}), od::Error);
}
