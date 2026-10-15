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
  auto o = parse({"--sta", "r.cif", "--mob", "a.sdf", "--mob=b.mol2", "--boxCenter", "1,-2.5,3",
                  "--boxSize", "20,20,30", "--outPrefix", "out", "--outCif", "true", "--outPdb=false",
                  "--mobmDockingConfSampleNum", "8"});
  CHECK(o.sta == "r.cif");
  REQUIRE(o.mob.size() == 2);
  CHECK(o.mob[1] == "b.mol2");
  CHECK((*o.boxCenter)[1] == doctest::Approx(-2.5));
  CHECK((*o.boxSize)[2] == doctest::Approx(30));
  CHECK(o.outCif);
  CHECK_FALSE(o.outPdb);
  CHECK(o.omegaTolerance == doctest::Approx(30));
  CHECK_FALSE(o.acceptNonPlanarOmega);
  CHECK(o.staFlexMargin == doctest::Approx(4));
  CHECK_NOTHROW(od::validateOptions(o));

  auto t = parse({"--omegaTolerance", "20", "--acceptNonPlanarOmega", "--staFlexMargin=6.5"});
  CHECK(t.omegaTolerance == doctest::Approx(20));
  CHECK(t.acceptNonPlanarOmega);
  CHECK(t.staFlexMargin == doctest::Approx(6.5));
}

TEST_CASE("cli: errors") {
  CHECK_THROWS_AS(parse({"--sta", "a", "--sta", "b"}), od::Error);
  CHECK_THROWS_AS(parse({"--unknown"}), od::Error);
  CHECK_THROWS_AS(parse({"--boxCenter", "1,2"}), od::Error);
  CHECK_THROWS_AS(parse({"--boxCenter", "1,2,x"}), od::Error);
  CHECK_THROWS_AS(parse({"--outPdb", "yes"}), od::Error);
  CHECK_THROWS_AS(parse({"--outPrefix"}), od::Error);
  CHECK_THROWS_AS(parse({"--help=1"}), od::Error);

  auto o = parse({"--sta", "r.cif", "--mob", "a.sdf", "--boxCenter", "0,0,0", "--boxSize", "10,10,10",
                  "--outPrefix", "o", "--outPdb", "false", "--mobmDockingConfSampleNum", "8"});
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);  // all output formats disabled
  o.outMol2 = true;
  CHECK_NOTHROW(od::validateOptions(o));
  o.omegaTolerance = 95;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.omegaTolerance = 30;
  o.staFlexMargin = -1;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.staFlexMargin = 4;
  o.boxSize = od::Vec3{10, 0, 10};
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
}

TEST_CASE("cli: 9a sample numbers and minimization") {
  auto o = parse({"--sta", "r.cif", "--mob", "a.sdf", "--boxCenter", "0,0,0", "--boxSize", "10,10,10", "--outPrefix", "o"});
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
  o.csgm1AngleStep = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.csgm1AngleStep = 1;
  o.confMaxDihedralStep = 0;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);

  auto t = parse({"--seed", "42", "--confMinMaxStep", "50", "--maxL1BatchNum=1000", "--ramaFloor", "0"});
  CHECK(*t.seed == 42);
  CHECK(t.confMinMaxStep == 50);
  CHECK(t.maxL1BatchNum == 1000);
  CHECK(t.ramaFloor == 0);
  CHECK_THROWS_AS(parse({"--seed", "-1"}), od::Error);
  CHECK_THROWS_AS(parse({"--confMinMaxStep", "1.5"}), od::Error);
}
