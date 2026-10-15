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
                  "--dockingConfSampleNum", "8"});
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
                  "--outPrefix", "o", "--outPdb", "false", "--dockingConfSampleNum", "8"});
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
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);  // --dockingConfSampleNum missing
  o.dockingConfSampleNum = 10;
  auto n = od::confSampleNums(o);
  CHECK(n.minL2 == 20);
  CHECK(n.initL1 == 20);
  CHECK(n.elite == 3);
  CHECK_NOTHROW(od::validateOptions(o));
  o.dockingConfSampleNum = 3;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.dockingConfSampleNum = 10;
  o.minL2confSampleNum = 19;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.minL2confSampleNum = 30;
  o.initL1confSampleNum = 29;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.initL1confSampleNum.reset();
  CHECK(od::confSampleNums(o).initL1 == 30);
  o.eliteConfSampleNum = 6;
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
  o.eliteConfSampleNum = 5;
  CHECK_NOTHROW(od::validateOptions(o));
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
