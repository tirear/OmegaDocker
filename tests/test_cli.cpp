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
                  "--boxSize", "20,20,30", "--outPrefix", "out", "--outCif", "true", "--outPdb=false"});
  CHECK(o.sta == "r.cif");
  REQUIRE(o.mob.size() == 2);
  CHECK(o.mob[1] == "b.mol2");
  CHECK((*o.boxCenter)[1] == doctest::Approx(-2.5));
  CHECK((*o.boxSize)[2] == doctest::Approx(30));
  CHECK(o.outCif);
  CHECK_FALSE(o.outPdb);
  CHECK_NOTHROW(od::validateOptions(o));
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
                  "--outPrefix", "o", "--outPdb", "false"});
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);  // all output formats disabled
  o.outMol2 = true;
  CHECK_NOTHROW(od::validateOptions(o));
  o.boxSize = od::Vec3{10, 0, 10};
  CHECK_THROWS_AS(od::validateOptions(o), od::Error);
}
