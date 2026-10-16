// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace od {

// Parameters of the binding free energy od_0100 (spec 10.2):
//   dG = c0 + sum_t inter[t] E_inter,t + sum_t intra[t] dE_intra,t
//        + entropyStam d(-TS)_StaM + entropyMobm d(-TS)_MobM + tors N_tors,
// t = vdW, hbond, elec, desolv. Defaults are not calibrated.
struct FeParams {
  static constexpr int kTerms = 4;  // vdw, hbond, elec, desolv
  double c0 = 0;
  double tors = 0;
  std::array<double, kTerms> inter{0.7564, 0.7564, 0.7564, 0.7564};
  std::array<double, kTerms> intra{0.7564, 0.7564, 0.7564, 0.7564};
  double entropyStam = 0.9;
  double entropyMobm = 0.9;
  double temperature = 298.15;  // K, RT of -TS
  int unboundCandidateNum = 20;
  double shakingRadius = 1.4;  // Å
  double shakingStep = 0.2;    // Å
  double rrtScaling = 4.0;     // mol/kcal, w_k = exp(-(E_k - E_min) rrtScaling)
};

inline constexpr std::array<const char*, FeParams::kTerms> kFeTermNames{"vdw", "hbond", "elec", "desolv"};

// Values given on the command line: (key, value text) in the order given,
// key as in the file (e.g. "od_0100.inter.vdw").
using FeOverrides = std::vector<std::pair<std::string, std::string>>;

// Built-in defaults, then the file (when path is not empty; keys not given keep
// their values), then the overrides. Throws od::Error listing every problem
// (file problems with their lines, overrides with their option).
FeParams loadFeParams(const std::string& path, const FeOverrides& overrides);

// The file with the values of fp.
std::string feParamsText(const FeParams& fp);

}  // namespace od
