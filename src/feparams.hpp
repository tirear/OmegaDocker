// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>

namespace od {

// Parameters of the binding free energy od_0100 (spec 10.2):
//   dG = c0 + sum_t inter[t] E_inter,t + sum_t intra[t] dE_intra,t
//        + entropyStam d(-TS)_StaM + entropyMobm d(-TS)_MobM + tors N_tors + ring dE_ring,
// t = vdW, hbond, elec, desolv. Read as the [od_0100] table of the parameter file (cli.hpp).
// Defaults (spec 17.3.1, 17.3.2): coefficient ratios fitted to pose ranking (Stage A Astex), then
// c0, tors and a common scale fitted to affinities (HiQBind subset); ring not calibrated (0).
struct FeParams {
  static constexpr int kTerms = 4;  // vdw, hbond, elec, desolv
  double c0 = -4.7559;
  double tors = 0.0075;
  double ring = 0;  // ring conformer strain (relative MMFF94, spec 4.9)
  std::array<double, kTerms> inter{0.4506, 0.1986, 0.2264, 0.2715};
  std::array<double, kTerms> intra{0.3003, 0.1123, 0.1742, 0.1860};
  double entropyStam = 0.0275;
  double entropyMobm = 0.0999;
  double temperature = 298.15;  // K, RT of -TS
  int unboundCandidateNum = 20;
  double shakingRadius = 1.4;  // Å
  double shakingStep = 0.2;    // Å
  double rrtScaling = 4.0;     // mol/kcal, w_k = exp(-(E_k - E_min) rrtScaling)
};

inline constexpr std::array<const char*, FeParams::kTerms> kFeTermNames{"vdw", "hbond", "elec", "desolv"};

}  // namespace od
