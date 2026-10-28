// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>

namespace od {

// Parameters of the binding free energy od_0100 (spec 10.2):
//   dG = c0 + sum_t inter[t] E_inter,t + sum_t intra[t] dE_intra,t
//        + entropyStam d(-TS)_StaM + entropyMobm d(-TS)_MobM + tors N_tors + ring dE_ring,
// t = vdW, hbond, elec, desolv. Defaults are not calibrated. Read as the
// [od_0100] table of the parameter file (cli.hpp).
struct FeParams {
  static constexpr int kTerms = 4;  // vdw, hbond, elec, desolv
  double c0 = 0;
  double tors = 0;
  double ring = 0;  // ring conformer strain (relative MMFF94, spec 4.9)
  // Ranking coefficients fitted on Stage A Astex (spec 17.3.1); c0, tors not calibrated.
  std::array<double, kTerms> inter{1.718, 0.757, 0.863, 1.035};
  std::array<double, kTerms> intra{1.145, 0.428, 0.664, 0.709};
  double entropyStam = 0.105;
  double entropyMobm = 0.381;
  double temperature = 298.15;  // K, RT of -TS
  int unboundCandidateNum = 20;
  double shakingRadius = 1.4;  // Å
  double shakingStep = 0.2;    // Å
  double rrtScaling = 4.0;     // mol/kcal, w_k = exp(-(E_k - E_min) rrtScaling)
};

inline constexpr std::array<const char*, FeParams::kTerms> kFeTermNames{"vdw", "hbond", "elec", "desolv"};

}  // namespace od
