// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <string>

namespace od {

// Coefficients of the binding free energy od_0100 (spec 10.2):
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
};

inline constexpr std::array<const char*, FeParams::kTerms> kFeTermNames{"vdw", "hbond", "elec", "desolv"};

// Reads a free energy parameter file (TOML); keys not given keep their defaults.
// Throws od::Error listing every problem with its line.
FeParams readFeParams(const std::string& path);

// The file with the values of fp.
std::string feParamsText(const FeParams& fp);

}  // namespace od
