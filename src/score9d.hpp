// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <span>
#include <vector>

#include "dock9b.hpp"
#include "shaking.hpp"

namespace od {

// 9d recalculation and binding free energies (spec 10.1, 10.2).

// World coordinates of all coordinates of sys.coordinates(s.angles), single
// precision as in 9b: R(float q)(y - float center) + float t.
std::vector<Vec3f> dockWorld(std::span<const Vec3f> local, const DockSample& s);

// StaM–MobM energy by atom pairs (single precision terms, double sums):
// AD4 terms; clash = penalty of every MobM heavy atom whose nearest StaM heavy
// atom is closer than clashDistance (as the 9b distance grid); metal atoms are
// not checked on either side.
struct InterTerms {
  double vdw = 0, hbond = 0, elec = 0, desolv = 0, clash = 0;
  double inter() const { return vdw + hbond + elec + desolv; }
};
InterTerms exactInter(const ShakeAtoms& stam, const ShakeAtoms& mob, const EnergyParams& ep);

// sum_i ln f_i of the complex over the StaM atoms and over the MobM atoms. StaM
// atoms farther than kElecCutoff + radius from every MobM atom are not shaken:
// their values equal those of the StaM conformation alone (stamAlone, per StaM atom).
struct ComplexLnF {
  double stam = 0, mobm = 0;
};
ComplexLnF complexSumLnF(const ShakeAtoms& stam, std::span<const double> stamAlone, const ShakeAtoms& mob,
                         const EnergyParams& ep, const ShakeParams& p, const ShakeFn& fn = {});

constexpr double kAdgpuTorsionCoeff = 0.2983;  // AutoDock4.2 torsional free energy per torsion (kcal/mol)

}  // namespace od
