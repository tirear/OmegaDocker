// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Pair energy terms shared by the CPU and GPU back ends (spec 15.1). Function
// forms follow AutoDock-GPU 1.6 cuda/calcMergeEneGra.cu (LGPL-2.1-or-later);
// derivatives are those of the evaluated forms.
// Copyright (C) 2017 TU Darmstadt, Embedded Systems and Applications Group, Germany.
// Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

#pragma once
#include <cmath>

#if defined(__CUDACC__)
#define OD_HD __host__ __device__
#else
#define OD_HD
#endif

namespace od {

constexpr double kVdwCutoff = 8.0;     // Å, vdW and H-bond
constexpr double kElecCutoff = 20.48;  // Å, electrostatics and desolvation

// 12-6 / 12-10 term with AutoDock smoothing: within ±smooth/2 of reqm the
// distance is replaced by reqm, otherwise shifted by smooth/2 towards reqm.
template <class T>
OD_HD inline void ad4VdwHb(T r, T AC, T BD, int m, int n, T reqm, T halfSmooth, T& e, T& dedr) {
  T s = reqm;
  T ds = 0;
  if (r <= reqm - halfSmooth) {
    s = r + halfSmooth;
    ds = 1;
  } else if (r >= reqm + halfSmooth) {
    s = r - halfSmooth;
    ds = 1;
  }
  const T rm = std::pow(s, T(-m));
  const T rn = std::pow(s, T(-n));
  e = AC * rm - BD * rn;
  dedr = ds * (n * BD * rn - m * AC * rm) / s;
}

// Electrostatics with the fitted Mehler-Solmajer dielectric; coeffQ1Q2 =
// coeff_elec * q1 * q2. Distances below rMin are evaluated at rMin.
template <class T>
OD_HD inline void ad4Elec(T r, T coeffQ1Q2, T rMin, T& e, T& dedr) {
  T dr = 1;
  if (r < rMin) {
    r = rMin;
    dr = 0;
  }
  const T sh = r + T(1.26366);
  const T diel = T(1.10859) / (sh * sh) + T(0.010358);
  const T es = coeffQ1Q2 / r;
  e = es * diel;
  dedr = dr * (-(es / r) * diel - es * T(2.21718) / (sh * sh * sh));
}

// Desolvation with the fitted Gaussian (sigma 3.6 Å); sv = coeff_desolv *
// (S1*V2 + S2*V1).
template <class T>
OD_HD inline void ad4Desolv(T r, T sv, T rMin, T& e, T& dedr) {
  T dr = 1;
  if (r < rMin) {
    r = rMin;
    dr = 0;
  }
  const T r2 = r * r;
  const T num = T(12.96) - T(0.1063) * r2 * (T(1) - T(0.001947) * r2);
  const T den = T(12.96) + r2 * (T(0.4137) + r2 * (T(0.00357) + T(0.000112) * r2));
  const T dnum = -T(0.2126) * r + T(4) * T(0.1063) * T(0.001947) * r2 * r;
  const T dden = r * (T(2) * T(0.4137) + r2 * (T(4) * T(0.00357) + T(6) * T(0.000112) * r2));
  e = sv * num / den;
  dedr = dr * sv * (dnum * den - num * dden) / (den * den);
}

// Clash penalty (spec 7.9) for d < clashDistance.
template <class T>
OD_HD inline void clashPenalty(T r, T clashDistance, T base, T k, T& e, T& dedr) {
  const T d = clashDistance - r;
  e = base + k * d * d;
  dedr = -2 * k * d;
}

}  // namespace od
