// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "grid.hpp"

namespace od {

// posm2_fft (spec 8.4.3): FFT correlation of the MobM probe atoms with the StaM maps.
// CPU reference implementation (pocketfft, single precision).

// Smallest n' >= n whose prime factors are 2, 3, 5 and 7.
int fftSize(int n);

// Spectra of the clamped maps of one StaM conformation: type maps, then elec, desolv.
struct FftStam {
  std::uint64_t id = 0;  // unique per buildFftStam call
  GridGeometry geom;
  std::array<int, 3> pad{};  // FFT points per axis (x, y, z), each >= geom.n
  double clamp = 10;
  std::vector<std::vector<std::complex<float>>> spectra;
  std::size_t realSize() const { return static_cast<std::size_t>(pad[0]) * pad[1] * pad[2]; }
  std::size_t complexSize() const { return static_cast<std::size_t>(pad[0] / 2 + 1) * pad[1] * pad[2]; }
  std::size_t index(int i, int j, int k) const { return (static_cast<std::size_t>(k) * pad[1] + j) * pad[0] + i; }
};
// Type maps are clamped to <= clamp, elec to [-clamp, clamp], desolv to <= clamp.
FftStam buildFftStam(const StamField& field, double clamp);

// Per-thread buffers of fftCorrelate.
struct FftWork {
  std::vector<float> energy;  // spread atoms of one channel; on return, energy per translation (index layout)
  std::vector<int> channel;   // spectra index of each probe spectrum (type maps with probe atoms, then elec, desolv)
  std::vector<std::vector<std::complex<float>>> spec;  // probe spectrum per channel
  std::vector<std::complex<float>> acc;                // sum over the channels
  std::vector<std::size_t> touched;
};

// Probe atom a sits at grid coordinates T + f[a] for the translation T (grid indices);
// map[a]: its type map; q[a]: its charge. On return, w.energy[index(T)] holds
// sum_a typeMap(T + f) + q elec(T + f) + |q| desolv(T + f) (trilinear, clamped maps) for every T
// with all T + f[a] inside the grid; other entries are not meaningful.
void fftCorrelate(const FftStam& s, std::span<const Vec3> f, std::span<const int> map, std::span<const double> q,
                  FftWork& w);
// fftCorrelate in two steps (spec 15.1.22): the probe spectra (w.channel, w.spec; only s.pad is
// used), then the energies against the spectra of s (w.energy). The probe spectra serve every
// StaM conformation of the same grid geometry; the result equals fftCorrelate bit for bit.
void fftProbeSpectra(const FftStam& s, std::span<const Vec3> f, std::span<const int> map, std::span<const double> q,
                     FftWork& w);
void fftCombine(const FftStam& s, FftWork& w);

// The top lowest (energy, index) of energy (index layout) over the translations lo..hi, ascending.
void fftSelect(const FftStam& s, const float* energy, const std::array<int, 3>& lo, const std::array<int, 3>& hi,
               std::size_t top, std::vector<std::pair<float, std::size_t>>& best);

}  // namespace od
