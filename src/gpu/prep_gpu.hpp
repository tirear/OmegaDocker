// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "energy_gpu.hpp"
#include "fft_dock.hpp"

namespace od {

// Device versions of the 9b quick filter (spec 8.5) and of the posm2_fft correlation
// (spec 8.4.3), shared by the MobMs (spec 15.1.5). Compiled without contracted
// multiply-adds: the filter follows MobmDocker::quickFilter operation by operation.
class GpuDockPrep {
 public:
  explicit GpuDockPrep(int device = 0);
  ~GpuDockPrep();
  GpuDockPrep(const GpuDockPrep&) = delete;
  GpuDockPrep& operator=(const GpuDockPrep&) = delete;

  struct FilterSetup {
    std::vector<float> y;   // probe atoms (torsion-tree frame), 3 per atom
    std::vector<int> heavy;  // per probe atom
    float center[3];
    float boxLo[3], boxHi[3];
    float clashDistance, nearDistance;
    int maxOverlap;
  };
  struct FilterPose {
    float t[3], q[4];  // q: w, x, y, z
    long long k;       // L1pos
  };
  struct FilterCounts {
    long long outsideBox = 0, overlap = 0, far = 0;
  };
  // Explicit poses; appends the L1pos of the passing ones to pass (any order).
  void filter(const GpuField& field, const FilterSetup& s, int n, const FilterPose* poses, std::vector<long long>& pass,
              FilterCounts& cnt);
  // posm2 L1pos [k0, k1), generated on the device as MobmDocker::l1Pose.
  struct Posm2Setup {
    std::vector<float> quats;              // 4 per orientation (float values of the Hopf grid)
    std::vector<int> range;                // lo[3], hi[3] per orientation
    std::vector<long long> offset;         // orientations + 1 prefix sums
    double base[3], step;                  // boxMin, posStep
  };
  void filterPosm2(const GpuField& field, const FilterSetup& s, const Posm2Setup& p2, long long k0, long long k1,
                   std::vector<long long>& pass, FilterCounts& cnt);

  // posm2_fft: uploads the spectra when s differs from the current ones.
  void setFft(const FftStam& s);
  // Orientations per fftTop call, limited by the free device memory.
  int fftBatch();
  // Probe spectra: direct (per-axis factors, all channels in one pass) or one FFT per channel;
  // Auto: direct up to 40 probe atoms per channel (spec 15.1.5).
  enum class Spectra { Auto, Direct, Fft };
  void setSpectra(Spectra s);
  // As PoseEvaluator::fftTop.
  void fftTop(int n, int probes, const double* f, const int* map, const double* q, const std::array<int, 3>* lo,
              const std::array<int, 3>* hi, int top, std::pair<float, std::size_t>* out);
  // Frees the work memory (spectra and plans are kept).
  void releaseWork();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
