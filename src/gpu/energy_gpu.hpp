// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstddef>
#include <memory>
#include <string>

#include "energy_kernel.hpp"

namespace od {

// CUDA back end (spec 15.1): batched evaluation of the energy of an
// EnergySystem (spec 7.9) in its flat form, one thread block per sample.

struct GpuDeviceInfo {
  std::string name;
  int major = 0, minor = 0;  // compute capability
  std::size_t totalMem = 0, freeMem = 0;
};

// Number of CUDA devices; 0 when there is no device or no driver.
int gpuDeviceCount();
GpuDeviceInfo gpuDeviceInfo(int device);

class GpuEnergy {
 public:
  // Uploads the model to the device. Throws od::Error on CUDA errors.
  GpuEnergy(const EnergyModel& model, int device = 0);
  ~GpuEnergy();
  GpuEnergy(const GpuEnergy&) = delete;
  GpuEnergy& operator=(const GpuEnergy&) = delete;

  // angles: n × torsionCount (degrees). terms: n. grad (optional): n ×
  // torsionCount, dE/dangle of the total energy in kcal/mol/rad. coords
  // (optional): n × 3 × coordCount.
  void evaluate(const double* angles, int n, FlatTerms* terms, double* grad, double* coords);
  // Single precision (spec 12.1): coordinates and per-item terms in float, sums in double.
  void evaluate(const float* angles, int n, FlatTerms* terms, double* grad, float* coords);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
