// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <memory>
#include <span>
#include <vector>

#include "closure.hpp"
#include "conf9a.hpp"

namespace od {

// RCD attempts of the eligible links on a GPU (spec 15.1.18): one warp per attempt, the points of
// the link in shared memory. The generators of 9a run as coroutines (GenBatch); every flush runs
// the attempts of the waiting slots in one kernel per device, the slots divided among the devices
// in turn (spec 15.1.37).
class GpuRcd : public GenBatch, public RcdRunner {
 public:
  GpuRcd(const LoopCloser& closer, std::span<const int> devices, int slots);
  ~GpuRcd() override;
  GpuRcd(const GpuRcd&) = delete;
  GpuRcd& operator=(const GpuRcd&) = delete;
  // At least one link can be closed on the GPU.
  bool active() const;
  int slots() const override;
  void flush(std::span<const int> slots) override;
  void run(int slot, RcdJob& job) override;
  // Statistics: kernels, attempts, RCD updates and kernel time (s).
  long long kernels() const;
  long long attempts() const;
  long long updates() const;
  double kernelSec() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
