// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <memory>
#include <span>
#include <vector>

#include "shaking.hpp"

namespace od {

// shakeLnF on a GPU (spec 10.2, 15.1.17): one block per shaken atom, the threads over the
// sphere points; per-pair terms in single precision as the CPU, the energy of a point summed in
// single precision per tile of partners and in double over the tiles.
class GpuShaker {
 public:
  GpuShaker(const EnergyParams& ep, const ShakeParams& p, int device);
  ~GpuShaker();
  GpuShaker(const GpuShaker&) = delete;
  GpuShaker& operator=(const GpuShaker&) = delete;
  std::vector<double> lnF(const ShakeAtoms& atoms, std::span<const char> shake = {});

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
