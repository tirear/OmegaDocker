// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <algorithm>
#include <vector>

#include "conf9a.hpp"
#include "energy_gpu.hpp"

namespace od {

// 9a energies on a GPU (spec 15.1.3). fp32 (gpu.conf9aPrecision = "float", spec 15.1.12): the
// evaluations with gradients during minimization in single precision; the final energies stay
// in double precision.
class GpuConfEvaluator : public ConfEvaluator {
 public:
  // batch: requested samples per call, limited to 90% of the free device memory.
  GpuConfEvaluator(const EnergyModel& model, int device, int batch, bool fp32 = false) : gpu_(model, device), fp32_(fp32) {
    const std::size_t per = gpu_.bytesPerSample(true, true);
    const std::size_t fit = static_cast<std::size_t>(0.9 * static_cast<double>(gpuDeviceInfo(device).freeMem)) / per;
    batch_ = static_cast<int>(std::clamp<std::size_t>(fit, 1, static_cast<std::size_t>(std::max(1, batch))));
  }
  int maxBatch() const override { return batch_; }
  void evaluate(const double* anglesDeg, int n, EnergyTerms* terms, double* gradRad) override {
    buf_.resize(n);
    gpu_.evaluate(anglesDeg, n, buf_.data(), gradRad, static_cast<double*>(nullptr));
    convert(n, terms);
  }
  int lanes() const override { return GpuEnergy::kLanes; }
  void submit(int lane, const double* anglesDeg, int n) override {
    size_[lane] = n;
    gpu_.submit(lane, anglesDeg, n, true, fp32_);
  }
  void complete(int lane, EnergyTerms* terms, double* gradRad) override {
    buf_.resize(size_[lane]);
    gpu_.complete(lane, buf_.data(), gradRad);
    convert(size_[lane], terms);
  }
  void idle() override { gpu_.releaseScratch(); }

 private:
  void convert(int n, EnergyTerms* terms) const {
    for (int i = 0; i < n; ++i) {
      const FlatTerms& f = buf_[i];
      EnergyTerms& e = terms[i];
      e = {};
      e.vdw = f.vdw;
      e.hbond = f.hbond;
      e.elec = f.elec;
      e.desolv = f.desolv;
      e.clash = f.clash;
      e.closure = f.closure;
      e.clashPairs = f.clashPairs;
      e.closureMaxDeviation = f.closureMaxDeviation;
    }
  }

  GpuEnergy gpu_;
  bool fp32_;
  int size_[GpuEnergy::kLanes] = {};
  int batch_ = 1;
  std::vector<FlatTerms> buf_;
};

}  // namespace od
