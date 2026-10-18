// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <algorithm>
#include <vector>

#include "conf9a.hpp"
#include "energy_gpu.hpp"

namespace od {

// 9a energies on a GPU in double precision (spec 15.1.3).
class GpuConfEvaluator : public ConfEvaluator {
 public:
  // batch: requested samples per call, limited to 90% of the free device memory.
  GpuConfEvaluator(const EnergyModel& model, int device, int batch) : gpu_(model, device) {
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
    gpu_.submit(lane, anglesDeg, n, true);
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
  int size_[GpuEnergy::kLanes] = {};
  int batch_ = 1;
  std::vector<FlatTerms> buf_;
};

}  // namespace od
