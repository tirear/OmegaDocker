// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <algorithm>
#include <memory>
#include <vector>

#include "dock9b.hpp"
#include "energy_gpu.hpp"

namespace od {

// 9b pose energies of one MobM on a GPU (spec 15.1.4); field is shared by the MobMs.
class GpuPoseEvaluator : public PoseEvaluator {
 public:
  // batch: requested poses per call, limited to 90% of the free device memory.
  GpuPoseEvaluator(const EnergyModel& model, std::shared_ptr<GpuField> field, int device, int batch)
      : gpu_(model, device), field_(std::move(field)), nt_(model.torsionCount), device_(device), want_(std::max(1, batch)) {
    batch_ = want_;
  }
  int maxBatch() const override { return batch_; }
  void bind(const MobmDocker& dk) override {
    if (field_->id() != dk.field().id()) field_->set(dk.field().model());
    const Dock9bParams& p = dk.params();
    gpu_.setProbes(dk.probes(), p.clashDistance, p.clashBase, p.clashK);
    if (sized_ == field_->id()) return;
    // Sized once per field, with this evaluator's work memory freed.
    gpu_.releaseScratch();
    const std::size_t per = gpu_.bytesPerSampleDock();
    const std::size_t fit = static_cast<std::size_t>(0.9 * static_cast<double>(gpuFreeMem(device_))) / per;
    batch_ = static_cast<int>(std::clamp<std::size_t>(fit, 1, static_cast<std::size_t>(want_)));
    sized_ = field_->id();
  }
  void evaluate(int n, const PoseRequest* req, MobmDocker::PoseEnergy* out, double* grad) override {
    const std::size_t nn = static_cast<std::size_t>(n), nt = static_cast<std::size_t>(nt_);
    angles_.resize(nn * nt);
    poses_.resize(nn);
    terms_.resize(nn);
    dterms_.resize(nn);
    for (std::size_t i = 0; i < nn; ++i) {
      const PoseRequest& r = req[i];
      for (std::size_t k = 0; k < nt; ++k) angles_[i * nt + k] = static_cast<float>(r.anglesDeg[k]);
      DockPose& ps = poses_[i];
      for (int c = 0; c < 3; ++c) {
        ps.t[c] = static_cast<float>(r.pose.t[c]);
        ps.center[c] = static_cast<float>(r.center[c]);
      }
      ps.q[0] = static_cast<float>(r.pose.q.w);
      ps.q[1] = static_cast<float>(r.pose.q.x);
      ps.q[2] = static_cast<float>(r.pose.q.y);
      ps.q[3] = static_cast<float>(r.pose.q.z);
    }
    gpu_.evaluateDock(*field_, angles_.data(), poses_.data(), n, terms_.data(), dterms_.data(), grad);
    for (std::size_t i = 0; i < nn; ++i) {
      const FlatTerms& f = terms_[i];
      MobmDocker::PoseEnergy& pe = out[i];
      pe = {};
      pe.intra.vdw = f.vdw;
      pe.intra.hbond = f.hbond;
      pe.intra.elec = f.elec;
      pe.intra.desolv = f.desolv;
      pe.intra.clash = f.clash;
      pe.intra.closure = f.closure;
      pe.intra.clashPairs = f.clashPairs;
      pe.intra.closureMaxDeviation = f.closureMaxDeviation;
      pe.inter = dterms_[i].inter;
      pe.interClash = dterms_[i].interClash;
      pe.intra.external = pe.inter + pe.interClash;
      pe.outside = dterms_[i].outside;
    }
  }
  void idle() override { gpu_.releaseScratch(); }

 private:
  GpuEnergy gpu_;
  std::shared_ptr<GpuField> field_;
  int nt_, device_, want_, batch_;
  std::uint64_t sized_ = 0;
  std::vector<float> angles_;
  std::vector<DockPose> poses_;
  std::vector<FlatTerms> terms_;
  std::vector<DockTerms> dterms_;
};

}  // namespace od
