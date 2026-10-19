// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstddef>
#include <memory>
#include <string>

#include "dock_kernel.hpp"
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
// Free device memory (bytes).
std::size_t gpuFreeMem(int device);

// StaM fields (spec 8.1, 8.2) on the device: up to slots() fields stay resident (spec 15.1.16);
// fields with the same grids (FieldModel::base) share the reference maps.
class GpuField {
 public:
  explicit GpuField(int device = 0);
  ~GpuField();
  GpuField(const GpuField&) = delete;
  std::uint64_t id() const;  // FieldModel::id of the current field; 0: none
  // Makes fm the current field, uploading it unless resident (replacing the least recently used).
  void set(const FieldModel& fm);
  // Makes the resident field id current (most recently used); false when it is not resident.
  bool touch(std::uint64_t id);
  int slots() const;
  void setSlots(int n);  // resident fields (>= 1); keeps the most recently used
  int slotOf(std::uint64_t id) const;  // slot of a resident field; -1: none
  std::size_t deviceBytes() const;
  const FieldView& deviceView() const;  // current field, arrays in device memory
  const FieldView* deviceViews() const;  // device array: the view of each slot

 private:
  friend class GpuEnergy;
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

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
  // 9b pose energies (dock::evaluate, spec 15.1.4) against field with the probe atoms of
  // setProbes; angles: n × torsionCount; grad (optional): n × (6 + torsionCount). Each pose is
  // evaluated against the resident field of slot DockPose::field (spec 15.1.16).
  void setProbes(const std::vector<DockProbe>& probes, double clashDistance, double clashBase, double clashK);
  void evaluateDock(const GpuField& field, const float* angles, const DockPose* poses, int n, FlatTerms* terms,
                    DockTerms* dterms, double* grad);
  // Asynchronous form of evaluateDock (spec 15.1.8): kLanes independent evaluations may be in
  // flight. submitDock copies the inputs and returns; completeDock waits and stores the results
  // (grad: when requested at submitDock). Lanes have their own work memory.
  static constexpr int kLanes = 2;
  void submitDock(int lane, const GpuField& field, const float* angles, const DockPose* poses, int n, bool grad);
  void completeDock(int lane, FlatTerms* terms, DockTerms* dterms, double* grad);
  // Asynchronous form of evaluate without coordinates (9a, spec 15.1.11); fp32: coordinates and
  // pair terms in single precision (as 9b, spec 15.1.12).
  void submit(int lane, const double* angles, int n, bool grad, bool fp32 = false);
  void complete(int lane, FlatTerms* terms, double* grad);
  std::size_t bytesPerSampleDock() const;
  // 9b poses (with gradients, several fields) evaluated at the same time on the device.
  int residentSamplesDock(int probeCount) const;
  // Device work memory per sample of evaluate (bytes).
  std::size_t bytesPerSample(bool fp64, bool grad) const;
  // Frees the work memory (it is allocated again by the next evaluate).
  void releaseScratch();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
