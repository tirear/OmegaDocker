// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "grid_kernel.hpp"

namespace od {

// buildEnergyGrids, buildDistanceGrid and the delta blocks of StamField::delta on a GPU
// (spec 15.1.10): the same values as the CPU (gridPoint, distPoint).
class GpuGridBuilder {
 public:
  // slabValues: output values per kernel launch (0: default).
  explicit GpuGridBuilder(int device, std::size_t slabValues = 0);
  ~GpuGridBuilder();
  GpuGridBuilder(const GpuGridBuilder&) = delete;
  GpuGridBuilder& operator=(const GpuGridBuilder&) = delete;

  // v: host arrays (GridBuildData::view). typeMaps (v.typeMapCount), elec, desolv: host
  // arrays of n[0] × n[1] × n[2] values (GridGeometry::index order), type maps clamped at clamp.
  void build(const GridBuildView& v, double clamp, float* const* typeMaps, float* elec, float* desolv);
  // d: host array of n[0] × n[1] × n[2] values.
  void buildDistance(const DistBuildView& v, float* d);

  // delta: the reference maps (type maps, elec, desolv, distance; host arrays of `points` values)
  // are kept on the device under id.
  void setReference(std::uint64_t id, const float* const* maps, int mapCount, std::size_t points);
  std::uint64_t referenceId() const;
  // Blocks of the sample's maps (g, d) differing from the reference by more than threshold, in the
  // layout of StamField::delta(ref, refDist, block, data).
  void buildDelta(const GridBuildView& g, const DistBuildView& d, double clamp, double threshold,
                  std::vector<int>& block, std::vector<float>& data);
  // The same (spec 15.1.19); changed (optional): changedCount positions (3 values each,
  // changedPositions) as StamField::delta; data (optional): host copy of the blocks. Returns the
  // blocks in device memory.
  std::shared_ptr<const DeviceBlocks> buildDelta(const GridBuildView& g, const DistBuildView& d, double clamp,
                                                 double threshold, const double* changed, int changedCount,
                                                 std::vector<int>& block, std::vector<float>* data);
  // buildDelta in two parts (spec 15.1.20): submitDelta uploads the inputs (g, d and changed may
  // be released when it returns) and starts the computation on the builder's stream; finishDelta
  // waits for it and returns the blocks. One delta at a time.
  void submitDelta(const GridBuildView& g, const DistBuildView& d, double clamp, double threshold,
                   const double* changed, int changedCount);
  std::shared_ptr<const DeviceBlocks> finishDelta(std::vector<int>& block, std::vector<float>* data);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace od
