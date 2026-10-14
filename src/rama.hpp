// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace od {

// Ramachandran classes of the Top8000 rama8000 grids (spec 7.5).
enum class RamaClass : std::uint8_t { General, Gly, TransPro, CisPro, PrePro, IleVal };

inline constexpr int kRamaClassCount = 6;
inline constexpr int kRamaBins = 180;  // per angle, 2 degrees each, from -180

// Percentile contour values of one rama8000 file; cells not listed are 0.
struct RamaGrid {
  std::vector<float> value;  // [phiBin * kRamaBins + psiBin]
  int listed = 0;            // cells listed in the file

  static int bin(double angleDeg);  // periodic
  float at(double phiDeg, double psiDeg) const { return value[bin(phiDeg) * kRamaBins + bin(psiDeg)]; }
};

// Parsed on first use from the data embedded in the executable.
const RamaGrid& ramaGrid(RamaClass c);

const char* ramaTypeName(RamaClass c);    // "RAMA_GENERAL", ...
const char* ramaClassLabel(RamaClass c);  // "general", "Gly", "trans-Pro", ...
std::optional<RamaClass> ramaClassFromTypeName(const std::string& name);

// Class of a residue from its name, the name of the next residue in the chain
// ("" when none) and whether the omega before it is cis.
RamaClass ramaClassOf(const std::string& resName, const std::string& nextResName, bool cisBefore);

}  // namespace od
