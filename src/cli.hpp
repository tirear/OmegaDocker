// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <optional>
#include <string>
#include <vector>

#include "molecule.hpp"

namespace od {

// Command-line parameters (spec section 14). Options are added as the
// corresponding stages are implemented.
struct Options {
  std::string sta;
  std::vector<std::string> mob;
  std::string molConfig;
  std::vector<std::string> ccd;
  std::optional<Vec3> boxCenter;
  std::optional<Vec3> boxSize;
  std::string outPrefix;
  bool outPdb = true;
  bool outCif = false;
  bool outSdf = false;
  bool outMol2 = false;

  bool showHelp = false;
  bool showVersion = false;
};

// Parses argv; throws od::Error on invalid input.
Options parseCommandLine(int argc, const char* const* argv);

// Checks required parameters and value ranges; throws od::Error.
void validateOptions(const Options& o);

std::string usageText();

}  // namespace od
