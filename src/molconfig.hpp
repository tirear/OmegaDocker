// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "molecule.hpp"
#include "topology.hpp"

namespace od {

// Molecule configuration file (spec section 5).

inline constexpr int kMolConfigFormatVersion = 1;

enum class TorsionMode : std::uint8_t { Fixed, Free, CisTrans };
enum class Csgm1 : std::uint8_t { A, B };

struct SsPoint {
  double phi = 0, psi = 0, prob = 0;
};

struct SsType {
  std::string name;
  bool builtin = false;
  std::vector<SsPoint> points;
};

struct ChainBreak {
  int bond = -1;       // index into Molecule::bonds
  std::string source;  // auto-loop | auto-crosslink | user
};

struct MolSettings {
  Topology topology;               // analyzed with the configured broken bonds
  std::vector<TorsionMode> mode;   // per Topology::torsions
  std::vector<std::string> ssType; // per Topology::units; protein backbones only
  Csgm1 csgm1 = Csgm1::A;
  std::vector<int> ringBreaks;     // indices into Molecule::bonds
  std::vector<ChainBreak> chainBreaks;
};

struct MolConfig {
  std::vector<SsType> ssTypes;
  std::vector<MolSettings> molecules;  // same order as the input molecules
};

// A molecule with its default topology (spec 4.1-4.7, DSSP assigned).
struct MolInput {
  const Molecule* mol = nullptr;
  const Topology* topo = nullptr;
};

// Built-in ssType names: R, O and the Ramachandran types of spec 7.5.
const std::vector<std::string>& builtinSsTypes();

std::string defaultMolConfigText(const std::vector<MolInput>& mols);

// Writes the default configuration to <outPrefix>.molconfig.default.toml, or to
// <outPrefix>.molconfig.default.<n>.toml when that file exists; returns the path.
std::string writeDefaultMolConfig(const std::string& outPrefix, const std::vector<MolInput>& mols);

// Reads and validates a configuration. Returns nothing and fills errors
// ("line <n>: <field>: <reason>") when the file is not valid. base holds the
// topology parameters of the run; broken bonds come from the file.
std::optional<MolConfig> readMolConfig(const std::string& path, const std::vector<MolInput>& mols,
                                       const TopologyParams& base, std::vector<std::string>& errors);

}  // namespace od
