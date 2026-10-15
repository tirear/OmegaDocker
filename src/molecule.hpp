// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace od {

template <class T>
using Vec3T = std::array<T, 3>;
using Vec3 = Vec3T<double>;
using Vec3f = Vec3T<float>;

// Kekulé bond orders. Aromaticity is stored separately (Bond::aromatic).
enum class BondOrder : std::uint8_t { Single = 1, Double = 2, Triple = 3, Metal = 4 };

enum class Role : std::uint8_t { StaM, MobM };

struct ResidueInfo {
  std::string chain;
  int seq = 0;
  std::string iCode;
  std::string name;
};

struct Atom {
  std::string name;          // empty when the format has no atom names
  int atomicNum = 0;
  int formalCharge = 0;
  double partialCharge = 0.0;
  bool aromatic = false;
  Vec3 xyz{};
  int residue = -1;          // index into Molecule::residues, -1 if none
  int fileSerial = 0;        // 1-based atom index within the input record
};

struct Bond {
  int a = 0;
  int b = 0;
  BondOrder order = BondOrder::Single;
  bool aromatic = false;
};

struct Residue {
  ResidueInfo info;
  std::vector<int> atoms;
};

// One connected component of the bond graph (spec 3.2).
struct Molecule {
  std::string id;
  Role role = Role::MobM;
  std::string sourceFile;
  bool hasResidueInfo = false;
  bool chargesFromFile = false;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Residue> residues;
};

// A whole input record before splitting into molecules.
struct RawStructure {
  std::string sourceFile;
  bool hasResidueInfo = false;
  bool hasPartialCharges = false;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Residue> residues;
};

const char* elementSymbol(int atomicNum);
int atomicNumber(const std::string& symbol);  // case-insensitive; 0 if unknown
bool isMetal(int atomicNum);
// Metal coordination bond: bond order Metal, or a metal atom at either end.
bool isMetalBond(const Molecule& m, const Bond& b);

// Splits a structure into connected components, preserving atom order.
std::vector<Molecule> splitMolecules(const RawStructure& s);

// "stam" / "mobm"
const char* roleTag(Role r);

}  // namespace od
