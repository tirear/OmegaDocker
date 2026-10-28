// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace od {

// One chemical component definition (wwPDB CCD).
struct CcdComp {
  struct Atom {
    std::string name;
    int atomicNum = 0;
    int charge = 0;
    std::array<double, 3> xyz{};  // model or ideal coordinates of the definition
    bool hasXyz = false;
    bool leaving = false;  // pdbx_leaving_atom_flag
  };
  struct Bond {
    int a = 0;
    int b = 0;
    int order = 1;  // 1, 2, 3 (Kekulé form of the definition)
    bool aromatic = false;
  };
  std::string id;
  std::string type;  // _chem_comp.type, e.g. "L-PEPTIDE LINKING"
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::unordered_map<std::string, int> atomIndex;

  int find(const std::string& atomName) const {
    auto it = atomIndex.find(atomName);
    return it == atomIndex.end() ? -1 : it->second;
  }
};

class CcdLibrary {
public:
  // Loads the CCD subset embedded in the executable.
  CcdLibrary();
  // Adds definitions from a user file (--ccd); they replace embedded ones with the same id.
  void addFile(const std::string& path);
  const CcdComp* find(const std::string& id) const;
  std::size_t size() const { return comps_.size(); }

private:
  void addDocument(const void* doc, const std::string& origin, bool fromUser);
  std::unordered_map<std::string, CcdComp> comps_;
};

}  // namespace od
