// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <vector>

#include "energy.hpp"
#include "molecule.hpp"
#include "protonate.hpp"
#include "torsion_tree.hpp"

namespace od {

// Protonation states of standard amino acids recalculated at other conformations (protonation.recompute,
// spec 3.4.1): PROPKA 3 at the conformation, states by the pKa, and protonation variants of the molecules.

enum class SiteType : std::uint8_t { Asp, Glu, Cterm, Nterm, Cys, Tyr, Lys, Arg, His };

struct TitratableSite {
  int residue = -1;
  SiteType type{};
  int key = -1;        // atom of the PROPKA group: ASP CG, GLU CD, C- OXT, N+ N, CYS SG, TYR OH, LYS NZ, ARG CZ, HIS CG
  std::string state;   // state in the molecule: ASP / ASH, GLU / GLH, COO- / COOH, NH3+ / NH2, CYM / CYS, TYM / TYR,
                       // LYN / LYS, ARN / ARG, HID / HIE / HIP
  bool locked = false; // metal-bonded or set by [[protonation]]
  std::string label;   // "<chain>:<seq><iCode>:<residue name>"
};

// Sites of the standard amino acids of m; overrides: [[protonation]] entries of its input file.
std::vector<TitratableSite> titratableSites(const Molecule& m, const std::vector<ProtonationOverride>& overrides);
bool protonatedState(const std::string& state);

// A molecule at one conformation (coordinates per atom) for PROPKA.
struct PkaMolecule {
  const Molecule* mol = nullptr;
  std::vector<Vec3> xyz;
};
// pKa values per molecule: (key atom, PROPKA residue type) -> pKa. All molecules form one structure.
using PkaMap = std::map<std::pair<int, std::string>, double>;
std::vector<PkaMap> calculatePka(const std::vector<PkaMolecule>& mols);

struct StateChange {
  int site = -1;  // index into the sites of the molecule
  std::string state;
  double pka = 0;
};
// State changes at the pKa values (protonated when pKa >= pH; a site without a value keeps its state;
// His becoming neutral: HIE, then chosen by ProtonationVariants::chooseTautomers).
std::vector<StateChange> stateChanges(const std::vector<TitratableSite>& sites, const PkaMap& pka, double pH);
// Protonation free energy of the changes (kcal/mol): per site ln(10) RT (pH - pKa) when a proton is gained,
// ln(10) RT (pKa - pH) when one is lost.
double protonationFreeEnergy(const std::vector<TitratableSite>& sites, const std::vector<StateChange>& changes, double pH,
                             double temperature);
std::string describeChanges(const std::vector<TitratableSite>& sites, const std::vector<StateChange>& changes);

// Protonation variant of a molecule: hydrogens, formal charges and bond orders of the changed sites,
// Gasteiger charges, and the torsion tree of the base molecule (added hydrogens move with their parent).
struct MoleculeVariant {
  Molecule mol;
  TorsionTree tree;
  std::vector<int> toBase;  // variant atom -> base atom; -1 for an added hydrogen
};
MoleculeVariant makeVariant(const Molecule& base, const TorsionTree& tree, const std::vector<TitratableSite>& sites,
                            const std::vector<StateChange>& changes);

// Protonation variants of a group of molecules (the StaM group or one MobM). Variant 0 is the base
// system; the others are created when first needed. Thread safe.
class ProtonationVariants {
 public:
  ProtonationVariants(const EnergySystem& base, std::vector<EnergyPart> parts,
                      const std::vector<std::vector<ProtonationOverride>>& overrides, const EnergyParams& ep);
  bool active() const { return active_; }  // a site that may change
  std::size_t partCount() const { return sites_.size(); }
  const std::vector<TitratableSite>& sites(int part) const { return sites_[part]; }
  // Molecules of the base system at coordinates x of the base system (all coordinate indices).
  std::vector<PkaMolecule> pkaMolecules(std::span<const Vec3f> x) const;
  // Variant of the changes of each part; 0 without changes. Ids are given in call order.
  int variantOf(const std::vector<std::vector<StateChange>>& changes);
  // Energy system and parts of the changes (created when first needed, without an id).
  using TautomerEnergy = std::function<double(const EnergySystem&, const std::vector<EnergyPart>&)>;
  const EnergySystem& sysOf(const std::vector<std::vector<StateChange>>& changes);
  const std::vector<EnergyPart>& partsOf(const std::vector<std::vector<StateChange>>& changes);
  // Tautomers (HID / HIE) of the neutral His sites not locked, by the lowest energy(system of the
  // changes, its parts): from the current states (a His becoming neutral: HIE), one site at a time with
  // the others fixed, a change only when the energy is lower, until no site changes (at most
  // kTautomerPasses passes). Returns the number of sites changed.
  static constexpr int kTautomerPasses = 4;
  int chooseTautomers(std::vector<std::vector<StateChange>>& changes, const TautomerEnergy& energy);
  std::size_t count() const;
  const EnergySystem& sys(int id) const;
  std::vector<EnergyPart> parts(int id) const;
  // Changes (site, state; pKa 0) of each part.
  std::vector<std::vector<StateChange>> changes(int id) const;
  std::string describe(int id) const;

 private:
  struct Variant {
    std::vector<std::unique_ptr<MoleculeVariant>> mols;  // null: part unchanged
    std::vector<EnergyPart> parts;
    std::unique_ptr<EnergySystem> sys;
    std::vector<std::vector<StateChange>> changes;
  };
  const EnergySystem& base_;
  std::vector<EnergyPart> parts_;
  EnergyParams ep_;
  std::vector<std::vector<TitratableSite>> sites_;
  bool active_ = false;
  using Key = std::vector<std::vector<std::pair<int, std::string>>>;
  Key keyOf(const std::vector<std::vector<StateChange>>& changes) const;
  const Variant& entry(const Key& key);  // nullptr sys: no changes
  mutable std::mutex mutex_;
  std::map<Key, std::unique_ptr<Variant>> cache_;
  std::vector<const Variant*> variants_;  // index id - 1
  std::map<Key, int> index_;
};

}  // namespace od
