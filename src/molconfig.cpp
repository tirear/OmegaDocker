// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "molconfig.hpp"

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <set>

#include <toml.hpp>

#include "log.hpp"
#include "rama.hpp"
#include "version.hpp"

namespace od {
namespace {

constexpr BackboneAngle kProteinAngles[] = {BackboneAngle::Phi, BackboneAngle::Psi, BackboneAngle::Omega};
constexpr BackboneAngle kNucleicAngles[] = {BackboneAngle::Alpha, BackboneAngle::Beta,    BackboneAngle::Gamma,
                                            BackboneAngle::Delta, BackboneAngle::Epsilon, BackboneAngle::Zeta};

std::string quote(const std::string& s) {
  std::string o = "\"";
  for (char c : s) {
    switch (c) {
      case '"': o += "\\\""; break;
      case '\\': o += "\\\\"; break;
      case '\n': o += "\\n"; break;
      case '\t': o += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20)
          o += std::format("\\u{:04X}", static_cast<int>(c));
        else
          o += c;
    }
  }
  return o + "\"";
}

std::string residueKey(const ResidueInfo& r) { return std::format("{}:{}{}", r.chain, r.seq, r.iCode); }

std::string chainList(const Molecule& m) {
  std::vector<std::string> chains;
  for (const Residue& r : m.residues)
    if (std::find(chains.begin(), chains.end(), r.info.chain) == chains.end()) chains.push_back(r.info.chain);
  std::string k;
  for (const auto& c : chains) k += (k.empty() ? "" : "+") + c;
  return k;
}

// Atom references of spec 5.2: "<name>" within the residue of the entry,
// "<chain>:<seq><iCode>:<name>" elsewhere, "#<serial>" when the atom has no
// residue or its name is empty or not unique in the residue.
class AtomRefs {
 public:
  explicit AtomRefs(const Molecule& m) : m_(m), unique_(m.atoms.size(), false) {
    for (std::size_t i = 0; i < m.atoms.size(); ++i) serial_[m.atoms[i].fileSerial] = static_cast<int>(i);
    for (std::size_t r = 0; r < m.residues.size(); ++r) {
      const std::string key = residueKey(m.residues[r].info);
      residueOf_[key] = residueOf_.contains(key) ? -2 : static_cast<int>(r);
      std::map<std::string, int> count;
      for (int a : m.residues[r].atoms) ++count[m.atoms[a].name];
      for (int a : m.residues[r].atoms) unique_[a] = !m.atoms[a].name.empty() && count[m.atoms[a].name] == 1;
    }
  }

  std::string ref(int a, int context) const {
    const Atom& x = m_.atoms[a];
    if (x.residue >= 0 && unique_[a] && residueOf_.at(residueKey(m_.residues[x.residue].info)) == x.residue)
      return x.residue == context ? x.name : residueKey(m_.residues[x.residue].info) + ":" + x.name;
    return std::format("#{}", x.fileSerial);
  }

  // Atom index, or -1 with err set.
  int parse(const std::string& s, int context, std::string& err) const {
    if (!s.empty() && s[0] == '#') {
      int v = 0;
      auto [p, ec] = std::from_chars(s.data() + 1, s.data() + s.size(), v);
      auto it = serial_.find(v);
      if (ec != std::errc() || p != s.data() + s.size() || it == serial_.end()) {
        err = std::format("atom '{}' not found", s);
        return -1;
      }
      return it->second;
    }
    int residue = context;
    std::string name = s;
    if (const auto c1 = s.find(':'); c1 != std::string::npos) {
      const auto c2 = s.find(':', c1 + 1);
      if (c2 == std::string::npos || s.find(':', c2 + 1) != std::string::npos) {
        err = std::format("atom reference '{}' must be <name>, <chain>:<seq><iCode>:<name> or #<serial>", s);
        return -1;
      }
      auto it = residueOf_.find(s.substr(0, c2));
      if (it == residueOf_.end() || it->second < 0) {
        err = std::format("residue '{}' not found{}", s.substr(0, c2),
                          it != residueOf_.end() ? " (not unique; use #<serial>)" : "");
        return -1;
      }
      residue = it->second;
      name = s.substr(c2 + 1);
    } else if (residue < 0) {
      err = std::format("atom reference '{}' must be <chain>:<seq><iCode>:<name> or #<serial> here", s);
      return -1;
    }
    for (int a : m_.residues[residue].atoms)
      if (m_.atoms[a].name == name) {
        if (!unique_[a]) {
          err = std::format("atom name '{}' is not unique in residue {}; use #<serial>", name,
                            residueKey(m_.residues[residue].info));
          return -1;
        }
        return a;
      }
    err = std::format("atom '{}' not found in residue {}", name, residueKey(m_.residues[residue].info));
    return -1;
  }

 private:
  const Molecule& m_;
  std::vector<bool> unique_;
  std::map<int, int> serial_;
  std::map<std::string, int> residueOf_;  // -2 when the key is not unique
};

// Order in which backbone units and torsions appear in the file.
struct Layout {
  std::vector<int> units;                               // unit indices along each chain
  std::vector<std::map<BackboneAngle, int>> backbone;   // per position: angle -> torsion index
  std::vector<std::vector<int>> torsions;               // per position: other torsions of the residue
  std::vector<int> molTorsions;                         // torsions outside backbone residues
};

Layout layoutOf(const Topology& t) {
  Layout l;
  std::vector<bool> seen(t.units.size(), false);
  auto walk = [&](int s) {
    for (int u = s; u >= 0 && !seen[u]; u = t.units[u].next) {
      seen[u] = true;
      l.units.push_back(u);
    }
  };
  for (std::size_t i = 0; i < t.units.size(); ++i)
    if (t.units[i].prev < 0) walk(static_cast<int>(i));
  for (std::size_t i = 0; i < t.units.size(); ++i)
    if (!seen[i]) walk(static_cast<int>(i));
  std::vector<int> posOfUnit(t.units.size(), -1);
  std::map<int, int> posOfResidue;
  for (std::size_t k = 0; k < l.units.size(); ++k) {
    posOfUnit[l.units[k]] = static_cast<int>(k);
    if (t.units[l.units[k]].residue >= 0) posOfResidue[t.units[l.units[k]].residue] = static_cast<int>(k);
  }
  l.backbone.resize(l.units.size());
  l.torsions.resize(l.units.size());
  for (std::size_t i = 0; i < t.torsions.size(); ++i) {
    const Torsion& tor = t.torsions[i];
    if (tor.backbone != BackboneAngle::None) {
      l.backbone[posOfUnit[tor.unit]][tor.backbone] = static_cast<int>(i);
    } else if (auto it = posOfResidue.find(tor.residue); tor.residue >= 0 && it != posOfResidue.end()) {
      l.torsions[it->second].push_back(static_cast<int>(i));
    } else {
      l.molTorsions.push_back(static_cast<int>(i));
    }
  }
  return l;
}

// Chain-end residues (without phi or psi) always use R (spec 7.5).
bool chainEnd(const std::map<BackboneAngle, int>& backbone) {
  return !backbone.contains(BackboneAngle::Phi) || !backbone.contains(BackboneAngle::Psi);
}

std::string defaultSsType(const Topology& t, int unit, bool end) {
  const int r = t.units[unit].residue;
  if (end || r < 0 || r >= static_cast<int>(t.dsspSS.size())) return "R";
  const char c = t.dsspSS[r];
  return c == 'H' || c == 'G' || c == 'I' ? "O" : "R";
}

// Empty when the torsion may be rotatable, otherwise the reason.
std::string notRotatableReason(const Torsion& tor) {
  if (tor.broken) return "the bond is broken (ringBreak / chainBreak)";
  if (tor.category == TorsionCategory::Single) return "";
  return std::format("category '{}' is not rotatable{}", torsionCategoryName(tor.category),
                     tor.category == TorsionCategory::Ring ? " (a ringBreak in the ring changes this)" : "");
}

// ---------- writing ----------

void writeTorsion(std::string& o, const std::string& table, const Torsion& tor, const AtomRefs& refs, int context,
                  const std::string& indent) {
  o += std::format("\n{}[[{}]]\n", indent, table);
  o += std::format("{}atoms     = [{}, {}, {}, {}]\n", indent, quote(refs.ref(tor.atoms[0], context)),
                   quote(refs.ref(tor.atoms[1], context)), quote(refs.ref(tor.atoms[2], context)),
                   quote(refs.ref(tor.atoms[3], context)));
  o += std::format("{}rotatable = {}\n", indent, tor.rotatable ? "true" : "false");
  o += std::format("{}category  = {}\n", indent, quote(torsionCategoryName(tor.category)));
}

void writeMolecule(std::string& o, const Molecule& m, const Topology& t) {
  const AtomRefs refs(m);
  const Layout l = layoutOf(t);
  const bool protein = t.backboneType == BackboneType::Protein;
  o += "\n[[molecule]]\n";
  o += std::format("id           = {}\n", quote(m.id));
  o += std::format("role         = {}\n", quote(m.role == Role::StaM ? "StaM" : "MobM"));
  o += std::format("file         = {}\n", quote(m.sourceFile));
  if (m.hasResidueInfo) o += std::format("chain        = {}\n", quote(chainList(m)));
  o += std::format("backboneType = {}\n", quote(backboneTypeName(t.backboneType)));
  if (protein) o += "csgm1        = \"1a\"\n";
  if (hasMetalRing(m)) o += "metalBond    = \"break\"           # \"break\" | \"keep\" (see the top of this file)\n";

  auto writeBreak = [&](int bi, const char* source) {
    o += "\n  [[molecule.chainBreak]]\n";
    o += std::format("  bond   = [{}, {}]\n", quote(refs.ref(m.bonds[bi].a, -1)), quote(refs.ref(m.bonds[bi].b, -1)));
    o += std::format("  source = \"{}\"\n", source);
  };
  for (int bi : t.crosslinkBreaks) writeBreak(bi, "auto-crosslink");
  if (hasMetalRing(m))
    for (int bi : t.metalBreaks) writeBreak(bi, "auto-metal");

  for (std::size_t k = 0; k < l.units.size(); ++k) {
    const int u = l.units[k];
    const int r = t.units[u].residue;
    o += "\n  [[molecule.residue]]\n";
    if (r >= 0) {
      const ResidueInfo& ri = m.residues[r].info;
      o += std::format("  chain   = {}\n  resName = {}\n  resSeq  = {}\n  iCode   = {}\n", quote(ri.chain), quote(ri.name),
                       ri.seq, quote(ri.iCode));
    } else {
      o += std::format("  unit    = {}\n", k + 1);
    }
    if (protein) {
      if (r >= 0 && r < static_cast<int>(t.dsspSS.size())) o += std::format("  dsspSS  = {}\n", quote(std::string(1, t.dsspSS[r])));
      o += std::format("  ssType  = {}\n", quote(defaultSsType(t, u, chainEnd(l.backbone[k]))));
    }
    for (const auto& [a, ti] : l.backbone[k]) {
      const char* name = backboneAngleName(a);
      if (a == BackboneAngle::Omega)
        o += std::format("  {:<7} = \"fixed\"\n", name);
      else
        o += std::format("  {:<7} = {}\n", name, t.torsions[ti].rotatable ? "true" : "false");
    }
    for (int ti : l.torsions[k]) writeTorsion(o, "molecule.residue.torsion", t.torsions[ti], refs, r, "    ");
  }
  for (int ti : l.molTorsions) writeTorsion(o, "molecule.torsion", t.torsions[ti], refs, -1, "  ");
}

// ---------- reading ----------

// Collects errors; flush() appends them to the output sorted by line, those
// without a line last.
class Errors {
 public:
  void add(int line, const std::string& field, const std::string& msg) {
    list_.push_back({line > 0 ? line : std::numeric_limits<int>::max(),
                     line > 0 ? std::format("line {}: {}: {}", line, field, msg) : std::format("{}: {}", field, msg)});
  }
  std::size_t count() const { return list_.size(); }
  void flush(std::vector<std::string>& out) {
    std::stable_sort(list_.begin(), list_.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& x : list_) out.push_back(std::move(x.second));
  }

 private:
  std::vector<std::pair<int, std::string>> list_;
};

int lineOf(const toml::node& n) { return static_cast<int>(n.source().begin.line); }

void checkKeys(const toml::table& t, const std::set<std::string>& allowed, const std::string& where, Errors& e) {
  for (auto&& [k, v] : t)
    if (!allowed.contains(std::string(k.str()))) e.add(lineOf(v), where + "." + std::string(k.str()), "unknown key");
}

template <class T>
bool hasType(const toml::node& n) {
  if constexpr (std::is_same_v<T, std::string>) return n.is_string();
  else if constexpr (std::is_same_v<T, bool>) return n.is_boolean();
  else if constexpr (std::is_same_v<T, std::int64_t>) return n.is_integer();
  else return n.is_number();
}

template <class T>
const char* typeName() {
  if constexpr (std::is_same_v<T, std::string>) return "a string";
  else if constexpr (std::is_same_v<T, bool>) return "true or false";
  else if constexpr (std::is_same_v<T, std::int64_t>) return "an integer";
  else return "a number";
}

template <class T>
std::optional<T> get(const toml::table& t, const std::string& key, const std::string& where, Errors& e,
                     bool required) {
  const toml::node* n = t.get(key);
  if (!n) {
    if (required) e.add(lineOf(t), where + "." + key, "missing");
    return std::nullopt;
  }
  if (!hasType<T>(*n)) {
    e.add(lineOf(*n), where + "." + key, std::format("must be {}", typeName<T>()));
    return std::nullopt;
  }
  return n->value<T>();
}

// Tables of an array of tables; nothing when absent.
std::vector<const toml::table*> tablesOf(const toml::table& t, const std::string& key, const std::string& where,
                                         Errors& e) {
  std::vector<const toml::table*> out;
  const toml::node* n = t.get(key);
  if (!n) return out;
  const toml::array* arr = n->as_array();
  bool ok = arr != nullptr;
  if (ok)
    for (const toml::node& x : *arr) ok &= x.is_table();
  if (!ok) {
    e.add(lineOf(*n), where + "." + key, "must be an array of tables ([[" + where + "." + key + "]])");
    return out;
  }
  for (const toml::node& x : *arr) out.push_back(x.as_table());
  return out;
}

// Atom indices of a string array of the given size; nothing on error.
std::optional<std::vector<int>> atomList(const toml::table& t, const std::string& key, std::size_t size,
                                         const AtomRefs& refs, int context, const std::string& where, Errors& e) {
  const toml::node* n = t.get(key);
  if (!n) {
    e.add(lineOf(t), where + "." + key, "missing");
    return std::nullopt;
  }
  const toml::array* arr = n->as_array();
  if (!arr || arr->size() != size || !std::all_of(arr->begin(), arr->end(), [](const toml::node& x) { return x.is_string(); })) {
    e.add(lineOf(*n), where + "." + key, std::format("must be an array of {} atom references", size));
    return std::nullopt;
  }
  std::vector<int> out;
  bool ok = true;
  for (const toml::node& x : *arr) {
    std::string err;
    const int a = refs.parse(*x.value<std::string>(), context, err);
    if (a < 0) {
      e.add(lineOf(x), where + "." + key, err);
      ok = false;
    }
    out.push_back(a);
  }
  if (!ok) return std::nullopt;
  return out;
}

std::vector<SsType> readSsTypes(const toml::table& root, Errors& e) {
  std::vector<SsType> out;
  const auto& builtins = builtinSsTypes();
  for (const toml::table* t : tablesOf(root, "ssType", "", e)) {
    checkKeys(*t, {"name", "builtin", "point"}, "ssType", e);
    auto name = get<std::string>(*t, "name", "ssType", e, true);
    auto builtin = get<bool>(*t, "builtin", "ssType", e, false);
    if (!name) continue;
    const std::string where = std::format("ssType '{}'", *name);
    SsType s;
    s.name = *name;
    s.builtin = builtin.value_or(false);
    const bool known = std::find(builtins.begin(), builtins.end(), s.name) != builtins.end();
    if (std::any_of(out.begin(), out.end(), [&](const SsType& x) { return x.name == s.name; }))
      e.add(lineOf(*t), where, "defined more than once");
    if (s.builtin) {
      if (!known) e.add(lineOf(*t), where, "is not a built-in type");
      if (t->contains("point")) e.add(lineOf(*t), where, "a built-in type cannot have points");
    } else {
      if (known || s.name.empty()) e.add(lineOf(*t), where, "this name is reserved for a built-in type (builtin = true)");
      const auto points = tablesOf(*t, "point", "ssType", e);
      if (points.empty()) e.add(lineOf(*t), where, "needs at least one [[ssType.point]]");
      double sum = 0;
      for (const toml::table* p : points) {
        checkKeys(*p, {"phi", "psi", "prob"}, "ssType.point", e);
        auto phi = get<double>(*p, "phi", where + ".point", e, true);
        auto psi = get<double>(*p, "psi", where + ".point", e, true);
        auto prob = get<double>(*p, "prob", where + ".point", e, true);
        if (!phi || !psi || !prob) continue;
        for (double v : {*phi, *psi})
          if (!(v > -180 && v <= 180)) e.add(lineOf(*p), where + ".point", "phi and psi must be within (-180, 180]");
        if (!(*prob > 0)) e.add(lineOf(*p), where + ".point.prob", "must be > 0");
        if (std::any_of(s.points.begin(), s.points.end(), [&](const SsPoint& x) { return x.phi == *phi && x.psi == *psi; }))
          e.add(lineOf(*p), where + ".point", std::format("(phi, psi) = ({}, {}) is repeated", *phi, *psi));
        sum += *prob;
        s.points.push_back({*phi, *psi, *prob});
      }
      if (sum > 1 + 1e-9) e.add(lineOf(*t), where, std::format("sum of prob is {} (must be <= 1)", sum));
    }
    out.push_back(std::move(s));
  }
  for (const char* req : {"R", "O"})
    if (std::none_of(out.begin(), out.end(), [&](const SsType& x) { return x.name == req && x.builtin; }))
      e.add(0, "ssType", std::format("built-in type '{}' is required", req));
  return out;
}

// Connected parts of the molecule without the broken bonds.
int componentCount(const Molecule& m, const std::vector<int>& broken) {
  std::vector<int> parent(m.atoms.size());
  for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
  auto find = [&](int x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  std::vector<bool> skip(m.bonds.size(), false);
  for (int b : broken) skip[b] = true;
  int n = static_cast<int>(m.atoms.size());
  for (std::size_t i = 0; i < m.bonds.size(); ++i) {
    if (skip[i]) continue;
    const int a = find(m.bonds[i].a), b = find(m.bonds[i].b);
    if (a != b) {
      parent[a] = b;
      --n;
    }
  }
  return n;
}

// Reads ringBreak and chainBreak; returns the broken bond indices.
std::vector<int> readBreaks(const toml::table& mt, const Molecule& m, const AtomRefs& refs, const std::string& where,
                            MetalBond metal, MolSettings& s, Errors& e) {
  std::map<std::pair<int, int>, int> bondOf;
  for (std::size_t i = 0; i < m.bonds.size(); ++i)
    bondOf[{std::min(m.bonds[i].a, m.bonds[i].b), std::max(m.bonds[i].a, m.bonds[i].b)}] = static_cast<int>(i);
  const std::vector<bool> ring = ringBondMask(m);
  std::set<int> used;
  auto bondIndex = [&](const toml::table& t, const std::string& field) -> int {
    auto atoms = atomList(t, "bond", 2, refs, -1, field, e);
    if (!atoms) return -1;
    auto it = bondOf.find({std::min((*atoms)[0], (*atoms)[1]), std::max((*atoms)[0], (*atoms)[1])});
    if (it == bondOf.end()) {
      e.add(lineOf(t), field + ".bond", "the two atoms are not bonded");
      return -1;
    }
    if (!used.insert(it->second).second) {
      e.add(lineOf(t), field + ".bond", "the bond is already broken by another entry");
      return -1;
    }
    return it->second;
  };
  auto metalBond = [&](int bi) { return isMetalBond(m, m.bonds[bi]); };
  for (const toml::table* t : tablesOf(mt, "ringBreak", where, e)) {
    const std::string field = where + ".ringBreak";
    checkKeys(*t, {"bond"}, field, e);
    const int bi = bondIndex(*t, field);
    if (bi < 0) continue;
    if (metalBond(bi))
      e.add(lineOf(*t), field + ".bond", "metal bonds are broken with [[molecule.chainBreak]] source = \"auto-metal\"");
    else if (!ring[bi])
      e.add(lineOf(*t), field + ".bond", "the bond is not in a ring");
    else if (m.bonds[bi].aromatic)
      e.add(lineOf(*t), field + ".bond", "aromatic bonds cannot be broken");
    else
      s.ringBreaks.push_back(bi);
  }
  for (const toml::table* t : tablesOf(mt, "chainBreak", where, e)) {
    const std::string field = where + ".chainBreak";
    checkKeys(*t, {"bond", "source"}, field, e);
    auto source = get<std::string>(*t, "source", field, e, true);
    if (source && *source != "auto-loop" && *source != "auto-crosslink" && *source != "auto-metal" &&
        *source != "user") {
      e.add(lineOf(*t), field + ".source", "must be auto-loop, auto-crosslink, auto-metal or user");
      source.reset();
    }
    const int bi = bondIndex(*t, field);
    if (bi < 0 || !source) continue;
    if (metalBond(bi) && metal == MetalBond::Keep)
      e.add(lineOf(*t), field + ".bond",
            "metalBond = \"keep\" does not break metal bonds; remove this entry or set metalBond = \"break\"");
    else if (metalBond(bi) && *source != "auto-metal")
      e.add(lineOf(*t), field + ".source", "must be auto-metal for a metal bond");
    else if (!metalBond(bi) && *source == "auto-metal")
      e.add(lineOf(*t), field + ".source", "auto-metal is only for metal bonds");
    else
      s.chainBreaks.push_back({bi, *source});
  }
  std::vector<int> broken(s.ringBreaks);
  for (const auto& c : s.chainBreaks) broken.push_back(c.bond);
  std::sort(broken.begin(), broken.end());
  return broken;
}

void readTorsions(const toml::table& owner, const std::string& key, const std::vector<int>& expected,
                  const AtomRefs& refs, int context, const std::string& where, MolSettings& s, Errors& e) {
  const Topology& t = s.topology;
  std::map<std::array<int, 4>, int> want;
  for (int i : expected) want[t.torsions[i].atoms] = i;
  std::set<int> seen;
  const std::string field = where + "." + key;
  for (const toml::table* tt : tablesOf(owner, key, where, e)) {
    checkKeys(*tt, {"atoms", "rotatable", "category"}, field, e);
    get<std::string>(*tt, "category", field, e, false);
    auto rot = get<bool>(*tt, "rotatable", field, e, true);
    auto atoms = atomList(*tt, "atoms", 4, refs, context, field, e);
    if (!atoms) continue;
    auto it = want.find({(*atoms)[0], (*atoms)[1], (*atoms)[2], (*atoms)[3]});
    if (it == want.end()) {
      e.add(lineOf(*tt), field + ".atoms", "not a torsion listed for this entry");
      continue;
    }
    if (!seen.insert(it->second).second) {
      e.add(lineOf(*tt), field + ".atoms", "listed more than once");
      continue;
    }
    if (!rot) continue;
    if (*rot) {
      if (auto why = notRotatableReason(t.torsions[it->second]); !why.empty()) {
        e.add(lineOf(*tt), field + ".rotatable", why);
        continue;
      }
    }
    s.mode[it->second] = *rot ? TorsionMode::Free : TorsionMode::Fixed;
  }
  for (int i : expected)
    if (!seen.contains(i)) {
      const auto& a = t.torsions[i].atoms;
      e.add(lineOf(owner), field,
            std::format("torsion [{}, {}, {}, {}] is missing", refs.ref(a[0], context), refs.ref(a[1], context),
                        refs.ref(a[2], context), refs.ref(a[3], context)));
    }
}

// RAMA_* types other than RAMA_AUTO must match the residue (spec 5.3, 7.5).
void checkRamaTypes(const Molecule& m, const Layout& l, const std::vector<std::string>& keyOf,
                    const std::vector<int>& ssLine, const std::string& where, const MolSettings& s, Errors& e) {
  const Topology& t = s.topology;
  std::vector<int> posOfUnit(t.units.size(), -1);
  for (std::size_t k = 0; k < l.units.size(); ++k) posOfUnit[l.units[k]] = static_cast<int>(k);
  auto resName = [&](int unit) { return unit >= 0 && t.units[unit].residue >= 0 ? m.residues[t.units[unit].residue].info.name : ""; };
  for (std::size_t k = 0; k < l.units.size(); ++k) {
    const int u = l.units[k];
    const auto given = ramaClassFromTypeName(s.ssType[u]);
    if (!given) continue;
    const std::string field = std::format("{}.residue {}.ssType", where, keyOf[k]);
    // The omega before the residue: cis when its input angle is within 90 degrees of 0; trans without one.
    bool cis = false, cisTrans = false;
    if (const int prev = t.units[u].prev; prev >= 0) {
      const auto& bb = l.backbone[posOfUnit[prev]];
      if (auto it = bb.find(BackboneAngle::Omega); it != bb.end()) {
        cis = std::abs(t.torsions[it->second].angle) < 90.0;
        cisTrans = s.mode[it->second] == TorsionMode::CisTrans;
      }
    }
    const RamaClass want = ramaClassOf(resName(u), resName(t.units[u].next), cis);
    if (cisTrans && (want == RamaClass::TransPro || want == RamaClass::CisPro) &&
        (*given == RamaClass::TransPro || *given == RamaClass::CisPro)) {
      e.add(ssLine[k], field,
            std::format("{} cannot be used: the omega before this Pro is cisTrans; use RAMA_AUTO", s.ssType[u]));
      continue;
    }
    if (*given != want)
      e.add(ssLine[k], field,
            std::format("{} is for {} residues; residue {} {} is {} (use {} or RAMA_AUTO)", s.ssType[u],
                        ramaClassLabel(*given), keyOf[k], resName(u).empty() ? "(no residue name)" : resName(u),
                        ramaClassLabel(want), ramaTypeName(want)));
  }
}

MolSettings readMolecule(const toml::table& mt, const MolInput& in, const TopologyParams& base,
                         const std::set<std::string>& ssNames, Errors& e) {
  const Molecule& m = *in.mol;
  const std::string where = std::format("molecule '{}'", m.id);
  const BackboneType bt = in.topo->backboneType;
  const bool protein = bt == BackboneType::Protein;
  MolSettings s;

  checkKeys(mt, {"id", "role", "file", "chain", "backboneType", "csgm1", "metalBond", "residue", "torsion", "ringBreak",
                 "chainBreak"},
            where, e);
  for (const char* k : {"role", "file", "chain"}) get<std::string>(mt, k, where, e, false);
  if (auto b = get<std::string>(mt, "backboneType", where, e, true); b && *b != backboneTypeName(bt))
    e.add(lineOf(*mt.get("backboneType")), where + ".backboneType",
          std::format("must be '{}' (detected; cannot be changed)", backboneTypeName(bt)));
  if (protein) {
    if (auto c = get<std::string>(mt, "csgm1", where, e, true)) {
      if (*c == "1a" || *c == "1b")
        s.csgm1 = *c == "1a" ? Csgm1::A : Csgm1::B;
      else
        e.add(lineOf(*mt.get("csgm1")), where + ".csgm1", "must be 1a or 1b");
    }
  } else if (mt.contains("csgm1")) {
    e.add(lineOf(*mt.get("csgm1")), where + ".csgm1", "only allowed for protein molecules");
  }

  if (hasMetalRing(m)) {
    if (auto v = get<std::string>(mt, "metalBond", where, e, true)) {
      if (*v == "keep" || *v == "break")
        s.metalBond = *v == "keep" ? MetalBond::Keep : MetalBond::Break;
      else
        e.add(lineOf(*mt.get("metalBond")), where + ".metalBond", "must be keep or break");
    }
  } else if (mt.contains("metalBond")) {
    e.add(lineOf(*mt.get("metalBond")), where + ".metalBond", "only allowed when a metal bond lies in a ring");
  }

  const AtomRefs refs(m);
  TopologyParams p = base;
  p.brokenBonds = readBreaks(mt, m, refs, where, s.metalBond, s, e);
  if (const int parts = componentCount(m, *p.brokenBonds) - componentCount(m, {}); parts > 0)
    e.add(lineOf(mt), where, std::format("the broken bonds (ringBreak / chainBreak) split the molecule into {} parts; "
                                         "only bonds in rings can be broken",
                                         parts + 1));
  if (s.metalBond == MetalBond::Break) {
    const std::vector<bool> ring = ringBondMask(m, *p.brokenBonds);
    for (std::size_t bi = 0; bi < m.bonds.size(); ++bi)
      if (ring[bi] && isMetalBond(m, m.bonds[bi])) {
        e.add(lineOf(*mt.get("metalBond")), where + ".metalBond",
              std::format("\"break\": metal bond [{}, {}] is still in a ring; break it or another metal bond of the "
                          "ring with [[molecule.chainBreak]] source = \"auto-metal\"",
                          refs.ref(m.bonds[bi].a, -1), refs.ref(m.bonds[bi].b, -1)));
        break;
      }
  }
  s.topology = analyzeTopology(m, p);
  s.topology.dsspSS = in.topo->dsspSS;
  s.mode.assign(s.topology.torsions.size(), TorsionMode::Fixed);
  s.ssType.assign(s.topology.units.size(), "");
  const Topology& t = s.topology;
  const Layout l = layoutOf(t);

  // Residue entries keyed by "<chain>:<seq><iCode>" or "unit <k>".
  std::map<std::string, int> posOf;
  std::vector<std::string> keyOf(l.units.size());
  for (std::size_t k = 0; k < l.units.size(); ++k) {
    const int r = t.units[l.units[k]].residue;
    keyOf[k] = r >= 0 ? residueKey(m.residues[r].info) : std::format("unit {}", k + 1);
    posOf[keyOf[k]] = static_cast<int>(k);
  }
  std::set<std::string> resKeys = {"chain", "resName", "resSeq", "iCode", "unit", "torsion"};
  if (protein) resKeys.insert({"dsspSS", "ssType"});
  if (bt != BackboneType::None)
    for (BackboneAngle a : protein ? std::vector<BackboneAngle>(std::begin(kProteinAngles), std::end(kProteinAngles))
                                   : std::vector<BackboneAngle>(std::begin(kNucleicAngles), std::end(kNucleicAngles)))
      resKeys.insert(backboneAngleName(a));

  const auto residues = tablesOf(mt, "residue", where, e);
  if (bt == BackboneType::None && !residues.empty())
    e.add(lineOf(*residues.front()), where + ".residue", "only allowed for molecules with a backbone");
  std::set<int> seen;
  std::vector<int> ssLine(l.units.size(), 0);
  for (const toml::table* rt : bt == BackboneType::None ? std::vector<const toml::table*>{} : residues) {
    const std::string rfield = where + ".residue";
    std::string key;
    if (rt->contains("unit")) {
      auto u = get<std::int64_t>(*rt, "unit", rfield, e, true);
      if (!u) continue;
      key = std::format("unit {}", *u);
    } else {
      auto chain = get<std::string>(*rt, "chain", rfield, e, true);
      auto seq = get<std::int64_t>(*rt, "resSeq", rfield, e, true);
      auto icode = get<std::string>(*rt, "iCode", rfield, e, false);
      if (!chain || !seq) continue;
      key = std::format("{}:{}{}", *chain, *seq, icode.value_or(""));
    }
    auto it = posOf.find(key);
    checkKeys(*rt, resKeys, std::format("{}.residue {}", where, key), e);
    if (it == posOf.end()) {
      e.add(lineOf(*rt), rfield, std::format("residue {} is not a backbone residue of this molecule", key));
      continue;
    }
    const int k = it->second;
    if (!seen.insert(k).second) {
      e.add(lineOf(*rt), rfield, std::format("residue {} is listed more than once", key));
      continue;
    }
    const std::string field = std::format("{}.residue {}", where, key);
    const int u = l.units[k];
    const int r = t.units[u].residue;
    if (auto name = get<std::string>(*rt, "resName", field, e, false); name && r >= 0 && *name != m.residues[r].info.name)
      e.add(lineOf(*rt->get("resName")), field + ".resName", std::format("must be '{}'", m.residues[r].info.name));
    get<std::string>(*rt, "dsspSS", field, e, false);
    if (protein) {
      if (auto ss = get<std::string>(*rt, "ssType", field, e, true)) {
        if (!ssNames.contains(*ss)) {
          e.add(lineOf(*rt->get("ssType")), field + ".ssType", std::format("ssType '{}' is not defined", *ss));
        } else {
          s.ssType[u] = *ss;
          ssLine[k] = lineOf(*rt->get("ssType"));
        }
      }
    }
    for (BackboneAngle a : protein ? std::vector<BackboneAngle>(std::begin(kProteinAngles), std::end(kProteinAngles))
                                   : std::vector<BackboneAngle>(std::begin(kNucleicAngles), std::end(kNucleicAngles))) {
      const std::string name = backboneAngleName(a);
      auto bit = l.backbone[k].find(a);
      if (bit == l.backbone[k].end()) {
        if (rt->contains(name)) e.add(lineOf(*rt->get(name)), field + "." + name, "this residue has no such backbone angle");
        continue;
      }
      const Torsion& tor = t.torsions[bit->second];
      if (a == BackboneAngle::Omega) {
        auto v = get<std::string>(*rt, name, field, e, true);
        if (!v) continue;
        if (*v == "fixed") {
          s.mode[bit->second] = TorsionMode::Fixed;
        } else if (*v == "cisTrans") {
          if (tor.broken)
            e.add(lineOf(*rt->get(name)), field + "." + name, "the bond is broken (chainBreak); must be fixed");
          else
            s.mode[bit->second] = TorsionMode::CisTrans;
        } else {
          e.add(lineOf(*rt->get(name)), field + "." + name, "must be fixed or cisTrans");
        }
        continue;
      }
      auto v = get<bool>(*rt, name, field, e, true);
      if (!v) continue;
      if (*v) {
        if (auto why = notRotatableReason(tor); !why.empty()) {
          e.add(lineOf(*rt->get(name)), field + "." + name, why);
          continue;
        }
      }
      s.mode[bit->second] = *v ? TorsionMode::Free : TorsionMode::Fixed;
    }
    readTorsions(*rt, "torsion", l.torsions[k], refs, r, field, s, e);
  }
  if (bt != BackboneType::None)
    for (std::size_t k = 0; k < l.units.size(); ++k)
      if (!seen.contains(static_cast<int>(k))) e.add(lineOf(mt), where + ".residue", std::format("residue {} is missing", keyOf[k]));
  if (protein) checkRamaTypes(m, l, keyOf, ssLine, where, s, e);
  readTorsions(mt, "torsion", l.molTorsions, refs, -1, where, s, e);
  return s;
}

}  // namespace

const std::vector<std::string>& builtinSsTypes() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> v = {"R", "O"};
    for (int i = 0; i < kRamaClassCount; ++i) v.push_back(ramaTypeName(static_cast<RamaClass>(i)));
    v.push_back("RAMA_AUTO");
    return v;
  }();
  return names;
}

bool hasMetalRing(const Molecule& m) {
  const std::vector<bool> ring = ringBondMask(m);
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi)
    if (ring[bi] && isMetalBond(m, m.bonds[bi])) return true;
  return false;
}

std::string defaultMolConfigText(const std::vector<MolInput>& mols) {
  std::string o;
  o += "# OmegaDocker molecule configuration (spec section 5).\n";
  o += "# Edit a copy and pass it with --molConfig. Entries must match the detected\n";
  o += "# molecules, residues and torsions; only the settings may be changed:\n";
  o += "#   csgm1 (\"1a\" | \"1b\"), ssType, phi / psi and alpha..zeta (true | false),\n";
  o += "#   omega (\"fixed\" | \"cisTrans\"), torsion rotatable, ringBreak, chainBreak, [[ssType]].\n";
  o += "# Atom references: <name> in the residue of the entry, <chain>:<seq><iCode>:<name>,\n";
  o += "# or #<atom serial in the input record>.\n\n";
  o += std::format("formatVersion = {}\n", kMolConfigFormatVersion);
  o += std::format("generatedBy   = {}\n", quote(std::format("OmegaDocker {}", kVersion)));
  o += "\n# ---------- secondary-structure sampling types (csgm1b only) ----------\n";
  o += "# Built-in types cannot be changed; R and O are required.\n";
  o += "#   R: uniform; O: keep the input angles; RAMA_*: Top8000 distributions (spec 7.5).\n";
  o += "#   RAMA_GLY: Gly; RAMA_TRANSPRO / RAMA_CISPRO: Pro after a trans / cis omega;\n";
  o += "#   RAMA_PREPRO: residue before Pro (not Gly, Pro); RAMA_ILEVAL: Ile, Val (not before Pro);\n";
  o += "#   RAMA_GENERAL: other residues; RAMA_AUTO: chosen per residue.\n";
  o += "# ssType is used only with csgm1 = \"1b\" for residues whose phi and psi are rotatable;\n";
  o += "# chain-end residues (without phi or psi) always use R.\n";
  o += "# Custom type: [[ssType]] name = \"A\" with [[ssType.point]] phi, psi, prob entries.\n";
  for (const auto& n : builtinSsTypes()) o += std::format("\n[[ssType]]\nname    = {}\nbuiltin = true\n", quote(n));
  o += "\n# ---------- molecules ----------\n";
  o += "# [[molecule.ringBreak]] bond = [\"<atom>\", \"<atom>\"] breaks a ring bond (spec 4.5).\n";
  o += "# [[molecule.chainBreak]] source = \"auto-loop\" | \"auto-crosslink\" | \"auto-metal\" | \"user\" (spec 4.8).\n";
  o += "#\n";
  o += "# metalBond: only for molecules with metal coordination bonds in rings (spec 4.7).\n";
  o += "#   \"break\" (default): the metal bonds listed as [[molecule.chainBreak]] source = \"auto-metal\"\n";
  o += "#            are broken and held by distance restraints (spec 4.8); the metal stays bonded\n";
  o += "#            to one ligand atom and moves with the molecule. Bonds in rings through the\n";
  o += "#            metal become rotatable again (e.g. the backbone between coordinating residues).\n";
  o += "#   \"keep\":  no metal bond is broken; every bond in a ring through the metal is not\n";
  o += "#            rotatable. Remove the auto-metal chainBreak entries when choosing \"keep\".\n";
  for (const MolInput& in : mols) writeMolecule(o, *in.mol, *in.topo);
  return o;
}

std::string writeDefaultMolConfig(const std::string& outPrefix, const std::vector<MolInput>& mols) {
  std::string path = outPrefix + ".molconfig.default.toml";
  for (int n = 1; std::filesystem::exists(path); ++n) path = std::format("{}.molconfig.default.{}.toml", outPrefix, n);
  std::ofstream f(path);
  if (!f) fail("Cannot write {}", path);
  f << defaultMolConfigText(mols);
  if (!f) fail("Cannot write {}", path);
  return path;
}

std::optional<MolConfig> readMolConfig(const std::string& path, const std::vector<MolInput>& mols,
                                       const TopologyParams& base, std::vector<std::string>& errors) {
  if (!std::filesystem::exists(path)) {
    errors.push_back(std::format("{}: file not found", path));
    return std::nullopt;
  }
  toml::table root;
  try {
    root = toml::parse_file(path);
  } catch (const toml::parse_error& ex) {
    errors.push_back(std::format("line {}: {}", ex.source().begin.line, ex.description()));
    return std::nullopt;
  }
  Errors e;
  checkKeys(root, {"formatVersion", "generatedBy", "ssType", "molecule"}, "", e);
  if (auto v = get<std::int64_t>(root, "formatVersion", "", e, true); v && *v != kMolConfigFormatVersion)
    e.add(lineOf(*root.get("formatVersion")), "formatVersion", std::format("must be {}", kMolConfigFormatVersion));
  get<std::string>(root, "generatedBy", "", e, false);

  MolConfig cfg;
  cfg.ssTypes = readSsTypes(root, e);
  std::set<std::string> ssNames;
  for (const auto& s : cfg.ssTypes) ssNames.insert(s.name);

  std::map<std::string, const toml::table*> byId;
  for (const toml::table* t : tablesOf(root, "molecule", "", e)) {
    auto id = get<std::string>(*t, "id", "molecule", e, true);
    if (!id) continue;
    if (std::none_of(mols.begin(), mols.end(), [&](const MolInput& in) { return in.mol->id == *id; }))
      e.add(lineOf(*t), "molecule", std::format("id '{}' does not match any input molecule", *id));
    else if (!byId.emplace(*id, t).second)
      e.add(lineOf(*t), "molecule", std::format("id '{}' is listed more than once", *id));
  }
  for (const MolInput& in : mols) {
    auto it = byId.find(in.mol->id);
    if (it == byId.end()) {
      e.add(0, "molecule", std::format("id '{}' is missing", in.mol->id));
      cfg.molecules.emplace_back();
      continue;
    }
    cfg.molecules.push_back(readMolecule(*it->second, in, base, ssNames, e));
  }
  if (e.count() > 0) {
    e.flush(errors);
    return std::nullopt;
  }
  return cfg;
}

}  // namespace od
