// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "feparams.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <map>
#include <set>
#include <variant>

#include <toml.hpp>

#include "log.hpp"

namespace od {

namespace {

// One parameter: key in the file, member, lower bound (exclusive when strict).
struct Key {
  std::string key;
  std::variant<double FeParams::*, int FeParams::*, std::pair<std::array<double, FeParams::kTerms> FeParams::*, int>> member;
  double min = -INFINITY;
  bool strict = false;
};

const std::vector<Key>& keys() {
  static const std::vector<Key> k = [] {
    std::vector<Key> v{{"od_0100.c0", &FeParams::c0},
                       {"od_0100.tors", &FeParams::tors},
                       {"od_0100.temperature", &FeParams::temperature, 0, true},
                       {"od_0100.unboundCandidateNum", &FeParams::unboundCandidateNum, 1, false}};
    for (int t = 0; t < FeParams::kTerms; ++t)
      v.push_back({std::string("od_0100.inter.") + kFeTermNames[t], std::pair{&FeParams::inter, t}});
    for (int t = 0; t < FeParams::kTerms; ++t)
      v.push_back({std::string("od_0100.intra.") + kFeTermNames[t], std::pair{&FeParams::intra, t}});
    v.push_back({"od_0100.entropy.stam", &FeParams::entropyStam});
    v.push_back({"od_0100.entropy.mobm", &FeParams::entropyMobm});
    v.push_back({"od_0100.shaking.radius", &FeParams::shakingRadius, 0, true});
    v.push_back({"od_0100.shaking.step", &FeParams::shakingStep, 0, true});
    v.push_back({"od_0100.shaking.rrtScaling", &FeParams::rrtScaling, 0, false});
    return v;
  }();
  return k;
}

const Key* findKey(const std::string& key) {
  for (const Key& k : keys())
    if (k.key == key) return &k;
  return nullptr;
}

bool isInt(const Key& k) { return std::holds_alternative<int FeParams::*>(k.member); }

// Checks v and stores it; the problem otherwise.
std::string set(FeParams& fp, const Key& k, double v) {
  if (!std::isfinite(v)) return "must be finite";
  if (isInt(k) && (v != std::floor(v) || std::abs(v) > 2147483647.0)) return "must be an integer";
  if (k.strict ? !(v > k.min) : !(v >= k.min)) return std::format("must be {} {}", k.strict ? ">" : ">=", k.min);
  std::visit(
      [&](auto m) {
        if constexpr (std::is_same_v<decltype(m), double FeParams::*>) fp.*m = v;
        else if constexpr (std::is_same_v<decltype(m), int FeParams::*>) fp.*m = static_cast<int>(v);
        else (fp.*(m.first))[m.second] = v;
      },
      k.member);
  return {};
}

double get(const FeParams& fp, const Key& k) {
  return std::visit(
      [&](auto m) -> double {
        if constexpr (std::is_same_v<decltype(m), double FeParams::*>) return fp.*m;
        else if constexpr (std::is_same_v<decltype(m), int FeParams::*>) return fp.*m;
        else return (fp.*(m.first))[m.second];
      },
      k.member);
}

int lineOf(const toml::node& n) { return static_cast<int>(n.source().begin.line); }

// Every leaf of a table checked against the keys; tables that hold keys are walked.
void walk(const toml::table& t, const std::string& where, FeParams& fp, std::multimap<int, std::string>& p) {
  for (const auto& [kk, v] : t) {
    const std::string key = where.empty() ? std::string(kk.str()) : where + "." + std::string(kk.str());
    auto add = [&](const std::string& msg) { p.emplace(lineOf(v), std::format("line {}: {}: {}", lineOf(v), key, msg)); };
    bool prefix = false;
    for (const Key& k : keys()) prefix |= k.key.starts_with(key + ".");
    if (prefix) {
      if (v.is_table()) walk(*v.as_table(), key, fp, p);
      else add("must be a table");
      continue;
    }
    const Key* k = findKey(key);
    if (!k) {
      add("unknown key");
      continue;
    }
    if (isInt(*k) ? !v.is_integer() : !v.is_number()) {
      add(isInt(*k) ? "must be an integer" : "must be a number");
      continue;
    }
    const double x = v.is_integer() ? static_cast<double>(*v.value<std::int64_t>()) : *v.value<double>();
    if (const std::string e = set(fp, *k, x); !e.empty()) add(e);
  }
}

}  // namespace

FeParams loadFeParams(const std::string& path, const FeOverrides& overrides) {
  FeParams fp;
  std::vector<std::string> problems;
  if (!path.empty()) {
    if (!std::filesystem::exists(path)) fail("{}: file not found", path);
    toml::table root;
    try {
      root = toml::parse_file(path);
    } catch (const toml::parse_error& ex) {
      fail("{}: line {}: {}", path, ex.source().begin.line, ex.description());
    }
    std::multimap<int, std::string> p;
    walk(root, "", fp, p);
    for (auto& [line, s] : p) problems.push_back(path + ": " + s);
  }
  for (const auto& [key, text] : overrides) {
    const Key* k = findKey(key);
    if (!k) fail("unknown free energy parameter '{}'", key);
    double v = 0;
    const auto r = std::from_chars(text.data(), text.data() + text.size(), v);
    std::string e = (r.ec != std::errc{} || r.ptr != text.data() + text.size()) ? "must be a number" : set(fp, *k, v);
    if (!e.empty()) problems.push_back(std::format("command line: {} ({}): {}", key, text, e));
  }
  if (!problems.empty()) {
    std::string msg = std::format("free energy parameters: {} problem(s) found:", problems.size());
    for (const auto& s : problems) msg += "\n  " + s;
    fail("{}", msg);
  }
  return fp;
}

std::string feParamsText(const FeParams& d) {
  std::string s =
      "# OmegaDocker free energy parameters (spec 10.2). Keys not given keep the defaults;\n"
      "# command-line options override this file. The defaults are not calibrated yet.\n"
      "#\n"
      "# od_0100: dG = c0 + sum_t inter.t * E_inter,t + sum_t intra.t * dE_intra,t\n"
      "#               + entropy.stam * d(-TS)_StaM + entropy.mobm * d(-TS)_MobM + tors * N_tors\n"
      "#   t = vdw, hbond, elec, desolv (AD4 terms, kcal/mol, AD4 weights included)\n"
      "#   dE_intra,t = intramolecular term t of StaM and MobM in the complex - in the unbound state\n"
      "#   d(-TS)_X   = -R * temperature * sum ln f_i over the heavy atoms of X, complex - unbound\n"
      "#   f_i        = sum_k exp(-(E_k - E_min) * shaking.rrtScaling) / points within shaking.radius\n"
      "#                (cubic lattice of spacing shaking.step)\n"
      "#   N_tors     = rotatable torsions of the MobM torsion tree\n"
      "#   unboundCandidateNum: 9a samples with the lowest energy evaluated for the unbound state\n";
  std::string table;
  for (const Key& k : keys()) {
    const auto dot = k.key.rfind('.');
    const std::string t = k.key.substr(0, dot);
    if (t != table) {
      s += std::format("\n[{}]\n", t);
      table = t;
    }
    s += std::format("{} = {}\n", k.key.substr(dot + 1), isInt(k) ? std::format("{}", static_cast<int>(get(d, k)))
                                                                   : std::format("{}", get(d, k)));
  }
  return s;
}

}  // namespace od
