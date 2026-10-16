// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "feparams.hpp"

#include <cmath>
#include <filesystem>
#include <format>
#include <map>
#include <set>
#include <vector>

#include <toml.hpp>

#include "log.hpp"

namespace od {

namespace {

int lineOf(const toml::node& n) { return static_cast<int>(n.source().begin.line); }

struct Problems {
  std::multimap<int, std::string> list;
  void add(int line, const std::string& field, const std::string& msg) {
    list.emplace(line, std::format("line {}: {}: {}", line, field, msg));
  }
};

void checkKeys(const toml::table& t, const std::set<std::string>& allowed, const std::string& where, Problems& p) {
  for (const auto& [k, v] : t) {
    const std::string key(k.str());
    if (!allowed.count(key)) p.add(lineOf(v), where.empty() ? key : where + "." + key, "unknown key");
  }
}

void number(const toml::table& t, const std::string& key, const std::string& where, double& out, Problems& p) {
  const toml::node* n = t.get(key);
  if (!n) return;
  const std::string field = where + "." + key;
  if (!n->is_number()) {
    p.add(lineOf(*n), field, "must be a number");
    return;
  }
  const double v = n->is_integer() ? static_cast<double>(*n->value<std::int64_t>()) : *n->value<double>();
  if (!std::isfinite(v)) {
    p.add(lineOf(*n), field, "must be finite");
    return;
  }
  out = v;
}

const toml::table* subTable(const toml::table& t, const std::string& key, const std::string& where, Problems& p) {
  const toml::node* n = t.get(key);
  if (!n) return nullptr;
  if (!n->is_table()) {
    p.add(lineOf(*n), where.empty() ? key : where + "." + key, "must be a table");
    return nullptr;
  }
  return n->as_table();
}

}  // namespace

FeParams readFeParams(const std::string& path) {
  if (!std::filesystem::exists(path)) fail("{}: file not found", path);
  toml::table root;
  try {
    root = toml::parse_file(path);
  } catch (const toml::parse_error& ex) {
    fail("{}: line {}: {}", path, ex.source().begin.line, ex.description());
  }
  FeParams fp;
  Problems p;
  checkKeys(root, {"od_0100"}, "", p);
  if (const toml::table* od = subTable(root, "od_0100", "", p)) {
    checkKeys(*od, {"c0", "tors", "inter", "intra", "entropy"}, "od_0100", p);
    number(*od, "c0", "od_0100", fp.c0, p);
    number(*od, "tors", "od_0100", fp.tors, p);
    const std::set<std::string> terms(kFeTermNames.begin(), kFeTermNames.end());
    for (auto [name, arr] : {std::pair{"inter", &fp.inter}, std::pair{"intra", &fp.intra}}) {
      const std::string where = std::string("od_0100.") + name;
      if (const toml::table* t = subTable(*od, name, "od_0100", p)) {
        checkKeys(*t, terms, where, p);
        for (int k = 0; k < FeParams::kTerms; ++k) number(*t, kFeTermNames[k], where, (*arr)[k], p);
      }
    }
    if (const toml::table* t = subTable(*od, "entropy", "od_0100", p)) {
      checkKeys(*t, {"stam", "mobm"}, "od_0100.entropy", p);
      number(*t, "stam", "od_0100.entropy", fp.entropyStam, p);
      number(*t, "mobm", "od_0100.entropy", fp.entropyMobm, p);
    }
  }
  if (!p.list.empty()) {
    std::string msg = std::format("{}: {} problem(s) found:", path, p.list.size());
    for (const auto& [line, s] : p.list) msg += "\n  " + s;
    fail("{}", msg);
  }
  return fp;
}

std::string feParamsText(const FeParams& d) {
  std::string s =
      "# OmegaDocker free energy parameters (spec 10.2). Keys not given keep the defaults.\n"
      "# The defaults are not calibrated yet.\n"
      "#\n"
      "# od_0100: dG = c0 + sum_t inter.t * E_inter,t + sum_t intra.t * dE_intra,t\n"
      "#               + entropy.stam * d(-TS)_StaM + entropy.mobm * d(-TS)_MobM + tors * N_tors\n"
      "#   t = vdw, hbond, elec, desolv (AD4 terms, kcal/mol, AD4 weights included)\n"
      "#   dE_intra,t = intramolecular term t of StaM and MobM in the complex - in the unbound state\n"
      "#   d(-TS)_X   = -RT sum ln f_i over the heavy atoms of X, complex - unbound\n"
      "#   N_tors     = rotatable torsions of the MobM torsion tree\n"
      "[od_0100]\n";
  s += std::format("c0 = {}\ntors = {}\n", d.c0, d.tors);
  for (auto [name, arr] : {std::pair{"inter", &d.inter}, std::pair{"intra", &d.intra}}) {
    s += std::format("\n[od_0100.{}]\n", name);
    for (int k = 0; k < FeParams::kTerms; ++k) s += std::format("{} = {}\n", kFeTermNames[k], (*arr)[k]);
  }
  s += std::format("\n[od_0100.entropy]\nstam = {}\nmobm = {}\n", d.entropyStam, d.entropyMobm);
  return s;
}

}  // namespace od
