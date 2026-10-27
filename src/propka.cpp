// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// C++ port of PROPKA 3 (LGPL-2.1, https://github.com/jensengroup/propka, commit 58ceb7b).
// The structure follows the Python modules (bonds, ligand, protonate, group, energy,
// determinants, iterative, conformation_container); orders of atoms, bonds and groups
// are kept so that the results match the Python program.

#include "propka.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <format>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <unordered_map>

#include "log.hpp"

extern "C" const char od_propka_cfg_begin[];
extern "C" const char od_propka_cfg_end[];

namespace od::propka {
namespace {

// ---------- parameters (propka.cfg) ----------

struct Params {
  std::vector<std::string> matrixKeys;
  std::map<std::string, std::map<std::string, char>> matrix;
  std::pair<double, double> cutoffDefault{0, 0};
  std::map<std::string, std::map<std::string, std::pair<double, double>>> cutoffs;
  std::map<std::string, double> vdwVolume, charge, modelPkas, ions, customModelPkas;
  std::map<std::string, std::vector<double>> bbNH, bbCO;
  std::map<std::string, std::string> proteinGroupMapping;
  std::vector<std::string> ignoreResidues, angular, acidList, baseList, bbReorg;
  double desolvCutoff = 20, buriedCutoff = 15, coulombCutoff1 = 4, coulombCutoff2 = 10;
  int nmin = 280, nmax = 560;
  double surfScale = 0.25, desolvPrefactor = -13, desolvAllowance = 0;
  double cooHis = 1.6, ocoHis = 1.6, cysHis = 1.6, cysCys = 3.6;
  int couplingMaxBonds = 3, minBondDistHb = 4;
  bool remove_penalised = true;
  double sidechainInteraction = 0.85;

  char matrixValue(const std::string& a, const std::string& b) const {
    auto i = matrix.find(a);
    if (i == matrix.end()) return 0;
    auto j = i->second.find(b);
    return j == i->second.end() ? 0 : j->second;
  }
  std::pair<double, double> cutoff(const std::string& a, const std::string& b) const {
    auto i = cutoffs.find(a);
    if (i == cutoffs.end()) return cutoffDefault;
    auto j = i->second.find(b);
    return j == i->second.end() ? cutoffDefault : j->second;
  }
  bool inList(const std::vector<std::string>& l, const std::string& v) const {
    return std::find(l.begin(), l.end(), v) != l.end();
  }
};

double toD(const std::string& s) { return std::stod(s); }

Params loadParams() {
  Params p;
  std::istringstream in(std::string(od_propka_cfg_begin, od_propka_cfg_end));
  std::string line;
  while (std::getline(in, line)) {
    if (auto c = line.find('#'); c != std::string::npos) line.resize(c);
    std::istringstream ls(line);
    std::vector<std::string> w;
    for (std::string x; ls >> x;) w.push_back(x);
    if (w.empty()) continue;
    const std::string& k = w[0];
    auto numDict = [&](std::map<std::string, double>& d) { d[w[1]] = toD(w[2]); };
    if (k == "interaction_matrix") {
      if (w.size() != p.matrixKeys.size() + 3) fail("propka.cfg: interaction_matrix row {} has {} values", w[1], w.size() - 2);
      const std::string g = w[1];
      p.matrixKeys.push_back(g);
      for (std::size_t i = 0; i < p.matrixKeys.size(); ++i) {
        p.matrix[p.matrixKeys[i]][g] = w[i + 2][0];
        p.matrix[g][p.matrixKeys[i]] = w[i + 2][0];
      }
    } else if (k == "sidechain_cutoffs") {
      if (w.size() == 4 && w[1] == "default") {
        p.cutoffDefault = {toD(w[2]), toD(w[3])};
      } else {
        p.cutoffs[w[1]][w[2]] = {toD(w[3]), toD(w[4])};
        p.cutoffs[w[2]][w[1]] = {toD(w[3]), toD(w[4])};
      }
    } else if (k == "VanDerWaalsVolume") numDict(p.vdwVolume);
    else if (k == "charge") numDict(p.charge);
    else if (k == "model_pkas") numDict(p.modelPkas);
    else if (k == "ions") numDict(p.ions);
    else if (k == "custom_model_pkas") numDict(p.customModelPkas);
    else if (k == "backbone_NH_hydrogen_bond" || k == "backbone_CO_hydrogen_bond") {
      auto& d = k == "backbone_NH_hydrogen_bond" ? p.bbNH : p.bbCO;
      for (std::size_t i = 2; i < w.size(); ++i) d[w[1]].push_back(toD(w[i]));
    } else if (k == "protein_group_mapping") p.proteinGroupMapping[w[1]] = w[2];
    else if (k == "ignore_residues") p.ignoreResidues.push_back(w[1]);
    else if (k == "angular_dependent_sidechain_interactions") p.angular.push_back(w[1]);
    else if (k == "acid_list") p.acidList.push_back(w[1]);
    else if (k == "base_list") p.baseList.push_back(w[1]);
    else if (k == "backbone_reorganisation_list") p.bbReorg.push_back(w[1]);
    else if (k == "desolv_cutoff") p.desolvCutoff = toD(w[1]);
    else if (k == "buried_cutoff") p.buriedCutoff = toD(w[1]);
    else if (k == "coulomb_cutoff1") p.coulombCutoff1 = toD(w[1]);
    else if (k == "coulomb_cutoff2") p.coulombCutoff2 = toD(w[1]);
    else if (k == "Nmin") p.nmin = std::stoi(w[1]);
    else if (k == "Nmax") p.nmax = std::stoi(w[1]);
    else if (k == "desolvationSurfaceScalingFactor") p.surfScale = toD(w[1]);
    else if (k == "desolvationPrefactor") p.desolvPrefactor = toD(w[1]);
    else if (k == "desolvationAllowance") p.desolvAllowance = toD(w[1]);
    else if (k == "COO_HIS_exception") p.cooHis = toD(w[1]);
    else if (k == "OCO_HIS_exception") p.ocoHis = toD(w[1]);
    else if (k == "CYS_HIS_exception") p.cysHis = toD(w[1]);
    else if (k == "CYS_CYS_exception") p.cysCys = toD(w[1]);
    else if (k == "coupling_max_number_of_bonds") p.couplingMaxBonds = std::stoi(w[1]);
    else if (k == "min_bond_distance_for_hydrogen_bonds") p.minBondDistHb = std::stoi(w[1]);
    else if (k == "remove_penalised_group") p.remove_penalised = std::stoi(w[1]) != 0;
    else if (k == "sidechain_interaction") p.sidechainInteraction = toD(w[1]);
    // Other keys (version, ligand_typing = groups, shared_determinants = 0, common_charge_centre = 0,
    // non-covalent coupling, write_out_order, ...) have the values this port implements.
  }
  return p;
}

const Params& params() {
  static const Params p = loadParams();
  return p;
}

// ---------- atoms ----------

struct PAtom {
  std::string name, resName /* padded to 3 */, chain, icode /* " " when empty */, element, type, terminal;
  int resNum = 0;
  double x = 0, y = 0, z = 0;
  std::vector<int> bonded;
  int group = -1;
  std::string groupType;
  bool cysteineBridge = false, isProtonated = false, stericSet = false, sybylAssigned = false;
  double charge = 0;
  int stericNumber = 0, lonePairs = 0, protonsToAdd = 0, piElec23 = 0, piElecConj23 = 0;
  std::string sybylType;
  int input = -1;  // index of the input atom, -1 for added hydrogens
};

struct Mol;

double sqDist(double ax, double ay, double az, double bx, double by, double bz) {
  const double dx = bx - ax, dy = by - ay, dz = bz - az;
  return dx * dx + dy * dy + dz * dz;
}

// Python round(x, 3).
double round3(double v) {
  const std::string s = std::format("{:.3f}", v);
  double r = 0;
  std::from_chars(s.data(), s.data() + s.size(), r);
  return r;
}

struct Vec {
  double x = 0, y = 0, z = 0;
  Vec operator+(const Vec& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec operator-(const Vec& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec operator-() const { return {-x, -y, -z}; }
  double dot(const Vec& o) const { return x * o.x + y * o.y + z * o.z; }
  Vec cross(const Vec& o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
  double sqLength() const { return x * x + y * y + z * z; }
  double length() const { return std::sqrt(sqLength()); }
  Vec rescale(double l) const {
    const double f = l / length();
    return {x * f, y * f, z * f};
  }
  Vec orthogonal() const {
    Vec r{y, -x, 0};
    if (std::abs(y) < std::abs(z)) r = {z, 0, -x};
    return r;
  }
};

struct Mat {  // 3x4 part of Matrix4x4
  double a[3][4] = {};
  Vec operator*(const Vec& v) const {
    return {a[0][0] * v.x + a[0][1] * v.y + a[0][2] * v.z + a[0][3], a[1][0] * v.x + a[1][1] * v.y + a[1][2] * v.z + a[1][3],
            a[2][0] * v.x + a[2][1] * v.y + a[2][2] * v.z + a[2][3]};
  }
};
Mat rotZ(double t) {
  Mat m;
  m.a[0][0] = std::cos(t); m.a[0][1] = -std::sin(t);
  m.a[1][0] = std::sin(t); m.a[1][1] = std::cos(t);
  m.a[2][2] = 1;
  return m;
}
Mat rotY(double t) {
  Mat m;
  m.a[0][0] = std::cos(t); m.a[0][2] = std::sin(t);
  m.a[1][1] = 1;
  m.a[2][0] = -std::sin(t); m.a[2][2] = std::cos(t);
  return m;
}
Vec rotateAroundAxis(double theta, Vec axis, Vec v) {
  double gamma = 0;
  if (axis.y != 0) {
    if (axis.x != 0)
      gamma = -axis.x / std::abs(axis.x) * std::asin(axis.y / std::sqrt(axis.x * axis.x + axis.y * axis.y));
    else
      gamma = M_PI / 2.0;
    const Mat r = rotZ(gamma);
    v = r * v;
    axis = r * axis;
  }
  double beta = 0;
  if (axis.x != 0) {
    beta = -axis.x / std::abs(axis.x) * std::acos(axis.z / std::sqrt(axis.x * axis.x + axis.z * axis.z));
    const Mat r = rotY(beta);
    v = r * v;
    axis = r * axis;
  }
  v = rotZ(theta) * v;
  v = rotY(-beta) * v;
  v = rotZ(-gamma) * v;
  return v;
}

// ---------- groups ----------

struct Det {
  int group;
  std::string label;
  double value;
};
enum DetType { kSidechain = 0, kBackbone = 1, kCoulomb = 2 };

struct Group {
  int atom = -1;
  std::string type, residueType, label;
  std::vector<Det> dets[3];
  double pka = 0, modelPka = 0, energyVolume = 0, energyLocal = 0, buried = 0;
  int numVolume = 0;
  double x = 0, y = 0, z = 0;
  double charge = 0;
  std::vector<int> acidAtoms, baseAtoms;
  bool modelPkaSet = false, titratable = false;
  std::vector<int> covCoupled;
  int coupledTitrating = -1;
};

struct Mol {
  std::vector<PAtom> atoms;
  std::vector<Group> groups;
  const Params& P = params();

  Vec pos(int i) const { return {atoms[i].x, atoms[i].y, atoms[i].z}; }
  std::vector<int> bondedElements(int a, const std::string& e) const {
    std::vector<int> r;
    for (int b : atoms[a].bonded)
      if (atoms[b].element == e) r.push_back(b);
    return r;
  }
  std::vector<int> bondedHeavy(int a) const {
    std::vector<int> r;
    for (int b : atoms[a].bonded)
      if (atoms[b].element != "H") r.push_back(b);
    return r;
  }
  bool withinBondDistance(int a, int other, int maxBonds, int cur) const {
    for (int b : atoms[a].bonded) {
      if (b == other) return true;
      if (maxBonds > cur && withinBondDistance(b, other, maxBonds, cur + 1)) return true;
    }
    return false;
  }
  std::string residueKey(int a) const {
    return std::format("{}|{}|{}", atoms[a].chain, atoms[a].resNum, atoms[a].icode);
  }
  bool groupEq(int g1, int g2) const {
    return groups[g1].label == groups[g2].label && residueKey(groups[g1].atom) == residueKey(groups[g2].atom);
  }
  bool groupIn(int g, const std::vector<int>& l) const {
    for (int x : l)
      if (groupEq(g, x)) return true;
    return false;
  }

  // ----- bonds (bonds.py) -----
  void makeBond(int a, int b) {
    if (a == b) return;
    if (std::find(atoms[b].bonded.begin(), atoms[b].bonded.end(), a) == atoms[b].bonded.end()) atoms[b].bonded.push_back(a);
    if (std::find(atoms[a].bonded.begin(), atoms[a].bonded.end(), b) == atoms[a].bonded.end()) atoms[a].bonded.push_back(b);
  }
  bool checkDistance(int a, int b) const {
    const double d = sqDist(atoms[a].x, atoms[a].y, atoms[a].z, atoms[b].x, atoms[b].y, atoms[b].z);
    if (d > 6.25) return false;
    const std::string key = atoms[a].element + "-" + atoms[b].element;
    const int h = static_cast<int>(std::count(key.begin(), key.end(), 'H'));
    if (d < 1.5 * 1.5 && h == 1) return true;
    if (d < 2.0 * 2.0 && h == 0) return true;
    if (key == "S-S" && d < 2.5 * 2.5) return true;
    if (key == "F-F" && d < 1.7 * 1.7) return true;
    return false;
  }
  void tryBond(int a, int b) {
    const auto& ba = atoms[b].bonded;
    if (std::find(ba.begin(), ba.end(), a) != ba.end()) return;
    if (checkDistance(a, b)) {
      makeBond(a, b);
      if (atoms[a].element == "S" && atoms[b].element == "S") atoms[a].cysteineBridge = atoms[b].cysteineBridge = true;
    }
  }
  void findBondsUsingBoxes() {
    const double box = std::max(2.5, std::sqrt(2.5 * 2.5) + 0.01);
    std::vector<std::array<long, 3>> keys;
    std::map<std::array<long, 3>, std::vector<int>> boxes;
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
      const std::array<long, 3> k{static_cast<long>(std::floor(atoms[i].x / box)), static_cast<long>(std::floor(atoms[i].y / box)),
                                  static_cast<long>(std::floor(atoms[i].z / box))};
      auto [it, added] = boxes.try_emplace(k);
      if (added) keys.push_back(k);
      it->second.push_back(i);
    }
    static const int nb[13][3] = {{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1},
                                  {-1, 1, 0},   {-1, 1, 1},  {0, -1, -1}, {0, -1, 0},  {0, -1, 1}, {0, 0, -1}};
    for (const auto& k : keys) {
      const std::vector<int>& v = boxes[k];
      for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = i + 1; j < v.size(); ++j) tryBond(v[i], v[j]);
      for (const auto& d : nb) {
        auto it = boxes.find({k[0] + d[0], k[1] + d[1], k[2] + d[2]});
        if (it == boxes.end()) continue;
        for (int a : v)
          for (int b : it->second) tryBond(a, b);
      }
    }
  }
  void addPiElectronInfo() {
    static const std::map<std::string, int> side{
        {"ARG-CZ", 1},  {"ARG-NH1", 1}, {"ASN-OD1", 1}, {"ASN-CG", 1},  {"ASP-OD1", 1}, {"ASP-CG", 1},  {"GLU-OE1", 1},
        {"GLU-CD", 1},  {"GLN-OE1", 1}, {"GLN-CD", 1},  {"HIS-CG", 1},  {"HIS-CD2", 1}, {"HIS-ND1", 1}, {"HIS-CE1", 1},
        {"PHE-CG", 1},  {"PHE-CD1", 1}, {"PHE-CE1", 1}, {"PHE-CZ", 1},  {"PHE-CE2", 1}, {"PHE-CD2", 1}, {"TRP-CG", 1},
        {"TRP-CD1", 1}, {"TRP-CE2", 1}, {"TRP-CD2", 1}, {"TRP-CE3", 1}, {"TRP-CZ3", 1}, {"TRP-CH2", 1}, {"TRP-CZ2", 1},
        {"TYR-CG", 1},  {"TYR-CD1", 1}, {"TYR-CE1", 1}, {"TYR-CZ", 1},  {"TYR-CE2", 1}, {"TYR-CD2", 1}};
    static const std::map<std::string, int> sideConj{{"ARG-NE", 1},  {"ARG-NH2", 1}, {"ASN-ND2", 1},
                                                     {"GLN-NE2", 1}, {"HIS-NE2", 1}, {"TRP-NE1", 1}};
    static const std::map<std::string, int> lig{{"C.ar", 1}, {"N.pl3", 0}, {"C.2", 1}, {"O.2", 1},
                                                {"O.co2", 1}, {"N.ar", 1}, {"C.1", 2}, {"N.1", 2}};
    static const std::map<std::string, int> ligConj{{"N.am", 1}, {"N.pl3", 1}};
    for (PAtom& a : atoms) {
      if (a.type == "hetatm") {
        if (auto it = lig.find(a.sybylType); it != lig.end()) a.piElec23 = it->second;
        if (auto it = ligConj.find(a.sybylType); it != ligConj.end()) a.piElecConj23 = it->second;
      }
      if (a.type == "atom") {
        const std::string key = a.resName + "-" + a.name;
        if (auto it = side.find(key); it != side.end()) a.piElec23 = it->second;
        if (auto it = sideConj.find(key); it != sideConj.end()) a.piElecConj23 = it->second;
        if (a.name == "C" || a.name == "O") a.piElec23 = 1;
        if (a.name == "N" && a.bonded.size() > 1) a.piElecConj23 = 1;
      }
    }
  }

  // ----- ligand typing (ligand.py) -----
  std::vector<int> identifyRing(int self, int original, int number, std::vector<int> past) const {
    ++number;
    past.push_back(self);
    std::vector<int> ret;
    if (number > 10) return ret;
    for (int a : bondedHeavy(self)) {
      if (a == original && number > 2) return past;
      if (std::find(past.begin(), past.end(), a) == past.end()) {
        std::vector<int> r = identifyRing(a, original, number, past);
        if (!r.empty() && (ret.size() > r.size() || ret.empty())) ret = r;
      }
    }
    return ret;
  }
  std::vector<int> ringMember(int a) const { return identifyRing(a, a, 0, {}); }
  bool atomsPlanar(const std::vector<int>& v) const {
    if (v.size() < 4) return false;
    const Vec v1 = pos(v[1]) - pos(v[0]), v2 = pos(v[2]) - pos(v[0]);
    Vec n = v1.cross(v2);
    if (n.sqLength() == 0.0) return false;
    n = n.rescale(1.0);
    for (std::size_t i = 3; i < v.size(); ++i) {
      Vec w = pos(v[i]) - pos(v[0]);
      if (w.sqLength() == 0.0) return false;
      w = w.rescale(1.0);
      if (std::abs(w.dot(n)) > 0.20) return false;
    }
    return true;
  }
  bool planar(int a) const {
    std::vector<int> v{a};
    v.insert(v.end(), atoms[a].bonded.begin(), atoms[a].bonded.end());
    return atomsPlanar(v);
  }
  bool aromaticRing(const std::vector<int>& r) const {
    if (r.size() < 5) return false;
    for (std::size_t i = 0; i < r.size(); ++i) {
      std::vector<int> v(r.begin() + static_cast<long>(i), r.end());
      v.insert(v.end(), r.begin(), r.begin() + static_cast<long>(i));
      if (!atomsPlanar(v)) return false;
    }
    return true;
  }
  void setType(int a, const std::string& t) {
    atoms[a].sybylType = t;
    atoms[a].sybylAssigned = true;
  }
  void assignSybylType(int a) {
    PAtom& at = atoms[a];
    if (at.sybylAssigned) return;
    const std::vector<int> ring = ringMember(a);
    const bool aromatic = aromaticRing(ring);
    const bool pl = planar(a);
    std::vector<std::string> bondedEl;
    for (int b : at.bonded) bondedEl.push_back(atoms[b].element);
    auto countEl = [&](const std::string& e) { return static_cast<int>(std::count(bondedEl.begin(), bondedEl.end(), e)); };
    auto indexEl = [&](const std::string& e, std::size_t from) {
      return static_cast<int>(std::find(bondedEl.begin() + static_cast<long>(from), bondedEl.end(), e) - bondedEl.begin());
    };
    if (aromatic) {
      for (int r : ring)
        if (atoms[r].element == "C" || atoms[r].element == "N") setType(r, atoms[r].element + ".ar");
      return;
    }
    const std::string el = at.element;
    if (el == "O" || el == "N" || el == "C") {
      int o = -1, n = -1, c = -1;
      if (el == "O" || el == "N") {
        for (int be : bondedElements(a, "C"))
          for (int ba : atoms[be].bonded) {
            if (atoms[ba].element == "N" && el == "O") { o = a; c = be; n = ba; }
            else if (atoms[ba].element == "O" && el == "N") { n = a; c = be; o = ba; }
          }
      }
      if (el == "C") {
        const auto ns = bondedElements(a, "N"), os = bondedElements(a, "O");
        if (ns.size() == 1 && os.size() == 1) { c = a; n = ns[0]; o = os[0]; }
      }
      if (c >= 0 && n >= 0 && o >= 0) {
        if (!aromaticRing(ringMember(n)) && bondedHeavy(n).size() == 2) {
          setType(n, "N.am");
          setType(c, "C.2");
          setType(o, "O.2");
          return;
        }
      }
    }
    if (el == "C") {
      if (at.bonded.size() == 3 && countEl("O") == 2) {
        const int i1 = indexEl("O", 0), i2 = indexEl("O", static_cast<std::size_t>(i1) + 1);
        if (atoms[at.bonded[i1]].bonded.size() == 1 && atoms[at.bonded[i2]].bonded.size() == 1) {
          setType(at.bonded[i1], "O.co2-");
          setType(at.bonded[i2], "O.co2");
          setType(a, "C.2");
          return;
        }
      }
      if (at.bonded.size() <= 2) {
        for (int b : std::vector<int>(at.bonded)) {
          if (sqDist(at.x, at.y, at.z, atoms[b].x, atoms[b].y, atoms[b].z) < 1.2 * 1.2) {
            setType(a, "C.1");
            setType(b, atoms[b].element + ".1");
          }
        }
        if (atoms[a].sybylAssigned) return;
      }
      if (pl) {
        setType(a, "C.2");
        for (int b : std::vector<int>(atoms[a].bonded))
          if (atoms[b].element == "N" && (atoms[b].bonded.size() < 3 || planar(b))) setType(b, "N.pl3");
        return;
      }
      setType(a, "C.3");
      return;
    }
    if (el == "N") {
      if (at.bonded.size() == 1 && planar(at.bonded[0])) {
        setType(a, "N.pl3");
        return;
      }
      if (pl) {
        setType(a, "N.pl3");
        return;
      }
      setType(a, "N.3");
      return;
    }
    if (el == "O") {
      setType(a, "O.3");
      if (at.bonded.size() == 1) {
        if (atoms[at.bonded[0]].element == "C") {
          const int c = at.bonded[0];
          if (atoms[c].bonded.size() == 3 && bondedElements(c, "O").size() == 2) {
            const auto os = bondedElements(c, "O");
            if (atoms[os[0]].bonded.size() == 1 && atoms[os[1]].bonded.size() == 1) {
              setType(os[0], "O.co2-");
              setType(os[1], "O.co2");
              setType(c, "C.2");
              return;
            }
          }
        }
        const int b = at.bonded[0];
        if (sqDist(at.x, at.y, at.z, atoms[b].x, atoms[b].y, atoms[b].z) < 1.3 * 1.3) {
          setType(a, "O.2");
          if (atoms[b].element == "C") setType(b, "C.2");
        }
      }
      return;
    }
    if (el == "S") {
      if (countEl("O") == 2) {
        const int i1 = indexEl("O", 0), i2 = indexEl("O", static_cast<std::size_t>(i1) + 1);
        setType(at.bonded[i1], "O.2");
        setType(at.bonded[i2], "O.2");
        setType(a, "S.o2");
        return;
      }
      if (countEl("O") == 4) {
        int no2 = 0;
        for (int b : at.bonded) {
          if (atoms[b].bonded.size() == 1 && no2 < 2) {
            setType(b, "O.2");
            ++no2;
          } else {
            setType(b, "O.3");
          }
        }
      }
      setType(a, "S.3");
      return;
    }
    if (el == "P") {
      setType(a, "P.3");
      for (int o : bondedElements(a, "O")) setType(o, "O.3");
      return;
    }
    std::string e = el;
    for (std::size_t i = 1; i < e.size(); ++i) e[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(e[i])));
    setType(a, e);
  }

  // ----- protonation (protonate.py) -----
  static int valenceElectrons(const std::string& e) {
    static const std::map<std::string, int> v{
        {"H", 1},  {"He", 2}, {"Li", 1}, {"Be", 2}, {"B", 3},  {"C", 4},  {"N", 5},  {"O", 6},  {"F", 7},  {"Ne", 8}, {"Na", 1},
        {"Mg", 2}, {"Al", 3}, {"Si", 4}, {"P", 5},  {"S", 6},  {"Cl", 7}, {"Ar", 8}, {"K", 1},  {"Ca", 2}, {"Sc", 2}, {"Ti", 2},
        {"V", 2},  {"Cr", 1}, {"Mn", 2}, {"Fe", 2}, {"Co", 2}, {"Ni", 2}, {"Cu", 1}, {"Zn", 2}, {"Ga", 3}, {"Ge", 4}, {"As", 5},
        {"Se", 6}, {"Br", 7}, {"Kr", 8}, {"Rb", 1}, {"Sr", 2}, {"Y", 2},  {"Zr", 2}, {"Nb", 1}, {"Mo", 1}, {"Tc", 2}, {"Ru", 1},
        {"Rh", 1}, {"Pd", 8}, {"Ag", 1}, {"Cd", 2}, {"In", 3}, {"Sn", 4}, {"Sb", 5}, {"Te", 6}, {"I", 7},  {"Xe", 8}, {"Cs", 1},
        {"Ba", 2}, {"La", 2}, {"Ce", 2}, {"Pr", 2}, {"Nd", 2}, {"Pm", 2}, {"Sm", 2}, {"Eu", 2}, {"Gd", 2}, {"Tb", 2}, {"Dy", 2},
        {"Ho", 2}, {"Er", 2}, {"Tm", 2}, {"Yb", 2}, {"Lu", 2}, {"Hf", 2}, {"Ta", 2}, {"W", 2},  {"Re", 2}, {"Os", 2}, {"Ir", 2},
        {"Pt", 1}, {"Au", 1}, {"Hg", 2}, {"Tl", 3}, {"Pb", 4}, {"Bi", 5}, {"Po", 6}, {"At", 7}, {"Rn", 8}, {"Fr", 1}, {"Ra", 2},
        {"Ac", 2}, {"Th", 2}, {"Pa", 2}, {"U", 2},  {"Np", 2}, {"Pu", 2}, {"Am", 2}, {"Cm", 2}, {"Bk", 2}, {"Cf", 2}, {"Es", 2},
        {"Fm", 2}, {"Md", 2}, {"No", 2}, {"Lr", 3}};
    auto it = v.find(e);
    return it == v.end() ? 4 : it->second;
  }
  void setCharge(int a) {
    PAtom& at = atoms[a];
    if (at.type == "atom") {
      static const std::map<std::string, double> std{{"ARG-NH1", 1.0}, {"ASP-OD2", -1.0}, {"GLU-OE2", -1.0}, {"HIS-ND1", 1.0},
                                                     {"LYS-NZ", 1.0},  {"N+", 1.0},      {"C-", -1.0}};
      std::string key = at.resName + "-" + at.name;
      if (!at.terminal.empty()) key = at.terminal;
      if (auto it = std.find(key); it != std.end()) at.charge = it->second;
    } else if (at.type == "hetatm") {
      static const std::map<std::string, double> syb{{"N.pl3", 1}, {"N.3", 1}, {"N.4", 1}, {"N.ar", 1}, {"O.co2-", 1}};
      if (auto it = syb.find(at.sybylType); it != syb.end()) {
        at.charge = it->second;
        at.sybylType.erase(std::remove(at.sybylType.begin(), at.sybylType.end(), '-'), at.sybylType.end());
      }
    }
  }
  void setProtonsToAdd(int a) {
    PAtom& at = atoms[a];
    at.protonsToAdd = 8 - valenceElectrons(at.element) - static_cast<int>(at.bonded.size()) - at.piElec23 +
                      static_cast<int>(at.charge);
  }
  void setStericNumber(int a) {
    PAtom& at = atoms[a];
    if (at.stericSet) return;
    int s = valenceElectrons(at.element) + static_cast<int>(at.bonded.size()) + at.protonsToAdd - at.piElec23 - at.piElecConj23;
    at.stericNumber = static_cast<int>(std::floor((s - at.charge) / 2));
    at.lonePairs = at.stericNumber - static_cast<int>(at.bonded.size()) - at.protonsToAdd;
    at.stericSet = true;
  }
  Vec setBondDistance(const Vec& v, const std::string& e) const {
    static const std::map<std::string, double> L{{"C", 1.09}, {"N", 1.01}, {"O", 0.96}, {"F", 0.92},
                                                 {"Cl", 1.27}, {"Br", 1.41}, {"I", 1.61}, {"S", 1.35}};
    auto it = L.find(e);
    return v.rescale(it == L.end() ? 1.0 : it->second);
  }
  void addProton(int a, const Vec& p) {
    PAtom h;
    h.name = "H" + atoms[a].name.substr(std::min<std::size_t>(1, atoms[a].name.size()));
    h.resName = atoms[a].resName;
    h.chain = atoms[a].chain;
    h.resNum = atoms[a].resNum;
    h.icode = atoms[a].icode;
    h.x = round3(p.x);
    h.y = round3(p.y);
    h.z = round3(p.z);
    h.element = "H";
    h.type = atoms[a].type;
    h.bonded = {a};
    const int hi = static_cast<int>(atoms.size());
    atoms.push_back(std::move(h));
    atoms[a].bonded.push_back(hi);
    atoms[a].protonsToAdd -= 1;
  }
  void trigonal(int a) {
    const double rot = 120.0 * M_PI / 180.0;
    const Vec c = pos(a);
    const auto& bd = atoms[a].bonded;
    if (bd.size() == 1 && atoms[a].protonsToAdd > 0) {
      const int n0 = bd[0];
      Vec av = pos(n0) - c;
      setStericNumber(n0);
      Vec axis;
      if (atoms[n0].stericNumber == 3 && atoms[n0].bonded.size() > 1) {
        std::vector<int> others;
        for (std::size_t i = 0; i < atoms[n0].bonded.size(); ++i)
          if (atoms[n0].bonded[i] != a) others.push_back(static_cast<int>(i));
        const Vec v1 = pos(n0) - c;
        const Vec v2 = pos(atoms[n0].bonded[others[0]]) - pos(n0);
        axis = v1.cross(v2);
        if (others.size() > 1) {
          const Vec v3 = pos(atoms[n0].bonded[others[1]]) - pos(n0);
          const Vec axis2 = v1.cross(v3);
          axis = axis.dot(axis2) > 0 ? axis + axis2 : axis - axis2;
        }
      } else {
        axis = av.orthogonal();
      }
      av = rotateAroundAxis(rot, axis, av);
      av = setBondDistance(av, atoms[a].element);
      addProton(a, c + av);
    }
    if (atoms[a].bonded.size() == 2 && atoms[a].protonsToAdd > 0) {
      const Vec a1 = (pos(atoms[a].bonded[0]) - c).rescale(1.0), a2 = (pos(atoms[a].bonded[1]) - c).rescale(1.0);
      addProton(a, c + setBondDistance(-a1 - a2, atoms[a].element));
    }
  }
  void tetrahedral(int a) {
    const double rot = 109.5 * M_PI / 180.0;
    const Vec c = pos(a);
    if (atoms[a].bonded.size() == 1 && atoms[a].protonsToAdd > 0) {
      Vec av = pos(atoms[a].bonded[0]) - c;
      const Vec axis = av.orthogonal();
      av = rotateAroundAxis(rot, axis, av);
      addProton(a, c + setBondDistance(av, atoms[a].element));
    }
    if (atoms[a].bonded.size() == 2 && atoms[a].protonsToAdd > 0) {
      const Vec a1 = (pos(atoms[a].bonded[0]) - c).rescale(1.0), a2 = (pos(atoms[a].bonded[1]) - c).rescale(1.0);
      const Vec na = rotateAroundAxis(90.0 * M_PI / 180.0, a1 + a2, -a1);
      addProton(a, c + setBondDistance(na, atoms[a].element));
    }
    if (atoms[a].bonded.size() == 3 && atoms[a].protonsToAdd > 0) {
      const Vec a1 = (pos(atoms[a].bonded[0]) - c).rescale(1.0), a2 = (pos(atoms[a].bonded[1]) - c).rescale(1.0),
                a3 = (pos(atoms[a].bonded[2]) - c).rescale(1.0);
      addProton(a, c + setBondDistance(-a1 - a2 - a3, atoms[a].element));
    }
  }
  void protonateAtom(int a) {
    if (atoms[a].isProtonated || atoms[a].element == "H") return;
    setCharge(a);
    setProtonsToAdd(a);
    setStericNumber(a);
    if (atoms[a].stericNumber == 4) tetrahedral(a);
    else if (atoms[a].stericNumber == 3) trigonal(a);
    atoms[a].isProtonated = true;
  }

  // ----- groups (group.py) -----
  int newGroup(int a, const std::string& type, const std::string& residueTypeOverride = "") {
    Group g;
    g.atom = a;
    const PAtom& at = atoms[a];
    g.residueType = at.terminal.empty() ? at.resName : at.terminal;
    if (at.type == "atom") {
      g.label = std::format("{:<3}{:>4}{:>2}", g.residueType, at.resNum, at.chain);
    } else if (at.resName == "DA " || at.resName == "DC " || at.resName == "DG " || at.resName == "DT ") {
      std::string nm = at.name;
      nm.erase(std::remove(nm.begin(), nm.end(), '\''), nm.end());
      g.label = std::format("{}{}{}{:>4}{:>2}", g.residueType.substr(1, 1), at.element,
                            nm.empty() ? "" : nm.substr(nm.size() - 1), at.resNum, at.chain);
    } else {
      g.label = std::format("{:<3}{:>4}{:>2}", g.residueType, at.name, at.chain);
    }
    if (at.icode != " " && !at.icode.empty()) g.label += ":" + at.icode;
    g.type = type;
    if (!residueTypeOverride.empty()) g.residueType = residueTypeOverride;
    groups.push_back(std::move(g));
    atoms[a].group = static_cast<int>(groups.size()) - 1;
    return static_cast<int>(groups.size()) - 1;
  }
  void setCenter(Group& g, const std::vector<int>& v) {
    g.x = g.y = g.z = 0;
    for (int a : v) {
      g.x += atoms[a].x;
      g.y += atoms[a].y;
      g.z += atoms[a].z;
    }
    g.x /= static_cast<double>(v.size());
    g.y /= static_cast<double>(v.size());
    g.z /= static_cast<double>(v.size());
  }
  void setInteraction(Group& g, const std::vector<int>& acids, const std::vector<int>& bases) {
    for (int a : acids) atoms[a].groupType = g.type;
    for (int a : bases) atoms[a].groupType = g.type;
    g.acidAtoms = acids;
    g.baseAtoms = bases;
  }
  std::vector<int> concat(std::vector<int> a, const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end());
    return a;
  }
  void setupAtoms(int gi) {
    Group& g = groups[gi];
    const int a = g.atom;
    const std::string& t = g.type;
    auto hOf = [&](int x) { return bondedElements(x, "H"); };
    if (t == "COO") {
      if (atoms[a].terminal == "C-") {  // CtermGroup
        const auto cs = bondedElements(a, "C");
        if (cs.empty()) {
          setCenter(g, {a});
        } else {
          auto other = bondedElements(cs[0], "O");
          other.erase(std::find(other.begin(), other.end(), a));
          const std::vector<int> os = concat({a}, other);
          setCenter(g, os);
          setInteraction(g, os, os);
        }
        return;
      }
      const auto os = bondedElements(a, "O");
      setCenter(g, os.empty() ? std::vector<int>{a} : os);
      setInteraction(g, os, os);
    } else if (t == "HIS") {
      const auto ring = ringMember(a);
      if (ring.size() != 5) logWarn(std::format("PROPKA: His group does not seem to contain a ring ({})", g.label));
      for (int r : ring) protonateAtom(r);
      setCenter(g, ring.empty() ? std::vector<int>{a} : ring);
      std::vector<int> hs, ns;
      for (int r : ring)
        if (atoms[r].element == "N") ns.push_back(r);
      for (int n : ns) {
        const auto h = hOf(n);
        hs.insert(hs.end(), h.begin(), h.end());
      }
      setInteraction(g, concat(hs, ns), ns);
    } else if (t == "ARG") {
      setCenter(g, {a});
      const auto ns = bondedElements(a, "N");
      for (int n : ns) protonateAtom(n);
      std::vector<int> hs;
      for (int n : ns) {
        const auto h = hOf(n);
        hs.insert(hs.end(), h.begin(), h.end());
      }
      setInteraction(g, concat(ns, hs), ns);
    } else if (t == "AMD") {
      const auto o = bondedElements(a, "O"), n = bondedElements(a, "N");
      if (o.empty() || n.empty()) {
        logWarn(std::format("PROPKA: missing N or O atom ({})", g.label));
        setCenter(g, {a});
        return;
      }
      protonateAtom(n[0]);
      const auto hs = hOf(n[0]);
      setCenter(g, concat(o, n));
      setInteraction(g, concat(n, hs), o);
    } else if (t == "TRP") {
      setCenter(g, {a});
      protonateAtom(a);
      setInteraction(g, concat(hOf(a), {a}), {a});
    } else if (t == "BBN" || t == "NAR" || t == "NAM" || t == "OH") {
      protonateAtom(a);
      const auto h = hOf(a);
      setCenter(g, {a});
      setInteraction(g, concat(h, {a}), concat(h, {a}));
    } else if (t == "BBC") {
      const auto o = bondedElements(a, "O");
      setCenter(g, {a});
      setInteraction(g, o, o);
    } else if (t == "OP") {
      protonateAtom(a);
      setCenter(g, {a});
      setInteraction(g, {a}, {a});
    } else if (t == "CG") {
      const auto ns = bondedElements(a, "N");
      setCenter(g, {a});
      std::vector<int> hs;
      for (int n : ns) {
        protonateAtom(n);
        const auto h = hOf(n);
        hs.insert(hs.end(), h.begin(), h.end());
      }
      setInteraction(g, concat(hs, ns), ns);
    } else if (t == "C2N") {
      std::vector<int> ns;
      for (int n : bondedElements(a, "N"))
        if (bondedHeavy(n).size() == 1) ns.push_back(n);
      setCenter(g, {a});
      std::vector<int> hs;
      for (int n : ns) {
        protonateAtom(n);
        const auto h = hOf(n);
        hs.insert(hs.end(), h.begin(), h.end());
      }
      setInteraction(g, concat(hs, ns), ns);
    } else if (t == "OCO") {
      const auto os = bondedElements(a, "O");
      setCenter(g, os);
      setInteraction(g, os, os);
    } else if (t == "N30" || t == "N31" || t == "N32" || t == "N33" || t == "NP1") {
      protonateAtom(a);
      const auto h = hOf(a);
      setCenter(g, {a});
      setInteraction(g, concat(h, {a}), {a});
    } else {  // Group.setup_atoms
      setCenter(g, {a});
      setInteraction(g, {a}, {a});
    }
  }
  void setupGroup(int gi) {
    Group& g = groups[gi];
    if (auto it = P.charge.find(g.type); it != P.charge.end()) g.charge = it->second;
    if (auto it = P.ions.find(g.residueType); it != P.ions.end()) g.charge = it->second;
    setupAtoms(gi);
    Group& h = groups[gi];
    h.titratable = false;
    if (auto it = P.modelPkas.find(h.residueType); it != P.modelPkas.end()) {
      if (!h.modelPkaSet) {
        h.modelPka = it->second;
        std::string rn = atoms[h.atom].resName;
        rn.erase(rn.find_last_not_of(' ') + 1);
        if (auto c = P.customModelPkas.find(rn + "-" + atoms[h.atom].name); c != P.customModelPkas.end()) h.modelPka = c->second;
        h.modelPkaSet = true;
      }
    }
    if (h.modelPkaSet && !atoms[h.atom].cysteineBridge) h.titratable = true;
  }
  std::optional<int> isProteinGroup(int a) {
    const PAtom& at = atoms[a];
    if (at.type != "atom") return std::nullopt;
    if (at.terminal == "N+") return newGroup(a, "N+");
    if (at.terminal == "C-") return newGroup(a, "COO");
    if (at.name == "N" && at.resName != "PRO") return newGroup(a, "BBN", "BBN");
    if (at.name == "C" && bondedElements(a, "O").size() == 1) return newGroup(a, "BBC", "BBC");
    if (auto it = P.proteinGroupMapping.find(at.resName + "-" + at.name); it != P.proteinGroupMapping.end())
      return newGroup(a, it->second);
    return std::nullopt;
  }
  std::optional<int> isLigandGroup(int a) {
    if (atoms[a].type != "hetatm") return std::nullopt;
    protonateAtom(a);
    const PAtom& at = atoms[a];
    const std::string& s = at.sybylType;
    if (s == "N.ar" && bondedHeavy(a).size() == 2) return newGroup(a, "NAR", "NAR");
    if (s == "N.am") return newGroup(a, "NAM", "NAM");
    if (s == "N.3" || s == "N.4") {
      static const char* t[] = {"N30", "N31", "N32", "N33"};
      const std::size_t h = bondedHeavy(a).size();
      if (h <= 3) return newGroup(a, t[h], t[h]);
    }
    if (s == "N.1") return newGroup(a, "N1", "N1");
    if (s == "N.pl3") {
      const auto cs = bondedElements(a, "C");
      if (cs.size() == 1 && bondedElements(cs[0], "N").size() == 1) return newGroup(a, "NP1", "NP1");
    }
    if (s == "C.2") {
      const auto ns = bondedElements(a, "N");
      int npls = 0, maxTwo = 0;
      for (int n : ns)
        if (atoms[n].sybylType == "N.pl3" && bondedHeavy(n).size() == 1) ++npls;
      if (npls == 2) {
        for (int n : ns)
          if (bondedHeavy(n).size() < 3) ++maxTwo;
        if (maxTwo == 2) return newGroup(a, "C2N", "C2N");
        if (ns.size() == 3) return newGroup(a, "CG", "CG");
      }
      int co2 = 0;
      for (int o : bondedElements(a, "O"))
        if (atoms[o].sybylType.find("O.co2") != std::string::npos) ++co2;
      if (co2 == 2) return newGroup(a, "OCO", "OCO");
    }
    if (s == "F") return newGroup(a, "F", "F");
    if (s == "Cl") return newGroup(a, "Cl", "Cl");
    if (s == "O.3") {
      if (bondedHeavy(a).size() == 1) {
        if (bondedElements(a, "P").size() == 1) return newGroup(a, "OP", "OP");
        return newGroup(a, "OH", "OH");
      }
      return newGroup(a, "O3", "O3");
    }
    if (s == "O.2") return newGroup(a, "O2", "O2");
    if (s == "S.3" && bondedHeavy(a).size() == 1) return newGroup(a, "SH", "SH");
    return std::nullopt;
  }
  void extractGroups(const std::vector<int>& heavy) {
    for (int a : heavy) {
      std::optional<int> g = isProteinGroup(a);
      if (!g) {
        std::string rn = atoms[a].resName;
        rn.erase(rn.find_last_not_of(' ') + 1);
        if (P.ions.contains(rn)) {
          g = newGroup(a, "ION", rn);
        }
      }
      if (!g) g = isLigandGroup(a);
      if (g) setupGroup(*g);
    }
  }

  // ----- calculations (energy.py, determinants.py, iterative.py) -----
  double pairWeight(int n1, int n2) const {
    const double w = static_cast<double>(n1 + n2 - 2 * P.nmin) / static_cast<double>(2 * P.nmax - 2 * P.nmin);
    return std::max(0.0, std::min(1.0, w));
  }
  void desolvation(Group& g, const std::vector<int>& heavySorted) {
    double volume = 0;
    g.numVolume = 0;
    const double minDist4 = std::pow(2.75, 4);
    const PAtom& ga = atoms[g.atom];
    for (int i : heavySorted) {
      const PAtom& a = atoms[i];
      if (a.resNum == ga.resNum && a.chain == ga.chain && a.icode == ga.icode) continue;
      const double sq = sqDist(g.x, g.y, g.z, a.x, a.y, a.z);
      if (sq < P.desolvCutoff * P.desolvCutoff) {
        double dvol;
        if (a.element == "C" && a.name != "CA" && a.name != "C") {
          dvol = P.vdwVolume.at("C4");
        } else {
          auto it = P.vdwVolume.find(a.element);
          dvol = it == P.vdwVolume.end() ? 1.0 : it->second;
        }
        volume += dvol / std::max(minDist4, sq * sq);
      }
      if (sq < P.buriedCutoff * P.buriedCutoff) g.numVolume += 1;
    }
    double w = static_cast<double>(g.numVolume - P.nmin) / static_cast<double>(P.nmax - P.nmin);
    g.buried = std::max(0.0, std::min(1.0, w));
    const double scale = 1.0 - (1.0 - P.surfScale) * (1.0 - g.buried);
    g.energyVolume = g.charge * P.desolvPrefactor * std::max(0.0, volume - P.desolvAllowance) * scale;
  }
  static double hbondEnergy(double dist, double dpkaMax, std::pair<double, double> c, double fAngle) {
    double v;
    if (dist < c.first) v = 1.0;
    else if (dist > c.second) v = 0.0;
    else v = 1.0 - (dist - c.first) / (c.second - c.first);
    return std::abs(dpkaMax * v * fAngle);
  }
  // angle_distance_factors(atom1 or center, atom2, atom3) -> dist12, f_angle
  std::pair<double, double> angleDistance(const Vec& p1, const Vec& p2, const Vec& p3) const {
    double dx32 = p2.x - p3.x, dy32 = p2.y - p3.y, dz32 = p2.z - p3.z;
    const double d23 = std::sqrt(dx32 * dx32 + dy32 * dy32 + dz32 * dz32);
    dx32 /= d23; dy32 /= d23; dz32 /= d23;
    double dx21 = p1.x - p2.x, dy21 = p1.y - p2.y, dz21 = p1.z - p2.z;
    const double d12 = std::sqrt(dx21 * dx21 + dy21 * dy21 + dz21 * dz21);
    dx21 /= d12; dy21 /= d12; dz21 /= d12;
    return {d12, dx21 * dx32 + dy21 * dy32 + dz21 * dz32};
  }
  struct Closest {
    int a = -1, b = -1;
    double dist = 1e6;
  };
  Closest smallestDistance(const std::vector<int>& v1, const std::vector<int>& v2) const {
    double best = 1e6;
    Closest c;
    for (int a : v1)
      for (int b : v2) {
        const double d = sqDist(atoms[a].x, atoms[a].y, atoms[a].z, atoms[b].x, atoms[b].y, atoms[b].z);
        if (d < best) {
          best = d;
          c.a = a;
          c.b = b;
        }
      }
    c.dist = std::sqrt(best);
    return c;
  }
  const std::vector<int>& interactionAtoms(int g, int other) const {
    return P.inList(P.baseList, groups[other].residueType) ? groups[g].baseAtoms : groups[g].acidAtoms;
  }
  std::pair<double, std::pair<double, double>> hbParams(int a1, int a2) const {
    return {P.sidechainInteraction, P.cutoff(atoms[a1].groupType, atoms[a2].groupType)};
  }
  static bool checkBuried(int n1, int n2) { return !(n1 + n2 <= 900 && (n1 <= 400 || n2 <= 400)); }
  // check_exceptions: {exception, value}
  std::pair<bool, double> exceptions(int g1, int g2) const {
    const std::string &t1 = groups[g1].type, &t2 = groups[g2].type;
    auto cooArg = [&](int gc, int ga) {
      double tot = 0;
      std::vector<int> ac = interactionAtoms(gc, ga), aa = interactionAtoms(ga, gc);
      for (int k = 0; k < 2; ++k) {
        const Closest c = smallestDistance(ac, aa);
        if (c.a < 0 || c.b < 0) continue;
        const auto [dmax, cut] = hbParams(c.a, c.b);
        double f = 1.0, dist = c.dist;
        if (P.inList(P.angular, groups[ga].type)) {
          const int a3 = atoms[c.b].bonded[0];
          std::tie(dist, f) = angleDistance(pos(c.a), pos(c.b), pos(a3));
        }
        tot += hbondEnergy(dist, dmax, cut, f);
        ac.erase(std::find(ac.begin(), ac.end(), c.a));
        aa.erase(std::find(aa.begin(), aa.end(), c.b));
      }
      return std::pair<bool, double>{true, tot};
    };
    if (t1 == "COO" && t2 == "ARG") return cooArg(g1, g2);
    if (t1 == "ARG" && t2 == "COO") return cooArg(g2, g1);
    if (t1 == "COO" && t2 == "COO") {
      const Closest c = smallestDistance(interactionAtoms(g1, g2), interactionAtoms(g2, g1));
      const auto [dmax, cut] = hbParams(c.a, c.b);
      double v = hbondEnergy(c.dist, dmax, cut, 1.0);
      v *= 1.0 + pairWeight(groups[g1].numVolume, groups[g2].numVolume);
      return {true, v};
    }
    const bool b = checkBuried(groups[g1].numVolume, groups[g2].numVolume);
    if (t1 == "CYS" && t2 == "CYS") return {b, P.cysCys};
    if ((t1 == "COO" && t2 == "HIS") || (t1 == "HIS" && t2 == "COO")) return {b, P.cooHis};
    if ((t1 == "OCO" && t2 == "HIS") || (t1 == "HIS" && t2 == "OCO")) return {b, P.ocoHis};
    if ((t1 == "CYS" && t2 == "HIS") || (t1 == "HIS" && t2 == "CYS")) return {b, P.cysHis};
    return {false, 0.0};
  }
  std::optional<double> hbondInteraction(int g1, int g2) const {
    const Closest c = smallestDistance(interactionAtoms(g1, g2), interactionAtoms(g2, g1));
    if (c.a < 0 || c.b < 0) {
      logWarn(std::format("PROPKA: side chain interaction failed for {} and {}", groups[g1].label, groups[g2].label));
      return std::nullopt;
    }
    const auto [dmax, cut] = hbParams(c.a, c.b);
    double dist = c.dist;
    if (dist >= cut.second) return std::nullopt;
    if (withinBondDistance(groups[g1].atom, groups[g2].atom, P.minBondDistHb, 1)) return std::nullopt;
    double f = 1.0;
    if (P.inList(P.angular, groups[g2].type)) {
      if (atoms[c.b].element == "H") std::tie(dist, f) = angleDistance(pos(c.a), pos(c.b), pos(atoms[c.b].bonded[0]));
      else f = 0.0;
    } else if (P.inList(P.angular, groups[g1].type)) {
      if (atoms[c.a].element == "H") std::tie(dist, f) = angleDistance(pos(c.b), pos(c.a), pos(atoms[c.a].bonded[0]));
      else f = 0.0;
    }
    const auto [exc, ev] = exceptions(g1, g2);
    return exc ? ev : hbondEnergy(dist, dmax, cut, f);
  }
  double coulombEnergy(double dist, double weight) const {
    const double diel = 160.0 - (160.0 - 30.0) * weight;
    dist = std::max(dist, P.coulombCutoff1);
    double scale = (dist - P.coulombCutoff2) / (P.coulombCutoff1 - P.coulombCutoff2);
    scale = std::max(0.0, std::min(1.0, scale));
    return std::abs(244.12 / (diel * dist) * scale);
  }
  std::optional<double> electrostatic(int g1, int g2, double dist) const {
    const Group &a = groups[g1], &b = groups[g2];
    if (!(a.titratable && b.titratable) || dist > P.coulombCutoff2 || a.numVolume + b.numVolume < P.nmin) return std::nullopt;
    return coulombEnergy(dist, pairWeight(a.numVolume, b.numVolume));
  }
  static bool truthy(const std::optional<double>& v) { return v && *v != 0.0; }
  void addDet(int g, DetType t, int other, double v) { groups[g].dets[t].push_back({other, groups[other].label, v}); }

  void addSidechainDeterminants(int g1, int g2) {
    const auto hb = hbondInteraction(g1, g2);
    if (!truthy(hb)) return;
    const double h = *hb;
    const Group &a = groups[g1], &b = groups[g2];
    double v1, v2;  // determinant on g1 (from g2), on g2 (from g1)
    if (a.charge == b.charge) {
      if (a.modelPka < b.modelPka) { v1 = -h; v2 = h; }
      else { v1 = h; v2 = -h; }
    } else {
      v1 = h * a.charge;
      v2 = h * b.charge;
    }
    addDet(g1, kSidechain, g2, v1);
    addDet(g2, kSidechain, g1, v2);
  }
  void addCoulombDeterminants(int g1, int g2, double dist) {
    const auto c = electrostatic(g1, g2, dist);
    if (!truthy(c)) return;
    const double v = *c;
    const double q1 = groups[g1].charge, q2 = groups[g2].charge;
    if (q1 < 0 && q2 < 0) {
      if (groups[g1].modelPka > groups[g2].modelPka) addDet(g1, kCoulomb, g2, v);
      else addDet(g2, kCoulomb, g1, v);
    } else if (q1 > 0 && q2 > 0) {
      if (groups[g1].modelPka < groups[g2].modelPka) addDet(g1, kCoulomb, g2, -v);
      else addDet(g2, kCoulomb, g1, -v);
    } else {
      addDet(g1, kCoulomb, g2, q1 * v);
      addDet(g2, kCoulomb, g1, q2 * v);
    }
  }

  struct Interaction {
    int g1, g2;
    double hb, coul;
    double ann[2] = {0, 0};
  };
  struct Iter {
    int group;
    double q, pkaOld, pkaNew, pkaNon;
    std::vector<std::pair<int, double>> dets[3];
  };
  void iterativeDeterminants(std::vector<Interaction>& inter) {
    std::vector<Iter> its;
    std::vector<int> done;
    for (const Interaction& in : inter)
      for (int g : {in.g1, in.g2}) {
        if (groupIn(g, done)) continue;
        Iter it{};
        it.group = g;
        it.q = groups[g].charge;
        double sc = 0, bb = 0, co = 0;
        for (const Det& d : groups[g].dets[kSidechain]) sc += d.value;
        for (const Det& d : groups[g].dets[kBackbone]) bb += d.value;
        for (const Det& d : groups[g].dets[kCoulomb]) co += d.value;
        it.pkaNon = groups[g].modelPka;
        it.pkaNon += groups[g].energyVolume;
        it.pkaNon += groups[g].energyLocal;
        it.pkaNon += sc;
        it.pkaNon += bb;
        it.pkaNon += co;
        it.pkaOld = it.pkaNon;
        its.push_back(std::move(it));
        done.push_back(g);
      }
    auto find = [&](int g) -> Iter& {
      for (Iter& it : its)
        if (groupEq(it.group, g)) return it;
      fail("PROPKA: iterative group not found");
    };
    for (int iteration = 1;; ++iteration) {
      for (Iter& it : its) {
        for (auto& d : it.dets) d.clear();
        it.pkaNew = it.pkaNon;
      }
      for (Interaction& in : inter) {
        Iter& o1 = find(in.g1);
        Iter& o2 = find(in.g2);
        const double hb = in.hb, co = in.coul;
        if (o1.q < 0 && o2.q < 0) {
          const double diff = co + 2 * hb;
          const double c1 = o1.pkaOld + in.ann[0] + diff, c2 = o2.pkaOld + in.ann[1] + diff;
          in.ann[0] = in.ann[1] = 0;
          if (c1 > c2) {
            o1.dets[kSidechain].push_back({o2.group, hb});
            o2.dets[kSidechain].push_back({o1.group, -hb});
            o1.dets[kCoulomb].push_back({o2.group, co});
            in.ann[0] = -diff;
          } else {
            o2.dets[kSidechain].push_back({o1.group, hb});
            o1.dets[kSidechain].push_back({o2.group, -hb});
            o2.dets[kCoulomb].push_back({o1.group, co});
            in.ann[1] = -diff;
          }
        } else if (o1.q > 0 && o2.q > 0) {
          const double diff = -(co + 2 * hb);
          const double c1 = o1.pkaOld + in.ann[0] + diff, c2 = o2.pkaOld + in.ann[1] + diff;
          in.ann[0] = in.ann[1] = 0;
          if (c1 < c2) {
            o1.dets[kSidechain].push_back({o2.group, -hb});
            o2.dets[kSidechain].push_back({o1.group, hb});
            o1.dets[kCoulomb].push_back({o2.group, -co});
            in.ann[0] = -diff;
          } else {
            o2.dets[kSidechain].push_back({o1.group, -hb});
            o1.dets[kSidechain].push_back({o2.group, hb});
            o2.dets[kCoulomb].push_back({o1.group, -co});
            in.ann[1] = -diff;
          }
        } else {
          const double q1 = o1.q, q2 = o2.q;
          double c1 = o1.pkaOld + in.ann[0] + q1 * co, c2 = o2.pkaOld + in.ann[1] + q2 * co;
          c1 += q1 * hb;  // exclude_sidechain_interactions is empty
          c2 += q2 * hb;
          const bool add = (q1 == -1.0 && c1 < c2) || (q1 == 1.0 && c1 > c2);
          in.ann[0] = in.ann[1] = 0;
          if (add) {
            if (co > 0.005) {
              o1.dets[kCoulomb].push_back({o2.group, q1 * co});
              in.ann[0] += -q1 * co;
              o2.dets[kCoulomb].push_back({o1.group, q2 * co});
              in.ann[1] += -q2 * co;
            }
            if (hb > 0.005) {
              o1.dets[kSidechain].push_back({o2.group, q1 * hb});
              in.ann[0] += -q1 * hb;
              o2.dets[kSidechain].push_back({o1.group, q2 * hb});
              in.ann[1] += -q2 * hb;
            }
          }
        }
      }
      for (Iter& it : its)
        for (auto& dl : it.dets)
          for (auto& d : dl) it.pkaNew += d.second;
      bool converged = true;
      for (Iter& it : its)
        if (it.pkaNew != it.pkaOld) converged = false;
      for (Iter& it : its) it.pkaOld = it.pkaNew;
      if (converged || iteration == 10) break;
    }
    for (Iter& it : its)
      for (int t = 0; t < 3; ++t)
        for (const auto& [g, v] : it.dets[t])
          if (v > 0.005 || v < -0.005) addDet(it.group, static_cast<DetType>(t), g, v);
  }

  void setDeterminants(const std::vector<int>& sc) {
    std::vector<Interaction> inter;
    for (int g1 : sc)
      for (int g2 : sc) {
        if (groupEq(g1, g2)) break;
        if (groupIn(g2, groups[g1].covCoupled)) break;
        const double dist = std::sqrt(sqDist(groups[g1].x, groups[g1].y, groups[g1].z, groups[g2].x, groups[g2].y, groups[g2].z));
        if (dist < P.coulombCutoff2) {
          const char t = P.matrixValue(groups[g1].type, groups[g2].type);
          if (t == 'I') {
            const auto hb = hbondInteraction(g1, g2);
            const auto co = electrostatic(g1, g2, dist);
            if (truthy(hb) || truthy(co)) inter.push_back({g1, g2, hb.value_or(0.0), co.value_or(0.0)});
          } else if (t == 'N') {
            addSidechainDeterminants(g1, g2);
            addCoulombDeterminants(g1, g2, dist);
          }
        }
      }
    iterativeDeterminants(inter);
  }

  void backboneDeterminants(const std::vector<int>& tit, const std::vector<int>& bb) {
    for (int tg : tit) {
      const std::vector<int>& ta = groups[tg].acidAtoms;
      if (ta.empty()) continue;
      for (int bg : bb) {
        const std::vector<int>& ba = interactionAtoms(bg, tg);
        if (ba.empty()) continue;
        const Closest c = smallestDistance(ba, ta);
        const int batom = c.a, tatom = c.b;
        const std::map<std::string, std::vector<double>>* tab = nullptr;
        if (atoms[batom].groupType == "BBC" && P.bbCO.contains(atoms[tatom].groupType)) tab = &P.bbCO;
        else if (atoms[batom].groupType == "BBN" && P.bbNH.contains(atoms[tatom].groupType)) tab = &P.bbNH;
        if (!tab) continue;
        const std::vector<double>& prm = tab->at(atoms[tatom].groupType);
        const double dmax = prm[0];
        const std::pair<double, double> cut{prm[1], prm[2]};
        if (c.dist < cut.second) {
          double f = 1.0;
          if (groups[bg].type == "BBC" && P.inList(P.angular, groups[tg].type)) {
            if (atoms[tatom].element == "H") f = angleDistance(pos(atoms[tatom].bonded[0]), pos(tatom), pos(batom)).second;
            else f = 0.0;
          }
          if (groups[bg].type == "BBN") {
            if (atoms[batom].element == "H") f = angleDistance(pos(tatom), pos(batom), pos(atoms[batom].bonded[0])).second;
            else f = 0.0;
          }
          if (f > 0.001) addDet(tg, kBackbone, bg, groups[tg].charge * hbondEnergy(c.dist, dmax, cut, f));
        }
      }
    }
  }

  std::vector<int> bondedTitratable(int a, int nb, int orig) const {
    std::vector<int> res;
    for (int b : atoms[a].bonded) {
      if (b == orig) continue;
      const int g = atoms[b].group;
      if (g >= 0 && groups[g].titratable && nb <= P.couplingMaxBonds && std::find(res.begin(), res.end(), g) == res.end())
        res.push_back(g);
      if (nb < P.couplingMaxBonds)
        for (int x : bondedTitratable(b, nb + 1, orig))
          if (std::find(res.begin(), res.end(), x) == res.end()) res.push_back(x);
    }
    return res;
  }
  void coupleCovalently() {
    for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
      if (!groups[g].titratable) continue;
      for (int b : bondedTitratable(groups[g].atom, 1, groups[g].atom)) {
        if (groupIn(b, groups[g].covCoupled)) continue;
        if (atoms[groups[b].atom].sybylType == atoms[groups[g].atom].sybylType) {
          if (!groupIn(b, groups[g].covCoupled)) groups[g].covCoupled.push_back(b);
          if (!groupIn(g, groups[b].covCoupled)) groups[b].covCoupled.push_back(g);
        }
      }
    }
  }
  void totalPka(Group& g) {
    if (atoms[g.atom].cysteineBridge) {
      g.pka = 99.99;
      return;
    }
    g.pka = g.modelPka + g.energyVolume + g.energyLocal;
    for (const auto& dl : g.dets)
      for (const Det& d : dl) g.pka += d.value;
  }
};

}  // namespace

std::vector<GroupResult> calculate(const std::vector<InputAtom>& input) {
  Mol m;
  const Params& P = m.P;
  // Reading (input.py read_mmcif): hydrogens and ignored residues dropped, terminals marked.
  std::map<std::string, std::string> firstResidue;
  for (std::size_t i = 0; i < input.size(); ++i) {
    const InputAtom& in = input[i];
    if (in.element == "H" || P.inList(P.ignoreResidues, in.resName)) continue;
    PAtom a;
    a.name = in.name;
    a.resName = std::format("{:<3}", in.resName);
    a.chain = in.chain;
    a.resNum = in.resNum;
    a.icode = in.iCode.empty() ? " " : in.iCode;
    a.element = in.element;
    a.type = in.hetatm ? "hetatm" : "atom";
    if (a.resName == "DA " || a.resName == "DC " || a.resName == "DG " || a.resName == "DT ") a.type = "hetatm";
    a.x = in.xyz[0];
    a.y = in.xyz[1];
    a.z = in.xyz[2];
    a.input = static_cast<int>(i);
    if (!in.hetatm) {
      const std::string rk = std::format("{}|{}|{}", a.chain, a.resNum, a.icode);
      auto [it, added] = firstResidue.try_emplace(a.chain, rk);
      if (a.name == "N" && it->second == rk) a.terminal = "N+";
      if (a.name == "OXT" || a.name == "O''") a.terminal = "C-";
    }
    m.atoms.push_back(std::move(a));
  }
  const int nHeavy = static_cast<int>(m.atoms.size());
  // setup_bonding_and_protonation (VersionA).
  m.findBondsUsingBoxes();
  for (int i = 0; i < nHeavy; ++i)
    if (m.atoms[i].type == "hetatm") m.assignSybylType(i);
  m.addPiElectronInfo();
  std::vector<int> heavy(nHeavy);
  for (int i = 0; i < nHeavy; ++i) heavy[i] = i;
  m.extractGroups(heavy);
  // sort_atoms: (chain, res_num, icode, element order), stable.
  std::vector<int> order(m.atoms.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
  auto key = [&](int i) {
    const PAtom& a = m.atoms[i];
    std::string eo;
    if (a.name.size() > a.element.size()) eo = a.name.substr(a.element.size(), 1);
    return std::tuple<const std::string&, int, const std::string&, std::string>(a.chain, a.resNum, a.icode, eo);
  };
  std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return key(x) < key(y); });
  std::vector<int> heavySorted;
  for (int i : order)
    if (m.atoms[i].element != "H") heavySorted.push_back(i);
  m.coupleCovalently();

  // calculate_pka
  std::vector<int> tit, ions, bb, sc;
  for (int g = 0; g < static_cast<int>(m.groups.size()); ++g) {
    const Group& G = m.groups[g];
    if (G.titratable) tit.push_back(g);
    if (P.ions.contains(G.residueType)) ions.push_back(g);
    if (G.type.find("BB") != std::string::npos) bb.push_back(g);
    else if (!m.atoms[G.atom].cysteineBridge) sc.push_back(g);
  }
  for (int g : tit) m.desolvation(m.groups[g], heavySorted);
  for (int g : ions)
    if (!m.groups[g].titratable) m.desolvation(m.groups[g], heavySorted);
  m.backboneDeterminants(tit, bb);
  for (int tg : tit)
    for (int ig : ions) {
      const Group &T = m.groups[tg], &I = m.groups[ig];
      const double d2 = sqDist(T.x, T.y, T.z, I.x, I.y, I.z);
      if (d2 < P.coulombCutoff2 * P.coulombCutoff2) {
        const double w = m.pairWeight(T.numVolume, I.numVolume);
        m.addDet(tg, kCoulomb, ig, -I.charge * m.coulombEnergy(std::sqrt(d2), w));
      }
    }
  // backbone_reorganization
  for (int g = 0; g < static_cast<int>(m.groups.size()); ++g) {
    Group& T = m.groups[g];
    if (!P.inList(P.bbReorg, T.residueType) || m.atoms[T.atom].cysteineBridge) continue;
    double dpka = 0;
    for (int b : bb) {
      if (m.groups[b].type != "BBC") continue;
      const int o = m.interactionAtoms(b, g)[0];
      const auto [dist, f] = m.angleDistance({T.x, T.y, T.z}, m.pos(o), m.pos(m.groups[b].atom));
      if (dist < 6.0 && f > 0.001) dpka += 0.8 * std::min(1.0, 1.0 - (dist - 3.0) / (6.0 - 3.0));
    }
    T.energyLocal = dpka * T.buried;
  }
  m.setDeterminants(sc);
  for (Group& g : m.groups) m.totalPka(g);
  // coupling_effects
  std::vector<std::string> penalised;
  {
    std::vector<int> coupled;
    for (int g = 0; g < static_cast<int>(m.groups.size()); ++g)
      if (!m.groups[g].covCoupled.empty()) coupled.push_back(g);
    std::vector<bool> used(m.groups.size(), false);
    for (int start : coupled) {
      if (used[start]) continue;
      std::vector<int> sys, stack{start};
      used[start] = true;
      while (!stack.empty()) {
        const int g = stack.back();
        stack.pop_back();
        sys.push_back(g);
        for (int c : m.groups[g].covCoupled)
          if (!used[c]) {
            used[c] = true;
            stack.push_back(c);
          }
      }
      int first = sys[0];
      for (int g : sys)
        if (m.groups[g].pka > m.groups[first].pka) first = g;
      if (m.groups[first].charge < 0) {
        int mn = sys[0];
        for (int g : sys)
          if (m.groups[g].pka < m.groups[mn].pka) mn = g;
        m.groups[first].coupledTitrating = mn;
        penalised.push_back(m.groups[first].label);
      } else {
        for (int g : sys) {
          if (g == first) continue;
          m.groups[g].coupledTitrating = first;
          penalised.push_back(m.groups[g].label);
        }
      }
    }
  }
  if (P.remove_penalised && !penalised.empty()) {
    for (int g : tit)
      for (auto& dl : m.groups[g].dets)
        std::erase_if(dl, [&](const Det& d) { return std::find(penalised.begin(), penalised.end(), d.label) != penalised.end(); });
    for (Group& g : m.groups) m.totalPka(g);
  }

  std::vector<GroupResult> out;
  for (int g : tit) {
    const Group& G = m.groups[g];
    out.push_back({m.atoms[G.atom].input, G.type, G.residueType, G.label, G.pka, G.modelPka, G.buried, G.coupledTitrating >= 0});
  }
  return out;
}

}  // namespace od::propka
