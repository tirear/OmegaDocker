// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// The DSSP algorithm (Kabsch & Sander, Biopolymers 22:2577, 1983), with the
// behavior of DSSP 4 (https://github.com/PDB-REDO/dssp):
//   Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute. All rights reserved.
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are met:
//   1. Redistributions of source code must retain the above copyright notice, this
//      list of conditions and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above copyright notice,
//      this list of conditions and the following disclaimer in the documentation
//      and/or other materials provided with the distribution.
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
//   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "dssp.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <numbers>

namespace od {
namespace {

constexpr double kCoupling = -27.888;  // -332 * 0.42 * 0.20 kcal/mol
constexpr double kMinEnergy = -9.9;
constexpr double kMaxHBondEnergy = -0.5;
constexpr double kMinDistance = 0.5;
constexpr double kMaxCADistance = 9.0;
constexpr double kBendAngle = 70.0;
constexpr double kPPPhi = -75.0, kPPPsi = 145.0, kPPRange = 29.0;
constexpr int kPPMinLength = 3;

// Single precision as in DSSP 4, so that energies near the thresholds round alike.
struct V {
  float x = 0, y = 0, z = 0;
};
V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
V operator/(V a, float s) { return {a.x / s, a.y / s, a.z / s}; }
float dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
float dist(V a, V b) { return std::sqrt(dot(a - b, a - b)); }
V toV(const Vec3& v) { return {static_cast<float>(v[0]), static_cast<float>(v[1]), static_cast<float>(v[2])}; }

struct HBond {
  int partner = -1;
  double energy = 0;
};

enum class HelixFlag : std::uint8_t { None, Start, Middle, End, StartAndEnd };

struct Res {
  int mol = -1, residue = -1;
  int chain = -1;            // chain index for ordering ladders
  V n, ca, c, o, h;
  bool donor = false;        // N bonded to no heavy atom besides CA and C(i-1) (excludes proline)
  bool breakBefore = false;  // chain break between this and the previous residue
  double phi = 360, psi = 360;
  HBond acceptor[2];         // best acceptors of this N-H
  HBond donorOf[2];          // best donors to this C=O
  HelixFlag flag[4] = {};    // strides 3, 4, 5 and PP (index 3)
  char ss = '.';
};

void keepBest(HBond (&b)[2], int partner, double e) {
  if (e < b[0].energy) {
    b[1] = b[0];
    b[0] = {partner, e};
  } else if (e < b[1].energy) {
    b[1] = {partner, e};
  }
}

// Energy of N-H(d) ... O=C(a).
double hbondEnergy(const Res& d, const Res& a) {
  if (!d.donor) return 0;
  const double dHO = dist(d.h, a.o), dHC = dist(d.h, a.c), dNC = dist(d.n, a.c), dNO = dist(d.n, a.o);
  double e;
  if (dHO < kMinDistance || dHC < kMinDistance || dNC < kMinDistance || dNO < kMinDistance)
    e = kMinEnergy;
  else
    e = kCoupling / dHO - kCoupling / dHC + kCoupling / dNC - kCoupling / dNO;
  e = std::round(e * 1000) / 1000;
  return std::max(e, kMinEnergy);
}

double dihedral(V a, V b, V c, V d) {
  Vec3 p{a.x, a.y, a.z}, q{b.x, b.y, b.z}, r{c.x, c.y, c.z}, s{d.x, d.y, d.z};
  return dihedralDeg(p, q, r, s);
}

class Dssp {
 public:
  explicit Dssp(std::vector<Res> r) : r_(std::move(r)) {}

  void run() {
    hbonds();
    sheets();
    helices();
    ppHelices();
  }
  const std::vector<Res>& residues() const { return r_; }

 private:
  std::vector<Res> r_;
  int n() const { return static_cast<int>(r_.size()); }

  bool noBreak(int a, int b) const {
    if (a > b) std::swap(a, b);
    for (int i = a + 1; i <= b; ++i)
      if (r_[i].breakBefore) return false;
    return true;
  }
  // N-H(d) ... O=C(a)
  bool bond(int d, int a) const {
    for (const HBond& h : r_[d].acceptor)
      if (h.partner == a && h.energy < kMaxHBondEnergy) return true;
    return false;
  }

  // Residue pairs with CA-CA distance <= kMaxCADistance, i < j.
  std::vector<std::pair<int, int>> near_;

  void hbonds() {
    for (int i = 0; i < n(); ++i)
      for (int j = i + 1; j < n(); ++j) {
        const V d = r_[i].ca - r_[j].ca;
        if (dot(d, d) > static_cast<float>(kMaxCADistance * kMaxCADistance)) continue;
        near_.emplace_back(i, j);
        double e = hbondEnergy(r_[i], r_[j]);
        if (e < 0) {
          keepBest(r_[i].acceptor, j, e);
          keepBest(r_[j].donorOf, i, e);
        }
        if (j != i + 1) {
          e = hbondEnergy(r_[j], r_[i]);
          if (e < 0) {
            keepBest(r_[j].acceptor, i, e);
            keepBest(r_[i].donorOf, j, e);
          }
        }
      }
  }

  enum class BridgeType { None, Parallel, Antiparallel };
  struct Ladder {
    BridgeType type;
    int chain;
    std::deque<int> i, j;
  };

  BridgeType bridge(int i, int j) const {
    const int a = i - 1, b = i, c = i + 1, d = j - 1, e = j, f = j + 1;
    if (a < 0 || f >= n() || !noBreak(a, c) || !noBreak(d, f)) return BridgeType::None;
    if ((bond(c, e) && bond(e, a)) || (bond(f, b) && bond(b, d))) return BridgeType::Parallel;
    if ((bond(c, d) && bond(f, a)) || (bond(e, b) && bond(b, e))) return BridgeType::Antiparallel;
    return BridgeType::None;
  }

  void sheets() {
    std::vector<Ladder> ladders;
    for (const auto& [i, j] : near_) {
      const BridgeType t = bridge(i, j);
      if (t == BridgeType::None) continue;
      bool found = false;
      for (Ladder& l : ladders) {
        if (l.type != t || i != l.i.back() + 1) continue;
        if (t == BridgeType::Parallel && l.j.back() + 1 == j) {
          l.i.push_back(i);
          l.j.push_back(j);
          found = true;
          break;
        }
        if (t == BridgeType::Antiparallel && l.j.front() - 1 == j) {
          l.i.push_back(i);
          l.j.push_front(j);
          found = true;
          break;
        }
      }
      if (!found) ladders.push_back({t, r_[i].chain, {i}, {j}});
    }

    // Bulge-linked ladders: at most one extra residue on one strand and four on the other.
    // Gaps are compared as unsigned values, so a ladder that does not follow on
    // both strands is never linked.
    std::stable_sort(ladders.begin(), ladders.end(), [](const Ladder& x, const Ladder& y) {
      return x.chain != y.chain ? x.chain < y.chain : x.i.front() < y.i.front();
    });
    for (std::size_t a = 0; a < ladders.size(); ++a)
      for (std::size_t b = a + 1; b < ladders.size(); ++b) {
        Ladder& x = ladders[a];
        const Ladder& y = ladders[b];
        const unsigned ibi = x.i.front(), iei = x.i.back(), jbi = x.j.front(), jei = x.j.back();
        const unsigned ibj = y.i.front(), iej = y.i.back(), jbj = y.j.front(), jej = y.j.back();
        if (x.type != y.type || !noBreak(static_cast<int>(std::min(ibi, ibj)), static_cast<int>(std::max(iei, iej))) ||
            !noBreak(static_cast<int>(std::min(jbi, jbj)), static_cast<int>(std::max(jei, jej))) || ibj - iei >= 6 ||
            (iei >= ibj && ibi <= iej))
          continue;
        bool bulge;
        if (x.type == BridgeType::Parallel)
          bulge = (jbj - jei < 6 && ibj - iei < 3) || jbj - jei < 3;
        else
          bulge = (jbi - jej < 6 && ibj - iei < 3) || jbi - jej < 3;
        if (!bulge) continue;
        x.i.insert(x.i.end(), y.i.begin(), y.i.end());
        if (x.type == BridgeType::Parallel)
          x.j.insert(x.j.end(), y.j.begin(), y.j.end());
        else
          x.j.insert(x.j.begin(), y.j.begin(), y.j.end());
        ladders.erase(ladders.begin() + static_cast<std::ptrdiff_t>(b));
        b = a;
      }

    for (const Ladder& l : ladders) {
      const char ss = l.i.size() > 1 ? 'E' : 'B';
      for (const auto* strand : {&l.i, &l.j})
        for (int k = strand->front(); k <= strand->back(); ++k)
          if (r_[k].ss != 'E') r_[k].ss = ss;
    }
  }

  bool helixStart(int i, int s) const {
    const HelixFlag f = r_[i].flag[s];
    return f == HelixFlag::Start || f == HelixFlag::StartAndEnd;
  }

  void helices() {
    for (int s = 0; s < 3; ++s) {
      const int stride = s + 3;
      for (int i = 0; i + stride < n(); ++i) {
        if (!noBreak(i, i + stride) || !bond(i + stride, i)) continue;
        r_[i + stride].flag[s] = HelixFlag::End;
        for (int j = i + 1; j < i + stride; ++j)
          if (r_[j].flag[s] == HelixFlag::None) r_[j].flag[s] = HelixFlag::Middle;
        r_[i].flag[s] = r_[i].flag[s] == HelixFlag::End ? HelixFlag::StartAndEnd : HelixFlag::Start;
      }
    }
    for (int i = 1; i + 4 < n(); ++i)
      if (helixStart(i, 1) && helixStart(i - 1, 1))
        for (int j = i; j <= i + 3; ++j) r_[j].ss = 'H';
    for (int i = 1; i + 3 < n(); ++i)
      if (helixStart(i, 0) && helixStart(i - 1, 0)) {
        bool empty = true;
        for (int j = i; j <= i + 2; ++j) empty &= r_[j].ss == '.' || r_[j].ss == 'G';
        if (empty)
          for (int j = i; j <= i + 2; ++j) r_[j].ss = 'G';
      }
    for (int i = 1; i + 5 < n(); ++i)
      if (helixStart(i, 2) && helixStart(i - 1, 2)) {
        bool empty = true;
        for (int j = i; j <= i + 4; ++j) empty &= r_[j].ss == '.' || r_[j].ss == 'I' || r_[j].ss == 'H';
        if (empty)
          for (int j = i; j <= i + 4; ++j) r_[j].ss = 'I';
      }
    for (int i = 1; i + 1 < n(); ++i) {
      if (r_[i].ss != '.') continue;
      bool turn = false;
      for (int s = 0; s < 3 && !turn; ++s)
        for (int l = 1; l < s + 3 && !turn; ++l) turn = i - l >= 0 && helixStart(i - l, s);
      if (turn)
        r_[i].ss = 'T';
      else if (bend(i))
        r_[i].ss = 'S';
    }
  }

  bool bend(int i) const {
    if (i < 2 || i + 2 >= n() || !noBreak(i - 2, i + 2)) return false;
    const V u = r_[i].ca - r_[i - 2].ca, v = r_[i + 2].ca - r_[i].ca;
    const double nu = std::sqrt(dot(u, u)), nv = std::sqrt(dot(v, v));
    if (nu == 0 || nv == 0) return false;
    return std::acos(std::clamp(static_cast<double>(dot(u, v)) / (nu * nv), -1.0, 1.0)) * 180.0 / std::numbers::pi >
           kBendAngle;
  }

  bool ppAngles(int i) const {
    return r_[i].phi != 360 && r_[i].psi != 360 && std::abs(r_[i].phi - kPPPhi) <= kPPRange &&
           std::abs(r_[i].psi - kPPPsi) <= kPPRange;
  }

  void ppHelices() {
    for (int i = 1; i + 3 < n(); ++i) {
      bool ok = true;
      for (int k = 0; k < kPPMinLength && ok; ++k) ok = ppAngles(i + k);
      if (!ok) continue;
      for (int k = 0; k < kPPMinLength; ++k)
        if (r_[i + k].flag[3] == HelixFlag::None) r_[i + k].flag[3] = HelixFlag::Middle;
    }
    for (Res& x : r_)
      if (x.flag[3] != HelixFlag::None && x.ss == '.') x.ss = 'P';
  }
};

}  // namespace

void assignDssp(const std::vector<const Molecule*>& mols, const std::vector<Topology*>& topos) {
  std::vector<Res> res;
  int nextChain = 0;
  for (std::size_t mi = 0; mi < mols.size(); ++mi) {
    const Molecule& m = *mols[mi];
    Topology& t = *topos[mi];
    if (t.backboneType != BackboneType::Protein) continue;
    t.dsspSS.assign(m.residues.size(), ' ');
    std::vector<int> heavyDegree(m.atoms.size(), 0);
    for (const Bond& b : m.bonds) {
      heavyDegree[b.a] += m.atoms[b.b].atomicNum > 1;
      heavyDegree[b.b] += m.atoms[b.a].atomicNum > 1;
    }
    std::vector<bool> seen(t.units.size(), false);
    auto walk = [&](int start) {
      const int chainIndex = nextChain++;
      bool first = true;
      for (int u = start; u >= 0 && !seen[u]; u = t.units[u].next) {
        seen[u] = true;
        const auto& a = t.units[u].atoms;
        Res r;
        r.mol = static_cast<int>(mi);
        r.residue = t.units[u].residue;
        r.n = toV(m.atoms[a[0]].xyz);
        r.ca = toV(m.atoms[a[1]].xyz);
        r.c = toV(m.atoms[a[2]].xyz);
        r.o = toV(m.atoms[a[3]].xyz);
        r.chain = chainIndex;
        r.breakBefore = first;
        r.donor = heavyDegree[a[0]] - (first ? 0 : 1) == 1;
        // As in DSSP, H is placed from the C=O of the preceding residue in the list,
        // also across chain breaks; the first residue keeps H at N.
        r.h = r.n;
        if (r.donor && !res.empty()) r.h = r.n + (res.back().c - res.back().o) / dist(res.back().c, res.back().o);
        if (!first) {
          const Res& p = res.back();
          r.phi = dihedral(p.c, r.n, r.ca, r.c);
          res.back().psi = dihedral(p.n, p.ca, p.c, r.n);
        }
        res.push_back(r);
        first = false;
      }
    };
    for (int u = 0; u < static_cast<int>(t.units.size()); ++u)
      if (t.units[u].prev < 0) walk(u);
    for (int u = 0; u < static_cast<int>(t.units.size()); ++u)
      if (!seen[u]) walk(u);  // cyclic chains
  }
  if (res.empty()) return;
  Dssp d(std::move(res));
  d.run();
  for (const Res& r : d.residues())
    if (r.residue >= 0) topos[r.mol]->dsspSS[r.residue] = r.ss;
}

}  // namespace od
