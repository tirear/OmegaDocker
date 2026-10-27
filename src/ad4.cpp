// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Parameter values from AutoDock-GPU 1.6 (host/src/processligand.cpp,
// host/inc/getparameters.h), LGPL-2.1-or-later:
// Copyright (C) 2017 TU Darmstadt, Embedded Systems and Applications Group, Germany.
// Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

#include "ad4.hpp"

#include <array>
#include <cmath>

namespace od {
namespace {

constexpr std::array<Ad4TypeParams, static_cast<int>(Ad4Type::Count)> kParams{{
    {"H", 2.00, 0.020, 0.0, 0.0, 0.0000, 0.00051},
    {"HD", 2.00, 0.020, 0.0, 1.0, 0.0000, 0.00051},
    {"HS", 2.00, 0.020, 0.0, 1.0, 0.0000, 0.00051},
    {"C", 4.00, 0.150, 0.0, 0.0, 33.5103, -0.00143},
    {"A", 4.00, 0.150, 0.0, 0.0, 33.5103, -0.00052},
    {"N", 3.50, 0.160, 0.0, 0.0, 22.4493, -0.00162},
    {"NA", 3.50, 0.160, 1.9, 5.0, 22.4493, -0.00162},
    {"NS", 3.50, 0.160, 1.9, 5.0, 22.4493, -0.00162},
    {"OA", 3.20, 0.200, 1.9, 5.0, 17.1573, -0.00251},
    {"OS", 3.20, 0.200, 1.9, 5.0, 17.1573, -0.00251},
    {"F", 3.09, 0.080, 0.0, 0.0, 15.4480, -0.00110},
    {"Mg", 1.30, 0.875, 0.0, 0.0, 1.5600, -0.00110},
    {"P", 4.20, 0.200, 0.0, 0.0, 38.7924, -0.00110},
    {"SA", 4.00, 0.200, 2.5, 1.0, 33.5103, -0.00214},
    {"S", 4.00, 0.200, 0.0, 0.0, 33.5103, -0.00214},
    {"Cl", 4.09, 0.276, 0.0, 0.0, 35.8235, -0.00110},
    {"Ca", 1.98, 0.550, 0.0, 0.0, 2.7700, -0.00110},
    {"Mn", 1.30, 0.875, 0.0, 0.0, 2.1400, -0.00110},
    {"Fe", 1.30, 0.010, 0.0, 0.0, 1.8400, -0.00110},
    {"Zn", 1.48, 0.550, 0.0, 0.0, 1.7000, -0.00110},
    {"Br", 4.33, 0.389, 0.0, 0.0, 42.5661, -0.00110},
    {"I", 4.72, 0.550, 0.0, 0.0, 55.0585, -0.00110},
    {"Si", 4.10, 0.200, 0.0, 0.0, 35.8235, -0.00143},
    {"B", 3.84, 0.155, 0.0, 0.0, 29.6478, -0.00152},
}};

}  // namespace

const Ad4TypeParams& ad4Params(Ad4Type t) { return kParams[static_cast<int>(t)]; }

const char* ad4TypeName(Ad4Type t) { return t == Ad4Type::Merged ? "-" : kParams[static_cast<int>(t)].name; }

bool ad4IsHbDonor(Ad4Type t) { return t == Ad4Type::HD || t == Ad4Type::HS; }

bool ad4IsHbAcceptor(Ad4Type t) {
  return t == Ad4Type::NA || t == Ad4Type::NS || t == Ad4Type::OA || t == Ad4Type::OS || t == Ad4Type::SA;
}

bool ad4IsMetal(Ad4Type t) {
  return t == Ad4Type::Mg || t == Ad4Type::Ca || t == Ad4Type::Mn || t == Ad4Type::Fe || t == Ad4Type::Zn;
}

bool ad4Clashes(Ad4Type t) {
  return t != Ad4Type::Merged && t != Ad4Type::H && t != Ad4Type::HD && t != Ad4Type::HS && !ad4IsMetal(t);
}

Ad4Typing assignAd4Types(const Molecule& m) {
  const int n = static_cast<int>(m.atoms.size());
  std::vector<std::vector<int>> nb(n);
  std::vector<int> valence(n, 0);
  for (const Bond& b : m.bonds) {
    if (isMetalBond(m, b)) continue;
    nb[b.a].push_back(b.b);
    nb[b.b].push_back(b.a);
    valence[b.a] += static_cast<int>(b.order);
    valence[b.b] += static_cast<int>(b.order);
  }
  auto elem = [&](int a) { return m.atoms[a].atomicNum; };

  Ad4Typing out;
  out.type.assign(n, Ad4Type::C);
  out.charge.resize(n);
  for (int a = 0; a < n; ++a) out.charge[a] = m.atoms[a].partialCharge;
  for (int a = 0; a < n; ++a) {
    const Atom& at = m.atoms[a];
    const int x = static_cast<int>(nb[a].size());
    Ad4Type t;
    switch (at.atomicNum) {
      case 1:
        if (x == 0) {
          t = Ad4Type::H;
        } else {
          const int z = elem(nb[a].front());
          t = (z == 7 || z == 8 || z == 9 || z == 15 || z == 16) ? Ad4Type::HD : Ad4Type::Merged;
        }
        break;
      case 6: t = at.aromatic ? Ad4Type::A : Ad4Type::C; break;
      case 7: {
        t = Ad4Type::NA;
        if (at.formalCharge == 1) {
          t = Ad4Type::N;
        } else if (x == 3 && valence[a] == 3) {
          for (int o : nb[a]) {
            const bool sp2C = elem(o) == 6 && nb[o].size() == 3 && valence[o] == 4;
            if (m.atoms[o].aromatic || sp2C) t = Ad4Type::N;
          }
        }
        break;
      }
      case 8: t = Ad4Type::OA; break;
      case 16: t = !at.aromatic && x == 2 ? Ad4Type::SA : Ad4Type::S; break;
      case 15: t = Ad4Type::P; break;
      case 9: t = Ad4Type::F; break;
      case 17: t = Ad4Type::Cl; break;
      case 35: t = Ad4Type::Br; break;
      case 53: t = Ad4Type::I; break;
      case 12: t = Ad4Type::Mg; break;
      case 20: t = Ad4Type::Ca; break;
      case 25: t = Ad4Type::Mn; break;
      case 26: t = Ad4Type::Fe; break;
      case 30: t = Ad4Type::Zn; break;
      case 5: t = Ad4Type::B; break;
      case 14: t = Ad4Type::Si; break;
      default:
        t = Ad4Type::C;
        out.unsupported.push_back(a);
    }
    out.type[a] = t;
  }
  for (int a = 0; a < n; ++a)
    if (out.type[a] == Ad4Type::Merged) {
      out.charge[nb[a].front()] += out.charge[a];
      out.charge[a] = 0.0;
    }
  return out;
}

Ad4PairParams ad4PairParams(Ad4Type a, Ad4Type b, const Ad4Coeffs& c) {
  const Ad4TypeParams& pa = ad4Params(a);
  const Ad4TypeParams& pb = ad4Params(b);
  Ad4PairParams p;
  if ((ad4IsHbDonor(a) && ad4IsHbAcceptor(b)) || (ad4IsHbDonor(b) && ad4IsHbAcceptor(a))) {
    const double eps = c.hb * pa.epsHb * pb.epsHb;
    p.reqm = pa.reqmHb + pb.reqmHb;
    p.m = 12;
    p.n = 10;
    p.AC = 5 * eps * std::pow(p.reqm, 12);
    p.BD = 6 * eps * std::pow(p.reqm, 10);
  } else {
    const double eps = c.vdW * std::sqrt(pa.eps * pb.eps);
    p.reqm = 0.5 * (pa.reqm + pb.reqm);
    p.AC = eps * std::pow(p.reqm, 12);
    p.BD = 2 * eps * std::pow(p.reqm, 6);
  }
  return p;
}

}  // namespace od
