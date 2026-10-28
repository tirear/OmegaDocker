// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "ringconf.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <format>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>

#include <Eigen/Dense>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>

#include "pool.hpp"
#include "rng.hpp"

namespace od {
namespace {

std::vector<std::vector<std::pair<int, int>>> adjacency(const Molecule& m, const std::vector<bool>& skipBond) {
  std::vector<std::vector<std::pair<int, int>>> adj(m.atoms.size());
  for (std::size_t i = 0; i < m.bonds.size(); ++i) {
    if (skipBond[i]) continue;
    adj[m.bonds[i].a].push_back({m.bonds[i].b, static_cast<int>(i)});
    adj[m.bonds[i].b].push_back({m.bonds[i].a, static_cast<int>(i)});
  }
  return adj;
}

std::vector<bool> brokenOrMetal(const Molecule& m, const std::vector<int>& broken) {
  std::vector<bool> skip(m.bonds.size(), false);
  for (int b : broken) skip[b] = true;
  for (std::size_t i = 0; i < m.bonds.size(); ++i)
    if (isMetalBond(m, m.bonds[i])) skip[i] = true;
  return skip;
}

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 unit(const Vec3& a) {
  const double n = std::sqrt(dot(a, a));
  return n > 0 ? Vec3{a[0] / n, a[1] / n, a[2] / n} : Vec3{0, 0, 0};
}

// Rigid transform x' = R x + t.
struct Rigid {
  double R[3][3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  Vec3 t{0, 0, 0};
  Vec3 apply(const Vec3& x) const {
    Vec3 o;
    for (int i = 0; i < 3; ++i) o[i] = R[i][0] * x[0] + R[i][1] * x[1] + R[i][2] * x[2] + t[i];
    return o;
  }
};

// Orthonormal frame (columns x, y, z) and origin.
struct Frame {
  Vec3 o, e[3];
};
// Origin o, x along p − o, z normal to the plane of p and q.
Frame frameAlong(const Vec3& o, const Vec3& p, const Vec3& q) {
  Frame f;
  f.o = o;
  f.e[0] = unit(sub(p, o));
  f.e[2] = unit(cross(f.e[0], sub(q, o)));
  f.e[1] = cross(f.e[2], f.e[0]);
  return f;
}
// Origin a, x along the bisector of the bonds to n1 and n2, z normal to their plane.
Frame frameBisector(const Vec3& a, const Vec3& n1, const Vec3& n2) {
  const Vec3 u1 = unit(sub(n1, a)), u2 = unit(sub(n2, a));
  Frame f;
  f.o = a;
  f.e[0] = unit({u1[0] + u2[0], u1[1] + u2[1], u1[2] + u2[2]});
  f.e[2] = unit(cross(u1, u2));
  f.e[1] = cross(f.e[2], f.e[0]);
  return f;
}
// The transform taking frame a onto frame b.
Rigid frameMap(const Frame& a, const Frame& b) {
  Rigid r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0;
      for (int k = 0; k < 3; ++k) s += b.e[k][i] * a.e[k][j];
      r.R[i][j] = s;
    }
  const Vec3 ro = r.apply(a.o);
  for (int i = 0; i < 3; ++i) r.t[i] = b.o[i] - ro[i];
  return r;
}

// Least-squares superposition of p onto q (Kabsch).
Rigid kabsch(const std::vector<Vec3>& p, const std::vector<Vec3>& q) {
  const std::size_t n = p.size();
  Vec3 cp{0, 0, 0}, cq{0, 0, 0};
  for (std::size_t i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k) {
      cp[k] += p[i][k] / static_cast<double>(n);
      cq[k] += q[i][k] / static_cast<double>(n);
    }
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  for (std::size_t i = 0; i < n; ++i)
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) H(a, b) += (p[i][a] - cp[a]) * (q[i][b] - cq[b]);
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
  if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0) D(2, 2) = -1;
  const Eigen::Matrix3d R = svd.matrixV() * D * svd.matrixU().transpose();
  Rigid r;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) r.R[i][j] = R(i, j);
  for (int i = 0; i < 3; ++i) r.t[i] = cq[i] - (R(i, 0) * cp[0] + R(i, 1) * cp[1] + R(i, 2) * cp[2]);
  return r;
}

double signedVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  return dot(sub(b, a), cross(sub(c, a), sub(d, a)));
}

double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) { return dihedralDeg(a, b, c, d); }

}  // namespace

int RingConformers::conformerCount() const {
  int n = 0;
  for (const auto& s : systems) n += static_cast<int>(s.conf.size());
  return n;
}

std::vector<int> smallestCycles(const Molecule& m, const std::vector<int>& broken, int limit) {
  const std::vector<bool> skip = brokenOrMetal(m, broken);
  const std::vector<bool> ring = ringBondMask(m, broken);
  const auto adj = adjacency(m, skip);
  std::vector<int> out(m.bonds.size(), 0);
  std::vector<int> depth(m.atoms.size(), -1);
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) {
    if (skip[bi] || !ring[bi]) continue;
    const int a = m.bonds[bi].a, b = m.bonds[bi].b;
    // Shortest path a -> b without this bond, up to limit - 1 bonds.
    std::vector<int> q{a};
    depth[a] = 0;
    int found = 0;
    for (std::size_t i = 0; i < q.size() && !found; ++i) {
      const int u = q[i];
      if (depth[u] >= limit - 2) continue;
      for (const auto& [v, e] : adj[u]) {
        if (e == static_cast<int>(bi) || depth[v] >= 0) continue;
        depth[v] = depth[u] + 1;
        if (v == b) {
          found = depth[v] + 1;
          break;
        }
        q.push_back(v);
      }
    }
    for (int u : q) depth[u] = -1;
    depth[b] = -1;
    out[bi] = found ? found : limit;
  }
  return out;
}

std::vector<RingSystem> findRingSystems(const Molecule& m, const std::vector<int>& broken, int macroMinSize) {
  const std::vector<int> cyc = smallestCycles(m, broken, macroMinSize);
  std::vector<int> parent(m.atoms.size());
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](int x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  std::vector<bool> inRing(m.atoms.size(), false);
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) {
    if (cyc[bi] == 0 || cyc[bi] >= macroMinSize) continue;
    const int a = m.bonds[bi].a, b = m.bonds[bi].b;
    inRing[a] = inRing[b] = true;
    parent[find(a)] = find(b);
  }
  std::map<int, int> sysOf;
  std::vector<RingSystem> out;
  for (int a = 0; a < static_cast<int>(m.atoms.size()); ++a) {
    if (!inRing[a]) continue;
    auto [it, fresh] = sysOf.try_emplace(find(a), static_cast<int>(out.size()));
    if (fresh) out.emplace_back();
    out[it->second].ring.push_back(a);
  }
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi)
    if (cyc[bi] > 0 && cyc[bi] < macroMinSize) {
      RingSystem& s = out[sysOf[find(m.bonds[bi].a)]];
      s.largestRing = std::max(s.largestRing, cyc[bi]);
    }
  for (RingSystem& s : out) {
    std::set<int> sh;
    const std::set<int> rs(s.ring.begin(), s.ring.end());
    for (const Bond& b : m.bonds) {
      if (isMetalBond(m, b)) continue;
      for (auto [x, y] : {std::pair{b.a, b.b}, std::pair{b.b, b.a}})
        if (rs.count(x) && !rs.count(y) && m.atoms[y].atomicNum > 1 && !isMetal(m.atoms[y].atomicNum)) sh.insert(y);
    }
    s.shell.assign(sh.begin(), sh.end());
    s.conf.push_back({});
    for (int a : s.ring) s.conf[0].push_back(m.atoms[a].xyz);
    for (int a : s.shell) s.conf[0].push_back(m.atoms[a].xyz);
  }
  return out;
}

std::vector<int> macroRingBreaks(const Molecule& m, const Topology& t, int macroMinSize) {
  std::map<int, const Torsion*> torOf;
  for (const Torsion& tor : t.torsions) torOf[tor.bond] = &tor;
  std::vector<int> broken = t.brokenBonds, out;
  std::vector<int> doubles(m.atoms.size(), 0);
  for (const Bond& b : m.bonds)
    if (b.order == BondOrder::Double || b.aromatic) ++doubles[b.a], ++doubles[b.b];
  auto carbonylToN = [&](int c, int n) {
    if (m.atoms[c].atomicNum != 6 || m.atoms[n].atomicNum != 7) return false;
    for (const Bond& b : m.bonds)
      if (b.order == BondOrder::Double && (b.a == c || b.b == c) && m.atoms[b.a == c ? b.b : b.a].atomicNum == 8) return true;
    return false;
  };
  for (;;) {
    const std::vector<int> cyc = smallestCycles(m, broken, macroMinSize);
    int best = -1, bestClass = INT_MAX;
    for (std::size_t bi = 0; bi < m.bonds.size(); ++bi) {
      const Bond& b = m.bonds[bi];
      if (cyc[bi] < macroMinSize || b.aromatic || b.order != BondOrder::Single || isMetalBond(m, b)) continue;
      if (m.role == Role::StaM) {
        auto it = torOf.find(static_cast<int>(bi));
        if (it == torOf.end() || !it->second->inFlexRange) continue;
      }
      const int cls = doubles[b.a] == 0 && doubles[b.b] == 0 ? 0 : carbonylToN(b.a, b.b) || carbonylToN(b.b, b.a) ? 1 : 2;
      if (cls < bestClass) {
        bestClass = cls;
        best = static_cast<int>(bi);
      }
    }
    if (best < 0) break;
    broken.push_back(best);
    out.push_back(best);
  }
  return out;
}

namespace {

// ETKDG candidates of one system, MMFF94 (UFF) optimized; fills s.conf (beyond the input),
// s.energy, s.generated, s.note.
void generateSystem(const Molecule& m, RingSystem& s, const RingParams& p,
                    std::uint64_t seed) {
  bool molHasH = false;
  for (const Atom& a : m.atoms) molHasH |= a.atomicNum == 1;
  const std::size_t nr = s.ring.size();
  // Fragment: ring atoms, shell atoms, and the heavy neighbors of the shell atoms (so that their
  // hybridization, e.g. of an amide, stays); other bonds are replaced by hydrogens.
  std::vector<int> frag(s.ring);
  frag.insert(frag.end(), s.shell.begin(), s.shell.end());
  std::map<int, int> idx;
  for (std::size_t i = 0; i < frag.size(); ++i) idx[frag[i]] = static_cast<int>(i);
  {
    std::set<int> second;
    for (const Bond& b : m.bonds) {
      if (isMetalBond(m, b)) continue;
      for (auto [x, y] : {std::pair{b.a, b.b}, std::pair{b.b, b.a}}) {
        auto ix = idx.find(x);
        if (ix != idx.end() && ix->second >= static_cast<int>(nr) && !idx.count(y) && m.atoms[y].atomicNum > 1 &&
            !isMetal(m.atoms[y].atomicNum))
          second.insert(y);
      }
    }
    for (int a : second) {
      idx[a] = static_cast<int>(frag.size());
      frag.push_back(a);
    }
  }
  const std::size_t nStore = nr + s.shell.size();  // ring and shell atoms are kept per conformer
  std::vector<int> hCount(m.atoms.size(), 0);
  for (const Bond& b : m.bonds) {
    if (m.atoms[b.a].atomicNum == 1) ++hCount[b.b];
    if (m.atoms[b.b].atomicNum == 1) ++hCount[b.a];
  }
  RDKit::RWMol mol;
  auto* conf = new RDKit::Conformer(static_cast<unsigned>(frag.size()));
  conf->set3D(true);
  for (std::size_t i = 0; i < frag.size(); ++i) {
    const Atom& a = m.atoms[frag[i]];
    RDKit::Atom at(a.atomicNum);
    at.setFormalCharge(a.formalCharge);
    if (i < nr && molHasH) {
      at.setNoImplicit(true);
      at.setNumExplicitHs(static_cast<unsigned>(hCount[frag[i]]));
    } else {
      at.setNoImplicit(false);
    }
    mol.addAtom(&at, true, false);
    conf->setAtomPos(static_cast<unsigned>(i), RDGeom::Point3D(a.xyz[0], a.xyz[1], a.xyz[2]));
  }
  mol.addConformer(conf, true);
  std::vector<std::vector<int>> nb(frag.size());
  for (const Bond& b : m.bonds) {
    auto ia = idx.find(b.a), ib = idx.find(b.b);
    if (ia == idx.end() || ib == idx.end() || isMetalBond(m, b)) continue;
    const RDKit::Bond::BondType bt = b.order == BondOrder::Double   ? RDKit::Bond::DOUBLE
                                     : b.order == BondOrder::Triple ? RDKit::Bond::TRIPLE
                                                                    : RDKit::Bond::SINGLE;
    mol.addBond(static_cast<unsigned>(ia->second), static_cast<unsigned>(ib->second), bt);
    nb[ia->second].push_back(ib->second);
    nb[ib->second].push_back(ia->second);
  }
  try {
    RDKit::MolOps::sanitizeMol(mol);
    RDKit::MolOps::assignStereochemistryFrom3D(mol);
    RDKit::MolOps::addHs(mol);
  } catch (const std::exception& e) {
    s.note = std::format("RDKit could not prepare the fragment ({})", e.what());
    return;
  }
  RDKit::DGeomHelpers::EmbedParameters ep = RDKit::DGeomHelpers::ETKDGv3;
  ep.randomSeed = static_cast<int>(seed & 0x7fffffff);
  ep.numThreads = 1;
  const std::vector<int> cids = RDKit::DGeomHelpers::EmbedMultipleConfs(mol, static_cast<unsigned>(p.confGenNum), ep);
  s.generated = static_cast<int>(cids.size());
  if (cids.empty()) {
    s.note = "ETKDG embedded no conformer";
    return;
  }
  std::vector<std::pair<int, double>> res;
  if (RDKit::MMFF::MMFFMolProperties(mol).isValid())
    RDKit::MMFF::MMFFOptimizeMoleculeConfs(mol, res, 1, 2000);
  else
    RDKit::UFF::UFFOptimizeMoleculeConfs(mol, res, 1, 2000);

  // Ring dihedrals (both bond atoms have another ring neighbor) and the centres checked for inversion.
  struct Dih {
    int x, a, b, y;
  };
  std::vector<Dih> dih;
  for (std::size_t i = 0; i < nr; ++i)
    for (int j : nb[i]) {
      if (j <= static_cast<int>(i) || j >= static_cast<int>(nr)) continue;
      int x = -1, y = -1;
      for (int k : nb[i])
        if (k != j && k < static_cast<int>(nr) && (x < 0 || k < x)) x = k;
      for (int k : nb[j])
        if (k != static_cast<int>(i) && k < static_cast<int>(nr) && (y < 0 || k < y)) y = k;
      if (x >= 0 && y >= 0) dih.push_back({x, static_cast<int>(i), j, y});
    }
  for (const Dih& d : dih) s.dih.push_back({d.x, d.a, d.b, d.y});
  auto dihedrals = [&](const std::vector<Vec3>& c) {
    std::vector<double> v;
    for (const Dih& d : dih) v.push_back(dihedral(c[d.x], c[d.a], c[d.b], c[d.y]));
    return v;
  };
  auto rmsd = [&](const std::vector<double>& u, const std::vector<double>& v) {
    double s = 0;
    for (std::size_t i = 0; i < u.size(); ++i) s += std::pow(wrapDeg(u[i] - v[i]), 2);
    return u.empty() ? 0.0 : std::sqrt(s / static_cast<double>(u.size()));
  };
  // Configuration of the tetrahedral ring atoms (|normalized signed volume| of the input > 0.3;
  // planar centres such as an amide N are not checked).
  std::vector<std::array<int, 4>> centres;
  for (std::size_t i = 0; i < nr; ++i)
    if (nb[i].size() >= 3) {
      std::vector<int> n3(nb[i]);
      std::sort(n3.begin(), n3.end());
      if (n3[2] >= static_cast<int>(nStore)) continue;
      const auto& c = s.conf[0];
      const Vec3 u = sub(c[n3[0]], c[i]), v = sub(c[n3[1]], c[i]), w = sub(c[n3[2]], c[i]);
      const double nv = std::sqrt(dot(u, u) * dot(v, v) * dot(w, w));
      if (nv > 0 && std::abs(dot(u, cross(v, w))) / nv > 0.3) centres.push_back({static_cast<int>(i), n3[0], n3[1], n3[2]});
    }
  auto chirality = [&](const std::vector<Vec3>& c) {
    std::vector<int> sg;
    for (const auto& q : centres) sg.push_back(signedVolume(c[q[0]], c[q[1]], c[q[2]], c[q[3]]) > 0 ? 1 : -1);
    return sg;
  };
  const std::vector<int> chi0 = chirality(s.conf[0]);
  struct Cand {
    double e;
    std::vector<Vec3> x;
    std::vector<double> d;
  };
  std::vector<Cand> cand;
  for (std::size_t c = 0; c < cids.size(); ++c) {
    const RDKit::Conformer& cf = mol.getConformer(cids[c]);
    std::vector<Vec3> x(nStore);
    for (std::size_t i = 0; i < nStore; ++i) {
      const auto& pt = cf.getAtomPos(static_cast<unsigned>(i));
      x[i] = {pt.x, pt.y, pt.z};
    }
    if (chirality(x) != chi0) continue;
    const double e = c < res.size() ? res[c].second : 0.0;
    cand.push_back({e, x, dihedrals(x)});
  }
  if (cand.empty()) {
    s.note = "no ETKDG conformer kept the configuration of the ring atoms";
    return;
  }
  std::stable_sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b) { return a.e < b.e; });
  const double eMin = cand[0].e;
  std::vector<std::vector<double>> kept{dihedrals(s.conf[0])};
  for (const Cand& c : cand) {
    if (static_cast<int>(s.conf.size()) >= p.maxConfNum) break;
    if (c.e - eMin > p.confEnergyWindow) break;
    bool fresh = true;
    for (const auto& k : kept) fresh &= rmsd(c.d, k) > p.confRmsd;
    if (!fresh) continue;
    kept.push_back(c.d);
    s.conf.push_back(c.x);
    s.energy.push_back(c.e - eMin);
  }
  if (s.conf.size() == 1) s.note = "every candidate is the input conformer";
}

}  // namespace

int nearestRingConformer(const RingSystem& s, std::span<const Vec3> xyz, double* rmsdDeg) {
  if (rmsdDeg) *rmsdDeg = 0;
  if (s.conf.size() < 2 || s.dih.empty()) return 0;
  std::vector<double> want;
  for (const auto& d : s.dih)
    want.push_back(dihedral(xyz[s.ring[d[0]]], xyz[s.ring[d[1]]], xyz[s.ring[d[2]]], xyz[s.ring[d[3]]]));
  int best = 0;
  double bestR = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < s.conf.size(); ++k) {
    const auto& c = s.conf[k];
    double sum = 0;
    for (std::size_t i = 0; i < s.dih.size(); ++i) {
      const auto& d = s.dih[i];
      sum += std::pow(wrapDeg(dihedral(c[d[0]], c[d[1]], c[d[2]], c[d[3]]) - want[i]), 2);
    }
    const double r = std::sqrt(sum / static_cast<double>(s.dih.size()));
    if (r < bestR) {
      bestR = r;
      best = static_cast<int>(k);
    }
  }
  if (rmsdDeg) *rmsdDeg = bestR;
  return best;
}

RingConformers generateRingConformers(const Molecule& m, const Topology& t, const RingParams& p, std::uint64_t stream) {
  RingConformers rc;
  rc.systems = findRingSystems(m, t.brokenBonds, p.macroMinSize);
  if (rc.systems.empty()) return rc;
  const std::vector<int> cyc = smallestCycles(m, t.brokenBonds, p.macroMinSize);
  std::set<int> small;
  for (std::size_t bi = 0; bi < m.bonds.size(); ++bi)
    if (cyc[bi] > 0 && cyc[bi] < p.macroMinSize) small.insert(static_cast<int>(bi));
  std::vector<bool> skip(m.bonds.size(), false);
  for (int b : t.brokenBonds) skip[b] = true;
  const auto adj = adjacency(m, skip);
  std::vector<int> todo;
  for (std::size_t si = 0; si < rc.systems.size(); ++si) {
    RingSystem& s = rc.systems[si];
    const std::set<int> rs(s.ring.begin(), s.ring.end());
    bool flex = false;
    for (const Torsion& tor : t.torsions)
      if (tor.inFlexRange && (rs.count(tor.atoms[1]) || rs.count(tor.atoms[2]))) flex = true;
    if (!flex) {
      s.note = "outside the flexible range";
      continue;
    }
    bool allAromatic = true;
    for (int bi : small)
      if (rs.count(m.bonds[bi].a) && rs.count(m.bonds[bi].b)) allAromatic &= m.bonds[bi].aromatic;
    if (allAromatic) {
      s.note = "aromatic";
      continue;
    }
    // Every exit must leave the system.
    bool bridges = true;
    std::vector<int> seen(m.atoms.size(), 0);
    for (int a : s.ring)
      for (const auto& [e, bi] : adj[a]) {
        if (rs.count(e) || !bridges) continue;
        std::vector<int> q{e};
        std::fill(seen.begin(), seen.end(), 0);
        seen[e] = seen[a] = 1;
        for (std::size_t i = 0; i < q.size() && bridges; ++i)
          for (const auto& [v, b2] : adj[q[i]]) {
            if (seen[v]) continue;
            if (rs.count(v)) {
              bridges = false;
              break;
            }
            seen[v] = 1;
            q.push_back(v);
          }
      }
    if (!bridges) {
      s.note = "part of a larger unbroken ring";
      continue;
    }
    todo.push_back(static_cast<int>(si));
  }
  parallelEach(static_cast<long long>(todo.size()), [&](long long k) {
    const int si = todo[k];
    Rng r(p.seed, stream, static_cast<std::uint64_t>(si));
    generateSystem(m, rc.systems[si], p, r.next());
  });
  return rc;
}

RingOps buildRingOps(const Molecule& m, const TorsionTree& tree, const RingConformers& rc, std::span<const int> toBase,
                     std::vector<std::string>* notes) {
  RingOps out;
  const int na = static_cast<int>(m.atoms.size());
  const int nc = na + static_cast<int>(tree.virtuals.size());
  std::vector<int> localOf;
  if (!toBase.empty()) {
    int maxBase = 0;
    for (int b : toBase) maxBase = std::max(maxBase, b);
    localOf.assign(static_cast<std::size_t>(maxBase) + 1, -1);
    for (int i = 0; i < na; ++i)
      if (toBase[i] >= 0) localOf[toBase[i]] = i;
  }
  auto local = [&](int base) { return toBase.empty() ? base : base < static_cast<int>(localOf.size()) ? localOf[base] : -1; };
  std::set<std::pair<int, int>> brokenPair;
  for (const ClosureLink& l : tree.links) brokenPair.insert({std::min(l.a, l.b), std::max(l.a, l.b)});
  std::vector<std::vector<int>> adj(na);
  for (const Bond& b : m.bonds) {
    if (brokenPair.count({std::min(b.a, b.b), std::max(b.a, b.b)})) continue;
    adj[b.a].push_back(b.b);
    adj[b.b].push_back(b.a);
  }
  std::vector<char> anchor(na, 0);
  for (int a : tree.anchors) anchor[a] = 1;
  auto note = [&](int si, const std::string& why) {
    if (notes) notes->push_back(std::format("ring system {}: {}", si + 1, why));
  };

  struct Built {
    RingOp op;
    std::vector<int> ringAtoms;  // local
  };
  std::vector<Built> built;
  for (std::size_t si = 0; si < rc.systems.size(); ++si) {
    const RingSystem& s = rc.systems[si];
    if (s.conf.size() < 2) continue;
    const int nr = static_cast<int>(s.ring.size());
    std::vector<int> S, shell;
    bool mapped = true;
    for (int b : s.ring) S.push_back(local(b)), mapped &= S.back() >= 0;
    for (int b : s.shell) shell.push_back(local(b)), mapped &= shell.back() >= 0;
    if (!mapped) {
      note(static_cast<int>(si), "atoms missing in this molecule");
      continue;
    }
    std::vector<int> fragIdx(na, -1);  // local atom -> index in the conformer arrays
    for (int i = 0; i < nr; ++i) fragIdx[S[i]] = i;
    for (std::size_t i = 0; i < shell.size(); ++i) fragIdx[shell[i]] = nr + static_cast<int>(i);
    std::vector<char> inS(na, 0);
    for (int a : S) inS[a] = 1;
    auto ringNeighbors = [&](int a) {
      std::vector<int> v;
      for (int o : adj[a])
        if (inS[o]) v.push_back(o);
      std::sort(v.begin(), v.end());
      return v;
    };
    // Exits and their branches.
    struct Exit {
      int a = -1, e = -1;      // e = -1: a virtual atom on a ring atom (local frame at a)
      std::vector<int> atoms;  // coordinate indices
      bool anchored = false;
    };
    std::vector<Exit> exits;
    std::vector<int> branchOf(nc, -1);
    bool ok = true;
    for (int a : S) {
      for (int e : adj[a]) {
        if (inS[e]) continue;
        Exit x{a, e, {}, false};
        std::vector<int> q{e};
        std::vector<char> seen(na, 0);
        seen[e] = seen[a] = 1;
        for (std::size_t i = 0; i < q.size(); ++i)
          for (int v : adj[q[i]]) {
            if (seen[v]) continue;
            if (inS[v]) ok = false;
            seen[v] = 1;
            q.push_back(v);
          }
        for (int v : q) {
          x.atoms.push_back(v);
          branchOf[v] = static_cast<int>(exits.size());
          x.anchored |= anchor[v] != 0;
        }
        exits.push_back(std::move(x));
      }
    }
    if (!ok) {
      note(static_cast<int>(si), "part of a larger unbroken ring in this molecule");
      continue;
    }
    for (std::size_t v = 0; v < tree.virtuals.size(); ++v) {
      const int p = tree.virtuals[v].parent, c = na + static_cast<int>(v);
      if (inS[p]) {
        exits.push_back({p, -1, {c}, false});
      } else if (branchOf[p] >= 0) {
        exits[branchOf[p]].atoms.push_back(c);
      }
    }
    std::vector<int> sAnchors;
    for (int a : S)
      if (anchor[a]) sAnchors.push_back(a);
    int anchoredExits = 0, parentExit = -1;
    for (std::size_t x = 0; x < exits.size(); ++x)
      if (exits[x].anchored) ++anchoredExits, parentExit = static_cast<int>(x);
    enum class Mode { Parent, All, Anchored } mode;
    if (sAnchors.empty() && anchoredExits == 1) mode = Mode::Parent;
    else if (m.role == Role::MobM || (sAnchors.empty() && anchoredExits == 0)) mode = Mode::All;
    else mode = Mode::Anchored;
    if (mode == Mode::Parent && fragIdx[exits[parentExit].e] < 0) mode = m.role == Role::MobM ? Mode::All : Mode::Anchored;
    if (mode == Mode::All)
      for (Exit& x : exits) x.anchored = false;
    std::vector<int> fitAtoms;  // local atoms with conformer positions, kept in place (Anchored)
    if (mode == Mode::Anchored) {
      fitAtoms = sAnchors;
      for (const Exit& x : exits)
        if (x.anchored && x.e >= 0 && fragIdx[x.e] >= 0) {
          fitAtoms.push_back(x.e);
          fitAtoms.push_back(x.a);
        }
      std::sort(fitAtoms.begin(), fitAtoms.end());
      fitAtoms.erase(std::unique(fitAtoms.begin(), fitAtoms.end()), fitAtoms.end());
      if (fitAtoms.size() < 3) {
        note(static_cast<int>(si), "too few fixed atoms to place its conformers");
        continue;
      }
    }
    std::vector<char> fixedS(na, 0);
    if (mode == Mode::Anchored)
      for (int a : sAnchors) fixedS[a] = 1;
    // Anchored: the frame of a fixed ring atom, a fixed ring neighbor and another fixed neighbor
    // (ring atom or atom of a fixed exit) when there is one; otherwise a fit of fitAtoms.
    int frameA = -1, frameB = -1, frameC = -1;
    if (mode == Mode::Anchored) {
      std::vector<char> fit(na, 0);
      for (int a : fitAtoms) fit[a] = 1;
      for (int a : sAnchors) {
        int b = -1, c = -1;
        for (int o : ringNeighbors(a))
          if (fixedS[o]) {
            b = o;
            break;
          }
        if (b < 0) continue;
        std::vector<int> nbs(adj[a]);
        std::sort(nbs.begin(), nbs.end());
        for (int o : nbs)
          if (o != b && fit[o] && fragIdx[o] >= 0) {
            c = o;
            break;
          }
        if (c < 0) continue;
        frameA = a, frameB = b, frameC = c;
        break;
      }
    }
    RingOp op;
    op.system = static_cast<int>(si);
    op.confCount = static_cast<int>(s.conf.size());
    for (int a : S)
      if (!fixedS[a]) op.move.push_back(a);
    std::vector<int> movingExit;
    for (std::size_t x = 0; x < exits.size(); ++x)
      if (!exits[x].anchored && !(mode == Mode::Parent && static_cast<int>(x) == parentExit)) {
        const int id = static_cast<int>(movingExit.size());
        movingExit.push_back(static_cast<int>(x));
        for (int c : exits[x].atoms) {
          op.exitAtom.push_back(c);
          op.exitOf.push_back(id);
        }
      }
    op.exitCount = static_cast<int>(movingExit.size());
    auto inputPos = [&](int c) { return c < na ? m.atoms[c].xyz : tree.virtuals[c - na].xyz; };
    auto ringX = [&](int a, int not1) {
      for (int o : ringNeighbors(a))
        if (o != not1) return o;
      return -1;
    };
    for (int k = 1; k < op.confCount; ++k) {
      const std::vector<Vec3>& C = s.conf[k];
      Rigid A;
      if (mode == Mode::Parent) {
        const Exit& x = exits[parentExit];
        const int X = ringX(x.a, -1);
        A = frameMap(frameAlong(C[fragIdx[x.a]], C[fragIdx[x.e]], C[fragIdx[X]]),
                     frameAlong(inputPos(x.a), inputPos(x.e), inputPos(X)));
      } else if (mode == Mode::Anchored && frameA >= 0) {
        A = frameMap(frameAlong(C[fragIdx[frameA]], C[fragIdx[frameB]], C[fragIdx[frameC]]),
                     frameAlong(inputPos(frameA), inputPos(frameB), inputPos(frameC)));
      } else {
        std::vector<Vec3> P, Q;
        for (int a : mode == Mode::All ? S : fitAtoms) {
          P.push_back(C[fragIdx[a]]);
          Q.push_back(inputPos(a));
        }
        A = kabsch(P, Q);
      }
      std::vector<Vec3> fin(na);  // final positions of the ring atoms (and aligned shell atoms)
      for (int a : S) fin[a] = fixedS[a] ? inputPos(a) : A.apply(C[fragIdx[a]]);
      for (int a : op.move)
        for (int c = 0; c < 3; ++c) op.pos.push_back(fin[a][c]);
      for (int id = 0; id < op.exitCount; ++id) {
        const Exit& x = exits[movingExit[id]];
        Rigid T;
        if (x.e >= 0 && fragIdx[x.e] >= 0) {
          const int X = ringX(x.a, -1);
          T = frameMap(frameAlong(inputPos(x.a), inputPos(x.e), inputPos(X)),
                       frameAlong(fin[x.a], A.apply(C[fragIdx[x.e]]), fin[X]));
        } else {
          const auto rn = ringNeighbors(x.a);
          T = frameMap(frameBisector(inputPos(x.a), inputPos(rn[0]), inputPos(rn[1])),
                       frameBisector(fin[x.a], fin[rn[0]], fin[rn[1]]));
        }
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) op.xf.push_back(T.R[i][j]);
        for (int i = 0; i < 3; ++i) op.xf.push_back(T.t[i]);
      }
    }
    op.token.assign(nc, 0);
    for (std::size_t i = 0; i < op.exitAtom.size(); ++i) op.token[op.exitAtom[i]] = 1 + op.exitOf[i];
    for (std::size_t i = 0; i < op.move.size(); ++i) op.token[op.move[i]] = 1 + op.exitCount + static_cast<int>(i);
    built.push_back({std::move(op), S});
  }
  // Order: a system moved by another one (its ring atoms carry a token of the other) comes first.
  const std::size_t nb = built.size();
  std::vector<std::vector<int>> before(nb);  // before[j]: systems to apply before j
  for (std::size_t i = 0; i < nb; ++i)
    for (std::size_t j = 0; j < nb; ++j)
      if (i != j && built[j].op.token[built[i].ringAtoms[0]] != 0) before[j].push_back(static_cast<int>(i));
  std::vector<int> state(nb, 0);
  std::function<void(int)> visit = [&](int j) {
    if (state[j]) return;
    state[j] = 1;
    for (int i : before[j]) visit(i);
    out.ops.push_back(built[j].op);
  };
  for (std::size_t j = 0; j < nb; ++j) visit(static_cast<int>(j));
  return out;
}

}  // namespace od
