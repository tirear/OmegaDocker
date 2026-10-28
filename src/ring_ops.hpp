// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cmath>
#include <span>
#include <type_traits>
#include <vector>

#include "molecule.hpp"

namespace od {

// Ring conformers of one small ring system in a molecule with its torsion tree (spec 4.9).
// Conformer k > 0 replaces the coordinates of the moving ring atoms and moves each exit branch
// (the atoms beyond a bond from a ring atom to an atom outside the ring system) rigidly; all in
// the input frame, applied before the torsions. Conformer 0 is the input.
struct RingOp {
  int system = -1;               // index into RingConformers::systems
  int confCount = 1;
  std::vector<int> move;         // coordinate indices (atoms of the molecule, then virtual atoms)
  std::vector<int> exitAtom;     // coordinate indices moved with an exit
  std::vector<int> exitOf;       // per exitAtom entry: its exit
  int exitCount = 0;
  std::vector<double> pos;       // (confCount − 1) × move.size() × 3
  std::vector<double> xf;        // (confCount − 1) × exitCount × 12: R (row-major), t; x' = R x + t
  // Per coordinate index: 0 on the side that does not move, otherwise an id of the rigid piece
  // (a moving ring atom or an exit). Two coordinates with different ids may change their distance.
  std::vector<int> token;
};

// Ring operations of one molecule, in the order they are applied (a ring system inside an exit
// branch of another one comes first).
struct RingOps {
  std::vector<RingOp> ops;
};

// Applies the ring conformers k[i] (conformer index of ops.ops[i], rounded) to xyz.
template <class V, class K>
inline void applyRingOps(const RingOps& ops, std::span<const K> k, std::vector<V>& xyz) {
  using T = std::remove_cvref_t<decltype(xyz[0][0])>;
  for (std::size_t r = 0; r < ops.ops.size(); ++r) {
    const RingOp& op = ops.ops[r];
    const long c = std::lround(static_cast<double>(k[r]));
    if (c <= 0 || c >= op.confCount) continue;
    const double* P = op.pos.data() + static_cast<std::size_t>(c - 1) * op.move.size() * 3;
    for (std::size_t i = 0; i < op.move.size(); ++i)
      for (int a = 0; a < 3; ++a) xyz[op.move[i]][a] = static_cast<T>(P[3 * i + a]);
    const double* X = op.xf.data() + static_cast<std::size_t>(c - 1) * op.exitCount * 12;
    for (std::size_t i = 0; i < op.exitAtom.size(); ++i) {
      const double* M = X + 12 * op.exitOf[i];
      V& p = xyz[op.exitAtom[i]];
      const T p0 = p[0], p1 = p[1], p2 = p[2];
      for (int a = 0; a < 3; ++a)
        p[a] = static_cast<T>(M[3 * a]) * p0 + static_cast<T>(M[3 * a + 1]) * p1 + static_cast<T>(M[3 * a + 2]) * p2 +
               static_cast<T>(M[9 + a]);
    }
  }
}

}  // namespace od
