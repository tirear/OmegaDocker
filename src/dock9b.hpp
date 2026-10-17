// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <vector>

#include "energy.hpp"
#include "grid.hpp"
#include "minimize.hpp"
#include "orient.hpp"

namespace od {

// 9b docking of one MobM against one StaM conformation (spec 8.4–8.6).
// CPU reference implementation. Coordinates, grid interpolation and per-atom /
// per-pair terms are single precision; sums, L-BFGS and pose updates are double
// (spec 12.1).

enum class Posm { Posm1, Posm2 };

struct Dock9bParams {
  Posm posm = Posm::Posm1;
  double posStep = 1.0;          // Å
  int hopfLevel = 3;             // posm2 orientations (rotStep)
  long long posm1SampleNum = 1000;
  Vec3 boxCenter{}, boxSize{};
  int maxOverlapAtomNum = 1;
  double maxNearDistance = 6.0;
  double clashDistance = 1.8;
  double clashBase = 10000.0;
  double clashK = 1000.0;
  // L-BFGS (spec 12)
  LbfgsParams lbfgs;
  double scaleTrans = 1, scaleRot = 1, scaleTors = 1;
  double maxTransStep = 0.5;     // Å
  double maxAngleStep = 10.0;    // degrees
  std::uint64_t seed = 0;
};

// MobM placement: world x = R(q) (y - center) + t, where y are the coordinates
// of the conformation (torsion-tree frame) and center its heavy-atom centre.
// The energy is evaluated at float(t), float(q).
struct Pose {
  Vec3 t{};
  Quat q{};
};

struct DockSample {
  Pose pose;                  // float values
  Vec3 center{};              // heavy-atom centre of the starting conformation (the Pose frame)
  std::vector<float> angles;  // degrees
  float energy = 0;           // inter + interClash + MobM intra + MobM clash; closure excluded (spec 4.5)
  float inter = 0;            // grid energy
  float interClash = 0;       // StaM–MobM clash penalty
  float intra = 0;            // MobM EnergyTerms::filterEnergy
  float closureMaxDeviation = 0;
  int steps = 0;
};

struct Dock9bStats {
  long long l1 = 0;          // L1poss generated
  long long overlap = 0;     // removed by 8.5
  long long far = 0;
  long long outsideBox = 0;  // posm1: sampled orientation does not fit in the box
  long long l2 = 0;          // minimized
  long long leftGrid = 0;    // removed during minimization (spec 8.6)
  long long steps = 0;
  long long evaluations = 0;
  Dock9bStats& operator+=(const Dock9bStats& o);
};

class MobmDocker {
 public:
  // sys: energy system of the single MobM (parts = {mob}); field: the StaM conformation (spec 8.1, 8.2).
  MobmDocker(const EnergySystem& sys, const StamField& field, const Dock9bParams& p);

  // Lattice indices of the positions keeping every probe atom kPosMargin inside
  // the box for one orientation (spec 8.4); empty when lo > hi on an axis.
  static constexpr double kPosMargin = 0.01;  // Å, absorbs float rounding
  struct PosRange {
    std::array<int, 3> lo{}, hi{};
    long long count() const;
  };

  // Conformation of a MobM sample in the torsion-tree frame.
  struct Conf {
    std::vector<double> anglesDeg;
    std::vector<Vec3f> coords;  // all coordinates of sys.evaluate, at float(anglesDeg)
    Vec3 center{};              // heavy-atom centre (unweighted)
    // posm2: position range per Hopf orientation and the prefix sums of their counts.
    std::vector<PosRange> hopfRange;
    std::vector<long long> hopfOffset;
  };
  Conf conformation(std::span<const double> anglesDeg) const;

  PosRange posRange(const Conf& c, const Quat& q) const;

  // L1pos k (0 <= k < l1Count(c)); stream: PRNG stream of the pair (posm1).
  // nullopt (posm1) when the sampled orientation does not fit in the box.
  std::optional<Pose> l1Pose(const Conf& c, std::uint64_t stream, long long k) const;
  // 8.5 quick filter; updates the counters of stats.
  bool quickFilter(const Conf& c, const Pose& pose, Dock9bStats& stats) const;

  // Energy of a pose; grad (optional): d/d(t, rotation vector, torsions in rad).
  struct PoseEnergy {
    EnergyTerms intra;
    double inter = 0, interClash = 0;
    int outside = 0;  // probe atoms outside the grid
    double objective() const { return inter + interClash + intra.total(); }
    double ranking() const { return inter + interClash + intra.filterEnergy(); }
  };
  // Single precision at float(pose), float(anglesDeg) (spec 8.6).
  PoseEnergy evaluate(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                      std::vector<double>* grad = nullptr) const;
  // Double-precision reference of evaluate (tests).
  PoseEnergy evaluateDouble(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                            std::vector<double>* grad = nullptr) const;

  // 9b minimization of one pose; false when the pose left the grid.
  // out holds the float pose, angles and energies of the last accepted point.
  bool minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats) const;

  // All of 9b for one conformation pair: L1poss -> 8.5 -> minimization.
  // Appends the docking samples to out (in L1poss order).
  Dock9bStats dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream, std::vector<DockSample>& out) const;
  // Same, minimizing the L2poss in chunks of chunkSize (in L1poss order); after each
  // chunk, sink receives the docking samples of the chunk.
  Dock9bStats dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream, std::size_t chunkSize,
                       const std::function<void(std::vector<DockSample>&)>& sink) const;

  long long l1Count(const Conf& c) const;  // L1poss per conformation pair

 private:
  template <class T>
  PoseEnergy evaluateT(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                       std::vector<double>* grad) const;

  const EnergySystem& sys_;
  const StamField& field_;
  Dock9bParams p_;
  std::vector<int> probe_;  // coordinate indices of the probe atoms
  std::vector<int> map_;    // type map per probe atom
  std::vector<double> charge_;
  std::vector<char> heavy_;
  std::array<int, 3> posN_{};  // positions per axis
  std::vector<Quat> hopf_;
};

// PRNG stream of a 9b conformation pair.
std::uint64_t dockStream(std::uint64_t round, std::uint64_t mobm, std::uint64_t stamSample, std::uint64_t mobmSample);

}  // namespace od
