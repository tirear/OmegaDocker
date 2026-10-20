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

#include "dock_kernel.hpp"
#include "energy.hpp"
#include "fft_dock.hpp"
#include "grid.hpp"
#include "minimize.hpp"
#include "orient.hpp"

namespace od {

// 9b docking of one MobM against one StaM conformation (spec 8.4–8.6).
// CPU reference implementation. Coordinates, grid interpolation and per-atom /
// per-pair terms are single precision; sums, L-BFGS and pose updates are double
// (spec 12.1).

enum class Posm { Posm1, Posm2, Posm3, Posm4, Posm2Fft };

class PoseEvaluator;

struct Dock9bParams {
  Posm posm = Posm::Posm1;
  double posStep = 1.0;          // Å
  int hopfLevel = 3;             // posm2, posm2_fft orientations (rotStep)
  int fftTopPerRot = 1;          // posm2_fft (spec 8.4.3)
  long long posm1SampleNum = 1000;
  // posm3 (spec 8.4.4)
  double sepaShellMin = 2.5, sepaShellMax = 6.0;  // Å
  double sepaPotPercentile = 5.0;                 // %
  std::optional<double> sepaPotThreshold;         // kcal/mol/e; replaces the percentile
  double sepaChargeThreshold = 0.3;               // e
  double sepaDistTolerance = 1.0;                 // Å
  double sepaPosStep = 2.0;                       // Å, posm3 translation spacing
  double sepaPosRange = 2.0;                      // Å
  // posm3 sampling (dock9b.sepaSampling): pair alignments, single placements, triple alignments.
  bool sepaPair = true, sepaSingle = false, sepaTriple = true;
  double rotStep = 7.5;                           // degrees, posm3 rotation about the pair axis
  Vec3 boxCenter{}, boxSize{};
  int maxOverlapAtomNum = 1;
  double maxNearDistance = 6.0;
  double clashDistance = 1.8;
  double clashBase = 10000.0;
  double clashK = 1000.0;
  // Minimization (spec 12)
  VecMinParams minimizer;
  double scaleTrans = 1, scaleRot = 1, scaleTors = 1;
  double maxTransStep = 0.5;     // Å
  double maxAngleStep = 10.0;    // degrees
  std::uint64_t seed = 0;
  PoseEvaluator* evaluator = nullptr;  // minimization energies; null: MobmDocker::evaluate on the CPU
};

// MobM placement: world x = R(q) (y - center) + t, where y are the coordinates
// of the conformation (torsion-tree frame) and center its heavy-atom centre.
// The energy is evaluated at float(t), float(q).
struct Pose {
  Vec3 t{};
  Quat q{};
};

// posm3: a StaM region of large |electrostatic potential| (spec 8.4.4).
struct StamSepa {
  Vec3 xyz{};     // centroid of the grid points
  int sign = 0;   // sign of the potential
  int points = 0;
};
struct StamSepaInfo {
  long long shellPoints = 0;  // grid points in the shell
  long long selected = 0;     // grid points passing the |potential| cut
  double potCut = 0;          // |potential| cut (kcal/mol/e)
};
// Regions of one StaM conformation; neighbours: 6 (face-sharing grid points) or 26.
std::vector<StamSepa> stamSepas(const StamField& field, const Dock9bParams& p, StamSepaInfo* info = nullptr,
                                int neighbours = 6);

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
  long long sepaPairs = 0;    // posm3: matched (StaM pair, MobM pair) alignments
  long long sepaSingles = 0;  // posm3: (StaM region, MobM atom) placements
  long long sepaTriples = 0;  // posm3: matched (StaM triple, MobM triple) alignments
  Dock9bStats& operator+=(const Dock9bStats& o);
};

class MobmDocker {
 public:
  // sys: energy system of the single MobM (parts = {mob}); field: the StaM conformation (spec 8.1, 8.2);
  // sepas: the regions of the StaM conformation (posm3, posm4); fft: its spectra (posm2_fft).
  MobmDocker(const EnergySystem& sys, const StamField& field, const Dock9bParams& p,
             const std::vector<StamSepa>* sepas = nullptr, const FftStam* fft = nullptr);

  // Lattice indices of the positions keeping every probe atom kPosMargin inside
  // the box for one orientation (spec 8.4); empty when lo > hi on an axis.
  static constexpr double kPosMargin = 0.01;  // Å, absorbs float rounding
  static constexpr double kTripleMinHeight = 0.5;  // Å, posm3 triples flatter than this are skipped
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
    // posm3: matched pairs (StaM regions s1, s2; probe atoms m1, m2), single placements (region, probe atom).
    struct PairMatch {
      int s1, s2, m1, m2;
    };
    // Triple: regions s, probe atoms m and the rotation superposing the m triangle on the s triangle.
    struct TripleMatch {
      std::array<int, 3> s, m;
      Quat q;
    };
    std::vector<PairMatch> pairs;
    std::vector<std::array<int, 2>> singles;
    std::vector<TripleMatch> triples;
    long long pairPoses = 0, singlePoses = 0, triplePoses = 0;
    long long sepaPoses() const { return pairPoses + singlePoses + triplePoses; }
    // posm2_fft: the kept translations, orientation-major, ascending FFT energy per orientation.
    std::vector<Pose> fftPoses;
    std::vector<float> fftEnergy;  // FFT (clamped-map) grid energy of each kept pose
  };
  Conf conformation(std::span<const double> anglesDeg) const;
  // posm2_fft (spec 15.1.22): conformation(anglesDeg) of each docker of dks (same MobM, parameters
  // and grid geometry, e.g. different StaM conformations); the probe spectra of each orientation
  // are computed once for all of them.
  static std::vector<Conf> conformations(std::span<const MobmDocker* const> dks, std::span<const double> anglesDeg);

  PosRange posRange(const Conf& c, const Quat& q) const;

  // L1pos k (0 <= k < l1Count(c)); stream: PRNG stream of the pair (posm1).
  // nullopt (posm1) when the sampled orientation does not fit in the box.
  std::optional<Pose> l1Pose(const Conf& c, std::uint64_t stream, long long k) const;
  // 8.5 quick filter; updates the counters of stats.
  bool quickFilter(const Conf& c, const Pose& pose, Dock9bStats& stats) const;
  // quickFilter of a pose whose probe atoms are all inside the box (fitsBox): decided at the
  // overlap limit (spec 15.1.21).
  bool quickFilterInBox(const Conf& c, const Pose& pose, Dock9bStats& stats) const;

  // posm3 L1pos k (k < c.sepaPoses()): rotation r = k / posm3Translations(), translation
  // k % posm3Translations(); the rotation gives the orientation, the region anchor and the
  // rotated conformation anchor rm, pose.t = anchorS + delta - rm (spec 8.4.4, 15.1.21).
  struct Posm3Rot {
    Quat q;
    Vec3 anchorS{}, rm{};
  };
  int posm3Translations() const { return transN_ * transN_ * transN_; }
  Posm3Rot posm3Rot(const Conf& c, long long r) const;
  Vec3 posm3Delta(int tr) const;
  Pose posm3Pose(const Posm3Rot& rot, int tr) const;
  // Box bounds of the L1pos box check (kPosMargin inside the box, float).
  void fitBounds(float lo[3], float hi[3]) const;

  // Energy of a pose; grad (optional): d/d(t, rotation vector, torsions in rad).
  struct PoseEnergy {
    EnergyTerms intra;
    double inter = 0, interClash = 0;
    int outside = 0;  // probe atoms outside the grid
    double objective() const { return inter + interClash + intra.total(); }
    double ranking() const { return inter + interClash + intra.filterEnergy(); }
  };
  // Energy of a pose at float(pose) and float(anglesDeg); grad: 6 + torsions, as evaluate.
  using PoseEval = std::function<PoseEnergy(const Pose&, std::span<const double> anglesDeg, std::vector<double>* grad)>;
  // Single precision at float(pose), float(anglesDeg) (spec 8.6).
  PoseEnergy evaluate(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                      std::vector<double>* grad = nullptr) const;
  // Double-precision reference of evaluate (tests).
  PoseEnergy evaluateDouble(const Conf& c, const Pose& pose, std::span<const double> anglesDeg,
                            std::vector<double>* grad = nullptr) const;

  // 9b minimization of one pose; false when the pose left the grid.
  // out holds the float pose, angles and energies of the last accepted point.
  bool minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats) const;
  // Same with the energies of ev.
  bool minimize(const Conf& c, const Pose& start, DockSample& out, Dock9bStats& stats, const PoseEval& ev) const;

  // All of 9b for one conformation pair: L1poss -> 8.5 -> minimization.
  // Appends the docking samples to out (in L1poss order).
  Dock9bStats dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream, std::vector<DockSample>& out) const;
  // Same, minimizing the L2poss in chunks of chunkSize (in L1poss order); after each
  // chunk, sink receives the docking samples of the chunk. posm1Only (posm4): only the posm1 part.
  Dock9bStats dockPair(std::span<const double> mobAnglesDeg, std::uint64_t stream, std::size_t chunkSize,
                       const std::function<void(std::vector<DockSample>&)>& sink, bool posm1Only = false) const;

  // Conformation, L1pos and quick-filter statistics of a pair.
  struct Prepared {
    Conf c;
    std::vector<long long> l2;
    Dock9bStats stats;
  };
  // Several conformation pairs: each job as dockPair(mobAnglesDeg, stream, chunkSize, sink,
  // posm1Only), the sinks receiving the same chunks in job order. With an evaluator the
  // L2poss of all jobs are minimized together (spec 15.1.4).
  struct PairJob {
    std::span<const double> mobAnglesDeg;
    std::uint64_t stream = 0;
    bool posm1Only = false;
    std::function<void(std::vector<DockSample>&)> sink;
    Dock9bStats stats;  // result
    std::optional<Prepared> prepared;  // set by prepareJobs; used (and cleared) by dockPairs
    std::optional<Conf> conf;          // optional: conformation(mobAnglesDeg), e.g. from conformations()
  };
  void dockPairs(std::vector<PairJob>& jobs, std::size_t chunkSize) const;
  // L1poss and 8.5 quick filter of the jobs ahead of dockPairs (e.g. to count their L2poss).
  void prepareJobs(std::vector<PairJob>& jobs) const;
  // Pair jobs of several dockers of the same MobM and parameters (e.g. different StaM
  // conformations): as dockPairs of each group in turn, the sinks receiving the same chunks in
  // group and job order. With an evaluator the L2poss of all groups are minimized together, at
  // most evaluator->fields() groups (spec 15.1.16).
  struct PairGroup {
    const MobmDocker* dk;
    std::vector<PairJob>* jobs;
  };
  static void dockPairs(std::span<const PairGroup> groups, std::size_t chunkSize);

  long long l1Count(const Conf& c) const;  // L1poss per conformation pair

  const EnergySystem& sys() const { return sys_; }
  const StamField& field() const { return field_; }
  const Dock9bParams& params() const { return p_; }
  const FftStam* fft() const { return fft_; }
  const std::vector<Quat>& hopf() const { return hopf_; }
  // Probe atoms for the flat pose energy (dock::evaluate, spec 15.1.4).
  std::vector<DockProbe> probes() const;
  DockView dockView(const std::vector<DockProbe>& probes, const FieldView& field) const;

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
  const std::vector<StamSepa>* sepas_ = nullptr;
  const FftStam* fft_ = nullptr;
  static void fftPoses(std::span<const MobmDocker* const> dks, std::span<Conf* const> cs);
  Conf baseConformation(std::span<const double> anglesDeg) const;  // conformation without the posm2_fft poses
  int axisAngles_ = 0;  // posm3: rotations about the pair axis
  int transN_ = 0;      // posm3: translations per axis
  std::optional<Pose> posm3Pose(const Conf& c, long long k) const;
  bool fitsBox(const Conf& c, const Pose& pose) const;
  // Extent of the rotated probe atoms (float, relative to the centre): fitsBox(box, pose) equals
  // fitsBox(c, pose) for any pose with this orientation (float addition is monotonic).
  struct RotBox {
    float lo[3], hi[3];
  };
  RotBox rotBox(const Conf& c, const Quat& q) const;
  template <bool InBox>
  bool quickFilterT(const Conf& c, const Pose& pose, Dock9bStats& stats) const;
  bool fitsBox(const RotBox& b, const Pose& pose) const;
  Prepared prepare(std::span<const double> mobAnglesDeg, std::uint64_t stream, bool posm1Only,
                   std::optional<Conf> conf = std::nullopt) const;
  Dock9bStats dockPrepared(Prepared pr, std::uint64_t stream, std::size_t chunkSize,
                           const std::function<void(std::vector<DockSample>&)>& sink) const;
  struct BatchItem {
    const MobmDocker* dk;
    int field;  // index of dk in the bound dockers
    const Conf* conf;
    std::uint64_t stream;
    long long k;  // L1pos
  };
  static void minimizeBatched(PoseEvaluator& pev, std::span<const MobmDocker* const> dks, std::span<const BatchItem> items,
                              std::vector<DockSample>& res, std::vector<char>& ok, std::vector<Dock9bStats>& stats);
};

// A pose energy request of the batched 9b minimization (spec 15.1.4).
struct PoseRequest {
  Pose pose;                 // float values
  Vec3 center{};             // the Pose frame
  const double* anglesDeg;   // torsion angles (float values)
  int field = 0;             // bound docker (PoseEvaluator::bindAll)
};

// Batched 9b pose energies, e.g. on a GPU (spec 15.1.4).
class PoseEvaluator {
 public:
  virtual ~PoseEvaluator() = default;
  virtual int maxBatch() const = 0;  // most poses per call
  // Selects the StaM conformation and probe atoms of dk for the following calls.
  virtual void bind(const MobmDocker& dk) = 0;
  // Several StaM conformations at once (spec 15.1.16): up to fields() dockers of the same MobM
  // bound together, PoseRequest::field selecting one of them.
  virtual int fields() const { return 1; }
  // Poses of dk evaluated at the same time (e.g. resident on a GPU); 0: unknown.
  virtual int residentPoses(const MobmDocker& /*dk*/) const { return 0; }
  virtual void bindAll(std::span<const MobmDocker* const> dks) { bind(*dks[0]); }
  // grad: n × (6 + torsions), as MobmDocker::evaluate.
  virtual void evaluate(int n, const PoseRequest* req, MobmDocker::PoseEnergy* out, double* grad) = 0;
  // Called when 9b of a round ends (e.g. to free work memory).
  virtual void idle() {}
  // Asynchronous form (spec 15.1.8): with lanes() > 1, up to lanes() evaluations are in flight;
  // submit copies the requests of a lane, complete waits and stores the results of its last submit.
  virtual int lanes() const { return 1; }
  virtual void submit(int /*lane*/, int /*n*/, const PoseRequest* /*req*/) {}
  virtual void complete(int /*lane*/, MobmDocker::PoseEnergy* /*out*/, double* /*grad*/) {}
  // Whole L-BFGS minimizations on the device (spec 15.1.26) from the poses and angles of req with
  // the parameters p; out and anglesOut (n × torsions). false: not supported (the minimization
  // then runs on the host with evaluate or submit and complete).
  virtual bool minimize(int /*n*/, const PoseRequest* /*req*/, const DockMinParams& /*p*/, DockMinResult* /*out*/,
                        float* /*anglesOut*/) {
    return false;
  }

  // Optional device versions of L1poss preparation (spec 15.1.5); the defaults leave it to the CPU.
  // 8.5 quick filter of the L1pos [k0, k1) of c: appends the passing ones to l2 (any order) and
  // adds to the outsideBox, overlap and far counters of stats. false: not supported.
  virtual bool quickFilter(const MobmDocker&, const MobmDocker::Conf&, std::uint64_t /*stream*/, long long /*k0*/,
                           long long /*k1*/, std::vector<long long>& /*l2*/, Dock9bStats& /*stats*/) {
    return false;
  }
  // posm2_fft: orientations per fftTop call against the spectra of the dockers dks (spec 15.1.22);
  // 0: not supported.
  virtual int fftBatch(std::span<const MobmDocker* const> /*dks*/) { return 0; }
  // n orientations: f (n × probes × 3, grid units, as fftCorrelate), map and q (per probe atom),
  // translation ranges lo, hi; out (dks × n × top, docker-major): fftSelect of each, entries
  // beyond the range count with index SIZE_MAX.
  virtual void fftTop(int /*n*/, int /*probes*/, const double* /*f*/, const int* /*map*/, const double* /*q*/,
                      const std::array<int, 3>* /*lo*/, const std::array<int, 3>* /*hi*/, int /*top*/,
                      std::pair<float, std::size_t>* /*out*/) {}
  // Called when the L1poss of a pair are prepared (e.g. to free the work memory of the above).
  virtual void prepared() {}
};

// PRNG stream of a 9b conformation pair.
std::uint64_t dockStream(std::uint64_t round, std::uint64_t mobm, std::uint64_t stamSample, std::uint64_t mobmSample);

}  // namespace od
