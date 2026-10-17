// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <algorithm>
#include <chrono>
#include <cmath>
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <format>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <string_view>
#include <unordered_map>

#include "ccd.hpp"
#include "checkpoint.hpp"
#include "closure.hpp"
#include "cli.hpp"
#include "log.hpp"
#include "read_input.hpp"
#include "conf9a.hpp"
#include "csgm1.hpp"
#include "csgm2.hpp"
#include "csgm3.hpp"
#include "dssp.hpp"
#include "energy.hpp"
#include "dock9b.hpp"
#include "dock9c.hpp"
#include "grid.hpp"
#include "grow.hpp"
#include "orient.hpp"
#include "feparams.hpp"
#include "score9d.hpp"
#include "shaking.hpp"
#include "write_struct.hpp"
#include "molconfig.hpp"
#include "topology.hpp"
#include "version.hpp"

namespace {

void printMolecule(const od::Molecule& m) {
  int netCharge = 0;
  double sumPartial = 0.0;
  for (const auto& a : m.atoms) {
    netCharge += a.formalCharge;
    sumPartial += a.partialCharge;
  }
  od::logInfo(std::format("  {:<40} atoms {:>7}  bonds {:>7}  residues {:>6}  formal charge {:>+4}  "
                          "sum partial charge {:>+8.3f}{}",
                          m.id, m.atoms.size(), m.bonds.size(), m.residues.size(), netCharge, sumPartial,
                          m.chargesFromFile ? " (file)" : " (Gasteiger)"));
}

od::TopologyParams topologyParams(const od::Options& opt) {
  od::TopologyParams p;
  p.box = od::Box{*opt.boxCenter, *opt.boxSize};
  p.stamFlexMargin = opt.stamFlexMargin;
  p.stamFixedResMovable = opt.stamFixedResMovable;
  return p;
}

od::Topology analyze(const od::Molecule& m, const od::Options& opt) {
  od::Topology t = od::analyzeTopology(m, topologyParams(opt));
  int rotatable = 0;
  for (const auto& tor : t.torsions) rotatable += tor.rotatable;
  od::logInfo(std::format("    backbone {} ({} units), torsions {}, rotatable by default {}",
                          od::backboneTypeName(t.backboneType), t.units.size(), t.torsions.size(), rotatable));
  return t;
}

// Spec 4.3, for omega angles set to fixed.
void checkOmegas(const od::Molecule& m, const od::MolSettings& s, const od::Options& opt) {
  const od::Topology& t = s.topology;
  std::vector<int> bad;
  for (int i : od::nonPlanarOmegas(t, opt.omegaTolerance))
    if (s.mode[i] == od::TorsionMode::Fixed) bad.push_back(i);
  if (bad.empty()) return;
  std::string list;
  for (int i : bad) {
    const od::Torsion& tor = t.torsions[i];
    const od::Atom& c = m.atoms[tor.atoms[1]];
    std::string where = c.residue >= 0 ? std::format("{}:{}{}:{}", m.residues[c.residue].info.chain,
                                                     m.residues[c.residue].info.seq,
                                                     m.residues[c.residue].info.iCode, m.residues[c.residue].info.name)
                                       : std::format("#{}", c.fileSerial);
    list += std::format("\n  {} omega {:.1f}", where, tor.angle);
  }
  if (!opt.acceptNonPlanarOmega)
    od::fail("{}: {} fixed omega angle(s) differ from 0 and 180 degrees by more than {} degrees "
             "(topology.acceptNonPlanarOmega = true continues):{}",
             m.id, bad.size(), opt.omegaTolerance, list);
  od::logWarn(std::format("{}: {} non-planar fixed omega angle(s):{}", m.id, bad.size(), list));
}

// Torsion tree, loops and added loop breaks (spec 4.7, 4.8).
void printTree(const od::Molecule& m, const od::MolSettings& s) {
  const od::Topology& t = s.topology;
  const od::TorsionTree& tree = s.tree;
  auto atomRef = [&](int a) {
    const od::Atom& x = m.atoms[a];
    if (x.residue < 0 || x.name.empty()) return std::format("#{}", x.fileSerial);
    const auto& r = m.residues[x.residue].info;
    return std::format("{}:{}{}:{}", r.chain, r.seq, r.iCode, x.name);
  };
  auto unitRef = [&](int u) {
    const int r = t.units[u].residue;
    return r >= 0 ? std::format("{}:{}{}", m.residues[r].info.chain, m.residues[r].info.seq, m.residues[r].info.iCode)
                  : std::format("unit {}", u + 1);
  };
  int depth = 0;
  std::vector<int> level(tree.torsions.size(), 0);
  for (std::size_t i = 0; i < tree.torsions.size(); ++i) {
    if (tree.torsions[i].parent >= 0) level[i] = level[tree.torsions[i].parent] + 1;
    depth = std::max(depth, level[i] + 1);
  }
  od::logInfo(std::format("    torsion tree: {} torsions, depth {}, anchor {}, closure pairs {}, distance restraints {}",
                          tree.torsions.size(), depth,
                          tree.anchors.size() == 1 ? atomRef(tree.anchors.front())
                                                   : std::format("{} fixed backbone atoms", tree.anchors.size()),
                          tree.virtuals.size(), tree.restraints.size()));
  for (const od::Loop& l : tree.loops) {
    std::string where = l.before >= 0 && l.after >= 0 ? "between fixed residues"
                        : l.before >= 0               ? "C-terminal, free end"
                        : l.after >= 0                ? "N-terminal, free end"
                                                      : "no fixed residue in the chain";
    std::string brk;
    if (l.breakBond >= 0)
      brk = std::format(", break [{}, {}]", atomRef(m.bonds[l.breakBond].a), atomRef(m.bonds[l.breakBond].b));
    od::logInfo(std::format("    loop {}-{} ({} residues, {}{})", unitRef(l.units.front()), unitRef(l.units.back()),
                            l.units.size(), where, brk));
  }
  for (const auto& c : s.chainBreaks)
    if (c.line == 0)
      od::logInfo(std::format("    added chainBreak [{}, {}] source = \"auto-loop\"", atomRef(m.bonds[c.bond].a),
                              atomRef(m.bonds[c.bond].b)));
}

od::EnergyParams energyParams(const od::Options& opt) {
  od::EnergyParams p;
  p.clashDistance = opt.clashDistance;
  p.clashBase = opt.EclashBase;
  p.clashK = opt.Kclash;
  p.closureG = opt.closureG;
  return p;
}

// 9a energy of the input conformation (spec 7.9).
void printInputEnergy(const std::string& name, const std::vector<od::EnergyPart>& parts, const od::Options& opt) {
  const od::EnergySystem sys(parts, energyParams(opt));
  int merged = 0, hd = 0;
  for (std::size_t k = 0; k < parts.size(); ++k)
    for (od::Ad4Type t : sys.typing(static_cast<int>(k)).type) {
      merged += t == od::Ad4Type::Merged;
      hd += t == od::Ad4Type::HD;
    }
  const auto e = sys.evaluate(sys.inputAngles());
  od::logInfo(std::format("  {:<40} torsions {:>5}  pairs {:>9}  merged H {:>5}  HD {:>5}", name, sys.torsionCount(),
                          sys.pairCount(), merged, hd));
  od::logInfo(std::format("    input: E {:.3f} (vdW {:.3f}  hbond {:.3f}  elec {:.3f}  desolv {:.3f}  clash {:.1f} "
                          "in {} pairs)  closure {:.3g}",
                          e.filterEnergy(), e.vdw, e.hbond, e.elec, e.desolv, e.clash, e.clashPairs, e.closure));
}

// Conformation samples of the StaM group or one MobM: the L2 pool, its 9a
// refills and the samples moved to a dcss (spec 7.3). Sample ids index samples().
class ConfSet {
 public:
  ConfSet(std::string name, std::vector<od::EnergyPart> parts, std::vector<const od::MolSettings*> settings,
          const od::MolConfig& cfg, const od::Options& opt, od::SampleSet set, std::uint64_t seed, std::uint64_t stream)
      : name_(std::move(name)), parts_(std::move(parts)), settings_(std::move(settings)), opt_(opt),
        csgm_(set == od::SampleSet::StaM ? opt.stamCsgm : opt.mobmCsgm), sys_(parts_, energyParams(opt)) {
    if (sys_.torsionCount() == 0) {
      samples_.push_back({});
      return;
    }
    od::Csgm1Params cp;
    cp.angleStep = csgm_.csgm1AngleStep;
    cp.ramaFloor = csgm_.ramaFloor;
    for (std::size_t k = 0; k < parts_.size(); ++k)
      sampler_.add(od::Csgm1Sampler(*parts_[k].mol, *settings_[k], cfg, cp));
    cisTrans_ = sampler_.cisTransMask();
    ga_ = {csgm_.tournamentRate / 100.0, csgm_.crossoverRate / 100.0, csgm_.mutationRate / 100.0, csgm_.absMaxDang};
    pe_ = {csgm_.peMaxOffset, csgm_.peCoverAngle, csgm_.peSampleNumPerLine};
    n_ = od::confSampleNums(opt, set);
    p_.minL2 = n_.minL2;
    p_.maxL1Batch = csgm_.maxL1BatchNum;
    p_.maxL1Total = csgm_.maxL1TotalNum;
    p_.Emax = opt.Emax;
    p_.maxBreakDist = opt.maxBreakDist;
    p_.min = {opt.confMinMaxStep, opt.confMinDeltaE, opt.confMaxDihedralStep};
    p_.seed = seed;
    p_.stream = stream;
    od::ClosureParams clp;
    clp.method = csgm_.closure == "kic" ? od::ClosureMethod::Kic
                 : csgm_.closure == "rcd" ? od::ClosureMethod::Rcd
                                          : od::ClosureMethod::None;
    clp.rcdTolerance = csgm_.rcdTolerance;
    clp.rcdMaxIter = csgm_.rcdMaxIter;
    clp.rcdMaxRestart = csgm_.rcdMaxRestart;
    clp.kicMaxTrial = csgm_.kicMaxTrial;
    clp.rcdClash = csgm_.rcdClash == "full" ? od::RcdClash::Full : csgm_.rcdClash == "repair" ? od::RcdClash::Repair : od::RcdClash::None;
    clp.rcdClashWeight = csgm_.rcdClashWeight;
    clp.clashDistance = opt.clashDistance;
    clp.growClose = csgm_.build == "growClose";
    std::vector<bool> freeMask(cisTrans_.size());
    for (std::size_t g = 0; g < freeMask.size(); ++g) freeMask[g] = !cisTrans_[g];
    closer_.emplace(sys_, parts_, freeMask, sampler_.segments(), clp);
    if (clp.method != od::ClosureMethod::None)
      for (const auto& line : closer_->describe()) od::logInfo(std::format("  {:<40} closure: {}", name_, line));
    const bool closing = closer_->active();
    if (csgm_.build != "whole") {
      // growClose without a loop closure: grow (spec 7.8).
      growClose_ = csgm_.build == "growClose" && closing;
      if (csgm_.build == "growClose" && !closing)
        od::logInfo(std::format("  {:<40} growClose: no loop closure; built as grow", name_));
      std::vector<const od::Topology*> topos;
      for (const auto* s : settings_) topos.push_back(&s->topology);
      od::GrowParams gp{csgm_.growMaxRetry, opt.clashDistance, opt.reachProtein, opt.reachNucleic};
      gp.maxBacktrack = csgm_.growMaxBacktrack;
      gp.reach = !closing;
      if (growClose_) gp.fixed = closer_->pathTorsions();
      grow_.emplace(sys_, parts_, topos, sampler_.segments(), gp);
      if (growClose_)
        od::logInfo(std::format("  {:<40} growClose: {} segments, {} applied before the loop closure", name_,
                                grow_->stepCount(), grow_->fixedStepCount()));
      else
        od::logInfo(std::format("  {:<40} grow: {} segments, {} reach checks", name_, grow_->stepCount(),
                                grow_->reachCheckCount()));
    }
  }
  ConfSet(const ConfSet&) = delete;
  ConfSet& operator=(const ConfSet&) = delete;

  const std::string& name() const { return name_; }
  bool movable() const { return sys_.torsionCount() > 0; }
  const od::EnergySystem& sys() const { return sys_; }
  const std::vector<od::EnergyPart>& parts() const { return parts_; }
  const std::vector<const od::MolSettings*>& settings() const { return settings_; }
  const od::ConfSample& sample(int id) const { return samples_[id]; }
  // Lowest 9a energy of the conformation samples made so far (L2 pool and dcss samples).
  float minEnergy() const {
    float e = std::numeric_limits<float>::infinity();
    for (const auto& s : samples_) e = std::min(e, s.energy);
    for (const auto& s : pool_) e = std::min(e, s.energy);
    return e;
  }
  // Angles of the n samples with the lowest 9a energy (L2 pool and dcss samples); the input conformation when not movable.
  std::vector<std::vector<float>> lowestAngles(int n) const {
    if (!movable()) return {{}};
    std::vector<const od::ConfSample*> v;
    for (const auto& s : samples_) v.push_back(&s);
    for (const auto& s : pool_) v.push_back(&s);
    std::stable_sort(v.begin(), v.end(), [](auto a, auto b) { return a->energy < b->energy; });
    std::vector<std::vector<float>> out;
    for (std::size_t k = 0; k < v.size() && k < static_cast<std::size_t>(n); ++k) out.push_back(v[k]->angles);
    return out;
  }
  // Checkpoint state (spec 16.3): random stream counter, L1 batch state, L2 pool and dcss samples.
  void save(od::BinWriter& w) const {
    w.put(counter_);
    w.put<std::uint8_t>(first_);
    w.put(nextL1_);
    for (const auto* v : {&pool_, &samples_}) {
      w.put<std::uint64_t>(v->size());
      for (const auto& s : *v) {
        w.putVec(s.angles);
        w.put(s.energy);
      }
    }
  }
  void load(od::BinReader& r) {
    counter_ = r.get<std::uint64_t>();
    first_ = r.get<std::uint8_t>() != 0;
    nextL1_ = r.get<long long>();
    for (auto* v : {&pool_, &samples_}) {
      const auto n = r.get<std::uint64_t>();
      v->clear();
      for (std::uint64_t k = 0; k < n; ++k) {
        od::ConfSample s;
        s.angles = r.getVec<float>();
        s.energy = r.get<float>();
        if (s.angles.size() != static_cast<std::size_t>(sys_.torsionCount()))
          od::fail("checkpoint: {}: {} torsions in a sample, {} expected", name_, s.angles.size(), sys_.torsionCount());
        v->push_back(std::move(s));
      }
    }
    if (samples_.empty() && !movable()) od::fail("checkpoint: {}: input conformation sample missing", name_);
  }
  std::size_t sampleCount() const { return samples_.size(); }
  int dockingNum() const { return movable() ? n_.docking : 1; }
  int eliteNum() const { return movable() ? n_.elite : 1; }

  // Moves need samples from the L2 pool into samples(), refilling the pool by
  // 9a first when it has fewer (spec 7.3); returns their ids. dcss: the samples
  // left in the dcss with their lowest ranking energy of the last round.
  std::vector<int> take(int need, const std::vector<std::pair<int, float>>& dcss) {
    if (static_cast<int>(pool_.size()) < need) refill(dcss);
    std::vector<int> ids;
    for (int k = 0; k < need; ++k) {
      ids.push_back(static_cast<int>(samples_.size()));
      samples_.push_back(std::move(pool_[k]));
    }
    pool_.erase(pool_.begin(), pool_.begin() + need);
    return ids;
  }

 private:
  // csgm1 while the dcss has fewer than eliteConfSampleNum samples, csgmGA (csgm2 or csgm3) otherwise.
  void refill(const std::vector<std::pair<int, float>>& dcss) {
    const bool ga = static_cast<int>(dcss.size()) >= n_.elite;
    const bool pe = ga && csgm_.pathExploring;
    std::vector<od::Csgm2Parent> parents;
    std::vector<const std::vector<float>*> ends;
    for (const auto& [id, e] : dcss) {
      parents.push_back({&samples_[id].angles, e});
      ends.push_back(&samples_[id].angles);
    }
    const od::SegmentSampler resample = [&](od::Rng& r, int seg, std::span<double> a) { sampler_.sampleSegment(r, seg, a); };
    const bool closing = closer_ && closer_->active();
    const od::ConfGenerator gen = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> a) {
      if (pe) {
        const int m = static_cast<int>(idx.local % pe_.perLine);
        od::csgm3Sample(od::Rng(p_.seed, p_.stream, idx.global - m), ends, idx.local / pe_.perLine, m, cisTrans_, pe_, a);
      } else if (ga) {
        od::csgm2Offspring(r, parents, cisTrans_, ga_, a);
      } else {
        sampler_.sample(r, a);
      }
      // csgm2 / csgm3: the torsions are given; building does not resample.
      const od::SegmentSampler& growResample = pe || ga ? od::SegmentSampler{} : resample;
      if (growClose_) {
        // growClose: closure paths as sampled, closure, then the rest built segment by segment (spec 7.8).
        if (!closer_->close(r, a, resample)) return od::GenResult::CloseFailed;
        if (!grow_->build(r, a, growResample).ok) return od::GenResult::SideGrowFailed;
        return od::GenResult::Built;
      }
      if (grow_ && !grow_->build(r, a, growResample).ok) return od::GenResult::GrowFailed;
      if (closing && !closer_->close(r, a, resample)) return od::GenResult::CloseFailed;
      return od::GenResult::Built;
    };
    const long long l1 = first_ ? n_.initL1 : nextL1_;
    first_ = false;
    p_.fixedL1 = pe ? od::csgm3SampleNum(ends.size(), pe_) : 0;
    const auto t0 = std::chrono::steady_clock::now();
    const auto st = od::fillL2Pool(sys_, gen, l1, p_, counter_, pool_);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    nextL1_ = st.nextL1;
    const long long minimized = st.l1 - st.growFailed - st.closeFailed - st.sideGrowFailed;
    od::logInfo(std::format("  9a {:<37} {} L1 {:>9} in {} batch(es)  -> L2 +{} = {}  (rejected: grow {}, closure {}, "
                            "grow after closure {}, E >= Emax {}, not closed {})  steps/sample {:.1f}  {:.2f} s",
                            name_, pe ? "csgm3" : ga ? "csgm2" : "csgm1", st.l1, st.batches, st.added, pool_.size(),
                            st.growFailed, st.closeFailed, st.sideGrowFailed, st.overEmax, st.notClosed,
                            minimized > 0 ? static_cast<double>(st.minSteps) / static_cast<double>(minimized) : 0.0, sec));
  }

  std::string name_;
  std::vector<od::EnergyPart> parts_;
  std::vector<const od::MolSettings*> settings_;
  const od::Options& opt_;
  const od::CsgmOpts& csgm_;
  od::EnergySystem sys_;
  od::Csgm1GroupSampler sampler_;
  std::vector<bool> cisTrans_;
  od::Csgm2Params ga_;
  od::Csgm3Params pe_;
  std::optional<od::GrowBuilder> grow_;
  bool growClose_ = false;
  std::optional<od::LoopCloser> closer_;
  od::ConfSampleNums n_;
  od::Fill9aParams p_;
  std::uint64_t counter_ = 0;
  bool first_ = true;
  long long nextL1_ = 0;
  std::vector<od::ConfSample> pool_;     // L2 pool, first in first out
  std::vector<od::ConfSample> samples_;  // samples moved to a dcss; the input conformation when not movable
};

od::Dock9bParams dockParams(const od::Options& opt, std::uint64_t seed) {
  od::Dock9bParams dp;
  dp.posm = opt.posm == "posm2"   ? od::Posm::Posm2
            : opt.posm == "posm2_fft" ? od::Posm::Posm2Fft
            : opt.posm == "posm3" ? od::Posm::Posm3
            : opt.posm == "posm4" ? od::Posm::Posm4
                                  : od::Posm::Posm1;
  dp.posStep = opt.posStep;
  dp.hopfLevel = od::hopfLevelForStep(opt.rotStep);
  dp.posm1SampleNum = opt.posm1SampleNum;
  dp.fftTopPerRot = static_cast<int>(opt.fftTopPerRot);
  dp.sepaShellMin = opt.sepaShellMin;
  dp.sepaShellMax = opt.sepaShellMax;
  dp.sepaPotPercentile = opt.sepaPotPercentile;
  dp.sepaPotThreshold = opt.sepaPotThreshold;
  dp.sepaChargeThreshold = opt.sepaChargeThreshold;
  dp.sepaDistTolerance = opt.sepaDistTolerance;
  dp.sepaPosStep = opt.sepaPosStep;
  dp.sepaPosRange = opt.sepaPosRange.value_or(opt.sepaPosStep);
  dp.sepaPair = opt.sepaSamplingPair;
  dp.sepaSingle = opt.sepaSamplingSingle;
  dp.sepaTriple = opt.sepaSamplingTriple;
  dp.rotStep = opt.rotStep;
  dp.boxCenter = *opt.boxCenter;
  dp.boxSize = *opt.boxSize;
  dp.maxOverlapAtomNum = opt.maxOverlapAtomNum;
  dp.maxNearDistance = opt.maxNearDistance;
  dp.clashDistance = opt.clashDistance;
  dp.clashBase = opt.EclashBase;
  dp.clashK = opt.Kclash;
  dp.lbfgs.maxSteps = opt.dockMinMaxStep;
  dp.lbfgs.convSteps = opt.dockMinConvSteps;
  dp.lbfgs.deltaE = opt.dockMinDeltaE;
  dp.lbfgs.memory = opt.lbfgsMemory;
  dp.scaleTrans = opt.scaleTrans;
  dp.scaleRot = opt.scaleRot;
  dp.scaleTors = opt.scaleTors;
  dp.maxTransStep = opt.maxTransStep;
  dp.maxAngleStep = opt.maxAngleStep;
  dp.seed = seed;
  return dp;
}

// Per-MobM state of the rounds (spec 6).
struct MobState {
  std::vector<int> stamDcss, mobDcss;  // sample ids
  std::unordered_map<std::uint64_t, float> pairLowest;  // (StaM id, MobM id) -> lowest docking energy, +inf if none
  std::unordered_map<int, float> stamLowest, mobLowest;  // dcss sample -> lowest ranking energy of the last round (csgm2 fitness)
  std::vector<od::RoundSample> top;
  int topNum = 0;
  std::size_t bufferNum = 0;  // sampleBufferNum
  std::unique_ptr<od::MobmHeavyAtoms> heavy;
  std::unique_ptr<od::SymmetryRmsd> sym;
};

std::uint64_t pairKey(int stam, int mobm) {
  return static_cast<std::uint64_t>(stam) << 32 | static_cast<std::uint32_t>(mobm);
}

// Maps are written sorted by key.
template <typename K>
void putMap(od::BinWriter& w, const std::unordered_map<K, float>& m) {
  std::vector<std::pair<K, float>> v(m.begin(), m.end());
  std::sort(v.begin(), v.end());
  std::vector<K> keys;
  std::vector<float> values;
  for (const auto& [k, e] : v) {
    keys.push_back(k);
    values.push_back(e);
  }
  w.putVec(keys);
  w.putVec(values);
}

template <typename K>
std::unordered_map<K, float> getMap(od::BinReader& r) {
  const auto keys = r.getVec<K>();
  const auto values = r.getVec<float>();
  if (keys.size() != values.size()) od::fail("checkpoint state: corrupt map");
  std::unordered_map<K, float> m;
  for (std::size_t k = 0; k < keys.size(); ++k) m.emplace(keys[k], values[k]);
  return m;
}

// Checkpoint state of one MobM (spec 16.3).
void saveMobState(od::BinWriter& w, const MobState& m) {
  w.putVec(m.stamDcss);
  w.putVec(m.mobDcss);
  putMap(w, m.pairLowest);
  putMap(w, m.stamLowest);
  putMap(w, m.mobLowest);
  w.put<std::uint64_t>(m.top.size());
  for (const auto& r : m.top) {
    w.put(r.stam);
    w.put(r.mobm);
    w.put(r.stamRel);
    w.put(r.s.pose);
    w.put(r.s.center);
    w.putVec(r.s.angles);
    w.put(r.s.energy);
    w.put(r.s.inter);
    w.put(r.s.interClash);
    w.put(r.s.intra);
    w.put(r.s.closureMaxDeviation);
    w.put(r.s.steps);
  }
}

void loadMobState(od::BinReader& r, MobState& m, std::size_t stamSamples, std::size_t mobSamples, int mobTorsions) {
  m.stamDcss = r.getVec<int>();
  m.mobDcss = r.getVec<int>();
  m.pairLowest = getMap<std::uint64_t>(r);
  m.stamLowest = getMap<int>(r);
  m.mobLowest = getMap<int>(r);
  const auto n = r.get<std::uint64_t>();
  m.top.clear();
  for (std::uint64_t k = 0; k < n; ++k) {
    od::RoundSample x;
    x.stam = r.get<int>();
    x.mobm = r.get<int>();
    x.stamRel = r.get<float>();
    x.s.pose = r.get<od::Pose>();
    x.s.center = r.get<od::Vec3>();
    x.s.angles = r.getVec<float>();
    x.s.energy = r.get<float>();
    x.s.inter = r.get<float>();
    x.s.interClash = r.get<float>();
    x.s.intra = r.get<float>();
    x.s.closureMaxDeviation = r.get<float>();
    x.s.steps = r.get<int>();
    m.top.push_back(std::move(x));
  }
  auto bad = [](int id, std::size_t count) { return id < 0 || static_cast<std::size_t>(id) >= count; };
  for (int id : m.stamDcss)
    if (bad(id, stamSamples)) od::fail("checkpoint state: StaM sample id {} out of range", id);
  for (int id : m.mobDcss)
    if (bad(id, mobSamples)) od::fail("checkpoint state: MobM sample id {} out of range", id);
  for (const auto& x : m.top)
    if (bad(x.stam, stamSamples) || bad(x.mobm, mobSamples) || x.s.angles.size() != static_cast<std::size_t>(mobTorsions))
      od::fail("checkpoint state: invalid topSampleSet entry");
}

// shift: added to RoundSample::rank so that the StaM term is relative to the current lowest StaM energy.
std::string energyList(const std::vector<od::RoundSample>& v, std::size_t n, float shift) {
  std::string s;
  for (std::size_t k = 0; k < v.size() && k < n; ++k)
    s += std::format("  {:.3f} (inter {:.3f}, intra {:.3f}, StaM {:.3f}; StaM {}, MobM {})", v[k].rank() + shift,
                     v[k].s.inter + v[k].s.interClash, v[k].s.intra, v[k].stamRel + shift, v[k].stam, v[k].mobm);
  return s;
}

// Output structure of the parts of sys at coordinates x (one entry per coordinate of sys).
std::vector<od::StructPart> structParts(const od::EnergySystem& sys, const std::vector<od::EnergyPart>& parts,
                                        const std::vector<od::Vec3f>& x) {
  std::vector<od::StructPart> out;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const int off = sys.coordOffset(static_cast<int>(k));
    const std::size_t n = parts[k].mol->atoms.size();
    out.push_back({parts[k].mol, std::vector<od::Vec3f>(x.begin() + off, x.begin() + off + static_cast<long>(n))});
  }
  return out;
}

// Directory name of a MobM: characters other than letters, digits and ._+- become _.
std::string safeName(std::string s) {
  for (char& c : s)
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '.' || c == '_' || c == '+' || c == '-')) c = '_';
  return s;
}

std::string joinNames(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : ",") + x;
  return s;
}

// 9d: recalculation by atom pairs, binding free energies, rankings and output (spec 10).
// AD4 terms (vdw, hbond, elec, desolv) as an array, in the order of od::kFeTermNames.
using Terms4 = std::array<double, od::FeParams::kTerms>;
Terms4 terms4(const od::EnergyTerms& e) { return {e.vdw, e.hbond, e.elec, e.desolv}; }
Terms4 terms4(const od::InterTerms& e) { return {e.vdw, e.hbond, e.elec, e.desolv}; }
double dot4(const Terms4& c, const Terms4& e) {
  double s = 0;
  for (int k = 0; k < od::FeParams::kTerms; ++k) s += c[k] * e[k];
  return s;
}

void run9d(ConfSet& stamSet, std::vector<std::unique_ptr<ConfSet>>& mobSets, const std::vector<MobState>& ms,
           const od::Options& opt, const od::FeParams& fp, float shift) {
  const od::EnergyParams ep = energyParams(opt);
  od::ShakeParams sp;
  sp.radius = fp.shakingRadius;
  sp.step = fp.shakingStep;
  sp.temperature = fp.temperature;
  sp.rrtScaling = fp.rrtScaling;
  const od::OutFormats fmt{opt.outPdb, opt.outCif, opt.outSdf, opt.outMol2};

  // A conformation evaluated alone: atoms, 9a terms (single precision) and, when needed, ln f of its atoms.
  struct Alone {
    std::vector<od::Vec3f> x;
    od::ShakeAtoms atoms;
    od::EnergyTerms terms;
    std::vector<double> lnF;
    double sumLnF = 0;
  };
  auto evalAlone = [&](const ConfSet& set, const std::vector<float>& angles, bool shake) {
    Alone a;
    a.x = set.sys().coordinates(angles);
    a.atoms = od::shakeAtoms(set.sys(), set.parts(), a.x);
    a.terms = set.sys().evaluate(std::span<const float>(angles), nullptr, nullptr, nullptr);
    if (shake) {
      a.lnF = od::shakeLnF(a.atoms, ep, sp);
      a.sumLnF = std::accumulate(a.lnF.begin(), a.lnF.end(), 0.0);
    }
    return a;
  };
  std::map<std::vector<float>, Alone> stamAlone;  // by angles
  auto stamOf = [&](const std::vector<float>& angles) -> const Alone& {
    auto it = stamAlone.find(angles);
    if (it == stamAlone.end()) it = stamAlone.emplace(angles, evalAlone(stamSet, angles, opt.feOd0100)).first;
    return it->second;
  };
  // Unbound state: the candidate with the lowest G = sum_t intra[t] E_intra,t + entropy (-TS).
  struct Unbound {
    double g = 0, mts = 0;
    Terms4 intra{};
  };
  auto unbound = [&](const ConfSet& set, bool stam) {
    Unbound u;
    u.g = std::numeric_limits<double>::infinity();
    const double cS = stam ? fp.entropyStam : fp.entropyMobm;
    for (const auto& ang : set.lowestAngles(fp.unboundCandidateNum)) {
      const Alone tmp = stam ? Alone{} : evalAlone(set, ang, true);
      const Alone& a = stam ? stamOf(ang) : tmp;
      const Terms4 e = terms4(a.terms);
      const double mts = od::minusTS(a.sumLnF, sp);
      const double g = dot4(fp.intra, e) + cS * mts;
      if (g < u.g) u = {g, mts, e};
    }
    return u;
  };
  Unbound ubStam;
  if (opt.feOd0100) {
    const auto t0 = std::chrono::steady_clock::now();
    ubStam = unbound(stamSet, true);
    od::logInfo(std::format("9d od_0100: G_unbound StaM group {:.3f}  ({} candidate(s), {:.2f} s)", ubStam.g,
                            stamSet.lowestAngles(fp.unboundCandidateNum).size(),
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
  }

  for (std::size_t i = 0; i < ms.size(); ++i) {
    const ConfSet& mob = *mobSets[i];
    const auto& top = ms[i].top;
    const auto t0 = std::chrono::steady_clock::now();
    Unbound ubMob;
    if (opt.feOd0100) ubMob = unbound(mob, false);
    const int nTors = mob.sys().torsionCount();
    const double unboundSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    struct Row {
      double e9b = 0, e9bGrid = 0, inter = 0, interClash = 0, mobIntra = 0, mobClash = 0, stamTerm = 0, stamIntra = 0,
             stamClash = 0;
      double dgAd = 0, dgOd = 0, gComplex = 0, mtsStam = 0, mtsMobm = 0;
      Terms4 inter4{}, dIntra4{};
      std::vector<od::Vec3f> world;
    };
    std::vector<Row> rows(top.size());
    for (std::size_t k = 0; k < top.size(); ++k) {
      const auto& r = top[k];
      const Alone& st = stamOf(stamSet.sample(r.stam).angles);
      Row& w = rows[k];
      const auto y = mob.sys().coordinates(std::span<const float>(r.s.angles));
      w.world = od::dockWorld(y, r.s);
      const od::ShakeAtoms mobAtoms = od::shakeAtoms(mob.sys(), mob.parts(), w.world);
      const od::EnergyTerms mt = mob.sys().evaluate(std::span<const float>(r.s.angles), nullptr, nullptr, nullptr);
      const od::InterTerms it = od::exactInter(st.atoms, mobAtoms, ep);
      w.inter = it.inter();
      w.interClash = it.clash;
      w.mobIntra = mt.intra();
      w.mobClash = mt.clash;
      w.stamTerm = r.stamRel + shift;
      w.stamIntra = st.terms.intra();
      w.stamClash = st.terms.clash;
      w.e9b = w.inter + w.interClash + w.mobIntra + w.mobClash + w.stamTerm;
      w.e9bGrid = r.rank() + shift;
      w.dgAd = w.inter + od::kAdgpuTorsionCoeff * nTors;
      if (opt.feOd0100) {
        const od::ComplexLnF l = od::complexSumLnF(st.atoms, st.lnF, mobAtoms, ep, sp);
        w.mtsStam = od::minusTS(l.stam, sp);
        w.mtsMobm = od::minusTS(l.mobm, sp);
        w.inter4 = terms4(it);
        const Terms4 si = terms4(st.terms), mi = terms4(mt);
        Terms4 intra4{};
        for (int t = 0; t < od::FeParams::kTerms; ++t) {
          intra4[t] = si[t] + mi[t];
          w.dIntra4[t] = intra4[t] - ubStam.intra[t] - ubMob.intra[t];
        }
        w.gComplex = dot4(fp.inter, w.inter4) + dot4(fp.intra, intra4) + fp.entropyStam * w.mtsStam +
                     fp.entropyMobm * w.mtsMobm;
        w.dgOd = fp.c0 + fp.tors * nTors + w.gComplex - ubStam.g - ubMob.g;
      }
    }
    const double evalSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() - unboundSec;

    // Rankings over the whole topSampleSet; output: union of the first outMinPerRank of each.
    struct Method {
      std::string name;
      double Row::*energy;
    };
    std::vector<Method> methods{{"9b", &Row::e9b}};
    if (opt.feAdgpu16) methods.push_back({"adgpu_16", &Row::dgAd});
    if (opt.feOd0100) methods.push_back({"od_0100", &Row::dgOd});
    std::vector<std::vector<std::size_t>> order(methods.size());
    std::vector<std::vector<int>> rankOf(methods.size(), std::vector<int>(top.size()));
    std::vector<char> out(top.size(), 0);
    for (std::size_t m = 0; m < methods.size(); ++m) {
      auto& o = order[m];
      o.resize(top.size());
      std::iota(o.begin(), o.end(), 0);
      std::stable_sort(o.begin(), o.end(), [&](auto a, auto b) { return rows[a].*methods[m].energy < rows[b].*methods[m].energy; });
      for (std::size_t z = 0; z < o.size(); ++z) {
        rankOf[m][o[z]] = static_cast<int>(z) + 1;
        if (z < static_cast<std::size_t>(opt.outMinPerRank)) out[o[z]] = 1;
      }
    }
    std::vector<int> pose(top.size(), 0);  // pose number, by the 9b ranking
    int poseNum = 0;
    for (std::size_t k : order[0])
      if (out[k]) pose[k] = ++poseNum;

    const std::string dir = opt.outPrefix + "/" + safeName(mob.name());
    std::filesystem::create_directories(dir);
    for (const auto& e : std::filesystem::directory_iterator(dir)) {  // output of an earlier run (e.g. before --resume)
      const std::string n = e.path().filename().string();
      if (e.is_regular_file() && (n.starts_with("pose_") || n.starts_with("stam_") || n == "summary.tsv" ||
                                  (n.starts_with("rank_") && n.ends_with(".tsv"))))
        std::filesystem::remove(e.path());
    }
    int fallback = 0;
    std::map<int, std::vector<std::string>> stamFiles;
    std::vector<std::vector<std::string>> poseFiles(top.size());
    for (std::size_t k : order[0]) {
      if (!out[k]) continue;
      const auto& r = top[k];
      if (!stamFiles.count(r.stam)) {
        const Alone& st = stamOf(stamSet.sample(r.stam).angles);
        bool fb = false;
        stamFiles[r.stam] = od::writeStructure(dir + std::format("/stam_{}", r.stam),
                                               structParts(stamSet.sys(), stamSet.parts(), st.x), fmt,
                                               std::format("StaM conformation sample {}", r.stam), fb);
        fallback += fb;
      }
      bool fb = false;
      poseFiles[k] = od::writeStructure(dir + std::format("/pose_{:04}", pose[k]),
                                        structParts(mob.sys(), mob.parts(), rows[k].world), fmt,
                                        std::format("{} pose {} (StaM sample {}, MobM sample {}) E_9b {:.3f}", mob.name(),
                                                    pose[k], r.stam, r.mobm, rows[k].e9b),
                                        fb);
      fallback += fb;
    }
    if (fallback > 0)
      od::logWarn(std::format("  {}: {} structure(s) exceed the PDB format limits and were written as mmCIF", mob.name(), fallback));

    auto openTsv = [&](const std::string& name) {
      std::ofstream f(dir + "/" + name);
      if (!f) od::fail("cannot write {}/{}", dir, name);
      return f;
    };
    {
      std::ofstream f = openTsv("summary.tsv");
      f << "pose\tmobm_files\tstam_files\tstam_sample\tmobm_sample";
      for (const auto& m : methods) f << "\trank_" << m.name;
      f << "\tE_9b\tE_9b_grid\tinter\tinter_clash\tmobm_intra\tmobm_clash\tstam_term\tstam_intra\tstam_clash\tclash";
      if (opt.feAdgpu16) f << "\tdG_adgpu_16";
      if (opt.feOd0100) {
        f << "\tdG_od_0100\tG_complex\tG_unbound_stam\tG_unbound_mobm";
        for (const char* t : od::kFeTermNames) f << "\tinter_" << t;
        for (const char* t : od::kFeTermNames) f << "\tdintra_" << t;
        f << "\tmTS_complex_stam\tmTS_complex_mobm\tmTS_unbound_stam\tmTS_unbound_mobm\tn_tors";
      }
      f << "\n";
      for (std::size_t k : order[0]) {
        if (!out[k]) continue;
        const Row& w = rows[k];
        f << std::format("{}\t{}\t{}\t{}\t{}", pose[k], joinNames(poseFiles[k]), joinNames(stamFiles[top[k].stam]),
                         top[k].stam, top[k].mobm);
        for (std::size_t m = 0; m < methods.size(); ++m) f << '\t' << rankOf[m][k];
        f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}", w.e9b, w.e9bGrid,
                         w.inter, w.interClash, w.mobIntra, w.mobClash, w.stamTerm, w.stamIntra, w.stamClash,
                         w.interClash + w.mobClash + w.stamClash > 0 ? "yes" : "no");
        if (opt.feAdgpu16) f << std::format("\t{:.4f}", w.dgAd);
        if (opt.feOd0100) {
          f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}", w.dgOd, w.gComplex, ubStam.g, ubMob.g);
          for (double v : w.inter4) f << std::format("\t{:.4f}", v);
          for (double v : w.dIntra4) f << std::format("\t{:.4f}", v);
          f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}", w.mtsStam, w.mtsMobm, ubStam.mts, ubMob.mts, nTors);
        }
        f << "\n";
      }
    }
    for (std::size_t m = 0; m < methods.size(); ++m) {
      std::ofstream f = openTsv("rank_" + methods[m].name + ".tsv");
      f << "rank\tpose\tmobm_files\tstam_files\tenergy\tclash\n";
      for (std::size_t k : order[m]) {
        if (!out[k]) continue;
        const Row& w = rows[k];
        f << std::format("{}\t{}\t{}\t{}\t{:.4f}\t{}\n", rankOf[m][k], pose[k], joinNames(poseFiles[k]),
                         joinNames(stamFiles[top[k].stam]), w.*methods[m].energy,
                         w.interClash + w.mobClash + w.stamClash > 0 ? "yes" : "no");
      }
    }

    od::logInfo(std::format("9d {:<38} recalculated {} sample(s) ({:.2f} s){}; output {} pose(s), {} StaM conformation(s) to {}",
                            mob.name(), top.size(), evalSec,
                            opt.feOd0100 ? std::format(", G_unbound MobM {:.3f} ({:.2f} s)", ubMob.g, unboundSec) : "",
                            poseNum, stamFiles.size(), dir));
    for (std::size_t m = 0; m < methods.size(); ++m) {
      std::string s;
      for (std::size_t z = 0; z < order[m].size() && z < 5; ++z) {
        const std::size_t k = order[m][z];
        s += std::format("  pose {} {:.3f}", pose[k], rows[k].*methods[m].energy);
      }
      od::logInfo(std::format("    lowest {:<9}{}", methods[m].name + ":", s));
    }
  }
}

// Checkpoints (spec 16).
struct Checkpointing {
  std::string params;                // params.toml of each checkpoint
  std::vector<od::InputHash> inputs;
  std::string resumeDir;             // empty: new run
  int resumeRound = 0;
  std::string stateSha256;
};

// topSampleSet of a checkpoint (spec 16.2): 9b ranking energies, no 9d recalculation.
void writeTopSnapshot(const std::string& dir, const ConfSet& stamSet, const std::vector<std::unique_ptr<ConfSet>>& mobSets,
                      const std::vector<MobState>& ms, const od::Options& opt, float shift) {
  const od::OutFormats fmt{opt.outPdb, opt.outCif, opt.outSdf, opt.outMol2};
  for (std::size_t i = 0; i < ms.size(); ++i) {
    const ConfSet& mob = *mobSets[i];
    const std::string d = dir + "/" + safeName(mob.name());
    std::filesystem::create_directories(d);
    std::ofstream f(d + "/topSampleSet.tsv");
    if (!f) od::fail("cannot write {}/topSampleSet.tsv", d);
    f << "rank\tpose_files\tstam_sample\tmobm_sample\tE_9b_grid\tinter\tinter_clash\tmobm_intra\tstam_term\n";
    std::map<int, bool> stamDone;
    for (std::size_t k = 0; k < ms[i].top.size(); ++k) {
      const auto& r = ms[i].top[k];
      std::string files;
      if (k < static_cast<std::size_t>(opt.outMinPerRank)) {
        bool fb = false;
        if (!stamDone[r.stam]) {
          stamDone[r.stam] = true;
          const auto& sa = stamSet.sample(r.stam).angles;
          od::writeStructure(d + std::format("/stam_{}", r.stam),
                             structParts(stamSet.sys(), stamSet.parts(), stamSet.sys().coordinates(std::span<const float>(sa))),
                             fmt, std::format("StaM conformation sample {}", r.stam), fb);
        }
        const auto y = mob.sys().coordinates(std::span<const float>(r.s.angles));
        files = joinNames(od::writeStructure(d + std::format("/top_{:04}", k + 1),
                                             structParts(mob.sys(), mob.parts(), od::dockWorld(y, r.s)), fmt,
                                             std::format("{} topSampleSet {} (StaM sample {}, MobM sample {}) E {:.3f}",
                                                         mob.name(), k + 1, r.stam, r.mobm, r.rank() + shift),
                                             fb));
      }
      f << std::format("{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n", k + 1, files.empty() ? "-" : files, r.stam,
                       r.mobm, r.rank() + shift, r.s.inter, r.s.interClash, r.s.intra, r.stamRel + shift);
    }
  }
}

void writeCheckpoint(int round, const ConfSet& stamSet, const std::vector<std::unique_ptr<ConfSet>>& mobSets,
                     const std::vector<MobState>& ms, float stamRef, const od::Options& opt, const Checkpointing& ck,
                     float shift) {
  const auto t0 = std::chrono::steady_clock::now();
  od::BinWriter w;
  od::putStateHeader(w);
  w.put<std::int32_t>(round);
  w.put<std::uint64_t>(mobSets.size());
  stamSet.save(w);
  for (const auto& m : mobSets) m->save(w);
  w.put(stamRef);
  for (const auto& m : ms) saveMobState(w, m);
  const od::Manifest man{od::kVersion, round, ck.inputs, od::sha256Hex(w.data())};
  const std::string dir = od::commitCheckpoint(
      opt.outPrefix, round, {{"state.bin", w.data()}, {"params.toml", ck.params}, {"manifest.toml", od::manifestText(man)}},
      [&](const std::string& tmp) { writeTopSnapshot(tmp, stamSet, mobSets, ms, opt, shift); });
  od::logInfo(std::format("  checkpoint {} ({:.1f} MB state)  {:.2f} s", dir, static_cast<double>(w.data().size()) / 1e6,
                          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
}

// Rounds of 9a refill, 9b docking and 9c selection (spec 6, 7.3, 8, 9).
void runRounds(ConfSet& stamSet, std::vector<std::unique_ptr<ConfSet>>& mobSets, const od::Options& opt,
               const od::FeParams& fp, std::uint64_t seed, const Checkpointing& ck) {
  std::vector<od::Ad4Type> probe;
  for (const auto& m : mobSets)
    for (od::Ad4Type t : m->sys().typing(0).type)
      if (t != od::Ad4Type::Merged && std::find(probe.begin(), probe.end(), t) == probe.end()) probe.push_back(t);
  const od::GridGeometry geom = od::boxGrid(*opt.boxCenter, *opt.boxSize, opt.gridSpacing);
  {
    const int level = od::hopfLevelForStep(opt.rotStep);
    od::logInfo(std::format("9b grids: {} x {} x {} points, spacing {} Å; orientations: Hopf level {}, {} "
                            "orientations, spacing {} deg (dock9b.rotStep {})",
                            geom.n[0], geom.n[1], geom.n[2], geom.spacing, level, 72LL << (3 * level),
                            od::hopfSpacingDeg(level), opt.rotStep));
  }
  const od::Dock9bParams dp = dockParams(opt, seed);
  const bool deterministic = dp.posm == od::Posm::Posm2 || dp.posm == od::Posm::Posm2Fft || dp.posm == od::Posm::Posm3;
  const bool usesFft = dp.posm == od::Posm::Posm2Fft;
  const bool usesSepas = dp.posm == od::Posm::Posm3 || dp.posm == od::Posm::Posm4;

  // grid.mode (spec 8.1); auto: delta when no StaM backbone torsion is rotatable.
  od::GridMode gridMode = od::GridMode::Full;
  bool backboneMoves = false;
  for (std::size_t k = 0; k < stamSet.parts().size(); ++k)
    for (const od::TreeTorsion& tt : stamSet.parts()[k].tree->torsions)
      backboneMoves |= stamSet.settings()[k]->topology.torsions[tt.torsion].backbone != od::BackboneAngle::None;
  if (opt.gridMode == "flexres") gridMode = od::GridMode::Flexres;
  else if (opt.gridMode == "delta") gridMode = od::GridMode::Delta;
  else if (opt.gridMode == "auto" && !backboneMoves) gridMode = od::GridMode::Delta;
  const std::vector<char> movedMask = od::movedGridAtoms(stamSet.sys(), stamSet.parts());
  const auto movedNum = std::count(movedMask.begin(), movedMask.end(), 1);
  od::logInfo(std::format("9b grid.mode {}{}: StaM grid atoms {}, moved by torsions {}, backbone {}", od::gridModeName(gridMode),
                          opt.gridMode == "auto" ? " (auto)" : "", movedMask.size(), movedNum,
                          backboneMoves ? "movable" : "fixed"));
  const double distCap = std::max({opt.maxNearDistance, opt.clashDistance,
                                   opt.posm == "posm3" || opt.posm == "posm4" ? opt.sepaShellMax : 0.0}) +
                         2 * opt.gridSpacing;
  // Grids built once and kept across rounds: full without StaM torsions (input conformation),
  // flexres (fixed atoms of the input conformation), delta (reference: input conformation).
  std::shared_ptr<const od::EnergyGrids> fixedGrids;
  std::shared_ptr<const od::FftStam> fixedFft;  // posm2_fft, StaM without torsions
  std::shared_ptr<const od::DistanceGrid> fixedDist;

  std::vector<MobState> ms(mobSets.size());
  for (std::size_t i = 0; i < ms.size(); ++i) {
    ms[i].topNum = opt.topSampleNum ? *opt.topSampleNum : stamSet.eliteNum() * mobSets[i]->eliteNum();
    ms[i].heavy = std::make_unique<od::MobmHeavyAtoms>(mobSets[i]->sys());
    ms[i].bufferNum = std::max<std::size_t>(opt.sampleBufferNum, static_cast<std::size_t>(std::max(4, od::midTopSampleFactor(opt) + 2)) *
                                                                       static_cast<std::size_t>(ms[i].topNum));
    const auto s0 = std::chrono::steady_clock::now();
    ms[i].sym = std::make_unique<od::SymmetryRmsd>(*mobSets[i]->parts()[0].mol, ms[i].heavy->atoms(), ms[i].heavy->types(),
                                                   static_cast<std::size_t>(opt.symMaxPerm));
    const double symSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - s0).count();
    od::logInfo(std::format("  {:<40} topSampleNum {}  sampleBufferNum {}  9d symmetry: {} ({:.3f} s)",
                            mobSets[i]->name(), ms[i].topNum, ms[i].bufferNum,
                            ms[i].sym->usesPermutations()
                                ? std::format("{} automorphism(s)", ms[i].sym->permutationCount())
                                : std::format("more than {} automorphisms, one-to-one pairing within each AD4 type",
                                              opt.symMaxPerm),
                            symSec));
  }

  float stamRef = std::numeric_limits<float>::infinity();  // fixed reference of RoundSample::stamRel (spec 9)
  int firstRound = 1;
  if (!ck.resumeDir.empty()) {
    const std::string path = od::checkpointDir(ck.resumeDir, ck.resumeRound) + "/state.bin";
    std::string data = od::readFile(path);
    if (od::sha256Hex(data) != ck.stateSha256) od::fail("{}: SHA-256 differs from manifest.toml (damaged checkpoint)", path);
    od::BinReader r(std::move(data));
    od::checkStateHeader(r);
    if (r.get<std::int32_t>() != ck.resumeRound || r.get<std::uint64_t>() != mobSets.size())
      od::fail("{}: round or MobM count differs from manifest.toml", path);
    stamSet.load(r);
    for (auto& m : mobSets) m->load(r);
    stamRef = r.get<float>();
    for (std::size_t i = 0; i < ms.size(); ++i)
      loadMobState(r, ms[i], stamSet.sampleCount(), mobSets[i]->sampleCount(), mobSets[i]->sys().torsionCount());
    if (!r.atEnd()) od::fail("{}: unexpected data at the end", path);
    firstRound = ck.resumeRound + 1;
    od::logInfo(std::format("Resumed from {} (round {} completed)", path, ck.resumeRound));
    if (firstRound > opt.maxRound)
      od::logInfo(std::format("  dock9c.maxRound {}: no further rounds", opt.maxRound));
  }
  for (int round = firstRound; round <= opt.maxRound; ++round) {
    od::logInfo(std::format("Round {} of {}:", round, opt.maxRound));
    // 9a: refill the dcss (spec 7.3).
    for (std::size_t i = 0; i < ms.size(); ++i) {
      auto fill = [](ConfSet& set, std::vector<int>& dcss, const std::unordered_map<int, float>& lowest) {
        if (!set.movable()) {
          dcss = {0};  // the input conformation (spec 7.1)
          return;
        }
        const int need = set.dockingNum() - static_cast<int>(dcss.size());
        if (need <= 0) return;
        std::vector<std::pair<int, float>> fit;
        for (int id : dcss) fit.emplace_back(id, lowest.at(id));
        const auto ids = set.take(need, fit);
        dcss.insert(dcss.end(), ids.begin(), ids.end());
      };
      fill(stamSet, ms[i].stamDcss, ms[i].stamLowest);
      fill(*mobSets[i], ms[i].mobDcss, ms[i].mobLowest);
    }
    if (stamSet.movable() && !std::isfinite(stamRef)) stamRef = stamSet.minEnergy();
    const float shift = stamSet.movable() ? stamRef - stamSet.minEnergy() : 0.0f;

    // 9b: one grid set per StaM conformation sample in a dcss, shared by all MobM.
    std::vector<int> stamIds;
    for (const auto& s : ms) stamIds.insert(stamIds.end(), s.stamDcss.begin(), s.stamDcss.end());
    std::sort(stamIds.begin(), stamIds.end());
    stamIds.erase(std::unique(stamIds.begin(), stamIds.end()), stamIds.end());
    std::vector<od::Dock9bStats> st(ms.size());
    std::vector<long long> docked(ms.size(), 0), skipped(ms.size(), 0), sampleNum(ms.size(), 0), partial(ms.size(), 0),
        posm3Docked(ms.size(), 0);
    std::vector<std::vector<od::RoundSample>> buffer(ms.size());
    std::vector<double> selectSec(ms.size(), 0);
    std::vector<float> droppedMin(ms.size(), std::numeric_limits<float>::infinity());  // --midKeepSkipped
    std::vector<std::size_t> poolMax(ms.size(), 0);
    // 9c steps 1-4 on the topSampleSet and the waiting docking samples (spec 9).
    auto select = [&](std::size_t i, bool last) {
      MobState& m = ms[i];
      std::vector<od::RoundSample> merged = std::move(m.top);
      merged.insert(merged.end(), std::make_move_iterator(buffer[i].begin()), std::make_move_iterator(buffer[i].end()));
      buffer[i].clear();
      const auto c0 = std::chrono::steady_clock::now();
      const int midNum = static_cast<int>(
          std::min<long long>(1LL * od::midTopSampleFactor(opt) * m.topNum, std::numeric_limits<int>::max()));
      if (last) {
        m.top = od::selectTopSamples(std::move(merged), *m.heavy, m.topNum, opt.diversityRmsd);
      } else if (opt.midKeepSkipped) {
        auto r = od::selectPool(std::move(merged), *m.heavy, static_cast<std::size_t>(midNum), opt.diversityRmsd, m.bufferNum);
        m.top = std::move(r.pool);
        droppedMin[i] = std::min(droppedMin[i], r.droppedMin);
        poolMax[i] = std::max(poolMax[i], m.top.size());
      } else {
        m.top = od::selectTopSamples(std::move(merged), *m.heavy, midNum, opt.diversityRmsd);
      }
      selectSec[i] += std::chrono::duration<double>(std::chrono::steady_clock::now() - c0).count();
    };
    double gridSec = 0, fftSec = 0;
    bool gridsReused = false;
    const auto d0 = std::chrono::steady_clock::now();
    long long movedSum = 0, blockSum = 0, blockTotal = 0;
    long long sepaShell = 0, sepaSelected = 0, sepaRegions = 0, sepaPos = 0;
    double sepaCutSum = 0;
    for (int sid : stamIds) {
      const auto g0 = std::chrono::steady_clock::now();
      auto stamAtoms = [&](bool input) {
        const auto& sa = stamSet.sample(sid).angles;
        std::vector<double> angles = !input && stamSet.movable() ? std::vector<double>(sa.begin(), sa.end())
                                                                 : stamSet.sys().inputAngles();
        std::vector<od::Vec3> x;
        stamSet.sys().evaluate(angles, nullptr, &x);
        return od::gridAtoms(stamSet.sys(), stamSet.parts(), x);
      };
      auto build = [&](const od::GridAtoms& a) {
        return std::pair{std::make_shared<const od::EnergyGrids>(od::buildEnergyGrids(geom, a, probe, energyParams(opt))),
                         std::make_shared<const od::DistanceGrid>(od::buildDistanceGrid(geom, a, distCap))};
      };
      const bool firstBuild = !fixedGrids;
      od::StamField field;
      if (gridMode == od::GridMode::Flexres) {
        if (!fixedGrids) std::tie(fixedGrids, fixedDist) = build(od::subset(stamAtoms(true), movedMask, false));
        od::GridAtoms moved = od::subset(stamAtoms(false), movedMask, true);
        movedSum += static_cast<long long>(moved.size());
        field = od::StamField::flexres(fixedGrids, fixedDist, std::move(moved), energyParams(opt));
      } else if (!stamSet.movable()) {
        if (!fixedGrids) std::tie(fixedGrids, fixedDist) = build(stamAtoms(true));
        field = od::StamField::full(fixedGrids, fixedDist);
      } else if (gridMode == od::GridMode::Delta) {
        if (!fixedGrids) std::tie(fixedGrids, fixedDist) = build(stamAtoms(true));
        const auto [g, d] = build(stamAtoms(false));
        field = od::StamField::delta(fixedGrids, fixedDist, *g, *d, opt.gridDeltaThreshold);
        blockSum += static_cast<long long>(field.deltaBlockCount());
        blockTotal += static_cast<long long>(field.blockCount());
      } else {
        const auto [g, d] = build(stamAtoms(false));
        field = od::StamField::full(g, d);
      }
      if (!firstBuild && (gridMode == od::GridMode::Flexres || !stamSet.movable())) gridsReused = true;
      std::vector<od::StamSepa> sepas;
      if (usesSepas) {
        od::StamSepaInfo si;
        sepas = od::stamSepas(field, dp, &si);
        int pos = 0;
        for (const auto& s : sepas) pos += s.sign > 0;
        sepaShell += si.shellPoints;
        sepaSelected += si.selected;
        sepaRegions += static_cast<long long>(sepas.size());
        sepaPos += pos;
        sepaCutSum += si.potCut;
      }
      std::shared_ptr<const od::FftStam> fft;
      if (usesFft) {
        const auto f0 = std::chrono::steady_clock::now();
        if (!stamSet.movable() && gridMode != od::GridMode::Flexres && fixedFft) {
          fft = fixedFft;
        } else {
          fft = std::make_shared<const od::FftStam>(od::buildFftStam(field, opt.fftClamp));
          if (!stamSet.movable() && gridMode != od::GridMode::Flexres) fixedFft = fft;
        }
        fftSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - f0).count();
      }
      gridSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - g0).count();
      for (std::size_t i = 0; i < ms.size(); ++i) {
        if (std::find(ms[i].stamDcss.begin(), ms[i].stamDcss.end(), sid) == ms[i].stamDcss.end()) continue;
        const ConfSet& mob = *mobSets[i];
        const od::MobmDocker dk(mob.sys(), field, dp, &sepas, fft.get());
        for (int mid : ms[i].mobDcss) {
          const std::uint64_t key = pairKey(sid, mid);
          if (deterministic && ms[i].pairLowest.contains(key)) {
            ++skipped[i];
            continue;
          }
          // posm4: the posm3 part of a pair docked before is not redone (spec 9).
          const bool posm1Only = dp.posm == od::Posm::Posm4 && ms[i].pairLowest.contains(key);
          const auto& ma = mob.sample(mid).angles;
          const std::vector<double> mobAngles =
              mob.movable() ? std::vector<double>(ma.begin(), ma.end()) : mob.sys().inputAngles();
          const float stamRel = stamSet.movable() ? stamSet.sample(sid).energy - stamRef : 0.0f;
          float lowest = std::numeric_limits<float>::infinity();
          st[i] += dk.dockPair(mobAngles, od::dockStream(round, i, sid, mid), ms[i].bufferNum,
                               [&](std::vector<od::DockSample>& chunk) {
                                 sampleNum[i] += static_cast<long long>(chunk.size());
                                 for (auto& s : chunk) {
                                   od::RoundSample r{sid, mid, std::move(s), stamRel};
                                   lowest = std::min(lowest, r.rank());
                                   buffer[i].push_back(std::move(r));
                                 }
                                 if (buffer[i].size() >= ms[i].bufferNum) {
                                   select(i, false);
                                   ++partial[i];
                                 }
                               },
                               posm1Only);
          ++docked[i];
          if (usesSepas && !posm1Only) ++posm3Docked[i];
          ms[i].pairLowest[key] = posm1Only ? std::min(lowest, ms[i].pairLowest[key]) : lowest;
        }
      }
    }
    double partialSec = 0;
    for (double s : selectSec) partialSec += s;
    const double dockSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - d0).count() - gridSec - partialSec;
    {
      std::string what = std::format("  9b grids for {} StaM conformation(s): {:.2f} s", stamIds.size(), gridSec);
      if (gridsReused)
        what += gridMode == od::GridMode::Flexres ? " (grids of the fixed atoms reused)" : " (grids of the input conformation reused)";
      if (gridMode == od::GridMode::Flexres && !stamIds.empty())
        what += std::format("; moved atoms per conformation {:.1f}", static_cast<double>(movedSum) / static_cast<double>(stamIds.size()));
      if (blockTotal > 0)
        what += std::format("; delta blocks stored {} of {} ({:.1f}%)", blockSum, blockTotal, 100.0 * static_cast<double>(blockSum) / static_cast<double>(blockTotal));
      if (usesFft) what += std::format("; posm2_fft StaM spectra {:.2f} s", fftSec);
      od::logInfo(what);
      if (usesSepas && !stamIds.empty()) {
        const double k = static_cast<double>(stamIds.size());
        od::logInfo(std::format("  9b posm3 StaM regions per conformation: shell points {:.0f}, selected {:.0f} "
                                "(|potential| >= {:.3f} kcal/mol/e), regions {:.1f} (+ {:.1f}, - {:.1f})",
                                static_cast<double>(sepaShell) / k, static_cast<double>(sepaSelected) / k, sepaCutSum / k,
                                static_cast<double>(sepaRegions) / k, static_cast<double>(sepaPos) / k,
                                static_cast<double>(sepaRegions - sepaPos) / k));
      }
    }

    // 9c (spec 9).
    for (std::size_t i = 0; i < ms.size(); ++i) {
      MobState& m = ms[i];
      const auto& s = st[i];
      od::logInfo(std::format("  9b {:<37} pairs docked {}  skipped {}  L1poss {}  removed (outside box {}, overlap {}, "
                              "far {})  L2poss {}  left grid {}  docking samples {}  steps/pose {:.1f}",
                              mobSets[i]->name(), docked[i], skipped[i], s.l1, s.outsideBox, s.overlap, s.far, s.l2,
                              s.leftGrid, sampleNum[i],
                              s.l2 > 0 ? static_cast<double>(s.steps) / static_cast<double>(s.l2) : 0.0));
      if (posm3Docked[i] > 0)
        od::logInfo(std::format("  9b {:<37} posm3 pairs {}: per pair matched pair alignments {:.1f}, single placements {:.1f}, "
                                "triple alignments {:.1f}",
                                mobSets[i]->name(), posm3Docked[i], static_cast<double>(s.sepaPairs) / posm3Docked[i],
                                static_cast<double>(s.sepaSingles) / posm3Docked[i],
                                static_cast<double>(s.sepaTriples) / posm3Docked[i]));
      if (s.l1 > 0 && s.l1 == s.outsideBox)
        od::logWarn(std::format("  {}: no sampled orientation fits in the docking box; consider a larger box.size",
                                mobSets[i]->name()));
      std::unordered_map<int, float> stamLowest, mobLowest;
      for (int a : m.stamDcss)
        for (int b : m.mobDcss) {
          const float e = m.pairLowest.at(pairKey(a, b));
          if (!std::isfinite(e)) continue;
          auto upd = [e](std::unordered_map<int, float>& mp, int id) {
            const auto [it, added] = mp.try_emplace(id, e);
            if (!added) it->second = std::min(it->second, e);
          };
          upd(stamLowest, a);
          upd(mobLowest, b);
        }
      select(i, true);
      const float eTop = od::topMaxEnergy(m.top);
      const int rmStam = stamSet.movable() ? od::filterDcss(m.stamDcss, stamLowest, eTop) : 0;
      const int rmMob = mobSets[i]->movable() ? od::filterDcss(m.mobDcss, mobLowest, eTop) : 0;
      m.stamLowest = std::move(stamLowest);
      m.mobLowest = std::move(mobLowest);
      od::logInfo(std::format("  9c {:<37} topSampleSet {}/{}  E_top_max {:.3f}  dcss removed: StaM {}, MobM {}  "
                              "9c runs {}  {:.2f} s",
                              mobSets[i]->name(), m.top.size(), m.topNum, eTop + shift, rmStam, rmMob, partial[i] + 1,
                              selectSec[i]));
      if (opt.midKeepSkipped) {
        const bool same = std::isinf(droppedMin[i]) ||
                          (m.top.size() == static_cast<std::size_t>(m.topNum) && eTop < droppedMin[i]);
        od::logInfo(std::format("    midKeepSkipped: largest pool {}  lowest removed energy {:.3f}  same as a single 9c: {}",
                                poolMax[i], droppedMin[i] + shift, same ? "yes" : "not guaranteed"));
      }
      if (!m.top.empty()) od::logInfo("    lowest energies:" + energyList(m.top, 5, shift));
    }
    od::logInfo(std::format("  9b docking: {:.2f} s", dockSec));
    writeCheckpoint(round, stamSet, mobSets, ms, stamRef, opt, ck, shift);
  }

  for (std::size_t i = 0; i < ms.size(); ++i) {
    const std::size_t before = ms[i].top.size();
    const auto f0 = std::chrono::steady_clock::now();
    ms[i].top = od::symmetryFilter(std::move(ms[i].top), *ms[i].heavy, *ms[i].sym, opt.diversityRmsd);
    od::logInfo(std::format("9d {:<38} symmetry-corrected RMSD < {} Å: removed {} of {}  {:.2f} s", mobSets[i]->name(),
                            opt.diversityRmsd, before - ms[i].top.size(), before,
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - f0).count()));
  }
  const float shift = stamSet.movable() && std::isfinite(stamRef) ? stamRef - stamSet.minEnergy() : 0.0f;
  run9d(stamSet, mobSets, ms, opt, fp, shift);
  od::logInfo("Final topSampleSet:");
  for (std::size_t i = 0; i < ms.size(); ++i) {
    od::logInfo(std::format("  {:<40} {} sample(s)", mobSets[i]->name(), ms[i].top.size()));
    for (std::size_t k = 0; k < ms[i].top.size(); ++k) {
      const auto& r = ms[i].top[k];
      od::logInfo(std::format("    {:>4}  E {:>10.3f}  E_9b {:>10.3f}  inter {:>10.3f}  intra {:>10.3f}  StaM intra {:>10.3f}  "
                              "StaM {:>5}  MobM {:>5}  t ({:.3f}, {:.3f}, {:.3f})",
                              k + 1, r.rank() + shift, r.s.energy, r.s.inter + r.s.interClash, r.s.intra, r.stamRel + shift,
                              r.stam, r.mobm, r.s.pose.t[0],
                              r.s.pose.t[1], r.s.pose.t[2]));
    }
  }
}

// DSSP residue counts, e.g. "H 12  E 10  . 5".
void printDssp(const od::Topology& t) {
  if (t.dsspSS.empty()) return;
  std::string line;
  for (char c : std::string("HGIPEBTS.")) {
    const auto n = std::count(t.dsspSS.begin(), t.dsspSS.end(), c);
    if (n > 0) line += std::format("  {} {}", c, n);
  }
  od::logInfo("    DSSP" + line);
}

// SHA-256 of the input files (spec 16.3).
std::vector<od::InputHash> inputHashes(const od::Options& opt) {
  std::vector<od::InputHash> v;
  auto add = [&](const std::string& role, const std::string& path) {
    v.push_back({role, std::filesystem::absolute(path).lexically_normal().string(), od::sha256File(path)});
  };
  add("stamFile", opt.stamFile);
  for (const auto& f : opt.mobmFiles) add("mobmFiles", f);
  add("molConfig", opt.molConfig);
  for (const auto& f : opt.ccd) add("ccd", f);
  return v;
}

// --resume (spec 16.4): parameters of the last checkpoint, then the command line
// (input files may be given again, e.g. when moved); only dock9c.maxRound may differ.
od::Options resumeOptions(int argc, char** argv, const od::Options& cli, Checkpointing& ck) {
  if (!cli.paramFile.empty()) od::fail("--paramFile cannot be used with --resume (the checkpoint parameters are used)");
  if (!cli.outPrefix.empty()) od::fail("--outPrefix cannot be used with --resume (the run continues in the --resume directory)");
  const std::string dir = cli.resume;
  const int round = od::latestCheckpoint(dir);
  if (round == 0) od::fail("--resume: no complete checkpoint in {}/checkpoint", dir);
  const std::string cdir = od::checkpointDir(dir, round);
  const od::Manifest man = od::readManifest(cdir + "/manifest.toml");
  if (man.version != od::kVersion)
    od::fail("{}: written by OmegaDocker {}, this is {}", cdir, man.version, od::kVersion);
  std::vector<std::string> args{argv[0], "--paramFile=" + cdir + "/params.toml", "--outPrefix=" + dir};
  for (int i = 1; i < argc; ++i)
    if (!std::string_view(argv[i]).starts_with("--resume=")) args.push_back(argv[i]);
  for (const auto& h : man.inputs) {
    if (h.role == "stamFile" && cli.stamFile.empty()) args.push_back("--stamFile=" + h.path);
    if (h.role == "mobmFiles" && cli.mobmFiles.empty()) args.push_back("--mobmFiles=" + h.path);
    if (h.role == "molConfig" && cli.molConfig.empty()) args.push_back("--molConfig=" + h.path);
  }
  auto parse = [&](const std::vector<std::string>& a) {
    std::vector<const char*> p;
    for (const auto& s : a) p.push_back(s.c_str());
    return od::parseCommandLine(static_cast<int>(p.size()), p.data());
  };
  od::Options opt = parse(args);
  opt.resume = dir;
  const od::Options stored = parse({argv[0], "--paramFile=" + cdir + "/params.toml"});
  od::Options cmp = opt;
  cmp.maxRound = stored.maxRound;
  const std::string a = od::paramsText(stored), b = od::paramsText(cmp);
  if (a != b) {
    auto lines = [](const std::string& s) {
      std::vector<std::string> v;
      std::string table;
      std::size_t p = 0;
      for (std::size_t q; (q = s.find('\n', p)) != std::string::npos; p = q + 1) {
        if (q == p || s[p] == '#') continue;
        if (s[p] == '[') table = s.substr(p + 1, q - p - 2) + ".";
        else v.push_back(table + s.substr(p, q - p));
      }
      return v;
    };
    const auto la = lines(a), lb = lines(b);
    std::string diff;
    for (std::size_t k = 0; k < std::max(la.size(), lb.size()); ++k) {
      const std::string x = k < la.size() ? la[k] : "", y = k < lb.size() ? lb[k] : "";
      if (x != y) diff += std::format("\n  checkpoint: {}\n  now:        {}", x, y);
    }
    od::fail("--resume: only dock9c.maxRound may differ from the checkpoint parameters:{}", diff);
  }
  ck.resumeDir = dir;
  ck.resumeRound = round;
  ck.stateSha256 = man.stateSha256;
  ck.inputs = man.inputs;
  return opt;
}

int run(int argc, char** argv) {
  od::Options opt = od::parseCommandLine(argc, argv);
  if (opt.showHelp) {
    std::fputs(od::usageText().c_str(), stdout);
    return 0;
  }
  if (opt.showVersion) {
    od::logInfo(std::format("OmegaDocker {}", od::kVersion));
    return 0;
  }
  Checkpointing ck;
  if (!opt.resume.empty()) opt = resumeOptions(argc, argv, opt, ck);
  od::validateOptions(opt);

  od::CcdLibrary ccd;
  for (const auto& f : opt.ccd) ccd.addFile(f);

  auto sta = od::loadMolecules(opt.stamFile, od::Role::StaM, ccd);
  od::logInfo(std::format("StaM group from {}: {} molecule(s)", opt.stamFile, sta.size()));
  std::vector<od::Topology> staTopo;
  for (const auto& m : sta) {
    printMolecule(m);
    staTopo.push_back(analyze(m, opt));
  }
  {
    std::vector<const od::Molecule*> mols;
    std::vector<od::Topology*> topos;
    for (std::size_t i = 0; i < sta.size(); ++i) {
      mols.push_back(&sta[i]);
      topos.push_back(&staTopo[i]);
    }
    od::assignDssp(mols, topos);
    for (const auto& t : staTopo) printDssp(t);
  }
  std::vector<od::Molecule> mobs;
  std::vector<od::Topology> mobTopo;
  for (const auto& f : opt.mobmFiles) {
    auto mob = od::loadMolecules(f, od::Role::MobM, ccd);
    if (mob.size() != 1) od::fail("{}: --mobmFiles file must contain exactly one molecule, found {}", f, mob.size());
    od::logInfo(std::format("MobM from {}:", f));
    printMolecule(mob.front());
    od::Topology t = analyze(mob.front(), opt);
    od::assignDssp({&mob.front()}, {&t});
    printDssp(t);
    mobs.push_back(std::move(mob.front()));
    mobTopo.push_back(std::move(t));
  }

  // Molecule configuration (spec 5.1).
  std::vector<od::MolInput> inputs;
  for (std::size_t i = 0; i < sta.size(); ++i) inputs.push_back({&sta[i], &staTopo[i]});
  for (std::size_t i = 0; i < mobs.size(); ++i) inputs.push_back({&mobs[i], &mobTopo[i]});
  if (opt.molConfig.empty()) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs);
    od::logInfo(std::format("Default molecule configuration written to {}\n"
                            "Review it (a copy may be edited) and run again with --molConfig <file>.",
                            path));
    return 0;
  }
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(opt.molConfig, inputs, topologyParams(opt), errors);
  if (!cfg) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs);
    std::string list;
    for (const auto& e : errors) list += "\n  " + e;
    od::fail("{}: {} problem(s) found:{}\nDefault molecule configuration written to {}", opt.molConfig, errors.size(),
             list, path);
  }
  {
    const auto now = inputHashes(opt);
    if (!ck.resumeDir.empty()) {
      if (now.size() != ck.inputs.size()) od::fail("--resume: the number of input files differs from the checkpoint");
      for (std::size_t k = 0; k < now.size(); ++k)
        if (now[k].role != ck.inputs[k].role || now[k].sha256 != ck.inputs[k].sha256)
          od::fail("--resume: {} {} differs from the checkpoint ({}; SHA-256 {} vs {})", now[k].role, now[k].path,
                   ck.inputs[k].path, now[k].sha256, ck.inputs[k].sha256);
    }
    ck.inputs = now;
  }
  od::logInfo(std::format("Molecule configuration {}:", opt.molConfig));
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    const od::MolSettings& s = cfg->molecules[i];
    checkOmegas(*inputs[i].mol, s, opt);
    const auto free = std::count(s.mode.begin(), s.mode.end(), od::TorsionMode::Free);
    const auto cisTrans = std::count(s.mode.begin(), s.mode.end(), od::TorsionMode::CisTrans);
    od::logInfo(std::format("  {:<40} rotatable {:>5}  cisTrans omega {:>4}  ringBreak {:>3}  chainBreak {:>3}",
                            inputs[i].mol->id, free, cisTrans, s.ringBreaks.size(), s.chainBreaks.size()));
    printTree(*inputs[i].mol, s);
  }

  {
    std::size_t staTorsions = 0, mobTorsions = 0;
    for (std::size_t i = 0; i < inputs.size(); ++i)
      (i < sta.size() ? staTorsions : mobTorsions) += cfg->molecules[i].tree.torsions.size();
    if (staTorsions > 0) od::requireConfSampleNums(opt, od::SampleSet::StaM);
    if (mobTorsions > 0) od::requireConfSampleNums(opt, od::SampleSet::MobM);
  }

  od::logInfo("9a energy:");
  {
    std::vector<od::EnergyPart> parts;
    for (std::size_t i = 0; i < sta.size(); ++i) parts.push_back({&sta[i], &cfg->molecules[i].tree});
    printInputEnergy("StaM group", parts, opt);
  }
  for (std::size_t i = 0; i < mobs.size(); ++i)
    printInputEnergy(mobs[i].id, {{&mobs[i], &cfg->molecules[sta.size() + i].tree}}, opt);

  // Random seeds fit in a TOML integer (63 bits).
  const std::uint64_t seed =
      opt.seed ? *opt.seed
               : ((static_cast<std::uint64_t>(std::random_device{}()) << 32) ^ std::random_device{}()) >> 1;
  od::logInfo(std::format("seed {}", seed));
  opt.seed = seed;
  {
    std::filesystem::create_directories(opt.outPrefix);
    const std::string path = opt.outPrefix + "/params.used.toml";
    std::ofstream f(path);
    if (!f) od::fail("cannot write {}", path);
    ck.params = od::paramsText(opt);
    f << ck.params;
    od::logInfo(std::format("parameters used (defaults < {} < command line) written to {}",
                            opt.paramFile.empty() ? "no --paramFile" : opt.paramFile, path));
  }
  std::unique_ptr<ConfSet> stamSet;
  {
    std::vector<od::EnergyPart> parts;
    std::vector<const od::MolSettings*> settings;
    for (std::size_t i = 0; i < sta.size(); ++i) {
      parts.push_back({&sta[i], &cfg->molecules[i].tree});
      settings.push_back(&cfg->molecules[i]);
    }
    stamSet = std::make_unique<ConfSet>("StaM group", parts, settings, *cfg, opt, od::SampleSet::StaM, seed, 0);
  }
  std::vector<std::unique_ptr<ConfSet>> mobSets;
  for (std::size_t i = 0; i < mobs.size(); ++i) {
    const auto& s = cfg->molecules[sta.size() + i];
    mobSets.push_back(std::make_unique<ConfSet>(mobs[i].id, std::vector<od::EnergyPart>{{&mobs[i], &s.tree}},
                                                std::vector<const od::MolSettings*>{&s}, *cfg, opt,
                                                od::SampleSet::MobM, seed, 1 + i));
  }
  for (const ConfSet* c : {stamSet.get()})
    if (!c->movable()) od::logInfo(std::format("  {}: no rotatable torsions; the input conformation is the only sample", c->name()));
  for (const auto& c : mobSets)
    if (!c->movable()) od::logInfo(std::format("  {}: no rotatable torsions; the input conformation is the only sample", c->name()));
  runRounds(*stamSet, mobSets, opt, opt.fe, seed, ck);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const od::Error& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  } catch (const std::exception& e) {
    std::fflush(stdout);
    std::fprintf(stderr, "ERROR: unexpected: %s\n", e.what());
    return 1;
  }
}
