// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include <algorithm>
#include <chrono>
#include <cmath>
#include <array>
#include <cstdio>
#include <deque>
#include <condition_variable>
#include <mutex>
#include <exception>
#include <atomic>
#include <thread>
#include <filesystem>
#include <fstream>
#include <functional>
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

#include "box.hpp"
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
#include "pool.hpp"
#include "reprotonate.hpp"
#include "ringconf.hpp"
#include "symrmsd.hpp"
#include "feparams.hpp"
#include "score9d.hpp"
#include "shaking.hpp"
#include "sysmem.hpp"
#include "write_struct.hpp"
#include "molconfig.hpp"
#include "topology.hpp"
#include "version.hpp"
#ifdef OD_HAVE_CUDA
#include "conf9a_gpu.hpp"
#include "rcd_gpu.hpp"
#include "dock9b_gpu.hpp"
#include "grid_gpu.hpp"
#include "shake_gpu.hpp"
#endif

namespace {

// Protonation of one input file (spec 3.4): log summary.
void logProtonation(const std::string& file, const od::ProtonationReport& r, double pH) {
  std::map<std::string, int> states;
  for (const auto& s : r.sites) ++states[s.state];
  std::string st;
  for (const auto& [k, n] : states) st += std::format("{}{} {}", st.empty() ? "" : ", ", k, n);
  od::logInfo(std::format("Protonation at pH {} ({}): {} input hydrogen(s) removed, {} added; {} titratable site(s){}{}", pH, file,
                          r.removedH, r.addedH, r.sites.size(), st.empty() ? "" : ": ", st));
  od::logInfo(std::format("  hydrogen optimization: {} rotor(s), {} flip group(s), {} flipped{}; {} component(s){}", r.rotors,
                          r.flipGroups, r.flipped.size(),
                          r.flipped.empty() ? "" : " (" + [&] {
                            std::string v;
                            for (const auto& x : r.flipped) v += (v.empty() ? "" : " ") + x;
                            return v;
                          }() + ")",
                          r.components,
                          r.greedyComponents > 0 ? std::format(", {} solved greedily (too many states)", r.greedyComponents) : ""));
  for (const auto& s : r.sites)
    if (s.method == "override" || s.method == "metal" ||
        (s.state != "ASP" && s.state != "GLU" && s.state != "COO-" && s.state != "NH3+" && s.state != "CYS" && s.state != "TYR" &&
         s.state != "LYS" && s.state != "ARG" && s.state != "HID/HIE" && s.method.starts_with("PROPKA")))
      od::logInfo(std::format("  {:<24} {:<20} {}{}", s.atom, s.state, s.method, std::isnan(s.pka) ? "" : std::format(" pKa {:.2f}", s.pka)));
  if (!r.addedOxt.empty()) {
    std::string v;
    for (const auto& x : r.addedOxt) v += (v.empty() ? "" : " ") + x;
    od::logInfo(std::format("  OXT added (C terminus without a link): {}", v));
  }
  if (!r.incomplete.empty()) {
    std::string v;
    for (const auto& x : r.incomplete) v += (v.empty() ? "" : " ") + x;
    od::logWarn(std::format("{}: residue(s) with missing heavy atoms; their open valences are filled with hydrogens: {}", file, v));
  }
}

// [[protonation]] section of the default molecule configuration (spec 3.4, 5).
std::string protonationText(const od::Options& opt, const std::vector<std::pair<std::string, od::ProtonationReport>>& reports,
                            const std::vector<od::ProtonationOverride>& overrides) {
  if (reports.empty()) return "";
  std::string o = "\n# ---------- protonation (stam.protonation / mobm.protonation = \"pH\", spec 3.4) ----------\n";
  o += std::format("# protonation.pH = {}. Predicted states: atom, state, method (pKa).\n", opt.protonationPH);
  o += "# Change a state with an entry (the default configuration may then need to be written again):\n";
  o += "#   [[protonation]]  file = \"<input file>\", site = \"<chain>:<seq><iCode>\", state = \"<state>\"\n";
  o += "#     ASP/ASH, GLU/GLH, HID/HIE/HIP, CYS/CYM, TYR/TYM, LYS/LYN, ARG/ARN\n";
  o += "#   [[protonation]]  file = \"<input file>\", site = \"<atom>\", charge = <formal charge>\n";
  o += "#     <atom>: <chain>:<seq><iCode>:<name> or #<atom serial in the input record>; hydrogens follow\n";
  o += "#     from the valence (e.g. termini: N charge 0 or 1, OXT charge -1 or 0)\n";
  for (const auto& [file, r] : reports) {
    o += std::format("# {}\n", file);
    for (const auto& s : r.sites)
      o += std::format("#   {:<24} {:<20} {}{}\n", s.atom, s.state, s.method, std::isnan(s.pka) ? "" : std::format(" ({:.2f})", s.pka));
  }
  for (const auto& v : overrides) {
    o += std::format("\n[[protonation]]\nfile = \"{}\"\nsite = \"{}\"\n", v.file, v.site);
    if (v.charge) o += std::format("charge = {}\n", *v.charge);
    else o += std::format("state = \"{}\"\n", v.state);
  }
  return o;
}

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
  p.box = od::Box{opt.box.center, opt.box.size};
  p.stamAllFixed = opt.boxMode == "fixAuto";
  p.stamFlexMargin = opt.stamFlexMargin;
  p.stamFixedResMovable = opt.stamFixedResMovable;
  return p;
}

// Macrocycles (spec 4.9): broken by default ([[molecule.ringBreak]] source = "auto-macro").
void addMacroBreaks(const od::Molecule& m, od::Topology& t, od::TopologyParams p, int macroMinSize) {
  if (const auto macro = od::macroRingBreaks(m, t, macroMinSize); !macro.empty()) {
    p.brokenBonds = t.brokenBonds;
    p.brokenBonds->insert(p.brokenBonds->end(), macro.begin(), macro.end());
    const std::string dssp = t.dsspSS;
    t = od::analyzeTopology(m, p);
    t.dsspSS = dssp;
    t.macroBreaks = macro;
  }
}

od::Topology analyze(const od::Molecule& m, const od::Options& opt) {
  od::Topology t = od::analyzeTopology(m, topologyParams(opt));
  addMacroBreaks(m, t, topologyParams(opt), opt.ringMacroMinSize);
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

// Closure settings of one group (spec 7.11).
od::ClosureParams closureParams(const od::CsgmOpts& c, const od::Options& opt) {
  od::ClosureParams clp;
  clp.method = c.closure == "kic" ? od::ClosureMethod::Kic : c.closure == "rcd" ? od::ClosureMethod::Rcd : od::ClosureMethod::None;
  clp.rcdTolerance = c.rcdTolerance;
  clp.rcdMaxIter = c.rcdMaxIter;
  clp.rcdStallTolerance = c.rcdStallTolerance;
  clp.rcdMaxRestart = c.rcdMaxRestart;
  clp.kicMaxTrial = c.kicMaxTrial;
  clp.rcdClash = c.rcdClash == "full" ? od::RcdClash::Full : c.rcdClash == "repair" ? od::RcdClash::Repair : od::RcdClash::None;
  clp.rcdClashWeight = c.rcdClashWeight;
  clp.clashDistance = opt.clashDistance;
  clp.growClose = c.build == "growClose";
  return clp;
}

od::EnergyParams energyParams(const od::Options& opt) {
  od::EnergyParams p;
  p.clashDistance = opt.clashDistance;
  p.clashBase = opt.EclashBase;
  p.clashK = opt.Kclash;
  p.closureG = opt.closureG;
  return p;
}

// 9a / 9b energies (energy.searchWeights): od_0100.inter.* for the 9b grids, od_0100.intra.* otherwise.
od::EnergyParams searchParams(const od::Options& opt, bool inter) {
  const od::EnergyParams p = energyParams(opt);
  if (opt.searchWeights != "od_0100") return p;
  return od::weighted(p, inter ? opt.fe.inter : opt.fe.intra);
}

// 9a energy of the input conformation (spec 7.9).
void logRingConformers(const od::Molecule& m, const od::RingConformers& rc) {
  if (rc.systems.empty()) return;
  int sampled = 0, conformers = 0;
  std::map<std::string, int> why;
  for (const auto& s : rc.systems) {
    if (s.conf.size() > 1) {
      ++sampled;
      conformers += static_cast<int>(s.conf.size());
    } else {
      ++why[s.note];
    }
  }
  std::string skipped;
  for (const auto& [w, n] : why) skipped += std::format("{}{} {}", skipped.empty() ? "" : ", ", n, w);
  od::logInfo(std::format("  {:<40} small ring systems {}: {} with conformers ({} conformers){}", m.id, rc.systems.size(),
                          sampled, conformers, skipped.empty() ? "" : "; input only: " + skipped));
  if (m.role != od::Role::MobM) return;
  for (std::size_t k = 0; k < rc.systems.size(); ++k) {
    const auto& s = rc.systems[k];
    std::string e;
    for (double v : s.energy) e += std::format(" {:.2f}", v);
    od::logInfo(std::format("    system {}: {} ring atoms (largest ring {}), {} conformer(s) of {} embedded{}{}", k + 1,
                            s.ring.size(), s.largestRing, s.conf.size(), s.generated,
                            e.empty() ? "" : ", relative MMFF energies (input first):" + e, s.note.empty() ? "" : " (" + s.note + ")"));
  }
}

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

// Whether 9a and 9b minimization run on the GPU (spec 15.1.3, 15.1.4); logs the choice.
bool selectGpu(const od::Options& opt) {
#ifdef OD_HAVE_CUDA
  if (!opt.gpuEnable) {
    od::logInfo("9a/9b minimization energies: CPU (gpu.enable = false)");
    return false;
  }
  const int n = od::gpuDeviceCount();
  if (n == 0) {
    od::logInfo("9a/9b minimization energies: CPU (no CUDA device)");
    return false;
  }
  for (int d : opt.gpuDevices)
    if (d >= n) od::fail("gpu.device {}: only {} CUDA device(s)", d, n);
  if (opt.gpuDock9bDeviceMin && opt.dockMinMethod == "lbfgs" && opt.lbfgsMemory > od::kDeviceLbfgsMaxMemory)
    od::fail("dock9b.lbfgsMemory {}: at most {} with gpu.dock9bDeviceMin = true (set gpu.dock9bDeviceMin = false to "
             "minimize on the CPU)",
             opt.lbfgsMemory, od::kDeviceLbfgsMaxMemory);
  for (std::size_t k = 0; k < opt.gpuDevices.size(); ++k) {
    const int d = opt.gpuDevices[k];
    const auto info = od::gpuDeviceInfo(d);
    od::logInfo(std::format("{} GPU {} {} (sm_{}{}, {} MiB free of {} MiB)",
                            k == 0 ? "9a/9b minimization energies:" : "  9a/9b minimization also on:", d, info.name, info.major,
                            info.minor, info.freeMem >> 20, info.totalMem >> 20));
  }
  return true;
#else
  (void)opt;
  od::logInfo("9a/9b minimization energies: CPU (built without CUDA)");
  return false;
#endif
}

// Conformation samples of the StaM group or one MobM: the L2 pool, its 9a
// refills and the samples moved to a dcss (spec 7.3). Sample ids index samples().
class ConfSet {
 public:
  // prot: the [[protonation]] entries of each part for protonation.recompute (spec 3.4.1); nullopt
  // when the hydrogens are those of the input file.
  ConfSet(std::string name, std::vector<od::EnergyPart> parts, std::vector<const od::MolSettings*> settings,
          const od::MolConfig& cfg, const od::Options& opt, od::SampleSet set, std::uint64_t seed, std::uint64_t stream,
          bool gpu, const std::optional<std::vector<std::vector<od::ProtonationOverride>>>& prot = std::nullopt)
      : name_(std::move(name)), parts_(std::move(parts)), settings_(std::move(settings)), opt_(opt),
        csgm_(set == od::SampleSet::StaM ? opt.stamCsgm : opt.mobmCsgm), sys_(parts_, searchParams(opt, false)) {
    for (const auto& n : sys_.ringNotes()) od::logInfo(std::format("  {:<40} ring conformers not used: {}", name_, n));
    if (sys_.ringSlotCount() > 0) {
      std::string k;
      for (int c : sys_.ringConfCounts()) k += std::format("{}{}", k.empty() ? "" : ", ", c);
      od::logInfo(std::format("  {:<40} {} ring conformer variable(s), conformers: {}", name_, sys_.ringSlotCount(), k));
    }
    if (opt.protonationRecompute != "none" && prot) {
      prot_ = std::make_unique<od::ProtonationVariants>(sys_, parts_, *prot, searchParams(opt, false));
      int n = 0, locked = 0;
      for (std::size_t k = 0; k < prot_->partCount(); ++k)
        for (const auto& s : prot_->sites(static_cast<int>(k))) {
          ++n;
          locked += s.locked;
        }
      od::logInfo(std::format("  {:<40} protonation.recompute {}: {} titratable site(s) of standard amino acids, {} fixed "
                              "(metal-bonded or [[protonation]])",
                              name_, opt.protonationRecompute, n, locked));
      if (!prot_->active()) prot_.reset();
    }
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
    p_.preCloseFactor = opt.confPreCloseFactor;
    p_.min = {opt.confMinMaxStep,
              opt.confMinDeltaE,
              opt.confMaxDihedralStep,
              od::parseMinMethod(opt.confMinMethod),
              od::parseLineSearchInit(opt.confLineSearchInit),
              od::parseLineSearchReduce(opt.confLineSearchReduce),
              opt.confMaxBacktrack};
    p_.seed = seed;
    p_.stream = stream;
#ifdef OD_HAVE_CUDA
    if (gpu) {
      // Several devices: the lanes of all of them (spec 15.1.33).
      std::vector<std::unique_ptr<od::GpuConfEvaluator>> ev;
      for (int dev : opt.gpuDevices)
        ev.push_back(std::make_unique<od::GpuConfEvaluator>(sys_.model(), dev, opt.gpuConf9aBatchNum,
                                                            opt.gpuConf9aPrecision == "float"));
      if (ev.size() == 1) gpu_ = std::move(ev[0]);
      else gpu_ = std::make_unique<od::MultiGpuConfEvaluator>(std::move(ev), sys_.torsionCount());
      p_.evaluator = gpu_.get();
      if (gpu_->maxBatch() < opt.gpuConf9aBatchNum * static_cast<int>(opt.gpuDevices.size()))
        od::logInfo(std::format("  {:<40} 9a GPU batch limited to {} samples by the free device memory", name_, gpu_->maxBatch()));
    }
#else
    (void)gpu;
#endif
    const od::ClosureParams clp = closureParams(csgm_, opt);
    std::vector<bool> freeMask(cisTrans_.size());
    for (std::size_t g = 0; g < freeMask.size(); ++g) freeMask[g] = !cisTrans_[g];
    closer_.emplace(sys_, parts_, freeMask, sampler_.segments(), clp);
    if (clp.method != od::ClosureMethod::None)
      for (const auto& line : closer_->describe()) od::logInfo(std::format("  {:<40} closure: {}", name_, line));
    const bool closing = closer_->active();
#ifdef OD_HAVE_CUDA
    if (gpu && opt.gpuConf9aRcd && closing && clp.method != od::ClosureMethod::None) {
      auto r = std::make_unique<od::GpuRcd>(*closer_, opt.gpuDevices, opt.gpuConf9aBatchNum);
      if (r->active()) {
        rcdGpu_ = std::move(r);
        p_.genBatch = rcdGpu_.get();
        p_.genBatchMin = opt.gpuConf9aRcdMinBatch;
        od::logInfo(std::format("  {:<40} closure: RCD attempts on the GPU for L1 batches of at least {} samples", name_,
                                opt.gpuConf9aRcdMinBatch));
      }
    }
#endif
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
  // Protonation variants (spec 3.4.1); null without protonation.recompute or titratable sites.
  od::ProtonationVariants* prot() const { return prot_.get(); }
  const od::EnergySystem& sysOf(int variant) const { return variant == 0 || !prot_ ? sys_ : prot_->sys(variant); }
  std::vector<od::EnergyPart> partsOf(int variant) const { return variant == 0 || !prot_ ? parts_ : prot_->parts(variant); }
  const std::vector<const od::MolSettings*>& settings() const { return settings_; }
  const od::ConfSample& sample(int id) const { return samples_[id]; }
  // Lowest 9a energy of the conformation samples made so far (L2 pool and dcss samples).
  float minEnergy() const {
    float e = std::numeric_limits<float>::infinity();
    if (!poolOnly_)
      for (const auto& s : samples_) e = std::min(e, s.energy);
    for (const auto& s : pool_) e = std::min(e, s.energy);
    return e;
  }
  // The n samples with the lowest 9a energy (L2 pool and dcss samples); the input conformation when not movable.
  std::vector<od::ConfSample> lowestSamples(int n) const {
    if (!movable()) return {od::ConfSample{}};
    std::vector<const od::ConfSample*> v;
    if (!poolOnly_)
      for (const auto& s : samples_) v.push_back(&s);
    for (const auto& s : pool_) v.push_back(&s);
    std::stable_sort(v.begin(), v.end(), [](auto a, auto b) { return a->energy < b->energy; });
    std::vector<od::ConfSample> out;
    for (std::size_t k = 0; k < v.size() && k < static_cast<std::size_t>(n); ++k) out.push_back(*v[k]);
    return out;
  }
  // Checkpoint state (spec 16.3): random stream counter, L1 batch state, L2 pool and dcss samples.
  void save(od::BinWriter& w) const {
    for (int v : {n_.docking, n_.minL2, n_.initL1, n_.elite}) w.put<std::int32_t>(v);
    w.put(counter_);
    w.put<std::uint8_t>(first_);
    w.put(nextL1_);
    const bool l2 = recomputeL2();
    if (l2) {
      // Protonation variants in creation order: per part, (site, state) of the changes.
      const auto n = static_cast<std::uint32_t>(prot_->count());
      w.put(n);
      for (std::uint32_t id = 1; id < n; ++id)
        for (const auto& part : prot_->changes(static_cast<int>(id))) {
          w.put<std::uint32_t>(static_cast<std::uint32_t>(part.size()));
          for (const auto& c : part) {
            w.put<std::int32_t>(c.site);
            w.putVec(std::vector<char>(c.state.begin(), c.state.end()));
          }
        }
    }
    for (const auto* v : {&pool_, &samples_}) {
      w.put<std::uint64_t>(v->size());
      for (const auto& s : *v) {
        w.putVec(s.angles);
        w.put(s.energy);
        if (l2) {
          w.put<std::int32_t>(s.variant);
          w.put(s.protG);
        }
      }
    }
  }
  void load(od::BinReader& r) {
    od::ConfSampleNums n;
    for (int* v : {&n.docking, &n.minL2, &n.initL1, &n.elite}) *v = r.get<std::int32_t>();
    if (movable()) {
      if (n.docking < n_.docking || n.minL2 < 2 * n.docking || n.elite < 2)
        od::fail("checkpoint: {}: invalid sample numbers (dcss size {})", name_, n.docking);
      setSampleNums(n);
    }
    counter_ = r.get<std::uint64_t>();
    first_ = r.get<std::uint8_t>() != 0;
    nextL1_ = r.get<long long>();
    const bool l2 = recomputeL2();
    if (l2) {
      const auto n = r.get<std::uint32_t>();
      for (std::uint32_t id = 1; id < n; ++id) {
        std::vector<std::vector<od::StateChange>> ch(prot_->partCount());
        for (auto& part : ch) {
          const auto m = r.get<std::uint32_t>();
          for (std::uint32_t k = 0; k < m; ++k) {
            od::StateChange c;
            c.site = r.get<std::int32_t>();
            const auto st = r.getVec<char>();
            c.state.assign(st.begin(), st.end());
            if (c.site < 0 || c.site >= static_cast<int>(prot_->sites(static_cast<int>(&part - ch.data())).size()))
              od::fail("checkpoint: {}: protonation variant with an unknown site", name_);
            part.push_back(std::move(c));
          }
        }
        if (prot_->variantOf(ch) != static_cast<int>(id)) od::fail("checkpoint: {}: protonation variants differ", name_);
      }
    }
    for (auto* v : {&pool_, &samples_}) {
      const auto n = r.get<std::uint64_t>();
      v->clear();
      for (std::uint64_t k = 0; k < n; ++k) {
        od::ConfSample s;
        s.angles = r.getVec<float>();
        s.energy = r.get<float>();
        if (l2) {
          s.variant = r.get<std::int32_t>();
          s.protG = r.get<float>();
          if (s.variant < 0 || s.variant >= static_cast<int>(prot_->count()))
            od::fail("checkpoint: {}: unknown protonation variant {}", name_, s.variant);
        }
        if (s.angles.size() != static_cast<std::size_t>(sys_.torsionCount()))
          od::fail("checkpoint: {}: {} torsions in a sample, {} expected", name_, s.angles.size(), sys_.torsionCount());
        v->push_back(std::move(s));
      }
    }
    if (samples_.empty() && !movable()) od::fail("checkpoint: {}: input conformation sample missing", name_);
  }
  std::size_t sampleCount() const { return samples_.size(); }

  // score / localMin (spec 6.1): the input conformation as a sample; returns its id. fill9a: the L2
  // pool filled by csgm1, and minEnergy / lowestAngles from the pool only (run.unbound = conf9a).
  int addInput(bool fill9a) {
    if (!movable()) return 0;
    od::ConfSample s;
    const auto in = sys_.inputAngles();
    s.angles.assign(in.begin(), in.end());
    s.energy = static_cast<float>(sys_.evaluate(std::span<const float>(s.angles), nullptr, nullptr, nullptr).filterEnergy());
    samples_.push_back(std::move(s));
    if (fill9a) {
      refill({});
      poolOnly_ = true;
    }
    return static_cast<int>(samples_.size()) - 1;
  }
  int dockingNum() const { return movable() ? n_.docking : 1; }
  int eliteNum() const { return movable() ? n_.elite : 1; }
  const od::ConfSampleNums& sampleNums() const { return n_; }
  // dcss size change (spec 9, step 6): minL2confSampleNum and the others recomputed.
  void setSampleNums(const od::ConfSampleNums& n) {
    n_ = n;
    p_.minL2 = n.minL2;
  }
  std::size_t poolSize() const { return pool_.size(); }

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
    od::RcdRunner* runner = nullptr;
#ifdef OD_HAVE_CUDA
    runner = rcdGpu_.get();
    const long long k0 = rcdGpu_ ? rcdGpu_->kernels() : 0, a0 = rcdGpu_ ? rcdGpu_->attempts() : 0;
    const double s0 = rcdGpu_ ? rcdGpu_->kernelSec() : 0;
#endif
    const od::ConfGenerator gen = [&](od::Rng& r, const od::GenIndex& idx, std::span<double> a) {
      const od::RcdSlot rs{idx.slot >= 0 ? runner : nullptr, idx.slot};
      if (pe) {
        const int m = static_cast<int>(idx.local % pe_.perLine);
        od::csgm3Sample(od::Rng(p_.seed, p_.stream, idx.global - m), ends, idx.local / pe_.perLine, m, cisTrans_, pe_, a);
      } else if (ga) {
        od::csgm2Offspring(r, parents, cisTrans_, ga_, a, sys_.ringConfCounts());
      } else {
        sampler_.sample(r, a);
        // Ring conformers (spec 4.9): uniform, from their own random stream.
        if (sys_.ringSlotCount() > 0) {
          od::Rng rr(p_.seed, (1ULL << 61) + p_.stream, idx.global);
          const int t0 = sys_.treeTorsionCount();
          for (int s = 0; s < sys_.ringSlotCount(); ++s)
            a[t0 + s] = static_cast<double>(rr.below(static_cast<std::uint64_t>(sys_.ringConfCounts()[s])));
        }
      }
      // csgm2 / csgm3: the torsions are given; building does not resample.
      const od::SegmentSampler& growResample = pe || ga ? od::SegmentSampler{} : resample;
      if (growClose_) {
        // growClose: closure paths as sampled, closure, then the rest built segment by segment (spec 7.8).
        if (!closer_->close(r, a, resample, rs)) return od::GenResult::CloseFailed;
        if (!grow_->build(r, a, growResample).ok) return od::GenResult::SideGrowFailed;
        return od::GenResult::Built;
      }
      if (grow_ && !grow_->build(r, a, growResample).ok) return od::GenResult::GrowFailed;
      if (closing && !closer_->close(r, a, resample, rs)) return od::GenResult::CloseFailed;
      return od::GenResult::Built;
    };
    const long long l1 = first_ ? n_.initL1 : nextL1_;
    first_ = false;
    p_.fixedL1 = pe ? od::csgm3SampleNum(ends.size(), pe_) : 0;
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t before = pool_.size();
    const auto st = od::fillL2Pool(sys_, gen, l1, p_, counter_, pool_);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (recomputeL2()) recomputePool(before);
    nextL1_ = st.nextL1;
    std::string rcdInfo;
#ifdef OD_HAVE_CUDA
    if (rcdGpu_)
      rcdInfo = std::format(", GPU RCD {} attempts in {} kernels {:.2f} s", rcdGpu_->attempts() - a0, rcdGpu_->kernels() - k0,
                            rcdGpu_->kernelSec() - s0);
#endif
    const long long minimized = st.l1 - st.growFailed - st.closeFailed - st.sideGrowFailed - st.preNotClosed;
    od::logInfo(std::format("  9a {:<37} {} L1 {:>9} in {} batch(es)  -> L2 +{} = {}  (rejected: grow {}, closure {}, "
                            "grow after closure {}, not closed before minimization {}, E >= Emax {}, not closed {})  "
                            "steps/sample {:.1f}  {:.2f} s{}",
                            name_, pe ? "csgm3" : ga ? "csgm2" : "csgm1", st.l1, st.batches, st.added, pool_.size(),
                            st.growFailed, st.closeFailed, st.sideGrowFailed, st.preNotClosed, st.overEmax, st.notClosed,
                            minimized > 0 ? static_cast<double>(st.minSteps) / static_cast<double>(minimized) : 0.0, sec,
                            p_.evaluator ? std::format(" (generate {:.2f} s{}, GPU minimize {:.2f} s)", st.genSec, rcdInfo, st.minSec) : ""));
  }

  bool recomputeL2() const { return prot_ && opt_.protonationRecompute == "l2" && movable(); }
  // protonation.recompute = l2 (spec 3.4.1): states of the samples pool_[begin..] from PROPKA at their
  // conformation; a sample with changes takes the energy of its variant plus the protonation free energy.
  void recomputePool(std::size_t begin) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t n = pool_.size() - begin;
    std::vector<char> drop(n, 0);
    std::vector<std::vector<std::vector<od::StateChange>>> ch(n);
    std::vector<double> g(n, 0.0);
    std::vector<int> tautomerFlips(n, 0);
    od::parallelEach(static_cast<long long>(n), [&](long long k) {
      const od::ConfSample& s = pool_[begin + static_cast<std::size_t>(k)];
      const auto x = sys_.coordinates(std::span<const float>(s.angles));
      const auto pka = od::calculatePka(prot_->pkaMolecules(x));
      ch[k].resize(prot_->partCount());
      for (std::size_t p = 0; p < ch[k].size(); ++p)
        ch[k][p] = od::stateChanges(prot_->sites(static_cast<int>(p)), pka[p], opt_.protonationPH);
      // Neutral His: the tautomer of the lowest energy of the sample.
      tautomerFlips[k] = prot_->chooseTautomers(ch[k], [&](const od::EnergySystem& sys, const std::vector<od::EnergyPart>&) {
        return sys.evaluate(std::span<const float>(s.angles), nullptr, nullptr, nullptr).filterEnergy();
      });
      for (std::size_t p = 0; p < ch[k].size(); ++p)
        g[k] += od::protonationFreeEnergy(prot_->sites(static_cast<int>(p)), ch[k][p], opt_.protonationPH, opt_.fe.temperature);
    });
    // Variant ids in sample order (reproducible).
    for (std::size_t k = 0; k < n; ++k) pool_[begin + k].variant = prot_->variantOf(ch[k]);
    od::parallelEach(static_cast<long long>(n), [&](long long k) {
      od::ConfSample& s = pool_[begin + static_cast<std::size_t>(k)];
      if (s.variant == 0) return;
      const double e = prot_->sys(s.variant).evaluate(std::span<const float>(s.angles), nullptr, nullptr, nullptr).filterEnergy();
      s.protG = static_cast<float>(g[k]);
      s.energy = static_cast<float>(e + g[k]);
      if (!(e < opt_.Emax)) drop[k] = 1;
    });
    std::size_t changed = 0, dropped = 0;
    std::vector<od::ConfSample> kept(pool_.begin(), pool_.begin() + static_cast<long>(begin));
    for (std::size_t k = 0; k < n; ++k) {
      od::ConfSample& s = pool_[begin + k];
      changed += s.variant != 0;
      if (drop[k]) ++dropped;
      else kept.push_back(std::move(s));
    }
    pool_ = std::move(kept);
    od::logInfo(std::format("  9a {:<37} protonation.recompute: {} of {} sample(s) changed (E >= Emax with the variant: {} "
                            "removed); {} variant(s), {} His tautomer change(s)  {:.2f} s",
                            name_, changed, n, dropped, prot_->count() - 1,
                            std::accumulate(tautomerFlips.begin(), tautomerFlips.end(), 0LL),
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
  }

  std::string name_;
  std::vector<od::EnergyPart> parts_;
  std::vector<const od::MolSettings*> settings_;
  const od::Options& opt_;
  const od::CsgmOpts& csgm_;
  od::EnergySystem sys_;
  std::unique_ptr<od::ProtonationVariants> prot_;
  od::Csgm1GroupSampler sampler_;
  std::vector<bool> cisTrans_;
  od::Csgm2Params ga_;
  od::Csgm3Params pe_;
  std::optional<od::GrowBuilder> grow_;
  bool growClose_ = false;
  std::optional<od::LoopCloser> closer_;
  od::ConfSampleNums n_;
  std::unique_ptr<od::ConfEvaluator> gpu_;
#ifdef OD_HAVE_CUDA
  std::unique_ptr<od::GpuRcd> rcdGpu_;
#endif
  od::Fill9aParams p_;
  std::uint64_t counter_ = 0;
  bool first_ = true;
  long long nextL1_ = 0;
  std::vector<od::ConfSample> pool_;     // L2 pool, first in first out
  std::vector<od::ConfSample> samples_;  // samples moved to a dcss; the input conformation when not movable
  bool poolOnly_ = false;                // addInput(true)
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
  dp.elecWeight = searchParams(opt, true).weights[2];
  dp.sepaChargeThreshold = opt.sepaChargeThreshold;
  dp.sepaDistTolerance = opt.sepaDistTolerance;
  dp.sepaPosStep = opt.sepaPosStep;
  dp.sepaPosRange = opt.sepaPosRange.value_or(opt.sepaPosStep);
  dp.sepaPair = opt.sepaSamplingPair;
  dp.sepaSingle = opt.sepaSamplingSingle;
  dp.sepaTriple = opt.sepaSamplingTriple;
  dp.rotStep = opt.rotStep;
  dp.boxCenter = opt.box.center;
  dp.boxSize = opt.box.size;
  dp.maxOverlapAtomNum = opt.maxOverlapAtomNum;
  dp.maxNearDistance = opt.maxNearDistance;
  dp.clashDistance = opt.clashDistance;
  dp.clashBase = opt.EclashBase;
  dp.clashK = opt.Kclash;
  dp.minimizer.method = od::parseMinMethod(opt.dockMinMethod);
  dp.minimizer.maxSteps = opt.dockMinMaxStep;
  dp.minimizer.convSteps = opt.dockMinConvSteps;
  dp.minimizer.deltaE = opt.dockMinDeltaE;
  dp.minimizer.memory = opt.lbfgsMemory;
  dp.minimizer.lineSearchInit = od::parseLineSearchInit(opt.dockLineSearchInit);
  dp.minimizer.lineSearchReduce = od::parseLineSearchReduce(opt.dockLineSearchReduce);
  dp.minimizer.maxBacktrack = opt.dockMaxBacktrack;
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
  int topNum = 0;  // from the initial eliteConfSampleNum values; not changed with the dcss size
  // dcss size (spec 9, step 6): host memory of a prepared pair, docking samples per docked pair (last round with docked pairs).
  double preparedPerPair = 0, samplesPerPair = 0;
  std::size_t bufferNum = 0;  // sampleBufferNum
  std::unique_ptr<od::MobmHeavyAtoms> heavy;
  std::unique_ptr<od::SymmetryRmsd> sym;
  std::vector<od::Vec3f> ref;  // reference heavy atoms (box frame) in MobmHeavyAtoms order; empty without --refFiles
  std::vector<int> refSel;     // MobmHeavyAtoms positions present in the reference; empty: all
  std::unique_ptr<od::SymmetryRmsd> refSym;  // SymmetryRmsd over refSel (refSel not empty)
  // score / localMin (spec 6.1), per entry of top: record (from 1), fit RMSD of the record, grid status.
  struct Input {
    int record = 0;
    float fitRmsd = 0;
    std::string gridStatus;
  };
  std::vector<Input> inputs;
};

// Reference heavy atoms (spec 10.3) in MobmHeavyAtoms order (ref in the box frame); MobM atoms
// absent from the reference (C-terminal OXT) are left out, with m.refSel and m.refSym over the rest.
void setReference(MobState& m, const ConfSet& mob, const od::Molecule& ref, const od::Options& opt) {
  const std::vector<int> map = od::referenceAtomMap(*mob.parts()[0].mol, m.heavy->atoms(), ref);
  m.ref.clear();
  m.refSel.clear();
  m.refSym.reset();
  std::vector<int> atoms, types;
  for (std::size_t k = 0; k < map.size(); ++k) {
    if (map[k] < 0) continue;
    const od::Vec3& p = ref.atoms[map[k]].xyz;
    m.ref.push_back({static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])});
    m.refSel.push_back(static_cast<int>(k));
    atoms.push_back(m.heavy->atoms()[k]);
    types.push_back(m.heavy->types()[k]);
  }
  if (m.refSel.size() == map.size()) {
    m.refSel.clear();
    return;
  }
  m.refSym = std::make_unique<od::SymmetryRmsd>(*mob.parts()[0].mol, atoms, types, static_cast<std::size_t>(opt.symMaxPerm));
  od::logInfo(std::format("{}: {} atom(s) absent from the reference (C-terminal OXT) are left out of rmsd_ref", mob.name(),
                          map.size() - m.refSel.size()));
}

double refRmsd(const MobState& m, const std::vector<od::Vec3f>& world) {
  if (m.refSel.empty()) return m.sym->rmsd(world, m.ref);
  std::vector<od::Vec3f> w;
  for (int k : m.refSel) w.push_back(world[k]);
  return m.refSym->rmsd(w, m.ref);
}

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
  w.put(m.preparedPerPair);
  w.put(m.samplesPerPair);
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
  m.preparedPerPair = r.get<double>();
  m.samplesPerPair = r.get<double>();
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
// Output coordinates in the input frame (spec 3.3).
std::vector<od::StructPart> structParts(const od::EnergySystem& sys, const std::vector<od::EnergyPart>& parts,
                                        const std::vector<od::Vec3f>& x, const od::BoxFrame& box) {
  std::vector<od::StructPart> out;
  for (std::size_t k = 0; k < parts.size(); ++k) {
    const int off = sys.coordOffset(static_cast<int>(k));
    const std::size_t n = parts[k].mol->atoms.size();
    out.push_back({parts[k].mol, std::vector<od::Vec3f>(x.begin() + off, x.begin() + off + static_cast<long>(n))});
    if (box.rotated)
      for (od::Vec3f& p : out.back().xyz) {
        const od::Vec3 q = box.toInput({p[0], p[1], p[2]});
        p = {static_cast<float>(q[0]), static_cast<float>(q[1]), static_cast<float>(q[2])};
      }
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
           const od::Options& opt, const od::FeParams& fp, float shift, bool gpu) {
  const od::EnergyParams ep = energyParams(opt);
  od::ShakeParams sp;
  sp.radius = fp.shakingRadius;
  sp.step = fp.shakingStep;
  sp.temperature = fp.temperature;
  sp.rrtScaling = fp.rrtScaling;
  // Atom shaking of od_0100 on the GPU (spec 15.1.17).
  od::ShakeFn shakeFn;
#ifdef OD_HAVE_CUDA
  std::shared_ptr<od::GpuShaker> gpuShaker;
  if (gpu && opt.feOd0100) {
    gpuShaker = std::make_shared<od::GpuShaker>(ep, sp, opt.gpuDevices[0]);
    shakeFn = [gpuShaker](const od::ShakeAtoms& a, std::span<const char> s) { return gpuShaker->lnF(a, s); };
  }
#else
  (void)gpu;
#endif
  const auto lnF = [&](const od::ShakeAtoms& a) { return shakeFn ? shakeFn(a, {}) : od::shakeLnF(a, ep, sp); };
  const od::OutFormats fmt{opt.outPdb, opt.outCif, opt.outSdf, opt.outMol2};

  // Protonation states (protonation.recompute, spec 3.4.1): l2, the variants of the samples; final,
  // recalculated for each structure of 9d (complex, and each unbound candidate alone).
  const bool recompute = opt.protonationRecompute != "none";
  const bool recFinal = opt.protonationRecompute == "final";
  // PROPKA molecules of a set at coordinates x of its base system.
  const auto pkaMols = [](const ConfSet& set, std::span<const od::Vec3f> x) {
    if (set.prot()) return set.prot()->pkaMolecules(x);
    std::vector<od::PkaMolecule> v;
    for (std::size_t k = 0; k < set.parts().size(); ++k) {
      od::PkaMolecule p;
      p.mol = set.parts()[k].mol;
      const int off = set.sys().coordOffset(static_cast<int>(k));
      for (std::size_t i = 0; i < p.mol->atoms.size(); ++i) p.xyz.push_back({x[off + i][0], x[off + i][1], x[off + i][2]});
      v.push_back(std::move(p));
    }
    return v;
  };
  // Variant and protonation free energy of a set from the pKa values of its parts.
  struct ProtState {
    int variant = 0;
    double g = 0;
  };
  using Changes = std::vector<std::vector<od::StateChange>>;
  const auto changesFrom = [&](const ConfSet& set, std::span<const od::PkaMap> pka) {
    Changes ch(set.prot() ? set.prot()->partCount() : 0);
    for (std::size_t p = 0; p < ch.size(); ++p) ch[p] = od::stateChanges(set.prot()->sites(static_cast<int>(p)), pka[p], opt.protonationPH);
    return ch;
  };
  const auto stateOf = [&](const ConfSet& set, const Changes& ch) {
    ProtState s;
    if (!set.prot()) return s;
    for (std::size_t p = 0; p < ch.size(); ++p)
      s.g += od::protonationFreeEnergy(set.prot()->sites(static_cast<int>(p)), ch[p], opt.protonationPH, fp.temperature);
    s.variant = set.prot()->variantOf(ch);
    return s;
  };
  // Energy alone (9a terms) of a conformation with the system of a protonation variant.
  const auto aloneEnergy = [](const std::vector<float>& angles) {
    return [&angles](const od::EnergySystem& sys, const std::vector<od::EnergyPart>&) {
      return sys.unweighted(sys.evaluate(std::span<const float>(angles), nullptr, nullptr, nullptr)).filterEnergy();
    };
  };
  // A sample alone: l2, its variant; final, from PROPKA of the set alone (neutral His: the tautomer
  // of the lowest energy alone).
  const auto aloneState = [&](const ConfSet& set, const od::ConfSample& s) {
    if (!recFinal || !set.prot()) return ProtState{s.variant, s.protG};
    const auto x = set.sys().coordinates(std::span<const float>(s.angles));
    const auto pka = od::calculatePka(pkaMols(set, x));
    Changes ch = changesFrom(set, pka);
    set.prot()->chooseTautomers(ch, aloneEnergy(s.angles));
    return stateOf(set, ch);
  };

  // A conformation evaluated alone: atoms, 9a terms (single precision) and, when needed, ln f of its atoms.
  struct Alone {
    std::vector<od::Vec3f> x;
    od::ShakeAtoms atoms;
    od::EnergyTerms terms;
    std::vector<double> lnF;
    double sumLnF = 0;
    ProtState prot;
  };
  auto evalAlone = [&](const ConfSet& set, const std::vector<float>& angles, ProtState ps, bool shake) {
    Alone a;
    a.prot = ps;
    const od::EnergySystem& sys = set.sysOf(ps.variant);
    a.x = sys.coordinates(angles);
    a.atoms = od::shakeAtoms(sys, set.partsOf(ps.variant), a.x);
    a.terms = sys.unweighted(sys.evaluate(std::span<const float>(angles), nullptr, nullptr, nullptr));
    if (shake) {
      a.lnF = lnF(a.atoms);
      a.sumLnF = std::accumulate(a.lnF.begin(), a.lnF.end(), 0.0);
    }
    return a;
  };
  std::map<std::pair<int, std::vector<float>>, Alone> stamAlone;  // by variant and angles
  auto stamOf = [&](const std::vector<float>& angles, ProtState ps) -> const Alone& {
    auto it = stamAlone.find({ps.variant, angles});
    if (it == stamAlone.end()) it = stamAlone.emplace(std::pair{ps.variant, angles}, evalAlone(stamSet, angles, ps, opt.feOd0100)).first;
    return it->second;
  };
  // Unbound state: the candidate with the lowest G = sum_t intra[t] E_intra,t + entropy (-TS) + protonation free energy.
  struct Unbound {
    double g = 0, mts = 0, protG = 0, ring = 0;
    Terms4 intra{};
  };
  auto unbound = [&](const ConfSet& set, bool stam) {
    Unbound u;
    u.g = std::numeric_limits<double>::infinity();
    const double cS = stam ? fp.entropyStam : fp.entropyMobm;
    for (const auto& smp : set.lowestSamples(fp.unboundCandidateNum)) {
      const ProtState ps = aloneState(set, smp);
      const Alone tmp = stam ? Alone{} : evalAlone(set, smp.angles, ps, true);
      const Alone& a = stam ? stamOf(smp.angles, ps) : tmp;
      const Terms4 e = terms4(a.terms);
      const double mts = od::minusTS(a.sumLnF, sp);
      const double ring = set.sysOf(ps.variant).ringEnergy(std::span<const float>(smp.angles));
      const double g = dot4(fp.intra, e) + cS * mts + ps.g + fp.ring * ring;
      if (g < u.g) u = {g, mts, ps.g, ring, e};
    }
    return u;
  };
  Unbound ubStam;
  if (opt.feOd0100) {
    const auto t0 = std::chrono::steady_clock::now();
    ubStam = unbound(stamSet, true);
    od::logInfo(std::format("9d od_0100: G_unbound StaM group {:.3f}  ({} candidate(s), {:.2f} s)", ubStam.g,
                            stamSet.lowestSamples(fp.unboundCandidateNum).size(),
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
      double dgAd = 0, dgOd = 0, gComplex = 0, mtsStam = 0, mtsMobm = 0, rmsdRef = 0, dRing = 0;
      Terms4 inter4{}, dIntra4{};
      std::vector<od::Vec3f> world;
      int stamVariant = 0, mobVariant = 0;
      double protG = 0;  // protonation free energy of the complex (relative to the input protonation)
    };
    std::vector<Row> rows(top.size());
    for (std::size_t k = 0; k < top.size(); ++k) {
      const auto& r = top[k];
      Row& w = rows[k];
      const od::ConfSample& ss = stamSet.sample(r.stam);
      const od::ConfSample& msmp = mob.sample(r.mobm);
      ProtState ps{ss.variant, ss.protG}, pm{msmp.variant, msmp.protG};
      if (recFinal && (stamSet.prot() || mob.prot())) {
        // The complex: StaM conformation and MobM pose with the input protonation.
        const auto xs = stamSet.sys().coordinates(std::span<const float>(ss.angles));
        const auto xm = od::dockWorld(mob.sys().coordinates(std::span<const float>(r.s.angles)), r.s);
        auto mols = pkaMols(stamSet, xs);
        const std::size_t ns = mols.size();
        for (auto& m : pkaMols(mob, xm)) mols.push_back(std::move(m));
        const auto pka = od::calculatePka(mols);
        Changes cs = changesFrom(stamSet, std::span<const od::PkaMap>(pka).first(ns));
        Changes cm = changesFrom(mob, std::span<const od::PkaMap>(pka).subspan(ns));
        // Neutral His: the tautomers of the lowest complex energy (energy alone of the group plus the
        // intermolecular energy), the StaM group and the MobM in turn.
        if (stamSet.prot() || mob.prot()) {
          od::ShakeAtoms sa = od::shakeAtoms(stamSet.sys(), stamSet.parts(), xs);
          od::ShakeAtoms ma = od::shakeAtoms(mob.sys(), mob.parts(), xm);
          for (int turn = 0; turn < od::ProtonationVariants::kTautomerPasses; ++turn) {
            int flips = 0;
            if (stamSet.prot())
              flips += stamSet.prot()->chooseTautomers(cs, [&](const od::EnergySystem& sys, const std::vector<od::EnergyPart>& parts) {
                const auto x = sys.coordinates(std::span<const float>(ss.angles));
                const od::InterTerms it = od::exactInter(od::shakeAtoms(sys, parts, x), ma, ep);
                return sys.unweighted(sys.evaluate(std::span<const float>(ss.angles), nullptr, nullptr, nullptr)).filterEnergy() + it.inter() + it.clash;
              });
            if (stamSet.prot()) {
              const od::EnergySystem& sys = stamSet.prot()->sysOf(cs);
              sa = od::shakeAtoms(sys, stamSet.prot()->partsOf(cs), sys.coordinates(std::span<const float>(ss.angles)));
            }
            if (mob.prot() && (turn == 0 || flips > 0))
              flips += mob.prot()->chooseTautomers(cm, [&](const od::EnergySystem& sys, const std::vector<od::EnergyPart>& parts) {
                const auto x = od::dockWorld(sys.coordinates(std::span<const float>(r.s.angles)), r.s);
                const od::InterTerms it = od::exactInter(sa, od::shakeAtoms(sys, parts, x), ep);
                return sys.unweighted(sys.evaluate(std::span<const float>(r.s.angles), nullptr, nullptr, nullptr)).filterEnergy() + it.inter() + it.clash;
              });
            if (mob.prot()) {
              const od::EnergySystem& sys = mob.prot()->sysOf(cm);
              ma = od::shakeAtoms(sys, mob.prot()->partsOf(cm), od::dockWorld(sys.coordinates(std::span<const float>(r.s.angles)), r.s));
            }
            if (flips == 0) break;
          }
        }
        ps = stateOf(stamSet, cs);
        pm = stateOf(mob, cm);
      }
      w.stamVariant = ps.variant;
      w.mobVariant = pm.variant;
      w.protG = ps.g + pm.g;
      const Alone& st = stamOf(ss.angles, ps);
      const od::EnergySystem& msys = mob.sysOf(pm.variant);
      const auto y = msys.coordinates(std::span<const float>(r.s.angles));
      w.world = od::dockWorld(y, r.s);
      if (!ms[i].ref.empty()) w.rmsdRef = refRmsd(ms[i], ms[i].heavy->world(r.s));
      const od::ShakeAtoms mobAtoms = od::shakeAtoms(msys, mob.partsOf(pm.variant), w.world);
      const od::EnergyTerms mt = msys.unweighted(msys.evaluate(std::span<const float>(r.s.angles), nullptr, nullptr, nullptr));
      const od::InterTerms it = od::exactInter(st.atoms, mobAtoms, ep);
      w.inter = it.inter();
      w.interClash = it.clash;
      w.mobIntra = mt.intra();
      w.mobClash = mt.clash;
      // final: the StaM term of 9b with the energy change of the StaM variant; l2: the StaM sample
      // energy (9b) includes its protonation free energy.
      w.stamTerm = r.stamRel + shift;
      if (recFinal && ps.variant != 0) w.stamTerm += st.terms.filterEnergy() - stamOf(ss.angles, {}).terms.filterEnergy();
      w.stamIntra = st.terms.intra();
      w.stamClash = st.terms.clash;
      w.e9b = w.inter + w.interClash + w.mobIntra + w.mobClash + w.stamTerm + (recFinal ? w.protG : pm.g);
      w.e9bGrid = r.rank() + shift;
      w.dgAd = w.inter + od::kAdgpuTorsionCoeff * nTors + w.protG;
      if (opt.feOd0100) {
        const od::ComplexLnF l = od::complexSumLnF(st.atoms, st.lnF, mobAtoms, ep, sp, shakeFn);
        w.mtsStam = od::minusTS(l.stam, sp);
        w.mtsMobm = od::minusTS(l.mobm, sp);
        w.inter4 = terms4(it);
        const Terms4 si = terms4(st.terms), mi = terms4(mt);
        Terms4 intra4{};
        for (int t = 0; t < od::FeParams::kTerms; ++t) {
          intra4[t] = si[t] + mi[t];
          w.dIntra4[t] = intra4[t] - ubStam.intra[t] - ubMob.intra[t];
        }
        const double ring = stamSet.sysOf(ps.variant).ringEnergy(std::span<const float>(ss.angles)) +
                            msys.ringEnergy(std::span<const float>(r.s.angles));
        w.dRing = ring - ubStam.ring - ubMob.ring;
        w.gComplex = dot4(fp.inter, w.inter4) + dot4(fp.intra, intra4) + fp.entropyStam * w.mtsStam +
                     fp.entropyMobm * w.mtsMobm + w.protG + fp.ring * ring;
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
    // score / localMin: every pose output, numbered and listed in record order (spec 6.1).
    const bool inputPoses = !ms[i].inputs.empty();
    std::vector<std::size_t> listOrder = order[0];
    if (inputPoses) {
      std::fill(out.begin(), out.end(), 1);
      std::iota(listOrder.begin(), listOrder.end(), 0);
    }
    std::vector<int> pose(top.size(), 0);  // pose number, by the 9b ranking
    int poseNum = 0;
    for (std::size_t k : listOrder)
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
    // StaM structures by sample and protonation variant (stam_<sample> or stam_<sample>_p<variant>).
    std::map<std::pair<int, int>, std::vector<std::string>> stamFiles;
    const auto stamKey = [&](std::size_t k) { return std::pair{top[k].stam, rows[k].stamVariant}; };
    std::vector<std::vector<std::string>> poseFiles(top.size());
    for (std::size_t k : listOrder) {
      if (!out[k]) continue;
      const auto& r = top[k];
      const int sv = rows[k].stamVariant;
      if (!stamFiles.count(stamKey(k))) {
        const Alone& st = stamOf(stamSet.sample(r.stam).angles, {sv, 0.0});
        bool fb = false;
        stamFiles[stamKey(k)] = od::writeStructure(
            dir + (sv == 0 ? std::format("/stam_{}", r.stam) : std::format("/stam_{}_p{}", r.stam, sv)),
            structParts(stamSet.sysOf(sv), stamSet.partsOf(sv), st.x, opt.box), fmt,
            sv == 0 ? std::format("StaM conformation sample {}", r.stam)
                    : std::format("StaM conformation sample {}, protonation {}", r.stam, stamSet.prot()->describe(sv)),
            fb);
        fallback += fb;
      }
      bool fb = false;
      const int mv = rows[k].mobVariant;
      poseFiles[k] = od::writeStructure(dir + std::format("/pose_{:04}", pose[k]),
                                        structParts(mob.sysOf(mv), mob.partsOf(mv), rows[k].world, opt.box), fmt,
                                        inputPoses
                                            ? std::format("{} pose {} (StaM sample {}, MobM record {}) E_9b {:.3f}",
                                                          mob.name(), pose[k], r.stam, ms[i].inputs[k].record,
                                                          rows[k].e9b)
                                            : std::format("{} pose {} (StaM sample {}, MobM sample {}) E_9b {:.3f}",
                                                          mob.name(), pose[k], r.stam, r.mobm, rows[k].e9b),
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
      if (inputPoses) f << "\tfit_rmsd\tgrid_status";
      if (!ms[i].ref.empty()) f << "\trmsd_ref";
      for (const auto& m : methods) f << "\trank_" << m.name;
      f << "\tE_9b\tE_9b_grid\tinter\tinter_clash\tmobm_intra\tmobm_clash\tstam_term\tstam_intra\tstam_clash\tclash";
      if (recompute) f << "\tdG_prot\tprot_stam\tprot_mobm";
      if (opt.feAdgpu16) f << "\tdG_adgpu_16";
      if (opt.feOd0100) {
        f << "\tdG_od_0100\tG_complex\tG_unbound_stam\tG_unbound_mobm";
        for (const char* t : od::kFeTermNames) f << "\tinter_" << t;
        for (const char* t : od::kFeTermNames) f << "\tdintra_" << t;
        f << "\tmTS_complex_stam\tmTS_complex_mobm\tmTS_unbound_stam\tmTS_unbound_mobm\tn_tors\tdring";
      }
      f << "\n";
      for (std::size_t k : listOrder) {
        if (!out[k] && !opt.outAllSummary) continue;
        const Row& w = rows[k];
        const bool written = out[k] != 0;
        f << std::format("{}\t{}\t{}\t{}\t{}", written ? std::to_string(pose[k]) : "-",
                         written ? joinNames(poseFiles[k]) : "-", written ? joinNames(stamFiles[stamKey(k)]) : "-",
                         top[k].stam, inputPoses ? ms[i].inputs[k].record : top[k].mobm);
        if (inputPoses) f << std::format("\t{:.4f}\t{}", ms[i].inputs[k].fitRmsd, ms[i].inputs[k].gridStatus);
        if (!ms[i].ref.empty()) f << std::format("\t{:.4f}", w.rmsdRef);
        for (std::size_t m = 0; m < methods.size(); ++m) f << '\t' << rankOf[m][k];
        f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}", w.e9b, w.e9bGrid,
                         w.inter, w.interClash, w.mobIntra, w.mobClash, w.stamTerm, w.stamIntra, w.stamClash,
                         w.interClash + w.mobClash + w.stamClash > 0 ? "yes" : "no");
        if (recompute) {
          const std::string a = w.stamVariant ? stamSet.prot()->describe(w.stamVariant) : "-";
          const std::string b = w.mobVariant ? mob.prot()->describe(w.mobVariant) : "-";
          f << std::format("\t{:.4f}\t{}\t{}", w.protG, a, b);
        }
        if (opt.feAdgpu16) f << std::format("\t{:.4f}", w.dgAd);
        if (opt.feOd0100) {
          f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}", w.dgOd, w.gComplex, ubStam.g, ubMob.g);
          for (double v : w.inter4) f << std::format("\t{:.4f}", v);
          for (double v : w.dIntra4) f << std::format("\t{:.4f}", v);
          f << std::format("\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\t{:.4f}", w.mtsStam, w.mtsMobm, ubStam.mts, ubMob.mts, nTors,
                           w.dRing);
        }
        f << "\n";
      }
    }
    for (std::size_t m = 0; m < methods.size(); ++m) {
      std::ofstream f = openTsv("rank_" + methods[m].name + ".tsv");
      f << "rank\tpose\tmobm_files\tstam_files\tenergy\tclash" << (ms[i].ref.empty() ? "" : "\trmsd_ref") << "\n";
      for (std::size_t k : order[m]) {
        if (!out[k]) continue;
        const Row& w = rows[k];
        f << std::format("{}\t{}\t{}\t{}\t{:.4f}\t{}", rankOf[m][k], pose[k], joinNames(poseFiles[k]),
                         joinNames(stamFiles[stamKey(k)]), w.*methods[m].energy,
                         w.interClash + w.mobClash + w.stamClash > 0 ? "yes" : "no");
        if (!ms[i].ref.empty()) f << std::format("\t{:.4f}", w.rmsdRef);
        f << "\n";
      }
    }

    od::logInfo(std::format("9d {:<38} recalculated {} sample(s) ({:.2f} s){}; output {} pose(s), {} StaM conformation(s) to {}",
                            mob.name(), top.size(), evalSec,
                            opt.feOd0100 ? std::format(", G_unbound MobM {:.3f} ({:.2f} s)", ubMob.g, unboundSec) : "",
                            poseNum, stamFiles.size(), dir));
    if (recompute) {
      const auto changed = std::count_if(rows.begin(), rows.end(), [](const Row& w) { return w.stamVariant || w.mobVariant; });
      od::logInfo(std::format("    protonation.recompute {}: {} of {} sample(s) with changed states", opt.protonationRecompute,
                              changed, rows.size()));
    }
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
          const od::ConfSample& ss = stamSet.sample(r.stam);
          const od::EnergySystem& sys = stamSet.sysOf(ss.variant);
          od::writeStructure(d + std::format("/stam_{}", r.stam),
                             structParts(sys, stamSet.partsOf(ss.variant), sys.coordinates(std::span<const float>(ss.angles)), opt.box),
                             fmt, std::format("StaM conformation sample {}", r.stam), fb);
        }
        const int mv = mob.sample(r.mobm).variant;
        const auto y = mob.sysOf(mv).coordinates(std::span<const float>(r.s.angles));
        files = joinNames(od::writeStructure(d + std::format("/top_{:04}", k + 1),
                                             structParts(mob.sysOf(mv), mob.partsOf(mv), od::dockWorld(y, r.s), opt.box), fmt,
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
                     const std::vector<MobState>& ms, float stamRef, bool limitStop, std::uint64_t memBase,
                     const od::Options& opt, const Checkpointing& ck, float shift) {
  const auto t0 = std::chrono::steady_clock::now();
  od::BinWriter w;
  od::putStateHeader(w);
  w.put<std::int32_t>(round);
  w.put<std::uint64_t>(mobSets.size());
  stamSet.save(w);
  for (const auto& m : mobSets) m->save(w);
  w.put(stamRef);
  w.put<std::uint8_t>(limitStop);
  w.put(memBase);
  for (const auto& m : ms) saveMobState(w, m);
  const od::Manifest man{od::kVersion, round, ck.inputs, od::sha256Hex(w.data())};
  const std::string dir = od::commitCheckpoint(
      opt.outPrefix, round, {{"state.bin", w.data()}, {"params.toml", ck.params}, {"manifest.toml", od::manifestText(man)}},
      [&](const std::string& tmp) { writeTopSnapshot(tmp, stamSet, mobSets, ms, opt, shift); });
  od::logInfo(std::format("  checkpoint {} ({:.1f} MB state)  {:.2f} s", dir, static_cast<double>(w.data().size()) / 1e6,
                          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
}

// Automatic box (spec 3.3): the grids of one StaM conformation (type maps of the MobM probe
// types, electrostatic, desolvation and distance grid, float) must fit in the free host memory
// and in the free memory of each device.
void checkGridMemory(const od::GridGeometry& geom, std::size_t probeTypes, const od::Options& opt, bool gpu) {
  const std::size_t need = geom.size() * (probeTypes + 3) * sizeof(float);
  const auto fail = [&](const std::string& where, std::size_t avail) {
    od::fail("the automatic docking box (box.mode = {}; size {:.1f} x {:.1f} x {:.1f} Å, {} x {} x {} = {} grid points, "
             "{} maps) needs {} MiB for the grids of one StaM conformation, more than the {} ({} MiB); increase "
             "grid.spacing or use box.mode = manual with a smaller box.size",
             opt.boxMode, opt.box.size[0], opt.box.size[1], opt.box.size[2], geom.n[0], geom.n[1], geom.n[2], geom.size(),
             probeTypes + 3, need >> 20, where, avail >> 20);
  };
  const od::HostMemLimit lim = od::hostMemLimit(opt.memHostLimit);
  const std::size_t resident = od::residentBytes();
  const std::size_t host = lim.bytes > resident ? lim.bytes - resident : 0;
  if (lim.bytes > 0 && need > host) fail(std::format("free host memory of memory.hostLimit ({})", lim.source), host);
#ifdef OD_HAVE_CUDA
  if (gpu)
    for (int dev : opt.gpuDevices)
      if (const std::size_t f = od::gpuFreeMem(dev); need > f) fail(std::format("free memory of GPU {}", dev), f);
#else
  (void)gpu;
#endif
}

// gpu.dock9bStamNum = auto (spec 15.1.16): L2poss of a window per pose evaluated at the same time on the GPU.
constexpr long long kWindowPosesPerResident = 48;
// posm2_fft: most StaM conformations whose translations are found together (spec 15.1.22).
constexpr std::size_t kFftStamMax = 16;

// Rounds of 9a refill, 9b docking and 9c selection (spec 6, 7.3, 8, 9).
void runRounds(ConfSet& stamSet, std::vector<std::unique_ptr<ConfSet>>& mobSets, const std::vector<od::Molecule>& refs,
               const od::Options& opt, const od::FeParams& fp, std::uint64_t seed, const Checkpointing& ck, bool gpu) {
  std::vector<od::Ad4Type> probe;
  for (const auto& m : mobSets)
    for (od::Ad4Type t : m->sys().typing(0).type)
      if (t != od::Ad4Type::Merged && std::find(probe.begin(), probe.end(), t) == probe.end()) probe.push_back(t);
  const od::GridGeometry geom = od::boxGrid(opt.box.center, opt.box.size, opt.gridSpacing);
  if (opt.boxMode != "manual") checkGridMemory(geom, probe.size(), opt, gpu);
  {
    const int level = od::hopfLevelForStep(opt.rotStep);
    od::logInfo(std::format("9b grids: {} x {} x {} points, spacing {} Å; orientations: Hopf level {}, {} "
                            "orientations, spacing {} deg (dock9b.rotStep {})",
                            geom.n[0], geom.n[1], geom.n[2], geom.spacing, level, 72LL << (3 * level),
                            od::hopfSpacingDeg(level), opt.rotStep));
  }
  const od::Dock9bParams dp = dockParams(opt, seed);
  // 9b on the GPU (spec 15.1.4, 15.1.5): one evaluator per MobM; field, quick filter and posm2_fft shared by them.
  std::vector<std::unique_ptr<od::PoseEvaluator>> poseEval(mobSets.size());
  // Several devices (spec 15.1.33): the windows of StaM conformations divided among them
  // (windowDevices); poseEvalOthers[d - 1]: the evaluators of device d.
  std::vector<std::vector<std::unique_ptr<od::PoseEvaluator>>> poseEvalOthers;
  int windowDevices = 1;
  const auto evalOf = [&](int d, std::size_t i) { return d == 0 ? poseEval[i].get() : poseEvalOthers[d - 1][i].get(); };
  const double distCap = std::max({opt.maxNearDistance, opt.clashDistance,
                                   opt.posm == "posm3" || opt.posm == "posm4" ? opt.sepaShellMax : 0.0}) +
                         2 * opt.gridSpacing;
  std::function<od::EnergyGrids(const od::GridAtoms&)> energyGrids = [&](const od::GridAtoms& a) {
    return od::buildEnergyGrids(geom, a, probe, searchParams(opt, true));
  };
  std::function<od::DistanceGrid(const od::GridAtoms&)> distanceGrid = [&](const od::GridAtoms& a) {
    return od::buildDistanceGrid(geom, a, distCap);
  };
  // StaM conformations docked together (spec 15.1.16) and the device memory for their fields.
  int stamWindow = 1;
  std::function<std::size_t()> gpuFieldBudget = [] { return std::numeric_limits<std::size_t>::max(); };
  std::function<void(int, int)> gpuFieldSlots = [](int, int) {};  // (device index, slots)
  std::function<std::size_t()> gpuFree = [] { return std::numeric_limits<std::size_t>::max(); };
#ifdef OD_HAVE_CUDA
  // Grids on the GPU (spec 15.1.10): the same values as the CPU.
  std::shared_ptr<od::GpuGridBuilder> gpuGrid;
  std::vector<std::unique_ptr<od::GpuGridBuilder>> gpuGridOthers;  // delta blocks on the other devices (spec 15.1.33)
  std::vector<std::shared_ptr<od::GpuField>> gpuFields;
  // Minimizations taken over by idle devices (spec 15.1.35); gpuEvals[d][i]: the evaluator of
  // MobM i on device d.
  std::shared_ptr<od::DockStealHub> stealHub;
  std::vector<std::vector<od::GpuPoseEvaluator*>> gpuEvals;
  od::DockStealHub::Stats stealBase;  // statistics at the start of the round
  std::vector<std::shared_ptr<od::GpuDockPrep>> gpuPreps;
  if (gpu) {
    gpuGrid = std::make_shared<od::GpuGridBuilder>(opt.gpuDevices[0]);
    for (std::size_t k = 1; k < opt.gpuDevices.size(); ++k)
      gpuGridOthers.push_back(std::make_unique<od::GpuGridBuilder>(opt.gpuDevices[k]));
    distanceGrid = [&](const od::GridAtoms& a) {
      od::DistanceGrid g;
      g.geom = geom;
      g.cap = distCap;
      g.d.resize(geom.size());
      gpuGrid->buildDistance(od::DistBuildData(geom, a, distCap).view, g.d.data());
      return g;
    };
    energyGrids = [&](const od::GridAtoms& a) {
      const od::GridBuildData b(geom, a, probe, searchParams(opt, true));
      od::EnergyGrids g;
      g.geom = geom;
      g.types = b.types;
      g.typeMaps.assign(g.types.size(), std::vector<float>(geom.size()));
      g.elecMap.resize(geom.size());
      g.desolvMap.resize(geom.size());
      std::vector<float*> maps;
      for (auto& m : g.typeMaps) maps.push_back(m.data());
      gpuGrid->build(b.view, od::kGridClamp, maps.data(), g.elecMap.data(), g.desolvMap.data());
      return g;
    };
    // One field per device (spec 15.1.33); the quick filter and posm2_fft on the first device.
    for (int d : opt.gpuDevices) gpuFields.push_back(std::make_shared<od::GpuField>(d));
    stamWindow = opt.gpuDock9bStamNum > 0 ? opt.gpuDock9bStamNum : std::numeric_limits<int>::max();
    // Resident fields: those of the window (the fields of the previous window are dropped); with
    // the windows divided among the devices, on the device of the window.
    gpuFieldSlots = [gpuFields, &windowDevices](int d, int n) {
      if (windowDevices > 1) gpuFields[d]->setSlots(n);
      else
        for (const auto& f : gpuFields) f->setSlots(n);
    };
    // Half of the free device memory and the memory of the resident fields (the smallest of the devices).
    gpuFieldBudget = [&opt, gpuFields] {
      std::size_t b = std::numeric_limits<std::size_t>::max();
      for (std::size_t k = 0; k < gpuFields.size(); ++k)
        b = std::min(b, (od::gpuFreeMem(opt.gpuDevices[k]) + gpuFields[k]->deviceBytes()) / 2);
      return b;
    };
    gpuFree = [&opt] { return od::gpuFreeMem(opt.gpuDevices[0]); };
    // Several devices: each window of StaM conformations on one device, with its own quick filter;
    // posm2_fft (spectra of several StaM conformations shared, spec 15.1.24): the translations and
    // the minimizations of every window divided among the devices (MultiGpuPoseEvaluator).
    const bool splitWindows = gpuFields.size() > 1 && dp.posm != od::Posm::Posm2Fft;
    for (std::size_t k = 0; k < gpuFields.size(); ++k)
      gpuPreps.push_back(std::make_shared<od::GpuDockPrep>(opt.gpuDevices[k]));
    if (splitWindows) {
      windowDevices = static_cast<int>(gpuFields.size());
      poseEvalOthers.resize(gpuFields.size() - 1);
      for (auto& v : poseEvalOthers) v.resize(mobSets.size());
    }
    if (gpuFields.size() > 1) stealHub = std::make_shared<od::DockStealHub>();
    gpuEvals.assign(gpuFields.size(), std::vector<od::GpuPoseEvaluator*>(mobSets.size(), nullptr));
    for (std::size_t i = 0; i < mobSets.size(); ++i) {
      std::vector<std::unique_ptr<od::GpuPoseEvaluator>> ev;
      for (std::size_t k = 0; k < gpuFields.size(); ++k) {
        ev.push_back(std::make_unique<od::GpuPoseEvaluator>(mobSets[i]->sys().model(), gpuFields[k], opt.gpuDevices[k],
                                                            opt.gpuDock9bBatchNum, gpuPreps[k], opt.gpuDock9bDeviceMin));
        if (stealHub) ev.back()->setSteal(stealHub, static_cast<int>(i));
        gpuEvals[k][i] = ev.back().get();
      }
      if (ev.size() == 1) poseEval[i] = std::move(ev[0]);
      else if (splitWindows) {
        poseEval[i] = std::move(ev[0]);
        for (std::size_t k = 1; k < ev.size(); ++k) poseEvalOthers[k - 1][i] = std::move(ev[k]);
      } else poseEval[i] = std::make_unique<od::MultiGpuPoseEvaluator>(std::move(ev), stealHub);
    }
  }
#else
  (void)gpu;
#endif
  // Evaluators of the protonation variants of MobM samples (protonation.recompute = l2, spec 3.4.1),
  // made when first needed; they do not take over minimizations of other devices.
  struct VariantKey {
    int dev = 0;
    std::size_t i = 0;
    int v = 0;
    auto operator<=>(const VariantKey&) const = default;
  };
  std::map<VariantKey, std::unique_ptr<od::PoseEvaluator>> variantEval;
  const auto evalFor = [&](int d, std::size_t i, int v) -> od::PoseEvaluator* {
    if (v == 0 || !evalOf(d, i)) return evalOf(d, i);
    auto& e = variantEval[{d, i, v}];
#ifdef OD_HAVE_CUDA
    if (!e) {
      const od::EnergyModel model = mobSets[i]->sysOf(v).model();
      if (windowDevices > 1 || gpuFields.size() == 1) {
        e = std::make_unique<od::GpuPoseEvaluator>(model, gpuFields[d], opt.gpuDevices[d], opt.gpuDock9bBatchNum, gpuPreps[d],
                                                   opt.gpuDock9bDeviceMin);
      } else {
        std::vector<std::unique_ptr<od::GpuPoseEvaluator>> ev;
        for (std::size_t k = 0; k < gpuFields.size(); ++k)
          ev.push_back(std::make_unique<od::GpuPoseEvaluator>(model, gpuFields[k], opt.gpuDevices[k], opt.gpuDock9bBatchNum,
                                                              gpuPreps[k], opt.gpuDock9bDeviceMin));
        e = std::make_unique<od::MultiGpuPoseEvaluator>(std::move(ev));
      }
    }
#endif
    return e.get();
  };
  // Host memory for the fields of a window (memory.hostLimit, spec 13, 15.1.16): 90% of the limit.
  const od::HostMemLimit hostLimit = od::hostMemLimit(opt.memHostLimit);
  const std::size_t hostBudget = hostLimit.bytes / 10 * 9;
  if (stamWindow > 1)
    od::logInfo(std::format("9b memory.hostLimit {} MiB ({}); StaM conformations docked together (gpu.dock9bStamNum {}) "
                            "use at most 90% of it and half of the free device memory",
                            hostLimit.bytes >> 20, hostLimit.source,
                            opt.gpuDock9bStamNum > 0 ? std::to_string(opt.gpuDock9bStamNum) : std::string("auto")));
  // Largest host (resident) and device memory of the field of a StaM conformation, over all rounds;
  // the host memory of the first field also holds the reference grids, so it counts only while no
  // other field has been built.
  std::size_t hostPerField = 0, gpuPerField = 0;
  long long fieldsBuilt = 0;
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
  // gpu.packedMaps (spec 15.1.32); auto: without delta blocks, when they fit in a quarter of the free device memory.
  const bool packedMaps = opt.gpuPackedMaps == "true" || (opt.gpuPackedMaps == "auto" && gridMode != od::GridMode::Delta);
#ifdef OD_HAVE_CUDA
  for (const auto& f : gpuFields)
    f->setPackedMaps(!packedMaps ? od::GpuField::Packed::Off
                     : opt.gpuPackedMaps == "true" ? od::GpuField::Packed::On
                                                   : od::GpuField::Packed::IfFits);
#endif
  // Grids built once and kept across rounds: full without StaM torsions (input conformation),
  // flexres (fixed atoms of the input conformation), delta (reference: input conformation).
  std::shared_ptr<const od::EnergyGrids> fixedGrids;
  std::shared_ptr<const od::FftStam> fixedFft;  // posm2_fft, StaM without torsions
  std::shared_ptr<const od::DistanceGrid> fixedDist;
  od::GridAtoms deltaRefAtoms;  // delta: StaM grid atoms of the reference
  // delta on the GPU: host copies of the blocks only for the CPU users of the field (posm3/posm4
  // regions, posm2_fft spectra; spec 15.1.19).
#ifdef OD_HAVE_CUDA
  const bool deltaHostData = !gpuGrid || !poseEval[0] || usesSepas || usesFft;
#endif

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
    if (!refs.empty()) setReference(ms[i], *mobSets[i], refs[i], opt);
    od::logInfo(std::format("  {:<40} topSampleNum {}  sampleBufferNum {}  9d symmetry: {} ({:.3f} s)",
                            mobSets[i]->name(), ms[i].topNum, ms[i].bufferNum,
                            ms[i].sym->usesPermutations()
                                ? std::format("{} automorphism(s)", ms[i].sym->permutationCount())
                                : std::format("more than {} automorphisms, one-to-one pairing within each AD4 type",
                                              opt.symMaxPerm),
                            symSec));
  }

  float stamRef = std::numeric_limits<float>::infinity();  // fixed reference of RoundSample::stamRel (spec 9)
  bool limitStop = false;  // spec 9, step 6: a dcss at the host memory limit did not change; no further rounds
  std::string limitText;
  std::uint64_t memBase = 0;  // spec 9, step 6: largest host memory of the rounds without the dcss-dependent parts
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
    limitStop = r.get<std::uint8_t>() != 0;
    memBase = r.get<std::uint64_t>();
    for (std::size_t i = 0; i < ms.size(); ++i)
      loadMobState(r, ms[i], stamSet.sampleCount(), mobSets[i]->sampleCount(), mobSets[i]->sys().torsionCount());
    if (!r.atEnd()) od::fail("{}: unexpected data at the end", path);
    firstRound = ck.resumeRound + 1;
    od::logInfo(std::format("Resumed from {} (round {} completed)", path, ck.resumeRound));
    if (limitStop) {
      limitText = std::format("the rounds ended after round {} because a dcss at the host memory limit did not change (spec 9, "
                              "step 6); a smaller dock9c.topSampleNum (now {}) in a new run lets the dcss change",
                              ck.resumeRound, ms.empty() ? 0 : ms[0].topNum);
      od::logWarn("  " + limitText);
      firstRound = opt.maxRound + 1;
    } else if (firstRound > opt.maxRound) {
      od::logInfo(std::format("  dock9c.maxRound {}: no further rounds", opt.maxRound));
    }
  }
  const bool peakReset = od::resetPeakResident();
  std::size_t roundPeak = 0;
  // Memory after a stage (spec 9, step 6): host resident and peak resident since the last stage, used device memory.
  const auto stageMemory = [&](const std::string& stage) {
    const std::size_t res = od::residentBytes(), peak = std::max(od::peakResidentBytes(), res);
    roundPeak = std::max(roundPeak, peak);
    od::resetPeakResident();
    std::string dev;
#ifdef OD_HAVE_CUDA
    if (gpu)
      for (int d : opt.gpuDevices) {
        const auto info = od::gpuDeviceInfo(d);
        dev += std::format(", GPU {} in use {} of {} MiB", d, (info.totalMem - info.freeMem) >> 20, info.totalMem >> 20);
      }
#endif
    od::logInfo(std::format("  {} memory: host resident {} MiB, peak {} MiB{}{}", stage, res >> 20, peak >> 20,
                            peakReset ? "" : " (peak since the start)", dev));
  };
  const auto sizeText = [](const ConfSet& set) {
    const od::ConfSampleNums& n = set.sampleNums();
    return std::format("{} (minL2confSampleNum {}, eliteConfSampleNum {}, L2 pool {})", n.docking, n.minL2, n.elite, set.poolSize());
  };
  for (int round = firstRound; round <= opt.maxRound; ++round) {
    od::logInfo(std::format("Round {} of {}:", round, opt.maxRound));
    const std::size_t roundBase = od::residentBytes();
    od::resetPeakResident();
    roundPeak = roundBase;
    {
      std::string t = stamSet.movable() ? "StaM " + sizeText(stamSet) : "StaM 1 (input conformation)";
      for (const auto& m : mobSets) t += ", " + m->name() + " " + (m->movable() ? sizeText(*m) : "1 (input conformation)");
      od::logInfo("  dcss size: " + t);
    }
    // 9a: refill the dcss (spec 7.3).
    const auto a0 = std::chrono::steady_clock::now();
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
    od::logInfo(std::format("  9a total {:.2f} s", std::chrono::duration<double>(std::chrono::steady_clock::now() - a0).count()));
    stageMemory("9a");
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
    // Prepared pairs (spec 9, step 6): their host memory per MobM and the largest of one StaM conformation.
    std::vector<double> preparedSum(ms.size(), 0);
    std::vector<long long> preparedCnt(ms.size(), 0);
    std::size_t itemMaxBytes = 0;
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
    double gridSec = 0, fftSec = 0, fftPassSec = 0;
    std::size_t fftPasses = 0, fftPassMax = 0;
    bool gridsReused = false;
    const auto d0 = std::chrono::steady_clock::now();
    long long movedSum = 0, blockSum = 0, blockTotal = 0;
    long long sepaShell = 0, sepaSelected = 0, sepaRegions = 0, sepaPos = 0;
    double sepaCutSum = 0;
    // 9b of a window of StaM conformations (spec 15.1.16): per MobM, the pairs of all of them
    // docked together; sinks run in StaM conformation and pair order.
    // Pairs of one StaM conformation and one MobM, prepared (L1poss, 8.5) when the conformation
    // joins the window.
    struct Group {
      int sid = 0;
      int variant = 0;  // protonation variant of the MobM samples (spec 3.4.1)
      float stamRel = 0;
      std::unique_ptr<od::MobmDocker> dk;
      std::vector<int> mids;
      std::vector<std::vector<double>> mobAngles;
      std::vector<char> posm1Only;
      std::vector<float> protG;  // protonation free energy of each MobM sample
      std::vector<float> lowest;
      std::vector<od::MobmDocker::PairJob> jobs;
      // Window docked on another thread (windowDevices > 1): the chunks of each job, delivered to
      // its sink (sinks) on this thread afterwards.
      std::vector<std::function<void(std::vector<od::DockSample>&)>> sinks;
      std::vector<std::vector<std::vector<od::DockSample>>> chunks;
    };
    // The members referenced by the dockers stay in place when the item moves.
    struct WindowItem {
      int sid;
      std::unique_ptr<od::StamField> field;
      std::unique_ptr<std::vector<od::StamSepa>> sepas;
      std::shared_ptr<const od::FftStam> fft;
      std::vector<std::vector<std::unique_ptr<Group>>> groups;  // per MobM, one per protonation variant; empty: not made yet
    };
    std::vector<WindowItem> window;
    int curDev = 0;  // device (index in gpu.device) of the window being built
    // The pairs of the StaM conformation that joined the window (its last item), prepared per MobM
    // (L1poss, 8.5); windowPoses: L2poss of the window per MobM.
    std::vector<long long> windowPoses(ms.size(), 0);
    // gpu.dock9bStamNum = auto: L2poss of a full window per MobM (0: unknown, e.g. without a GPU).
    std::vector<long long> poseTarget(ms.size(), -1);
    // The pair jobs of a StaM conformation per MobM (without the L1poss).
    const auto makeGroups = [&](WindowItem& w) {
      w.groups.resize(ms.size());
      const int sid = w.sid;
      for (std::size_t i = 0; i < ms.size(); ++i) {
        const ConfSet& mob = *mobSets[i];
        if (std::find(ms[i].stamDcss.begin(), ms[i].stamDcss.end(), sid) == ms[i].stamDcss.end()) continue;
        std::map<int, std::unique_ptr<Group>> byVariant;
        for (int mid : ms[i].mobDcss) {
          const std::uint64_t key = pairKey(sid, mid);
          if (deterministic && ms[i].pairLowest.contains(key)) {
            ++skipped[i];
            continue;
          }
          const od::ConfSample& smp = mob.sample(mid);
          auto& gp = byVariant[smp.variant];
          if (!gp) {
            gp = std::make_unique<Group>();
            gp->sid = sid;
            gp->variant = smp.variant;
          }
          Group& g = *gp;
          // posm4: the posm3 part of a pair docked before is not redone (spec 9).
          g.posm1Only.push_back(dp.posm == od::Posm::Posm4 && ms[i].pairLowest.contains(key));
          g.mobAngles.push_back(mob.movable() ? std::vector<double>(smp.angles.begin(), smp.angles.end()) : mob.sys().inputAngles());
          g.protG.push_back(smp.protG);
          g.mids.push_back(mid);
        }
        for (auto& [variant, gp] : byVariant) {
          Group& g = *gp;
          od::Dock9bParams dpi = dp;
          dpi.evaluator = evalFor(curDev, i, variant);
          g.dk = std::make_unique<od::MobmDocker>(mob.sysOf(variant), *w.field, dpi, w.sepas.get(), w.fft.get());
          g.stamRel = stamSet.movable() ? stamSet.sample(sid).energy - stamRef : 0.0f;
          g.lowest.assign(g.mids.size(), std::numeric_limits<float>::infinity());
          g.jobs.resize(g.mids.size());
          for (std::size_t j = 0; j < g.mids.size(); ++j) {
            g.jobs[j].mobAnglesDeg = g.mobAngles[j];
            g.jobs[j].stream = od::dockStream(round, i, sid, g.mids[j]);
            g.jobs[j].posm1Only = g.posm1Only[j];
            g.jobs[j].sink = [&, i, j, gr = &g](std::vector<od::DockSample>& chunk) {
              sampleNum[i] += static_cast<long long>(chunk.size());
              for (auto& s : chunk) {
                s.energy += gr->protG[j];  // l2: protonation free energy of the MobM sample
                od::RoundSample r{gr->sid, gr->mids[j], std::move(s), gr->stamRel};
                gr->lowest[j] = std::min(gr->lowest[j], r.rank());
                buffer[i].push_back(std::move(r));
              }
              if (buffer[i].size() >= ms[i].bufferNum) {
                select(i, false);
                ++partial[i];
              }
            };
          }
          w.groups[i].push_back(std::move(gp));
        }
      }
    };
    const auto prepareItem = [&] {
      std::size_t itemBytes = 0;
      WindowItem& w = window.back();
      gpuFieldSlots(curDev, static_cast<int>(window.size()));
      if (w.groups.empty()) makeGroups(w);
      for (std::size_t i = 0; i < ms.size(); ++i)
        for (const auto& gp : w.groups[i]) {
          Group& g = *gp;
          g.dk->prepareJobs(g.jobs);
          for (const auto& j : g.jobs) {
            windowPoses[i] += static_cast<long long>(j.prepared->l2.size());
            const std::size_t b = od::MobmDocker::preparedBytes(*j.prepared);
            itemBytes += b;
            preparedSum[i] += static_cast<double>(b);
            ++preparedCnt[i];
          }
          od::PoseEvaluator* ev = evalFor(curDev, i, g.variant);
          if (poseTarget[i] < 0) poseTarget[i] = ev ? kWindowPosesPerResident * static_cast<long long>(ev->residentPoses(*g.dk)) : 0;
        }
      itemMaxBytes = std::max(itemMaxBytes, itemBytes);
    };

    std::size_t windowBytes = 0, windows = 0, largestWindow = 0;
    const std::size_t fieldBudget = gpuFieldBudget();
    const std::size_t roundResident = od::residentBytes();
    long long endByPoses = 0, endByCount = 0, endByGpu = 0, endByHost = 0;
    // Device memory of a field (estimate): the stored delta blocks, or the maps of a full field of
    // a StaM with torsions (the other fields share the reference maps).
    const auto fieldBytes = [&](const od::StamField& f) -> std::size_t {
      const std::size_t maps = static_cast<std::size_t>(f.typeMapCount()) + 4;
      if (f.mode() == od::GridMode::Delta)
        return f.deltaBlockCount() * maps * od::kDeltaBlock * od::kDeltaBlock * od::kDeltaBlock * sizeof(float) +
               f.blockCount() * sizeof(int);
      if (f.mode() == od::GridMode::Full && stamSet.movable()) return f.geom().size() * maps * sizeof(float) * (packedMaps ? 5 : 1);
      return 0;
    };
    const auto groupsOf = [&](std::vector<WindowItem>& win, std::size_t i) {
      std::vector<Group*> groups;
      for (WindowItem& w : win)
        for (const auto& g : w.groups[i]) groups.push_back(g.get());
      return groups;
    };
    // The pairs of MobM i in the window; keep: the chunks kept for the sinks (finishWindow).
    const auto dockGroups = [&](std::vector<WindowItem>& win, std::size_t i, bool keep) {
      const std::vector<Group*> groups = groupsOf(win, i);
      std::map<int, std::vector<od::MobmDocker::PairGroup>> pgOf;  // per protonation variant (one evaluator each)
      for (Group* g : groups) {
        pgOf[g->variant].push_back({g->dk.get(), &g->jobs});
        if (!keep) continue;
        g->chunks.assign(g->jobs.size(), {});
        g->sinks.resize(g->jobs.size());
        for (std::size_t j = 0; j < g->jobs.size(); ++j) {
          g->sinks[j] = std::move(g->jobs[j].sink);
          g->jobs[j].sink = [c = &g->chunks[j]](std::vector<od::DockSample>& chunk) { c->push_back(std::move(chunk)); };
        }
      }
      for (auto& [variant, pg] : pgOf) od::MobmDocker::dockPairs(pg, ms[i].bufferNum);
    };
    // Statistics and lowest energies of the docked pairs of MobM i; kept: the kept chunks are
    // delivered first, in StaM conformation, pair and chunk order (as without keeping them).
    const auto finishWindow = [&](std::vector<WindowItem>& win, std::size_t i, bool kept) {
      const std::vector<Group*> groups = groupsOf(win, i);
      if (kept)
        for (Group* g : groups)
          for (std::size_t j = 0; j < g->jobs.size(); ++j) {
            for (auto& chunk : g->chunks[j]) g->sinks[j](chunk);
            g->chunks[j] = {};
          }
      {
        for (const Group* gp : groups) {
          const Group& g = *gp;
          for (std::size_t j = 0; j < g.mids.size(); ++j) {
            const std::uint64_t key = pairKey(g.sid, g.mids[j]);
            st[i] += g.jobs[j].stats;
            ++docked[i];
            if (usesSepas && !g.posm1Only[j]) ++posm3Docked[i];
            ms[i].pairLowest[key] = g.posm1Only[j] ? std::min(g.lowest[j], ms[i].pairLowest[key]) : g.lowest[j];
          }
        }
      }
    };
    // Several devices: a window docked on a host thread (its device's evaluators) while the next
    // window is built and prepared on the first device that is free; the sinks and statistics in
    // window order.
    struct Flight {
      std::vector<WindowItem> items;
      int dev = 0;
      bool done = false;  // flightMutex
      std::exception_ptr error;
      std::thread th;
      ~Flight() {
        if (th.joinable()) th.join();
      }
    };
    std::deque<std::unique_ptr<Flight>> flights;
    std::mutex flightMutex;
    std::condition_variable flightDone;
    const auto landFront = [&] {
      // Joined before it leaves flights: the helpers (spec 15.1.35) see its device busy until then.
      flights.front()->th.join();
      std::unique_ptr<Flight> f;
      {
        std::lock_guard lk(flightMutex);
        f = std::move(flights.front());
        flights.pop_front();
      }
      if (f->error) std::rethrow_exception(f->error);
      for (std::size_t i = 0; i < ms.size(); ++i) finishWindow(f->items, i, true);
    };
    // The docked windows at the front delivered; a free device for the next window.
    const auto nextDevice = [&] {
      for (;;) {
        bool front;
        {
          std::lock_guard lk(flightMutex);
          front = !flights.empty() && flights.front()->done;
        }
        if (!front) break;
        landFront();
      }
      std::unique_lock lk(flightMutex);
      for (;;) {
        for (int d = 0; d < windowDevices; ++d)
          if (std::none_of(flights.begin(), flights.end(), [&](const auto& f) { return f->dev == d && !f->done; })) return d;
        flightDone.wait(lk);
      }
    };
    const auto dockWindow = [&] {
      if (window.empty()) return;
      gpuFieldSlots(curDev, static_cast<int>(window.size()));
      ++windows;
      largestWindow = std::max(largestWindow, window.size());
      if (windowDevices == 1) {
        for (std::size_t i = 0; i < ms.size(); ++i) {
          dockGroups(window, i, false);
          finishWindow(window, i, false);
        }
        return;
      }
      auto f = std::make_unique<Flight>();
      f->items = std::move(window);
      f->dev = curDev;
      window.clear();
      Flight* fp = f.get();
      {
        std::lock_guard lk(flightMutex);
        flights.push_back(std::move(f));
      }
      fp->th = std::thread([&, fp] {
        try {
          for (std::size_t i = 0; i < ms.size(); ++i) dockGroups(fp->items, i, true);
        } catch (...) {
          fp->error = std::current_exception();
        }
        {
          std::lock_guard lk(flightMutex);
          fp->done = true;
        }
        flightDone.notify_all();
      });
    };
    // Grid atoms of a StaM conformation sample, with its protonation variant (spec 3.4.1).
    auto stamAtomsOf = [&](int sid, bool input) {
      const auto& sa = stamSet.sample(sid).angles;
      const int v = input ? 0 : stamSet.sample(sid).variant;
      const od::EnergySystem& sys = stamSet.sysOf(v);
      std::vector<double> angles = !input && stamSet.movable() ? std::vector<double>(sa.begin(), sa.end()) : sys.inputAngles();
      std::vector<od::Vec3> x;
      sys.evaluate(angles, nullptr, &x);
      return od::gridAtoms(sys, stamSet.partsOf(v), x);
    };
    // delta: positions that may change grid values; a protonation variant differs in atoms and
    // charges, so all atom positions of both.
    auto deltaChanged = [&](int sid, const od::GridAtoms& a) {
      if (stamSet.sample(sid).variant == 0) return od::changedPositions(deltaRefAtoms, a);
      std::vector<od::Vec3> v = deltaRefAtoms.xyz;
      v.insert(v.end(), a.xyz.begin(), a.xyz.end());
      return v;
    };
#ifdef OD_HAVE_CUDA
    // delta on the GPU (spec 15.1.20): the blocks of the next StaM conformation are computed while
    // the current one joins the window (and the window is docked).
    // Windows divided among the devices: the blocks only on the device of the window. deltaSid: the
    // StaM conformation submitted per device (-1: none); one submitted for a device whose window was
    // docked meanwhile is dropped when the device builds again (no wait behind the docking).
    std::vector<int> deltaSid(static_cast<std::size_t>(windowDevices), -1);
    const auto builderOf = [&](int d) { return d == 0 ? gpuGrid.get() : gpuGridOthers[d - 1].get(); };
    const auto submitDelta = [&](int sid) {
      const od::GridAtoms a = stamAtomsOf(sid, false);
      const std::vector<od::Vec3> changed = deltaChanged(sid, a);
      const od::GridBuildData gb(geom, a, probe, searchParams(opt, true));
      const od::DistBuildData db(geom, a, distCap);
      deltaSid[curDev] = sid;
      for (auto* b : [&] {
             std::vector<od::GpuGridBuilder*> all{builderOf(curDev)};
             if (windowDevices == 1)
               for (const auto& o : gpuGridOthers) all.push_back(o.get());
             return all;
           }())
        b->submitDelta(gb.view, db.view, od::kGridClamp, opt.gridDeltaThreshold,
                       changed.empty() ? nullptr : changed[0].data(), static_cast<int>(changed.size()));
    };
#endif
    // A StaM conformation with its field, regions and spectra (spec 8.1–8.4), in stamIds order.
    struct Built {
      WindowItem item;
      std::size_t bytes = 0, grown = 0;  // fieldBytes, resident memory added
    };
    const auto buildItem = [&](std::size_t si) {
      const int sid = stamIds[si];
      const auto g0 = std::chrono::steady_clock::now();
      const std::size_t resident0 = od::residentBytes();
      auto stamAtoms = [&](bool input) { return stamAtomsOf(sid, input); };
      auto build = [&](const od::GridAtoms& a) {
        return std::pair{std::make_shared<const od::EnergyGrids>(energyGrids(a)),
                         std::make_shared<const od::DistanceGrid>(distanceGrid(a))};
      };
      const bool firstBuild = !fixedGrids;
      auto fieldPtr = std::make_unique<od::StamField>();
      od::StamField& field = *fieldPtr;
      if (gridMode == od::GridMode::Flexres && stamSet.sample(sid).variant != 0) {
        // A protonation variant: other atoms than the fixed grids (spec 3.4.1).
        const auto [g, d] = build(stamAtoms(false));
        field = od::StamField::full(g, d);
      } else if (gridMode == od::GridMode::Flexres) {
        if (!fixedGrids) std::tie(fixedGrids, fixedDist) = build(od::subset(stamAtoms(true), movedMask, false));
        od::GridAtoms moved = od::subset(stamAtoms(false), movedMask, true);
        movedSum += static_cast<long long>(moved.size());
        field = od::StamField::flexres(fixedGrids, fixedDist, std::move(moved), searchParams(opt, true));
      } else if (!stamSet.movable()) {
        if (!fixedGrids) std::tie(fixedGrids, fixedDist) = build(stamAtoms(true));
        field = od::StamField::full(fixedGrids, fixedDist);
      } else if (gridMode == od::GridMode::Delta) {
        if (!fixedGrids) {
          deltaRefAtoms = stamAtoms(true);
          std::tie(fixedGrids, fixedDist) = build(deltaRefAtoms);
        }
#ifdef OD_HAVE_CUDA
        if (gpuGrid) {
          // The reference maps stay on the device; only the differing blocks are copied back.
          if (gpuGrid->referenceId() != 1) {
            std::vector<const float*> maps;
            for (const auto& m : fixedGrids->typeMaps) maps.push_back(m.data());
            maps.push_back(fixedGrids->elecMap.data());
            maps.push_back(fixedGrids->desolvMap.data());
            maps.push_back(fixedDist->d.data());
            gpuGrid->setReference(1, maps.data(), static_cast<int>(maps.size()), geom.size());
            for (const auto& o : gpuGridOthers) o->setReference(1, maps.data(), static_cast<int>(maps.size()), geom.size());
          }
          if (deltaSid[curDev] >= 0 && deltaSid[curDev] != sid) {
            std::vector<int> bk;
            builderOf(curDev)->finishDelta(bk, nullptr);
            deltaSid[curDev] = -1;
          }
          if (deltaSid[curDev] < 0) submitDelta(sid);
          std::vector<int> block;
          std::vector<float> data;
          deltaSid[curDev] = -1;
          auto dev = builderOf(curDev)->finishDelta(block, deltaHostData ? &data : nullptr);
          // The same blocks on the other devices (otherwise GpuField copies them from the first).
          if (windowDevices == 1)
            for (const auto& o : gpuGridOthers) {
              std::vector<int> bk;
              auto other = o->finishDelta(bk, nullptr);
              if (dev && other && bk == block) std::const_pointer_cast<od::DeviceBlocks>(dev)->others.push_back(std::move(other));
            }
          if (si + 1 < stamIds.size()) submitDelta(stamIds[si + 1]);
          field = od::StamField::delta(fixedGrids, fixedDist, std::move(block), std::move(data), std::move(dev));
        } else
#endif
        {
          const od::GridAtoms a = stamAtoms(false);
          const std::vector<od::Vec3> changed = deltaChanged(sid, a);
          const auto [g, d] = build(a);
          field = od::StamField::delta(fixedGrids, fixedDist, *g, *d, opt.gridDeltaThreshold, &changed);
        }
        blockSum += static_cast<long long>(field.deltaBlockCount());
        blockTotal += static_cast<long long>(field.blockCount());
      } else {
#ifdef OD_HAVE_CUDA
        if (gpuGrid && !deltaHostData) {
          // Built on the device of the window and kept there (spec 15.1.34).
          const od::GridAtoms a = stamAtoms(false);
          const od::GridBuildData gb(geom, a, probe, searchParams(opt, true));
          auto g = std::make_shared<od::EnergyGrids>();
          g->geom = geom;
          g->types = gb.types;
          auto d = std::make_shared<od::DistanceGrid>();
          d->geom = geom;
          d->cap = distCap;
          auto dev = builderOf(curDev)->buildDevice(gb.view, od::DistBuildData(geom, a, distCap).view, od::kGridClamp);
          field = od::StamField::full(std::move(g), std::move(d), std::move(dev));
        } else
#endif
        {
          const auto [g, d] = build(stamAtoms(false));
          field = od::StamField::full(g, d);
        }
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
          fft = std::make_shared<const od::FftStam>(od::buildFftStam(field, opt.fftClamp, static_cast<int>(opt.fftStride)));
          if (!stamSet.movable() && gridMode != od::GridMode::Flexres) fixedFft = fft;
        }
        fftSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - f0).count();
      }
      gridSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - g0).count();
      const std::size_t resident = od::residentBytes();
      Built b{{sid, std::move(fieldPtr), std::make_unique<std::vector<od::StamSepa>>(std::move(sepas)), std::move(fft), {}},
              fieldBytes(field), resident > resident0 ? resident - resident0 : 0};
      hostPerField = ++fieldsBuilt == 2 ? b.grown : std::max(hostPerField, b.grown);
      return b;
    };
    // posm2_fft (spec 15.1.22): the translations of the pairs of the built StaM conformations,
    // each MobM conformation against all of them at once.
    const auto fftShare = [&](std::deque<Built>& items) {
      const auto f0 = std::chrono::steady_clock::now();
      for (Built& b : items)
        if (b.item.groups.empty()) makeGroups(b.item);
      for (std::size_t i = 0; i < ms.size(); ++i) {
        std::map<int, std::vector<std::pair<Group*, std::size_t>>> byMid;
        for (Built& b : items)
          for (const auto& gp : b.item.groups[i])
            for (std::size_t j = 0; j < gp->mids.size(); ++j) byMid[gp->mids[j]].emplace_back(gp.get(), j);
        for (auto& [mid, list] : byMid) {
          std::vector<const od::MobmDocker*> dks;
          for (const auto& [g, j] : list) dks.push_back(g->dk.get());
          auto confs = od::MobmDocker::conformations(dks, list[0].first->mobAngles[list[0].second]);
          for (std::size_t k = 0; k < list.size(); ++k) list[k].first->jobs[list[k].second].conf = std::move(confs[k]);
        }
        if (poseEval[i]) poseEval[i]->prepared();
        for (auto& [key, ev] : variantEval)
          if (key.i == i && key.dev == 0) ev->prepared();
      }
      ++fftPasses;
      fftPassMax = std::max(fftPassMax, items.size());
      fftPassSec += std::chrono::duration<double>(std::chrono::steady_clock::now() - f0).count();
    };
    std::deque<Built> ahead;
    for (std::size_t si = 0; si < stamIds.size(); ++si) {
      if (window.empty() && windowDevices > 1) curDev = nextDevice();
      if (ahead.empty()) {
        ahead.push_back(buildItem(si));
        if (usesFft) {
          // Built ahead while the spectra take at most a quarter of the free device memory and the
          // fields fit in the host memory budget.
          std::size_t specBytes = 0;
          for (const auto& s : ahead.front().item.fft->spectra) specBytes += s.size() * sizeof(s[0]);
          while (si + ahead.size() < stamIds.size() && ahead.size() < kFftStamMax &&
                 (ahead.size() + 1) * specBytes <= gpuFree() / 4 &&
                 od::residentBytes() + std::max(hostPerField, ahead.back().grown) <= hostBudget)
            ahead.push_back(buildItem(si + ahead.size()));
          fftShare(ahead);
        }
      }
      Built built = std::move(ahead.front());
      ahead.pop_front();
      // Docked together with the following StaM conformations (spec 15.1.16): the window is docked
      // when it is full or when the field of one more conformation could exceed the memory limits.
      const std::size_t bytes = built.bytes, resident = od::residentBytes();
      gpuPerField = std::max(gpuPerField, bytes);
      windowBytes += bytes;
      window.push_back(std::move(built.item));
      prepareItem();
      // auto: full when every MobM with L2poss in the window has poseTarget of them.
      bool byPoses = opt.gpuDock9bStamNum == 0, anyPoses = false;
      for (std::size_t i = 0; i < ms.size() && byPoses; ++i) {
        if (windowPoses[i] == 0) continue;
        anyPoses = true;
        byPoses = poseTarget[i] > 0 && windowPoses[i] >= poseTarget[i];
      }
      byPoses = byPoses && anyPoses;
      const bool byCount = !byPoses && window.size() >= static_cast<std::size_t>(stamWindow);
      const bool byGpu = !byPoses && !byCount && windowBytes + gpuPerField > fieldBudget;
      const bool byHost = !byPoses && !byCount && !byGpu && resident + hostPerField > hostBudget;
      if (byPoses || byCount || byGpu || byHost) {
        endByPoses += byPoses;
        endByCount += byCount;
        endByGpu += byGpu;
        endByHost += byHost;
        dockWindow();
        window.clear();
        windowBytes = 0;
        std::fill(windowPoses.begin(), windowPoses.end(), 0);
        od::trimHeap();
      }
    }
    dockWindow();
    window.clear();
#ifdef OD_HAVE_CUDA
    // The last windows of the round: devices without a window take over minimizations of the
    // others (spec 15.1.35).
    std::atomic<bool> helpStop{false};
    std::vector<std::thread> helpers;
    std::vector<std::exception_ptr> helpErr(static_cast<std::size_t>(windowDevices));
    if (stealHub && windowDevices > 1 && !flights.empty())
      for (int d = 0; d < windowDevices; ++d)
        helpers.emplace_back([&, d] {
          try {
            const auto evOf = [&](int key) { return gpuEvals[static_cast<std::size_t>(d)][static_cast<std::size_t>(key)]; };
            while (!helpStop.load()) {
              bool busy;
              {
                std::lock_guard lk(flightMutex);
                busy = std::any_of(flights.begin(), flights.end(), [&](const auto& f) { return f->dev == d && !f->done; });
              }
              if (busy || !stealHub->helpOnce(d, evOf)) stealHub->wait(2);
            }
          } catch (...) {
            helpErr[static_cast<std::size_t>(d)] = std::current_exception();
          }
        });
    const auto stopHelpers = [&] {
      helpStop = true;
      for (auto& h : helpers) h.join();
      helpers.clear();
    };
    try {
      while (!flights.empty()) landFront();
    } catch (...) {
      stopHelpers();
      throw;
    }
    stopHelpers();
    for (const auto& e : helpErr)
      if (e) std::rethrow_exception(e);
#else
    while (!flights.empty()) landFront();
#endif
    od::trimHeap();
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
      if (usesFft && fftPasses > 0)
        od::logInfo(std::format("  9b posm2_fft translations: {:.2f} s in {} pass(es), up to {} StaM conformation(s) each",
                                fftPassSec, fftPasses, fftPassMax));
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
    std::vector<int> rmStamOf(ms.size(), -1), rmMobOf(ms.size(), -1);
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
      const int rmStam = stamSet.movable() ? od::filterDcss(m.stamDcss, m.top, &od::RoundSample::stam) : 0;
      const int rmMob = mobSets[i]->movable() ? od::filterDcss(m.mobDcss, m.top, &od::RoundSample::mobm) : 0;
      rmStamOf[i] = rmStam;
      rmMobOf[i] = rmMob;
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
    for (auto& e : poseEval)
      if (e) e->idle();
    for (auto& v : poseEvalOthers)
      for (auto& e : v)
        if (e) e->idle();
    {
      std::string steal;
#ifdef OD_HAVE_CUDA
      if (stealHub) {
        const od::DockStealHub::Stats s = stealHub->stats();
        steal = std::format("; poses taken over by other devices {} in {} part(s), binding {:.2f} s, {:.0f} MiB", s.poses - stealBase.poses,
                            s.steals - stealBase.steals, s.bindSeconds - stealBase.bindSeconds,
                            static_cast<double>(s.bindBytes - stealBase.bindBytes) / (1 << 20));
        stealBase = s;
      }
#endif
      od::logInfo(std::format("  9b docking: {:.2f} s ({} window(s) of StaM conformations docked together, largest {}{})",
                              dockSec, windows, largestWindow, steal));
    }
    if (stamWindow > 1 && windows > 0) {
      // Memory of one window with all StaM conformations of the round (a hint for memory.hostLimit).
      const std::size_t n = stamIds.size();
      std::string target;
      if (opt.gpuDock9bStamNum == 0)
        for (std::size_t i = 0; i < ms.size(); ++i)
          if (poseTarget[i] > 0) target += std::format("{}{}", target.empty() ? "" : ", ", poseTarget[i]);
      od::logInfo(std::format("  9b windows ended by: L2poss {}{}, gpu.dock9bStamNum {}, device memory {}, host memory {}, "
                              "last {}; per StaM conformation: host {} MiB, device {} MiB",
                              endByPoses, target.empty() ? "" : " (" + target + " per window)", endByCount, endByGpu,
                              endByHost, windows - static_cast<std::size_t>(endByPoses + endByCount + endByGpu + endByHost),
                              hostPerField >> 20, gpuPerField >> 20));
      if (endByGpu + endByHost > 0)
        od::logInfo(std::format("  9b one window for all {} StaM conformations would need memory.hostLimit >= about {} MiB "
                                "and about {} MiB of free device memory",
                                n, (roundResident + n * hostPerField) / 9 * 10 >> 20, 2 * n * gpuPerField >> 20));
    }
    stageMemory("9b/9c");

    // dcss size (spec 9, step 6): doubled when no dcss sample was removed, within the host memory.
    {
      std::vector<int> sizes{stamSet.dockingNum()};
      std::vector<bool> grow{false};
      for (std::size_t i = 0; i < ms.size(); ++i) {
        sizes.push_back(mobSets[i]->dockingNum());
        grow.push_back(mobSets[i]->movable() && rmMobOf[i] == 0);
        if (stamSet.movable() && rmStamOf[i] == 0) grow[0] = true;
      }
      for (std::size_t i = 0; i < ms.size(); ++i) {
        if (preparedCnt[i] > 0) ms[i].preparedPerPair = preparedSum[i] / static_cast<double>(preparedCnt[i]);
        if (docked[i] > 0) ms[i].samplesPerPair = static_cast<double>(sampleNum[i]) / static_cast<double>(docked[i]);
      }
      {
        // Host memory of a round (spec 9, step 6): the measured peak of this round, with the parts that depend on the
        // dcss sizes replaced by their values at the new sizes (all pairs docked).
        enum { kPrepared, kSamples, kPairs, kPool, kBatch, kTerms };
        const auto terms = [&](const std::vector<int>& s, bool actual) {
          std::array<double, kTerms> t{};
          for (std::size_t i = 0; i < ms.size(); ++i) {
            const double pairs = static_cast<double>(s[0]) * s[i + 1];
            const double sampleBytes =
                sizeof(od::RoundSample) + 4.0 * mobSets[i]->sys().torsionCount() + 16.0;
            if (!actual) t[kPrepared] += s[i + 1] * ms[i].preparedPerPair;
            const double n = actual ? static_cast<double>(sampleNum[i]) : ms[i].samplesPerPair * pairs;
            t[kSamples] += 2.0 * std::min(static_cast<double>(ms[i].bufferNum), n) * sampleBytes;
            t[kPairs] += 64.0 * (static_cast<double>(ms[i].pairLowest.size()) + (actual ? 0.0 : pairs));
          }
          if (actual) t[kPrepared] = static_cast<double>(itemMaxBytes);
          const auto group = [&](const ConfSet& set, od::SampleSet which, int size, const od::CsgmOpts& c) {
            if (!set.movable()) return;
            const int minL2 = size == set.dockingNum() ? set.sampleNums().minL2 : od::scaledConfSampleNums(opt, which, size).minL2;
            const double tors = set.sys().torsionCount();
            t[kPool] += (minL2 + size) * (sizeof(od::ConfSample) + 4.0 * tors + 16.0);
            t[kBatch] += static_cast<double>(std::min<long long>(c.maxL1BatchNum, 2LL * minL2)) * (32.0 * tors + 64.0);
          };
          group(stamSet, od::SampleSet::StaM, s[0], opt.stamCsgm);
          for (std::size_t i = 0; i < ms.size(); ++i) group(*mobSets[i], od::SampleSet::MobM, s[i + 1], opt.mobmCsgm);
          return t;
        };
        const auto sum = [](const std::array<double, kTerms>& t) { return std::accumulate(t.begin(), t.end(), 0.0); };
        const double now = sum(terms(sizes, true));
        memBase = std::max<std::uint64_t>({memBase, roundBase, roundPeak > now ? roundPeak - static_cast<std::size_t>(now) : 0});
        const auto base = static_cast<std::size_t>(memBase);
        if (std::find(grow.begin(), grow.end(), true) != grow.end()) {
          const auto estimate = [&](const std::vector<int>& s) { return base + static_cast<std::size_t>(sum(terms(s, false))); };
          const od::DcssGrowth g = od::growDcssSizes(sizes, grow, estimate, hostBudget);
          const auto next = terms(g.sizes, false);
          const auto mib = [](double b) { return static_cast<long long>(b / (1 << 20)); };
          const auto nameOf = [&](std::size_t k) { return k == 0 ? std::string("StaM") : mobSets[k - 1]->name(); };
          std::string which, change;
          for (std::size_t k = 0; k < sizes.size(); ++k)
            if (grow[k]) {
              which += (which.empty() ? "" : ", ") + nameOf(k);
              change += std::format("{}{} {} -> {}", change.empty() ? "" : ", ", nameOf(k), sizes[k], g.atLimit ? sizes[k] : g.sizes[k]);
            }
          od::logInfo(std::format("  dcss: no sample removed ({}); host memory: round peak {} MiB (of which dcss-dependent {} MiB); "
                                  "estimate with the new sizes {} MiB = {} MiB + 9b prepared pairs of one StaM conformation {} MiB "
                                  "+ docking samples {} MiB + pair records {} MiB + L2 pool and dcss {} MiB + 9a batch {} MiB; budget {}",
                                  which, roundPeak >> 20, mib(now), g.estimate >> 20, base >> 20, mib(next[kPrepared]),
                                  mib(next[kSamples]), mib(next[kPairs]), mib(next[kPool]), mib(next[kBatch]),
                                  hostBudget > 0 ? std::format("{} MiB (90% of memory.hostLimit)", hostBudget >> 20)
                                                 : std::string("unknown (no limit)")));
          if (g.atLimit) {
            limitStop = true;
            limitText = std::format("round {}: no dcss sample of {} was removed and the dcss size cannot grow within the host memory "
                                    "budget{}. A smaller dock9c.topSampleNum (now {}) in a new run "
                                    "lets the dcss change",
                                    round, which, round < opt.maxRound ? "; the remaining rounds are skipped" : "",
                                    ms.empty() ? 0 : ms[0].topNum);
            od::logWarn("  " + limitText);
          } else if (round == opt.maxRound) {
            od::logInfo(std::format("  dcss size: {} not applied (last round)", change));
          } else {
            if (grow[0] && g.sizes[0] > sizes[0]) stamSet.setSampleNums(od::scaledConfSampleNums(opt, od::SampleSet::StaM, g.sizes[0]));
            for (std::size_t i = 0; i < ms.size(); ++i)
              if (grow[i + 1] && g.sizes[i + 1] > sizes[i + 1])
                mobSets[i]->setSampleNums(od::scaledConfSampleNums(opt, od::SampleSet::MobM, g.sizes[i + 1]));
            od::logInfo(std::format("  dcss size: {}", change));
          }
        }
      }
    }
    writeCheckpoint(round, stamSet, mobSets, ms, stamRef, limitStop, memBase, opt, ck, shift);
    stageMemory("checkpoint");
    if (limitStop) break;
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
  run9d(stamSet, mobSets, ms, opt, fp, shift, gpu);
  stageMemory("9d");
  if (limitStop) od::logWarn(limitText);
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
// (input files may be given again, e.g. when moved); only dock9c.maxRound, gpu.*, cpu.threadNum, memory.hostLimit and output.allSummary may differ.
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
  cmp.gpuEnable = stored.gpuEnable;
  cmp.gpuDevices = stored.gpuDevices;
  cmp.gpuConf9aBatchNum = stored.gpuConf9aBatchNum;
  cmp.gpuConf9aPrecision = stored.gpuConf9aPrecision;
  cmp.gpuConf9aRcd = stored.gpuConf9aRcd;
  cmp.gpuConf9aRcdMinBatch = stored.gpuConf9aRcdMinBatch;
  cmp.gpuDock9bBatchNum = stored.gpuDock9bBatchNum;
  cmp.gpuDock9bStamNum = stored.gpuDock9bStamNum;
  cmp.gpuDock9bDeviceMin = stored.gpuDock9bDeviceMin;
  cmp.gpuPackedMaps = stored.gpuPackedMaps;
  cmp.memHostLimit = stored.memHostLimit;
  cmp.cpuThreadNum = stored.cpuThreadNum;
  cmp.outAllSummary = stored.outAllSummary;
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
    od::fail("--resume: only dock9c.maxRound, gpu.*, cpu.threadNum, memory.hostLimit and output.allSummary may differ from the checkpoint parameters:{}", diff);
  }
  ck.resumeDir = dir;
  ck.resumeRound = round;
  ck.stateSha256 = man.stateSha256;
  ck.inputs = man.inputs;
  return opt;
}

std::string vecText(const od::Vec3& v) { return std::format("{:.3f}, {:.3f}, {:.3f}", v[0], v[1], v[2]); }

// Settings of the MobMs for box.mobmRadius = sampled / max (spec 3.3): from --molConfig when it is
// valid for them, otherwise from their default configuration; nullopt for a MobM without either.
std::vector<std::optional<od::MolSettings>> boxMobmSettings(const od::Options& opt, const std::vector<od::Molecule>& mobs,
                                                            od::MolConfig& cfgOut) {
  std::vector<od::Topology> topo;
  for (const auto& m : mobs) {
    topo.push_back(od::analyzeTopology(m, {}));
    addMacroBreaks(m, topo.back(), {}, opt.ringMacroMinSize);
    od::assignDssp({&m}, {&topo.back()});
  }
  std::vector<od::MolInput> in;
  for (std::size_t i = 0; i < mobs.size(); ++i) in.push_back({&mobs[i], &topo[i]});
  std::vector<std::optional<od::MolSettings>> out(mobs.size());
  std::vector<std::string> errors;
  if (!opt.molConfig.empty())
    if (auto cfg = od::readMolConfig(opt.molConfig, in, {}, errors, true)) {
      for (std::size_t i = 0; i < mobs.size(); ++i) out[i] = std::move(cfg->molecules[i]);
      cfg->molecules.clear();
      cfgOut = std::move(*cfg);
      return out;
    }
  bool haveCfg = false;
  for (std::size_t i = 0; i < mobs.size(); ++i) {
    errors.clear();
    if (auto cfg = od::parseMolConfig(od::defaultMolConfigText({in[i]}), {in[i]}, {}, errors)) {
      out[i] = std::move(cfg->molecules[0]);
      if (!haveCfg) {
        cfg->molecules.clear();
        cfgOut = std::move(*cfg);
        haveCfg = true;
      }
    }
  }
  return out;
}

// Inputs of generateRingConformers from the topology: the broken bonds and the central bonds of the
// torsions in the flexible range.
std::vector<std::array<int, 2>> ringKey(const od::Topology& t) {
  std::vector<std::array<int, 2>> k;
  for (int b : t.brokenBonds) k.push_back({-1, b});
  for (const auto& tor : t.torsions)
    if (tor.inFlexRange) k.push_back({std::min(tor.atoms[1], tor.atoms[2]), std::max(tor.atoms[1], tor.atoms[2])});
  std::sort(k.begin(), k.end());
  return k;
}

// Ring conformers of a MobM made for the box (input coordinates), reused by 9a when the key matches.
struct BoxRings {
  std::vector<std::array<int, 2>> key;
  od::RingConformers rc;
};

od::RingParams ringParams(const od::Options& opt, std::uint64_t seed) {
  return {opt.ringMacroMinSize, opt.ringConfGenNum, opt.ringConfRmsd, opt.ringConfEnergyWindow, opt.ringMaxConfNum, seed};
}

// Docking box of the run (spec 3.3). Automatic modes: from the StaM heavy atoms and the MobM radii.
// boxRings: per MobM, the ring conformers used for box.mobmRadius = max / sampled (spec 4.9).
od::BoxFrame resolveBox(const od::Options& opt, const std::vector<od::Molecule>& sta, const std::vector<od::Molecule>& mobs,
                        std::uint64_t seed, std::vector<std::optional<BoxRings>>& boxRings) {
  boxRings.assign(mobs.size(), std::nullopt);
  od::BoxFrame b;
  if (opt.boxMode == "manual") {
    b.center = *opt.boxCenter;
    b.size = *opt.boxSize;
    od::logInfo(std::format("Docking box (box.mode = manual): center {} Å (input coordinates), size {} Å (along the input axes)",
                            vecText(b.center), vecText(b.size)));
    return b;
  }
  if (opt.boxCenter)
    od::logWarn(std::format("box.center is given with box.mode = {}: it replaces the computed center", opt.boxMode));
  if (opt.boxSize)
    od::logWarn(std::format("box.size is given with box.mode = {}: it replaces the computed size", opt.boxMode));
  const auto t0 = std::chrono::steady_clock::now();
  double rmax = 0;
  std::string rInfo;
  {
    od::MolConfig cfg;
    std::vector<std::optional<od::MolSettings>> settings;
    if (opt.boxMobmRadius != "input") settings = boxMobmSettings(opt, mobs, cfg);
    for (std::size_t i = 0; i < mobs.size(); ++i) {
      const od::Molecule& m = mobs[i];
      const std::vector<od::Vec3> heavy = od::heavyCoords(m);
      double r = od::radiusOf(heavy);
      std::string how = "input conformation";
      // Ring conformers (spec 4.9) as 9a will use them: same parameters and random stream.
      const od::RingConformers* rc = nullptr;
      std::string ringHow;
      if (opt.boxMobmRadius != "input" && settings[i] && settings[i]->ringConformers &&
          (opt.runMode == "dock" || opt.runUnbound == "conf9a") && opt.ringMaxConfNum > 1) {
        boxRings[i] = BoxRings{ringKey(settings[i]->topology),
                               od::generateRingConformers(m, settings[i]->topology, ringParams(opt, seed),
                                                          (1ULL << 60) + sta.size() + i)};
        rc = &boxRings[i]->rc;
        int vars = 0, confs = 0;
        for (const auto& s : rc->systems)
          if (s.conf.size() > 1) {
            ++vars;
            confs += static_cast<int>(s.conf.size());
          }
        if (vars > 0) ringHow = std::format(", {} ring system(s) with {} conformers", vars, confs);
      }
      if (opt.boxMobmRadius != "input" && !settings[i]) {
        od::logWarn(std::format("{}: no valid molecule configuration for box.mobmRadius = {}; the radius of the input "
                                "conformation is used",
                                m.id, opt.boxMobmRadius));
      } else if (opt.boxMobmRadius == "max") {
        const od::RingOps ops = rc ? od::buildRingOps(m, settings[i]->tree, *rc, {}) : od::RingOps{};
        int overflow = 0;
        r = od::mobmRadiusMax(m, settings[i]->tree, &ops, &overflow);
        how = std::format("upper bound, {} rotatable torsions{}", settings[i]->tree.torsions.size(), ringHow);
        if (overflow > 0)
          od::logWarn(std::format("{}: box.mobmRadius = max: {} rigid part(s) with too many ring conformer combinations; "
                                  "each ring system was varied alone there, the radius may not be an upper bound",
                                  m.id, overflow));
      } else if (opt.boxMobmRadius == "sampled") {
        const od::EnergySystem sys({{&m, &settings[i]->tree, rc, nullptr}}, energyParams(opt));
        od::Csgm1GroupSampler smp;
        smp.add(od::Csgm1Sampler(m, *settings[i], cfg, {opt.mobmCsgm.csgm1AngleStep, opt.mobmCsgm.ramaFloor}));
        // Breaks are closed as in 9a (spec 7.11) with the MobM closure settings.
        od::ClosureParams clp = closureParams(opt.mobmCsgm, opt);
        clp.growClose = false;
        const std::vector<bool> cis = smp.cisTransMask();
        std::vector<bool> freeMask(cis.size());
        for (std::size_t g = 0; g < freeMask.size(); ++g) freeMask[g] = !cis[g];
        const od::LoopCloser closer(sys, {{&m, &settings[i]->tree}}, freeMask, smp.segments(), clp);
        const double maxDev = opt.confPreCloseFactor > 0 && sys.hasClosures()
                                  ? opt.confPreCloseFactor * opt.maxBreakDist
                                  : 0.0;
        // Random streams of the box samples, apart from those of 9a.
        const od::SampledRadius s = od::mobmRadiusSampled(m, sys, smp, opt.boxMobmSampleNum, seed, (1ULL << 62) + i,
                                                          &closer, maxDev);
        r = s.radius;
        how = sys.torsionCount() == 0 ? std::string("no rotatable torsions")
                                      : std::format("input and {} of {} csgm1 samples without clashes{}", s.accepted,
                                                    opt.boxMobmSampleNum,
                                                    closer.active() ? std::format(", closed ({} not closed)", s.notClosed)
                                                                    : std::string()) +
                  ringHow;
      }
      od::logInfo(std::format("  box: MobM {:<40} radius {:.3f} Å ({})", m.id, r, how));
      rmax = std::max(rmax, r);
    }
  }
  std::vector<od::Vec3> heavy;
  for (const auto& m : sta) {
    const auto h = od::heavyCoords(m);
    heavy.insert(heavy.end(), h.begin(), h.end());
  }
  od::AutoBoxParams ap;
  ap.stamMargin = opt.boxStamMargin;
  ap.mobmRadius = rmax;
  ap.autoRotate = opt.boxAutoRotate;
  ap.center = opt.boxCenter;
  ap.size = opt.boxSize;
  b = od::autoBox(heavy, ap);
  const char* axes = b.rotated ? "box axes" : "input axes";
  od::logInfo(std::format("Docking box (box.mode = {}, box.autoRotate = {}; r_max {:.3f} Å, box.stamMargin {} Å): "
                          "center {} Å (input coordinates, {}), size {} Å (along the {}, {})  {:.2f} s",
                          opt.boxMode, opt.boxAutoRotate ? "true" : "false", rmax, opt.boxStamMargin,
                          vecText(b.toInput(b.center)), opt.boxCenter ? "given" : "computed", vecText(b.size), axes,
                          opt.boxSize ? "given" : "computed",
                          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
  if (b.rotated)
    od::logInfo(std::format("  box axes in input coordinates: x' ({}), y' ({}), z' ({}); the molecules are turned to "
                            "them internally, output coordinates are in the input frame",
                            vecText(b.axes[0]), vecText(b.axes[1]), vecText(b.axes[2])));
  return b;
}

// score / localMin (spec 6.1): the poses of every MobM file (box frame) on the 9b grid of the
// input StaM conformation, minimized on it (localMin), then 9d.
void runScore(ConfSet& stamSet, std::vector<std::unique_ptr<ConfSet>>& mobSets,
              const std::vector<std::vector<std::vector<od::Vec3>>>& poses, const std::vector<od::Molecule>& refs,
              const od::Options& opt, std::uint64_t seed, bool gpu) {
  const bool localMin = opt.runMode == "localMin";
  const bool fill9a = opt.runUnbound == "conf9a";
  const int stamId = stamSet.addInput(fill9a);
  std::vector<int> mobIds;
  for (auto& m : mobSets) mobIds.push_back(m->addInput(fill9a));
  float stamRel = stamSet.movable() ? stamSet.sample(stamId).energy - stamSet.minEnergy() : 0.0f;
  if (!std::isfinite(stamRel)) {
    od::logWarn("StaM group: no 9a conformation sample; the StaM term is 0");
    stamRel = 0;
  }

  std::vector<od::Ad4Type> probe;
  for (const auto& m : mobSets)
    for (od::Ad4Type t : m->sys().typing(0).type)
      if (t != od::Ad4Type::Merged && std::find(probe.begin(), probe.end(), t) == probe.end()) probe.push_back(t);
  const od::GridGeometry geom = od::boxGrid(opt.box.center, opt.box.size, opt.gridSpacing);
  if (opt.boxMode != "manual") checkGridMemory(geom, probe.size(), opt, false);
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<od::Vec3> x;
  stamSet.sys().evaluate(stamSet.sys().inputAngles(), nullptr, &x);
  const od::GridAtoms ga = od::gridAtoms(stamSet.sys(), stamSet.parts(), x);
  const double distCap = std::max(opt.maxNearDistance, opt.clashDistance) + 2 * opt.gridSpacing;
  const od::StamField field =
      od::StamField::full(std::make_shared<od::EnergyGrids>(od::buildEnergyGrids(geom, ga, probe, searchParams(opt, true))),
                          std::make_shared<od::DistanceGrid>(od::buildDistanceGrid(geom, ga, distCap)));
  od::logInfo(std::format("{}: 9b grid of the input StaM conformation, {} x {} x {} points, spacing {} Å ({:.2f} s)",
                          opt.runMode, geom.n[0], geom.n[1], geom.n[2], geom.spacing,
                          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
  od::Dock9bParams dp = dockParams(opt, seed);
  dp.posm = od::Posm::Posm1;

  std::vector<MobState> ms(mobSets.size());
  for (std::size_t i = 0; i < mobSets.size(); ++i) {
    const ConfSet& mob = *mobSets[i];
    MobState& m = ms[i];
    m.heavy = std::make_unique<od::MobmHeavyAtoms>(mob.sys());
    if (!refs.empty()) {
      m.sym = std::make_unique<od::SymmetryRmsd>(*mob.parts()[0].mol, m.heavy->atoms(), m.heavy->types(),
                                                 static_cast<std::size_t>(opt.symMaxPerm));
      setReference(m, mob, refs[i], opt);
    }
    const od::MobmDocker dk(mob.sys(), field, dp);
    const od::TorsionTree& tree = *mob.parts()[0].tree;
    const od::Topology& topo = mob.settings()[0]->topology;
    const std::size_t nAtoms = mob.parts()[0].mol->atoms.size();
    const std::size_t n = poses[i].size();
    m.top.resize(n);
    m.inputs.resize(n);
    std::vector<int> steps(n, 0);
    const od::RingConformers* rings = mob.parts()[0].rings;
    const od::RingOps& ringOps = mob.sys().ringOps(0);
    const int ringOff = mob.sys().ringSlotOffset(0);
    std::vector<int> ringChanged(n, 0);
    const auto t1 = std::chrono::steady_clock::now();
    od::parallelEach(static_cast<long long>(n), [&](long long r) {
      const std::vector<od::Vec3>& X = poses[i][r];
      std::vector<double> ang;
      for (const auto& tt : tree.torsions) {
        const auto& a = topo.torsions[tt.torsion].atoms;
        ang.push_back(od::dihedralDeg(X[a[0]], X[a[1]], X[a[2]], X[a[3]]));
      }
      ang.resize(static_cast<std::size_t>(mob.sys().torsionCount()), 0.0);
      const std::span<const od::Vec3> to(X.data(), nAtoms);
      std::vector<od::Vec3> from(nAtoms);
      auto fitOf = [&](const std::vector<double>& g) {
        const od::MobmDocker::Conf cg = dk.conformation(g);
        for (std::size_t k = 0; k < nAtoms; ++k) from[k] = {cg.coords[k][0], cg.coords[k][1], cg.coords[k][2]};
        double f = 0;
        od::kabschQuat(from, to, f);
        return f;
      };
      // Ring slots (spec 6.1, 4.9): the conformers with the smallest fit RMSD; every combination when
      // there are at most 1024, otherwise from the conformers with the nearest ring dihedrals, one slot
      // at a time until the fit RMSD no longer drops.
      const std::size_t nSlot = rings ? ringOps.ops.size() : 0;
      if (nSlot > 0) {
        const std::size_t o = static_cast<std::size_t>(ringOff);
        double combos = 1;
        for (const auto& op : ringOps.ops) combos *= op.confCount;
        double best = std::numeric_limits<double>::infinity();
        if (combos <= 1024) {
          std::vector<double> g = ang, keep = ang;
          for (long long idx = 0; idx < static_cast<long long>(combos); ++idx) {
            long long v = idx;
            for (std::size_t k = 0; k < nSlot; ++k) {
              g[o + k] = static_cast<double>(v % ringOps.ops[k].confCount);
              v /= ringOps.ops[k].confCount;
            }
            if (const double f = fitOf(g); f < best - 1e-9) {
              best = f;
              keep = g;
            }
          }
          ang = keep;
        } else {
          for (std::size_t k = 0; k < nSlot; ++k)
            ang[o + k] = od::nearestRingConformer(rings->systems[ringOps.ops[k].system], X);
          best = fitOf(ang);
          for (bool better = true; better;) {
            better = false;
            for (std::size_t k = 0; k < nSlot; ++k) {
              std::vector<double> g = ang;
              for (int c = 0; c < ringOps.ops[k].confCount; ++c) {
                if (c == static_cast<int>(ang[o + k])) continue;
                g[o + k] = c;
                if (const double f = fitOf(g); f < best - 1e-9) {
                  best = f;
                  ang = g;
                  better = true;
                }
              }
            }
          }
        }
        for (std::size_t k = 0; k < nSlot; ++k) ringChanged[r] += ang[o + k] > 0;
      }
      const od::MobmDocker::Conf c = dk.conformation(ang);
      // Pose: the record superposed by the conformation at its torsion angles (bond lengths and angles of record 1).
      for (std::size_t k = 0; k < nAtoms; ++k) from[k] = {c.coords[k][0], c.coords[k][1], c.coords[k][2]};
      double fit = 0;
      const od::Quat q = od::kabschQuat(from, to, fit);
      od::Vec3 cf{}, ct{};
      for (std::size_t k = 0; k < nAtoms; ++k)
        for (int d = 0; d < 3; ++d) {
          cf[d] += from[k][d] / static_cast<double>(nAtoms);
          ct[d] += to[k][d] / static_cast<double>(nAtoms);
        }
      const od::Vec3 rc = od::rotate(q, {c.center[0] - cf[0], c.center[1] - cf[1], c.center[2] - cf[2]});
      od::Pose pose;
      for (int d = 0; d < 3; ++d) pose.t[d] = static_cast<float>(ct[d] + rc[d]);
      pose.q = {static_cast<float>(q.w), static_cast<float>(q.x), static_cast<float>(q.y), static_cast<float>(q.z)};
      const od::MobmDocker::PoseEnergy pe = dk.evaluate(c, pose, ang);
      od::DockSample s;
      s.pose = pose;
      s.center = c.center;
      s.angles.assign(ang.begin(), ang.end());
      s.energy = static_cast<float>(pe.ranking());
      s.inter = static_cast<float>(pe.inter);
      s.interClash = static_cast<float>(pe.interClash);
      s.intra = static_cast<float>(pe.intra.filterEnergy());
      s.closureMaxDeviation = static_cast<float>(pe.intra.closureMaxDeviation);
      std::string status = "ok";
      if (pe.outside > 0) {
        status = "outsideGrid";
        s.energy = s.inter = std::numeric_limits<float>::quiet_NaN();
      } else if (localMin) {
        od::Dock9bStats st;
        od::DockSample out;
        if (dk.minimize(c, pose, out, st)) {
          s = std::move(out);
          steps[r] = s.steps;
        } else {
          status = "leftGrid";
        }
      }
      od::RoundSample& rs = m.top[r];
      rs.stam = stamId;
      rs.mobm = mobIds[i];
      rs.s = std::move(s);
      rs.stamRel = stamRel;
      m.inputs[r] = {static_cast<int>(r) + 1, static_cast<float>(fit), status};
    });
    float fitMax = 0;
    long long outside = 0, left = 0, stepSum = 0;
    for (std::size_t r = 0; r < n; ++r) {
      fitMax = std::max(fitMax, m.inputs[r].fitRmsd);
      outside += m.inputs[r].gridStatus == "outsideGrid";
      left += m.inputs[r].gridStatus == "leftGrid";
      stepSum += steps[r];
    }
    const long long ringPoses = std::count_if(ringChanged.begin(), ringChanged.end(), [](int v) { return v > 0; });
    od::logInfo(std::format("{} {:<38} {} pose(s): largest fit RMSD {:.4f} Å, outside the grid {}{}{} ({:.2f} s)", opt.runMode,
                            mob.name(), n, fitMax, outside,
                            ringOps.ops.empty() ? "" : std::format(", poses with a ring conformer other than the input {}", ringPoses),
                            localMin ? std::format(", left the grid during minimization {}, steps/pose {:.1f}", left,
                                                   n > static_cast<std::size_t>(outside + left)
                                                       ? static_cast<double>(stepSum) / static_cast<double>(n - outside - left)
                                                       : 0.0)
                                     : "",
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count()));
  }
  run9d(stamSet, mobSets, ms, opt, opt.fe, 0.0f, gpu);
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
  {
    const od::CpuThreads ct = od::applyCpuThreads(opt.cpuThreadNum);
    od::logInfo(std::format("CPU threads {} ({}{}), CPUs {}", ct.threads, ct.source, opt.cpuThreadNum ? "" : ", cpu.threadNum auto",
                            od::cpuListText(ct.cpus)));
  }

  od::CcdLibrary ccd;
  for (const auto& f : opt.ccd) ccd.addFile(f);

  // Protonation (spec 3.4).
  const std::vector<od::ProtonationOverride> overrides =
      opt.molConfig.empty() ? std::vector<od::ProtonationOverride>{} : od::readProtonationOverrides(opt.molConfig);
  for (const auto& v : overrides) {
    const bool stam = v.file == opt.stamFile;
    const bool mobm = std::find(opt.mobmFiles.begin(), opt.mobmFiles.end(), v.file) != opt.mobmFiles.end();
    if (!stam && !mobm) od::fail("{}: [[protonation]] file \"{}\" is neither --stamFile nor a --mobmFiles file", v.origin, v.file);
    if ((stam && opt.stamProtonation != "pH") || (!stam && opt.mobmProtonation != "pH"))
      od::fail("{}: [[protonation]] entries need {}.protonation = \"pH\"", v.origin, stam ? "stam" : "mobm");
  }
  auto protOf = [&](const std::string& file, bool stam) -> std::optional<od::ProtonationOptions> {
    if ((stam ? opt.stamProtonation : opt.mobmProtonation) != "pH") return std::nullopt;
    od::ProtonationOptions p;
    p.pH = opt.protonationPH;
    for (const auto& v : overrides)
      if (v.file == file) p.overrides.push_back(v);
    return p;
  };
  if (opt.protonationRecompute != "none" && opt.stamProtonation != "pH" && opt.mobmProtonation != "pH")
    od::logWarn("protonation.recompute: no effect, stam.protonation and mobm.protonation are input");
  std::vector<std::pair<std::string, od::ProtonationReport>> protReports;
  auto load = [&](const std::string& file, od::Role role, bool poseRecords, od::PoseRecords* pr) {
    const auto p = protOf(file, role == od::Role::StaM);
    od::ProtonationReport rep;
    std::vector<od::Molecule> v;
    if (poseRecords) *pr = od::loadPoseRecords(file, role, ccd, p ? &*p : nullptr, &rep);
    else v = od::loadMolecules(file, role, ccd, p ? &*p : nullptr, &rep);
    if (p) {
      logProtonation(file, rep, p->pH);
      protReports.push_back({file, std::move(rep)});
    }
    return v;
  };
  auto sta = load(opt.stamFile, od::Role::StaM, false, nullptr);
  std::vector<od::Molecule> mobs;
  std::vector<std::vector<std::vector<od::Vec3>>> poses;  // score / localMin: coordinates of every record (spec 6.1)
  for (const auto& f : opt.mobmFiles) {
    if (opt.runMode != "dock") {
      od::PoseRecords p;
      load(f, od::Role::MobM, true, &p);
      od::logInfo(std::format("{}: {} pose record(s)", f, p.xyz.size()));
      mobs.push_back(std::move(p.mol));
      poses.push_back(std::move(p.xyz));
      continue;
    }
    auto mob = load(f, od::Role::MobM, false, nullptr);
    if (mob.size() != 1) od::fail("{}: --mobmFiles file must contain exactly one molecule, found {}", f, mob.size());
    mobs.push_back(std::move(mob.front()));
  }
  std::vector<od::Molecule> refs;
  for (const auto& f : opt.refFiles) refs.push_back(od::loadReference(f, ccd));
  // Random seeds fit in a TOML integer (63 bits).
  const std::uint64_t seed =
      opt.seed ? *opt.seed
               : ((static_cast<std::uint64_t>(std::random_device{}()) << 32) ^ std::random_device{}()) >> 1;
  od::logInfo(std::format("seed {}", seed));
  opt.seed = seed;
  std::vector<std::optional<BoxRings>> boxRings;
  opt.box = resolveBox(opt, sta, mobs, seed, boxRings);
  if (opt.box.rotated) {
    for (auto* v : {&sta, &mobs, &refs})
      for (auto& m : *v)
        for (auto& a : m.atoms) a.xyz = opt.box.toBox(a.xyz);
    for (auto& p : poses)
      for (auto& rec : p)
        for (auto& xyz : rec) xyz = opt.box.toBox(xyz);
  }
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
  std::vector<od::Topology> mobTopo;
  for (std::size_t i = 0; i < mobs.size(); ++i) {
    od::logInfo(std::format("MobM from {}:", opt.mobmFiles[i]));
    printMolecule(mobs[i]);
    od::Topology t = analyze(mobs[i], opt);
    od::assignDssp({&mobs[i]}, {&t});
    printDssp(t);
    mobTopo.push_back(std::move(t));
  }

  // Molecule configuration (spec 5.1).
  std::vector<od::MolInput> inputs;
  for (std::size_t i = 0; i < sta.size(); ++i) inputs.push_back({&sta[i], &staTopo[i]});
  for (std::size_t i = 0; i < mobs.size(); ++i) inputs.push_back({&mobs[i], &mobTopo[i]});
  const std::string protText = protonationText(opt, protReports, overrides);
  if (opt.molConfig.empty()) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs, protText, opt.ringMacroMinSize);
    od::logInfo(std::format("Default molecule configuration written to {}\n"
                            "Review it (a copy may be edited) and run again with --molConfig <file>.",
                            path));
    return 0;
  }
  std::vector<std::string> errors;
  auto cfg = od::readMolConfig(opt.molConfig, inputs, topologyParams(opt), errors);
  if (!cfg) {
    const std::string path = od::writeDefaultMolConfig(opt.outPrefix, inputs, protText, opt.ringMacroMinSize);
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

  // score / localMin with run.unbound = input: no 9a (spec 6.1).
  const bool conf9a = opt.runMode == "dock" || opt.runUnbound == "conf9a";
  // Ring conformers of the small ring systems (spec 4.9).
  std::vector<od::RingConformers> ringConfs(inputs.size());
  if (conf9a && opt.ringMaxConfNum > 1) {
    const auto t0 = std::chrono::steady_clock::now();
    const od::RingParams rp = ringParams(opt, seed);
    od::logInfo("Ring conformers (spec 4.9):");
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      if (!cfg->molecules[i].ringConformers) continue;
      std::optional<BoxRings>* bp = i >= sta.size() ? &boxRings[i - sta.size()] : nullptr;
      if (bp && *bp && (*bp)->key == ringKey(cfg->molecules[i].topology)) {
        std::optional<BoxRings>& br = *bp;
        // Made for the box (input coordinates): conformer 0 follows the molecule into the box frame.
        ringConfs[i] = std::move(br->rc);
        br.reset();
        if (opt.box.rotated)
          for (auto& rs : ringConfs[i].systems)
            for (auto& x : rs.conf[0]) x = opt.box.toBox(x);
        od::logInfo(std::format("  {}: from the box (box.mobmRadius)", inputs[i].mol->id));
      } else {
        ringConfs[i] = od::generateRingConformers(*inputs[i].mol, cfg->molecules[i].topology, rp, (1ULL << 60) + i);
      }
      logRingConformers(*inputs[i].mol, ringConfs[i]);
    }
    od::logInfo(std::format("  ({:.2f} s)", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
  }

  {
    std::size_t staTorsions = 0, mobTorsions = 0;
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      std::size_t ringVars = 0;
      for (const auto& rs : ringConfs[i].systems) ringVars += rs.conf.size() > 1;
      (i < sta.size() ? staTorsions : mobTorsions) += cfg->molecules[i].tree.torsions.size() + ringVars;
    }
    if (staTorsions > 0 && conf9a) od::requireConfSampleNums(opt, od::SampleSet::StaM);
    if (mobTorsions > 0 && conf9a) od::requireConfSampleNums(opt, od::SampleSet::MobM);
  }

  od::logInfo("9a energy:");
  {
    std::vector<od::EnergyPart> parts;
    for (std::size_t i = 0; i < sta.size(); ++i) parts.push_back({&sta[i], &cfg->molecules[i].tree});
    printInputEnergy("StaM group", parts, opt);
  }
  for (std::size_t i = 0; i < mobs.size(); ++i)
    printInputEnergy(mobs[i].id, {{&mobs[i], &cfg->molecules[sta.size() + i].tree}}, opt);

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
  const bool gpu = selectGpu(opt);
  // protonation.recompute (spec 3.4.1): the [[protonation]] entries of each part of a group whose
  // hydrogens were placed at protonation.pH.
  const auto recomputeOverrides = [&](const std::string& file, bool stam, std::size_t parts)
      -> std::optional<std::vector<std::vector<od::ProtonationOverride>>> {
    const auto p = protOf(file, stam);
    if (!p) return std::nullopt;
    return std::vector<std::vector<od::ProtonationOverride>>(parts, p->overrides);
  };
  std::unique_ptr<ConfSet> stamSet;
  {
    std::vector<od::EnergyPart> parts;
    std::vector<const od::MolSettings*> settings;
    for (std::size_t i = 0; i < sta.size(); ++i) {
      parts.push_back({&sta[i], &cfg->molecules[i].tree, &ringConfs[i]});
      settings.push_back(&cfg->molecules[i]);
    }
    stamSet = std::make_unique<ConfSet>("StaM group", parts, settings, *cfg, opt, od::SampleSet::StaM, seed, 0, gpu,
                                        recomputeOverrides(opt.stamFile, true, sta.size()));
  }
  std::vector<std::unique_ptr<ConfSet>> mobSets;
  for (std::size_t i = 0; i < mobs.size(); ++i) {
    const auto& s = cfg->molecules[sta.size() + i];
    mobSets.push_back(std::make_unique<ConfSet>(mobs[i].id, std::vector<od::EnergyPart>{{&mobs[i], &s.tree, &ringConfs[sta.size() + i]}},
                                                std::vector<const od::MolSettings*>{&s}, *cfg, opt,
                                                od::SampleSet::MobM, seed, 1 + i, gpu,
                                                recomputeOverrides(opt.mobmFiles[i], false, 1)));
  }
  for (const ConfSet* c : {stamSet.get()})
    if (!c->movable()) od::logInfo(std::format("  {}: no rotatable torsions; the input conformation is the only sample", c->name()));
  for (const auto& c : mobSets)
    if (!c->movable()) od::logInfo(std::format("  {}: no rotatable torsions; the input conformation is the only sample", c->name()));
  if (opt.runMode == "dock") runRounds(*stamSet, mobSets, refs, opt, opt.fe, seed, ck, gpu);
  else runScore(*stamSet, mobSets, poses, refs, opt, seed, gpu);
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
