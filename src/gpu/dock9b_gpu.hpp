// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <thread>
#include <vector>

#include "dock9b.hpp"
#include "error.hpp"
#include "pool.hpp"
#include "energy_gpu.hpp"
#include "prep_gpu.hpp"

namespace od {

class GpuPoseEvaluator;

// Whole 9b minimizations taken over by other devices (spec 15.1.35): each running minimize of an
// evaluator with a hub is a call; an idle device (helpOnce) takes the last untaken poses of the
// call with the longest expected remaining time from its kernel and minimizes them with its own
// evaluator of the same MobM (key), writing the results into the call's output.
class DockStealHub {
 public:
  struct Call {
    GpuPoseEvaluator* ev = nullptr;
    int device = 0, key = 0;
    const PoseRequest* req = nullptr;
    DockMinParams p{};
    DockMinResult* out = nullptr;
    float* anglesOut = nullptr;
    std::size_t nt = 0;
    int n = 0;
    std::vector<const MobmDocker*> dks;
    int pending = 0;            // steals in progress
    bool closed = false;        // the kernel ended
    bool dry = false;           // a steal got nothing
    std::exception_ptr error;   // of a helper
  };
  void add(Call& c) {
    {
      std::lock_guard lk(m_);
      calls_.push_back(&c);
    }
    cv_.notify_all();
  }
  // After the call's kernel: waits for its steals; rethrows a helper's error.
  void finish(Call& c) {
    std::unique_lock lk(m_);
    c.closed = true;
    cv_.wait(lk, [&] { return c.pending == 0; });
    calls_.erase(std::find(calls_.begin(), calls_.end(), &c));
    lk.unlock();
    cv_.notify_all();
    if (c.error) std::rethrow_exception(c.error);
  }
  // One steal for device (evaluators ev(key)); false when no call is worth it.
  bool helpOnce(int device, const std::function<GpuPoseEvaluator*(int)>& ev);
  // Waits for a change of the calls, at most ms.
  void wait(int ms) {
    std::unique_lock lk(m_);
    cv_.wait_for(lk, std::chrono::milliseconds(ms));
  }
  bool anyCall() {
    std::lock_guard lk(m_);
    return !calls_.empty();
  }
  struct Stats {
    long long steals = 0, poses = 0;
    double bindSeconds = 0;
    std::size_t bindBytes = 0;
  };
  Stats stats() {
    std::lock_guard lk(m_);
    return stats_;
  }

 private:
  std::mutex m_;
  std::condition_variable cv_;
  std::vector<Call*> calls_;
  Stats stats_;
  double bandwidth_ = 6e9;  // bytes/s of binding the fields on another device (measured)
};

// 9b pose energies of one MobM on a GPU (spec 15.1.4); field is shared by the MobMs.
// prep (optional, shared by the MobMs): the quick filter and posm2_fft on the GPU (spec 15.1.5).
class GpuPoseEvaluator : public PoseEvaluator {
 public:
  // batch: requested poses per call, limited to 90% of the free device memory.
  // deviceMin: minimizations entirely on the GPU (minimize, spec 15.1.26, 15.1.28).
  GpuPoseEvaluator(const EnergyModel& model, std::shared_ptr<GpuField> field, int device, int batch,
                   std::shared_ptr<GpuDockPrep> prep = nullptr, bool deviceMin = true)
      : gpu_(model, device),
        field_(std::move(field)),
        prep_(std::move(prep)),
        nt_(model.torsionCount),
        device_(device),
        want_(std::max(1, batch)),
        deviceMin_(deviceMin) {
    batch_ = want_;
  }
  int maxBatch() const override { return batch_; }
  int torsionCount() const { return nt_; }
  int residentPoses(const MobmDocker& dk) const override {
    return gpu_.residentSamplesDock(static_cast<int>(dk.probes().size()));
  }
  void bind(const MobmDocker& dk) override {
    if (!field_->touch(dk.field().id())) field_->set(dk.field().model());
    slot_.assign(1, std::max(0, field_->slotOf(dk.field().id())));
    bound_.assign(1, &dk);
    setup(dk, dk.field().id());
  }
  int device() const { return device_; }
  // Minimizations taken over by other devices (spec 15.1.35); key: the evaluators of the same MobM.
  void setSteal(std::shared_ptr<DockStealHub> hub, int key) {
    hub_ = std::move(hub);
    stealKey_ = key;
  }
  // Device memory to copy when binding dks here (fields not resident; estimate).
  std::size_t bindBytes(const std::vector<const MobmDocker*>& dks) const {
    std::size_t b = 0;
    for (const MobmDocker* dk : dks) {
      const StamField& f = dk->field();
      if (field_->slotOf(f.id()) >= 0) continue;
      const std::size_t maps = static_cast<std::size_t>(f.typeMapCount()) + 4;
      if (f.mode() == GridMode::Delta) b += f.deltaBlockCount() * maps * kDeltaBlock * kDeltaBlock * kDeltaBlock * sizeof(float);
      else if (f.mode() == GridMode::Full && !f.hasHostData()) b += f.geom().size() * maps * sizeof(float);
    }
    return b;
  }
  // Binds dks for poses taken from another device (bindAll with at least dks.size() slots).
  void bindStolen(const std::vector<const MobmDocker*>& dks) {
    if (field_->slots() < static_cast<int>(dks.size())) field_->setSlots(static_cast<int>(dks.size()));
    bindAll(dks);
  }
  GpuEnergy& energy() { return gpu_; }
  // Resident fields: GpuField::slots() (spec 15.1.16).
  int fields() const override { return field_->slots(); }
  void bindAll(std::span<const MobmDocker* const> dks) override {
    if (static_cast<int>(dks.size()) > field_->slots()) throw Error("bindAll: more fields than GpuField slots");
    std::uint64_t key = 0;
    std::vector<char> have(dks.size());
    for (std::size_t i = 0; i < dks.size(); ++i) have[i] = field_->touch(dks[i]->field().id());
    for (std::size_t i = 0; i < dks.size(); ++i)
      if (!have[i]) field_->set(dks[i]->field().model());
    slot_.clear();
    for (const MobmDocker* dk : dks) {
      const int s = field_->slotOf(dk->field().id());
      if (s < 0) throw Error("bindAll: field not resident");
      slot_.push_back(s);
      key = key * 0x9e3779b97f4a7c15ULL + dk->field().id();
    }
    bound_.assign(dks.begin(), dks.end());
    setup(*dks[0], key);
  }
  void evaluate(int n, const PoseRequest* req, MobmDocker::PoseEnergy* out, double* grad) override {
    submit(0, n, req);
    complete(0, out, grad);
  }
  int lanes() const override { return GpuEnergy::kLanes; }
  void submit(int lane, int n, const PoseRequest* req) override {
    Work& w = work_[lane];
    const std::size_t nn = static_cast<std::size_t>(n), nt = static_cast<std::size_t>(nt_);
    w.n = n;
    w.angles.resize(nn * nt);
    w.poses.resize(nn);
    w.terms.resize(nn);
    w.dterms.resize(nn);
    parallelRange(n, [&](long long i) {
      const PoseRequest& r = req[i];
      for (std::size_t k = 0; k < nt; ++k) w.angles[i * nt + k] = static_cast<float>(r.anglesDeg[k]);
      DockPose& ps = w.poses[i];
      for (int c = 0; c < 3; ++c) {
        ps.t[c] = static_cast<float>(r.pose.t[c]);
        ps.center[c] = static_cast<float>(r.center[c]);
      }
      ps.q[0] = static_cast<float>(r.pose.q.w);
      ps.q[1] = static_cast<float>(r.pose.q.x);
      ps.q[2] = static_cast<float>(r.pose.q.y);
      ps.q[3] = static_cast<float>(r.pose.q.z);
      ps.field = slot_[static_cast<std::size_t>(r.field)];
    });
    gpu_.submitDock(lane, *field_, w.angles.data(), w.poses.data(), n, true);
  }
  void complete(int lane, MobmDocker::PoseEnergy* out, double* grad) override {
    Work& w = work_[lane];
    gpu_.completeDock(lane, w.terms.data(), w.dterms.data(), grad);
    parallelRange(w.n, [&](long long i) {
      const FlatTerms& f = w.terms[i];
      MobmDocker::PoseEnergy& pe = out[i];
      pe = {};
      pe.intra.vdw = f.vdw;
      pe.intra.hbond = f.hbond;
      pe.intra.elec = f.elec;
      pe.intra.desolv = f.desolv;
      pe.intra.clash = f.clash;
      pe.intra.closure = f.closure;
      pe.intra.clashPairs = f.clashPairs;
      pe.intra.closureMaxDeviation = f.closureMaxDeviation;
      pe.inter = w.dterms[i].inter;
      pe.interClash = w.dterms[i].interClash;
      pe.intra.external = pe.inter + pe.interClash;
      pe.outside = w.dterms[i].outside;
    });
  }
  bool minimize(int n, const PoseRequest* req, const DockMinParams& p, DockMinResult* out, float* anglesOut) override {
    if (!deviceMin_ || (p.method == static_cast<int>(MinMethod::Lbfgs) && p.memory > kDeviceLbfgsMaxMemory)) return false;
    const std::size_t nn = static_cast<std::size_t>(n), nt = static_cast<std::size_t>(nt_);
    minStarts_.resize(nn);
    minAngles_.resize(nn * nt);
    parallelRange(n, [&](long long i) {
      const PoseRequest& r = req[i];
      DockMinStart& s = minStarts_[i];
      for (int c = 0; c < 3; ++c) {
        s.t[c] = r.pose.t[c];
        s.center[c] = static_cast<float>(r.center[c]);
      }
      s.q[0] = r.pose.q.w;
      s.q[1] = r.pose.q.x;
      s.q[2] = r.pose.q.y;
      s.q[3] = r.pose.q.z;
      s.field = slot_[static_cast<std::size_t>(r.field)];
      std::copy_n(r.anglesDeg, nt, minAngles_.begin() + static_cast<std::ptrdiff_t>(i * nt));
    });
    if (!hub_) {
      gpu_.minimizeDock(*field_, p, minStarts_.data(), minAngles_.data(), n, out, anglesOut);
      return true;
    }
    DockStealHub::Call c;
    c.ev = this;
    c.device = device_;
    c.key = stealKey_;
    c.req = req;
    c.p = p;
    c.out = out;
    c.anglesOut = anglesOut;
    c.nt = nt;
    c.n = n;
    c.dks = bound_;
    hub_->add(c);
    try {
      gpu_.minimizeDock(*field_, p, minStarts_.data(), minAngles_.data(), n, out, anglesOut, true);
    } catch (...) {
      hub_->finish(c);
      throw;
    }
    hub_->finish(c);
    return true;
  }
  void idle() override {
    gpu_.releaseScratch();
    if (prep_) prep_->releaseWork();
  }

  bool quickFilter(const MobmDocker& dk, const MobmDocker::Conf& c, std::uint64_t stream, long long k0, long long k1,
                   std::vector<long long>& l2, Dock9bStats& stats) override {
    if (!prep_) return false;
    bind(dk);
    const Dock9bParams& p = dk.params();
    GpuDockPrep::FilterSetup s;
    for (const DockProbe& pr : dk.probes()) {
      for (int d = 0; d < 3; ++d) s.y.push_back(c.coords[pr.coord][d]);
      s.clash.push_back(pr.clash);
    }
    for (int d = 0; d < 3; ++d) {
      s.center[d] = static_cast<float>(c.center[d]);
      s.boxLo[d] = static_cast<float>(p.boxCenter[d] - 0.5 * p.boxSize[d]);
      s.boxHi[d] = static_cast<float>(p.boxCenter[d] + 0.5 * p.boxSize[d]);
    }
    s.clashDistance = static_cast<float>(p.clashDistance);
    s.nearDistance = static_cast<float>(p.maxNearDistance);
    s.maxOverlap = p.maxOverlapAtomNum;
    GpuDockPrep::FilterCounts cnt;
    if (p.posm == Posm::Posm2) {
      GpuDockPrep::Posm2Setup p2;
      for (const Quat& q : dk.hopf()) {
        const QuatF f = toFloat(q);
        p2.quats.insert(p2.quats.end(), {f.w, f.x, f.y, f.z});
      }
      for (const auto& r : c.hopfRange) {
        p2.range.insert(p2.range.end(), r.lo.begin(), r.lo.end());
        p2.range.insert(p2.range.end(), r.hi.begin(), r.hi.end());
      }
      p2.offset = c.hopfOffset;
      for (int d = 0; d < 3; ++d) p2.base[d] = p.boxCenter[d] - 0.5 * p.boxSize[d];
      p2.step = p.posStep;
      prep_->filterPosm2(*field_, s, p2, k0, k1, l2, cnt);
    } else {
      if (p.posm == Posm::Posm3 || p.posm == Posm::Posm4) {
        // posm3 part: the rotations from the CPU, the poses generated on the device (spec 15.1.21).
        const long long n3 = std::min(k1, c.sepaPoses());
        const long long nt = dk.posm3Translations();
        if (k0 < n3) {
          GpuDockPrep::Posm3Setup p3;
          for (int tr = 0; tr < nt; ++tr) {
            const Vec3 dl = dk.posm3Delta(tr);
            p3.delta.insert(p3.delta.end(), dl.begin(), dl.end());
          }
          dk.fitBounds(p3.fitLo, p3.fitHi);
          const long long rotChunk = std::max<long long>(1, (1LL << 22) / nt);
          for (long long r0 = k0 / nt; r0 * nt < n3; r0 += rotChunk) {
            const long long nr = std::min(rotChunk, (n3 + nt - 1) / nt - r0);
            p3.r0 = r0;
            p3.quats.resize(4 * static_cast<std::size_t>(nr));
            p3.base.resize(6 * static_cast<std::size_t>(nr));
            parallelRange(nr, [&](long long j) {
              const MobmDocker::Posm3Rot rot = dk.posm3Rot(c, r0 + j);
              const QuatF q = toFloat(rot.q);
              float* qd = p3.quats.data() + 4 * j;
              qd[0] = q.w;
              qd[1] = q.x;
              qd[2] = q.y;
              qd[3] = q.z;
              double* bd = p3.base.data() + 6 * j;
              for (int d = 0; d < 3; ++d) {
                bd[d] = rot.anchorS[d];
                bd[3 + d] = rot.rm[d];
              }
            });
            prep_->filterPosm3(*field_, s, p3, std::max(k0, r0 * nt), std::min(n3, (r0 + nr) * nt), l2, cnt);
          }
          k0 = n3;
        }
      }
      // Poses from the CPU (MobmDocker::l1Pose), in chunks.
      constexpr long long kChunk = 1LL << 20;
      std::vector<GpuDockPrep::FilterPose> poses;
      std::vector<char> valid;
      for (long long b = k0; b < k1; b += kChunk) {
        const long long m = std::min(kChunk, k1 - b);
        poses.resize(static_cast<std::size_t>(m));
        valid.assign(static_cast<std::size_t>(m), 0);
        parallelRange(m, [&](long long i) {
          const auto pose = dk.l1Pose(c, stream, b + i);
          if (!pose) return;
          GpuDockPrep::FilterPose& fp = poses[static_cast<std::size_t>(i)];
          for (int d = 0; d < 3; ++d) fp.t[d] = static_cast<float>(pose->t[d]);
          const QuatF q = toFloat(pose->q);
          fp.q[0] = q.w;
          fp.q[1] = q.x;
          fp.q[2] = q.y;
          fp.q[3] = q.z;
          fp.k = b + i;
          valid[static_cast<std::size_t>(i)] = 1;
        });
        std::size_t n = 0;
        for (std::size_t i = 0; i < poses.size(); ++i)
          if (valid[i]) poses[n++] = poses[i];
        cnt.outsideBox += m - static_cast<long long>(n);
        prep_->filter(*field_, s, static_cast<int>(n), poses.data(), l2, cnt);
      }
    }
    stats.outsideBox += cnt.outsideBox;
    stats.overlap += cnt.overlap;
    stats.far += cnt.far;
    return true;
  }
  int fftBatch(std::span<const MobmDocker* const> dks) override {
    if (!prep_) return 0;
    std::vector<const FftStam*> s;
    for (const MobmDocker* dk : dks) {
      if (!dk->fft()) return 0;
      s.push_back(dk->fft());
    }
    prep_->setFft(s);
    return prep_->fftBatch();
  }
  void fftTop(int n, int probes, const double* f, const int* map, const double* q, const std::array<int, 3>* lo,
              const std::array<int, 3>* hi, int top, std::pair<float, std::size_t>* out) override {
    prep_->fftTop(n, probes, f, map, q, lo, hi, top, out);
  }
  void prepared() override {
    if (prep_) prep_->releaseWork();
  }

 private:
  // Probe atoms of dk; the batch is sized once per set of fields (key), with this evaluator's
  // work memory freed.
  void setup(const MobmDocker& dk, std::uint64_t key) {
    const Dock9bParams& p = dk.params();
    gpu_.setProbes(dk.probes(), p.clashDistance, p.clashBase, p.clashK);
    if (sized_ == key) return;
    gpu_.releaseScratch();
    const std::size_t per = gpu_.bytesPerSampleDock();
    const std::size_t fit = static_cast<std::size_t>(0.9 * static_cast<double>(gpuFreeMem(device_))) / per;
    batch_ = static_cast<int>(std::clamp<std::size_t>(fit, 1, static_cast<std::size_t>(want_)));
    sized_ = key;
  }

  GpuEnergy gpu_;
  std::shared_ptr<GpuField> field_;
  std::shared_ptr<GpuDockPrep> prep_;
  int nt_, device_, want_, batch_;
  bool deviceMin_;
  std::vector<DockMinStart> minStarts_;
  std::vector<double> minAngles_;
  std::uint64_t sized_ = 0;
  std::vector<int> slot_{0};  // GpuField slot of each bound docker
  std::vector<const MobmDocker*> bound_;
  std::shared_ptr<DockStealHub> hub_;
  int stealKey_ = 0;
  // Host data of each lane.
  struct Work {
    int n = 0;
    std::vector<float> angles;
    std::vector<DockPose> poses;
    std::vector<FlatTerms> terms;
    std::vector<DockTerms> dterms;
  };
  Work work_[GpuEnergy::kLanes];
};

inline bool DockStealHub::helpOnce(int device, const std::function<GpuPoseEvaluator*(int)>& evOf) {
  std::unique_lock lk(m_);
  Call* best = nullptr;
  GpuPoseEvaluator* helper = nullptr;
  double bestGive = 0, bestFrac = 0;
  int bestMin = 0;
  std::size_t bestBytes = 0;
  for (Call* c : calls_) {
    if (c->closed || c->dry || c->device == device || c->dks.empty()) continue;
    int taken = 0, end = 0;
    double sec = 0;
    if (!c->ev->energy().dockProgress(taken, end, sec) || sec <= 0) continue;
    GpuPoseEvaluator* h = evOf(c->key);
    if (!h) continue;
    const double untaken = end - taken;
    // The kernel takes a wave of poses at its start: the rate from the poses after it.
    const int wave = h->residentPoses(*c->dks[0]);
    if (untaken < 2.0 * wave || taken <= wave) continue;
    const double rate = static_cast<double>(taken - wave) / sec;
    const std::size_t bytes = h->bindBytes(c->dks);
    const double delay = static_cast<double>(bytes) / bandwidth_ + 0.005;
    // Both devices end together: the helper starts after delay (same rate).
    const double give = (untaken - delay * rate) / 2;
    if (give < wave || give <= bestGive) continue;
    best = c;
    helper = h;
    bestGive = give;
    bestFrac = give / untaken;
    bestMin = wave;
    bestBytes = bytes;
  }
  if (!best) return false;
  ++best->pending;
  lk.unlock();
  const auto [lo, hi] = best->ev->energy().stealDock(bestFrac, bestMin);
  if (lo >= hi) {
    lk.lock();
    best->dry = true;
    --best->pending;
    lk.unlock();
    cv_.notify_all();
    return true;
  }
  std::exception_ptr err;
  double bindSec = 0;
  try {
    const auto t0 = std::chrono::steady_clock::now();
    helper->bindStolen(best->dks);
    bindSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!helper->minimize(hi - lo, best->req + lo, best->p, best->out + lo, best->anglesOut + static_cast<std::size_t>(lo) * best->nt))
      throw Error("DockStealHub: minimize on the device refused");
  } catch (...) {
    err = std::current_exception();
  }
  lk.lock();
  if (err && !best->error) best->error = err;
  ++stats_.steals;
  stats_.poses += hi - lo;
  stats_.bindSeconds += bindSec;
  stats_.bindBytes += bestBytes;
  if (bestBytes > (1u << 24) && bindSec > 0) bandwidth_ = 0.5 * bandwidth_ + 0.5 * static_cast<double>(bestBytes) / bindSec;
  --best->pending;
  lk.unlock();
  cv_.notify_all();
  return true;
}

// 9b pose energies of one MobM on several GPUs (spec 15.1.33): the whole minimizations (minimize)
// and the posm2_fft translations (fftTop) are divided among the evaluators, each on its own host
// thread; the other work (evaluate, quick filter) runs on the first one.
class MultiGpuPoseEvaluator : public PoseEvaluator {
 public:
  // hub (optional): an evaluator done with its part takes over poses of the others (spec 15.1.35);
  // the evaluators then have it (setSteal).
  explicit MultiGpuPoseEvaluator(std::vector<std::unique_ptr<GpuPoseEvaluator>> ev,
                                 std::shared_ptr<DockStealHub> hub = nullptr)
      : ev_(std::move(ev)), hub_(std::move(hub)) {
    if (ev_.empty()) throw Error("MultiGpuPoseEvaluator: no evaluator");
  }
  int maxBatch() const override { return ev_[0]->maxBatch(); }
  int residentPoses(const MobmDocker& dk) const override {
    int n = 0;
    for (const auto& e : ev_) n += e->residentPoses(dk);
    return n;
  }
  void bind(const MobmDocker& dk) override {
    for (auto& e : ev_) e->bind(dk);
  }
  int fields() const override {
    int n = ev_[0]->fields();
    for (const auto& e : ev_) n = std::min(n, e->fields());
    return n;
  }
  void bindAll(std::span<const MobmDocker* const> dks) override {
    for (auto& e : ev_) e->bindAll(dks);
  }
  void evaluate(int n, const PoseRequest* req, MobmDocker::PoseEnergy* out, double* grad) override {
    ev_[0]->evaluate(n, req, out, grad);
  }
  int lanes() const override { return ev_[0]->lanes(); }
  void submit(int lane, int n, const PoseRequest* req) override { ev_[0]->submit(lane, n, req); }
  void complete(int lane, MobmDocker::PoseEnergy* out, double* grad) override { ev_[0]->complete(lane, out, grad); }
  bool minimize(int n, const PoseRequest* req, const DockMinParams& p, DockMinResult* out, float* anglesOut) override {
    const std::size_t k = ev_.size();
    if (k == 1 || n < 2) return ev_[0]->minimize(n, req, p, out, anglesOut);
    std::vector<int> begin(k + 1, 0);
    for (std::size_t i = 1; i <= k; ++i) begin[i] = static_cast<int>(static_cast<long long>(n) * static_cast<long long>(i) / static_cast<long long>(k));
    const std::size_t nt = static_cast<std::size_t>(ev_[0]->torsionCount());
    std::vector<char> ok(k, 0);
    std::vector<std::exception_ptr> err(k);
    std::atomic<std::size_t> left{k};
    const auto run = [&](std::size_t i) {
      try {
        const int b = begin[i], m = begin[i + 1] - b;
        ok[i] = m == 0 || ev_[i]->minimize(m, req + b, p, out + b, anglesOut + static_cast<std::size_t>(b) * nt);
      } catch (...) {
        err[i] = std::current_exception();
      }
      --left;
      // Poses of the parts still running.
      if (hub_ && ok[i] && !err[i])
        try {
          const int dev = ev_[i]->device();
          const auto evOf = [&](int) { return ev_[i].get(); };
          while (left.load() > 0)
            if (!hub_->helpOnce(dev, evOf)) hub_->wait(1);
        } catch (...) {
          err[i] = std::current_exception();
        }
    };
    std::vector<std::thread> th;
    for (std::size_t i = 1; i < k; ++i) th.emplace_back(run, i);
    run(0);
    for (auto& x : th) x.join();
    for (const auto& e : err)
      if (e) std::rethrow_exception(e);
    for (char o : ok)
      if (!o) return false;
    return true;
  }
  void idle() override {
    for (auto& e : ev_) e->idle();
  }
  bool quickFilter(const MobmDocker& dk, const MobmDocker::Conf& c, std::uint64_t stream, long long k0, long long k1,
                   std::vector<long long>& l2, Dock9bStats& stats) override {
    return ev_[0]->quickFilter(dk, c, stream, k0, k1, l2, stats);
  }
  // posm2_fft: the orientations of each fftTop divided among the evaluators in proportion to their
  // batches (spec 15.1.35); all on the first when one has no quick filter on its device.
  int fftBatch(std::span<const MobmDocker* const> dks) override {
    fftParts_.clear();
    fftSets_ = static_cast<int>(dks.size());
    int total = 0;
    for (auto& e : ev_) {
      const int b = e->fftBatch(dks);
      if (b <= 0) {
        fftParts_.clear();
        return ev_[0]->fftBatch(dks);
      }
      fftParts_.push_back(b);
      total += b;
    }
    return total;
  }
  void fftTop(int n, int probes, const double* f, const int* map, const double* q, const std::array<int, 3>* lo,
              const std::array<int, 3>* hi, int top, std::pair<float, std::size_t>* out) override {
    const std::size_t k = fftParts_.size();
    if (k < 2 || n < 2) {
      ev_[0]->fftTop(n, probes, f, map, q, lo, hi, top, out);
      return;
    }
    long long sum = 0;
    for (int b : fftParts_) sum += b;
    std::vector<int> begin(k + 1, 0);
    long long acc = 0;
    for (std::size_t i = 0; i < k; ++i) {
      acc += fftParts_[i];
      begin[i + 1] = static_cast<int>(static_cast<long long>(n) * acc / sum);
    }
    const std::size_t sets = static_cast<std::size_t>(fftSets_), tp = static_cast<std::size_t>(top);
    std::vector<std::vector<std::pair<float, std::size_t>>> part(k);
    std::vector<std::exception_ptr> err(k);
    const auto run = [&](std::size_t i) {
      try {
        const int b = begin[i], m = begin[i + 1] - b;
        if (m == 0) return;
        part[i].resize(sets * static_cast<std::size_t>(m) * tp);
        ev_[i]->fftTop(m, probes, f + static_cast<std::size_t>(b) * probes * 3, map, q, lo + b, hi + b, top, part[i].data());
      } catch (...) {
        err[i] = std::current_exception();
      }
    };
    std::vector<std::thread> th;
    for (std::size_t i = 1; i < k; ++i) th.emplace_back(run, i);
    run(0);
    for (auto& x : th) x.join();
    for (const auto& e : err)
      if (e) std::rethrow_exception(e);
    // out: sets × n × top.
    for (std::size_t i = 0; i < k; ++i) {
      const std::size_t b = static_cast<std::size_t>(begin[i]), m = static_cast<std::size_t>(begin[i + 1]) - b;
      for (std::size_t s = 0; s < sets; ++s)
        std::copy_n(part[i].begin() + static_cast<std::ptrdiff_t>(s * m * tp), m * tp,
                    out + (s * static_cast<std::size_t>(n) + b) * tp);
    }
  }
  void prepared() override {
    for (auto& e : ev_) e->prepared();
  }

 private:
  std::vector<std::unique_ptr<GpuPoseEvaluator>> ev_;
  std::shared_ptr<DockStealHub> hub_;
  std::vector<int> fftParts_;  // fftBatch of each evaluator (empty: all on the first)
  int fftSets_ = 1;
};

}  // namespace od
