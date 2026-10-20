// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "pool.hpp"

#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace od {
namespace {

thread_local bool tInPool = false;

// Runs fn with OpenMP regions limited to one thread (restored afterwards).
template <class F>
void serialOmp(F&& fn) {
  const int saved = omp_get_max_threads();
  omp_set_num_threads(1);
  try {
    fn();
  } catch (...) {
    omp_set_num_threads(saved);
    throw;
  }
  omp_set_num_threads(saved);
}

class Pool {
 public:
  explicit Pool(int threads) : workers_(std::max(0, threads - 1)) {
    for (int w = 0; w < static_cast<int>(workers_.size()); ++w)
      workers_[w].thread = std::thread([this, w] { loop(w); });
  }
  ~Pool() {
    stop_ = true;
    for (Worker& w : workers_) {
      w.gen.fetch_add(1);
      w.gen.notify_one();
    }
    for (Worker& w : workers_) w.thread.join();
  }
  int threads() const { return static_cast<int>(workers_.size()) + 1; }

  void run(long long n, long long chunk, const std::function<void(long long, long long)>& fn) {
    if (n <= 0) return;
    if (chunk <= 0) chunk = std::max(1LL, (n + 4LL * threads() - 1) / (4LL * threads()));
    const long long ranges = (n + chunk - 1) / chunk;
    std::unique_lock lock(mutex_, std::try_to_lock);
    if (tInPool || !lock.owns_lock() || ranges <= 1 || workers_.empty()) {
      serialOmp([&] {
        for (long long lo = 0; lo < n; lo += chunk) fn(lo, std::min(n, lo + chunk));
      });
      return;
    }
    fn_ = &fn;
    n_ = n;
    chunk_ = chunk;
    next_.store(0);
    error_ = nullptr;
    const int k = static_cast<int>(std::min<long long>(static_cast<long long>(workers_.size()), ranges - 1));
    busy_.store(k);
    for (int w = 0; w < k; ++w) {
      workers_[w].gen.fetch_add(1);
      workers_[w].gen.notify_one();
    }
    tInPool = true;
    serialOmp([&] { work(); });
    tInPool = false;
    for (int b = busy_.load(); b != 0; b = busy_.load()) busy_.wait(b);
    fn_ = nullptr;
    if (error_) std::rethrow_exception(error_);
  }

 private:
  struct alignas(64) Worker {
    std::atomic<std::uint32_t> gen{0};
    std::thread thread;
  };

  void loop(int w) {
    tInPool = true;
    omp_set_num_threads(1);
    std::uint32_t seen = 0;
    for (;;) {
      workers_[w].gen.wait(seen);
      seen = workers_[w].gen.load();
      if (stop_) return;
      work();
      if (busy_.fetch_sub(1) == 1) busy_.notify_one();
    }
  }
  void work() {
    for (;;) {
      const long long lo = next_.fetch_add(chunk_);
      if (lo >= n_) return;
      try {
        (*fn_)(lo, std::min(n_, lo + chunk_));
      } catch (...) {
        std::lock_guard g(errorMutex_);
        if (!error_) error_ = std::current_exception();
        next_.store(n_);
      }
    }
  }

  std::vector<Worker> workers_;
  std::mutex mutex_, errorMutex_;
  const std::function<void(long long, long long)>* fn_ = nullptr;
  long long n_ = 0, chunk_ = 1;
  std::atomic<long long> next_{0};
  std::atomic<int> busy_{0};
  std::atomic<bool> stop_{false};
  std::exception_ptr error_;
};

Pool& pool() {
  static Pool p(omp_get_max_threads());
  return p;
}

}  // namespace

int poolThreads() { return pool().threads(); }

void parallelFor(long long n, long long chunk, const std::function<void(long long, long long)>& fn) {
  pool().run(n, chunk, fn);
}

}  // namespace od
