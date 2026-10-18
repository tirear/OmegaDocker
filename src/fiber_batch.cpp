// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "fiber_batch.hpp"

#include <sys/mman.h>
#include <ucontext.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <memory>
#include <vector>

#include "log.hpp"

namespace od {
namespace {

constexpr std::size_t kFiberStack = 256 * 1024;

struct Fiber {
  ucontext_t ctx{};
  ucontext_t* sched = nullptr;  // context of the thread that resumed it last
  char* mem = nullptr;  // guard page + stack
  std::size_t memSize = 0;
  const std::function<void(int, int)>* task = nullptr;
  int slot = 0, index = -1;
  bool running = false, waiting = false;
  std::exception_ptr error;

  Fiber() = default;
  Fiber(const Fiber&) = delete;
  Fiber& operator=(const Fiber&) = delete;
  ~Fiber() {
    if (mem) munmap(mem, memSize);
  }
  void allocate() {
    const std::size_t page = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    memSize = kFiberStack + page;
    void* m = mmap(nullptr, memSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (m == MAP_FAILED) fail("cannot allocate a coroutine stack");
    mem = static_cast<char*>(m);
    mprotect(mem, page, PROT_NONE);
  }
};

thread_local Fiber* tCurrent = nullptr;

// Never returns: a finished task switches back to the thread that resumed it last
// (uc_link would name the thread that started it).
void fiberMain(unsigned lo, unsigned hi) {
  Fiber& fb = *reinterpret_cast<Fiber*>((static_cast<std::uintptr_t>(hi) << 32) | lo);
  try {
    (*fb.task)(fb.index, fb.slot);
  } catch (...) {
    fb.error = std::current_exception();
  }
  fb.running = false;
  setcontext(fb.sched);
}

void resume(Fiber& fb) {
  ucontext_t sched;
  Fiber* const outer = tCurrent;
  tCurrent = &fb;
  fb.sched = &sched;
  swapcontext(&sched, &fb.ctx);
  tCurrent = outer;
}

}  // namespace

void batchWait() {
  Fiber* fb = tCurrent;
  if (!fb) fail("batchWait outside runBatched");
  fb->waiting = true;
  swapcontext(&fb->ctx, fb->sched);
}

void runBatched(int count, int maxActive, const std::function<void(int, int)>& task,
                const std::function<void(std::span<const int>)>& flush) {
  if (count <= 0) return;
  const int nf = std::max(1, std::min(maxActive, count));
  std::vector<std::unique_ptr<Fiber>> fibers;
  for (int k = 0; k < nf; ++k) {
    fibers.push_back(std::make_unique<Fiber>());
    Fiber& fb = *fibers.back();
    fb.allocate();
    fb.task = &task;
    fb.slot = k;
  }
  std::atomic<int> next{0};
  std::atomic<bool> failed{false};
  // Starts indices on fb until one waits for a flush, none are left or a task failed.
  auto launch = [&](Fiber& fb) {
    while (!fb.running && !failed.load(std::memory_order_relaxed)) {
      const int i = next.fetch_add(1);
      if (i >= count) return;
      fb.index = i;
      fb.running = true;
      fb.waiting = false;
      getcontext(&fb.ctx);
      fb.ctx.uc_stack.ss_sp = fb.mem + (fb.memSize - kFiberStack);
      fb.ctx.uc_stack.ss_size = kFiberStack;
      fb.ctx.uc_link = nullptr;
      const auto ptr = reinterpret_cast<std::uintptr_t>(&fb);
      makecontext(&fb.ctx, reinterpret_cast<void (*)()>(fiberMain), 2, static_cast<unsigned>(ptr & 0xffffffffu),
                  static_cast<unsigned>(ptr >> 32));
      resume(fb);
      if (fb.error) failed = true;
    }
  };
  auto rethrow = [&] {
    for (auto& fb : fibers)
      if (fb->error) std::rethrow_exception(fb->error);
  };
  // The tasks of different slots are independent: they are resumed on several threads.
#pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < nf; ++k) launch(*fibers[k]);
  rethrow();
  std::vector<int> wait;
  for (;;) {
    wait.clear();
    for (int k = 0; k < nf; ++k)
      if (fibers[k]->running && fibers[k]->waiting) wait.push_back(k);
    if (wait.empty()) break;
    flush(wait);
    const int nw = static_cast<int>(wait.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int w = 0; w < nw; ++w) {
      Fiber& fb = *fibers[wait[w]];
      fb.waiting = false;
      resume(fb);
      if (fb.error) failed = true;
      else launch(fb);
    }
    rethrow();
  }
}

}  // namespace od
