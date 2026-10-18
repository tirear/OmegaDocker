// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "fiber_batch.hpp"

#include <sys/mman.h>
#include <unistd.h>
#if !defined(__x86_64__)
#include <ucontext.h>
#endif

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

#if defined(__x86_64__)
// Context switch without system calls (spec 15.1.9): saves the callee-saved registers, MXCSR and
// the x87 control word on the current stack, stores the stack pointer in *save and continues on
// the stack `to` (System V ABI: everything else is caller-saved).
extern "C" void odFiberSwitch(void** save, void* to);
extern "C" void odFiberEntry();
asm(R"(
  .text
  .globl odFiberSwitch
  .type odFiberSwitch, @function
odFiberSwitch:
  .cfi_startproc
  pushq %rbp
  pushq %rbx
  pushq %r12
  pushq %r13
  pushq %r14
  pushq %r15
  subq $16, %rsp
  stmxcsr 8(%rsp)
  fnstcw (%rsp)
  movq %rsp, (%rdi)
  movq %rsi, %rsp
  ldmxcsr 8(%rsp)
  fldcw (%rsp)
  addq $16, %rsp
  popq %r15
  popq %r14
  popq %r13
  popq %r12
  popq %rbx
  popq %rbp
  ret
  .cfi_endproc
  .size odFiberSwitch, .-odFiberSwitch

  .globl odFiberEntry
  .type odFiberEntry, @function
odFiberEntry:
  .cfi_startproc
  .cfi_undefined rip
  movq %r12, %rdi
  call odFiberMain@PLT
  ud2
  .cfi_endproc
  .size odFiberEntry, .-odFiberEntry
)");
#endif

struct Fiber {
#if defined(__x86_64__)
  void* sp = nullptr;
  void** sched = nullptr;  // saved stack pointer of the thread that resumed it last
#else
  ucontext_t ctx{};
  ucontext_t* sched = nullptr;  // context of the thread that resumed it last
#endif
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

void runTask(Fiber& fb) {
  try {
    (*fb.task)(fb.index, fb.slot);
  } catch (...) {
    fb.error = std::current_exception();
  }
  fb.running = false;
}

}  // namespace

#if defined(__x86_64__)
// Never returns: a finished task switches back to the thread that resumed it last.
extern "C" [[noreturn]] void odFiberMain(Fiber* fb) {
  runTask(*fb);
  void* unused;
  odFiberSwitch(&unused, *fb->sched);
  __builtin_unreachable();
}
#endif

namespace {

#if defined(__x86_64__)
// Initial stack of a new task: odFiberSwitch "returns" to odFiberEntry with r12 = &fb.
void prepare(Fiber& fb) {
  auto top = reinterpret_cast<std::uintptr_t>(fb.mem + fb.memSize) & ~std::uintptr_t(15);
  auto* s = reinterpret_cast<std::uint64_t*>(top);
  s[-1] = reinterpret_cast<std::uint64_t>(&odFiberEntry);  // return address
  for (int k = 2; k <= 7; ++k) s[-k] = 0;                    // rbp rbx r12 r13 r14 r15
  s[-4] = reinterpret_cast<std::uint64_t>(&fb);              // r12
  std::uint32_t mxcsr;
  std::uint16_t fpcw;
  asm volatile("stmxcsr %0" : "=m"(mxcsr));
  asm volatile("fnstcw %0" : "=m"(fpcw));
  s[-9] = fpcw;
  s[-8] = mxcsr;
  fb.sp = s - 9;
}
#else
// Never returns: a finished task switches back to the thread that resumed it last
// (uc_link would name the thread that started it).
void fiberMain(unsigned lo, unsigned hi) {
  Fiber& fb = *reinterpret_cast<Fiber*>((static_cast<std::uintptr_t>(hi) << 32) | lo);
  runTask(fb);
  setcontext(fb.sched);
}

void prepare(Fiber& fb) {
  getcontext(&fb.ctx);
  fb.ctx.uc_stack.ss_sp = fb.mem + (fb.memSize - kFiberStack);
  fb.ctx.uc_stack.ss_size = kFiberStack;
  fb.ctx.uc_link = nullptr;
  const auto ptr = reinterpret_cast<std::uintptr_t>(&fb);
  makecontext(&fb.ctx, reinterpret_cast<void (*)()>(fiberMain), 2, static_cast<unsigned>(ptr & 0xffffffffu),
              static_cast<unsigned>(ptr >> 32));
}
#endif

void resume(Fiber& fb) {
  Fiber* const outer = tCurrent;
  tCurrent = &fb;
#if defined(__x86_64__)
  void* sched;
  fb.sched = &sched;
  odFiberSwitch(&sched, fb.sp);
#else
  ucontext_t sched;
  fb.sched = &sched;
  swapcontext(&sched, &fb.ctx);
#endif
  tCurrent = outer;
}

}  // namespace

void batchWait() {
  Fiber* fb = tCurrent;
  if (!fb) fail("batchWait outside runBatched");
  fb->waiting = true;
#if defined(__x86_64__)
  odFiberSwitch(&fb->sp, *fb->sched);
#else
  swapcontext(&fb->ctx, fb->sched);
#endif
}

namespace {

// Fibers and the start of new indices shared by runBatched and runBatchedAsync.
struct Batch {
  int count;
  const std::function<void(int, int)>& task;
  std::vector<std::unique_ptr<Fiber>> fibers;
  std::atomic<int> next{0};
  std::atomic<bool> failed{false};

  Batch(int c, int nf, const std::function<void(int, int)>& t) : count(c), task(t) {
    for (int k = 0; k < nf; ++k) {
      fibers.push_back(std::make_unique<Fiber>());
      Fiber& fb = *fibers.back();
      fb.allocate();
      fb.task = &task;
      fb.slot = k;
    }
  }
  // Starts indices on fb until one waits for a flush, none are left or a task failed.
  void launch(Fiber& fb) {
    while (!fb.running && !failed.load(std::memory_order_relaxed)) {
      const int i = next.fetch_add(1);
      if (i >= count) return;
      fb.index = i;
      fb.running = true;
      fb.waiting = false;
      prepare(fb);
      resume(fb);
      if (fb.error) failed = true;
    }
  }
  // Resumes the waiting fibers of slots (they continue with new indices when done).
  void resumeAll(std::span<const int> slots) {
    const int nw = static_cast<int>(slots.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int w = 0; w < nw; ++w) {
      Fiber& fb = *fibers[slots[w]];
      fb.waiting = false;
      resume(fb);
      if (fb.error) failed = true;
      else launch(fb);
    }
  }
  void launchAll(int k0, int k1) {
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = k0; k < k1; ++k) launch(*fibers[k]);
  }
  void waiting(int k0, int k1, std::vector<int>& out) const {
    out.clear();
    for (int k = k0; k < k1; ++k)
      if (fibers[k]->running && fibers[k]->waiting) out.push_back(k);
  }
  void rethrow() const {
    for (auto& fb : fibers)
      if (fb->error) std::rethrow_exception(fb->error);
  }
};

}  // namespace

void runBatched(int count, int maxActive, const std::function<void(int, int)>& task,
                const std::function<void(std::span<const int>)>& flush) {
  if (count <= 0) return;
  const int nf = std::max(1, std::min(maxActive, count));
  Batch b(count, nf, task);
  // The tasks of different slots are independent: they are resumed on several threads.
  b.launchAll(0, nf);
  b.rethrow();
  std::vector<int> wait;
  for (;;) {
    b.waiting(0, nf, wait);
    if (wait.empty()) break;
    flush(wait);
    b.resumeAll(wait);
    b.rethrow();
  }
}

void runBatchedAsync(int count, int maxActive, int groups, const std::function<void(int, int)>& task,
                     const std::function<void(int, std::span<const int>)>& submit,
                     const std::function<void(int)>& complete) {
  if (count <= 0) return;
  const int nf = std::max(1, std::min(maxActive, count));
  const int ng = std::clamp(groups, 1, nf);
  Batch b(count, nf, task);
  auto lo = [&](int g) { return static_cast<int>(static_cast<long long>(nf) * g / ng); };
  std::vector<char> pending(ng, 0);
  std::vector<std::vector<int>> slots(ng);
  auto drain = [&] {
    for (int g = 0; g < ng; ++g)
      if (pending[g]) {
        pending[g] = 0;
        complete(g);
      }
  };
  try {
    b.launchAll(0, nf);
    b.rethrow();
    for (int g = 0; g < ng; ++g) {
      b.waiting(lo(g), lo(g + 1), slots[g]);
      if (slots[g].empty()) continue;
      submit(g, slots[g]);
      pending[g] = 1;
    }
    for (bool any = true; any;) {
      any = false;
      for (int g = 0; g < ng; ++g) {
        if (!pending[g]) continue;
        pending[g] = 0;
        complete(g);
        b.resumeAll(slots[g]);
        b.rethrow();
        b.waiting(lo(g), lo(g + 1), slots[g]);
        if (slots[g].empty()) continue;
        submit(g, slots[g]);
        pending[g] = 1;
        any = true;
      }
    }
  } catch (...) {
    drain();
    throw;
  }
}

}  // namespace od
