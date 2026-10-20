// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>

namespace od {

// Thread pool for work repeated at every minimization step (spec 15.1.25). Idle threads block
// (futex) instead of spinning, so a thread that loses its CPU to another program does not make
// the others wait actively. The pool has omp_get_max_threads() threads, the calling thread
// included; they are started at the first use.
int poolThreads();
// Calls fn(lo, hi) for ranges covering [0, n) of at most `chunk` indices (0: about four ranges per
// thread) on the pool and returns when all are done; the first exception of fn is rethrown. Runs
// on the calling thread alone when called from a pool task, or while another thread uses the pool.
// OpenMP regions inside fn run with one thread.
void parallelFor(long long n, long long chunk, const std::function<void(long long lo, long long hi)>& fn);
// fn(i) for every i in [0, n), one index per range (dynamic schedule).
template <class F>
void parallelEach(long long n, F&& fn) {
  parallelFor(n, 1, [&](long long lo, long long hi) {
    for (long long i = lo; i < hi; ++i) fn(i);
  });
}
// fn(i) for every i in [0, n), in ranges of about n / (4 × poolThreads()) indices.
template <class F>
void parallelRange(long long n, F&& fn) {
  parallelFor(n, 0, [&](long long lo, long long hi) {
    for (long long i = lo; i < hi; ++i) fn(i);
  });
}

}  // namespace od
