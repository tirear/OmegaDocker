// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <functional>
#include <span>

namespace od {

// Batched evaluation for sequential algorithms (spec 15.1.3): task(index, slot) runs for every
// index in [0, count) as a coroutine on the calling thread, at most maxActive at a time; slot
// (0 <= slot < maxActive) identifies the running coroutine. A task calls batchWait() after
// storing a request for its slot; when every running task waits, flush(slots) is called with
// the waiting slots (ascending) and must store their results, then the tasks continue.
// Tasks of different slots may run concurrently on OpenMP threads (between flushes); a task
// must only write data of its own index or slot. An exception of a task is rethrown by runBatched.
void runBatched(int count, int maxActive, const std::function<void(int index, int slot)>& task,
                const std::function<void(std::span<const int> slots)>& flush);
// Suspends the calling task until its request has been flushed.
void batchWait();

}  // namespace od
