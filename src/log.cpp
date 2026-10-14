// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "log.hpp"

#include <atomic>
#include <cstdio>
#include <mutex>

namespace od {
namespace {
std::mutex gMutex;
std::atomic<int> gWarnings{0};
}  // namespace

void logInfo(const std::string& msg) {
  std::lock_guard lock(gMutex);
  std::fprintf(stdout, "%s\n", msg.c_str());
  std::fflush(stdout);
}

void logWarn(const std::string& msg) {
  ++gWarnings;
  std::lock_guard lock(gMutex);
  std::fflush(stdout);
  std::fprintf(stderr, "WARNING: %s\n", msg.c_str());
  std::fflush(stderr);
}

int warningCount() { return gWarnings.load(); }

}  // namespace od
