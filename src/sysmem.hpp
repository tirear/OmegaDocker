// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

// Host memory of the process and its limit (memory.hostLimit, spec 13).
#pragma once
#include <cstddef>
#include <optional>
#include <string>

namespace od {

// Resident memory of this process (bytes); 0 when unknown.
std::size_t residentBytes(const std::string& procRoot = "/proc");

// Memory limit of the cgroup of this process: the smallest limit of its cgroup and of its
// ancestors (cgroup v2 memory.max, or cgroup v1 memory.limit_in_bytes); nullopt when unlimited or
// unknown. source: the file of the limit.
std::optional<std::size_t> cgroupMemLimit(std::string* source = nullptr, const std::string& procRoot = "/proc",
                                          const std::string& cgroupRoot = "/sys/fs/cgroup");

// MemAvailable of /proc/meminfo (bytes); 0 when unknown.
std::size_t memAvailableBytes(const std::string& procRoot = "/proc");

struct HostMemLimit {
  std::size_t bytes = 0;
  std::string source;
};

// memory.hostLimit: mib > 0, otherwise the cgroup limit, otherwise the resident memory plus
// MemAvailable.
HostMemLimit hostMemLimit(long long mib);

// Returns the freed heap memory to the system (glibc malloc_trim).
void trimHeap();

}  // namespace od
