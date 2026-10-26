// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

// Host memory of the process and its limit (memory.hostLimit, spec 13); CPU threads
// (cpu.threadNum, spec 13).
#pragma once
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

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

// CPU quota of the cgroup of this process in CPUs, rounded down and at least 1: the smallest of its
// cgroup and of its ancestors (cgroup v2 cpu.max, or cgroup v1 cpu.cfs_quota_us / cpu.cfs_period_us);
// nullopt when unlimited or unknown. source: the file of the quota.
std::optional<int> cgroupCpuLimit(std::string* source = nullptr, const std::string& procRoot = "/proc",
                                  const std::string& cgroupRoot = "/sys/fs/cgroup");

// cpus ordered by their index among the hardware threads of their core that are in cpus
// (thread_siblings_list under sysCpuRoot), then by number: one thread of every core before the
// second ones.
std::vector<int> orderCpusByCore(std::vector<int> cpus, const std::string& sysCpuRoot = "/sys/devices/system/cpu");

struct CpuThreads {
  int threads = 1;
  std::string source;
  std::vector<int> cpus;  // the CPUs this process is bound to (sched_setaffinity)
};

// The CPUs this process may use: those of the affinity mask, at most the cgroup CPU quota.
struct CpuLimit {
  int count = 0;          // 0: unknown
  std::vector<int> cpus;  // affinity mask
  std::optional<int> quota;
  std::string quotaSource;
};
CpuLimit cpuLimit();

// Error text when n threads exceed the limit; empty when within it (or the limit is unknown).
std::string cpuThreadsError(long long n, const CpuLimit& l);

// cpu.threadNum: n > 0 (more than cpuLimit() → Error), otherwise (auto) omp_get_max_threads()
// (OMP_NUM_THREADS or the CPUs of the affinity mask) limited by cpuLimit(). Sets the OpenMP and
// thread pool thread count and, when fewer than the CPUs of the affinity mask, binds the process to
// the first of them (orderCpusByCore). Call before any parallel work.
CpuThreads applyCpuThreads(long long n);

// "0-3,8": the CPUs as ranges.
std::string cpuListText(const std::vector<int>& cpus);

}  // namespace od
