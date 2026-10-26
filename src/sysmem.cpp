// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "sysmem.hpp"
#include "error.hpp"
#include <format>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>

#include <malloc.h>
#include <omp.h>
#include <sched.h>
#include <unistd.h>

namespace od {

namespace {

std::optional<std::size_t> readLimit(const std::string& file) {
  std::ifstream in(file);
  std::string s;
  if (!(in >> s) || s == "max") return std::nullopt;
  try {
    const unsigned long long v = std::stoull(s);
    // cgroup v1 writes a huge page-aligned number for "unlimited".
    if (v >= (1ULL << 60)) return std::nullopt;
    return static_cast<std::size_t>(v);
  } catch (...) {
    return std::nullopt;
  }
}

std::string parentDir(const std::string& p) {
  const std::size_t k = p.find_last_of('/');
  return k == std::string::npos || k == 0 ? std::string("/") : p.substr(0, k);
}

// Paths of the cgroup of this process: cgroup v2 ("0::<path>") and cgroup v1 of the controller ctl.
struct CgroupPaths {
  std::string v2, v1;
  bool hasV1 = false;
};

CgroupPaths cgroupPaths(const std::string& procRoot, const char* ctlName) {
  std::ifstream in(procRoot + "/self/cgroup");
  CgroupPaths r;
  std::string line;
  // Lines "<id>:<controllers>:<path>".
  while (std::getline(in, line)) {
    const std::size_t a = line.find(':'), b = a == std::string::npos ? a : line.find(':', a + 1);
    if (b == std::string::npos) continue;
    const std::string ctl = line.substr(a + 1, b - a - 1), path = line.substr(b + 1);
    if (line.compare(0, a, "0") == 0 && ctl.empty()) r.v2 = path;
    std::stringstream ss(ctl);
    for (std::string c; std::getline(ss, c, ',');)
      if (c == ctlName) {
        r.v1 = path;
        r.hasV1 = true;
      }
  }
  return r;
}

// Calls fn(dir) for the cgroup directory of path under base and for those of its ancestors.
template <class Fn>
void walkCgroup(const std::string& base, std::string path, Fn&& fn) {
  if (path.empty()) path = "/";
  while (true) {
    fn(path == "/" ? base : base + path);
    if (path == "/") break;
    path = parentDir(path);
  }
}

std::optional<long long> readInt(const std::string& file) {
  std::ifstream in(file);
  long long v = 0;
  if (!(in >> v)) return std::nullopt;
  return v;
}

// "0,8", "0-1", "4": the listed numbers.
std::vector<int> parseCpuList(const std::string& s) {
  std::vector<int> r;
  std::stringstream ss(s);
  for (std::string part; std::getline(ss, part, ',');) {
    try {
      const std::size_t d = part.find('-');
      const int a = std::stoi(part.substr(0, d)), b = d == std::string::npos ? a : std::stoi(part.substr(d + 1));
      for (int c = a; c <= b; ++c) r.push_back(c);
    } catch (...) {
    }
  }
  return r;
}

}  // namespace

std::size_t residentBytes(const std::string& procRoot) {
  std::ifstream in(procRoot + "/self/statm");
  unsigned long long size = 0, resident = 0;
  if (!(in >> size >> resident)) return 0;
  return static_cast<std::size_t>(resident) * static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
}

std::optional<std::size_t> cgroupMemLimit(std::string* source, const std::string& procRoot,
                                          const std::string& cgroupRoot) {
  const CgroupPaths cp = cgroupPaths(procRoot, "memory");
  std::optional<std::size_t> best;
  const auto check = [&](const std::string& file) {
    if (const auto v = readLimit(file); v && (!best || *v < *best)) {
      best = v;
      if (source) *source = file;
    }
  };
  if (cp.hasV1) walkCgroup(cgroupRoot + "/memory", cp.v1, [&](const std::string& d) { check(d + "/memory.limit_in_bytes"); });
  else if (!cp.v2.empty()) walkCgroup(cgroupRoot, cp.v2, [&](const std::string& d) { check(d + "/memory.max"); });
  return best;
}

std::optional<int> cgroupCpuLimit(std::string* source, const std::string& procRoot, const std::string& cgroupRoot) {
  const CgroupPaths cp = cgroupPaths(procRoot, "cpu");
  std::optional<double> best;
  const auto take = [&](double cpus, const std::string& file) {
    if (cpus > 0 && (!best || cpus < *best)) {
      best = cpus;
      if (source) *source = file;
    }
  };
  if (cp.hasV1)
    walkCgroup(cgroupRoot + "/cpu", cp.v1, [&](const std::string& d) {
      const auto q = readInt(d + "/cpu.cfs_quota_us"), p = readInt(d + "/cpu.cfs_period_us");
      if (q && p && *q > 0 && *p > 0) take(static_cast<double>(*q) / static_cast<double>(*p), d + "/cpu.cfs_quota_us");
    });
  else if (!cp.v2.empty())
    walkCgroup(cgroupRoot, cp.v2, [&](const std::string& d) {
      // "max <period>" or "<quota> <period>".
      std::ifstream in(d + "/cpu.max");
      std::string q;
      long long p = 0;
      if (!(in >> q >> p) || q == "max" || p <= 0) return;
      try {
        take(static_cast<double>(std::stoll(q)) / static_cast<double>(p), d + "/cpu.max");
      } catch (...) {
      }
    });
  if (!best) return std::nullopt;
  return std::max(1, static_cast<int>(std::floor(*best + 1e-9)));
}

std::vector<int> orderCpusByCore(std::vector<int> cpus, const std::string& sysCpuRoot) {
  std::vector<std::pair<int, int>> key;  // (index among the threads of the core, cpu)
  for (int c : cpus) {
    std::ifstream in(sysCpuRoot + "/cpu" + std::to_string(c) + "/topology/thread_siblings_list");
    std::string s;
    std::getline(in, s);
    std::vector<int> sib;
    for (int x : parseCpuList(s))
      if (std::find(cpus.begin(), cpus.end(), x) != cpus.end()) sib.push_back(x);
    std::sort(sib.begin(), sib.end());
    const auto it = std::find(sib.begin(), sib.end(), c);
    key.emplace_back(it == sib.end() ? 0 : static_cast<int>(it - sib.begin()), c);
  }
  std::sort(key.begin(), key.end());
  for (std::size_t i = 0; i < key.size(); ++i) cpus[i] = key[i].second;
  return cpus;
}

CpuLimit cpuLimit() {
  CpuLimit l;
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) == 0)
    for (int c = 0; c < CPU_SETSIZE; ++c)
      if (CPU_ISSET(c, &set)) l.cpus.push_back(c);
  l.quota = cgroupCpuLimit(&l.quotaSource);
  l.count = static_cast<int>(l.cpus.size());
  if (l.quota && (l.count == 0 || *l.quota < l.count)) l.count = *l.quota;
  return l;
}

std::string cpuThreadsError(long long n, const CpuLimit& l) {
  if (l.count <= 0 || n <= l.count) return {};
  std::string s = std::format("cpu.threadNum {} is more than the {} CPU(s) this process may use (affinity mask: ", n, l.count);
  s += l.cpus.empty() ? std::string("unknown") : std::format("{} CPU(s) {}", l.cpus.size(), cpuListText(l.cpus));
  s += l.quota ? std::format("; cgroup CPU quota {} ({})", *l.quota, l.quotaSource) : std::string("; no cgroup CPU quota");
  return s + ")";
}

CpuThreads applyCpuThreads(long long n) {
  const CpuLimit l = cpuLimit();
  const std::vector<int>& avail = l.cpus;
  cpu_set_t set;
  CpuThreads r;
  if (n > 0) {
    if (const std::string e = cpuThreadsError(n, l); !e.empty()) throw Error(e);
    r.threads = static_cast<int>(n);
    r.source = "cpu.threadNum";
  } else {
    r.threads = omp_get_max_threads();
    r.source = std::getenv("OMP_NUM_THREADS") ? "OMP_NUM_THREADS" : "affinity mask";
    if (l.quota && *l.quota < r.threads) {
      r.threads = *l.quota;
      r.source = "cgroup " + l.quotaSource;
    }
    if (!avail.empty() && r.threads > static_cast<int>(avail.size())) {
      r.threads = static_cast<int>(avail.size());
      r.source = "affinity mask";
    }
  }
  r.cpus = avail;
  if (!avail.empty() && r.threads < static_cast<int>(avail.size())) {
    std::vector<int> pick = orderCpusByCore(avail);
    pick.resize(r.threads);
    std::sort(pick.begin(), pick.end());
    CPU_ZERO(&set);
    for (int c : pick) CPU_SET(c, &set);
    if (sched_setaffinity(0, sizeof(set), &set) == 0) r.cpus = pick;
  }
  omp_set_num_threads(r.threads);
  return r;
}

std::size_t memAvailableBytes(const std::string& procRoot) {
  std::ifstream in(procRoot + "/meminfo");
  std::string key, unit;
  unsigned long long v = 0;
  while (in >> key >> v) {
    std::getline(in, unit);
    if (key == "MemAvailable:") return static_cast<std::size_t>(v) * 1024;
  }
  return 0;
}

HostMemLimit hostMemLimit(long long mib) {
  if (mib > 0) return {static_cast<std::size_t>(mib) << 20, "memory.hostLimit"};
  std::string src;
  if (const auto v = cgroupMemLimit(&src)) return {*v, "cgroup " + src};
  const std::size_t avail = memAvailableBytes();
  if (avail == 0) return {std::numeric_limits<std::size_t>::max(), "unknown (no limit)"};
  return {residentBytes() + avail, "resident + MemAvailable"};
}

void trimHeap() { malloc_trim(0); }

std::string cpuListText(const std::vector<int>& cpus) {
  std::string s;
  for (std::size_t i = 0; i < cpus.size();) {
    std::size_t j = i;
    while (j + 1 < cpus.size() && cpus[j + 1] == cpus[j] + 1) ++j;
    if (!s.empty()) s += ',';
    s += std::to_string(cpus[i]);
    if (j > i) s += '-' + std::to_string(cpus[j]);
    i = j + 1;
  }
  return s.empty() ? "unknown" : s;
}

}  // namespace od
