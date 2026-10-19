// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "sysmem.hpp"

#include <fstream>
#include <limits>
#include <sstream>

#include <malloc.h>
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

}  // namespace

std::size_t residentBytes(const std::string& procRoot) {
  std::ifstream in(procRoot + "/self/statm");
  unsigned long long size = 0, resident = 0;
  if (!(in >> size >> resident)) return 0;
  return static_cast<std::size_t>(resident) * static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
}

std::optional<std::size_t> cgroupMemLimit(std::string* source, const std::string& procRoot,
                                          const std::string& cgroupRoot) {
  std::ifstream in(procRoot + "/self/cgroup");
  std::string line, v2, v1;
  bool hasV1 = false;
  // Lines "<id>:<controllers>:<path>": cgroup v2 "0::<path>", cgroup v1 with "memory" in the controllers.
  while (std::getline(in, line)) {
    const std::size_t a = line.find(':'), b = a == std::string::npos ? a : line.find(':', a + 1);
    if (b == std::string::npos) continue;
    const std::string ctl = line.substr(a + 1, b - a - 1), path = line.substr(b + 1);
    if (line.compare(0, a, "0") == 0 && ctl.empty()) v2 = path;
    std::stringstream ss(ctl);
    for (std::string c; std::getline(ss, c, ',');)
      if (c == "memory") {
        v1 = path;
        hasV1 = true;
      }
  }
  std::optional<std::size_t> best;
  const auto walk = [&](const std::string& base, std::string path, const char* file) {
    while (true) {
      const std::string dir = path == "/" ? base : base + path;
      if (const auto v = readLimit(dir + "/" + file); v && (!best || *v < *best)) {
        best = v;
        if (source) *source = dir + "/" + file;
      }
      if (path == "/" || path.empty()) break;
      path = parentDir(path);
    }
  };
  if (hasV1) walk(cgroupRoot + "/memory", v1.empty() ? "/" : v1, "memory.limit_in_bytes");
  else if (!v2.empty()) walk(cgroupRoot, v2, "memory.max");
  return best;
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

}  // namespace od
