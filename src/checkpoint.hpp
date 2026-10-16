// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "log.hpp"

namespace od {

// Checkpoints at the end of each round and resuming (spec 16).

// SHA-256 as 64 lowercase hexadecimal digits.
std::string sha256Hex(const void* data, std::size_t size);
std::string sha256Hex(const std::string& s);
std::string sha256File(const std::string& path);

// Binary state: values in host byte order (the header records it).
class BinWriter {
 public:
  template <typename T>
  void put(const T& v) {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto* p = reinterpret_cast<const char*>(&v);
    buf_.append(p, sizeof(T));
  }
  template <typename T>
  void putVec(const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>);
    put<std::uint64_t>(v.size());
    buf_.append(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
  }
  const std::string& data() const { return buf_; }

 private:
  std::string buf_;
};

class BinReader {
 public:
  explicit BinReader(std::string data) : buf_(std::move(data)) {}
  template <typename T>
  T get() {
    static_assert(std::is_trivially_copyable_v<T>);
    need(sizeof(T));
    T v;
    std::memcpy(&v, buf_.data() + pos_, sizeof(T));
    pos_ += sizeof(T);
    return v;
  }
  template <typename T>
  std::vector<T> getVec() {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto n = get<std::uint64_t>();
    if (n > (buf_.size() - pos_) / sizeof(T)) fail("checkpoint state: truncated data");
    std::vector<T> v(n);
    std::memcpy(v.data(), buf_.data() + pos_, n * sizeof(T));
    pos_ += n * sizeof(T);
    return v;
  }
  bool atEnd() const { return pos_ == buf_.size(); }

 private:
  void need(std::size_t n) const {
    if (buf_.size() - pos_ < n) fail("checkpoint state: truncated data");
  }
  std::string buf_;
  std::size_t pos_ = 0;
};

// Header of state.bin: magic, byte-order mark and format version.
void putStateHeader(BinWriter& w);
void checkStateHeader(BinReader& r);

struct InputHash {
  std::string role;  // stamFile, mobmFiles, molConfig, ccd
  std::string path;  // absolute
  std::string sha256;
};

struct Manifest {
  std::string version;
  int round = 0;  // completed rounds
  std::vector<InputHash> inputs;
  std::string stateSha256;
};

std::string manifestText(const Manifest& m);
Manifest readManifest(const std::string& path);

// <outPrefix>/checkpoint/round_<r>
std::string checkpointDir(const std::string& outPrefix, int round);

// Writes files (name, content) to round_<r>.tmp, calls extra(temporary
// directory) for further files, renames it to round_<r> (replacing an existing
// one) and points `last` to it. Returns the directory.
std::string commitCheckpoint(const std::string& outPrefix, int round,
                                  const std::vector<std::pair<std::string, std::string>>& files,
                                  const std::function<void(const std::string& dir)>& extra);

// Largest r with a complete round_<r> (manifest.toml present) under
// <outPrefix>/checkpoint; 0 when there is none.
int latestCheckpoint(const std::string& outPrefix);

std::string readFile(const std::string& path);

}  // namespace od
