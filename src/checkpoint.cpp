// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "checkpoint.hpp"

#include <array>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>

#include <toml.hpp>

namespace od {

namespace {

namespace fs = std::filesystem;

// FIPS 180-4.
class Sha256 {
 public:
  void update(const unsigned char* p, std::size_t n) {
    len_ += n;
    while (n > 0) {
      const std::size_t k = std::min(n, 64 - fill_);
      std::memcpy(block_.data() + fill_, p, k);
      fill_ += k;
      p += k;
      n -= k;
      if (fill_ == 64) {
        compress();
        fill_ = 0;
      }
    }
  }
  std::string hex() {
    const std::uint64_t bits = len_ * 8;
    const unsigned char one = 0x80, zero = 0;
    update(&one, 1);
    while (fill_ != 56) update(&zero, 1);
    for (int i = 7; i >= 0; --i) {
      const auto b = static_cast<unsigned char>(bits >> (8 * i));
      update(&b, 1);
    }
    std::string s;
    for (std::uint32_t v : h_) s += std::format("{:08x}", v);
    return s;
  }

 private:
  static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
  void compress() {
    static constexpr std::array<std::uint32_t, 64> k{
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
    std::array<std::uint32_t, 64> w;
    for (int i = 0; i < 16; ++i)
      w[i] = static_cast<std::uint32_t>(block_[4 * i]) << 24 | static_cast<std::uint32_t>(block_[4 * i + 1]) << 16 |
             static_cast<std::uint32_t>(block_[4 * i + 2]) << 8 | block_[4 * i + 3];
    for (int i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = h_;
    for (int i = 0; i < 64; ++i) {
      const std::uint32_t t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g)) + k[i] + w[i];
      const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    const std::array<std::uint32_t, 8> v{a, b, c, d, e, f, g, h};
    for (int i = 0; i < 8; ++i) h_[i] += v[i];
  }
  std::array<std::uint32_t, 8> h_{0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                  0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::array<unsigned char, 64> block_{};
  std::size_t fill_ = 0;
  std::uint64_t len_ = 0;
};

constexpr char kMagic[8] = {'O', 'D', 'C', 'K', 'P', 'T', 0, 0};
constexpr std::uint32_t kStateFormat = 2;
constexpr std::uint32_t kByteOrder = 0x01020304;

void writeFile(const std::string& path, const std::string& content) {
  std::ofstream f(path, std::ios::binary);
  if (!f) fail("cannot write {}", path);
  f.write(content.data(), static_cast<std::streamsize>(content.size()));
  f.close();
  if (!f) fail("cannot write {}", path);
}

}  // namespace

std::string sha256Hex(const void* data, std::size_t size) {
  Sha256 s;
  s.update(static_cast<const unsigned char*>(data), size);
  return s.hex();
}

std::string sha256Hex(const std::string& s) { return sha256Hex(s.data(), s.size()); }

std::string sha256File(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) fail("{}: cannot read", path);
  Sha256 s;
  std::vector<char> buf(1 << 20);
  while (f) {
    f.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    s.update(reinterpret_cast<const unsigned char*>(buf.data()), static_cast<std::size_t>(f.gcount()));
  }
  return s.hex();
}

std::string readFile(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) fail("{}: cannot read", path);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

void putStateHeader(BinWriter& w) {
  w.put(kMagic);
  w.put(kByteOrder);
  w.put(kStateFormat);
}

void checkStateHeader(BinReader& r) {
  const auto m = r.get<std::array<char, 8>>();
  if (std::memcmp(m.data(), kMagic, 8) != 0) fail("checkpoint state: not an OmegaDocker state file");
  if (r.get<std::uint32_t>() != kByteOrder) fail("checkpoint state: written on a machine with another byte order");
  const auto v = r.get<std::uint32_t>();
  if (v != kStateFormat) fail("checkpoint state: format {} (expected {})", v, kStateFormat);
}

std::string manifestText(const Manifest& m) {
  toml::table t;
  t.insert("version", m.version);
  t.insert("round", m.round);
  t.insert("state_sha256", m.stateSha256);
  toml::array in;
  for (const auto& h : m.inputs) in.push_back(toml::table{{"role", h.role}, {"path", h.path}, {"sha256", h.sha256}});
  t.insert("inputs", in);
  std::ostringstream s;
  s << "# OmegaDocker checkpoint (spec 16)\n" << toml::toml_formatter(t) << "\n";
  return s.str();
}

Manifest readManifest(const std::string& path) {
  toml::table t;
  try {
    t = toml::parse_file(path);
  } catch (const toml::parse_error& ex) {
    fail("{}: line {}: {}", path, ex.source().begin.line, ex.description());
  }
  Manifest m;
  auto str = [&](const toml::node_view<toml::node>& v, const char* key) {
    const auto s = v.value<std::string>();
    if (!s) fail("{}: {} missing", path, key);
    return *s;
  };
  m.version = str(t["version"], "version");
  m.stateSha256 = str(t["state_sha256"], "state_sha256");
  const auto r = t["round"].value<std::int64_t>();
  if (!r || *r < 1) fail("{}: round missing", path);
  m.round = static_cast<int>(*r);
  if (const auto* a = t["inputs"].as_array())
    for (auto& e : *a) {
      const auto* x = e.as_table();
      if (!x) fail("{}: inputs: tables expected", path);
      toml::table& tt = const_cast<toml::table&>(*x);
      m.inputs.push_back({str(tt["role"], "inputs.role"), str(tt["path"], "inputs.path"), str(tt["sha256"], "inputs.sha256")});
    }
  return m;
}

std::string checkpointDir(const std::string& outPrefix, int round) {
  return std::format("{}/checkpoint/round_{}", outPrefix, round);
}

std::string commitCheckpoint(const std::string& outPrefix, int round,
                             const std::vector<std::pair<std::string, std::string>>& files,
                             const std::function<void(const std::string& dir)>& extra) {
  const std::string dir = checkpointDir(outPrefix, round);
  const std::string tmp = dir + ".tmp";
  std::error_code ec;
  fs::remove_all(tmp, ec);
  fs::create_directories(tmp);
  if (extra) extra(tmp);
  for (const auto& [name, content] : files) writeFile(tmp + "/" + name, content);
  fs::remove_all(dir, ec);
  fs::rename(tmp, dir);
  const std::string base = outPrefix + "/checkpoint/";
  fs::remove(base + "last.tmp", ec);
  fs::create_directory_symlink(std::format("round_{}", round), base + "last.tmp");
  fs::rename(base + "last.tmp", base + "last");
  return dir;
}

int latestCheckpoint(const std::string& outPrefix) {
  const fs::path base = fs::path(outPrefix) / "checkpoint";
  int best = 0;
  std::error_code ec;
  if (!fs::is_directory(base, ec)) return 0;
  for (const auto& e : fs::directory_iterator(base)) {
    const std::string n = e.path().filename().string();
    if (!n.starts_with("round_") || !e.is_directory() || e.is_symlink()) continue;
    const std::string num = n.substr(6);
    if (num.empty() || num.find_first_not_of("0123456789") != std::string::npos || num.size() > 9) continue;
    if (!fs::exists(e.path() / "manifest.toml")) continue;
    best = std::max(best, std::stoi(num));
  }
  return best;
}

}  // namespace od
