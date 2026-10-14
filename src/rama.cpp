// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "rama.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <mutex>
#include <string_view>

#include "log.hpp"

extern "C" {
#define OD_RAMA_SYMBOLS(name) extern const char od_rama_##name##_begin[], od_rama_##name##_end[];
OD_RAMA_SYMBOLS(general)
OD_RAMA_SYMBOLS(gly)
OD_RAMA_SYMBOLS(transpro)
OD_RAMA_SYMBOLS(cispro)
OD_RAMA_SYMBOLS(prepro)
OD_RAMA_SYMBOLS(ileval)
#undef OD_RAMA_SYMBOLS
}

namespace od {
namespace {

struct Source {
  const char* begin;
  const char* end;
  const char* file;
};

const std::array<Source, kRamaClassCount> kSources = {{
    {od_rama_general_begin, od_rama_general_end, "rama8000-general-noGPIVpreP.data"},
    {od_rama_gly_begin, od_rama_gly_end, "rama8000-gly-sym.data"},
    {od_rama_transpro_begin, od_rama_transpro_end, "rama8000-transpro.data"},
    {od_rama_cispro_begin, od_rama_cispro_end, "rama8000-cispro.data"},
    {od_rama_prepro_begin, od_rama_prepro_end, "rama8000-prepro-noGP.data"},
    {od_rama_ileval_begin, od_rama_ileval_end, "rama8000-ileval-nopreP.data"},
}};

bool parseNumber(std::string_view& s, double& v) {
  while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.remove_prefix(1);
  auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
  if (ec != std::errc()) return false;
  s.remove_prefix(static_cast<std::size_t>(p - s.data()));
  return true;
}

// Lines are "<phi> <psi> <value>" at cell centers (odd degrees); '#' starts a comment line.
RamaGrid parse(const Source& src) {
  RamaGrid g;
  g.value.assign(kRamaBins * kRamaBins, 0.0f);
  std::string_view text(src.begin, static_cast<std::size_t>(src.end - src.begin));
  int lineNo = 0;
  while (!text.empty()) {
    const auto nl = text.find('\n');
    std::string_view line = text.substr(0, nl);
    text.remove_prefix(nl == std::string_view::npos ? text.size() : nl + 1);
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    if (line.empty() || line.front() == '#') continue;
    double phi = 0, psi = 0, v = 0;
    if (!parseNumber(line, phi) || !parseNumber(line, psi) || !parseNumber(line, v) || !(v >= 0) ||
        std::fmod(phi + 181.0, 2.0) != 0 || std::fmod(psi + 181.0, 2.0) != 0 || phi < -179 || phi > 179 ||
        psi < -179 || psi > 179)
      fail("embedded {}: line {}: not a valid grid entry", src.file, lineNo);
    float& cell = g.value[RamaGrid::bin(phi) * kRamaBins + RamaGrid::bin(psi)];
    if (cell != 0.0f) fail("embedded {}: line {}: cell listed twice", src.file, lineNo);
    cell = static_cast<float>(v);
    ++g.listed;
  }
  if (g.listed == 0) fail("embedded {}: no grid entries", src.file);
  return g;
}

}  // namespace

int RamaGrid::bin(double angleDeg) {
  double a = std::fmod(angleDeg + 180.0, 360.0);
  if (a < 0) a += 360.0;
  const int b = static_cast<int>(a / 2.0);
  return b >= kRamaBins ? b - kRamaBins : b;
}

const RamaGrid& ramaGrid(RamaClass c) {
  static std::array<RamaGrid, kRamaClassCount> grids;
  static std::array<std::once_flag, kRamaClassCount> once;
  const auto i = static_cast<std::size_t>(c);
  std::call_once(once[i], [&] { grids[i] = parse(kSources[i]); });
  return grids[i];
}

const char* ramaTypeName(RamaClass c) {
  switch (c) {
    case RamaClass::General: return "RAMA_GENERAL";
    case RamaClass::Gly: return "RAMA_GLY";
    case RamaClass::TransPro: return "RAMA_TRANSPRO";
    case RamaClass::CisPro: return "RAMA_CISPRO";
    case RamaClass::PrePro: return "RAMA_PREPRO";
    case RamaClass::IleVal: return "RAMA_ILEVAL";
  }
  return "";
}

const char* ramaClassLabel(RamaClass c) {
  switch (c) {
    case RamaClass::General: return "general (not Gly, Pro, Ile, Val or pre-Pro)";
    case RamaClass::Gly: return "Gly";
    case RamaClass::TransPro: return "trans-Pro";
    case RamaClass::CisPro: return "cis-Pro";
    case RamaClass::PrePro: return "pre-Pro (not Gly or Pro)";
    case RamaClass::IleVal: return "Ile / Val (not pre-Pro)";
  }
  return "";
}

std::optional<RamaClass> ramaClassFromTypeName(const std::string& name) {
  for (int i = 0; i < kRamaClassCount; ++i)
    if (name == ramaTypeName(static_cast<RamaClass>(i))) return static_cast<RamaClass>(i);
  return std::nullopt;
}

RamaClass ramaClassOf(const std::string& resName, const std::string& nextResName, bool cisBefore) {
  if (resName == "GLY") return RamaClass::Gly;
  if (resName == "PRO") return cisBefore ? RamaClass::CisPro : RamaClass::TransPro;
  if (nextResName == "PRO") return RamaClass::PrePro;
  if (resName == "ILE" || resName == "VAL") return RamaClass::IleVal;
  return RamaClass::General;
}

}  // namespace od
