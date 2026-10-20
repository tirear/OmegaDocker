// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "fft_dock.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>

#include "pocketfft_hdronly.h"

namespace od {

namespace {

// C-order shape (z, y, x) and byte strides of the real and half-complex arrays.
struct Layout {
  pocketfft::shape_t shape;
  pocketfft::stride_t sr, sc;
  pocketfft::shape_t axes{0, 1, 2};
};
Layout layout(const FftStam& s) {
  const std::size_t nx = s.pad[0], ny = s.pad[1], nz = s.pad[2], hx = nx / 2 + 1;
  Layout l;
  l.shape = {nz, ny, nx};
  l.sr = {static_cast<std::ptrdiff_t>(ny * nx * sizeof(float)), static_cast<std::ptrdiff_t>(nx * sizeof(float)),
          static_cast<std::ptrdiff_t>(sizeof(float))};
  l.sc = {static_cast<std::ptrdiff_t>(ny * hx * sizeof(std::complex<float>)),
          static_cast<std::ptrdiff_t>(hx * sizeof(std::complex<float>)),
          static_cast<std::ptrdiff_t>(sizeof(std::complex<float>))};
  return l;
}

}  // namespace

int fftSize(int n) {
  for (int m = std::max(n, 1);; ++m) {
    int r = m;
    for (int p : {2, 3, 5, 7})
      while (r % p == 0) r /= p;
    if (r == 1) return m;
  }
}

FftStam buildFftStam(const StamField& field, double clamp) {
  static std::atomic<std::uint64_t> next{1};
  FftStam s;
  s.id = next.fetch_add(1);
  s.geom = field.geom();
  s.clamp = clamp;
  for (int d = 0; d < 3; ++d) s.pad[d] = fftSize(s.geom.n[d]);
  const int nt = field.typeMapCount();
  const Layout l = layout(s);
  const float c = static_cast<float>(clamp);
  s.spectra.resize(static_cast<std::size_t>(nt) + 2);
  for (int m = 0; m < nt + 2; ++m) {
    const std::vector<float> v = field.pointMap(m);
    std::vector<float> in(s.realSize(), 0.0f);
    const GridGeometry& g = s.geom;
    for (int k = 0; k < g.n[2]; ++k)
      for (int j = 0; j < g.n[1]; ++j)
        for (int i = 0; i < g.n[0]; ++i) {
          const float x = v[g.index(i, j, k)];
          in[s.index(i, j, k)] = m == nt ? std::clamp(x, -c, c) : std::min(x, c);
        }
    s.spectra[m].resize(s.complexSize());
    pocketfft::r2c(l.shape, l.sr, l.sc, l.axes, pocketfft::FORWARD, in.data(), s.spectra[m].data(), 1.0f, 0);
  }
  return s;
}

void fftCorrelate(const FftStam& s, std::span<const Vec3> f, std::span<const int> map, std::span<const double> q,
                  FftWork& w) {
  fftProbeSpectra(s, f, map, q, w);
  fftCombine(s, w);
}

void fftProbeSpectra(const FftStam& s, std::span<const Vec3> f, std::span<const int> map, std::span<const double> q,
                     FftWork& w) {
  const Layout l = layout(s);
  // The energy buffer doubles as the input of the forward FFTs (zero outside the spread atoms).
  w.energy.assign(s.realSize(), 0.0f);
  w.channel.clear();
  const int nt = static_cast<int>(s.spectra.size()) - 2;
  // Base corner and trilinear weights per atom.
  struct Corner {
    std::array<int, 3> b;
    std::array<float, 3> u;
  };
  std::vector<Corner> cn(f.size());
  for (std::size_t a = 0; a < f.size(); ++a)
    for (int d = 0; d < 3; ++d) {
      const double fl = std::floor(f[a][d]);
      const int p = s.pad[d];
      cn[a].b[d] = ((static_cast<int>(fl) % p) + p) % p;
      cn[a].u[d] = static_cast<float>(f[a][d] - fl);
    }
  auto spread = [&](std::size_t a, float wt) {
    const Corner& c = cn[a];
    for (int dz = 0; dz < 2; ++dz) {
      const int k = (c.b[2] + dz) % s.pad[2];
      const float wz = dz ? c.u[2] : 1.0f - c.u[2];
      for (int dy = 0; dy < 2; ++dy) {
        const int j = (c.b[1] + dy) % s.pad[1];
        const float wy = dy ? c.u[1] : 1.0f - c.u[1];
        for (int dx = 0; dx < 2; ++dx) {
          const int i = (c.b[0] + dx) % s.pad[0];
          const float wx = dx ? c.u[0] : 1.0f - c.u[0];
          const std::size_t idx = s.index(i, j, k);
          w.energy[idx] += wt * wx * wy * wz;
          w.touched.push_back(idx);
        }
      }
    }
  };
  auto channel = [&](int m) {
    const std::size_t c = w.channel.size();
    w.channel.push_back(m);
    if (w.spec.size() <= c) w.spec.resize(c + 1);
    w.spec[c].resize(s.complexSize());
    pocketfft::r2c(l.shape, l.sr, l.sc, l.axes, pocketfft::FORWARD, w.energy.data(), w.spec[c].data(), 1.0f, 1);
    for (std::size_t idx : w.touched) w.energy[idx] = 0.0f;
    w.touched.clear();
  };
  for (int m = 0; m < nt; ++m) {
    bool any = false;
    for (std::size_t a = 0; a < f.size(); ++a)
      if (map[a] == m) {
        spread(a, 1.0f);
        any = true;
      }
    if (any) channel(m);
  }
  bool charged = false;
  for (std::size_t a = 0; a < f.size(); ++a)
    if (q[a] != 0) {
      spread(a, static_cast<float>(q[a]));
      charged = true;
    }
  if (charged) {
    channel(nt);
    for (std::size_t a = 0; a < f.size(); ++a)
      if (q[a] != 0) spread(a, static_cast<float>(std::abs(q[a])));
    channel(nt + 1);
  }
}

void fftCombine(const FftStam& s, FftWork& w) {
  if (w.channel.empty()) {  // no probe atoms
    w.energy.assign(s.realSize(), 0.0f);
    return;
  }
  w.energy.resize(s.realSize());
  const Layout l = layout(s);
  w.acc.resize(s.complexSize());
  for (std::size_t c = 0; c < w.channel.size(); ++c) {
    const auto& sp = s.spectra[w.channel[c]];
    const auto& pc = w.spec[c];
    if (c == 0)
      for (std::size_t k = 0; k < w.acc.size(); ++k) w.acc[k] = std::conj(pc[k]) * sp[k];
    else
      for (std::size_t k = 0; k < w.acc.size(); ++k) w.acc[k] += std::conj(pc[k]) * sp[k];
  }
  pocketfft::c2r(l.shape, l.sc, l.sr, l.axes, pocketfft::BACKWARD, w.acc.data(), w.energy.data(),
                 1.0f / static_cast<float>(s.realSize()), 1);
}

void fftSelect(const FftStam& s, const float* energy, const std::array<int, 3>& lo, const std::array<int, 3>& hi,
               std::size_t top, std::vector<std::pair<float, std::size_t>>& best) {
  best.clear();  // max-heap on (energy, index)
  for (int k = lo[2]; k <= hi[2]; ++k)
    for (int j = lo[1]; j <= hi[1]; ++j)
      for (int i = lo[0]; i <= hi[0]; ++i) {
        const std::size_t idx = s.index(i, j, k);
        const std::pair<float, std::size_t> e{energy[idx], idx};
        if (best.size() < top) {
          best.push_back(e);
          std::push_heap(best.begin(), best.end());
        } else if (e < best.front()) {
          std::pop_heap(best.begin(), best.end());
          best.back() = e;
          std::push_heap(best.begin(), best.end());
        }
      }
  std::sort_heap(best.begin(), best.end());
}

}  // namespace od
