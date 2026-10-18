// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <format>
#include <string>

#include "error.hpp"

namespace od {

void logInfo(const std::string& msg);
void logWarn(const std::string& msg);
int warningCount();

template <class... Args>
[[noreturn]] void fail(std::format_string<Args...> fmt, Args&&... args) {
  throw Error(std::format(fmt, std::forward<Args>(args)...));
}

}  // namespace od
