# OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
# Copyright (C) 2026 Gao Quan-Ze.
# SPDX-License-Identifier: LGPL-2.1-or-later

add_executable(od_tests test_main.cpp test_cli.cpp test_input.cpp)
target_include_directories(od_tests PRIVATE ${PROJECT_SOURCE_DIR}/third_party/gemmi/third_party)
target_compile_definitions(od_tests PRIVATE OD_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(od_tests PRIVATE od_core)
add_test(NAME od_tests COMMAND od_tests)
