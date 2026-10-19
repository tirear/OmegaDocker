# OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
# Copyright (C) 2026 Gao Quan-Ze.
# SPDX-License-Identifier: LGPL-2.1-or-later

cmake_minimum_required(VERSION 3.20)
project(OmegaDocker VERSION 0.40.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(OD_STATIC_DEPS "Link third-party libraries statically" OFF)
option(OD_BUILD_TESTS "Build tests" ON)
option(OD_CUDA "Build the CUDA back end" OFF)

# ---------- dependencies ----------
if(OD_STATIC_DEPS)
  # Must be set before find_package(RDKit), which calls find_package(Boost) internally.
  set(Boost_USE_STATIC_LIBS ON)
  set(OD_RDKIT_SUFFIX "_static")
else()
  set(OD_RDKIT_SUFFIX "")
endif()
find_package(RDKit REQUIRED CONFIG)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(OD_RDKIT_LIBS)
foreach(c PartialCharges FileParsers SmilesParse SubstructMatch GraphMol RDGeneral)
  list(APPEND OD_RDKIT_LIBS RDKit::${c}${OD_RDKIT_SUFFIX})
endforeach()

# gemmi (compiled from third_party/gemmi; gz support not used)
file(GLOB OD_GEMMI_SRC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/src/*.cpp)
list(FILTER OD_GEMMI_SRC EXCLUDE REGEX "/(gz|mmread_gz)\\.cpp$")
add_library(od_gemmi STATIC ${OD_GEMMI_SRC})
target_include_directories(od_gemmi
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gemmi/third_party)
target_compile_definitions(od_gemmi PRIVATE GEMMI_BUILD)
set_target_properties(od_gemmi PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(od_tomlpp INTERFACE)
target_include_directories(od_tomlpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/tomlplusplus)

# ---------- embedded data ----------
file(GLOB OD_CCD_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/ccd/*.cif)
list(SORT OD_CCD_FILES)
set(OD_CCD_EMBED ${CMAKE_CURRENT_BINARY_DIR}/generated/ccd_embedded.cpp)
set(OD_CCD_LIST ${CMAKE_CURRENT_BINARY_DIR}/generated/ccd_files.txt)
file(WRITE ${OD_CCD_LIST} "${OD_CCD_FILES}")
add_custom_command(
  OUTPUT ${OD_CCD_EMBED}
  COMMAND ${CMAKE_COMMAND} -DOUT=${OD_CCD_EMBED} -DLIST_FILE=${OD_CCD_LIST}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_ccd.cmake
  DEPENDS ${OD_CCD_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_ccd.cmake
  COMMENT "Embedding CCD subset")

set(OD_RAMA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/top8000)
file(GLOB OD_RAMA_FILES CONFIGURE_DEPENDS ${OD_RAMA_DIR}/rama8000-*.data)
set(OD_RAMA_EMBED ${CMAKE_CURRENT_BINARY_DIR}/generated/rama_embedded.cpp)
configure_file(cmake/rama_embedded.cpp.in ${OD_RAMA_EMBED} @ONLY)
set_source_files_properties(${OD_RAMA_EMBED} PROPERTIES OBJECT_DEPENDS "${OD_RAMA_FILES}")

configure_file(src/version.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/generated/version.hpp @ONLY)

# ---------- core library ----------
add_library(od_core STATIC
  src/log.cpp
  src/cli.cpp
  src/molecule.cpp
  src/ccd.cpp
  src/bond_orders.cpp
  src/chem.cpp
  src/read_mmcif.cpp
  src/read_rdkit.cpp
  src/read_input.cpp
  src/topology.cpp
  src/dssp.cpp
  src/molconfig.cpp
  src/rama.cpp
  src/torsion_tree.cpp
  src/ad4.cpp
  src/energy.cpp
  src/csgm1.cpp
  src/csgm2.cpp
  src/csgm3.cpp
  src/minimize.cpp
  src/fiber_batch.cpp
  src/conf9a.cpp
  src/grow.cpp
  src/grid.cpp
  src/shaking.cpp
  src/sysmem.cpp
  src/score9d.cpp
  src/write_struct.cpp
  src/orient.cpp
  src/dock9b.cpp
  src/fft_dock.cpp
  src/dock9c.cpp
  src/symrmsd.cpp
  src/checkpoint.cpp
  src/closure.cpp
  ${OD_CCD_EMBED}
  ${OD_RAMA_EMBED})
target_include_directories(od_core PUBLIC src ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_include_directories(od_core SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pocketfft)
target_link_libraries(od_core PUBLIC od_gemmi od_tomlpp ${OD_RDKIT_LIBS} Eigen3::Eigen)
if(OD_STATIC_DEPS)
  target_compile_options(od_core PUBLIC ${OpenMP_CXX_FLAGS})
  target_link_libraries(od_core PUBLIC -l:libgomp.a Threads::Threads ${CMAKE_DL_LIBS})
else()
  target_link_libraries(od_core PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(od_core PRIVATE -Wall -Wextra)

# ---------- CUDA back end ----------
if(OD_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70)
  endif()
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  add_library(od_gpu STATIC src/gpu/energy_gpu.cu src/gpu/prep_gpu.cu src/gpu/grid_gpu.cu src/gpu/shake_gpu.cu)
  # The quick filter and the energy grids follow the CPU operations exactly (spec 15.1.5, 15.1.10).
  set_source_files_properties(src/gpu/prep_gpu.cu src/gpu/grid_gpu.cu PROPERTIES COMPILE_OPTIONS "-fmad=false")
  target_include_directories(od_gpu PUBLIC src src/gpu)
  target_compile_definitions(od_gpu PUBLIC OD_HAVE_CUDA)
  if(OD_STATIC_DEPS)
    target_link_libraries(od_gpu PUBLIC CUDA::cudart_static CUDA::cufft_static_nocallback)
  else()
    target_link_libraries(od_gpu PUBLIC CUDA::cudart CUDA::cufft)
  endif()
  target_link_libraries(od_core PUBLIC od_gpu)
endif()

add_executable(omegadocker src/main.cpp)
target_link_libraries(omegadocker PRIVATE od_core)
if(OD_STATIC_DEPS)
  target_link_options(omegadocker PRIVATE -static-libstdc++ -static-libgcc)
endif()

if(OD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
