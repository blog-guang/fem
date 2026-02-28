cmake_minimum_required(VERSION 3.20)
project(fem_multiphysics VERSION 0.1.0 LANGUAGES CXX)

# ── C++17 ──
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── 默认 Release ──
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ── 编译选项 ──
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0 -Wall -Wextra -Wpedantic)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O2 -DNDEBUG)
endif()

# ── 输出目录 ──
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ── AMGCL (header-only library) ──
set(AMGCL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/third_party/amgcl)
include_directories(${AMGCL_INCLUDE_DIR})

# ── 子目录 ──
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(examples)

# ── 测试 ──
enable_testing()
