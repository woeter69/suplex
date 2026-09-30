cmake_minimum_required(VERSION 3.24)
project(suplex_person3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ─────────────────────────────────────────────────────────────
add_compile_options(-Wall -Wextra -Wpedantic)
# Enable optimisations in release, debug info in debug
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0 -fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
else()
    add_compile_options(-O3 -march=native)
endif()

# ── Core stubs (shared by all targets) ────────────────────────────────────────
add_library(suplex_core STATIC
    src/core/core_stubs.cpp
)
target_include_directories(suplex_core PUBLIC src)

# ── MILP library ───────────────────────────────────────────────────────────────
add_library(suplex_milp STATIC
    src/milp/presolve.cpp
    src/milp/node.cpp
    src/milp/branching.cpp
    src/milp/cuts.cpp
    src/milp/heuristics.cpp
    src/milp/branch_bound.cpp
    src/milp/conflict.cpp
    src/milp/milp_solver.cpp
)
target_include_directories(suplex_milp PUBLIC src)
target_link_libraries(suplex_milp PUBLIC suplex_core)

# ── Unit tests ─────────────────────────────────────────────────────────────────
enable_testing()

# Fetch Google Test
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_presolve
    tests/unit/milp/test_presolve.cpp
)
target_link_libraries(test_presolve PRIVATE suplex_milp GTest::gtest_main)

add_executable(test_branching
    tests/unit/milp/test_branching.cpp
)
target_link_libraries(test_branching PRIVATE suplex_milp GTest::gtest_main)

add_executable(test_cuts
    tests/unit/milp/test_cuts.cpp
)
target_link_libraries(test_cuts PRIVATE suplex_milp GTest::gtest_main)

add_executable(test_heuristics
    tests/unit/milp/test_heuristics.cpp
)
target_link_libraries(test_heuristics PRIVATE suplex_milp GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_presolve)
gtest_discover_tests(test_branching)
gtest_discover_tests(test_cuts)
gtest_discover_tests(test_heuristics)
