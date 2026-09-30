cmake_minimum_required(VERSION 3.24)
project(suplex LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Strict warnings and optimization flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Werror -O3 -fPIC)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------
# Core Library
# ---------------------------------------------------------
add_library(suplex_core
    src/core/sparse_matrix.cpp
    src/core/problem.cpp
    src/core/solution.cpp
)
target_include_directories(suplex_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------
# Simplex Library (Person 1)
# ---------------------------------------------------------
add_library(suplex_simplex
    src/simplex/basis.cpp
    src/simplex/lu_factor.cpp
    src/simplex/lu_update.cpp
    src/simplex/pricing.cpp
    src/simplex/primal_simplex.cpp
    src/simplex/dual_simplex.cpp
)
target_link_libraries(suplex_simplex PUBLIC suplex_core)
target_include_directories(suplex_simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ---------------------------------------------------------
# Unit Tests
# ---------------------------------------------------------
enable_testing()

add_executable(suplex_test_runner
    tests/unit/test_runner.cpp
    tests/unit/test_sparse_matrix.cpp
    tests/unit/test_lu.cpp
    tests/unit/test_simplex.cpp
)
target_link_libraries(suplex_test_runner PRIVATE suplex_simplex suplex_core)
target_include_directories(suplex_test_runner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME suplex_unit_tests COMMAND suplex_test_runner)
