# Tech Stack & Runtime Environment

## Language & Compiler
- Language Standard: C++20
- Supported Compilers: GCC 11+, Clang 14+ (Host: GCC 16.2.1)
- Build System: CMake 3.24+
- Build Flags: `-Wall -Wextra -Werror -O3 -fPIC` (Debug: `-g -O0`)

## External Dependencies
- Solver Core (`src/core/`, `src/simplex/`): Zero external dependencies. Self-contained sparse linear algebra.
- Testing: Custom test runner / lightweight assertion harness without mandatory external dependencies, compatible with ctest.

## Hardware & Architecture
- Platform: Linux x86_64
- Memory Model: 64-bit address space, 32-bit indices (`int32_t`) for cache efficiency, 64-bit double precision (`double`) for numerical stability.
