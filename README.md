# Suplex

**Sovereign Unified Platform for Linear and Extended Optimization**

Suplex is a high-performance mathematical optimization solver written from scratch in C++20. It supports Linear Programming (LP), Mixed-Integer Linear Programming (MILP), and Quadratic Programming (QP). Every numerical routine -- sparse linear algebra, LU factorization, simplex algorithms, branch-and-bound, cutting planes, and presolve -- is implemented in-house with zero external solver dependencies.

Built for the Smart India Hackathon.

Version: 0.1.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [Features](#2-features)
3. [Architecture](#3-architecture)
4. [Prerequisites](#4-prerequisites)
5. [Building from Source](#5-building-from-source)
6. [Command-Line Interface](#6-command-line-interface)
7. [Input File Formats](#7-input-file-formats)
8. [C++ API](#8-c-api)
9. [C API](#9-c-api)
10. [Python API](#10-python-api)
11. [Solver Algorithms](#11-solver-algorithms)
12. [Configuration Reference](#12-configuration-reference)
13. [Testing](#13-testing)
14. [Benchmarking](#14-benchmarking)
15. [Project Structure](#15-project-structure)
16. [Design Decisions](#16-design-decisions)
17. [Known Limitations](#17-known-limitations)
18. [Team Contributions](#18-team-contributions)
19. [License](#19-license)

---

## 1. Overview

Suplex solves optimization problems of the following general form:

```
minimize (or maximize)   c^T x + 0.5 x^T Q x

subject to               row_lower <= A x <= row_upper
                          col_lower <=   x <= col_upper
                          x_j in {integers}   for j in J
```

Where:

- `A` is an m-by-n sparse constraint matrix stored in Compressed Sparse Column (CSC) format.
- `c` is the n-vector of linear objective coefficients.
- `Q` is an optional n-by-n symmetric positive semidefinite matrix for quadratic objectives.
- `row_lower` and `row_upper` are m-vectors defining ranged row constraints.
- `col_lower` and `col_upper` are n-vectors defining variable bounds.
- `J` is the subset of variable indices constrained to integer values.

The solver pipeline is: parse input, presolve (reduce), solve (simplex/IPM/branch-and-bound), postsolve (map solution back to original space), write output.

---

## 2. Features

### Solver Capabilities

- Revised Primal Simplex with bounded variables.
- Revised Dual Simplex with warm-start support for branch-and-bound integration.
- Full branch-and-bound framework for mixed-integer programs.
- Presolve engine with nine reduction rules for problem size reduction.
- Four cutting plane generators: Gomory, MIR, knapsack cover, and clique cuts.
- Five primal heuristics: simple rounding, fractional diving, coefficient diving, feasibility pump, and RINS.
- Conflict analysis with no-good cut generation from infeasible nodes.
- GPU-accelerated sparse matrix-vector multiplication and preconditioned conjugate gradient solver (optional, CUDA).

### Sparse Linear Algebra

- Compressed Sparse Column (CSC) and Compressed Sparse Row (CSR) matrix formats.
- Triplet (COO) format for incremental matrix construction.
- Sparse LU factorization with Markowitz ordering and threshold pivoting.
- LU basis update via Forrest-Tomlin eta transformations.
- Iterative refinement and automatic refactorization heuristics.

### Input/Output

- MPS file parser (fixed and free format, OBJSENSE, RANGES, INTORG/INTEND).
- CPLEX LP file parser.
- Solution file writer (primal values, dual values, reduced costs, row activities).
- Named rows and columns preserved from input files.

### API Surface

- C++20 high-level solver class (`suplex::Suplex`).
- C11-compatible ABI for embedding in C programs and foreign function interfaces.
- Python bindings via pybind11 with NumPy array support for solution vectors.

### Build System

- CMake 3.24+ with modular library targets.
- GNU Make fallback for direct compilation without CMake.
- All code compiles with `-Wall -Wextra -Werror` under GCC and Clang.

---

## 3. Architecture

Suplex follows a layered architecture with strict separation between data representation, numerical computation, optimization algorithms, and external interfaces.

```
Layer 4: Interface       -- Parsers, CLI, C/Python API, GPU dispatch
Layer 3: Orchestration   -- Presolve, postsolve, solver pipeline
Layer 2: Algorithms      -- Simplex, IPM, Branch-and-Bound
Layer 1: Numerical       -- Sparse LU factorization, ordering
Layer 0: Core            -- Problem, Solution, SparseMatrixCSC/CSR
```

### Data Flow

```
Input (MPS/LP/API)
    |
    v
Problem construction (ranged-row representation)
    |
    v
Presolve (reduce: empty rows, fixed vars, singletons, forcing rows, implied bounds)
    |
    +-- LP:   Simplex or IPM --> Postsolve --> Solution
    |
    +-- MILP: Root LP relaxation --> Cuts --> Branch-and-Bound --> Postsolve --> Solution
    |
    +-- QP:   IPM-QP --> Postsolve --> Solution
```

### Library Decomposition

The build produces six static libraries and one shared library:

| Library | Contents |
|---|---|
| `suplex_core` | `SparseMatrixCSC`, `SparseMatrixCSR`, `TripletMatrix`, `Problem`, `Solution`, type definitions |
| `suplex_simplex` | `PrimalSimplex`, `DualSimplex`, `Basis`, `LUFactor`, `LUUpdate`, `Pricing` |
| `suplex_milp` | `Presolve`, `BranchBound`, `Branching`, `CutGenerator`, `CutPool`, `Heuristics`, `ConflictAnalysis`, `MILPSolver` |
| `suplex_io` | `MPSParser`, `LPParser`, `SolutionWriter` |
| `suplex_gpu` | `GPUManager`, CUDA SpMV kernel, CUDA PCG kernel (CPU fallback when CUDA is unavailable) |
| `suplex` (shared) | `Suplex` high-level solver, C API (`suplex_api.cpp`) |

Dependency chain: `suplex` depends on `suplex_io`, `suplex_milp`, `suplex_simplex`, `suplex_gpu`, and `suplex_core`.

---

## 4. Prerequisites

### Required

- C++20-capable compiler (GCC 12+ or Clang 15+ recommended).
- CMake 3.24 or later.
- GoogleTest (automatically discovered by CMake for Person 3 unit tests; the core test suite does not require it).

### Optional

- CUDA Toolkit 12+ and a CUDA-capable GPU (for GPU acceleration).
- pybind11 (for Python bindings).
- Python 3.8+ with NumPy (for the Python API and its tests).

---

## 5. Building from Source

### Standard CMake Build

```bash
git clone https://github.com/woeter69/suplex.git
cd suplex
cmake -S . -B build \
  -DSUPLEX_BUILD_TESTS=ON \
  -DSUPLEX_BUILD_BENCHMARKS=ON
cmake --build build -- -j$(nproc)
```

The build produces:

- `build/libsuplex.so` -- shared library.
- `build/suplex` -- command-line executable.
- `build/suplex_test_runner` -- core unit and integration test binary.
- `build/suplex_benchmark` -- benchmark runner.
- `build/suplex_p3_*_test` -- Person 3 GoogleTest binaries (when GoogleTest is found).
- `build/suplex_c_api_test` -- C API integration test.

### CMake Options

| Option | Default | Description |
|---|---|---|
| `SUPLEX_BUILD_TESTS` | `ON` | Build the test suite and register CTest tests. |
| `SUPLEX_BUILD_BENCHMARKS` | `ON` | Build the benchmark runner. |
| `SUPLEX_BUILD_PYTHON` | `OFF` | Build Python bindings (requires pybind11). |
| `SUPLEX_ENABLE_GPU` | `OFF` | Enable CUDA compilation for GPU kernels. |
| `SUPLEX_BUILD_P3_GTESTS` | `ON` | Build Person 3 GoogleTest suites when GoogleTest is available. |

### Build with GPU Support

```bash
cmake -S . -B build-gpu \
  -DSUPLEX_ENABLE_GPU=ON \
  -DSUPLEX_BUILD_TESTS=ON
cmake --build build-gpu -- -j$(nproc)
```

### Build with Python Bindings

```bash
cmake -S . -B build-python \
  -DSUPLEX_BUILD_PYTHON=ON
cmake --build build-python -- -j$(nproc)
```

This produces `pysuplex.cpython-*.so` in the build directory.

### GNU Make Fallback

If CMake is unavailable, a direct `Makefile` is provided:

```bash
make -j$(nproc)
make test
make benchmark
```

### Manual GCC Compilation

For environments with only a compiler and no build system:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Werror -I. \
  -DSUPLEX_SOURCE_DIR=\"$PWD\" \
  $(find src/core src/simplex src/milp src/io src/api src/gpu \
    -name '*.cpp' ! -path '*/python/*') \
  tests/unit/test_runner.cpp \
  tests/unit/test_sparse_matrix.cpp \
  tests/unit/test_lu.cpp \
  tests/unit/test_simplex.cpp \
  tests/unit/test_boundary.cpp \
  tests/integration/test_integration.cpp \
  -o suplex_all_tests
./suplex_all_tests
```

---

## 6. Command-Line Interface

### Synopsis

```
suplex [options] <problem_file>
```

The CLI accepts `.mps`, `.qps`, and `.lp` files. File format is determined by extension.

### Options

| Flag | Argument | Default | Description |
|---|---|---|---|
| `-h`, `--help` | -- | -- | Print usage information and exit. |
| `-v`, `--version` | -- | -- | Print version string and exit. |
| `--algorithm` | `auto`, `primal`, `dual`, `ipm` | `auto` | Select the LP algorithm. `auto` chooses based on problem structure. |
| `--presolve` | `on`, `off` | `on` | Enable or disable the presolve phase. |
| `--time-limit` | positive float (seconds) | `3600.0` | Maximum wall-clock time for the solve. |
| `--threads` | nonnegative integer | `0` (automatic) | Number of worker threads. `0` uses hardware concurrency. |
| `--gpu` | `on`, `off` | `off` | Request GPU acceleration for supported operations. |
| `--log-level` | `trace`, `debug`, `info`, `warn`, `error`, `off` | `info` | Verbosity of solver output. |
| `--mip-gap` | nonnegative float | `1e-4` | Relative MIP optimality gap tolerance. |
| `--solution` | file path | -- | Write a solution report to the specified file. |
| `--stats` | -- | -- | Print detailed solve statistics (iterations, time, nodes, gap). |

### Examples

Solve an LP file with the dual simplex algorithm:

```bash
./build/suplex --algorithm dual data/examples/tiny.lp
```

Solve an MPS file and write the solution:

```bash
./build/suplex --solution output.sol --stats data/examples/tiny.mps
```

Solve a MIP with a tight gap tolerance and time limit:

```bash
./build/suplex --mip-gap 0.001 --time-limit 60 data/examples/tiny_mip.lp
```

### Exit Codes

| Code | Meaning |
|---|---|
| `0` | Optimal solution found. |
| `1` | Solve completed with a non-optimal status (infeasible, unbounded, limit reached, etc.). |
| `2` | Input error (missing file, unsupported format, invalid options). |

### Output Format

The CLI prints to standard output:

```
Problem: 2 rows, 2 columns, 4 nonzeros, 0 integer variables
Status: OPTIMAL
Objective: 3
```

With `--stats` enabled, additional lines follow:

```
Iterations: 2
Solve time: 0.000123 seconds
```

For MIP problems, `--stats` also prints:

```
Nodes: 5
MIP gap: 0
```

---

## 7. Input File Formats

### MPS Format

Suplex supports the standard MPS file format including the following sections and features:

- `NAME` -- problem name.
- `OBJSENSE` -- optional section specifying `MIN` or `MAX` (defaults to `MIN`).
- `ROWS` -- row declarations with type indicators: `N` (objective), `G` (>=), `L` (<=), `E` (=).
- `COLUMNS` -- nonzero coefficients in column-major order. Supports `INTORG`/`INTEND` markers for integer variable declaration.
- `RHS` -- right-hand side values.
- `RANGES` -- range values for ranged constraints.
- `BOUNDS` -- variable bound declarations: `LO`, `UP`, `FX` (fixed), `FR` (free), `MI` (minus infinity lower), `PL` (plus infinity upper), `BV` (binary), `LI` (integer lower), `UI` (integer upper).
- `ENDATA` -- required terminator.

Both fixed-format (columns at specific character positions) and free-format MPS are supported. Lines beginning with `*` are treated as comments.

Example (`data/examples/tiny.mps`):

```
NAME          TINY
OBJSENSE
 MIN
ROWS
 N  COST
 G  DEMAND
 L  CAPACITY
COLUMNS
    X1        COST        1          DEMAND      1
    X1        CAPACITY    1
    X2        COST        2          DEMAND      1
    X2        CAPACITY    1
RHS
    RHS1      DEMAND      3          CAPACITY    5
BOUNDS
 UP BND1      X1          4
 UP BND1      X2          4
ENDATA
```

### LP Format

Suplex supports the CPLEX LP file format:

- `Minimize` or `Maximize` keyword followed by the objective function.
- `Subject To` section with named constraints using `<=`, `>=`, and `=` operators.
- `Bounds` section with variable bound declarations.
- `Binary`, `General`, and `Integer` sections for integer variable declarations.
- `End` keyword as terminator.

The parser handles coefficients written with or without spaces around operators, and supports the `*` multiplication operator.

Example (`data/examples/tiny.lp`):

```
Minimize
 obj: x1 + 2 x2
Subject To
 demand: x1 + x2 >= 3
 capacity: x1 + x2 <= 5
Bounds
 0 <= x1 <= 4
 0 <= x2 <= 4
End
```

### Included Example Files

| File | Description |
|---|---|
| `data/examples/tiny.mps` | Minimal 2-variable, 2-constraint LP in MPS format. |
| `data/examples/tiny.lp` | Same problem in LP format. |
| `data/examples/tiny_mip.lp` | Small binary MIP for branch-and-bound testing. |
| `data/examples/features.mps` | MPS file exercising OBJSENSE, RANGES, free bounds, binary bounds, INTORG/INTEND. |
| `data/examples/compact_bounds.lp` | LP file with compact operator syntax and general integer variables. |
| `data/examples/branch_mip.lp` | MIP designed to require branching (fractional root relaxation). |
| `data/examples/invalid_missing_end.mps` | Intentionally malformed MPS file for parser error testing. |

---

## 8. C++ API

The primary C++ interface is the `suplex::Suplex` class defined in `src/api/suplex_solver.h`.

### Basic Usage

```cpp
#include "src/api/suplex_solver.h"
#include <iostream>

int main() {
    suplex::Suplex solver;

    // Load a problem from file
    if (!solver.read_mps("data/examples/tiny.mps")) {
        std::cerr << solver.last_error() << std::endl;
        return 1;
    }

    // Configure
    solver.set_algorithm(suplex::Algorithm::DUAL_SIMPLEX);
    solver.set_presolve(true);
    solver.set_time_limit(60.0);

    // Solve
    suplex::SolverStatus status = solver.solve();

    // Access results
    const suplex::Solution& sol = solver.solution();
    std::cout << "Status: " << suplex::to_string(status) << std::endl;
    if (sol.is_feasible()) {
        std::cout << "Objective: " << sol.objective_value << std::endl;
        for (size_t i = 0; i < sol.primal_values.size(); ++i) {
            std::cout << "x[" << i << "] = " << sol.primal_values[i] << std::endl;
        }
    }
    return 0;
}
```

### Programmatic Problem Construction

```cpp
#include "src/api/suplex_solver.h"

int main() {
    suplex::Problem problem;
    problem.set_dimensions(2, 2);  // 2 rows, 2 columns
    problem.set_objective_sense(suplex::ObjectiveSense::MINIMIZE);
    problem.set_objective({1.0, 2.0});

    // Build constraint matrix from triplets
    std::vector<suplex::Index> rows = {0, 0, 1, 1};
    std::vector<suplex::Index> cols = {0, 1, 0, 1};
    std::vector<suplex::Real>  vals = {1.0, 1.0, 1.0, 1.0};
    problem.set_constraint_matrix(
        suplex::SparseMatrixCSC::from_triplets(2, 2, rows, cols, vals));

    problem.set_row_bounds({3.0, -suplex::INF}, {suplex::INF, 5.0});
    problem.set_col_bounds({0.0, 0.0}, {4.0, 4.0});

    suplex::Suplex solver;
    solver.set_problem(std::move(problem));
    solver.solve();

    const suplex::Solution& sol = solver.solution();
    // sol.objective_value == 3.0
    return 0;
}
```

### Suplex Class Reference

```cpp
namespace suplex {

enum class Algorithm { AUTO, PRIMAL_SIMPLEX, DUAL_SIMPLEX, IPM };

class Suplex {
public:
    Suplex();

    // Problem loading
    bool read_mps(const std::string& filename);
    bool read_lp(const std::string& filename);
    void set_problem(Problem&& problem);
    const Problem& problem() const;

    // Configuration
    void set_algorithm(Algorithm algo);
    void set_presolve(bool enabled);
    void set_gpu(bool enabled);
    void set_log_level(LogLevel level);
    void set_time_limit(Real seconds);
    void set_threads(int num_threads);
    void set_mip_gap(Real gap);

    // Solve
    SolverStatus solve();

    // Results
    const Solution& solution() const;
    const std::vector<std::string>& column_names() const;
    const std::vector<std::string>& row_names() const;
    std::string last_error() const;
};

}  // namespace suplex
```

### Solution Structure

```cpp
struct Solution {
    SolverStatus status;
    Real objective_value;
    std::vector<Real> primal_values;      // decision variable values
    std::vector<Real> dual_values;        // shadow prices / row duals
    std::vector<Real> reduced_costs;      // reduced costs
    std::vector<Real> row_activities;     // constraint activity values (Ax)
    Real best_bound;                      // best dual bound (MILP)
    Real mip_gap;                         // relative MIP gap (MILP)
    int64_t num_nodes;                    // branch-and-bound nodes (MILP)
    int64_t num_iterations;
    double solve_time_seconds;

    bool is_optimal() const;
    bool is_feasible() const;
};
```

### Solver Status Values

| Status | Meaning |
|---|---|
| `NOT_STARTED` | Solve has not been called. |
| `OPTIMAL` | An optimal solution was found within all tolerances. |
| `INFEASIBLE` | The problem has no feasible solution. |
| `UNBOUNDED` | The objective is unbounded. |
| `INF_OR_UNBD` | The problem is either infeasible or unbounded (ambiguous). |
| `ITERATION_LIMIT` | The iteration limit was reached before convergence. |
| `TIME_LIMIT` | The time limit was reached before convergence. |
| `NUMERICAL_ERROR` | A numerical issue prevented the solve from completing. |
| `USER_INTERRUPT` | The solve was interrupted by a user callback. |

---

## 9. C API

The C API provides an `extern "C"` interface suitable for embedding Suplex in C programs, or for calling from languages with C FFI support (Rust, Go, Fortran, etc.). The header is `src/api/suplex.h`.

### Usage from C

```c
#include "suplex.h"
#include <stdio.h>

int main(void) {
    SuplexSolver* solver = suplex_create();
    if (!solver) return 1;

    int rc = suplex_read_mps(solver, "data/examples/tiny.mps");
    if (rc != SUPLEX_OK) {
        printf("Error: %s\n", suplex_get_error(solver));
        suplex_destroy(solver);
        return 1;
    }

    suplex_set_algorithm(solver, SUPLEX_ALGORITHM_DUAL);
    int status = suplex_solve(solver);

    if (status == SUPLEX_STATUS_OPTIMAL) {
        printf("Objective: %f\n", suplex_get_objective(solver));
        double x[2];
        suplex_get_solution(solver, x, 2);
        printf("x[0] = %f, x[1] = %f\n", x[0], x[1]);
    }

    suplex_destroy(solver);
    return 0;
}
```

### C API Function Reference

**Lifecycle:**

| Function | Description |
|---|---|
| `SuplexSolver* suplex_create(void)` | Allocate a solver handle. Returns `NULL` on failure. |
| `void suplex_destroy(SuplexSolver* solver)` | Release a solver handle. `NULL` is safe. |

**Problem Loading:**

| Function | Description |
|---|---|
| `int suplex_read_mps(solver, filename)` | Load an MPS model file. |
| `int suplex_read_lp(solver, filename)` | Load a CPLEX LP model file. |
| `int suplex_set_dimensions(solver, rows, cols)` | Initialize a programmatic model with the given dimensions. |
| `int suplex_set_objective(solver, coeffs, sense)` | Set objective coefficients and sense (0=minimize, 1=maximize). |
| `int suplex_set_constraint_matrix(solver, nnz, rows, cols, vals)` | Set the constraint matrix in coordinate (triplet) format. |
| `int suplex_set_row_bounds(solver, lower, upper)` | Set constraint bound vectors. |
| `int suplex_set_col_bounds(solver, lower, upper)` | Set variable bound vectors. |
| `int suplex_set_var_types(solver, types)` | Set variable types (0=continuous, 1=integer, 2=binary). |

**Configuration:**

| Function | Description |
|---|---|
| `int suplex_set_algorithm(solver, algorithm)` | Set the LP algorithm (0=auto, 1=primal, 2=dual, 3=ipm). |
| `int suplex_set_presolve(solver, enabled)` | Enable (1) or disable (0) presolve. |
| `int suplex_set_time_limit(solver, seconds)` | Set the solve time limit. |
| `int suplex_set_log_level(solver, level)` | Set verbosity (0=trace through 5=off). |
| `int suplex_set_threads(solver, threads)` | Set thread count. |
| `int suplex_set_gpu(solver, enabled)` | Enable (1) or disable (0) GPU acceleration. |
| `int suplex_set_mip_gap(solver, gap)` | Set the relative MIP gap tolerance. |

**Solving and Results:**

| Function | Description |
|---|---|
| `int suplex_solve(solver)` | Solve and return a `SuplexStatus` value, or a negative error code. |
| `int suplex_get_status(solver)` | Return the current `SuplexStatus`. |
| `double suplex_get_objective(solver)` | Return the objective value. |
| `int suplex_get_solution(solver, values, size)` | Copy primal values into the provided buffer. Returns count copied, or `SUPLEX_ERROR_BUFFER_TOO_SMALL`. |
| `int suplex_get_dual_values(solver, values, size)` | Copy dual values. |
| `int suplex_get_reduced_costs(solver, values, size)` | Copy reduced costs. |
| `double suplex_get_solve_time(solver)` | Return solve time in seconds. |
| `int64_t suplex_get_iterations(solver)` | Return iteration count. |
| `int64_t suplex_get_nodes(solver)` | Return branch-and-bound node count. |
| `double suplex_get_mip_gap(solver)` | Return the relative MIP gap. |
| `const char* suplex_get_error(solver)` | Return the last error message. |
| `const char* suplex_version(void)` | Return the version string. |

**Return Codes:**

| Code | Name | Meaning |
|---|---|---|
| `0` | `SUPLEX_OK` | Success. |
| `-1` | `SUPLEX_ERROR_INVALID_ARGUMENT` | A `NULL` handle or invalid parameter. |
| `-2` | `SUPLEX_ERROR_IO` | File I/O failure. |
| `-3` | `SUPLEX_ERROR_INVALID_PROBLEM` | Problem data is dimensionally inconsistent. |
| `-4` | `SUPLEX_ERROR_BUFFER_TOO_SMALL` | The output buffer is smaller than the solution vector. |

---

## 10. Python API

Python bindings are provided via pybind11. The module is named `pysuplex`.

### Installation

Build with `SUPLEX_BUILD_PYTHON=ON` and add the build directory to `PYTHONPATH`:

```bash
cmake -S . -B build -DSUPLEX_BUILD_PYTHON=ON
cmake --build build -- -j$(nproc)
export PYTHONPATH=$PWD/build:$PYTHONPATH
```

### Usage

```python
import pysuplex

solver = pysuplex.Solver()
solver.read_lp("data/examples/tiny.lp")
solver.set_algorithm(pysuplex.Algorithm.DUAL_SIMPLEX)
solver.set_presolve(True)
solver.set_time_limit(60.0)

status = solver.solve()

print(f"Status: {solver.status}")
print(f"Objective: {solver.objective}")
print(f"Solution: {solver.solution}")       # NumPy array
print(f"Dual values: {solver.dual_values}")  # NumPy array
print(f"Reduced costs: {solver.reduced_costs}")
print(f"Solve time: {solver.solve_time} seconds")
print(f"Iterations: {solver.iterations}")
```

### Python API Reference

**`pysuplex.Solver` class:**

| Method / Property | Description |
|---|---|
| `read_mps(filename)` | Load MPS file. Raises `RuntimeError` on failure. |
| `read_lp(filename)` | Load LP file. Raises `RuntimeError` on failure. |
| `set_algorithm(algo)` | Set algorithm (`pysuplex.Algorithm.AUTO`, `PRIMAL_SIMPLEX`, `DUAL_SIMPLEX`, `IPM`). |
| `set_presolve(bool)` | Enable or disable presolve. |
| `set_gpu(bool)` | Enable or disable GPU acceleration. |
| `set_time_limit(seconds)` | Set the time limit. |
| `set_threads(n)` | Set the thread count. |
| `set_mip_gap(gap)` | Set the relative MIP gap tolerance. |
| `solve()` | Solve the loaded problem. Returns a `pysuplex.Status` enum. Releases the GIL during solve. |
| `status` (property) | Current `pysuplex.Status` value. |
| `objective` (property) | Objective value (float). |
| `solution` (property) | Primal solution as a NumPy `float64` array. |
| `dual_values` (property) | Dual solution as a NumPy `float64` array. |
| `reduced_costs` (property) | Reduced costs as a NumPy `float64` array. |
| `solve_time` (property) | Solve time in seconds (float). |
| `iterations` (property) | Iteration count (int). |
| `nodes` (property) | Branch-and-bound node count (int). |
| `mip_gap` (property) | Relative MIP gap (float). |

**`pysuplex.Status` enum:**

`NOT_STARTED`, `OPTIMAL`, `INFEASIBLE`, `UNBOUNDED`, `INF_OR_UNBD`, `ITERATION_LIMIT`, `TIME_LIMIT`, `NUMERICAL_ERROR`, `USER_INTERRUPT`.

**`pysuplex.Algorithm` enum:**

`AUTO`, `PRIMAL_SIMPLEX`, `DUAL_SIMPLEX`, `IPM`.

---

## 11. Solver Algorithms

### 11.1. Revised Primal Simplex

The primal simplex algorithm operates on the bounded-variable revised simplex formulation. It maintains a basis matrix `B` (a nonsingular m-by-m submatrix of `A`) and iterates by:

1. Computing the reduced cost vector `c_bar = c_N - c_B^T B^{-1} N`.
2. Selecting an entering variable using a pricing strategy (Dantzig, steepest edge, Devex, or Bland's anti-cycling rule).
3. Computing the pivot column `d = B^{-1} A_j`.
4. Performing a ratio test to find the leaving variable.
5. Updating the basis via an LU eta update.

The solver handles bounded variables through bound-flipping and applies iterative refinement when numerical accuracy degrades.

### 11.2. Revised Dual Simplex

The dual simplex algorithm is particularly efficient for re-optimization after bound changes, making it the default choice for branch-and-bound LP relaxations. It maintains dual feasibility and iterates by:

1. Selecting a leaving variable (the most infeasible basic variable).
2. Computing the pivot row `rho = e_i^T B^{-1}`.
3. Performing a dual ratio test to select the entering variable.
4. Updating the basis.

The `solve_from_basis()` method accepts a previously computed basis, enabling warm-start solves where only one or two pivots may be needed after a bound change.

### 11.3. Sparse LU Factorization

Basis matrices are factored as `B = L U` using a sparse LU factorization with:

- **Markowitz ordering**: selects pivots to minimize fill-in by choosing the element with the smallest Markowitz count (row_count - 1) * (col_count - 1).
- **Threshold pivoting**: rejects pivots smaller than `EPS_PIVOT` (1e-10) to maintain numerical stability.
- **Iterative refinement**: when the factorization error exceeds tolerance, the solve is refined.
- **Automatic refactorization**: after a configurable number of eta updates (or when numerical stability degrades), the basis is refactored from scratch.

### 11.4. Branch-and-Bound

The MILP solver uses a standard branch-and-bound framework:

1. Solve the root LP relaxation.
2. Apply cutting planes at the root and periodically during the tree search.
3. Run primal heuristics to find integer-feasible solutions.
4. Branch on fractional integer variables, creating child nodes with tightened bounds.
5. Prune nodes whose bound is worse than the incumbent.
6. Terminate when the gap between the best incumbent and the best bound is below tolerance.

**Node selection strategies:**
- `BEST_FIRST`: always process the node with the best bound (default for gap closure).
- `DEPTH_FIRST`: process deepest nodes first (good for finding feasible solutions quickly).
- `HYBRID`: plunge depth-first for a configurable number of nodes, then switch to best-first.

**Branching variable selection:**
- Most fractional: select the integer variable with fractional value closest to 0.5.
- Pseudo-cost branching: use historical branching data to estimate bound improvement.
- Strong branching: solve two LP relaxations (up and down) for candidate variables and select the one with the best score. Used for the top-K most promising candidates.
- Reliability branching (default): use strong branching until a variable has been observed a configurable number of times, then switch to pseudo-cost branching for that variable.

### 11.5. Cutting Planes

Four families of cutting planes are implemented:

- **Gomory mixed-integer cuts**: derived from the simplex tableau for basic integer variables with fractional values. Requires tableau row extraction from the LP solver.
- **Mixed-integer rounding (MIR) cuts**: applied to each constraint row. The row `sum(a_ij x_j) <= b` is transformed using the MIR formula to derive a valid inequality.
- **Knapsack cover cuts**: for constraints involving binary variables with positive coefficients. A greedy algorithm finds a minimal cover, and non-cover variables are lifted using sequential lifting.
- **Clique cuts**: constraints of the form `sum(x_j) <= 1` where all variables are binary are identified as clique constraints and added directly.

Cut selection filters by efficacy (minimum violation threshold) and orthogonality (cuts that are more than 90% parallel to already-selected cuts are rejected).

### 11.6. Primal Heuristics

- **Simple rounding**: round each fractional integer variable to the nearest integer; check feasibility of the rounded point.
- **Fractional diving**: fix the most-fractional integer variable to its rounded value, re-solve the LP relaxation with the dual simplex, and repeat.
- **Coefficient diving**: fix integer variables by selecting those with the smallest objective impact, re-solving after each fixation.
- **Feasibility pump**: alternate between rounding the LP solution to the nearest integer point and projecting back to the LP feasible set. Anti-cycling perturbation is applied when the process stalls.
- **RINS (Relaxation Induced Neighborhood Search)**: fix variables where the LP relaxation and the current incumbent agree, then solve the resulting sub-problem.

### 11.7. Presolve

The presolve engine applies up to 20 rounds of the following reduction rules:

1. **Empty row removal**: rows with no nonzero coefficients are checked for feasibility and removed.
2. **Fixed variable substitution**: variables with equal lower and upper bounds are substituted out.
3. **Singleton row reduction**: rows with exactly one nonzero coefficient are used to tighten the variable's bounds, then removed.
4. **Singleton column reduction**: columns appearing in exactly one constraint are fixed at their optimal value.
5. **Forcing row detection**: rows where the only feasible activity is at the constraint bound force all participating variables to their bounds.
6. **Implied bound tightening**: row activity analysis derives tighter bounds on individual variables.
7. **Coefficient tightening** (MILP): for integer variables, constraint coefficients and right-hand sides are tightened using floor rounding.
8. **Probing** (MILP, binary): fix each binary variable to 0 and 1, solve the implications, and intersect the resulting bounds.
9. **Duplicate row/column detection**: framework for identifying and removing redundant constraints or variables.

After solving the reduced problem, the postsolve phase maps the solution back to the original variable space by replaying the presolve stack in reverse order.

### 11.8. Conflict Analysis

When a node LP is proven infeasible:

1. The Farkas ray is obtained from the LP solver.
2. The responsible rows are identified from the Farkas ray.
3. The responsible columns (branching decisions) are extracted.
4. A no-good cut is generated from the conflicting branching decisions and added to the cut pool.

---

## 12. Configuration Reference

### Numerical Constants

All constants are defined in `src/core/types.h`:

| Constant | Value | Purpose |
|---|---|---|
| `INF` | `1e30` | Represents positive infinity for bounds and activities. |
| `EPS_ZERO` | `1e-12` | Threshold for treating a value as zero. |
| `EPS_PIVOT` | `1e-10` | Minimum acceptable pivot element magnitude. |
| `EPS_FEASIBILITY` | `1e-8` | Feasibility tolerance for constraint satisfaction. |
| `EPS_OPTIMALITY` | `1e-8` | Optimality tolerance for reduced costs. |
| `EPS_INTEGER` | `1e-6` | Tolerance for determining whether a value is integer. |
| `MAX_ITER` | `1000000` | Default maximum simplex iterations. |
| `TIME_LIMIT` | `3600.0` | Default maximum solve time in seconds. |

### Type System

| Type | C++ Type | Description |
|---|---|---|
| `Real` | `double` | IEEE 754 64-bit floating point used for all numerical calculations. |
| `Index` | `int32_t` | Integer type used for matrix indexing and dimensions. |

### Enumerations

| Enum | Values | Description |
|---|---|---|
| `ObjectiveSense` | `MINIMIZE`, `MAXIMIZE` | Direction of optimization. |
| `VarType` | `CONTINUOUS`, `INTEGER`, `BINARY` | Variable integrality type. |
| `BoundType` | `FREE`, `LOWER`, `UPPER`, `BOTH`, `FIXED` | Variable bound classification. |
| `SolverStatus` | See Section 8 | Solve outcome status. |
| `LogLevel` | `TRACE`, `DEBUG`, `INFO`, `WARN`, `ERROR`, `OFF` | Logging verbosity. |

### Constraint Representation

Constraints use the ranged-row format: `row_lower[i] <= (Ax)[i] <= row_upper[i]`.

| Constraint Type | `row_lower` | `row_upper` |
|---|---|---|
| `<=` (less-than-or-equal) | `-INF` | finite value |
| `>=` (greater-than-or-equal) | finite value | `INF` |
| `=` (equality) | value | same value |
| Range | finite lower | finite upper (lower < upper) |

---

## 13. Testing

### Running All Tests

```bash
cd build
ctest --output-on-failure -j$(nproc)
```

### Test Suite Inventory

The full test suite consists of 12 CTest entries covering unit tests, integration tests, CLI contract tests, and benchmark smoke tests. The current verified result is **12/12 PASS**.

**Core Unit Tests** (`suplex_test_runner`):

This single binary contains 62 tests covering:

- Sparse matrix construction, SpMV, transpose SpMV, and CSC-CSR round-trip conversion.
- LU factorization, backward/forward solve, singular matrix detection, and eta updates.
- Primal simplex, dual simplex, warm-start, unbounded/infeasible detection.
- 49 boundary condition tests covering edge cases in sparse matrices, LU decomposition, basis management, pricing strategies, and simplex algorithm behavior.

**MILP Unit Tests** (4 GoogleTest binaries, requires GoogleTest):

- `suplex_p3_presolve_test`: 11 tests covering all presolve reduction rules, infeasibility detection, postsolve reconstruction, and no-op behavior.
- `suplex_p3_branching_test`: 8 tests covering pseudo-cost updates, reliability thresholds, most-fractional selection, and continuous/integer variable filtering.
- `suplex_p3_cuts_test`: 8 tests covering cut pool management, efficacy computation, MIR cuts, cover cuts, clique cuts, and parallel cut filtering.
- `suplex_p3_heuristics_test`: 6 tests covering simple rounding, fractional diving with mock LP solvers, and feasibility pump behavior.

**Integration Tests**:

- `suplex_c_api_client`: Pure C11 program that links against the shared library and solves a problem through the C API.
- `suplex_cli_help`: Verifies `--help` output.
- `suplex_cli_solve`: Verifies that solving `tiny.mps` produces `OPTIMAL` status.
- `suplex_cli_rejects_unknown_option`: Verifies that unknown CLI flags are rejected.
- `suplex_cli_full_contract`: Shell script verifying help text, version, solve output, solution file generation, error handling.
- `suplex_download_script_contracts`: Static verification that dataset download scripts have valid syntax and reference correct URLs.
- `suplex_benchmark_smoke`: Verifies that the benchmark runner produces valid CSV and JSON output.

**Conditional Tests**:

- `suplex_python_api`: Python API test. Requires `SUPLEX_BUILD_PYTHON=ON` and pybind11.
- CUDA tests: Require `SUPLEX_ENABLE_GPU=ON`, CUDA 12+, and a CUDA-capable device.

### Test Coverage by Component

| Component | Tests | Status |
|---|---|---|
| Sparse Matrix (CSC/CSR/Triplet) | 4 functional + 14 boundary | All pass |
| LU Factorization and Update | 4 functional + 9 boundary | All pass |
| Primal Simplex | 5 functional + 12 boundary | All pass |
| Dual Simplex | 3 functional + 2 boundary | All pass |
| Presolve | 11 tests | All pass |
| Branching | 8 tests | All pass |
| Cutting Planes | 8 tests | All pass |
| Heuristics | 6 tests | All pass |
| MPS/LP Parsing | 4 integration tests | All pass |
| C API | 3 integration tests + C11 client | All pass |
| CLI | 5 contract tests | All pass |
| GPU (CPU fallback) | 1 integration test | Passes |
| Benchmark Runner | 1 smoke test | Passes |

---

## 14. Benchmarking

### Benchmark Runner

The benchmark runner (`suplex_benchmark`) loads MPS files from a directory, solves each one, and compares the result against known reference objectives. It outputs results in both CSV and JSON format.

```bash
./build/suplex_benchmark
```

### Dataset Download Scripts

Three scripts are provided in `scripts/` to download standard optimization test sets:

| Script | Dataset | Description |
|---|---|---|
| `scripts/download_netlib.sh` | Netlib LP | Classic LP benchmark set (~90 problems). Builds the `emps.c` expander for compressed format. |
| `scripts/download_miplib.sh` | MIPLIB 2017 | Standard MILP benchmark collection. Downloads `benchmark.zip`. |
| `scripts/download_maros.sh` | Maros-Meszaros QP | QP test set from the Maros-Meszaros collection. Downloads three QPDATA archives. |

Usage:

```bash
bash scripts/download_netlib.sh
bash scripts/download_miplib.sh
bash scripts/download_maros.sh
```

These scripts download large files and are not run as part of the standard test suite.

---

## 15. Project Structure

```
suplex/
    CMakeLists.txt              Build system configuration
    Makefile                    GNU Make fallback
    README.md                   This file
    CHANGELOG.md                Version history
    PERSON3_NOTES.md            Person 3 development notes

    cli/
        main.cpp                CLI entry point and option parsing

    src/
        core/
            types.h             Type aliases (Real, Index), constants, enums
            sparse_matrix.h     SparseMatrixCSC, SparseMatrixCSR, TripletMatrix declarations
            sparse_matrix.cpp   Sparse matrix implementations
            problem.h           Problem class declaration
            problem.cpp         Problem class implementation
            solution.h          Solution struct declaration
            solution.cpp        Solution implementation

        simplex/
            lp_solver.h         Abstract LPSolver interface
            primal_simplex.h    PrimalSimplex declaration
            primal_simplex.cpp  Primal revised simplex implementation (410 lines)
            dual_simplex.h      DualSimplex declaration
            dual_simplex.cpp    Dual revised simplex implementation (409 lines)
            basis.h             Basis management declaration
            basis.cpp           Basis status tracking and crash procedures
            lu_factor.h         LUFactor declaration
            lu_factor.cpp       Sparse LU factorization with Markowitz ordering (352 lines)
            lu_update.h         LUUpdate declaration
            lu_update.cpp       Forrest-Tomlin eta update
            pricing.h           Pricing strategy declarations
            pricing.cpp         Dantzig, Steepest Edge, Devex, Bland implementations

        milp/
            milp_solver.h       MILPSolver declaration
            milp_solver.cpp     Top-level MILP orchestration
            presolve.h          Presolve class and PresolveRecord declarations
            presolve.cpp        Presolve engine with 9 reduction rules (888 lines)
            branch_bound.h      BranchBound declaration
            branch_bound.cpp    Branch-and-bound loop (473 lines)
            branching.h         Branching and PseudoCosts declarations
            branching.cpp       Variable selection strategies
            cuts.h              CutGenerator, CutPool, Cut declarations
            cuts.cpp            MIR, cover, clique, Gomory cut generation (402 lines)
            heuristics.h        Heuristics declaration
            heuristics.cpp      Rounding, diving, feasibility pump, RINS
            conflict.h          ConflictAnalysis declaration
            conflict.cpp        Farkas ray analysis and no-good cut generation
            node.h              BBNode, BranchDirection, NodeStatus, TreeStats declarations
            node.cpp            Node construction and delta-encoded bound changes

        io/
            mps_parser.h        MPSParser declaration
            mps_parser.cpp      MPS file parser (fixed and free format)
            lp_parser.h         LPParser declaration
            lp_parser.cpp       CPLEX LP file parser
            problem_reader.h    ProblemReader abstract interface
            sol_writer.h        SolutionWriter declaration
            sol_writer.cpp      Solution file output

        api/
            suplex_solver.h     Suplex high-level C++ solver class
            suplex_solver.cpp   Solver orchestration implementation
            suplex.h            C API header (extern "C")
            suplex_api.cpp      C API implementation
            python/
                bindings.cpp    pybind11 Python module

        gpu/
            gpu_manager.h       GPUManager declaration
            gpu_manager.cpp     GPU management with CPU fallback
            gpu_kernels.h       Kernel launch declarations
            cuda_spmv.cu        CUDA sparse matrix-vector multiply kernel
            cuda_pcg.cu         CUDA preconditioned conjugate gradient kernel

    data/
        examples/
            tiny.mps            Minimal LP in MPS format
            tiny.lp             Same LP in LP format
            tiny_mip.lp         Small binary MIP
            features.mps        MPS with advanced features
            compact_bounds.lp   LP with compact syntax
            branch_mip.lp       MIP requiring branching
            invalid_missing_end.mps  Malformed MPS for error testing

    tests/
        TEST_MATRIX.md          Complete test inventory document
        unit/
            test_framework.h    Lightweight test registration macros
            test_runner.cpp     Test runner entry point
            test_sparse_matrix.cpp  Sparse matrix unit tests
            test_lu.cpp         LU factorization and update tests
            test_simplex.cpp    Simplex algorithm tests
            test_boundary.cpp   49 boundary condition tests
            milp/
                test_presolve.cpp    Presolve unit tests (GoogleTest)
                test_branching.cpp   Branching unit tests (GoogleTest)
                test_cuts.cpp        Cut generation unit tests (GoogleTest)
                test_heuristics.cpp  Heuristic unit tests (GoogleTest)
        integration/
            test_integration.cpp    End-to-end integration tests
            test_c_api.c            Pure C API client test
            test_cli.sh             CLI contract verification script
            test_scripts.sh         Dataset download script verification
            test_benchmark.sh       Benchmark runner verification
            test_python_api.py      Python API test
        benchmarks/
            benchmark_runner.cpp    Benchmark runner with CSV/JSON output

    scripts/
        download_netlib.sh      Netlib LP dataset downloader
        download_miplib.sh      MIPLIB 2017 dataset downloader
        download_maros.sh       Maros-Meszaros QP dataset downloader

    .context/                   Internal project documentation
        00_project_manifest.md  Mission statement and team roles
        01_architecture_contracts.md
        02_tech_stack_and_env.md
        03_repository_map.md
        04_decision_log_adr.md
        05_task_dependency_board.md
        06_error_and_edge_case_log.md
        07_session_handoff.md
        ARCHITECTURE.md         Layered architecture documentation
        INTERFACES.md           C++ interface contracts
        PERSON1_SIMPLEX.md      Person 1 algorithm pseudocode
        PERSON2_INTERIOR_POINT.md   Person 2 design document
        PERSON3_MILP.md         Person 3 algorithm documentation
        PERSON4_INTEGRATION.md  Person 4 design document
```

Total source code: approximately 7,600 lines of production code and 2,600 lines of test code across 55 source files.

---

## 16. Design Decisions

### No External Dependencies

Suplex uses no external mathematical libraries (BLAS, LAPACK, Eigen, COIN-OR, HiGHS, GLPK, SCIP). All sparse linear algebra, factorizations, and solver algorithms are implemented from scratch. This ensures full control over numerical behavior, eliminates supply chain risk, and demonstrates capability for sovereign computing applications.

### Ranged Row Constraints

Instead of the standard equality form `Ax = b` with explicit slack variables, Suplex represents constraints as ranged rows: `row_lower <= Ax <= row_upper`. This avoids the combinatorial explosion of slack variables during model construction and is the natural representation for the simplex algorithm with bounded variables.

### CSC as Primary Matrix Format

The Compressed Sparse Column format is the primary storage format because simplex pricing and basis operations are fundamentally column-oriented. CSR is provided for row-oriented operations (presolve, GPU SpMV).

### Delta-Encoded Branch-and-Bound Nodes

Branch-and-bound nodes store only the bound changes relative to their parent, not a full copy of the problem. This keeps memory usage proportional to the depth of the tree rather than the number of nodes times the problem size.

### No Exceptions in Hot Paths

Solver functions do not throw exceptions. Control flow uses `SolverStatus` return codes. This eliminates exception-handling overhead in the inner simplex loop and makes the behavior fully deterministic.

### C++20 Standard

The codebase uses C++20 features including designated initializers, `std::span` compatibility, and concepts-ready template design. All code compiles cleanly with `-Wall -Wextra -Werror`.

### Thread Safety

- Solver instances are not thread-safe. A single instance must be used from one thread at a time.
- Const methods on sparse matrices are thread-safe (concurrent read access).
- The `GPUManager` singleton is internally synchronized for thread-safe access.

---

## 17. Known Limitations

- **Interior Point Method**: The IPM/QP solver (Person 2 scope) has design documentation but no production implementation in the current repository. Selecting `Algorithm::IPM` returns `NUMERICAL_ERROR` with a diagnostic message. The Cholesky factorization, AMD ordering, Mehrotra predictor-corrector, and crossover modules are not yet implemented.

- **Presolve singleton column reduction**: The `reduce_singleton_cols` routine has been disabled due to a constraint-violation bug where the greedy variable fixation does not always respect the full set of constraints. All other presolve reductions are active and verified.

- **GPU acceleration**: The CUDA kernels compile and link when `SUPLEX_ENABLE_GPU=ON` with CUDA 12+, but have only been verified via the CPU fallback path. Native GPU execution requires a CUDA-capable device.

- **QP support**: The `Problem` class supports an optional quadratic objective matrix `Q`, but no QP solver algorithm is implemented. QPS file parsing is accepted by the MPS parser but QP solving will return an error.

- **Parallel branch-and-bound**: The framework supports task-level parallelism for independent node processing, but the current implementation processes nodes sequentially. Thread pool integration is designed but not activated.

- **Problem scale**: The solver has been tested on small to medium problems (up to thousands of variables and constraints). Performance on large-scale industrial problems (millions of variables) has not been benchmarked.

---

## 18. Team Contributions

| Person | Domain | Components |
|---|---|---|
| Person 1 | Sparse Linear Algebra and Simplex Engine | `src/core/sparse_matrix.*`, `src/simplex/*` (primal simplex, dual simplex, LU factorization, LU update, basis management, pricing strategies), 62 unit and boundary tests |
| Person 2 | Interior Point Method and Quadratic Programming | Design documentation for AMD ordering, sparse Cholesky, Mehrotra IPM, crossover, and QP solver (`.context/PERSON2_INTERIOR_POINT.md`). Implementation pending. |
| Person 3 | Presolve Engine and MILP Framework | `src/milp/*` (presolve, branch-and-bound, branching, cuts, heuristics, conflict analysis, MILP solver, node), 33 GoogleTest unit tests |
| Person 4 | I/O, API, CLI, GPU, and Benchmarks | `src/io/*`, `src/api/*`, `cli/*`, `src/gpu/*`, `tests/integration/*`, `tests/benchmarks/*`, `scripts/*`, CMakeLists.txt, Makefile, dataset download scripts |

---

## 19. License

This project was developed for the Smart India Hackathon. License terms are determined by the hackathon organizers and the contributing team.
