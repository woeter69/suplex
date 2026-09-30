# PERSON 4: I/O, API, GPU Acceleration & Benchmarking

## Your Mission
You are the glue that holds Suplex together. You build: (1) file parsers that read industry-standard problem formats, (2) the command-line interface, (3) the C and Python APIs, (4) GPU-accelerated kernels for performance-critical operations, and (5) the benchmark infrastructure that validates the solver against recognized optimization benchmarks. Without your work, no one can test their code on real problems.

## Deliverables Checklist

- [x] `src/io/mps_parser.h/cpp` — MPS file parser (fixed + free format) (5 days)
- [x] `src/io/lp_parser.h/cpp` — LP file parser (3 days)
- [x] `src/io/sol_writer.h/cpp` — Solution output writer (1 day)
- [x] `cli/main.cpp` — Command-line interface (2 days)
- [x] `src/api/suplex.h` — C API header (2 days)
- [x] `src/api/suplex_api.cpp` — C API implementation (2 days)
- [x] `src/api/python/bindings.cpp` — Python pybind11 bindings (3 days)
- [x] `src/gpu/gpu_manager.h/cpp` — GPU device management (2 days)
- [x] `src/gpu/cuda_spmv.cu` — GPU sparse matrix-vector multiply (3 days)
- [x] `src/gpu/cuda_pcg.cu` — GPU preconditioned conjugate gradient (3 days)
- [x] `tests/benchmarks/` — Benchmark runner infrastructure (3 days)
- [x] `scripts/` — Download scripts for benchmark sets (2 days)
- [x] `CMakeLists.txt` — Root build system (2 days)
- [x] Integration testing (4 days)

Total: ~37 working days

## PART A: FILE PARSERS

## 1. MPS File Parser

The MPS (Mathematical Programming System) format is THE industry standard for optimization problems. Every benchmark set (Netlib, MIPLIB, Maros-Mészáros) uses it.

### 1.1 Fixed-Format MPS

The original fixed-column format:
```
NAME          problem_name
ROWS
 N  COST          <- objective (N = no constraint)
 L  ROW1          <- less-than-or-equal
 G  ROW2          <- greater-than-or-equal  
 E  ROW3          <- equality
COLUMNS
    X1        ROW1        1.0        COST        2.0
    X1        ROW2        3.0        ROW3        1.0
    X2        ROW1        4.0        COST        1.0
RHS
    RHS       ROW1        10.0       ROW2        5.0
RANGES
    RANGE     ROW1        2.0
BOUNDS
 UP BOUND     X1          100.0
 LO BOUND     X1          0.0
 BV BOUND     X3                     <- binary variable
 FX BOUND     X4          5.0        <- fixed variable
 FR BOUND     X5                     <- free variable
 MI BOUND     X6                     <- minus infinity lower bound
ENDATA
```

Field positions (1-indexed):
- Field 1: columns 2-3 (indicator/type)
- Field 2: columns 5-12 (name)
- Field 3: columns 15-22 (name)
- Field 4: columns 25-36 (value)
- Field 5: columns 40-47 (name)
- Field 6: columns 50-61 (value)

Provide FULL parsing pseudocode for each section.

### 1.2 Free-Format MPS
Modern extension: fields separated by whitespace, not column positions. Must auto-detect.

### 1.3 ROWS Section
```
parse_rows():
  for each line:
    type = first token (N, L, G, E)
    name = second token
    if type == 'N': set as objective row
    if type == 'L': add row with bounds (-INF, 0)  // after RHS processing
    if type == 'G': add row with bounds (0, INF)
    if type == 'E': add row with bounds (0, 0)
    store row_name → row_index mapping
```

### 1.4 COLUMNS Section
The trickiest section:
```
parse_columns():
  for each line:
    col_name = first token
    if col_name not seen: create new column
    // Handle MARKER lines for integer variables:
    if second token == "'MARKER'":
      if third token == "'INTORG'": start integer section
      if third token == "'INTEND'": end integer section
      continue
    // One or two entries per line:
    row_name1 = second token
    value1 = third token (as double)
    add entry (row_index[row_name1], col_index[col_name], value1) to triplet matrix
    if fourth token exists:
      row_name2 = fourth token
      value2 = fifth token
      add entry
```

### 1.5 RHS Section
```
parse_rhs():
  for each line:
    rhs_name = first token (usually "RHS" or similar)
    row_name1 = second token, value1 = third token
    // For L rows: row_upper = value1
    // For G rows: row_lower = value1  
    // For E rows: row_lower = row_upper = value1
    if fourth token exists: similar
```

### 1.6 RANGES Section
Ranges create ranged constraints:
```
parse_ranges():
  for each line:
    range_name = first token
    row_name = second token, value = third token
    // For L rows: row_lower = row_upper - |value|
    // For G rows: row_upper = row_lower + |value|
    // For E rows: depends on sign of value
```

Document the exact semantics for each row type.

### 1.7 BOUNDS Section
```
parse_bounds():
  for each line:
    bound_type = first token (UP, LO, FX, FR, MI, PL, BV, LI, UI)
    bound_name = second token
    col_name = third token
    value = fourth token (if present)
    
    switch bound_type:
      UP: col_upper[col] = value
      LO: col_lower[col] = value
      FX: col_lower[col] = col_upper[col] = value
      FR: col_lower[col] = -INF, col_upper[col] = INF
      MI: col_lower[col] = -INF
      PL: col_upper[col] = INF
      BV: col_lower[col] = 0, col_upper[col] = 1, var_type = BINARY
      LI: col_lower[col] = value, var_type = INTEGER
      UI: col_upper[col] = value, var_type = INTEGER
```

### 1.8 Error Handling
- Report line number and content for parse errors
- Handle: missing sections, duplicate names, unknown row/col references
- Gracefully handle: empty lines, comments (starting with *), extra whitespace
- Return meaningful error messages via `last_error()`

### 1.9 Performance
- Parse files using memory-mapped I/O for speed on large files
- Hash maps for name → index lookup
- Pre-allocate triplet matrix based on file size estimate
- Target: parse a 100MB MPS file in < 5 seconds

## 2. LP File Parser

The CPLEX LP format is more human-readable:
```
Minimize
 obj: 2 x1 + 3 x2 - x3
Subject To
 c1: x1 + x2 <= 10
 c2: x1 + x3 >= 5
 c3: x2 - x3 = 3
Bounds
 0 <= x1 <= 100
 x2 Free
Generals
 x3
End
```

Provide full parsing specification and pseudocode.

## 3. Solution Writer

Output format:
```
Suplex Solution File
Status: OPTIMAL
Objective: 12345.6789

Variables:
  x1    1.23456789e+02
  x2    0.00000000e+00
  ...

Constraints:
  c1    Activity: 9.99999999e+00  Dual: 1.23456789e+00
  c2    Activity: 5.00000000e+00  Dual: 0.00000000e+00
  ...

Statistics:
  Iterations: 1234
  Solve Time: 1.234 seconds
  Nodes: 567 (MILP only)
  Gap: 0.00% (MILP only)
```

## PART B: COMMAND-LINE INTERFACE

## 4. CLI Design

```
Usage: suplex [options] <problem_file>

Options:
  -h, --help              Show help message
  -v, --version           Show version
  --algorithm <algo>      LP algorithm: primal, dual, ipm, auto (default: auto)
  --presolve <on|off>     Enable/disable presolve (default: on)
  --time-limit <seconds>  Time limit (default: 3600)
  --threads <n>           Number of threads (default: auto)
  --gpu <on|off>          Enable GPU acceleration (default: off)
  --log-level <level>     Logging: trace, debug, info, warn, error, off (default: info)
  --mip-gap <value>       MIP gap tolerance (default: 1e-4)
  --solution <file>       Write solution to file
  --stats                 Print detailed statistics

Examples:
  suplex problem.mps
  suplex --algorithm ipm --gpu on large_problem.mps
  suplex --time-limit 600 --mip-gap 0.01 --solution out.sol mip_problem.mps
```

Implement argument parsing (use a simple hand-rolled parser or a header-only lib like cxxopts).

Provide the FULL main.cpp pseudocode:
```
main(argc, argv):
  1. Parse command-line arguments
  2. Detect file format from extension (.mps or .lp)
  3. Parse problem file → Problem object
  4. Print problem statistics (rows, cols, nonzeros, integers)
  5. Create Suplex solver, set options
  6. Call solver.solve()
  7. Print solution summary
  8. If --solution specified: write solution file
  9. Return 0 if optimal, 1 otherwise
```

## PART C: API DESIGN

## 5. C API (`src/api/suplex.h`)

Design a clean C API for FFI compatibility:
```c
#ifndef SUPLEX_H
#define SUPLEX_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SuplexSolver SuplexSolver;

// Lifecycle
SuplexSolver* suplex_create(void);
void suplex_destroy(SuplexSolver* solver);

// Problem loading
int suplex_read_mps(SuplexSolver* solver, const char* filename);
int suplex_read_lp(SuplexSolver* solver, const char* filename);

// Problem building (programmatic)
int suplex_set_dimensions(SuplexSolver* solver, int rows, int cols);
int suplex_set_objective(SuplexSolver* solver, const double* coeffs, int sense); // 0=min, 1=max
int suplex_set_constraint_matrix(SuplexSolver* solver, int nnz, const int* row_indices, const int* col_indices, const double* values);
int suplex_set_row_bounds(SuplexSolver* solver, const double* lower, const double* upper);
int suplex_set_col_bounds(SuplexSolver* solver, const double* lower, const double* upper);
int suplex_set_var_types(SuplexSolver* solver, const int* types); // 0=continuous, 1=integer, 2=binary

// Configuration
int suplex_set_algorithm(SuplexSolver* solver, int algo); // 0=auto, 1=primal, 2=dual, 3=ipm
int suplex_set_presolve(SuplexSolver* solver, int enabled);
int suplex_set_time_limit(SuplexSolver* solver, double seconds);
int suplex_set_log_level(SuplexSolver* solver, int level);
int suplex_set_threads(SuplexSolver* solver, int threads);
int suplex_set_gpu(SuplexSolver* solver, int enabled);
int suplex_set_mip_gap(SuplexSolver* solver, double gap);

// Solve
int suplex_solve(SuplexSolver* solver);

// Results
int suplex_get_status(const SuplexSolver* solver);
double suplex_get_objective(const SuplexSolver* solver);
int suplex_get_solution(const SuplexSolver* solver, double* values, int size);
int suplex_get_dual_values(const SuplexSolver* solver, double* values, int size);
int suplex_get_reduced_costs(const SuplexSolver* solver, double* values, int size);
double suplex_get_solve_time(const SuplexSolver* solver);
int64_t suplex_get_iterations(const SuplexSolver* solver);
int64_t suplex_get_nodes(const SuplexSolver* solver);
double suplex_get_mip_gap(const SuplexSolver* solver);

// Error handling
const char* suplex_get_error(const SuplexSolver* solver);
const char* suplex_version(void);

#ifdef __cplusplus
}
#endif
#endif
```

Document every function: parameters, return values, error codes.

## 6. Python Bindings

Using pybind11:
```python
import suplex

solver = suplex.Solver()
solver.read_mps("problem.mps")
solver.set_algorithm("dual")
solver.set_time_limit(600)
status = solver.solve()

if status == suplex.Status.OPTIMAL:
    print(f"Objective: {solver.objective}")
    x = solver.solution  # numpy array
    y = solver.dual_values  # numpy array
```

Provide the full pybind11 binding code with:
- Class `Solver` wrapping `Suplex`
- Enum `Status` wrapping `SolverStatus`
- Enum `Algorithm`
- NumPy integration for solution vectors
- Proper error handling (translate C++ errors to Python exceptions)

## PART D: GPU ACCELERATION

## 7. GPU Manager

```cpp
class GPUManager {
    bool gpu_available_;
    int device_id_;
    cudaStream_t stream_;
    
public:
    static GPUManager& instance();
    bool initialize(int device_id = 0);
    bool is_available() const;
    void synchronize();
    
    // Memory management
    void* device_malloc(size_t bytes);
    void device_free(void* ptr);
    void copy_to_device(void* dst, const void* src, size_t bytes);
    void copy_to_host(void* dst, const void* src, size_t bytes);
};
```

Handle:
- CUDA not installed → GPU features disabled gracefully
- Multiple GPUs → select via config
- CUDA errors → log and fall back to CPU

## 8. GPU Sparse Matrix-Vector Multiply (SpMV)

### 8.1 CSR-based SpMV Kernel
```cuda
__global__ void spmv_csr_kernel(
    int num_rows,
    const int* row_start,
    const int* col_index,
    const double* values,
    const double* x,
    double* y)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < num_rows) {
        double sum = 0.0;
        for (int j = row_start[row]; j < row_start[row+1]; j++) {
            sum += values[j] * x[col_index[j]];
        }
        y[row] = sum;
    }
}
```

### 8.2 Advanced: Warp-based SpMV
For rows with many nonzeros, use warp-level reduction:
- Each warp (32 threads) processes one row
- Threads within the warp cooperatively process nonzeros
- Warp shuffle for reduction

Provide the warp-based kernel code.

### 8.3 When to Use GPU SpMV
- Only beneficial for large matrices (>10,000 rows, >100,000 nonzeros)
- Data transfer overhead must be amortized
- Keep matrix on GPU, only transfer vectors
- Decision function: `use_gpu = (nnz > 100000) && gpu_available`

## 9. GPU Preconditioned Conjugate Gradient (PCG)

For Person 2's IPM: the normal equations can be solved iteratively instead of directly. PCG is GPU-friendly.

```
PCG on GPU:
  r = b - A*x          (GPU SpMV)
  z = M^{-1} r          (GPU — Jacobi preconditioner: z = D^{-1} r)
  p = z
  rz = r^T z            (GPU dot product with reduction)
  
  for k = 1, ..., max_iter:
    q = A*p              (GPU SpMV)
    alpha = rz / (p^T q) (GPU dot product)
    x += alpha * p       (GPU AXPY)
    r -= alpha * q       (GPU AXPY)
    if ||r|| < tol: break
    z = M^{-1} r         (GPU)
    rz_new = r^T z       (GPU dot product)
    beta = rz_new / rz
    p = z + beta * p     (GPU)
    rz = rz_new
```

Provide CUDA kernels for: dot product (with block reduction), AXPY, norm.

## PART E: BUILD SYSTEM

## 10. CMake Configuration

Provide the FULL root CMakeLists.txt:
```cmake
cmake_minimum_required(VERSION 3.24)
project(suplex VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(SUPLEX_BUILD_TESTS "Build tests" ON)
option(SUPLEX_BUILD_BENCHMARKS "Build benchmarks" ON)
option(SUPLEX_BUILD_PYTHON "Build Python bindings" OFF)
option(SUPLEX_ENABLE_GPU "Enable CUDA GPU support" OFF)

# Compiler flags
add_compile_options(-Wall -Wextra -Werror)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -DNDEBUG)
endif()

# CUDA
if(SUPLEX_ENABLE_GPU)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# Core library
add_library(suplex_core
    src/core/sparse_matrix.cpp
    src/core/problem.cpp
    src/core/solution.cpp
    src/core/timer.cpp
    src/simplex/lu_factor.cpp
    src/simplex/lu_update.cpp
    src/simplex/basis.cpp
    src/simplex/pricing.cpp
    src/simplex/primal_simplex.cpp
    src/simplex/dual_simplex.cpp
    src/ipm/ordering.cpp
    src/ipm/cholesky.cpp
    src/ipm/ipm_solver.cpp
    src/ipm/qp_solver.cpp
    src/ipm/crossover.cpp
    src/milp/presolve.cpp
    src/milp/branch_bound.cpp
    src/milp/node.cpp
    src/milp/branching.cpp
    src/milp/cuts.cpp
    src/milp/heuristics.cpp
    src/milp/conflict.cpp
    src/io/mps_parser.cpp
    src/io/lp_parser.cpp
    src/io/sol_writer.cpp
    src/api/suplex_api.cpp
)

target_include_directories(suplex_core PUBLIC src/)

# GPU
if(SUPLEX_ENABLE_GPU)
    target_sources(suplex_core PRIVATE
        src/gpu/cuda_spmv.cu
        src/gpu/cuda_pcg.cu
        src/gpu/gpu_manager.cpp
    )
    target_link_libraries(suplex_core CUDA::cudart)
    target_compile_definitions(suplex_core PUBLIC SUPLEX_GPU_ENABLED)
endif()

# CLI
add_executable(suplex_cli cli/main.cpp)
target_link_libraries(suplex_cli suplex_core)

# Tests
if(SUPLEX_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    # ... add test targets
endif()

# Python
if(SUPLEX_BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(pysuplex src/api/python/bindings.cpp)
    target_link_libraries(pysuplex PRIVATE suplex_core)
endif()
```

## PART F: BENCHMARK INFRASTRUCTURE

## 11. Benchmark Test Sets

### 11.1 Netlib LP
- URL: https://www.netlib.org/lp/data/
- ~95 LP problems in compressed MPS format
- Known optimal values for all problems
- Write a download script: `scripts/download_netlib.sh`

### 11.2 MIPLIB 2017
- URL: https://miplib.zib.de/
- ~240 MIP problems, categorized by difficulty
- Focus on "easy" subset for initial validation
- Write a download script: `scripts/download_miplib.sh`

### 11.3 Maros-Mészáros QP
- URL: https://www.doc.ic.ac.uk/~im/
- ~138 QP problems in QPS format (MPS extension with QUADOBJ section)
- Write a download script

## 12. Benchmark Runner

```cpp
struct BenchmarkResult {
    std::string problem_name;
    int num_rows, num_cols, nnz;
    SolverStatus status;
    double objective_value;
    double known_optimal;     // from reference file
    double gap_to_optimal;    // |ours - known| / |known|
    double solve_time;
    int64_t iterations;
    int64_t nodes;           // for MILP
    bool passed;             // gap < tolerance
};
```

Provide pseudocode for the benchmark runner:
```
run_benchmark(directory, reference_file):
  results = []
  for each .mps file in directory:
    problem = parse(file)
    solver = Suplex()
    solver.set_time_limit(3600)
    solver.solve()
    result = compare_with_reference(solver.solution(), reference_file)
    results.append(result)
  
  print_summary_table(results)
  print_pass_rate(results)
  save_csv(results)
```

Create reference files with known optimal values for each benchmark set.

## 13. Performance Profiling
- Time breakdown: parse time, presolve time, factorization time, solve time, postsolve time
- Memory usage: peak RSS
- For MILP: nodes/second, cuts added, heuristic solutions found
- Output as JSON for programmatic analysis

## 14. Top-Level `Suplex` Class (Wiring Everything Together)

Provide detailed pseudocode for the orchestration:
```
Suplex::solve():
  1. Validate problem
  2. Start timer
  3. Determine problem type (LP, MILP, QP)
  4. If presolve enabled:
     presolved = Presolve::apply(problem)
     if presolved.is_infeasible → return INFEASIBLE
  5. Select solver:
     if MILP:
       lp_solver = create_lp_solver(algorithm)
       milp_solver = MILPSolver()
       milp_solver.set_lp_solver(lp_solver)
       status = milp_solver.solve(presolved.reduced_problem, solution)
     elif QP:
       status = QPSolver().solve(presolved.reduced_problem, solution)
     else (LP):
       if algorithm == AUTO:
         if problem.num_rows > 10000: use IPM
         else: use Dual Simplex
       status = lp_solver.solve(presolved.reduced_problem, solution)
  6. If presolve enabled:
     solution = Presolve::postsolve(solution)
  7. Stop timer, record statistics
  8. Return status
```

## 15. Dependencies on Other Persons
- You CONSUME Person 1's `PrimalSimplexSolver`, `DualSimplexSolver`
- You CONSUME Person 2's `IPMSolver`, `QPSolver`
- You CONSUME Person 3's `Presolve`, `MILPSolver`
- You PROVIDE `Problem` objects to everyone (via parsers)
- You PROVIDE GPU acceleration hooks that Person 2 can optionally use
- You PROVIDE the top-level `Suplex` class that orchestrates everything

## 16. Testing Strategy

### Parser Tests
- Parse `afiro.mps` (Netlib) — verify dimensions (28 rows, 32 cols), objective name, variable bounds
- Parse free-format MPS — verify same result as fixed-format
- Parse LP format — verify equivalent to MPS parse
- Error handling: corrupt file, missing ENDATA, unknown section
- Round-trip: parse → write → parse → compare

### CLI Tests
- `suplex --help` — verify output
- `suplex afiro.mps` — verify OPTIMAL status
- `suplex --algorithm ipm afiro.mps` — verify IPM path
- `suplex --solution /tmp/test.sol afiro.mps` — verify solution file created

### API Tests
- C API: create solver, build problem programmatically, solve, get results
- Python API: equivalent to C API test
- Python: verify NumPy integration for solution vectors

### GPU Tests
- SpMV: compare GPU result with CPU result (tolerance 1e-12)
- PCG: solve SPD system, compare with direct solve
- Graceful fallback when GPU not available

## 17. Key References
- IBM ILOG CPLEX MPS format documentation
- COIN-OR MPS reader (for format edge cases)
- CUDA C++ Programming Guide (NVIDIA)
- pybind11 documentation
- Google Test documentation

## 18. Definition of Done
- MPS parser reads all 95 Netlib files without errors
- LP parser reads CPLEX-format LP files
- CLI works end-to-end for LP, MILP, QP
- C API functional for all operations
- Python bindings work with NumPy
- GPU SpMV matches CPU within 1e-12
- Benchmark runner produces results table for all test sets
- CMake builds successfully on Linux (GCC 12+, Clang 15+)
- No memory leaks
- Doxygen documentation complete
