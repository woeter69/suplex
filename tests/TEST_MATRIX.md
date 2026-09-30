# Suplex Complete Test Matrix

## 1. Purpose and Status Vocabulary

This document is the test inventory for all four workstreams and the combined
solver. It records what each test proves, where it lives, how it runs, and its
current verification state.

| Status | Meaning |
|---|---|
| **PASS** | Executed successfully in the current GCC/Linux workspace. |
| **SOURCE-READY** | Test source and build registration exist, but an optional dependency is unavailable here. |
| **BLOCKED** | The production module required by the test does not exist in this repository. |
| **CONDITIONAL** | Runs only when its build option/toolchain/hardware is enabled. |

Current verified result:

- Native custom suite: **74/74 PASS**.
- Pure C11 ABI client: **PASS**.
- CLI shell contract: **PASS**.
- Benchmark runner contract: **PASS**.
- Dataset-script static contracts: **PASS**.
- Person 3 GoogleTest suites: **SOURCE-READY**; GoogleTest is not installed here.
- Person 2 IPM/QP algorithm tests: **BLOCKED**; `src/ipm/` is absent.
- Python runtime test: **CONDITIONAL** on `SUPLEX_BUILD_PYTHON=ON` and pybind11.
- Native CUDA execution: **CONDITIONAL** on `SUPLEX_ENABLE_GPU=ON`, CUDA 12+, and a CUDA device.

## 2. Build and Run Commands

### CMake/CTest

```bash
cmake -S . -B build \
  -DSUPLEX_BUILD_TESTS=ON \
  -DSUPLEX_BUILD_BENCHMARKS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Optional Python and CUDA coverage:

```bash
cmake -S . -B build-full \
  -DSUPLEX_BUILD_TESTS=ON \
  -DSUPLEX_BUILD_BENCHMARKS=ON \
  -DSUPLEX_BUILD_PYTHON=ON \
  -DSUPLEX_ENABLE_GPU=ON
cmake --build build-full -j
ctest --test-dir build-full --output-on-failure
```

### Make fallback

```bash
make -j test
make -j benchmark
```

### Direct GCC command used in this workspace

This environment has GCC but no CMake, Make, Ninja, GoogleTest, CUDA, or
pybind11. The native suite was therefore compiled directly:

```bash
g++ -std=c++20 -O0 -Wall -Wextra -Werror -I. \
  -DSUPLEX_SOURCE_DIR=\"$PWD\" \
  $(find src/core src/simplex src/milp src/io src/api src/gpu \
    -name '*.cpp' ! -path '*/python/*') \
  tests/unit/test_runner.cpp \
  tests/unit/test_sparse_matrix.cpp \
  tests/unit/test_lu.cpp \
  tests/unit/test_simplex.cpp \
  tests/unit/test_boundary.cpp \
  tests/integration/test_integration.cpp \
  -o /tmp/suplex_all_tests
/tmp/suplex_all_tests
```

## 3. Person 1 — Core, Sparse Linear Algebra, and Simplex

Person 1 tests use the lightweight registry in
`tests/unit/test_framework.h`. They are built into `suplex_test_runner`.

### 3.1 Sparse matrix tests

Source: `tests/unit/test_sparse_matrix.cpp`

| Test | What it verifies | Expected result |
|---|---|---|
| `test_csc_from_triplets` | COO/triplet assembly produces valid sorted CSC storage and preserves coefficients. | Dimensions, nonzeros, and values match. |
| `test_csc_spmv` | CSC matrix-vector multiplication computes `y = Ax`. | Dense reference values match within tolerance. |
| `test_csc_transpose_spmv` | Transpose multiplication computes `y = A^T x`. | Dense reference values match. |
| `test_csc_to_csr_roundtrip` | CSC→CSR→CSC conversion preserves the matrix. | Every entry matches the original. |

### 3.2 Sparse LU and update tests

Source: `tests/unit/test_lu.cpp`

| Test | What it verifies | Expected result |
|---|---|---|
| `test_lu_small` | Factorization and forward solve on a small nonsingular basis. | Recovered solution satisfies `Bx=b`. |
| `test_lu_btran` | Backward-transpose solve. | Recovered solution satisfies `B^T x=b`. |
| `test_lu_singular` | Singular basis detection. | Factorization rejects the matrix. |
| `test_lu_update` | Eta/product-form basis update followed by solve. | Updated solve matches the modified basis. |

### 3.3 Primary simplex tests

Source: `tests/unit/test_simplex.cpp`

| Test | What it verifies | Expected result |
|---|---|---|
| `test_primal_simplex_tiny` | Primal revised simplex on a bounded LP. | `OPTIMAL` with known objective/primal values. |
| `test_dual_simplex_tiny` | Dual revised simplex on the same LP. | Same optimum as primal simplex. |
| `test_unbounded` | Unbounded-ray/status detection. | `UNBOUNDED`. |
| `test_infeasible` | Contradictory row/bound detection. | `INFEASIBLE`. |
| `test_warm_start` | Dual re-optimization from a supplied basis. | Warm solve reaches the expected optimum. |

### 3.4 Core and sparse boundary tests

Source: `tests/unit/test_boundary.cpp`

| Test | What it verifies |
|---|---|
| `test_empty_csc_from_triplets` | Zero-by-zero triplet conversion is valid. |
| `test_csc_from_triplets_no_entries` | Nonzero dimensions with zero stored entries are valid. |
| `test_csc_1x1_matrix` | Smallest nonempty CSC matrix. |
| `test_csc_identity_matrix` | Identity assembly and access. |
| `test_csc_out_of_bounds_get` | Invalid CSC lookup is handled safely. |
| `test_csc_spmv_zero_vector` | Multiplication by zero produces zero. |
| `test_csc_dot_col_out_of_range` | Invalid column-dot request is safe. |
| `test_csc_multiply_col_out_of_range` | Invalid column AXPY request is safe. |
| `test_csc_multiply_col_zero_scalar` | Zero-scaled column update is a no-op. |
| `test_csc_multiple_duplicates_same_entry` | Duplicate triplets are accumulated correctly. |
| `test_csr_empty_multiply` | Empty CSR multiplication is valid. |
| `test_csr_out_of_bounds_get` | Invalid CSR lookup is safe. |
| `test_triplet_clear` | Triplet reset clears dimensions and storage. |
| `test_csc_clear` | CSC reset clears dimensions and storage. |
| `test_problem_zero_dimensions` | Empty `Problem` is dimensionally consistent. |
| `test_problem_validate_mismatched` | Array/matrix dimension mismatch is rejected. |
| `test_problem_validate_inverted_bounds` | Lower bound above upper bound is rejected. |
| `test_problem_set_bounds_out_of_range` | Invalid bound setter index does not corrupt state. |
| `test_problem_is_mip` | Integer/binary type detection and counting. |
| `test_problem_quadratic` | Optional quadratic matrix storage and QP detection. |
| `test_solution_reset` | Solution reset restores default status and clears vectors/statistics. |
| `test_solution_is_feasible` | Feasible-status classification. |

### 3.5 LU, basis, pricing, and simplex boundary tests

Source: `tests/unit/test_boundary.cpp`

| Test | What it verifies |
|---|---|
| `test_lu_1x1_factorize` | One-element basis factorization. |
| `test_lu_identity_factorize` | Identity factorization and solve. |
| `test_lu_near_singular` | Pivot tolerance catches near singularity. |
| `test_lu_zero_dimension` | Empty factorization does not crash. |
| `test_lu_dimension_mismatch` | RHS/basis dimension mismatch is rejected. |
| `test_lu_btran_after_update` | BTRAN remains correct after eta update. |
| `test_lu_update_near_zero_pivot` | Unsafe eta pivot is rejected. |
| `test_lu_update_out_of_range` | Invalid leaving-row index is rejected. |
| `test_lu_should_refactorize` | Update-count/stability refactorization trigger. |
| `test_basis_pivot_out_of_range` | Invalid basis pivot is safe. |
| `test_basis_set_status_out_of_range` | Invalid variable-status write is safe. |
| `test_pricing_bland_rule` | Bland selection returns the lowest eligible index. |
| `test_pricing_no_eligible` | Pricing reports no entering variable when optimal. |
| `test_pricing_free_variable` | Free-variable reduced-cost handling. |
| `test_pricing_fixed_variable_skipped` | Fixed variables cannot enter the basis. |
| `test_simplex_maximization` | Objective-sense conversion for maximization. |
| `test_simplex_equality_constraints` | Equality-row handling. |
| `test_simplex_fixed_variable` | Fixed column bound handling. |
| `test_simplex_single_variable` | One-column LP solve. |
| `test_simplex_upper_bounded_variables` | Finite upper-bound pivots and optimum. |
| `test_simplex_dual_infeasible_detected` | Dual simplex infeasibility status. |
| `test_simplex_iteration_limit` | Iteration-limit termination. |
| `test_simplex_larger_problem` | Multirow/multicolumn regression solve. |
| `test_simplex_primal_dual_agree` | Primal and dual algorithms return the same objective. |
| `test_simplex_zero_objective` | Feasible zero-objective model. |
| `test_dual_simplex_warm_start_tighten_both` | Warm start after lower/upper bound tightening. |
| `test_to_string_coverage` | Every `SolverStatus` has a stable string representation. |

Person 1 current result: **62/62 PASS**.

## 4. Person 2 — IPM, Cholesky, Ordering, Crossover, and QP

### 4.1 Current repository state

The integration history contains no `src/ipm/` directory. `origin/p2` contains
context documentation but no ordering, Cholesky, IPM, crossover, or QP source.
Algorithm-level P2 tests cannot be compiled honestly until those production
classes exist.

### 4.2 Executable P2 boundary test

Source: `tests/integration/test_integration.cpp`

| Test | What it verifies | Status |
|---|---|---|
| `test_missing_ipm_reports_error` | Selecting `Algorithm::IPM` without P2 never silently falls back to a different algorithm; it returns `NUMERICAL_ERROR` and a diagnostic naming P2. | **PASS** |

### 4.3 Required P2 algorithm tests once source is merged

These cases are mandatory activation criteria, currently **BLOCKED**:

| Planned test | Production component | Required assertion |
|---|---|---|
| AMD permutation validity | `ordering.*` | Permutation is bijective and reduces/does not unexpectedly explode fill on reference matrices. |
| Cholesky SPD solve | `cholesky.*` | `LL^T x=b` residual is below `1e-10`. |
| Cholesky non-SPD rejection | `cholesky.*` | Indefinite matrix is rejected without NaNs. |
| IPM feasible LP | `ipm_solver.*` | Matches simplex objective within `1e-7`. |
| IPM infeasible LP | `ipm_solver.*` | Returns `INFEASIBLE` or a valid certificate status. |
| IPM unbounded LP | `ipm_solver.*` | Returns `UNBOUNDED`/`INF_OR_UNBD`. |
| IPM iteration/time limits | `ipm_solver.*` | Stops with the configured status and finite iterates. |
| Predictor-corrector residual reduction | `ipm_solver.*` | Primal, dual, and complementarity residuals decrease. |
| Crossover basis validity | `crossover.*` | Returned basis has exactly `m` nonsingular basic columns. |
| Convex QP optimum | `qp_solver.*` | Known QP objective/primal values match within `1e-6`. |
| QP rejects indefinite Q | `qp_solver.*` | Nonconvex Q is diagnosed safely. |
| Maros QPS parse/solve | QPS parser + QP solver | Reference optimum matches within `1e-5`. |

P2 is not reported as passing; it is explicitly source-blocked.

## 5. Person 3 — Presolve and MILP

The standalone P3 tests use GoogleTest and are registered automatically when
`GTest` is discoverable. CMake option: `SUPLEX_BUILD_P3_GTESTS=ON`.

### 5.1 Presolve suite

Source: `tests/unit/milp/test_presolve.cpp`

| Test | What it verifies |
|---|---|
| `Presolve.EmptyRowSatisfied` | A redundant empty row is removed. |
| `Presolve.EmptyRowInfeasible` | An empty row excluding zero proves infeasibility. |
| `Presolve.FixedVariableRemoved` | Fixed column substitution/removal. |
| `Presolve.SingletonRowTightensUpperBound` | Positive singleton row tightens a column upper bound. |
| `Presolve.SingletonRowInfeasible` | Singleton contradiction is detected. |
| `Presolve.SingletonRowNegativeCoefficient` | Bound propagation handles coefficient sign. |
| `Presolve.ForcingRow` | Forcing constraint fixes participating variables. |
| `Presolve.ImpliedBoundTightening` | Row activity derives tighter implied bounds. |
| `Presolve.ImpliedBoundInfeasible` | Implied-bound contradiction is detected. |
| `Presolve.PostsolveSingletonRow` | Reduced solution is reconstructed through singleton records. |
| `Presolve.NoReductionOnGenericProblem` | Generic model survives a no-op presolve pass. |

Important integration finding: the current P3 presolve/postsolve mapping is not
safe for all-column elimination and a simple binary singleton model. Top-level
orchestration therefore preserves correctness by solving in original space.
Presolve tests remain source-ready, but presolve is not declared production-ready.

### 5.2 Branching suite

Source: `tests/unit/milp/test_branching.cpp`

| Test | What it verifies |
|---|---|
| `PseudoCosts.UpdateDown` | Down-branch pseudo-cost averaging. |
| `PseudoCosts.UpdateUp` | Up-branch pseudo-cost averaging. |
| `PseudoCosts.ReliabilityThreshold` | Reliability requires the configured observation count. |
| `PseudoCosts.ScorePositive` | Fractional candidate receives a positive score. |
| `Branching.MostFractionalSelected` | Most-fractional integer variable selection. |
| `Branching.IntegerVariableSkipped` | Already-integral integer variables are ignored. |
| `Branching.ContinuousSkipped` | Continuous columns are never branching candidates. |
| `Branching.PseudoCostUpdateRecorded` | Branch observations update stored statistics. |

### 5.3 Cut suite

Source: `tests/unit/milp/test_cuts.cpp`

| Test | What it verifies |
|---|---|
| `CutPool.AddAndAge` | Cut insertion and age tracking. |
| `CutPool.RemoveOld` | Stale cuts are removed. |
| `CutGenerator.EfficacyPositiveWhenViolated` | Violated cut has positive efficacy. |
| `CutGenerator.CoverCutFoundForKnapsack` | Cover separation finds a violated binary knapsack cover. |
| `CutGenerator.MinimalCoverCorrect` | Extracted cover is minimal. |
| `CutGenerator.MIRCutsNoSEGFAULT` | MIR separation is safe on representative rows. |
| `CutGenerator.CliqueCutFromClique` | Clique inequality construction. |
| `CutGenerator.SelectionFiltersParallelCuts` | Near-parallel/duplicate cuts are filtered. |

### 5.4 Heuristic suite

Source: `tests/unit/milp/test_heuristics.cpp`

| Test | What it verifies |
|---|---|
| `Heuristics.SimpleRoundingSuccess` | Rounding returns a feasible integer solution. |
| `Heuristics.SimpleRoundingFailsWhenInfeasible` | Infeasible rounded point is rejected. |
| `Heuristics.SimpleRoundingRespectsBounds` | Rounding never violates column bounds. |
| `Heuristics.FractionalDivingWithMockSolver` | Diving integrates with the LP-solver contract. |
| `Heuristics.FractionalDivingInfeasibleMock` | Infeasible dive terminates safely. |
| `Heuristics.FeasibilityPumpAlreadyInteger` | Pump immediately accepts an integral point. |

### 5.5 Executed MILP integration tests

Source: `tests/integration/test_integration.cpp`

| Test | What it verifies | Result |
|---|---|---|
| `test_milp_end_to_end` | LP parser → binary types → MILP dispatch → integer optimum on `tiny_mip.lp`. | **PASS** |
| `test_milp_fractional_root_branches` | Root relaxation is fractional; B&B must branch and return integer `x=0`, objective `0`, with at least one processed node. | **PASS** |

## 6. Person 4 — I/O, API, CLI, GPU, and Benchmarks

### 6.1 MPS and LP parsing

Source: `tests/integration/test_integration.cpp`

| Test | What it verifies | Fixture | Result |
|---|---|---|---|
| `test_mps_and_lp_equivalence` | Equivalent MPS and LP files create equal dimensions/nonzero counts and preserve names. | `tiny.mps`, `tiny.lp` | **PASS** |
| `test_mps_features` | Inline `OBJSENSE`, `RANGES`, free bounds, binary bounds, and `INTORG/INTEND`. | `features.mps` | **PASS** |
| `test_lp_compact_bounds_and_types` | Operators without surrounding spaces, `*` multiplication, free variables, and general integers. | `compact_bounds.lp` | **PASS** |
| `test_parser_error_diagnostic` | Missing `ENDATA` is rejected with line-aware diagnostic text. | `invalid_missing_end.mps` | **PASS** |

### 6.2 Solution writer

| Test | What it verifies | Result |
|---|---|---|
| `test_solution_writer` | Writes status, objective, named variables, activities, duals, and statistics; output is nonempty. | **PASS** |

### 6.3 C++ and C APIs

| Test | Source | What it verifies | Result |
|---|---|---|---|
| `test_top_level_file_solve` | `test_integration.cpp` | C++ facade loads MPS, selects primal simplex, and returns objective `3`. | **PASS** |
| `test_c_api_programmatic_problem` | `test_integration.cpp` | C API dimensions/objective/COO matrix/bounds → solve → primal copy. | **PASS** |
| `test_c_api_validation_and_buffer_size` | `test_integration.cpp` | Invalid dimensions set an error; undersized result buffer returns `SUPLEX_ERROR_BUFFER_TOO_SMALL`. | **PASS** |
| Pure C11 client | `tests/integration/test_c_api.c` | Header compiles as C11, links to C++ shared/static library, solves, and reads `x=2`. | **PASS** |

### 6.4 GPU manager, SpMV, and PCG

| Test | What it verifies | Result |
|---|---|---|
| `test_gpu_cpu_fallback` | CSR SpMV computes `[6,7]`; Jacobi-PCG solves a 2×2 SPD system within `1e-10`; works without CUDA. | **PASS** |
| CUDA-enabled build | Compiles `cuda_spmv.cu` and `cuda_pcg.cu`, initializes a CUDA device, and makes large matrices eligible for device execution. | **CONDITIONAL** |

The CPU fallback test is always required. Native CUDA execution additionally
requires CUDA 12+ and a supported device and was not available here.

### 6.5 CLI contract

Source: `tests/integration/test_cli.sh`

| Check | Expected result |
|---|---|
| `suplex --help` | Contains `Usage: suplex`. |
| `suplex --version` | Contains `Suplex 0.1.0`. |
| Solve `tiny.mps` | Objective is `3`. |
| Solve `tiny.lp` | Status is `OPTIMAL`. |
| `--solution` | File exists and contains `Status: OPTIMAL`. |
| Unknown option | Nonzero exit. |
| Unsupported extension | Nonzero exit. |

Current result: **PASS**.

### 6.6 Benchmark runner contract

Source: `tests/integration/test_benchmark.sh`

| Check | Expected result |
|---|---|
| Parse reference CSV | `tiny` reference objective `3` is loaded. |
| Run corpus directory | The copied MPS fixture is solved. |
| CSV output | Dimensions/status/objective row exists and `passed=1`. |
| JSON output | Contains `problem=tiny` and `passed=true`. |

Current result: **PASS**.

### 6.7 Dataset download scripts

Source: `tests/integration/test_scripts.sh`

| Check | Expected result |
|---|---|
| Shell syntax | All three scripts pass `bash -n`. |
| Netlib contract | Script references and builds the official `emps.c` expander. |
| MIPLIB contract | Script resolves/downloads `benchmark.zip`. |
| Maros contract | Script resolves all three `QPDATA` archives. |

This is deliberately a static/offline contract. Full corpus downloads are not
run in ordinary CI because they are large and network-dependent.

### 6.8 Python binding

Source: `tests/integration/test_python_api.py`

Enabled only by `SUPLEX_BUILD_PYTHON=ON`.

| Check | Expected result |
|---|---|
| Import `pysuplex` | Module and enums load. |
| Read LP and solve | Status `OPTIMAL`, objective `3`. |
| NumPy properties | Primal and dual arrays both have shape `(2,)`. |

Current status: **CONDITIONAL**, because pybind11 is not installed here.

## 7. Combined End-to-End Tests

These tests cross ownership boundaries and are the release gate for integration.

| Flow | Components crossed | Tests/checks | Current result |
|---|---|---|---|
| Text model parity | P4 MPS + P4 LP + P1 core | `test_mps_and_lp_equivalence` | **PASS** |
| LP solve | P4 parser → P4 facade → P1 simplex → P4 solution | `test_top_level_file_solve`, CLI MPS/LP checks | **PASS** |
| Programmatic solve | C client → P4 C ABI → P1 core/simplex | C11 client, `test_c_api_programmatic_problem` | **PASS** |
| MILP integral root | P4 LP parser → P3 MILP → P1 LP relaxation | `test_milp_end_to_end` | **PASS** |
| MILP branching | P4 LP parser → P3 B&B/warm start → P1 simplex | `test_milp_fractional_root_branches` | **PASS** |
| GPU-independent linear solve | P4 GPU manager → P1 sparse formats | `test_gpu_cpu_fallback` | **PASS** |
| Benchmark reporting | P4 parser/facade/solver → CSV + JSON | `test_benchmark.sh` | **PASS** |
| IPM selection without P2 | P4 facade → missing P2 boundary | `test_missing_ipm_reports_error` | **PASS** diagnostic behavior; algorithm remains **BLOCKED** |

## 8. Release Gates

A commit is ready for integration when all applicable gates hold:

1. The native registry reports zero failures.
2. The pure C client compiles and exits zero.
3. CLI and benchmark shell contracts exit zero.
4. `git diff --check` reports no whitespace errors.
5. Optional Python/CUDA tests pass when those options are enabled.
6. P3 presolve must remain disabled at orchestration level until its documented
   postsolve/binary reductions are repaired and regression-tested.
7. P2 must not be marked complete until `src/ipm/` exists and every planned P2
   algorithm test in Section 4.3 is implemented and passing.
