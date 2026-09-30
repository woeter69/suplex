# CHANGELOG — Suplex

All notable changes to this project will be documented here.
Format: [Unreleased] / [version] — date, then Added / Changed / Fixed.

---

## [Unreleased] — 2026-09-30

### Added — Person 3: Presolve Engine & MILP Framework (branch: person3/milp-engine)

#### Core stubs (`src/core/`)
- `types.h` — Real=double, Index=int32_t, all enums (SolverStatus, VarType, LogLevel, etc.)
- `sparse_matrix.h` — SparseMatrixCSC/CSR/TripletMatrix interface
- `problem.h` — Problem class with ranged-row constraint representation
- `solution.h` — Solution struct (primal, dual, rc, gap, nodes, time)
- `core_stubs.cpp` — Temporary implementations of sparse matrix ops and Problem
  methods so Person 3 can compile standalone (replaced by Person 1's code on integration)

#### LP solver interface stub (`src/simplex/lp_solver.h`)
- Abstract `LPSolver` with `solve()`, `solve_from_basis()` (warm start), `get_basis()`, `get_farkas_ray()`

#### Presolve engine (`src/milp/presolve.h/.cpp`)
- `PresolveRuleType` enum — 12 reduction types
- `PresolveRecord` / `PresolveStack` — reversible LIFO undo stack
- `Presolve::apply()` — runs up to 20 rounds of all reduction rules
- `Presolve::postsolve()` — maps reduced-space Solution back to original variable space
- Reduction rules implemented:
  1. **Empty rows** — trivially satisfied or infeasibility detected
  2. **Fixed variables** — substitute x_j = lb_j, update row bounds, remove column
  3. **Singleton rows** — tighten variable bounds from single-nonzero constraint
  4. **Singleton columns** — substitute optimal x_j from single-nonzero column
  5. **Forcing rows** — fix all variables when max/min activity meets constraint bound
  6. **Implied bound tightening** — derive tighter bounds from row activities
  7. **Coefficient tightening** (MILP) — floor-round RHS for integer variables
  8. **Probing** (MILP, binary) — probe x_j=0/1, intersect implied bounds
  9. **Duplicate row/col** detection — stub, framework ready
- `compact()` helper — renumbers rows and cols after removals, rebuilds CSC matrix

#### B&B node (`src/milp/node.h/.cpp`)
- `BBNode` — delta-encoded bound changes from parent (no full problem copy)
- `BranchDirection`, `NodeStatus` enums
- `BoundChange` struct
- `TreeStats` — node counters, cuts added, heuristic solutions, max depth

#### Branch-and-bound manager (`src/milp/branch_bound.h/.cpp`)
- `MILPCallback` — on_new_incumbent, on_node_solved, should_terminate hooks
- `BranchBound::solve()` — full B&B loop:
  - Root LP → cuts → heuristics → B&B loop
  - Warm-start LP solves via `solve_from_basis()`
  - Prune by bound, infeasibility, and gap
  - Gap checking, time/node limits
  - Conflict analysis hook on infeasible nodes
- Node selection: BEST_FIRST, DEPTH_FIRST, HYBRID (plunge then best-first)
- Progress logging every 500 nodes

#### Branching (`src/milp/branching.h/.cpp`)
- `PseudoCosts` — running averages for down/up cost per unit
- Product scoring (`μ = 1/6`, SCIP default)
- Variable selection strategies:
  - Most fractional
  - Pseudo-cost branching
  - Strong branching (actual LP solves, top-K candidates)
  - **Reliability branching** (default) — strong branch until reliable, then pseudo-cost

#### Cutting planes (`src/milp/cuts.h/.cpp`)
- `Cut` struct — indices, coefficients, rhs, type, efficacy, age
- `CutPool` — add, age_all, remove_old
- `CutGenerator`:
  - **Gomory cuts** — from simplex tableau rows (tableau extraction hook for Person 1)
  - **MIR cuts** — mixed-integer rounding from each constraint row
  - **Cover cuts** — greedy minimal cover + sequential lifting for non-cover variables
  - **Clique cuts** — detect rows that are already sum-of-binaries <= 1 constraints
  - **select_cuts** — filter by efficacy > threshold, orthogonality < 0.9

#### Primal heuristics (`src/milp/heuristics.h/.cpp`)
- `Heuristics::simple_rounding` — round to nearest integer, check feasibility
- `Heuristics::fractional_diving` — fix most-fractional, re-solve dual simplex
- `Heuristics::coefficient_diving` — fix by smallest objective impact
- `Heuristics::feasibility_pump` — alternate round/project with anti-cycling perturbation
- `Heuristics::rins` — fix variables where LP and incumbent agree, solve sub-LP

#### Conflict analysis (`src/milp/conflict.h/.cpp`)
- `ConflictAnalysis::analyse()` — Farkas ray → responsible rows → conflict set → no-good cut
- `responsible_cols()` — identify columns in active Farkas rows
- `minimize_conflict()` — framework for redundancy removal
- `build_nogood_cut()` — generate no-good clause from branching decisions

#### Top-level MILP solver (`src/milp/milp_solver.h/.cpp`)
- `MILPSolver::solve()` — presolve → B&B → postsolve pipeline
- Falls back to LP solver for continuous problems
- Full configuration: gap, nodes, time, presolve toggle, log level, callback

#### Build system
- `CMakeLists.txt` — CMake 3.24+, FetchContent GoogleTest, 4 test executables
- `Makefile` — Direct g++ fallback (no cmake required), `make check` for syntax verification

#### Unit tests (`tests/unit/milp/`)
- `test_presolve.cpp` — 9 tests: empty rows, fixed vars, singleton row/col, forcing row,
  implied bounds (normal + infeasible), postsolve round-trip, no-op
- `test_branching.cpp` — 5 tests: pseudo-cost update, reliability, most-fractional,
  continuous skip, all-integer detection, observation recording
- `test_cuts.cpp` — 7 tests: CutPool add/age/remove, efficacy, MIR, cover, clique,
  parallel cut filtering
- `test_heuristics.cpp` — 7 tests: rounding (success/fail/bounds), diving with mock
  solver (feasible/infeasible), feasibility pump (already-integer)

#### Verification
- All 9 `.cpp` source files pass `g++ -std=c++20 -Wall -Wextra -fsyntax-only` with **0 errors**
- Warnings only: unused parameters in virtual default implementations (intentional)

---

## [Context] — 2026-09-30

### Added
- `.context/` folder with 7 documentation files (~3,900 lines total):
  - `README.md` — master project context, tech stack, team table, benchmark targets
  - `ARCHITECTURE.md` — 5-layer architecture with Mermaid diagrams
  - `INTERFACES.md` — complete C++ interface contracts for all modules
  - `PERSON1_SIMPLEX.md` — sparse LU, simplex algorithms (~1,114 lines, full pseudocode)
  - `PERSON2_INTERIOR_POINT.md` — AMD ordering, Cholesky, Mehrotra IPM, QP (~531 lines)
  - `PERSON3_MILP.md` — presolve, B&B, cuts, heuristics, conflict analysis (~870 lines)
  - `PERSON4_INTEGRATION.md` — parsers, CLI, C/Python API, CUDA, benchmarks (~671 lines)
- Git worktree at `/mnt/d_drive/programs/suplex-person3` on branch `person3/milp-engine`
