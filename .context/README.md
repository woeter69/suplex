# Suplex Master Context

## 1. Project Identity
- **Name:** **Suplex** (Sovereign Unified Platform for Linear and Extended Optimization)
- **Mission:** Build a production-grade, from-scratch mathematical optimization solver supporting LP, MILP, QP (extensible to MIQP, NLP, MINLP).
- **Independence:** NOT built upon any existing solver library (no COIN-OR, no HiGHS, no GLPK, no SCIP internals).
- **Target:** Indian industrial optimization — refinery scheduling, crude blending, production planning, logistics, power dispatch, supply chain.
- **License:** Dual (AGPL-3.0 for open-source, commercial license for enterprise).

## 2. Problem Statement
India's refining, petrochemical, power, logistics, manufacturing and planning sectors depend on foreign solvers (CPLEX, Gurobi, Xpress). These have high license costs, restrictive models, no visibility into internals. Open-source alternatives (CBC, HiGHS, GLPK, SCIP) lag behind for large-scale MILP. The goal is a sovereign solver core — not a modeling environment.

## 3. Scope & Deliverables
- **Phase 1 (MVP):** LP solver (revised simplex + interior point), MPS/LP parser, CLI, Netlib benchmark passing
- **Phase 2:** MILP solver (branch-and-bound/cut), presolve, MIPLIB benchmarks
- **Phase 3:** QP solver, GPU acceleration, Python API
- **Phase 4:** Performance tuning, industrial validation, advanced heuristics

## 4. Technology Stack
- **Language:** C++20 (core solver), CUDA 12+ (GPU kernels), Python 3.10+ (bindings)
- **Build:** CMake 3.24+, Ninja
- **Testing:** Google Test, Google Benchmark
- **CI:** GitHub Actions
- **Dependencies:** NONE for solver core (no BLAS/LAPACK — we implement our own sparse LA). Only pybind11 for Python bindings, CUDA toolkit for GPU.
- **Platform:** Linux primary (x86_64, aarch64), Windows secondary

## 5. Repository Structure
```text
suplex/
├── .context/                    # Project documentation & work assignments
│   ├── README.md                # This file — master context
│   ├── ARCHITECTURE.md          # System architecture & design
│   ├── INTERFACES.md            # Shared interfaces & data structures
│   ├── PERSON1_SIMPLEX.md       # Work: Sparse LA & Simplex methods
│   ├── PERSON2_INTERIOR_POINT.md # Work: Interior Point & QP
│   ├── PERSON3_MILP.md          # Work: Presolve & MILP engine
│   └── PERSON4_INTEGRATION.md   # Work: I/O, API, GPU & Benchmarks
├── CMakeLists.txt               # Root CMake
├── src/
│   ├── core/                    # Shared data structures
│   │   ├── sparse_matrix.h/cpp  # CSC/CSR/Triplet sparse formats
│   │   ├── problem.h/cpp        # Optimization problem representation
│   │   ├── solution.h/cpp       # Solution container
│   │   ├── timer.h/cpp          # Performance timing
│   │   └── types.h              # Typedefs, constants, enums
│   ├── simplex/                 # Person 1's domain
│   │   ├── basis.h/cpp          # Basis management
│   │   ├── lu_factor.h/cpp      # Sparse LU factorization
│   │   ├── lu_update.h/cpp      # LU update (Forrest-Tomlin)
│   │   ├── pricing.h/cpp        # Pricing strategies
│   │   ├── primal_simplex.h/cpp
│   │   └── dual_simplex.h/cpp
│   ├── ipm/                     # Person 2's domain
│   │   ├── cholesky.h/cpp       # Sparse Cholesky
│   │   ├── ordering.h/cpp       # Fill-reducing ordering (AMD)
│   │   ├── ipm_solver.h/cpp     # Mehrotra predictor-corrector
│   │   ├── qp_solver.h/cpp      # QP interior point
│   │   └── crossover.h/cpp      # IPM → basic solution
│   ├── milp/                    # Person 3's domain
│   │   ├── presolve.h/cpp       # Presolve engine
│   │   ├── branch_bound.h/cpp   # B&B tree
│   │   ├── node.h/cpp           # Search tree node
│   │   ├── branching.h/cpp      # Variable selection
│   │   ├── cuts.h/cpp           # Cutting plane generators
│   │   ├── heuristics.h/cpp     # Primal heuristics
│   │   └── conflict.h/cpp       # Conflict analysis
│   ├── gpu/                     # Person 4's domain (GPU)
│   │   ├── cuda_spmv.cu         # Sparse matrix-vector multiply
│   │   ├── cuda_pcg.cu          # Preconditioned conjugate gradient
│   │   └── gpu_manager.h/cpp    # Device management
│   ├── io/                      # Person 4's domain (I/O)
│   │   ├── mps_parser.h/cpp     # MPS file reader
│   │   ├── lp_parser.h/cpp      # LP file reader
│   │   └── sol_writer.h/cpp     # Solution output
│   └── api/                     # Person 4's domain (API)
│       ├── suplex.h              # C API header
│       ├── suplex_api.cpp        # C API implementation
│       └── python/
│           └── bindings.cpp      # pybind11 bindings
├── cli/
│   └── main.cpp                 # CLI entry point
├── tests/
│   ├── unit/                    # Unit tests per module
│   ├── integration/             # Cross-module tests
│   └── benchmarks/              # Performance benchmarks
│       ├── netlib/              # Netlib LP test set
│       ├── miplib/              # MIPLIB 2017 test set
│       └── maros/               # Maros-Mészáros QP test set
├── data/                        # Sample problem files
├── docs/                        # Additional documentation
└── scripts/                     # Build/benchmark/CI scripts
```

## 6. Team Structure & Work Division
| Person | Role | Directory Ownership | Key Deliverables |
|--------|------|---------------------|------------------|
| **Person 1** | Sparse Linear Algebra & Simplex Engine | `src/simplex/`, `src/core/sparse_matrix.*` | CSC/CSR matrices, sparse LU, LU update, primal/dual simplex |
| **Person 2** | Interior Point Methods & QP Solver | `src/ipm/` | Sparse Cholesky, Mehrotra IPM, QP solver, crossover |
| **Person 3** | Presolve Engine & MILP Framework | `src/milp/` | Presolve, B&B tree, branching, cuts, heuristics |
| **Person 4** | I/O, API, GPU Acceleration & Benchmarks | `src/io/`, `src/api/`, `src/gpu/`, `cli/`, `tests/benchmarks/` | Parsers, CLI, C/Python API, CUDA kernels, benchmark suite |

## 7. Coding Standards
- C++20, `-Wall -Wextra -Werror`, no exceptions in hot paths (use error codes/`std::expected`).
- Naming: `snake_case` for functions/variables, `PascalCase` for classes/structs, `UPPER_SNAKE` for constants.
- All public headers use `#pragma once`.
- Every function has a Doxygen docstring.
- No raw `new`/`delete` — use `std::unique_ptr`, `std::vector`, arena allocators.
- Floating point: `double` everywhere (no `float` in solver core). Tolerances defined in `types.h`.
- Git: feature branches, conventional commits (`feat:`, `fix:`, `perf:`, `refactor:`, `test:`, `docs:`).

## 8. Shared Constants & Tolerances
| Constant | Value | Purpose |
|----------|-------|---------|
| `SUPLEX_EPS_ZERO` | 1e-12 | General zero tolerance |
| `SUPLEX_EPS_PIVOT` | 1e-10 | Minimum pivot element |
| `SUPLEX_EPS_FEASIBILITY` | 1e-8 | Primal/dual feasibility |
| `SUPLEX_EPS_OPTIMALITY` | 1e-8 | Reduced cost optimality |
| `SUPLEX_EPS_INTEGER` | 1e-6 | Integrality tolerance |
| `SUPLEX_INF` | 1e30 | Numerical infinity |
| `SUPLEX_MAX_ITER` | 1000000 | Default iteration limit |
| `SUPLEX_TIME_LIMIT` | 3600.0 | Default time limit (seconds) |

## 9. Integration Protocol
- **Week 1-2:** All persons set up skeleton code with interface headers (no implementation). Compile must succeed.
- **Week 3-6:** Independent implementation against interfaces. Unit tests.
- **Week 7-8:** Integration. Person 4 wires parsers → Person 2/3 presolve → Person 1/2 LP solve → Person 3 MILP solve → Person 4 output.
- **Week 9-10:** Benchmark validation against Netlib (LP), MIPLIB (MILP), Maros-Mészáros (QP).
- **Week 11-12:** Performance tuning, GPU integration, documentation.

## 10. Benchmark Targets
| Benchmark Set | # Problems | Target | Metric |
|---------------|-----------|--------|--------|
| Netlib LP | 95 | Solve 90+ correctly | Optimal within 1e-6 of known optimal |
| MIPLIB 2017 (easy) | 50 | Solve 35+ to optimality | Within 1% gap in 1 hour |
| Maros-Mészáros QP | 50 | Solve 40+ correctly | Optimal within 1e-5 |
