# Project Manifest & Mission Objectives

## Core Mission
Build Suplex (Sovereign Unified Platform for Linear and Extended Optimization), a high-performance, from-scratch mathematical optimization solver supporting LP, MILP, and QP for the Smart India Hackathon and sovereign industrial optimization needs.

## Non-Negotiable Domain Rules
- No external solver dependencies (no COIN-OR, HiGHS, GLPK, SCIP, Eigen, BLAS/LAPACK). All sparse linear algebra and solver cores are implemented in-house.
- Strict precision: IEEE 754 64-bit double precision (`Real`) for all numerical calculations.
- Clean C++20 standards, compiled with `-Wall -Wextra -Werror`. No exceptions in solver hot paths; return `SolverStatus` or status codes.
- Modular architecture with clean contracts between Person 1 (Linear Algebra & Simplex), Person 2 (Interior Point & QP), Person 3 (Presolve & MILP), and Person 4 (I/O, API, GPU & Benchmarks).
- Zero emojis across all code, commits, and documentation.

## Person 1 Role & Ownership
- Owner: Person 1
- Domain: Sparse Linear Algebra & Simplex Engine (`src/core/sparse_matrix.*`, `src/simplex/*`)
- Deliverables:
  - Compressed Sparse Column (CSC), Compressed Sparse Row (CSR), and Triplet matrix formats with efficient arithmetic.
  - Sparse LU Factorization with Markowitz ordering, threshold pivoting, and iterative refinement.
  - LU Update using Forrest-Tomlin / eta transformations and refactorization heuristics.
  - Basis management, status tracking, and triangular crash procedures.
  - Pricing strategies: Dantzig, Steepest Edge, Devex, and Bland anti-cycling rule.
  - Revised Primal Simplex Solver for bounded variables.
  - Revised Dual Simplex Solver with warm start (`solve_from_basis`) for Branch-and-Bound integration.
  - Complete unit testing suite and benchmark readiness.
