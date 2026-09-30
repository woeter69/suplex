# Architecture Contracts & System Diagrams

## System Layers
- Layer 0: Core Data Structures (`types.h`, `sparse_matrix.h`, `problem.h`, `solution.h`)
- Layer 1: Numerical Engines (`lu_factor.h`, `lu_update.h`, Cholesky, Ordering)
- Layer 2: Solver Algorithms (`primal_simplex.h`, `dual_simplex.h`, `ipm_solver.h`, `milp_solver.h`)
- Layer 3: Orchestration (`presolve.h`, `postsolve.h`, pipeline)
- Layer 4: Interfaces (`mps_parser.h`, `lp_parser.h`, C API, CLI)

## Component Boundaries for Person 1
- `SparseMatrixCSC` and `SparseMatrixCSR`: Used by Person 1 (Simplex), Person 2 (IPM normal equations), Person 3 (Presolve row operations), Person 4 (GPU SpMV).
- `LPSolver`: Base abstract class implemented by `PrimalSimplexSolver` (Person 1), `DualSimplexSolver` (Person 1), and `IPMSolver` (Person 2). Consumed by Person 3 in MILP B&B relaxations.
- `DualSimplexSolver::solve_from_basis`: Critical hot-path interface for Person 3 MILP engine. Must support resolving LPs under bound changes in 1 to 5 iterations.

## Data Flow for Person 1
```
Problem (CSC A, bounds, objective)
   │
   ▼
Basis Initialization (Crash / Warm Basis)
   │
   ▼
Sparse LU Factorization (Markowitz Ordering, Threshold Pivoting)
   │
   ├──> FTRAN (Bx = b) via L, U and Update transformations
   └──> BTRAN (B^T y = c) via L^T, U^T and Update transformations
   │
   ▼
Pricing (Dantzig / Steepest Edge / Devex)
   │
   ▼
Ratio Test (Primal: Harris/Bounded; Dual: Harris/Stability)
   │
   ▼
Basis Update & Forrest-Tomlin LU Update (or Refactorization trigger)
   │
   ▼
Convergence Check (Optimal / Infeasible / Unbounded)
   │
   ▼
Solution (x, y, reduced costs, status)
```
