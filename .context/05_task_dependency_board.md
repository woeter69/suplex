# Unified Task Matrix

## Epics & Subtasks for Person 1

### Epic 1: Core Foundation & Sparse Linear Algebra
- [x] Task 1.1: `src/core/types.h` (Aliases, tolerances, solver statuses)
- [x] Task 1.2: `src/core/sparse_matrix.h` & `src/core/sparse_matrix.cpp` (CSC, CSR, Triplet, SpMV, transpose, dot products)
- [x] Task 1.3: `src/core/problem.h` & `src/core/problem.cpp` (Problem container, row ranges, variable types)
- [x] Task 1.4: `src/core/solution.h` & `src/core/solution.cpp` (Solution container, metrics)

### Epic 2: Numerical Engine (Sparse LU Factorization & Update)
- [x] Task 2.1: `src/simplex/lu_factor.h` & `src/simplex/lu_factor.cpp` (Markowitz ordering, threshold pivoting, FTRAN, BTRAN, iterative refinement)
- [x] Task 2.2: `src/simplex/lu_update.h` & `src/simplex/lu_update.cpp` (Eta updates, refactorization heuristic)

### Epic 3: Simplex Infrastructure (Basis & Pricing)
- [x] Task 3.1: `src/simplex/basis.h` & `src/simplex/basis.cpp` (Basis management, triangular crash procedure)
- [x] Task 3.2: `src/simplex/pricing.h` & `src/simplex/pricing.cpp` (Dantzig, Steepest Edge, Devex, Bland)

### Epic 4: Simplex Solvers
- [x] Task 4.1: `src/simplex/lp_solver.h` (LPSolver abstract base class)
- [x] Task 4.2: `src/simplex/primal_simplex.h` & `src/simplex/primal_simplex.cpp` (Revised Primal Simplex with bounded variables, ratio test, Phase 1/Phase 2)
- [x] Task 4.3: `src/simplex/dual_simplex.h` & `src/simplex/dual_simplex.cpp` (Revised Dual Simplex with Harris ratio test, warm start `solve_from_basis`)

### Epic 5: Build System, Unit Tests & Graphify Setup
- [x] Task 5.1: `CMakeLists.txt` for clean C++20 build with `-Wall -Wextra -Werror`
- [x] Task 5.2: Unit test suite covering all modules and edge cases (13/13 passing)
- [x] Task 5.3: Graphify setup and generation of interactive visualizer and reports
