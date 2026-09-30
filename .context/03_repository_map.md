# Repository Map

## Directory Structure & Module Descriptions
- `.context/`: Project documentation, architecture guidelines, ADRs, state matrices, and task descriptions.
- `CMakeLists.txt`: Root build configuration for core libraries and test suite.
- `src/core/`: Common foundational data structures.
  - `types.h`: Universal numerical aliases (`Real`, `Index`), tolerances, enums (`SolverStatus`, `VarType`).
  - `sparse_matrix.h`, `sparse_matrix.cpp`: CSC, CSR, and Triplet matrix representations with SpMV, transpose, and column access.
  - `problem.h`, `problem.cpp`: General bounded/ranged linear and quadratic optimization problem container.
  - `solution.h`, `solution.cpp`: Optimization solution container with primal/dual values, reduced costs, status, and metrics.
- `src/simplex/`: Simplex solver engine (Person 1 ownership).
  - `lp_solver.h`: Abstract base class for all LP solvers.
  - `basis.h`, `basis.cpp`: Basis tracking, variable statuses, and triangular crash procedures.
  - `lu_factor.h`, `lu_factor.cpp`: Sparse LU factorization with Markowitz ordering, threshold pivoting, and FTRAN/BTRAN.
  - `lu_update.h`, `lu_update.cpp`: Forrest-Tomlin / product-form eta updates with refactorization thresholds.
  - `pricing.h`, `pricing.cpp`: Pricing strategies (Dantzig, Steepest Edge, Devex, Bland).
  - `primal_simplex.h`, `primal_simplex.cpp`: Revised Primal Simplex algorithm with bounded variables and Phase 1/Phase 2.
  - `dual_simplex.h`, `dual_simplex.cpp`: Revised Dual Simplex algorithm with Harris ratio test and warm-start support.
- `tests/unit/`: Comprehensive test suite for sparse matrices, LU factorization, basis, pricing, primal and dual simplex.
