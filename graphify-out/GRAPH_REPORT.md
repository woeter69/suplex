# Graph Report - suplex  (2026-09-30)

## Corpus Check
- Corpus is ~31,556 words - fits in a single context window. You may not need a graph.

## Summary
- 331 nodes · 656 edges · 10 communities
- Extraction: 86% EXTRACTED · 14% INFERRED · 0% AMBIGUOUS · INFERRED: 90 edges (avg confidence: 0.81)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Sparse Matrix Formats
- Optimization Problem Representation
- Primal Simplex & Solver Interface
- Dual Simplex Algorithm
- Sparse LU Factorization
- LU Updates & Eta Management
- Simplex Basis & Crash Procedure
- Core Types & Problem Validation
- Pricing Strategies & Weights
- Solution Representation & Metrics

## God Nodes (most connected - your core abstractions)
1. `Problem` - 42 edges
2. `PrimalSimplexSolver` - 40 edges
3. `SparseMatrixCSC` - 35 edges
4. `DualSimplexSolver` - 33 edges
5. `Solution` - 25 edges
6. `SparseLU` - 24 edges
7. `LUManager` - 22 edges
8. `Basis` - 21 edges
9. `SparseMatrixCSR` - 19 edges
10. `Pricing` - 17 edges

## Surprising Connections (you probably didn't know these)
- `test_lu_btran()` --calls--> `add_entry`  [INFERRED]
  tests/unit/test_lu.cpp → src/core/sparse_matrix.h
- `test_lu_singular()` --calls--> `add_entry`  [INFERRED]
  tests/unit/test_lu.cpp → src/core/sparse_matrix.h
- `test_lu_small()` --calls--> `add_entry`  [INFERRED]
  tests/unit/test_lu.cpp → src/core/sparse_matrix.h
- `test_lu_update()` --calls--> `add_entry`  [INFERRED]
  tests/unit/test_lu.cpp → src/core/sparse_matrix.h
- `test_dual_simplex_tiny()` --calls--> `add_entry`  [INFERRED]
  tests/unit/test_simplex.cpp → src/core/sparse_matrix.h

## Import Cycles
- None detected.

## Communities (10 total, 0 thin omitted)

### Community 0 - "Sparse Matrix Formats"
Cohesion: 0.07
Nodes (46): Index, Real, vector, Index, Real, vector, SparseMatrixCSC, clear (+38 more)

### Community 1 - "Optimization Problem Representation"
Cohesion: 0.11
Nodes (35): ObjectiveSense, optional, Index, Real, VarType, vector, Index, Real (+27 more)

### Community 2 - "Primal Simplex & Solver Interface"
Cohesion: 0.06
Nodes (37): LPSolver, set_iteration_limit, set_log_level, set_time_limit, solve, solve_from_basis, Index, Real (+29 more)

### Community 3 - "Dual Simplex Algorithm"
Cohesion: 0.09
Nodes (28): Index, SolverStatus, vector, DualSimplexSolver, A_aug_, basis_, c_aug_, compute_basic_primal (+20 more)

### Community 4 - "Sparse LU Factorization"
Cohesion: 0.10
Nodes (26): pair, Index, Real, SolverStatus, vector, Index, Real, vector (+18 more)

### Community 5 - "LU Updates & Eta Management"
Cohesion: 0.12
Nodes (23): Index, Real, SolverStatus, vector, EtaMatrix, indices, leaving_pos, pivot_val (+15 more)

### Community 6 - "Simplex Basis & Crash Procedure"
Cohesion: 0.14
Nodes (17): Basis, basic_row_, basic_vars_, crash, init, pivot, set_logical_basis, set_status (+9 more)

### Community 7 - "Core Types & Problem Validation"
Cohesion: 0.16
Nodes (10): function, vector, SolverStatus, string, to_string(), string, register_test(), TestCase (+2 more)

### Community 8 - "Pricing Strategies & Weights"
Cohesion: 0.14
Nodes (19): BasisStatus, Index, PricingStrategy, Real, vector, Index, PricingStrategy, Real (+11 more)

### Community 9 - "Solution Representation & Metrics"
Cohesion: 0.12
Nodes (17): Real, SolverStatus, vector, Solution, best_bound, dual_values, is_feasible, mip_gap (+9 more)

## Knowledge Gaps
- **95 isolated node(s):** `A_`, `c_`, `col_lower_`, `col_upper_`, `row_lower_` (+90 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 122 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PrimalSimplexSolver` connect `Primal Simplex & Solver Interface` to `Sparse Matrix Formats`, `Optimization Problem Representation`, `Dual Simplex Algorithm`, `LU Updates & Eta Management`, `Simplex Basis & Crash Procedure`, `Core Types & Problem Validation`, `Pricing Strategies & Weights`?**
  _High betweenness centrality (0.340) - this node is a cross-community bridge._
- **Why does `SparseMatrixCSC` connect `Sparse Matrix Formats` to `Optimization Problem Representation`, `Primal Simplex & Solver Interface`, `Dual Simplex Algorithm`, `Sparse LU Factorization`, `LU Updates & Eta Management`, `Simplex Basis & Crash Procedure`, `Core Types & Problem Validation`?**
  _High betweenness centrality (0.308) - this node is a cross-community bridge._
- **Why does `DualSimplexSolver` connect `Dual Simplex Algorithm` to `Sparse Matrix Formats`, `Primal Simplex & Solver Interface`, `LU Updates & Eta Management`, `Simplex Basis & Crash Procedure`, `Core Types & Problem Validation`?**
  _High betweenness centrality (0.203) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `PrimalSimplexSolver` (e.g. with `solve` and `test_warm_start()`) actually correct?**
  _`PrimalSimplexSolver` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `A_`, `c_`, `col_lower_` to the rest of the system?**
  _95 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Sparse Matrix Formats` be split into smaller, more focused modules?**
  _Cohesion score 0.07138535995160314 - nodes in this community are weakly interconnected._
- **Should `Optimization Problem Representation` be split into smaller, more focused modules?**
  _Cohesion score 0.10782241014799154 - nodes in this community are weakly interconnected._