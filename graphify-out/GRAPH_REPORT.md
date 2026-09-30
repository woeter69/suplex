# Graph Report - suplex  (2026-09-30)

## Corpus Check
- 104 files · ~73,808 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 13 file(s) not represented in the graph (top: .lp 4, (none) 3, .mps 3)

## Summary
- 1105 nodes · 2639 edges · 53 communities (40 shown, 13 thin omitted)
- Extraction: 76% EXTRACTED · 24% INFERRED · 0% AMBIGUOUS · INFERRED: 633 edges (avg confidence: 0.81)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Core Problem Model
- Shared Runtime Dependencies
- Branching Strategy Tests
- Cut Generation Tests
- C API Layer
- GPU PCG Kernels
- Conflict Analysis
- Command Line Interface
- Sparse Matrix Interface
- Sparse LU Decomposition
- LP Format Parser
- Project Architecture
- Presolve Engine
- Pricing Strategies
- Basis Management
- Primal Simplex State
- Team Workstream Design
- Branch and Bound
- MILP Solver Interface
- Dual Simplex State
- C++ Solver API
- LU Boundary Tests
- Primal Heuristics
- Graphify Workflows
- Solver Integration Layer
- Sparse Matrix Tests
- Solution Data Model
- LU Factorization
- LU Update Operations
- Sparse Matrix Operations
- Presolve Test Doubles
- MPS Parser Implementation
- Presolve Rule Tests
- Dual Simplex Operations
- Problem Reader Interfaces
- Primal Simplex Operations
- MILP Callback API
- MILP Unit Harness
- Heuristic Unit Tests
- Branch Bound Solve Flow
- Search Tree Statistics
- Solution Writer
- Algorithm Dispatch
- Solution Model Tests
- Node Pool Management
- Cross Repository Merge
- Media Transcription
- Maros Dataset Download
- MIPLIB Dataset Download
- Netlib Dataset Download
- Benchmark Integration Tests
- CLI Integration Tests
- Script Integration Tests

## God Nodes (most connected - your core abstractions)
1. `Problem` - 110 edges
2. `Solution` - 69 edges
3. `add_entry` - 44 edges
4. `to_csc` - 43 edges
5. `BranchBound` - 42 edges
6. `PrimalSimplexSolver` - 41 edges
7. `SuplexSolver` - 38 edges
8. `set_col_bounds` - 38 edges
9. `Suplex` - 37 edges
10. `SparseMatrixCSC` - 37 edges

## Surprising Connections (you probably didn't know these)
- `test_csr_out_of_bounds_get()` --calls--> `get`  [INFERRED]
  tests/unit/test_boundary.cpp → src/core/sparse_matrix.h
- `test_csr_empty_multiply()` --calls--> `multiply_vec`  [INFERRED]
  tests/unit/test_boundary.cpp → src/core/sparse_matrix.h
- `test_mps_and_lp_equivalence()` --calls--> `row_names_`  [INFERRED]
  tests/integration/test_integration.cpp → src/io/lp_parser.h
- `test_lu_zero_dimension()` --calls--> `factorize`  [INFERRED]
  tests/unit/test_boundary.cpp → src/simplex/lu_factor.h
- `Workstream Verification State` --semantically_similar_to--> `Person 3 Integration Readiness Gaps`  [INFERRED] [semantically similar]
  tests/TEST_MATRIX.md → PERSON3_NOTES.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Graphify Persistent Workflow** — _codex_skills_graphify_skill_full_pipeline, _codex_skills_graphify_references_update_incremental_graph_update, _codex_skills_graphify_references_query_graph_work_memory [INFERRED 0.85]
- **Suplex Solver Layers** — _context_architecture_suplex_architecture, _context_interfaces_shared_interfaces, _context_person1_simplex_sparse_linear_algebra_and_simplex [EXTRACTED 1.00]
- **Simplex Numerical Core** — _context_person1_simplex_sparse_lu_factorization, _context_person1_simplex_forrest_tomlin_update, _context_person1_simplex_revised_primal_simplex, _context_person1_simplex_revised_dual_simplex [EXTRACTED 1.00]
- **Solver Workstream Integration** — _context_person2_interior_point_person_2_workstream, _context_person3_milp_person_3_workstream, _context_person4_integration_person_4_workstream [EXTRACTED 1.00]
- **Release Readiness Evidence** — changelog_person_3_implementation_status, person3_notes_integration_readiness_gaps, tests_test_matrix_verification_state, tests_test_matrix_combined_release_gates [INFERRED 0.85]

## Communities (53 total, 13 thin omitted)

### Community 0 - "Core Problem Model"
Cohesion: 0.07
Nodes (74): Problem, A_, c_, col_lower_, col_upper_, is_mip, num_integers, Q_ (+66 more)

### Community 1 - "Shared Runtime Dependencies"
Cohesion: 0.05
Nodes (23): to_string(), BranchBound::BranchBound(), to_string(), BenchmarkResult, columns, gap, iterations, name (+15 more)

### Community 2 - "Branching Strategy Tests"
Cohesion: 0.07
Nodes (30): Branching, Branching::Branching(), most_fractional, pseudo_cost_branch, pseudo_costs_, record_observation, reliability_branch, reliability_threshold_ (+22 more)

### Community 3 - "Cut Generation Tests"
Cohesion: 0.10
Nodes (30): separate_cuts, Cut, age, coefficients, efficacy, indices, name, rhs (+22 more)

### Community 4 - "C API Layer"
Cohesion: 0.11
Nodes (38): copy_result(), invalid(), suplex_create(), suplex_destroy(), suplex_get_dual_values(), suplex_get_error(), suplex_get_iterations(), suplex_get_mip_gap() (+30 more)

### Community 5 - "GPU PCG Kernels"
Cohesion: 0.07
Nodes (22): apply_jacobi(), cuda_pcg_host(), device_dot(), dot_kernel(), initialize_residual(), pcg_spmv(), update_direction(), update_x_r() (+14 more)

### Community 6 - "Conflict Analysis"
Cohesion: 0.09
Nodes (24): ConflictAnalysis, analyse, build_nogood_cut, minimize_conflict, responsible_cols, BBNode, basis, bound_changes (+16 more)

### Community 7 - "Command Line Interface"
Cohesion: 0.06
Nodes (15): Options, algorithm, gpu, log_level, mip_gap, presolve, problem_file, solution_file (+7 more)

### Community 8 - "Sparse Matrix Interface"
Cohesion: 0.07
Nodes (27): SparseMatrixCSC, clear, col_start, nnz, num_cols, num_rows, row_index, to_csr (+19 more)

### Community 9 - "Sparse LU Decomposition"
Cohesion: 0.07
Nodes (20): SparseLU, condition_est_, L_rows_, m_, p_, p_inv_, q_, q_inv_ (+12 more)

### Community 10 - "LP Format Parser"
Cohesion: 0.12
Nodes (20): bound_words(), lower(), LPReader, col_names_, error_, read, row_names_, parse_linear_expression() (+12 more)

### Community 11 - "Project Architecture"
Cohesion: 0.08
Nodes (27): Sovereign Optimization Solver, Suplex Mission, Layered Solver Architecture, Dual Simplex Warm Start Contract, C++20 Build Environment, Double Precision and Sparse Indexing, Suplex Repository Map, Simplex Module (+19 more)

### Community 12 - "Presolve Engine"
Cohesion: 0.07
Nodes (24): PresolveRecord, a_ij, affected_cols, c_j, col_idx, fix_value, obj_offset, row_idx (+16 more)

### Community 13 - "Pricing Strategies"
Cohesion: 0.12
Nodes (14): Pricing, init, iterations_since_reset_, num_rows_, reset_weights, select_entering, total_vars_, update (+6 more)

### Community 14 - "Basis Management"
Cohesion: 0.14
Nodes (11): Basis, basic_row_, basic_vars_, crash, init, pivot, set_logical_basis, set_status (+3 more)

### Community 15 - "Primal Simplex State"
Cohesion: 0.09
Nodes (16): PrimalSimplexSolver, A_aug_, basis_, c_aug_, d_, log_level_, lower_aug_, lu_manager_ (+8 more)

### Community 16 - "Team Workstream Design"
Cohesion: 0.13
Nodes (19): Mehrotra Predictor Corrector IPM, Person 2 Interior Point Workstream, Sparse Cholesky and Crossover, Branch Cut MILP Engine, Person 3 Presolve and MILP Workstream, GPU and Benchmark Infrastructure, I O API and Solver Orchestration, Person 4 Integration Workstream (+11 more)

### Community 17 - "Branch and Bound"
Cohesion: 0.09
Nodes (16): BranchBound, callback_, cut_rounds_node_, cut_rounds_root_, gap_tol_, incumbent_obj_, incumbent_sol_, is_integer_feasible (+8 more)

### Community 18 - "MILP Solver Interface"
Cohesion: 0.13
Nodes (11): solve, MILPSolver, callback_, gap_tol_, log_level_, lp_solver_, node_limit_, set_callback (+3 more)

### Community 19 - "Dual Simplex State"
Cohesion: 0.11
Nodes (14): DualSimplexSolver, A_aug_, basis_, c_aug_, d_, log_level_, lower_aug_, lu_manager_ (+6 more)

### Community 20 - "C++ Solver API"
Cohesion: 0.17
Nodes (14): main(), Suplex, algorithm_, column_names_, error_, gpu_, log_level_, mip_gap_ (+6 more)

### Community 21 - "LU Boundary Tests"
Cohesion: 0.16
Nodes (12): from_triplets, get, test_csc_1x1_matrix(), test_csc_from_triplets_no_entries(), test_csc_identity_matrix(), test_csc_multiple_duplicates_same_entry(), test_csc_out_of_bounds_get(), test_csr_empty_multiply() (+4 more)

### Community 22 - "Primal Heuristics"
Cohesion: 0.25
Nodes (13): run_heuristics, Heuristics, coefficient_diving, feasibility_pump, fractional_diving, is_feasible, is_integer_feasible, log_level_ (+5 more)

### Community 23 - "Graphify Workflows"
Cohesion: 0.12
Nodes (15): Incremental Folder Watcher, URL Ingestion and Folder Watch, Graph Export Formats, Token Reduction Benchmark, Extraction Confidence Rubric, Semantic Extraction Contract, Automatic AST Graph Refresh, Graphify Commit Hook Integration (+7 more)

### Community 24 - "Solver Integration Layer"
Cohesion: 0.22
Nodes (7): solution_, fixture(), test_milp_end_to_end(), test_milp_fractional_root_branches(), test_missing_ipm_reports_error(), test_mps_and_lp_equivalence(), test_top_level_file_solve()

### Community 25 - "Sparse Matrix Tests"
Cohesion: 0.22
Nodes (14): add_entry, to_csc, test_csc_clear(), test_csc_dot_col_out_of_range(), test_csc_multiply_col_out_of_range(), test_csc_multiply_col_zero_scalar(), test_csc_spmv_zero_vector(), test_lu_dimension_mismatch() (+6 more)

### Community 26 - "Solution Data Model"
Cohesion: 0.13
Nodes (12): Solution, best_bound, dual_values, mip_gap, num_iterations, num_nodes, objective_value, primal_values (+4 more)

### Community 27 - "LU Factorization"
Cohesion: 0.29
Nodes (9): btran, factorize, factorize_matrix, ftran, iterative_refine, test_lu_1x1_factorize(), test_lu_identity_factorize(), test_lu_btran() (+1 more)

### Community 28 - "LU Update Operations"
Cohesion: 0.29
Nodes (10): btran, ftran, refactorize, should_refactorize, update, test_lu_btran_after_update(), test_lu_should_refactorize(), test_lu_update_near_zero_pivot() (+2 more)

### Community 29 - "Sparse Matrix Operations"
Cohesion: 0.26
Nodes (8): dot_col, multiply_col, multiply_transpose_vec, multiply_vec, get, multiply_vec, reserve, pcg_solve

### Community 30 - "Presolve Test Doubles"
Cohesion: 0.21
Nodes (3): MockLPSolver, fixed_sol, fixed_status

### Community 31 - "MPS Parser Implementation"
Cohesion: 0.29
Nodes (5): fields(), read, number(), unquote(), upper()

### Community 33 - "Dual Simplex Operations"
Cohesion: 0.24
Nodes (7): compute_basic_primal, compute_duals_and_reduced_costs, ensure_dual_feasibility, extract_solution, setup_augmented_system, solve, solve_from_basis

### Community 34 - "Problem Reader Interfaces"
Cohesion: 0.29
Nodes (7): read_lp, read_mps, MPSReader, col_names_, error_, row_names_, test_parser_error_diagnostic()

### Community 35 - "Primal Simplex Operations"
Cohesion: 0.22
Nodes (4): check_primal_feasibility, compute_basic_solution, compute_duals_and_reduced_costs, setup_augmented_system

### Community 39 - "Branch Bound Solve Flow"
Cohesion: 0.33
Nodes (6): compute_gap, log_progress, solve, solve_root_lp, try_update_incumbent, now_seconds()

### Community 40 - "Search Tree Statistics"
Cohesion: 0.22
Nodes (9): TreeStats, cuts_added, heuristic_sols, max_depth, nodes_created, nodes_infeasible, nodes_per_second, nodes_pruned (+1 more)

### Community 41 - "Solution Writer"
Cohesion: 0.29
Nodes (4): SolutionWriter, error_, write, test_solution_writer()

### Community 43 - "Solution Model Tests"
Cohesion: 0.40
Nodes (4): is_feasible, reset, test_solution_is_feasible(), test_solution_reset()

### Community 44 - "Node Pool Management"
Cohesion: 0.50
Nodes (3): apply_bound_changes, remove_from_pool, select_node

## Knowledge Gaps
- **254 isolated node(s):** `algorithm`, `log_level`, `presolve`, `gpu`, `stats` (+249 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 411 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **13 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Problem` connect `Core Problem Model` to `Shared Runtime Dependencies`, `Branching Strategy Tests`, `Cut Generation Tests`, `C API Layer`, `Conflict Analysis`, `Sparse Matrix Interface`, `LP Format Parser`, `Presolve Engine`, `Branch and Bound`, `C++ Solver API`, `LU Boundary Tests`, `Primal Heuristics`, `Sparse Matrix Tests`, `Presolve Test Doubles`, `MPS Parser Implementation`, `Dual Simplex Operations`, `Primal Simplex Operations`, `MILP Unit Harness`, `Heuristic Unit Tests`, `Branch Bound Solve Flow`, `Algorithm Dispatch`, `Node Pool Management`?**
  _High betweenness centrality (0.252) - this node is a cross-community bridge._
- **Why does `Solution` connect `Solution Data Model` to `Core Problem Model`, `Shared Runtime Dependencies`, `Branching Strategy Tests`, `Cut Generation Tests`, `MILP Callback API`, `Dual Simplex Operations`, `MILP Unit Harness`, `Branch Bound Solve Flow`, `Heuristic Unit Tests`, `Solution Writer`, `Solution Model Tests`, `Branch and Bound`, `C++ Solver API`, `LU Boundary Tests`, `Primal Heuristics`, `Presolve Test Doubles`?**
  _High betweenness centrality (0.110) - this node is a cross-community bridge._
- **Why does `SparseMatrixCSC` connect `Sparse Matrix Interface` to `Core Problem Model`, `Shared Runtime Dependencies`, `GPU PCG Kernels`, `Basis Management`, `Primal Simplex State`, `Dual Simplex State`, `LU Boundary Tests`, `Sparse Matrix Tests`, `LU Factorization`, `LU Update Operations`, `Sparse Matrix Operations`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `Problem` (e.g. with `main()` and `test_lp_compact_bounds_and_types()`) actually correct?**
  _`Problem` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 40 inferred relationships involving `add_entry` (e.g. with `test_csc_clear()` and `test_csc_dot_col_out_of_range()`) actually correct?**
  _`add_entry` has 40 INFERRED edges - model-reasoned connections that need verification._
- **Are the 39 inferred relationships involving `to_csc` (e.g. with `test_csc_clear()` and `test_csc_dot_col_out_of_range()`) actually correct?**
  _`to_csc` has 39 INFERRED edges - model-reasoned connections that need verification._
- **What connects `algorithm`, `log_level`, `presolve` to the rest of the system?**
  _254 weakly-connected nodes found - possible documentation gaps or missing edges._