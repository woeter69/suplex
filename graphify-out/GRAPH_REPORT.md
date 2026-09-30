# Graph Report - suplex  (2026-09-30)

## Corpus Check
- 103 files · ~73,953 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 13 file(s) not represented in the graph (top: .lp 4, (none) 3, .mps 3)

## Summary
- 1106 nodes · 2647 edges · 49 communities (39 shown, 10 thin omitted)
- Extraction: 76% EXTRACTED · 24% INFERRED · 0% AMBIGUOUS · INFERRED: 633 edges (avg confidence: 0.81)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `b7bab068`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- Problem
- vector
- TEST
- TEST
- suplex_api.cpp
- GPUManager
- BBNode
- Options
- SparseMatrixCSR
- LUManager
- lp_parser.cpp
- Sparse Linear Algebra and Simplex Engine
- set_col_bounds
- Pricing
- Basis
- PrimalSimplexSolver
- Branch Cut MILP Engine
- BranchBound
- MILPSolver
- DualSimplexSolver
- Suplex
- problem.cpp
- Solution
- Graphify Full Pipeline
- test_integration.cpp
- test_boundary.cpp
- LPReader
- SparseLU
- BenchmarkResult
- SparseMatrixCSC
- Algorithm
- mps_parser.cpp
- bindings.cpp
- dual_simplex.cpp
- suplex_solver.cpp
- primal_simplex.cpp
- MILPCallback
- solve
- TreeStats
- write
- apply_bound_changes
- Cross Repository Graph Merge
- Domain Prompted Whisper Transcription
- download_maros.sh
- download_miplib.sh
- download_netlib.sh
- test_benchmark.sh
- test_cli.sh
- test_scripts.sh

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
- `Workstream Verification State` --semantically_similar_to--> `Person 3 Integration Readiness Gaps`  [INFERRED] [semantically similar]
  tests/TEST_MATRIX.md → PERSON3_NOTES.md
- `main()` --calls--> `column_names_`  [INFERRED]
  cli/main.cpp → src/api/suplex_solver.h
- `main()` --calls--> `read_lp`  [INFERRED]
  cli/main.cpp → src/api/suplex_solver.h
- `main()` --calls--> `read_mps`  [INFERRED]
  cli/main.cpp → src/api/suplex_solver.h
- `main()` --calls--> `row_names_`  [INFERRED]
  cli/main.cpp → src/api/suplex_solver.h

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Simplex Numerical Core** — _context_person1_simplex_sparse_lu_factorization, _context_person1_simplex_forrest_tomlin_update, _context_person1_simplex_revised_primal_simplex, _context_person1_simplex_revised_dual_simplex [EXTRACTED 1.00]
- **Solver Workstream Integration** — _context_person2_interior_point_person_2_workstream, _context_person3_milp_person_3_workstream, _context_person4_integration_person_4_workstream [EXTRACTED 1.00]
- **Suplex Solver Layers** — _context_architecture_suplex_architecture, _context_interfaces_shared_interfaces, _context_person1_simplex_sparse_linear_algebra_and_simplex [EXTRACTED 1.00]
- **Graphify Persistent Workflow** — _codex_skills_graphify_skill_full_pipeline, _codex_skills_graphify_references_update_incremental_graph_update, _codex_skills_graphify_references_query_graph_work_memory [INFERRED 0.85]
- **Release Readiness Evidence** — changelog_person_3_implementation_status, person3_notes_integration_readiness_gaps, tests_test_matrix_verification_state, tests_test_matrix_combined_release_gates [INFERRED 0.85]

## Communities (49 total, 10 thin omitted)

### Community 0 - "Problem"
Cohesion: 0.05
Nodes (63): Problem, A_, c_, col_lower_, col_upper_, Q_, row_lower_, row_upper_ (+55 more)

### Community 1 - "vector"
Cohesion: 0.10
Nodes (7): to_string(), BranchBound::BranchBound(), to_string(), register_test(), TestCase, func, name

### Community 2 - "TEST"
Cohesion: 0.06
Nodes (30): Branching, Branching::Branching(), most_fractional, pseudo_cost_branch, pseudo_costs_, record_observation, reliability_branch, reliability_threshold_ (+22 more)

### Community 3 - "TEST"
Cohesion: 0.10
Nodes (30): separate_cuts, Cut, age, coefficients, efficacy, indices, name, rhs (+22 more)

### Community 4 - "suplex_api.cpp"
Cohesion: 0.11
Nodes (38): copy_result(), invalid(), suplex_create(), suplex_destroy(), suplex_get_dual_values(), suplex_get_error(), suplex_get_iterations(), suplex_get_mip_gap() (+30 more)

### Community 5 - "GPUManager"
Cohesion: 0.07
Nodes (23): apply_jacobi(), cuda_pcg_host(), device_dot(), dot_kernel(), initialize_residual(), pcg_spmv(), update_direction(), update_x_r() (+15 more)

### Community 6 - "BBNode"
Cohesion: 0.09
Nodes (24): ConflictAnalysis, analyse, build_nogood_cut, minimize_conflict, responsible_cols, BBNode, basis, bound_changes (+16 more)

### Community 7 - "Options"
Cohesion: 0.15
Nodes (11): Options, algorithm, gpu, log_level, mip_gap, presolve, problem_file, solution_file (+3 more)

### Community 8 - "SparseMatrixCSR"
Cohesion: 0.10
Nodes (16): SparseMatrixCSR, clear, col_index, nnz, num_cols, num_rows, row_start, to_csc (+8 more)

### Community 9 - "LUManager"
Cohesion: 0.12
Nodes (18): EtaMatrix, indices, leaving_pos, pivot_val, values, LUManager, base_lu_, btran (+10 more)

### Community 10 - "lp_parser.cpp"
Cohesion: 0.28
Nodes (13): bound_words(), lower(), read, parse_linear_expression(), parse_number(), ParsedRow, name, relation (+5 more)

### Community 11 - "Sparse Linear Algebra and Simplex Engine"
Cohesion: 0.08
Nodes (27): Sovereign Optimization Solver, Suplex Mission, Layered Solver Architecture, Dual Simplex Warm Start Contract, C++20 Build Environment, Double Precision and Sparse Indexing, Suplex Repository Map, Simplex Module (+19 more)

### Community 12 - "set_col_bounds"
Cohesion: 0.27
Nodes (25): set_col_bounds, set_constraint_matrix, set_objective, set_objective_coeff, set_row_bounds, solve, make_simple_mip(), make_problem() (+17 more)

### Community 13 - "Pricing"
Cohesion: 0.12
Nodes (14): Pricing, init, iterations_since_reset_, num_rows_, reset_weights, select_entering, total_vars_, update (+6 more)

### Community 14 - "Basis"
Cohesion: 0.14
Nodes (11): Basis, basic_row_, basic_vars_, crash, init, pivot, set_logical_basis, set_status (+3 more)

### Community 15 - "PrimalSimplexSolver"
Cohesion: 0.09
Nodes (16): PrimalSimplexSolver, A_aug_, basis_, c_aug_, d_, log_level_, lower_aug_, lu_manager_ (+8 more)

### Community 16 - "Branch Cut MILP Engine"
Cohesion: 0.13
Nodes (19): Mehrotra Predictor Corrector IPM, Person 2 Interior Point Workstream, Sparse Cholesky and Crossover, Branch Cut MILP Engine, Person 3 Presolve and MILP Workstream, GPU and Benchmark Infrastructure, I O API and Solver Orchestration, Person 4 Integration Workstream (+11 more)

### Community 17 - "BranchBound"
Cohesion: 0.09
Nodes (16): BranchBound, callback_, cut_rounds_node_, cut_rounds_root_, gap_tol_, incumbent_obj_, incumbent_sol_, is_integer_feasible (+8 more)

### Community 18 - "MILPSolver"
Cohesion: 0.13
Nodes (11): solve, MILPSolver, callback_, gap_tol_, log_level_, lp_solver_, node_limit_, set_callback (+3 more)

### Community 19 - "DualSimplexSolver"
Cohesion: 0.11
Nodes (14): DualSimplexSolver, A_aug_, basis_, c_aug_, d_, log_level_, lower_aug_, lu_manager_ (+6 more)

### Community 20 - "Suplex"
Cohesion: 0.13
Nodes (18): main(), Suplex, algorithm_, column_names_, error_, gpu_, log_level_, mip_gap_ (+10 more)

### Community 21 - "problem.cpp"
Cohesion: 0.15
Nodes (11): is_mip, num_integers, set_dimensions, set_quadratic_objective, set_var_types, validate, make_knapsack(), test_problem_is_mip() (+3 more)

### Community 22 - "Solution"
Cohesion: 0.06
Nodes (34): Solution, best_bound, dual_values, is_feasible, mip_gap, num_iterations, num_nodes, objective_value (+26 more)

### Community 23 - "Graphify Full Pipeline"
Cohesion: 0.12
Nodes (15): Incremental Folder Watcher, URL Ingestion and Folder Watch, Graph Export Formats, Token Reduction Benchmark, Extraction Confidence Rubric, Semantic Extraction Contract, Automatic AST Graph Refresh, Graphify Commit Hook Integration (+7 more)

### Community 24 - "test_integration.cpp"
Cohesion: 0.24
Nodes (3): parse_options(), parse_switch(), usage()

### Community 25 - "test_boundary.cpp"
Cohesion: 0.15
Nodes (23): add_entry, to_csc, factorize_matrix, test_csc_clear(), test_csc_dot_col_out_of_range(), test_csc_identity_matrix(), test_csc_multiply_col_out_of_range(), test_csc_multiply_col_zero_scalar() (+15 more)

### Community 26 - "LPReader"
Cohesion: 0.15
Nodes (8): LPReader, col_names_, error_, row_names_, ProblemReader, last_error, read, test_mps_and_lp_equivalence()

### Community 27 - "SparseLU"
Cohesion: 0.11
Nodes (15): SparseLU, btran, condition_est_, factorize, ftran, iterative_refine, L_rows_, m_ (+7 more)

### Community 28 - "BenchmarkResult"
Cohesion: 0.13
Nodes (13): BenchmarkResult, columns, gap, iterations, name, nodes, nonzeros, objective (+5 more)

### Community 29 - "SparseMatrixCSC"
Cohesion: 0.12
Nodes (26): SparseMatrixCSC, clear, col_start, dot_col, from_triplets, get, multiply_col, multiply_transpose_vec (+18 more)

### Community 30 - "Algorithm"
Cohesion: 0.24
Nodes (4): main(), read_references(), write_csv(), write_json()

### Community 31 - "mps_parser.cpp"
Cohesion: 0.23
Nodes (5): fields(), read, number(), unquote(), upper()

### Community 33 - "dual_simplex.cpp"
Cohesion: 0.24
Nodes (7): compute_basic_primal, compute_duals_and_reduced_costs, ensure_dual_feasibility, extract_solution, setup_augmented_system, solve, solve_from_basis

### Community 34 - "suplex_solver.cpp"
Cohesion: 0.17
Nodes (8): read_lp, read_mps, set_problem, MPSReader, col_names_, error_, row_names_, test_parser_error_diagnostic()

### Community 35 - "primal_simplex.cpp"
Cohesion: 0.21
Nodes (6): check_primal_feasibility, compute_basic_solution, compute_duals_and_reduced_costs, extract_solution, setup_augmented_system, solve_from_basis

### Community 39 - "solve"
Cohesion: 0.33
Nodes (6): compute_gap, log_progress, solve, solve_root_lp, try_update_incumbent, now_seconds()

### Community 40 - "TreeStats"
Cohesion: 0.22
Nodes (9): TreeStats, cuts_added, heuristic_sols, max_depth, nodes_created, nodes_infeasible, nodes_per_second, nodes_pruned (+1 more)

### Community 41 - "write"
Cohesion: 0.29
Nodes (4): SolutionWriter, error_, write, test_solution_writer()

### Community 44 - "apply_bound_changes"
Cohesion: 0.50
Nodes (3): apply_bound_changes, remove_from_pool, select_node

## Knowledge Gaps
- **254 isolated node(s):** `algorithm`, `log_level`, `presolve`, `gpu`, `stats` (+249 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 409 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **10 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Problem` connect `Problem` to `vector`, `TEST`, `TEST`, `suplex_api.cpp`, `BBNode`, `lp_parser.cpp`, `set_col_bounds`, `BranchBound`, `Suplex`, `problem.cpp`, `Solution`, `test_boundary.cpp`, `LPReader`, `SparseMatrixCSC`, `mps_parser.cpp`, `dual_simplex.cpp`, `suplex_solver.cpp`, `primal_simplex.cpp`, `solve`, `apply_bound_changes`?**
  _High betweenness centrality (0.251) - this node is a cross-community bridge._
- **Why does `Solution` connect `Solution` to `Problem`, `vector`, `TEST`, `TEST`, `MILPCallback`, `dual_simplex.cpp`, `primal_simplex.cpp`, `solve`, `write`, `set_col_bounds`, `BranchBound`, `Suplex`, `test_boundary.cpp`?**
  _High betweenness centrality (0.109) - this node is a cross-community bridge._
- **Why does `SparseMatrixCSC` connect `SparseMatrixCSC` to `Problem`, `vector`, `GPUManager`, `SparseMatrixCSR`, `LUManager`, `set_col_bounds`, `Basis`, `PrimalSimplexSolver`, `DualSimplexSolver`, `problem.cpp`, `test_boundary.cpp`, `SparseLU`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `Problem` (e.g. with `main()` and `test_lp_compact_bounds_and_types()`) actually correct?**
  _`Problem` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 40 inferred relationships involving `add_entry` (e.g. with `test_csc_clear()` and `test_csc_dot_col_out_of_range()`) actually correct?**
  _`add_entry` has 40 INFERRED edges - model-reasoned connections that need verification._
- **Are the 39 inferred relationships involving `to_csc` (e.g. with `test_csc_clear()` and `test_csc_dot_col_out_of_range()`) actually correct?**
  _`to_csc` has 39 INFERRED edges - model-reasoned connections that need verification._
- **What connects `algorithm`, `log_level`, `presolve` to the rest of the system?**
  _254 weakly-connected nodes found - possible documentation gaps or missing edges._