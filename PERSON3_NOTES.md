# Person 3 — Working Notes

> Worktree: `/mnt/d_drive/programs/suplex-person3`
> Branch: `person3/milp-engine`
> Spec: `.context/PERSON3_MILP.md`

---

## What's Done ✅

| File | Status | Notes |
|------|--------|-------|
| `src/core/types.h` | ✅ Done | All enums, Real/Index typedefs |
| `src/core/sparse_matrix.h` | ✅ Done | Interface stub |
| `src/core/problem.h` | ✅ Done | Full interface |
| `src/core/solution.h` | ✅ Done | Full struct |
| `src/core/core_stubs.cpp` | ✅ Done | Replace with Person 1's code at integration |
| `src/simplex/lp_solver.h` | ✅ Done | Interface stub (Person 1/2 implement) |
| `src/milp/presolve.h/.cpp` | ✅ Done | 9 rules, undo stack, postsolve, compact() |
| `src/milp/node.h/.cpp` | ✅ Done | BBNode, BoundChange, TreeStats |
| `src/milp/branch_bound.h/.cpp` | ✅ Done | Full B&B loop, 3 node selection strategies |
| `src/milp/branching.h/.cpp` | ✅ Done | Pseudo-costs, strong, reliability branching |
| `src/milp/cuts.h/.cpp` | ✅ Done | Gomory, MIR, cover (with lifting), clique |
| `src/milp/heuristics.h/.cpp` | ✅ Done | Rounding, diving, feas-pump, RINS |
| `src/milp/conflict.h/.cpp` | ✅ Done | Farkas-based conflict analysis |
| `src/milp/milp_solver.h/.cpp` | ✅ Done | Public API, presolve+B&B+postsolve |
| `CMakeLists.txt` | ✅ Done | cmake + GoogleTest |
| `Makefile` | ✅ Done | Direct g++ fallback |
| `tests/unit/milp/test_presolve.cpp` | ✅ Done | 9 tests |
| `tests/unit/milp/test_branching.cpp` | ✅ Done | 5 tests |
| `tests/unit/milp/test_cuts.cpp` | ✅ Done | 7 tests |
| `tests/unit/milp/test_heuristics.cpp` | ✅ Done | 7 tests |

**All 9 source files: 0 compiler errors** ✅

---

## What's Next 🔨

### High Priority (Implement now)

1. **Wire GoogleTest and run tests** — install cmake via `nix-env` or download binary,
   then `cmake -B build && cmake --build build && ctest`

2. **Presolve: duplicate row/col detection** — currently stubs returning `false`
   - Strategy: hash each row's sparsity pattern + coefficient ratios, find matches

3. **Branch-and-bound: cut insertion** — `separate_cuts()` currently counts but does
   NOT actually add cuts as new rows to the problem. Add:
   ```cpp
   // After generating selected cuts, rebuild problem_with_cuts:
   // Add each cut as a new row (lower=-INF, upper=rhs)
   ```

4. **Gomory cuts: tableau row extraction** — `extract_tableau_row()` returns `false`
   (stub). Needs Person 1 to expose `btran(e_i)` → dense vector. Add to `LPSolver`:
   ```cpp
   virtual void btran_row(Index basic_pos, std::vector<Real>& row) const = 0;
   ```

5. **Integration test: small MIP** — once Person 1 sends DualSimplexSolver,
   solve `p0033` from MIPLIB (17 rows, 33 cols, known optimal = 3089)

### Medium Priority (Week 3-6)

6. **Postsolve: exact dual recovery** — current postsolve approximates row duals;
   implement proper dual completion from KKT conditions

7. **Free column substitution** — presolve rule for free variables (lb=-INF):
   substitute from equality row, eliminate variable + row

8. **Coefficient aggregation** — aggregate rows to strengthen cuts further

9. **RINS: proper sub-MIP** — current RINS solves LP relaxation only; upgrade to
   call `MILPSolver::solve()` recursively with tight node limit (e.g., 500 nodes)

10. **Parallel node evaluation** — `BranchBound::process_node()` can be parallelised
    with OpenMP for independent nodes at same depth

### Integration Week (Week 7-8)

11. **Replace `core_stubs.cpp`** with Person 1's `sparse_matrix.cpp` + `problem.cpp`
12. **Set up LPSolver** — inject `DualSimplexSolver` from Person 1 into `MILPSolver`
13. **Run MIPLIB benchmark** against Person 3's targets (see `PERSON3_MILP.md §11`)

---

## Key Interfaces You Depend On

### From Person 1 (DualSimplexSolver)
```cpp
// Called thousands of times during B&B
SolverStatus solve_from_basis(const Problem&, const std::vector<Index>& basis, Solution&);

// Farkas ray for conflict analysis (when LP is infeasible)
const std::vector<Real>& get_farkas_ray() const;

// For Gomory cuts — expose tableau row (BTRAN)
virtual void btran_row(Index basic_pos, std::vector<Real>& dense_row) const;
```

### From Person 4 (Problem construction)
```cpp
// Person 4 parses MPS/LP files and hands you a Problem object
Problem problem;
mps_reader.read("instance.mps", problem);
milp_solver.solve(problem, solution);
```

---

## Numerical Notes

- **EPS_INTEGER = 1e-6** — a variable is "integer" if `|frac| < 1e-6`
- **EPS_FEASIBILITY = 1e-8** — constraint violation tolerance
- **Pseudo-cost initialisation = 1.0** — prevents division by zero; converges quickly
- **Probing is expensive** — enable only for problems with binary variables, limit to
  first few rounds of presolve
- **Gomory cuts can be numerically unstable** — skip rows with `f0 < 0.05 or f0 > 0.95`
- **Feasibility pump cycling** — detect when `x_tilde == prev_tilde` and perturb;
  reduce α geometrically to push toward feasibility

---

## Build Commands

```bash
# Syntax check all sources (no cmake needed):
cd /mnt/d_drive/programs/suplex-person3
for f in src/**/*.cpp; do
  g++ -std=c++20 -Wall -Wextra -Isrc -fsyntax-only "$f"
done

# When cmake is available:
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```
