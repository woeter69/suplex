# PERSON 3: Presolve Engine & MILP Framework

## Your Mission
You build two critical subsystems: (1) the presolve engine that simplifies problems before solving, and (2) the mixed-integer programming (MILP) engine that handles integer variables via branch-and-bound/branch-and-cut. Presolve is used for ALL problem types (LP, MILP, QP), while B&B is MILP-specific.

## Deliverables Checklist

- [ ] `src/milp/presolve.h/cpp` — Presolve engine with 10+ reduction rules (7 days)
- [ ] `src/milp/postsolve.h/cpp` — Postsolve (undo) mechanism (3 days)
- [ ] `src/milp/node.h/cpp` — B&B tree node data structure (2 days)
- [ ] `src/milp/branch_bound.h/cpp` — Branch-and-bound tree manager (5 days)
- [ ] `src/milp/branching.h/cpp` — Variable selection strategies (4 days)
- [ ] `src/milp/cuts.h/cpp` — Cutting plane generators (7 days)
- [ ] `src/milp/heuristics.h/cpp` — Primal heuristics (5 days)
- [ ] `src/milp/conflict.h/cpp` — Conflict analysis (3 days)
- [ ] Unit tests (5 days)
- [ ] MIPLIB benchmark validation (4 days)

Total: ~45 working days

## PART A: PRESOLVE ENGINE

Presolve is one of the most impactful components. A good presolve can reduce problem size by 50-90% and dramatically improve solver performance. It runs BEFORE the LP/MILP/QP solver.

## 1. Presolve Architecture

### 1.1 Undo Stack
Every presolve reduction must be REVERSIBLE. The postsolve mechanism undoes reductions in reverse order to recover the solution to the original problem.

Design:
```cpp
enum class PresolveRuleType {
    EMPTY_ROW,
    SINGLETON_ROW,
    SINGLETON_COL,
    FORCING_ROW,
    DOMINATED_COL,
    DUPLICATE_ROW,
    DUPLICATE_COL,
    FREE_COL_SUBSTITUTION,
    IMPLIED_BOUND,
    FIXED_VARIABLE,
    PROBING,
    COEFFICIENT_TIGHTENING
};

struct PresolveRecord {
    PresolveRuleType rule;
    // Rule-specific data to undo the transformation
    Index row_or_col;
    std::vector<Index> affected_indices;
    std::vector<Real> saved_values;
    Real saved_bound;
};

class PresolveStack {
    std::vector<PresolveRecord> records;
public:
    void push(PresolveRecord record);
    PresolveRecord pop();
    bool empty() const;
};
```

### 1.2 Presolve Loop
Run reductions in a loop until no more reductions apply:
```
PRESOLVE(Problem P):
  changed = true
  while changed:
    changed = false
    changed |= remove_empty_rows(P)
    changed |= remove_fixed_variables(P)
    changed |= singleton_row_reduction(P)
    changed |= singleton_col_reduction(P) 
    changed |= forcing_row_reduction(P)
    changed |= dominated_column_reduction(P)
    changed |= duplicate_row_detection(P)
    changed |= duplicate_col_detection(P)
    changed |= implied_bound_tightening(P)
    changed |= coefficient_tightening(P)  // MILP only
    changed |= probing(P)                  // MILP only
  return reduced problem + undo stack
```

## 2. Presolve Reductions (Detail Each)

### 2.1 Empty Row
- If row i has no nonzero coefficients: check if 0 ∈ [row_lower[i], row_upper[i]]. If yes, remove row. If no, problem is INFEASIBLE.

### 2.2 Fixed Variable
- If col_lower[j] == col_upper[j] (within tolerance): substitute x_j = col_lower[j], update b = b - A[:,j] * x_j, remove column j.

### 2.3 Singleton Row
Row with exactly one nonzero a_ij: implies a direct bound on x_j.

```
SINGLETON_ROW(P, row i, nonzero a_ij at column j):
  // The constraint is: row_lower[i] <= a_ij * x_j <= row_upper[i]
  if a_ij > 0:
    implied_lower = row_lower[i] / a_ij   // -INF/positive = -INF (ok)
    implied_upper = row_upper[i] / a_ij
  else:  // a_ij < 0 — division flips inequality
    implied_lower = row_upper[i] / a_ij   // positive / negative = negative
    implied_upper = row_lower[i] / a_ij

  // Tighten variable bounds
  new_lower = max(col_lower[j], implied_lower)
  new_upper = min(col_upper[j], implied_upper)

  if new_lower > new_upper + EPS_FEASIBILITY:
    return INFEASIBLE

  // Record for postsolve
  stack.push({SINGLETON_ROW, row=i, col=j, a_ij,
              old_lower=col_lower[j], old_upper=col_upper[j],
              row_lower[i], row_upper[i]})

  col_lower[j] = new_lower
  col_upper[j] = new_upper
  remove row i from problem
  return CHANGED
```

### 2.4 Singleton Column
Column j with exactly one nonzero a_ij — x_j can be determined from the constraint and objective:

```
SINGLETON_COL(P, col j, nonzero a_ij at row i):
  // x_j appears only in row i: row_lower[i] <= ... + a_ij*x_j + ... <= row_upper[i]
  // And in objective: c_j * x_j

  // Determine optimal value of x_j based on objective direction
  if c_j > EPS_ZERO:       // minimizing: want x_j as small as possible
    x_j_opt = col_lower[j]
  elif c_j < -EPS_ZERO:    // minimizing: want x_j as large as possible
    x_j_opt = col_upper[j]
  else:                     // c_j ≈ 0: x_j free to take any value
    // Choose value that gives most slack in constraint
    // Compute feasible range from row bounds
    // x_j ∈ [(row_lower[i] - activity_without_j) / a_ij,
    //         (row_upper[i] - activity_without_j) / a_ij]
    // Pick midpoint or lower bound
    x_j_opt = col_lower[j]

  // Substitute x_j = x_j_opt into the constraint
  // Update row bounds: row_lower[i] -= a_ij * x_j_opt
  //                    row_upper[i] -= a_ij * x_j_opt
  // Update objective constant: obj_offset += c_j * x_j_opt

  stack.push({SINGLETON_COL, row=i, col=j, a_ij, c_j, x_j_opt,
              saved row bounds, saved col bounds})

  row_lower[i] -= a_ij * x_j_opt
  row_upper[i] -= a_ij * x_j_opt
  remove column j and remove a_ij from row i
  return CHANGED
```

### 2.5 Forcing Row
A row i where the constraint forces ALL variables in that row to their bounds:

```
FORCING_ROW(P, row i):
  // Compute min and max activity for this row
  min_act = 0.0;  max_act = 0.0
  for each (j, a_ij) in row i:
    if a_ij > 0:
      min_act += a_ij * col_lower[j]
      max_act += a_ij * col_upper[j]
    else:
      min_act += a_ij * col_upper[j]
      max_act += a_ij * col_lower[j]

  // Check if row forces variables to achieve max activity
  if max_act <= row_upper[i] + EPS_FEASIBILITY:
    // ALL variables must be at the bound that achieves their max contribution
    for each (j, a_ij) in row i:
      if a_ij > 0:
        fix_value = col_upper[j]  // max contribution when at upper
      else:
        fix_value = col_lower[j]  // max contribution when at lower (neg coeff)
      stack.push({FORCING_ROW, col=j, fix_value, old bounds})
      col_lower[j] = fix_value
      col_upper[j] = fix_value
    remove row i
    return CHANGED

  // Check if row forces variables to achieve min activity
  if min_act >= row_lower[i] - EPS_FEASIBILITY:
    for each (j, a_ij) in row i:
      if a_ij > 0:
        fix_value = col_lower[j]
      else:
        fix_value = col_upper[j]
      stack.push({FORCING_ROW, col=j, fix_value, old bounds})
      col_lower[j] = fix_value
      col_upper[j] = fix_value
    remove row i
    return CHANGED

  return NO_CHANGE
```

### 2.6 Dominated Column
- Column j is dominated by column k if: same constraint coefficients (up to scaling) and c_j/a_ij >= c_k/a_ik for all rows i.
- In this case, x_j can be fixed at its lower bound.

### 2.7 Implied Bound Tightening
- For each row i, compute min/max activity:
  - max_act_without_j = max_activity - contribution_of_j
  - This implies an upper bound on a_ij * x_j: row_upper[i] - max_act_without_j
  - → implied bound on x_j
- If implied bound is tighter than current bound: tighten it.
- Iterate until no more tightening.

### 2.8 Coefficient Tightening (MILP only)
- For integer variables, round coefficients and bounds to exploit integrality
- Example: 3.7 x ≤ 10 with x integer → 3 x ≤ 8 (since x ≤ 2.7 → x ≤ 2)

### 2.9 Probing (MILP only)
For binary variables, temporarily fix and propagate to discover implied bounds:

```
PROBING(P):
  changed = false
  for each binary variable x_j (col_lower[j]=0, col_upper[j]=1):
    // Probe x_j = 0
    saved_bounds_0 = copy(col_lower, col_upper)
    col_upper[j] = 0
    status_0 = propagate_bounds(P)  // run implied bound tightening
    bounds_at_0 = copy(col_lower, col_upper)
    restore(col_lower, col_upper, saved_bounds_0)

    // Probe x_j = 1
    saved_bounds_1 = copy(col_lower, col_upper)
    col_lower[j] = 1
    status_1 = propagate_bounds(P)
    bounds_at_1 = copy(col_lower, col_upper)
    restore(col_lower, col_upper, saved_bounds_1)

    // Case 1: One fixing is infeasible → fix to other value
    if status_0 == INFEASIBLE and status_1 == INFEASIBLE:
      return INFEASIBLE  // problem infeasible
    if status_0 == INFEASIBLE:
      col_lower[j] = 1; col_upper[j] = 1  // must be 1
      stack.push({PROBING, col=j, fixed_value=1})
      changed = true; continue
    if status_1 == INFEASIBLE:
      col_lower[j] = 0; col_upper[j] = 0  // must be 0
      stack.push({PROBING, col=j, fixed_value=0})
      changed = true; continue

    // Case 2: Both fixings agree on a bound for another variable
    for each variable k ≠ j:
      agreed_lower = max(bounds_at_0.lower[k], bounds_at_1.lower[k])
      agreed_upper = min(bounds_at_0.upper[k], bounds_at_1.upper[k])
      if agreed_lower > col_lower[k] + EPS_FEASIBILITY:
        col_lower[k] = agreed_lower
        changed = true
      if agreed_upper < col_upper[k] - EPS_FEASIBILITY:
        col_upper[k] = agreed_upper
        changed = true

  return changed
```

## 3. Postsolve

Recover the original-space solution from the reduced-space solution. Process the undo stack in **reverse order**:

```
POSTSOLVE(Solution reduced_sol, PresolveStack stack):
  x = reduced_sol.primal_values   // indexed by reduced-space columns
  y = reduced_sol.dual_values     // indexed by reduced-space rows

  while not stack.empty():
    record = stack.pop()
    switch record.rule:

      case FIXED_VARIABLE:
        // x_j was fixed to value v and removed
        // Restore: x_j = v, dual contribution from all rows containing j
        insert x[j] = record.fix_value into solution
        // Compute reduced cost for optimality certificate
        rc[j] = c[j] - sum(y[i] * a_ij for each row i containing j)

      case SINGLETON_ROW:
        // Row i had one nonzero a_ij. Row was removed, bounds tightened on x_j.
        // Restore: compute row dual y[i] from optimality
        // If x_j is at its (tightened) lower bound: y[i] = (c_j - rc_j) / a_ij
        // Simpler: y[i] is determined by the bound that was active
        if x[j] == record.new_lower:
          y[i] = rc_j_from_original / a_ij  // adjust sign based on bound type
        else:
          y[i] = 0.0  // constraint not tight → dual is zero

      case SINGLETON_COL:
        // Column j had one nonzero a_ij in row i. Both removed.
        // x_j was set to optimal value. Now compute y[i]:
        // Optimality: c_j - y[i] * a_ij = 0 (if x_j is between bounds)
        // So: y[i] = c_j / a_ij
        insert x[j] = record.x_j_opt
        if col_lower[j] < x[j] < col_upper[j]:
          y[i] = record.c_j / record.a_ij
        else:
          // x_j at bound — y[i] free, compute from constraint
          y[i] = 0.0  // or derive from dual feasibility

      case FORCING_ROW:
        // Variables were fixed at bounds, row removed
        // Restore: x[j] = fix_value (already set during fixed_variable postsolve)
        // Compute y[i] for the removed row:
        // Activity = sum(a_ij * x_j) for j in row
        // If activity == row_upper[i]: y[i] >= 0 (for <= constraint)
        // If activity == row_lower[i]: y[i] <= 0 (for >= constraint)
        // Choose y[i] to satisfy dual feasibility of all fixed variables

      case PROBING:
        // Variable was fixed by probing — same as FIXED_VARIABLE
        insert x[j] = record.fixed_value

  return Solution(x, y, reduced_costs)
```

## PART B: MILP ENGINE

## 4. Branch-and-Bound Tree

### 4.1 Node Data Structure
```cpp
struct BBNode {
    int64_t id;
    int64_t parent_id;          // -1 for root
    int depth;
    
    // Branching info
    Index branch_var;            // which variable was branched on
    Real branch_value;           // the fractional value
    BranchDirection direction;   // UP (>= ceil) or DOWN (<= floor)
    
    // Bound changes from parent
    std::vector<BoundChange> bound_changes;
    
    // LP relaxation result
    Real lp_bound;               // LP optimal value at this node
    SolverStatus lp_status;
    
    // Basis for warm-starting LP solve at child nodes
    std::vector<Index> basis;
    
    // State
    NodeStatus status;           // OPEN, SOLVED, PRUNED, BRANCHED, INFEASIBLE
};

enum class BranchDirection { DOWN, UP };
enum class NodeStatus { OPEN, SOLVED, PRUNED, BRANCHED, INFEASIBLE };

struct BoundChange {
    Index var;
    Real new_lower;
    Real new_upper;
};
```

### 4.2 Node Pool
- Priority queue of open nodes
- Support multiple selection strategies (see section 5)
- Memory-efficient: only store CHANGES from parent, not full problem copy

### 4.3 Tree Statistics
Track and log:
- Total nodes created/explored/pruned
- Current best incumbent (best integer-feasible solution found so far)
- Current best bound (best LP relaxation among open nodes)
- Gap = |incumbent - best_bound| / |incumbent|
- Nodes per second
- Tree depth distribution

## 5. Node Selection Strategies

### 5.1 Best-First Search
- Select node with best (lowest for min) LP bound
- Pro: Tightest bounds, proves optimality efficiently
- Con: Large memory (many open nodes)

### 5.2 Depth-First Search
- Select deepest node (LIFO)
- Pro: Small memory, warm start benefit (parent basis similar)
- Con: May explore many nodes before proving optimality

### 5.3 Best-Estimate Search
- Estimate the integer objective at each node using pseudo-costs
- Select node with best estimate

### 5.4 Hybrid Strategy (RECOMMENDED DEFAULT)

```
HYBRID_NODE_SELECT(node_pool, last_node, incumbent):
  // PLUNGING: after branching, immediately dive depth-first
  if last_node != NULL and last_node.status == BRANCHED:
    // Pick the child with the better (lower) LP bound
    child = best_child_of(last_node)
    if child.lp_bound < incumbent - EPS_OPTIMALITY:
      return child  // continue plunging

  // BACKTRACK: plunge ended (infeasible, pruned, or integer-feasible)
  // Switch to best-first: pick the open node with best LP bound
  if node_pool.empty():
    return NULL  // no more nodes

  return node_pool.pop_best_bound()  // priority queue by LP bound
```

## 6. Branching Variable Selection

### 6.1 Most Fractional
- Select the integer variable closest to 0.5 (most fractional)
- Simple but poor performance

### 6.2 Pseudo-Cost Branching
Maintain pseudo-costs: average objective change per unit change in variable value when branching.

```
PSEUDOCOST_SCORE(j, frac_val, pseudocost_down, pseudocost_up):
  f = frac_val - floor(frac_val)      // fractional part
  down_est = pseudocost_down[j] * f    // estimated obj change branching down
  up_est   = pseudocost_up[j] * (1-f)  // estimated obj change branching up

  // Product scoring (used in SCIP — better than sum):
  score = (1 - mu) * min(down_est, up_est) + mu * max(down_est, up_est)
  // where mu = 1/6 (SCIP default)
  return score

UPDATE_PSEUDOCOSTS(j, direction, frac_val, parent_bound, child_bound):
  delta_obj = child_bound - parent_bound   // LP bound change
  if direction == DOWN:
    delta_var = frac_val - floor(frac_val)
    pseudocost_down[j] = update_average(pseudocost_down[j], delta_obj / delta_var)
    pseudocost_down_count[j] += 1
  else:
    delta_var = ceil(frac_val) - frac_val
    pseudocost_up[j] = update_average(pseudocost_up[j], delta_obj / delta_var)
    pseudocost_up_count[j] += 1
```

### 6.3 Strong Branching
- For the top-K most fractional variables (K ≈ 10-20):
  - Actually solve the LP relaxation for both child nodes (x_j <= floor, x_j >= ceil)
  - Score based on actual LP bound change
- Very expensive but best quality decisions
- Use only at root node or first few levels

### 6.4 Reliability Branching (RECOMMENDED DEFAULT)

```
RELIABILITY_BRANCH(lp_solution, fractional_vars, pseudocosts, rel_threshold=8):
  candidates = []
  strong_branch_needed = []

  for each fractional variable j in fractional_vars:
    if pseudocost_down_count[j] >= rel_threshold
       AND pseudocost_up_count[j] >= rel_threshold:
      // Pseudo-costs are reliable — use them
      score = PSEUDOCOST_SCORE(j, lp_solution[j], ...)
      candidates.append((j, score))
    else:
      // Not enough observations — need strong branching
      strong_branch_needed.append(j)

  // Strong branch the unreliable candidates (limit to top-20 by fractionality)
  sort strong_branch_needed by |frac - 0.5| ascending (most fractional first)
  for j in strong_branch_needed[:20]:
    // Solve LP with x_j <= floor(val)
    down_bound = dual_simplex_solve_with_bound_change(j, upper=floor(val))
    // Solve LP with x_j >= ceil(val)
    up_bound = dual_simplex_solve_with_bound_change(j, lower=ceil(val))
    UPDATE_PSEUDOCOSTS(j, DOWN, val, current_bound, down_bound)
    UPDATE_PSEUDOCOSTS(j, UP, val, current_bound, up_bound)
    score = PSEUDOCOST_SCORE(j, val, ...)
    candidates.append((j, score))

  // Select variable with highest score
  best = max(candidates, key=lambda x: x.score)
  return best.j
```

## 7. Cutting Planes

Cutting planes tighten the LP relaxation, reducing the gap and the number of B&B nodes.

### 7.1 Cut Management
```cpp
struct Cut {
    std::vector<Index> indices;   // variable indices
    std::vector<Real> coefficients;
    Real rhs;
    CutType type;                  // GOMORY, MIR, CLIQUE, COVER
    Real efficacy;                 // how much it cuts off the current LP solution
    int age;                       // rounds since last active in basis
};
```

### 7.2 Gomory Mixed-Integer Cuts

```
GENERATE_GOMORY_CUTS(basis, lp_solution, A, max_cuts=50):
  cuts = []
  for each basic variable i where var_type[basis[i]] == INTEGER:
    val = lp_solution[basis[i]]
    f_0 = val - floor(val)

    // Skip near-integer rows (numerically unsafe)
    if f_0 < 0.05 or f_0 > 0.95:
      continue

    // Extract tableau row: x_i + Σ_j a_bar_ij x_j = b_bar_i
    // a_bar_ij = (B^{-1} A)_{i,j} for non-basic j
    // Use BTRAN: compute row i of B^{-1}, then multiply by each A column
    e_i = unit_vector(i)
    beta = BTRAN(e_i)   // row i of B^{-1}

    cut_indices = []
    cut_coeffs = []

    for each non-basic variable j:
      a_bar = DOT_COL(A, j, beta)    // tableau coefficient

      // Skip tiny coefficients
      if |a_bar| < EPS_ZERO: continue

      f_j = a_bar - floor(a_bar)   // fractional part (always in [0,1))

      // Gomory coefficient
      if f_j <= f_0:
        gomory_coeff = f_j / f_0
      else:
        gomory_coeff = (1 - f_j) / (1 - f_0)

      // Strengthening: if x_j is integer, use stronger coefficient
      if var_type[j] == INTEGER:
        if f_j <= f_0:
          gomory_coeff = f_j / f_0
        else:
          gomory_coeff = (1 - f_j) / (1 - f_0) * f_0 / (1 - f_0)
          // Actually: for integer vars, the MIR strengthening applies

      if gomory_coeff > EPS_ZERO:
        cut_indices.append(j)
        cut_coeffs.append(gomory_coeff)

    // Gomory cut: Σ gomory_coeff_j * x_j >= 1
    // Compute efficacy (distance from LP solution to cut)
    lhs_at_lp = sum(cut_coeffs[k] * lp_solution[cut_indices[k]] for k)
    efficacy = (1.0 - lhs_at_lp) / norm(cut_coeffs)

    if efficacy > 0.01:  // only add effective cuts
      cuts.append(Cut(cut_indices, cut_coeffs, rhs=1.0, GOMORY, efficacy))

    if len(cuts) >= max_cuts: break

  return cuts
```

### 7.3 Mixed-Integer Rounding (MIR) Cuts

```
GENERATE_MIR_CUT(row_coeffs, row_rhs, var_types):
  // Given constraint: Σ a_j x_j <= b
  // Partition into integer (I) and continuous (C) variables
  // MIR: floor(b) + (Σ_{j∈I: f_j > f_0} (f_j - f_0)/(1-f_0) * x_j
  //       + Σ_{j∈C: a_j > 0} a_j/(1-f_0) * x_j) >= Σ_{j∈I} floor(a_j) * x_j

  f_0 = row_rhs - floor(row_rhs)
  if f_0 < 0.05 or f_0 > 0.95: return NULL  // skip

  cut_lhs = []
  for each variable j:
    if var_types[j] == INTEGER:
      f_j = row_coeffs[j] - floor(row_coeffs[j])
      if f_j <= f_0:
        mir_coeff = floor(row_coeffs[j])
      else:
        mir_coeff = floor(row_coeffs[j]) + (f_j - f_0) / (1 - f_0)
    else:  // continuous
      if row_coeffs[j] >= 0:
        mir_coeff = row_coeffs[j] / (1 - f_0)
      else:
        mir_coeff = row_coeffs[j]   // negative coefficients unchanged
    cut_lhs.append((j, mir_coeff))

  cut_rhs = floor(row_rhs)
  return Cut(cut_lhs, cut_rhs, type=MIR)
```

### 7.4 Clique Cuts
- For binary variables: if x_i + x_j <= 1 is implied (from constraint or bound), this is a clique
- Build conflict graph: edge between binary vars that can't both be 1
- Find maximal cliques → clique cuts: Σ_{j in clique} x_j <= 1

### 7.5 Knapsack Cover Cuts

```
GENERATE_COVER_CUT(knapsack_row, b):
  // knapsack_row: a_j x_j for binary variables, Σ a_j x_j <= b
  // Find minimal cover C: Σ_{j∈C} a_j > b

  // Greedy cover: sort by coefficient descending, add until sum > b
  sorted_vars = sort(knapsack_row, by=a_j, descending)
  cover = []
  cover_sum = 0
  for (j, a_j) in sorted_vars:
    cover.append(j)
    cover_sum += a_j
    if cover_sum > b:
      break

  if cover_sum <= b: return NULL  // no cover exists

  // Basic cover inequality: Σ_{j∈C} x_j <= |C| - 1
  // Lifted cover: extend to non-cover variables
  cut_coeffs = {}
  for j in cover:
    cut_coeffs[j] = 1.0
  cut_rhs = len(cover) - 1

  // Sequential lifting: for each non-cover variable k
  for (k, a_k) in knapsack_row:
    if k in cover: continue
    // Find largest alpha such that cover + alpha*x_k is valid
    // Solve: max alpha s.t. Σ_{j∈C} x_j + alpha*x_k <= |C|-1
    //        whenever Σ_{j∈C} a_j x_j + a_k * x_k <= b
    // Use the lifting procedure (solve small knapsack)
    alpha = compute_lifting_coefficient(cover, a_k, b)
    if alpha > EPS_ZERO:
      cut_coeffs[k] = alpha

  return Cut(cut_coeffs, cut_rhs, type=COVER)
```

### 7.6 Cut Selection & Management
- Generate many candidate cuts, select the best:
  - Efficacy: distance from current LP solution to the cut hyperplane
  - Orthogonality: prefer cuts that are not parallel to existing cuts
  - Sparsity: prefer sparser cuts (fewer nonzeros)
- Cut pool: store generated cuts, age them, remove inactive cuts
- Limit: add at most 100 cuts per round, remove cuts with age > 10

## 8. Primal Heuristics

Heuristics find good integer-feasible solutions quickly, providing tight incumbent bounds that help prune the B&B tree.

### 8.1 Simple Rounding
- Round each fractional integer variable to nearest integer
- Check feasibility of rounded solution
- If feasible: new incumbent!

### 8.2 Diving Heuristics

```
FRACTIONAL_DIVING(P, lp_solution, max_dives=100):
  current_bounds = copy(col_lower, col_upper)

  for dive = 1 to max_dives:
    // Find most fractional integer variable
    best_j = -1; best_frac = 0
    for each integer var j:
      f = lp_solution[j] - floor(lp_solution[j])
      frac_distance = min(f, 1 - f)
      if frac_distance > EPS_INTEGER and frac_distance > best_frac:
        best_frac = frac_distance
        best_j = j

    if best_j == -1:
      // All integer vars are integral → feasible solution found!
      return lp_solution  // new incumbent

    // Fix to nearest integer
    val = lp_solution[best_j]
    if val - floor(val) <= 0.5:
      col_upper[best_j] = floor(val)
    else:
      col_lower[best_j] = ceil(val)

    // Re-solve LP with new bound
    status = dual_simplex.solve_from_basis(P_modified, current_basis)

    if status == INFEASIBLE:
      break  // dive failed
    lp_solution = dual_simplex.solution()
    current_basis = dual_simplex.basis()

  // Restore original bounds
  restore(col_lower, col_upper, current_bounds)
  return NULL  // no feasible solution found
```

### 8.3 Feasibility Pump

```
FEASIBILITY_PUMP(P, lp_solution, max_iter=100):
  x_star = lp_solution   // LP relaxation solution
  alpha = 0.85

  for iter = 1 to max_iter:
    // Round to nearest integer
    x_tilde = copy(x_star)
    for each integer variable j:
      x_tilde[j] = round(x_star[j])

    // Check if rounded solution is feasible
    if is_feasible(P, x_tilde):
      return x_tilde  // found feasible solution!

    // Check for cycling: if x_tilde == previous x_tilde
    if iter > 1 and x_tilde == prev_tilde:
      // Perturb: randomly flip some roundings
      for each integer var j with probability 0.1:
        if x_tilde[j] == floor(x_star[j]):
          x_tilde[j] = ceil(x_star[j])
        else:
          x_tilde[j] = floor(x_star[j])

    prev_tilde = x_tilde

    // Solve LP: minimize distance to rounded solution
    // min Σ_{j ∈ integers} |x_j - x_tilde_j|
    // subject to: Ax ∈ [row_lower, row_upper], x ∈ [col_lower, col_upper]
    // Linearize absolute value with auxiliary variables:
    // min Σ d_j where d_j >= x_j - x_tilde_j, d_j >= x_tilde_j - x_j
    // Or use a weighted combination with original objective:
    // min alpha * Σ|x_j - x_tilde_j| + (1-alpha) * c^T x
    x_star = solve_distance_lp(P, x_tilde, alpha)

    // Gradually reduce alpha to push toward feasibility
    alpha *= 0.95

  return NULL  // failed to find feasible solution
```

### 8.4 RINS (Relaxation Induced Neighborhood Search)
- Fix variables where LP relaxation and incumbent agree (both at same integer value)
- Solve restricted MIP on remaining variables (with node limit)

## 9. Conflict Analysis & Learning

When a node is pruned due to infeasibility, analyze WHY:

```
CONFLICT_ANALYSIS(infeasible_node, branching_path):
  // Step 1: Get infeasibility proof from dual simplex
  // The Farkas ray y satisfies: y^T A >= 0, y^T b < 0 (proving infeasibility)
  farkas_ray = dual_simplex.get_farkas_ray()

  // Step 2: Identify which bound changes are responsible
  // The Farkas ray uses specific constraints (rows with y[i] != 0)
  responsible_rows = {i : |farkas_ray[i]| > EPS_ZERO}

  // Step 3: Trace back through branching decisions
  conflict_set = []
  for each branching decision (var_j, bound_change) in branching_path:
    // Check if this bound change affects any responsible row
    if variable j appears in any responsible row:
      conflict_set.append((var_j, bound_change))

  // Step 4: Minimize conflict set (remove redundant decisions)
  // Use the resolvent: if a bound change was implied by others, remove it
  minimal_conflict = minimize_conflict(conflict_set, farkas_ray)

  // Step 5: Generate conflict constraint (no-good cut)
  // For binary branching: conflict is a clause
  // "Not all these bound changes can be active"
  // Σ (x_j if branched down) + Σ (1-x_j if branched up) >= 1
  cut = generate_nogood_cut(minimal_conflict)

  return cut  // add to global cut pool
```

## 10. Overall MILP Solve Algorithm

```
MILP_SOLVE(Problem P):
  1. Presolve P → P_reduced
  2. Solve root LP relaxation of P_reduced
  3. If root LP is infeasible → INFEASIBLE
  4. If root LP solution is integer-feasible → OPTIMAL (lucky!)
  5. Run cut separation at root node (multiple rounds)
  6. Run root heuristics (rounding, feasibility pump)
  7. Initialize B&B tree with root node
  
  WHILE open nodes exist AND gap > tolerance AND within limits:
    8. Select next node (hybrid strategy)
    9. Apply bound changes from node's branch
    10. Solve LP relaxation (warm start from parent basis)
    11. If LP infeasible → prune node, conflict analysis
    12. If LP bound >= incumbent → prune node (bound)
    13. If LP solution is integer-feasible → update incumbent
    14. Else:
        a. Separate cuts (every N nodes)
        b. Run heuristics (every M nodes)
        c. Select branching variable
        d. Create two child nodes
        e. Add children to node pool
    15. Update best bound (min LP bound among open nodes)
    16. Log progress
  
  9. Postsolve incumbent solution
  10. Return solution
```

## 11. Testing Strategy

### Presolve Tests
- Each reduction rule: test on hand-crafted 5x5 problem where that rule applies
- Test postsolve: verify original-space solution is feasible and optimal
- Test that presolve detects infeasibility when appropriate
- Test reduction count: on Netlib `afiro`, verify specific reductions happen

### MILP Tests
- Trivial: 2-variable IP with known solution
- Small: knapsack problem (10 items) — known optimal
- Medium: set covering problem (50 rows, 200 cols)
- Prove infeasibility on infeasible MIP
- Verify cutting planes cut off LP solution but not IP optimum

### MIPLIB Benchmark

| Problem | Type | Rows | Cols | Difficulty |
|---------|------|------|------|-----------|
| `air04` | Set partitioning | 823 | 8904 | Easy |
| `bell5` | Network | 92 | 104 | Easy |
| `blend2` | Blending | 274 | 353 | Easy |
| `dcmulti` | Network | 291 | 548 | Easy |
| `egout` | Network | 99 | 141 | Easy |
| `enigma` | Puzzle | 22 | 100 | Easy |
| `flugpl` | Planning | 19 | 18 | Easy |
| `gen` | General | 781 | 870 | Easy |
| `p0033` | Covering | 17 | 33 | Easy |
| `p0201` | Covering | 134 | 201 | Moderate |
| `bell3a` | Network | 124 | 133 | Moderate |
| `mod008` | General | 7 | 319 | Moderate |
| `mod010` | General | 147 | 2655 | Moderate |

**Target:** solve 35+ of 50 easy MIPLIB instances to optimality within 1 hour each

## 12. Dependencies on Other Persons

| Who | What You Need | What You Provide | Priority |
|-----|--------------|-----------------|----------|
| **Person 1** | `DualSimplexSolver::solve_from_basis()` for every B&B node LP | — | CRITICAL |
| **Person 2** | `IPMSolver` optionally for root node relaxation | — | Nice-to-have |
| **Person 4** | `Problem` object (parsed from MPS/LP) | `Presolve`, `MILPSolver` | Week 4 |

## 13. Key References
- Achterberg, T. "Constraint Integer Programming" (PhD thesis, 2007) — SCIP's foundation
- Achterberg, T. et al. "SCIP: Solving Constraint Integer Programs" (Math. Prog. Computation, 2009)
- Wolsey, L.A. "Integer Programming" (Wiley, 1998) — chapters on cutting planes
- Savelsbergh, M. "Preprocessing and Probing Techniques for MIP" (ORSA, 1994)
- Berthold, T. "Primal Heuristics for Mixed Integer Programs" (PhD thesis, 2006)
- Witzig, J. et al. "Conflict-Driven Heuristics for MIP" (2019)

## 14. Definition of Done
- [ ] All presolve rules implemented and tested with postsolve
- [ ] B&B framework solves trivial MIPs correctly
- [ ] At least 3 cutting plane types working (Gomory, MIR, knapsack cover)
- [ ] At least 2 heuristics working (rounding, feasibility pump)
- [ ] Solves 35+ of 50 easy MIPLIB problems within 1 hour each
- [ ] Gap reporting is accurate
- [ ] No memory leaks
- [ ] Doxygen documentation complete
