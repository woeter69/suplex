# PERSON 1: Sparse Linear Algebra & Simplex Engine

## Your Mission
You build the numerical bedrock of Suplex. Every other solver component depends on your sparse matrix infrastructure and LP solver. Your code must be numerically bulletproof.

**Ownership:** `src/core/sparse_matrix.h/cpp`, `src/simplex/` (entire directory)

---

## Deliverables Checklist

| # | File | Description | Effort |
|---|------|-------------|--------|
| 1 | `src/core/sparse_matrix.h/cpp` | CSC, CSR, Triplet formats | 3 days |
| 2 | `src/simplex/lu_factor.h/cpp` | Sparse LU factorization with Markowitz ordering | 5 days |
| 3 | `src/simplex/lu_update.h/cpp` | Forrest-Tomlin LU update | 4 days |
| 4 | `src/simplex/basis.h/cpp` | Basis management, crash procedures | 2 days |
| 5 | `src/simplex/pricing.h/cpp` | Dantzig, steepest edge, Devex pricing | 3 days |
| 6 | `src/simplex/primal_simplex.h/cpp` | Full revised primal simplex | 5 days |
| 7 | `src/simplex/dual_simplex.h/cpp` | Full revised dual simplex | 5 days |
| 8 | Unit tests for all above | Ongoing | 5 days |
| 9 | Netlib LP benchmark validation | End-to-end | 3 days |

**Total: ~35 working days**

---

## 1. Sparse Matrix Implementation

### 1.1 CSC (Compressed Sparse Column) Format

The CSC format stores a sparse m×n matrix using three contiguous arrays:
- `col_start[num_cols + 1]` — column j's nonzeros span indices `col_start[j]` to `col_start[j+1] - 1`
- `row_index[nnz]` — row indices of each nonzero, sorted within each column
- `values[nnz]` — corresponding nonzero values

**Why CSC:** The simplex method is column-oriented — pricing computes `y^T a_j` (column dot), FTRAN solves `B^{-1} a_q` (column operation). CSC makes column access O(nnz_in_col).

#### Construction from Triplets

```
FROM_TRIPLETS(nrows, ncols, trip_rows[], trip_cols[], trip_vals[], num_entries):
  // Step 1: Count entries per column
  col_count = array[ncols], initialized to 0
  for k = 0 to num_entries - 1:
    col_count[trip_cols[k]] += 1

  // Step 2: Build col_start via prefix sum
  col_start = array[ncols + 1]
  col_start[0] = 0
  for j = 0 to ncols - 1:
    col_start[j + 1] = col_start[j] + col_count[j]
  nnz = col_start[ncols]

  // Step 3: Populate row_index and values using col_start as write cursors
  cursor = copy of col_start[0..ncols-1]
  row_index = array[nnz]
  values = array[nnz]
  for k = 0 to num_entries - 1:
    j = trip_cols[k]
    pos = cursor[j]
    row_index[pos] = trip_rows[k]
    values[pos] = trip_vals[k]
    cursor[j] += 1

  // Step 4: Sort each column by row index and sum duplicates
  for j = 0 to ncols - 1:
    start = col_start[j]
    end = col_start[j + 1]
    sort (row_index[start..end-1], values[start..end-1]) by row_index

    // Sum duplicates in-place
    write_pos = start
    for k = start to end - 1:
      if write_pos > start AND row_index[k] == row_index[write_pos - 1]:
        values[write_pos - 1] += values[k]   // sum duplicate
      else:
        row_index[write_pos] = row_index[k]
        values[write_pos] = values[k]
        write_pos += 1
    // Adjust col_start if duplicates were found (compact storage)
    // Shift subsequent data if nnz decreased

  return CSC(nrows, ncols, col_start, row_index, values)
```

#### Key Operations

```
MULTIPLY_VEC(A_csc, x, y):        // y = A * x
  fill y[0..num_rows-1] with 0.0
  for j = 0 to num_cols - 1:
    if x[j] == 0.0: continue       // skip zero entries
    for k = col_start[j] to col_start[j+1] - 1:
      y[row_index[k]] += values[k] * x[j]

MULTIPLY_TRANSPOSE_VEC(A_csc, x, y):  // y = A^T * x
  for j = 0 to num_cols - 1:
    sum = 0.0
    for k = col_start[j] to col_start[j+1] - 1:
      sum += values[k] * x[row_index[k]]
    y[j] = sum

DOT_COL(A_csc, j, x):             // returns a_j^T * x
  sum = 0.0
  for k = col_start[j] to col_start[j+1] - 1:
    sum += values[k] * x[row_index[k]]
  return sum

MULTIPLY_COL(A_csc, j, scalar, y):  // y += scalar * A[:,j]
  for k = col_start[j] to col_start[j+1] - 1:
    y[row_index[k]] += scalar * values[k]
```

### 1.2 CSR (Compressed Sparse Row) Format

Identical structure but row-oriented: `row_start[num_rows+1]`, `col_index[nnz]`, `values[nnz]`.

**Why CSR:** Row-wise operations in presolve (Person 3), GPU SpMV (Person 4), and constraint evaluation.

**Conversion CSC → CSR:** This is equivalent to a sparse matrix transpose.
```
CSC_TO_CSR(csc):
  // Count entries per row
  row_count = array[num_rows], init to 0
  for k = 0 to nnz - 1:
    row_count[csc.row_index[k]] += 1

  // Prefix sum → row_start
  row_start[0] = 0
  for i = 0 to num_rows - 1:
    row_start[i + 1] = row_start[i] + row_count[i]

  // Fill CSR arrays
  cursor = copy of row_start[0..num_rows-1]
  for j = 0 to num_cols - 1:
    for k = csc.col_start[j] to csc.col_start[j+1] - 1:
      i = csc.row_index[k]
      pos = cursor[i]
      col_index[pos] = j
      csr_values[pos] = csc.values[k]
      cursor[i] += 1

  return CSR(num_rows, num_cols, row_start, col_index, csr_values)
```

### 1.3 Performance Requirements
- SpMV must achieve >1 GFLOP/s on a single core for typical LP matrices (~10 nnz/row)
- No unnecessary copies — all operations work on raw array pointers
- Prefer `restrict` pointers and loop vectorization hints where applicable

---

## 2. Sparse LU Factorization

**This is the MOST CRITICAL component.** The simplex method's inner loop is dominated by LU solve time. Get this wrong and the solver is unusable.

### 2.1 Algorithm: Sparse LU with Markowitz Ordering and Threshold Pivoting

Given basis matrix B (m × m, sparse), compute P, Q, L, U such that **P B Q = L U** where:
- P, Q are permutation matrices (row and column pivoting)
- L is unit lower triangular (stored without diagonal)
- U is upper triangular

```
SPARSE_LU_FACTORIZE(B, threshold = 0.1):
  m = B.num_rows
  // Work on a mutable copy of B in linked-list column + row format
  // for efficient dynamic element insertion/deletion

  // Initialize row counts and column counts
  row_nnz[i] = number of nonzeros in row i of B,  for i = 0..m-1
  col_nnz[j] = number of nonzeros in col j of B,  for j = 0..m-1

  row_perm = identity(m)     // P: maps pivot step → original row
  col_perm = identity(m)     // Q: maps pivot step → original column

  L_entries = empty list      // will store (step, row, value) triples
  U_entries = empty list      // will store (step, col, value) triples

  for step = 0 to m - 1:
    // ── PIVOT SELECTION (Markowitz + threshold) ──
    best_cost = INFINITY
    pivot_row = -1
    pivot_col = -1

    // Search candidates: iterate over active submatrix
    for each active column j (not yet pivoted):
      // Find maximum absolute value in column j (active part)
      max_abs_col = max(|a_ij|) for active i in column j
      if max_abs_col < EPS_PIVOT:
        continue  // singular or near-singular column

      for each active row i with nonzero a_ij in column j:
        // Threshold test: pivot must be large enough relative to column max
        if |a_ij| < threshold * max_abs_col:
          continue

        // Markowitz cost
        cost = (row_nnz[i] - 1) * (col_nnz[j] - 1)
        if cost < best_cost:
          best_cost = cost
          pivot_row = i
          pivot_col = j

        if cost == 0:
          break  // can't do better than 0

      if best_cost == 0:
        break

    if pivot_row == -1:
      // No valid pivot found → matrix is numerically singular
      return SINGULAR_ERROR

    // Record permutations
    row_perm[step] = pivot_row
    col_perm[step] = pivot_col

    pivot_value = B[pivot_row][pivot_col]

    // ── STORE U ROW (pivot row entries in active columns) ──
    for each active column j with nonzero B[pivot_row][j]:
      U_entries.append(step, j, B[pivot_row][j])

    // ── COMPUTE L COLUMN AND ELIMINATE ──
    for each active row i ≠ pivot_row with nonzero B[i][pivot_col]:
      multiplier = B[i][pivot_col] / pivot_value
      L_entries.append(step, i, multiplier)

      // Row elimination: row_i -= multiplier * pivot_row
      for each active column j with nonzero B[pivot_row][j]:
        B[i][j] -= multiplier * B[pivot_row][j]
        if |B[i][j]| < EPS_ZERO:
          // Drop to zero — remove from sparse structure
          remove B[i][j]
          row_nnz[i] -= 1
          col_nnz[j] -= 1
        elif B[i][j] was zero before:
          // Fill-in — new nonzero created
          row_nnz[i] += 1
          col_nnz[j] += 1

      // Remove eliminated entry
      remove B[i][pivot_col]
      row_nnz[i] -= 1

    // Mark pivot_row and pivot_col as eliminated
    col_nnz[pivot_col] = 0
    row_nnz[pivot_row] = 0

  // Build compressed L and U from entries
  L = build_unit_lower_triangular(L_entries, m)
  U = build_upper_triangular(U_entries, m)

  return (L, U, row_perm, col_perm)
```

**Critical implementation notes:**
- Use a **doubly-linked list** or **hash map** per row and column for O(1) insertion/deletion during elimination
- Pre-allocate memory for fill-in (estimate 5-10× original nnz for typical LP bases)
- The `threshold` parameter trades sparsity vs. stability: lower = sparser but less stable. Start at 0.1, retry at 0.5 or 1.0 if growth factor is too large.

### 2.2 Solving with LU Factors

Two critical operations used every simplex iteration:

**FTRAN: Solve Bx = b (used to compute pivot column α = B⁻¹aₑ)**
```
FTRAN(L, U, row_perm, col_perm, b):
  // Apply row permutation: b' = P * b
  for k = 0 to m - 1:
    work[k] = b[row_perm[k]]

  // Forward solve: L * z = b'  (L is unit lower triangular)
  for k = 0 to m - 1:
    // z[k] = work[k] is already correct (unit diagonal)
    for each (row_i, l_val) in L column k where row_i > k:
      work[row_i] -= l_val * work[k]

  // Backward solve: U * x' = z
  for k = m - 1 downto 0:
    if |U[k][k]| < EPS_PIVOT:
      return NUMERICAL_ERROR
    work[k] /= U[k][k]
    for each (col_j, u_val) in U row k where col_j < k:
      work[col_j] -= u_val * work[k]

  // Apply column permutation: x = Q * x'
  for k = 0 to m - 1:
    result[col_perm[k]] = work[k]

  return result
```

**BTRAN: Solve Bᵀy = c (used to compute dual variables y = B⁻ᵀcB)**
```
BTRAN(L, U, row_perm, col_perm, c):
  // Apply column permutation inverse: c' = Q^T * c
  for k = 0 to m - 1:
    work[k] = c[col_perm[k]]

  // Forward solve: U^T * z = c'
  for k = 0 to m - 1:
    for each (col_j, u_val) in U row k where col_j < k:
      work[k] -= u_val * work[col_j]   // Note: U^T, so column ops become row ops
    work[k] /= U[k][k]

  // Backward solve: L^T * y' = z  (L^T is unit upper triangular)
  for k = m - 1 downto 0:
    for each (row_i, l_val) in L column k where row_i > k:
      work[k] -= l_val * work[row_i]

  // Apply row permutation inverse: y = P^T * y'
  for k = 0 to m - 1:
    result[row_perm[k]] = work[k]

  return result
```

**Sparse FTRAN/BTRAN:** When the RHS vector `b` is sparse (e.g., a single column of A), exploit this:
```
SPARSE_FTRAN(L, U, row_perm, col_perm, b_sparse):
  // b_sparse has only a few nonzeros
  // 1. Compute "reach" — the set of L columns that are affected
  reach = topological_sort(nonzero_indices_of_Pb, elimination_tree_of_L)
  // 2. Only process L columns in the reach set (in topological order)
  // 3. Only process U rows that have nonzero entries after L solve
  // This reduces O(m) work to O(nnz_in_solution) — crucial for steepest edge
```

### 2.3 Numerical Stability Measures

1. **Growth factor monitoring:**
   ```
   growth = max(|U_ij|) / max(|B_ij|)
   if growth > 1e8:
     refactorize with threshold = 0.5  // tighter pivoting
     if still > 1e8: threshold = 1.0   // full partial pivoting
   ```

2. **Iterative refinement:**
   ```
   ITERATIVE_REFINE(B, x, b, max_rounds = 3):
     for round = 1 to max_rounds:
       r = b - B * x                    // compute residual
       if ||r||_inf < EPS_FEASIBILITY * ||b||_inf:
         break                           // accurate enough
       d = FTRAN(L, U, perms, r)        // solve B*d = r
       x += d                            // correct solution
   ```

3. **Pivot rejection:** If `|pivot| < EPS_PIVOT`, widen the search to accept any sufficiently large element in the active submatrix, even if Markowitz cost is higher.

---

## 3. LU Update (Forrest-Tomlin)

When the basis changes by one column (column `p` leaves, column `q` enters), update LU instead of refactorizing from scratch.

### 3.1 Forrest-Tomlin Update Algorithm

```
FORREST_TOMLIN_UPDATE(L, U, row_perm, col_perm, leaving_pos, entering_col):
  // leaving_pos = position in basis of leaving variable
  // entering_col = a_q (the entering column from the constraint matrix)

  // Step 1: Compute spike = B^{-1} * entering_col via FTRAN
  spike = FTRAN(L, U, row_perm, col_perm, entering_col)

  // Step 2: spike replaces column 'leaving_pos' of the identity in B^{-1}
  // In U, this means column leaving_pos is replaced with spike values
  // We need to restore upper triangular form.

  // Step 3: Permute spike so leaving_pos moves to the last position
  // Shift columns leaving_pos+1, ..., m-1 one position left
  // The spike column is now in position m-1

  // Rearrange: save the spike values
  eta = spike[leaving_pos + 1 : m - 1]   // the subdiagonal elements
  spike_diag = spike[leaving_pos]          // new diagonal

  if |spike_diag| < EPS_PIVOT:
    return NUMERICAL_ERROR  // degenerate pivot

  // Step 4: Store the eta vector (update transformation)
  // The transformation is: U_new = E_k * U_shifted
  // where E_k is identity except column leaving_pos has the eta vector
  eta_vector = spike[0 : m-1]
  eta_vector[leaving_pos] = spike_diag

  // Step 5: Store update transformation
  update_etas.append({
    position: leaving_pos,
    eta: eta_vector,
    permutation: <the column shift performed>
  })

  num_updates += 1
```

**Applying updates during FTRAN/BTRAN:**
After the base LU solve, apply each update transformation in sequence:
```
FTRAN_WITH_UPDATES(b):
  x = FTRAN(L, U, perms, b)       // base factorization
  for each update eta in update_etas (oldest to newest):
    // Apply permutation
    apply eta.permutation to x
    // Apply eta transformation
    pivot_val = x[eta.position]
    for i ≠ eta.position:
      x[i] -= eta.eta[i] * pivot_val / eta.eta[eta.position]
    x[eta.position] = pivot_val / eta.eta[eta.position]
  return x
```

### 3.2 Refactorization Trigger

Refactorize the basis from scratch when ANY of:
- `num_updates > max(50, sqrt(num_cols))` — too many update transformations
- Fill-in in update etas exceeds 5× the original LU fill-in
- Residual check fails: `||Bx - b||_inf / ||b||_inf > 1e-8` after an FTRAN
- Growth factor (max update eta value) exceeds 1e10

```
SHOULD_REFACTORIZE():
  if num_updates > refactor_interval:
    return true
  if total_update_nnz > 5 * base_lu_nnz:
    return true
  // Periodic accuracy check (every 20 iterations)
  if iter_count % 20 == 0:
    x = FTRAN(L, U, updates, b)
    residual = ||B*x - b||_inf / max(1.0, ||b||_inf)
    if residual > 1e-8:
      return true
  return false
```

---

## 4. Basis Management

### 4.1 Basis Representation

```cpp
// For each variable (column), track its status:
enum class BasisStatus { BASIC, AT_LOWER, AT_UPPER, FIXED, FREE_ZERO };

struct Basis {
    std::vector<Index> basic_vars;      // size m: basic_vars[i] = column index of i-th basic variable
    std::vector<BasisStatus> var_status; // size n: status of every variable
    std::vector<Index> basic_row;       // size n: if var j is basic, basic_row[j] = its position in basis
};
```

### 4.2 Crash Procedure (Triangularity-Based)

Find a good initial basis without artificial variables when possible:

```
CRASH(A, row_lower, row_upper, col_lower, col_upper):
  m = num_rows, n = num_cols
  basis = empty list (capacity m)
  assigned_row = array[m], all false
  assigned_col = array[n], all false

  // Pass 1: Singleton columns — columns with exactly 1 nonzero in unassigned rows
  for j = 0 to n - 1:
    if col is slack/surplus: continue
    count = 0; last_row = -1
    for each (i, val) in column j:
      if not assigned_row[i]:
        count += 1; last_row = i
    if count == 1 and |val| > EPS_PIVOT:
      basis.append(j)
      assigned_row[last_row] = true
      assigned_col[j] = true

  // Pass 2: Columns with 2 nonzeros in unassigned rows, etc.
  for pass = 2 to min(5, m):
    for j = 0 to n - 1:
      if assigned_col[j]: continue
      count = 0; best_row = -1; best_val = 0
      for each (i, val) in column j:
        if not assigned_row[i]:
          count += 1
          if |val| > |best_val|:
            best_row = i; best_val = val
      if count == pass and |best_val| > EPS_PIVOT:
        basis.append(j)
        assigned_row[best_row] = true
        assigned_col[j] = true

  // Fill remaining basis positions with logical (slack) variables
  for i = 0 to m - 1:
    if not assigned_row[i]:
      // Add slack variable for row i (column index n + i)
      basis.append(n + i)
      assigned_row[i] = true

  return basis
```

---

## 5. Pricing Strategies

### 5.1 Dantzig Pricing (Full Pricing)
```
DANTZIG_PRICE(reduced_costs, var_status, sense):
  best_j = -1
  best_rc = -EPS_OPTIMALITY   // for minimization

  for j = 0 to n - 1:
    if var_status[j] == BASIC: continue
    rc = reduced_costs[j]
    // For AT_LOWER: enters if rc < 0 (min), variable increases
    // For AT_UPPER: enters if rc > 0 (min), variable decreases
    if var_status[j] == AT_LOWER and rc < best_rc:
      best_rc = rc
      best_j = j
    elif var_status[j] == AT_UPPER and -rc < best_rc:
      best_rc = -rc
      best_j = j

  return best_j   // -1 means optimal
```

### 5.2 Steepest Edge Pricing

Maintains weights `w[j] = ||B⁻¹ a_j||²` for each non-basic variable.

```
STEEPEST_EDGE_INIT(B_inv, A):
  // Expensive initialization — compute exact weights
  for each non-basic j:
    col_j = A.column(j)
    eta_j = FTRAN(col_j)     // B^{-1} a_j
    weights[j] = dot(eta_j, eta_j)

STEEPEST_EDGE_PRICE(reduced_costs, weights, var_status):
  best_j = -1
  best_score = -EPS_OPTIMALITY

  for each non-basic j:
    rc = reduced_costs[j]
    if eligible_to_enter(j, rc, var_status):
      score = (rc * rc) / weights[j]    // normalized: most negative rc per unit move
      if score > best_score:
        best_score = score
        best_j = j

  return best_j

STEEPEST_EDGE_UPDATE(weights, pivot_col_alpha, leaving_pos, A):
  // After a pivot, update all non-basic weights
  // pivot_col_alpha = B^{-1} a_q (the FTRAN result of the entering column)
  // tau = B^{-T} e_p (BTRAN of the leaving row unit vector)

  tau = BTRAN(e_{leaving_pos})   // row of B^{-1}
  pivot_elem = pivot_col_alpha[leaving_pos]

  for each non-basic j (including the variable that just left):
    // tau_j = tau^T a_j (compute via sparse dot product)
    tau_j = dot_col(A, j, tau)
    alpha_p = pivot_col_alpha[leaving_pos]

    // DSE update formula:
    // w_j_new = w_j_old - 2*(tau_j/alpha_p)*alpha_j + (tau_j/alpha_p)^2 * w_entering
    // where alpha_j = dot_col(pivot_col_alpha, j)... simplified:
    ratio = tau_j / alpha_p
    weights[j] = max(EPS_ZERO, weights[j] - 2.0 * ratio * pivot_col_alpha_j + ratio * ratio * weights_entering)
```

### 5.3 Devex Pricing

Approximate steepest edge — cheaper per iteration, slightly more iterations total:

```
DEVEX_INIT():
  weights[j] = 1.0 for all non-basic j   // start with uniform weights
  reference_frame = 0

DEVEX_PRICE():
  // Same selection as steepest edge but with approximate weights
  // Score = rc^2 / devex_weight[j]

DEVEX_UPDATE(pivot_col_alpha, leaving_pos):
  pivot_elem = pivot_col_alpha[leaving_pos]
  for each non-basic j:
    alpha_j = pivot_col_alpha[j]  // FTRAN component
    // Devex update: simple approximate formula
    candidate_weight = 1.0 + (alpha_j / pivot_elem)^2
    weights[j] = max(0.01, max(candidate_weight, 0.5 * weights[j]))

  // Every ~100 iterations, reset weights to 1.0 (new reference frame)
  if iterations_since_reset > max(100, num_cols / 5):
    weights = all 1.0
    reference_frame += 1
```

---

## 6. Revised Primal Simplex Method

### 6.1 Complete Algorithm with Bounded Variables

```
PRIMAL_SIMPLEX(Problem P) -> Solution:
  m = P.num_rows
  n = P.num_cols

  // ── PHASE 0: SETUP ──
  scale = compute_scaling(P)               // Section 9
  apply_scaling(P, scale)
  basis = CRASH(P.A, P.row_lower, ...)     // Section 4.2
  var_status = initialize_status(P, basis)

  // Set non-basic variable values
  for j = 0 to n - 1:
    if var_status[j] == AT_LOWER: x[j] = col_lower[j]
    elif var_status[j] == AT_UPPER: x[j] = col_upper[j]
    elif var_status[j] == FIXED: x[j] = col_lower[j]
    elif var_status[j] == FREE_ZERO: x[j] = 0.0

  // ── PHASE 1: Find feasible basis (if crash basis is infeasible) ──
  // Factor basis
  (L, U, rperm, cperm) = SPARSE_LU_FACTORIZE(B)

  // Compute basic variable values: x_B = B^{-1}(b - N*x_N)
  rhs = P.b - A_N * x_N       // b minus non-basic contributions
  x_B = FTRAN(L, U, rperm, cperm, rhs)

  // Check primal feasibility
  if any x_B[i] < col_lower[basis[i]] - EPS_FEASIBILITY
     or x_B[i] > col_upper[basis[i]] + EPS_FEASIBILITY:
    // Use Big-M: add penalty for infeasibilities to objective
    // Alternatively: switch to Phase I objective (minimize sum of infeasibilities)
    phase1 = true
    save original objective c_orig = c
    c = phase1_objective(x_B, col_lower, col_upper, basis)
  else:
    phase1 = false

  // Initialize pricing (steepest edge by default)
  pricing = STEEPEST_EDGE_INIT() or DEVEX_INIT()
  num_updates = 0
  degenerate_count = 0

  // ── MAIN LOOP ──
  for iter = 0 to MAX_ITER - 1:
    if elapsed_time() > TIME_LIMIT:
      return Solution(TIME_LIMIT)

    // ── STEP 1: Compute dual variables ──
    c_B = [c[basis[i]] for i = 0..m-1]
    y = BTRAN(L, U, rperm, cperm, c_B)    // y = B^{-T} c_B

    // ── STEP 2: Compute reduced costs and price ──
    for each non-basic j:
      reduced_cost[j] = c[j] - DOT_COL(A, j, y)   // d_j = c_j - y^T a_j

    // ── STEP 3: Select entering variable ──
    q = pricing.select(reduced_cost, var_status)

    if q == -1:
      // All reduced costs non-negative → optimal (or Phase 1 complete)
      if phase1:
        if objective_value > EPS_FEASIBILITY:
          return Solution(INFEASIBLE)   // can't reach feasibility
        // Switch to Phase 2 with original objective
        c = c_orig
        phase1 = false
        continue   // restart main loop with original objective
      else:
        // OPTIMAL
        break

    // ── STEP 4: Compute pivot column ──
    a_q = extract_column(A, q)
    alpha = FTRAN(L, U, rperm, cperm, a_q)    // α = B^{-1} a_q

    // ── STEP 5: Ratio test (with bound flipping) ──
    // Determine entering direction
    if var_status[q] == AT_LOWER:
      direction = +1.0   // variable increases
    else:  // AT_UPPER
      direction = -1.0   // variable decreases
      alpha = -alpha      // flip direction

    theta_max = INF
    leaving_pos = -1

    for i = 0 to m - 1:
      if alpha[i] > EPS_PIVOT:
        // Basic variable i will decrease
        slack = x_B[i] - col_lower[basis[i]]
        ratio = slack / alpha[i]
      elif alpha[i] < -EPS_PIVOT:
        // Basic variable i will increase
        slack = col_upper[basis[i]] - x_B[i]
        ratio = slack / (-alpha[i])
      else:
        continue

      if ratio < theta_max:
        theta_max = ratio
        leaving_pos = i

    // ── BOUND FLIPPING ──
    // If theta_max is limited by the entering variable's opposite bound:
    entering_range = col_upper[q] - col_lower[q]
    if entering_range < theta_max:
      // Flip the entering variable to its other bound (no basis change!)
      if var_status[q] == AT_LOWER:
        x[q] = col_upper[q]
        var_status[q] = AT_UPPER
      else:
        x[q] = col_lower[q]
        var_status[q] = AT_LOWER
      // Update x_B -= entering_range * alpha
      x_B -= entering_range * alpha
      degenerate_count = 0
      continue  // no basis change needed

    if leaving_pos == -1:
      return Solution(UNBOUNDED)

    // ── STEP 6: Update solution ──
    x[q] = (var_status[q] == AT_LOWER ? col_lower[q] : col_upper[q]) + direction * theta_max
    old_x_B = x_B[leaving_pos]
    x_B -= theta_max * alpha
    x_B[leaving_pos] = direction * theta_max   // entering variable's basic value... 
    // Actually: the entering variable replaces the leaving variable in the basis

    // Check for degeneracy (zero pivot step)
    if theta_max < EPS_ZERO:
      degenerate_count += 1
    else:
      degenerate_count = 0

    // ── STEP 7: Anti-cycling ──
    if degenerate_count > 100:
      apply_perturbation(col_lower, col_upper, basis)
      degenerate_count = 0

    // ── STEP 8: Update basis ──
    leaving_var = basis[leaving_pos]
    basis[leaving_pos] = q
    // Determine leaving variable's new status
    if x_B_old ≈ col_lower[leaving_var]:
      var_status[leaving_var] = AT_LOWER
    else:
      var_status[leaving_var] = AT_UPPER
    var_status[q] = BASIC

    // ── STEP 9: Update LU factorization ──
    status = FORREST_TOMLIN_UPDATE(L, U, rperm, cperm, leaving_pos, a_q)
    num_updates += 1

    if status == NUMERICAL_ERROR or SHOULD_REFACTORIZE():
      (L, U, rperm, cperm) = SPARSE_LU_FACTORIZE(B_current)
      num_updates = 0
      // Re-compute x_B for numerical accuracy
      x_B = FTRAN(L, U, rperm, cperm, rhs)

    // ── STEP 10: Update pricing weights ──
    pricing.update(alpha, leaving_pos, ...)

  // ── POST-SOLVE ──
  unscale_solution(x, y, scale)
  return Solution(OPTIMAL, x, y, reduced_cost, ...)
```

---

## 7. Revised Dual Simplex Method

**The dual simplex is CRITICAL for MILP.** Person 3 calls `solve_from_basis()` at every B&B node. The child node differs from the parent by one bound change — dual simplex resolves this in very few iterations.

### 7.1 Complete Algorithm

```
DUAL_SIMPLEX(Problem P, optional warm_basis) -> Solution:
  m = P.num_rows
  n = P.num_cols

  // ── SETUP ──
  if warm_basis provided:
    basis = warm_basis
  else:
    basis = CRASH(...)

  (L, U, rperm, cperm) = SPARSE_LU_FACTORIZE(B)
  x_B = FTRAN(b_adjusted)

  // Compute initial dual variables and reduced costs
  y = BTRAN(c_B)
  for each non-basic j:
    rc[j] = c[j] - DOT_COL(A, j, y)

  // Ensure dual feasibility (all reduced costs correct sign)
  // For AT_LOWER: rc[j] >= 0 (for min)
  // For AT_UPPER: rc[j] <= 0 (for min)
  // If not dual feasible: fix by flipping non-basic variable bounds

  shift_for_dual_feasibility(rc, var_status)

  // ── MAIN LOOP ──
  for iter = 0 to MAX_ITER - 1:
    if elapsed_time() > TIME_LIMIT:
      return Solution(TIME_LIMIT)

    // ── STEP 1: Check primal feasibility → select leaving variable ──
    leaving_pos = -1
    max_infeas = EPS_FEASIBILITY

    for i = 0 to m - 1:
      bi = basis[i]
      below = col_lower[bi] - x_B[i]   // how much below lower bound
      above = x_B[i] - col_upper[bi]    // how much above upper bound
      infeas = max(below, above)
      if infeas > max_infeas:
        max_infeas = infeas
        leaving_pos = i
        if below > above:
          leaving_direction = -1   // variable needs to increase (bound from below)
        else:
          leaving_direction = +1   // variable needs to decrease (bound from above)

    if leaving_pos == -1:
      // All basic variables feasible → OPTIMAL
      return Solution(OPTIMAL, x_B, y, rc)

    // ── STEP 2: Compute pivot row ──
    // e_p = unit vector at leaving_pos
    e_p = zeros(m)
    e_p[leaving_pos] = 1.0
    tau = BTRAN(L, U, rperm, cperm, e_p)   // τ = B^{-T} e_p (row of B^{-1})

    // Compute pivot row: ρ_j = τ^T a_j for each non-basic j
    for each non-basic j:
      rho[j] = DOT_COL(A, j, tau)

    // ── STEP 3: Dual ratio test (Harris) ──
    // Select entering variable q that maintains dual feasibility
    best_q = -1
    best_ratio = INF
    HARRIS_TOL = 1e-7

    // Pass 1: Harris — find the maximum step that keeps all duals feasible
    // within a tolerance
    harris_theta = INF
    for each non-basic j:
      if (leaving_direction == -1 and rho[j] > EPS_PIVOT):
        // Variable at lower bound, rho positive
        ratio = (rc[j] + HARRIS_TOL) / rho[j]
        harris_theta = min(harris_theta, ratio)
      elif (leaving_direction == -1 and rho[j] < -EPS_PIVOT):
        // Variable at upper bound, rho negative
        ratio = (rc[j] - HARRIS_TOL) / rho[j]
        harris_theta = min(harris_theta, ratio)
      // Similar for leaving_direction == +1 with flipped signs
      elif (leaving_direction == +1 and rho[j] < -EPS_PIVOT):
        ratio = (-rc[j] + HARRIS_TOL) / (-rho[j])
        harris_theta = min(harris_theta, ratio)
      elif (leaving_direction == +1 and rho[j] > EPS_PIVOT):
        ratio = (-rc[j] - HARRIS_TOL) / (-rho[j])  // sign handling
        harris_theta = min(harris_theta, ratio)

    // Pass 2: Among candidates within Harris theta, pick largest |rho[j]|
    // (most numerically stable pivot)
    best_pivot = 0.0
    for each non-basic j:
      if is_candidate(j, rho[j], rc[j], harris_theta, leaving_direction):
        if |rho[j]| > best_pivot:
          best_pivot = |rho[j]|
          best_q = j

    if best_q == -1:
      return Solution(INFEASIBLE)  // dual unbounded → primal infeasible

    // ── STEP 4: Compute pivot column (for basis update) ──
    a_q = extract_column(A, best_q)
    alpha = FTRAN(L, U, rperm, cperm, a_q)
    pivot_element = alpha[leaving_pos]

    if |pivot_element| < EPS_PIVOT:
      // Numerical trouble — try another candidate
      return NUMERICAL_ERROR

    // ── STEP 5: Compute step sizes ──
    theta_dual = rc[best_q] / rho[best_q]   // dual step
    theta_primal = (x_B[leaving_pos] - bound) / pivot_element  // primal step

    // ── STEP 6: Update ──
    // Update primal (basic variable values)
    for i = 0 to m - 1:
      x_B[i] -= theta_primal * alpha[i]
    x_B[leaving_pos] = bound_of_entering   // entering variable's new value

    // Update dual (reduced costs)
    for each non-basic j:
      rc[j] -= theta_dual * rho[j]
    rc[best_q] = -theta_dual   // leaving variable gets the dual step... 
    // (exact formula depends on leaving direction)

    // Update dual variables
    y -= theta_dual * tau

    // ── STEP 7: Basis swap ──
    leaving_var = basis[leaving_pos]
    basis[leaving_pos] = best_q
    var_status[leaving_var] = (leaving_direction == -1) ? AT_LOWER : AT_UPPER
    var_status[best_q] = BASIC

    // ── STEP 8: Update LU ──
    FORREST_TOMLIN_UPDATE(L, U, rperm, cperm, leaving_pos, a_q)
    num_updates += 1

    if SHOULD_REFACTORIZE():
      REFACTORIZE(...)
      RECOMPUTE_x_B_and_rc(...)

  return Solution(ITERATION_LIMIT)
```

### 7.2 Warm Start for B&B (Critical for Person 3)

```
SOLVE_FROM_BASIS(Problem P, basis_indices, bound_changes) -> Solution:
  // This is called thousands of times during branch-and-bound
  // The basis is from the parent node; only 1-2 bounds have changed

  // 1. Set basis from parent
  basis = basis_indices

  // 2. Apply bound changes (e.g., x_j >= ceil(frac_val))
  for each (var, new_lb, new_ub) in bound_changes:
    col_lower[var] = new_lb
    col_upper[var] = new_ub
    // If var is non-basic at a bound that's now violated, update its value
    if var_status[var] == AT_LOWER and new_lb > x[var]:
      x[var] = new_lb
    if var_status[var] == AT_UPPER and new_ub < x[var]:
      x[var] = new_ub

  // 3. The basis may now be primal infeasible (but is still dual feasible)
  //    → Run dual simplex — should resolve in very few iterations (typically 1-5)
  return DUAL_SIMPLEX(P, warm_basis = basis)
```

---

## 8. Anti-Cycling & Degeneracy Handling

### 8.1 Perturbation
```
APPLY_PERTURBATION(col_lower, col_upper, basis):
  // Shift lower bounds of basic variables by small random amounts
  for i = 0 to m - 1:
    j = basis[i]
    eps = EPS_FEASIBILITY * (1.0 + random(0, 1)) * (1 + i)
    col_lower[j] -= eps
    // After solving, round solution back to unperturbed bounds
```

### 8.2 Bland's Rule
When `degenerate_count > 200` and perturbation doesn't help:
```
BLAND_PRICE(reduced_costs, var_status):
  // Select the lowest-index eligible variable (guaranteed finite)
  for j = 0 to n - 1:
    if eligible_to_enter(j, reduced_costs[j], var_status[j]):
      return j
  return -1  // optimal
```

---

## 9. Scaling

Apply before solving, unscale after:

```
GEOMETRIC_MEAN_SCALING(A, max_rounds = 10):
  row_scale = array[m], all 1.0
  col_scale = array[n], all 1.0

  for round = 1 to max_rounds:
    max_change = 0.0

    // Row scaling: for each row, scale so geometric mean of |elements| ≈ 1
    for i = 0 to m - 1:
      min_abs = INF; max_abs = 0
      for each (j, val) in row i:
        abs_val = |val * col_scale[j] * row_scale[i]|
        if abs_val > EPS_ZERO:
          min_abs = min(min_abs, abs_val)
          max_abs = max(max_abs, abs_val)
      if max_abs > EPS_ZERO:
        factor = 1.0 / sqrt(min_abs * max_abs)
        row_scale[i] *= factor
        max_change = max(max_change, |factor - 1.0|)

    // Column scaling: same approach
    for j = 0 to n - 1:
      min_abs = INF; max_abs = 0
      for each (i, val) in column j:
        abs_val = |val * col_scale[j] * row_scale[i]|
        if abs_val > EPS_ZERO:
          min_abs = min(min_abs, abs_val)
          max_abs = max(max_abs, abs_val)
      if max_abs > EPS_ZERO:
        factor = 1.0 / sqrt(min_abs * max_abs)
        col_scale[j] *= factor
        max_change = max(max_change, |factor - 1.0|)

    if max_change < 0.01:
      break  // converged

  // Apply: A_scaled[i][j] = row_scale[i] * A[i][j] * col_scale[j]
  // Scale objective: c_scaled[j] = c[j] * col_scale[j]
  // Scale bounds: b_scaled[i] = b[i] * row_scale[i]
  // Scale variable bounds: lb_scaled[j] = lb[j] / col_scale[j], etc.
  return (row_scale, col_scale)

UNSCALE_SOLUTION(x, y, rc, row_scale, col_scale):
  x[j] *= col_scale[j]       // unscale primal
  y[i] *= row_scale[i]       // unscale dual
  rc[j] /= col_scale[j]      // unscale reduced costs
```

---

## 10. Testing Strategy

### Unit Tests (write these FIRST — TDD)

| Test | What to Verify | Input |
|------|---------------|-------|
| `test_csc_from_triplets` | Correct CSC arrays, duplicate summing | 5×5 matrix with duplicates |
| `test_csc_spmv` | `y = Ax` matches dense multiply | Random 10×10 sparse matrix |
| `test_csc_transpose_spmv` | `y = Aᵀx` matches dense | Same matrix |
| `test_csc_to_csr_roundtrip` | CSC→CSR→CSC gives same matrix | Various sizes |
| `test_lu_small` | `||PAQ - LU||_F < 1e-12` | 5×5, 10×10 known matrices |
| `test_lu_ftran` | `||Bx - b|| < 1e-10` | Random 20×20 basis |
| `test_lu_btran` | `||Bᵀy - c|| < 1e-10` | Same basis |
| `test_lu_update` | After 50 updates, FTRAN still accurate | Sequential pivots |
| `test_lu_singular` | Detect singular matrix | Rank-deficient 5×5 |
| `test_primal_simplex_tiny` | Known optimal for 2×3 LP | Hand-crafted |
| `test_dual_simplex_tiny` | Same LP, same optimal | Hand-crafted |
| `test_cycling` | Terminates on Beale's example | Beale's 3×7 cycling LP |
| `test_unbounded` | Returns UNBOUNDED | min -x, x ≥ 0 |
| `test_infeasible` | Returns INFEASIBLE | x ≤ 1 AND x ≥ 2 |
| `test_warm_start` | Few iterations from warm basis | Add bound, re-solve |

### Integration Tests (Netlib)

| Problem | Rows | Cols | Difficulty | Expected Optimal |
|---------|------|------|-----------|-----------------|
| `afiro` | 28 | 32 | Trivial | -464.7531 |
| `adlittle` | 57 | 97 | Easy | 225494.9 |
| `blend` | 75 | 83 | Moderate | -30.8121 |
| `sc205` | 206 | 203 | Moderate | -52.2021 |
| `share2b` | 97 | 79 | Numerical | -415.7325 |
| `25fv47` | 822 | 1571 | Performance | 5501.846 |
| `fit2p` | 3000 | 13525 | Large | 68464.293 |

**Pass criteria:** Objective within 1e-6 relative error of known optimal.

### Performance Targets

- `afiro`: < 0.01s
- `25fv47`: < 5s
- `fit2p`: < 60s
- Refactorization: < 2× cold factorization time

---

## 11. Key References

1. Maros, I. *Computational Techniques of the Simplex Method* (Springer, 2003) — **primary reference**
2. Suhl, U. & Suhl, L. "Computing Sparse LU Factorizations for Large-Scale LP" (ORSA J. Computing, 1990)
3. Koberstein, A. *The Dual Simplex Method, Techniques for a Fast and Stable Implementation* (PhD thesis, Paderborn, 2005) — **best dual simplex reference**
4. Goldfarb, D. & Forrest, J. "Steepest-edge simplex algorithms for LP" (Math. Programming, 1992)
5. Forrest, J. & Goldfarb, D. "Steepest-edge simplex algorithms" (Math. Programming, 1992)
6. Davis, T.A. *Direct Methods for Sparse Linear Systems* (SIAM, 2006) — sparse LU details
7. Huangfu, Q. & Hall, J.A.J. "Parallelizing the dual revised simplex method" (Math. Prog. Comp., 2018) — HiGHS implementation insights

---

## 12. Dependencies on Other Persons

| Who | What They Need From You | Priority |
|-----|------------------------|----------|
| **Person 2** | `SparseMatrixCSC`, `SparseMatrixCSR`, `TripletMatrix` | Week 1 (interface header) |
| **Person 3** | `DualSimplexSolver::solve_from_basis()` — called thousands of times in B&B. MUST be fast from warm start. | Week 6 (functional) |
| **Person 4** | `PrimalSimplexSolver::solve()`, `DualSimplexSolver::solve()` via `LPSolver` interface | Week 6 (functional) |
| **Everyone** | Sparse matrix correctness and performance | Week 2 (tested) |

---

## 13. Definition of Done

- [ ] All unit tests pass (100% of tests listed above)
- [ ] Solves **90+ of 95** Netlib LP problems correctly (optimal within 1e-6 of known)
- [ ] Dual simplex warm start resolves in ≤ 10 iterations for single bound change
- [ ] No memory leaks (Valgrind clean)
- [ ] Doxygen documentation for all public APIs
- [ ] Code compiles with `-Wall -Wextra -Werror` with zero warnings
- [ ] Performance targets met (see Section 10)
