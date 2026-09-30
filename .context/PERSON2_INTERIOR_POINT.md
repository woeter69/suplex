# PERSON 2: Interior Point Methods & QP Solver

## Your Mission
You build the barrier/interior-point solver, which provides an alternative LP solving path (often faster for very large problems) and the QP solver for quadratic objectives. Your methods are numerically intensive — sparse Cholesky factorization is your bottleneck.

## Deliverables Checklist

- [ ] `src/ipm/ordering.h/cpp` — Fill-reducing orderings (AMD, minimum degree) (3 days)
- [ ] `src/ipm/cholesky.h/cpp` — Sparse Cholesky factorization (5 days)
- [ ] `src/ipm/ipm_solver.h/cpp` — Mehrotra predictor-corrector IPM for LP (7 days)
- [ ] `src/ipm/qp_solver.h/cpp` — Interior point method for convex QP (5 days)
- [ ] `src/ipm/crossover.h/cpp` — IPM solution → basic feasible solution (4 days)
- [ ] Unit tests (5 days)
- [ ] Benchmark validation (3 days)

Total: ~32 working days

---

## 1. Fill-Reducing Ordering

### 1.1 Why Ordering Matters
In IPM, we solve the normal equations: `(A D^2 A^T) Δy = rhs` where D is a diagonal scaling. The matrix `A D^2 A^T` is symmetric positive definite. Direct Cholesky factorization produces fill-in (new nonzeros). The ordering of rows/columns DRAMATICALLY affects fill-in.

Bad ordering → 100x more fill-in → 100x slower factorization → solver is useless.

### 1.2 Approximate Minimum Degree (AMD)
Implement the AMD algorithm:
1. Build the adjacency graph of `A*A^T` (don't form the matrix, just the graph)
2. Repeatedly eliminate the node with minimum degree
3. Use mass elimination and element absorption to speed up
4. Output: permutation vector P such that `P^T (A D^2 A^T) P` has low fill-in

**AMD Algorithm Pseudocode:**
```text
function AMD(A_pattern):
    // Input: A_pattern represents the sparsity pattern of the symmetric matrix M = A A^T
    n = A_pattern.size()
    
    // Initialize Quotient Graph structures
    nodes = {0, 1, ..., n-1} // Supervariables
    degrees = compute_initial_degrees(A_pattern)
    eliminated = array of size n, initialized to false
    permutation = array of size n
    perm_idx = 0
    
    // Elements represent eliminated supervariables (cliques)
    elements = empty set
    
    while perm_idx < n:
        // 1. Minimum Degree Selection
        min_deg = infinity
        u = -1
        for i in nodes if not eliminated[i]:
            if degrees[i] < min_deg:
                min_deg = degrees[i]
                u = i
                
        // 2. Mass Elimination (Indistinguishable nodes)
        // Find all nodes sharing exactly the same set of neighbors as u
        mass_nodes = find_indistinguishable_nodes(u)
        
        for node in mass_nodes:
            permutation[perm_idx] = node
            perm_idx += 1
            eliminated[node] = true
            
        // 3. Quotient Graph Update (Element Creation & Absorption)
        // Treat u (and mass_nodes) as a new element E_new
        E_new = create_element(u, mass_nodes)
        
        // Find all adjacent elements to E_new
        adjacent_elements = find_adjacent_elements(E_new)
        
        // Absorb adjacent elements into E_new to keep graph small
        for E in adjacent_elements:
            absorb_element(E_new, E)
            
        // 4. Update Degrees
        // Compute approximate degree for all uneliminated neighbors of E_new
        // Exact degree is costly, so we bound it:
        // deg_approx(x) = |A_x| + | \cup_{e \in E_x} (e \setminus {x}) |
        uneliminated_neighbors = get_uneliminated_boundary(E_new)
        for v in uneliminated_neighbors:
            degrees[v] = compute_approx_degree(v, E_new)
            
    return permutation
```

### 1.3 Interface
```cpp
class Ordering {
public:
    // Compute fill-reducing ordering for symmetric matrix defined by its lower triangle
    // Input: CSC representation of lower triangle
    // Output: permutation vector perm[i] = new index of original row/col i
    static std::vector<Index> amd(const SparseMatrixCSC& lower_triangle);
    
    // For the normal equations, we don't form A*A^T explicitly
    // Instead, compute ordering from the sparsity pattern of A
    static std::vector<Index> amd_ata(const SparseMatrixCSC& A);
};
```

---

## 2. Sparse Cholesky Factorization

This is YOUR most critical component. IPM performance is dominated by Cholesky solve time.

### 2.1 Algorithm: Left-Looking Sparse Cholesky

**Symbolic and Numerical Phase Pseudocode:**
```text
function SparseCholesky_Symbolic(M_lower, perm):
    // M_lower is the lower triangle of the permuted matrix P^T M P
    n = M_lower.rows
    parent = array of size n, initialized to -1
    col_counts = array of size n, initialized to 1 (diagonal)
    L_pattern = array of lists for each column
    
    // 1. Elimination Tree
    for j = 0 to n-1:
        for i in M_lower.col(j).indices:
            if i > j:
                p = j
                while parent[p] != -1 and parent[p] != i:
                    p = parent[p]
                if parent[p] == -1:
                    parent[p] = i
                    
    // 2. Symbolic Factorization (compute sparsity pattern of L)
    for j = 0 to n-1:
        L_pattern[j].add(j) // diagonal
        // To find pattern of L(:, j), traverse elimination tree from structural nonzeros
        visited = array of size n, initialized to false
        for i in M_lower.col(j).indices:
            if i > j:
                curr = i
                while curr != -1 and not visited[curr]:
                    visited[curr] = true
                    L_pattern[j].add(curr)
                    col_counts[j] += 1
                    curr = parent[curr]
                    
    // Sort indices in L_pattern for CSC
    for j = 0 to n-1: sort(L_pattern[j])
    
    return L_pattern, parent, col_counts

function SparseCholesky_Numeric(M_lower, L_pattern, delta):
    n = M_lower.rows
    L_values = allocate array based on sum of col_counts
    
    for j = 0 to n-1:
        // Initialize column j of L
        for i in L_pattern[j]:
            L[i, j] = (i == j) ? M_lower[j, j] + delta : M_lower[i, j]
            
        // Left-looking updates: subtract contributions from previously computed columns k < j
        for k in columns_updating_j(j, L_pattern): // Requires inverse pattern mapping
            L_j_k = L[j, k]
            for i in L_pattern[k]:
                if i >= j:
                    L[i, j] -= L[i, k] * L_j_k
                    
        // Check SPD / Apply Regularization
        if L[j, j] <= 1e-12:
            return FAILED_NOT_SPD
            
        // Scale column j
        L_jj_sqrt = sqrt(L[j, j])
        L[j, j] = L_jj_sqrt
        for i in L_pattern[j]:
            if i > j:
                L[i, j] /= L_jj_sqrt
                
    return L
```

### 2.2 Numerical Stability
- If diagonal element becomes negative or too small → matrix is not SPD
- For IPM: this can happen due to poor scaling. Apply regularization: add small delta to diagonal
- Regularization: `M' = M + δI` where `δ` starts at `1e-12`, doubles if factorization fails
- After successful factorization with regularization, perform iterative refinement.

### 2.3 Solving with Cholesky
- Forward solve: `Lz = Pb` (sparse triangular solve)
- Backward solve: `L^T y = z`
- Solution: `x = P^T y`
- Iterative refinement: compute residual `r = b - Mx`, solve `Md = r`, update `x += d`

### 2.4 Performance Targets
- For a 10,000 × 10,000 sparse SPD matrix with ~100,000 nonzeros:
  - Symbolic analysis: < 0.1 seconds
  - Numerical factorization: < 1 second
  - Each solve (forward + backward): < 0.01 seconds

---

## 3. Mehrotra Predictor-Corrector IPM for LP

### 3.1 Standard Form Conversion

**Standard Form Conversion Pseudocode:**
```text
function ConvertToAugmentedForm(Problem P):
    // min c^T x
    // s.t. row_lower <= Ax <= row_upper
    //      col_lower <= x <= col_upper
    
    A_new = []
    c_new = []
    
    // Variables mapping
    for j in 0 to P.num_cols-1:
        // Apply variable shifts: x_new = x - col_lower
        // Handling free variables: split into x+ and x-
        if P.col_lower[j] == -infinity and P.col_upper[j] == infinity:
            add_col(A_new, P.A[:, j])  // x+
            add_col(A_new, -P.A[:, j]) // x-
            c_new.push(P.c[j], -P.c[j])
        else:
            add_col(A_new, P.A[:, j])
            c_new.push(P.c[j])
            
    // Slack variables for constraints
    for i in 0 to P.num_rows-1:
        if P.row_lower[i] == P.row_upper[i]: // Equality
            // No slack needed
        elif P.row_lower[i] > -infinity and P.row_upper[i] == infinity: // >= constraint
            // Ax - s = b -> add slack with coeff -1
            add_col(A_new, -I[:, i])
            c_new.push(0)
        elif P.row_lower[i] == -infinity and P.row_upper[i] < infinity: // <= constraint
            // Ax + s = b -> add slack with coeff 1
            add_col(A_new, I[:, i])
            c_new.push(0)
        else: // Ranged constraint
            // Transform into equality + bounded slack
            add_col(A_new, I[:, i])
            c_new.push(0)
            // Note: ranged constraints also introduce upper bounds on slacks!
            
    return A_new, c_new
```

### 3.2 KKT System
The optimality conditions (KKT) for LP:
```
A^T y + s = c          (dual feasibility)
Ax = b                 (primal feasibility)
X S e = 0              (complementarity)
x, s >= 0              (non-negativity)
```
where `X = diag(x)`, `S = diag(s)`.

### 3.3 Newton Direction
Linearize the KKT conditions around current `(x, y, s)`:
```
[0    A^T   I ] [Δx]   [rc]     rc = c - A^T y - s
[A    0     0 ] [Δy] = [rb]     rb = b - Ax
[S    0     X ] [Δs]   [rμ]     rμ = σμe - XSe
```

### 3.4 Normal Equations
Eliminate Δx and Δs to get:
```
(A D^2 A^T) Δy = A D^2 (rc - X^{-1} rμ) + rb
where D^2 = X / S = diag(x_i / s_i)
```
Then:
```
Δs = rc - A^T Δy
Δx = S^{-1}(rμ - X Δs)
```

### 3.5 Mehrotra Predictor-Corrector

**Full Mehrotra IPM Pseudocode:**
```text
function MehrotraIPM(A, b, c):
    n = A.cols
    m = A.rows
    
    // 1. Initial Starting Point
    (x, y, s) = ComputeStartingPoint(A, b, c)
    
    for iter = 1 to MAX_ITER:
        // Compute residuals and duality gap
        rc = c - A^T * y - s
        rb = b - A * x
        μ = (x^T * s) / n
        
        // Convergence Check
        prim_feas = ||rb|| / (1 + ||b||)
        dual_feas = ||rc|| / (1 + ||c||)
        gap = abs(c^T * x - b^T * y) / (1 + abs(c^T * x))
        if prim_feas < EPS and dual_feas < EPS and gap < EPS:
            return OPTIMAL, x, y, s
            
        // Predictor Step (Affine Scaling, sigma = 0)
        D2 = diag(x) / diag(s)
        N = A * D2 * A^T
        
        // Factorize Normal Equations
        L = SparseCholesky(N, ordering, delta)
        if L == FAILED:
            increase delta, refactorize
            
        // Solve predictor normal equations
        rhs_aff = A * (D2 * rc) + rb // since rμ = -XSe here, term X^-1 rμ = -s, but properly formulated:
        // rhs_aff = A * D2 * rc + rb + A * (X^-1 * X * s) = A * (D2*rc + x) + rb
        // Wait, correct rhs for affine: 
        // Δy_aff = (A D^2 A^T)^{-1} (rb + A(S^{-1}(X rc + x .* s)))
        dy_aff = CholeskySolve(L, rhs_aff)
        ds_aff = rc - A^T * dy_aff
        dx_aff = -x - (x / s) .* ds_aff
        
        // Compute affine step length
        alpha_p_aff = max_step(x, dx_aff)
        alpha_d_aff = max_step(s, ds_aff)
        
        // Centering Parameter
        μ_aff = ((x + alpha_p_aff * dx_aff)^T * (s + alpha_d_aff * ds_aff)) / n
        sigma = (μ_aff / μ)^3
        
        // Corrector Step (Combined direction)
        // rμ = σμe - dx_aff .* ds_aff
        // Form new right hand side
        rhs_corr = A * (x / s) .* (rc - X^{-1} rμ) + rb
        
        dy = CholeskySolve(L, rhs_corr)
        ds = rc - A^T * dy
        dx = s^{-1} .* (rμ - x .* ds)
        
        // Final Step Length
        alpha_p = 0.9995 * max_step(x, dx)
        alpha_d = 0.9995 * max_step(s, ds)
        
        // Update Variables
        x = x + alpha_p * dx
        y = y + alpha_d * dy
        s = s + alpha_d * ds
        
    return UNCONVERGED

function max_step(v, dv):
    alpha = 1.0
    for i = 0 to v.size - 1:
        if dv[i] < 0:
            alpha = min(alpha, -v[i] / dv[i])
    return alpha
```

### 3.6 Starting Point
The starting point heuristic (Mehrotra's):
1. Solve `A A^T y = A c` to get `y`, then `s = c - A^T y`
2. Solve `A A^T x_hat = A b` to get an approximate `x`
3. Shift: `δx = max(-1.5 * min(x), 0)`, `δs = max(-1.5 * min(s), 0)`
4. `x += δx`, `s += δs`
5. Additional shift to balance `x^T s`

---

## 4. QP Interior Point Method

For convex QP: `min 0.5 x^T Q x + c^T x`, subject to `Ax ∈ [row_lower, row_upper]`, `x ∈ [col_lower, col_upper]`.

### 4.1 KKT System for QP
```
[Q    A^T   I ] [Δx]   [rQ]    rQ = Qx + c - A^T y - s
[A    0     0 ] [Δy] = [rb]    rb = b - Ax
[S    0     X ] [Δs]   [rμ]    rμ = σμe - XSe
```

### 4.2 Normal Equations for QP
Eliminate Δs, Δx:
```
(A (Q + Θ^{-1})^{-1} A^T) Δy = rhs
where Θ = S^{-1} X
```
Note: `(Q + Θ^{-1})` may not be easy to invert unless Q is diagonal. Alternative: use the augmented system via LDLT factorization, or conjugate gradient for inner solves. 

### 4.3 Algorithm
Same Mehrotra predictor-corrector structure as LP-IPM, but with Q in the KKT system.

**QP IPM Pseudocode Adaptation:**
```text
function QPIPM(Q, A, b, c):
    // Similar initialization
    // Q MUST be symmetric positive semidefinite
    
    for iter = 1 to MAX_ITER:
        rQ = Q*x + c - A^T * y - s
        rb = b - A * x
        μ = (x^T * s) / n
        
        // Predictor step
        Θ_inv = s / x
        M_qp = Q + diag(Θ_inv)
        
        // If Q is diagonal, M_qp is diagonal, easy inversion
        // Normal equations: A * M_qp^{-1} * A^T * dy = ...
        
        // If Q is not diagonal, solve Augmented System:
        // [ -M_qp    A^T ] [ dx ] = [ rhs_x ]
        // [   A       0  ] [ dy ] = [ rhs_y ]
        // Using sparse LDL^T
        
        // ... remainder of Mehrotra predictor-corrector is identical, substituting rQ for rc
```

### 4.4 Positive Semidefinite Check
Q must be PSD for convex QP. If not, the problem is non-convex (out of scope for now).
Simple check: attempt Cholesky of Q. If it fails, report error.

---

## 5. Crossover: IPM → Basic Feasible Solution

IPM produces a solution in the interior of the feasible region. For the simplex method (and for MILP branch-and-bound), we need a basic feasible solution (vertex).

### 5.1 Why Crossover
- IPM solutions have ALL variables strictly between bounds (no variable is exactly at a bound)
- We need to identify which variables are basic and push non-basic variables to their bounds
- This gives Person 1's simplex an excellent starting basis for warm-starting

### 5.2 Algorithm

**Crossover Pseudocode:**
```text
function Crossover(P, x_ipm):
    n = P.cols
    m = P.rows
    
    basis = empty array
    fixed_vars = empty array
    
    // 1. Classification
    for j = 0 to n-1:
        if abs(x_ipm[j] - P.col_lower[j]) < 1e-6:
            fix_variable(j, P.col_lower[j])
            fixed_vars.add(j)
        elif abs(x_ipm[j] - P.col_upper[j]) < 1e-6:
            fix_variable(j, P.col_upper[j])
            fixed_vars.add(j)
        else:
            basis.add(j)
            
    // A basis needs exactly m variables.
    // If |basis| > m, we have too many candidates. (Degeneracy or loose bounds)
    // If |basis| < m, we have too few.
    
    // Phase 1: Ensure structurally independent basis
    if basis.size() > m:
        // Perform LU factorization with pivoting on columns of basis
        // Keep columns that form a full rank matrix, fix the rest
        basis = extract_linearly_independent(A[:, basis])
        
    if basis.size() < m:
        // Add artificial variables or slack variables to complete basis
        basis.pad_to_m()
        
    // Phase 2: Simplex Refinement
    // Initialize Dual Simplex with `basis`
    simplex_solver = new DualSimplex(P)
    simplex_solver.set_basis(basis)
    simplex_solver.solve()
    
    return simplex_solver.get_solution() // Exact basic feasible solution
```

---

## 6. Iterative Refinement

After every solve (Cholesky or the full IPM), apply iterative refinement:
```text
for k = 1, 2, 3:
    r = b - A*x  // compute residual in extended precision if possible
    solve A*d = r
    x += d
    if ||r|| < ε * ||b||: break
```

---

## 7. Testing Strategy

### Unit Tests
- AMD ordering: verify on 10x10 matrix that fill-in is reduced vs natural ordering
- Cholesky: verify `L L^T = P M P^T` on small SPD matrices (5x5, 20x20)
- Cholesky: verify solve accuracy `||Mx - b|| < 1e-12`
- Cholesky: verify regularization activates on near-singular matrix
- IPM: solve 3x3 LP by hand, verify matches
- IPM: solve simple LP with equality, inequality, bound constraints
- QP: solve 2D QP (quadratic in 2 vars) — verify against analytical solution
- Crossover: verify that crossover produces a vertex solution

### Integration Tests (Netlib LP)
- `afiro`, `adlittle`, `blend`, `sc205`, `share2b` (same as Person 1, to verify IPM matches simplex)
- Larger: `stocfor2`, `ship08l`, `d2q06c`

### QP Tests (Maros-Mészáros)
- `QAFIRO` — trivial QP
- `QBRANDY` — moderate QP
- `QSCAGR7` — larger QP

---

## 8. Dependencies on Other Persons
- You implement the `LPSolver` interface (same as Person 1). Person 3 may use your IPM for root node relaxation.
- Person 4 provides the `Problem` object (parsed from MPS/LP files).
- Person 4's GPU SpMV could accelerate your Cholesky if the matrix is large enough. Expose hooks.
- Person 1's sparse matrix code is your foundation — you use `SparseMatrixCSC` throughout.

## 9. Key References
- Wright, S.J. "Primal-Dual Interior-Point Methods" (SIAM, 1997)
- Mehrotra, S. "On the Implementation of a Primal-Dual Interior Point Method" (SIAM J. Optimization, 1992)
- Andersen, E.D. & Andersen, K.D. "The MOSEK Interior Point Optimizer for LP" (2000)
- Davis, T.A. "Direct Methods for Sparse Linear Systems" (SIAM, 2006) — for Cholesky and AMD
- Gondzio, J. "Interior Point Methods 25 Years Later" (European J. of OR, 2012)

## 10. Definition of Done
- All unit tests pass
- IPM solves same Netlib problems as Person 1's simplex (cross-validation)
- QP solver handles at least 40 Maros-Mészáros problems
- Crossover produces valid basic solution
- No memory leaks
- Doxygen documentation complete
