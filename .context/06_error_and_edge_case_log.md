# Error & Edge Case Debugging Ledger

## Known Numerical & Edge Cases to Handle

1. **Singular Basis Matrix**:
   - Detection: No pivot found exceeding `EPS_PIVOT = 1e-10` during LU factorization.
   - Guardrail: Return `NUMERICAL_ERROR` / trigger recovery basis perturbation.

2. **Cycling & Stalling in Degenerate Pivots**:
   - Detection: Zero-length step sizes in ratio test (`theta < EPS_ZERO`) persisting over 100 iterations.
   - Guardrail: Apply bound perturbation and switch to Bland's minimum index rule.

3. **Infeasible Problems**:
   - Detection: Primal Phase 1 terminates with strictly positive artificial sum / dual simplex unboundedness.
   - Return: `SolverStatus::INFEASIBLE`.

4. **Unbounded Problems**:
   - Detection: Entering variable selected with all negative/zero pivot column coefficients (or entering ray unbounded).
   - Return: `SolverStatus::UNBOUNDED`.

5. **LU Fill-in & Factor Drift**:
   - Detection: Number of updates > 50 or residual `||Bx - b||_inf / ||b||_inf > 1e-8`.
   - Guardrail: Automatically trigger refactorization from scratch.
