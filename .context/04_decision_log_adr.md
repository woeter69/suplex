# Architectural Decision Records (ADR Log)

## ADR-001: Autonomous Person 1 Simplex Architecture & Clean Separation
- **Status**: Accepted
- **Context**: Person 1 is assigned Sparse Linear Algebra and Simplex Methods. While Person 4 will provide parsers and CLI later, Person 1's work must be self-contained, rigorously testable, and independently compilable.
- **Decision**:
  1. Implement complete core types (`types.h`, `sparse_matrix.h/cpp`, `problem.h/cpp`, `solution.h/cpp`) adhering strictly to `INTERFACES.md`.
  2. Implement sparse LU with Markowitz minimum-degree strategy and threshold pivoting (`0.1` default, upgradable to `0.5`/`1.0` if growth factor exceeds safety limits).
  3. Implement product-form / eta update mechanism with Forrest-Tomlin principles, including refactorization triggers based on iteration count and residual checks.
  4. Provide modular pricing hierarchy: `DantzigPricing`, `SteepestEdgePricing`, `DevexPricing`, `BlandPricing`.
  5. Provide `PrimalSimplexSolver` and `DualSimplexSolver` inheriting from `LPSolver`. Dual simplex explicitly supports `solve_from_basis` for subsequent MILP branch-and-bound integration.
- **Consequences**:
  - Person 1 deliverables can be verified with standalone C++ unit tests immediately without waiting for Person 4's parser or Person 3's MILP engine.
  - Zero external dependencies ensures maximum portability and compile speed.
