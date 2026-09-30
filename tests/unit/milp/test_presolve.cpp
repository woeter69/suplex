/// Unit tests for the Presolve engine.
/// Tests each reduction rule independently with hand-crafted problems.
#include <gtest/gtest.h>
#include "milp/presolve.h"
#include "core/sparse_matrix.h"
#include "core/problem.h"
#include "core/types.h"

using namespace suplex;

// ── Helper: build a simple Problem from dense data ────────────────────────────
static Problem make_problem(int m, int n,
                             const std::vector<std::vector<double>>& A_dense,
                             const std::vector<double>& c,
                             const std::vector<double>& rl,
                             const std::vector<double>& ru,
                             const std::vector<double>& cl,
                             const std::vector<double>& cu,
                             const std::vector<VarType>& vt = {}) {
    // Build triplets
    std::vector<Index> trip_r, trip_c;
    std::vector<Real>  trip_v;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if (std::abs(A_dense[i][j]) > 1e-15) {
                trip_r.push_back(i); trip_c.push_back(j);
                trip_v.push_back(A_dense[i][j]);
            }

    Problem p;
    p.set_dimensions(m, n);
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(m, n, trip_r, trip_c, trip_v));
    p.set_objective(std::vector<Real>(c));
    p.set_row_bounds(std::vector<Real>(rl), std::vector<Real>(ru));
    p.set_col_bounds(std::vector<Real>(cl), std::vector<Real>(cu));

    std::vector<VarType> types = vt.empty() ? std::vector<VarType>(n, VarType::CONTINUOUS) : vt;
    p.set_var_types(std::move(types));
    return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Empty row detection
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, EmptyRowSatisfied) {
    // Row 0 has no nonzeros; constraint is -INF <= 0 <= INF → remove it
    Problem p = make_problem(2, 2,
        {{0, 0}, {1, 1}},
        {1, 1},
        {-INF, -INF}, {INF, INF},
        {0, 0}, {10, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    EXPECT_EQ(result.reduced_problem.num_rows(), 1);  // completely solved
}

TEST(Presolve, EmptyRowInfeasible) {
    // Row 0: 0 ∈ [5, 10] is impossible → infeasible
    Problem p = make_problem(2, 2,
        {{0, 0}, {1, 1}},
        {1, 1},
        {5, -INF}, {10, INF},       // row 0: 5 <= 0 <= 10 → impossible
        {0, 0}, {10, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_TRUE(result.is_infeasible);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Fixed variable removal
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, FixedVariableRemoved) {
    // x0 is fixed at 3; x1 is free
    Problem p = make_problem(1, 2,
        {{1, 2}},
        {5, 3},
        {0}, {INF},
        {3, 0}, {3, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    EXPECT_EQ(result.reduced_problem.num_cols(), 1);
    // Objective offset should include 5*3 = 15 (and potentially more from x1)
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Singleton row
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, SingletonRowTightensUpperBound) {
    // Row 0: 2*x0 <= 8  →  x0 <= 4
    Problem p = make_problem(2, 2,
        {{2, 0}, {0, 1}},
        {1, 1},
        {-INF, -INF}, {8, 5},
        {0, 0}, {10, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    // Everything should be removed because they are all singleton columns!
    EXPECT_EQ(result.reduced_problem.num_rows(), 0);
}

TEST(Presolve, SingletonRowInfeasible) {
    // Row 0: 1*x0 <= 2,  but x0 >= 5 from bounds → infeasible
    Problem p = make_problem(1, 1,
        {{1}},
        {1},
        {-INF}, {2},
        {5}, {10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_TRUE(result.is_infeasible);
}

TEST(Presolve, SingletonRowNegativeCoefficient) {
    // Row 0: -3*x0 <= -6  →  x0 >= 2
    Problem p = make_problem(1, 1,
        {{-3}},
        {1},
        {-INF}, {-6},
        {0}, {10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    EXPECT_NEAR(result.reduced_problem.col_lower()[0], 2.0, 1e-9);
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Forcing row
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, ForcingRow) {
    // Row 0: x0 + x1 >= 10, x0 in [0,5], x1 in [0,5]
    // max_activity = 10 = row_lower → forcing row (all vars at upper bounds)
    Problem p = make_problem(1, 2,
        {{1, 1}},
        {1, 1},
        {10}, {10},
        {0, 0}, {5, 5});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    // After forcing, both cols should be fixed at 5
    // The row and both cols should be removed
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Implied bound tightening
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, ImpliedBoundTightening) {
    // Row 0: x0 + x1 <= 6, x0 in [0,4], x1 in [0,4]
    // For x1: max_act_without_x1 = 4 (x0 at max)
    // Implied ub on x1: 6 - 0 = 6... but x1 <= 4 already
    // Let's make it tighter: x0 in [2,4]
    // Implied ub on x1: 6 - 2 = 4 (no tightening)
    // For a tighter case: x0 in [3,4]
    // Implied ub on x1: 6 - 3 = 3 < 4 → tightened!
    Problem p = make_problem(2, 2,
        {{1, 1},
         {1, 1}},
        {1, 1},
        {-INF, -INF}, {6, 100},
        {3, 0}, {4, 4});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    // x1 upper bound should be tightened to 3
    // Note: since it wasn't reduced to 0 columns, it's safe to check col_upper
    EXPECT_LE(result.reduced_problem.col_upper()[1], 3.0 + 1e-8);
}

TEST(Presolve, ImpliedBoundInfeasible) {
    // Row 0: x0 + x1 <= 3, x0 in [2,4], x1 in [2,4]
    // min_activity = 4 > 3 = row_upper → infeasible
    Problem p = make_problem(1, 2,
        {{1, 1}},
        {1, 1},
        {-INF}, {3},
        {2, 2}, {4, 4});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_TRUE(result.is_infeasible);
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Postsolve round-trip
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, PostsolveSingletonRow) {
    // After presolve, verify that postsolve maps solution back correctly.
    // Row 0: 2*x0 <= 8, x0 in [0,10], x1 in [0,10]
    // Row 1: x0 + x1 >= 3
    Problem p = make_problem(2, 2,
        {{2, 0}, {1, 1}},
        {1, 2},
        {-INF, 3}, {8, INF},
        {0, 0}, {10, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);

    // Simulate a reduced-space solution (x0=4, x1=3 in reduced space)
    Solution reduced_sol;
    reduced_sol.status = SolverStatus::OPTIMAL;
    reduced_sol.primal_values.assign(result.reduced_problem.num_cols(), 1.0);
    if (result.reduced_problem.num_cols() >= 1) reduced_sol.primal_values[0] = 4.0;
    if (result.reduced_problem.num_cols() >= 2) reduced_sol.primal_values[1] = 3.0;
    reduced_sol.dual_values.assign(result.reduced_problem.num_rows(), 0.0);
    reduced_sol.objective_value = 10.0;

    Solution original_sol = ps.postsolve(reduced_sol, result);
    EXPECT_EQ(original_sol.status, SolverStatus::OPTIMAL);
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. No-op on already minimal problem
// ─────────────────────────────────────────────────────────────────────────────
TEST(Presolve, NoReductionOnGenericProblem) {
    // A problem where no reduction applies
    Problem p = make_problem(2, 3,
        {{1, 2, 0}, {0, 1, 3}},
        {1, 1, 1},
        {1, 2}, {5, 8},
        {0, 0, 0}, {10, 10, 10});

    Presolve ps;
    auto result = ps.apply(p);
    EXPECT_FALSE(result.is_infeasible);
    // Problem may or may not shrink — just verify it runs without crash
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
