/// Unit tests for Primal Heuristics.
#include <gtest/gtest.h>
#include "milp/heuristics.h"
#include "core/problem.h"
#include "core/solution.h"
#include "core/types.h"

using namespace suplex;

// ── Mock LP solver for testing ────────────────────────────────────────────────
class MockLPSolver : public LPSolver {
public:
    // Fixed solution to return
    Solution fixed_sol;
    SolverStatus fixed_status = SolverStatus::OPTIMAL;

    SolverStatus solve(const Problem&, Solution& sol) override {
        sol = fixed_sol;
        return fixed_status;
    }
    SolverStatus solve_from_basis(const Problem& p,
                                   const std::vector<Index>&,
                                   Solution& sol) override {
        return solve(p, sol);
    }
    void set_iteration_limit(Index) override {}
    void set_time_limit(Real)      override {}
    void set_log_level(LogLevel)   override {}
    const std::vector<Index>& get_basis() const override {
        static std::vector<Index> empty;
        return empty;
    }
    const std::vector<Real>& get_farkas_ray() const override {
        static std::vector<Real> empty;
        return empty;
    }
};

// ── Helper: build simple MIP ──────────────────────────────────────────────────
static Problem make_simple_mip() {
    // min x0 + x1
    // s.t. x0 + x1 >= 3
    //      x0, x1 in {0,1,...,5}
    Problem p;
    p.set_dimensions(1, 2);
    std::vector<Index> r={0,0}, c={0,1};
    std::vector<Real>  v={1,1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 2, r, c, v));
    p.set_objective({1, 1});
    p.set_row_bounds({3}, {INF});
    p.set_col_bounds({0, 0}, {5, 5});
    p.set_var_types({VarType::INTEGER, VarType::INTEGER});
    return p;
}

// ── Simple rounding ───────────────────────────────────────────────────────────
TEST(Heuristics, SimpleRoundingSuccess) {
    Problem p = make_simple_mip();

    // LP solution: x0=1.5, x1=1.6 → rounds to 2,2 which satisfies x0+x1>=3
    Solution lp_sol;
    lp_sol.primal_values = {1.5, 1.6};
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(nullptr);
    auto result = h.simple_rounding(p, lp_sol);

    ASSERT_TRUE(result.has_value());
    // Rounded: x0=2, x1=2, sum=4>=3 ✓
    EXPECT_NEAR(result->primal_values[0], 2.0, 0.5);
    EXPECT_NEAR(result->primal_values[1], 2.0, 0.5);
}

TEST(Heuristics, SimpleRoundingFailsWhenInfeasible) {
    // Rounding leads to infeasible: x0=0.4→0, x1=0.4→0, sum=0 < 3
    Problem p = make_simple_mip();

    Solution lp_sol;
    lp_sol.primal_values = {0.4, 0.4};
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(nullptr);
    auto result = h.simple_rounding(p, lp_sol);
    // Either fails or produces infeasible solution
    if (result.has_value()) {
        // If it claims feasible, verify manually
        Real sum = result->primal_values[0] + result->primal_values[1];
        EXPECT_GE(sum, 3.0 - 1e-8);
    }
}

TEST(Heuristics, SimpleRoundingRespectsBounds) {
    // x0 in [0,3], value=3.7 → rounds to 4 > ub → clamped to 3
    Problem p;
    p.set_dimensions(1, 1);
    std::vector<Index> r={0}, c={0};
    std::vector<Real>  v={1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 1, r, c, v));
    p.set_objective({1});
    p.set_row_bounds({-INF}, {INF});
    p.set_col_bounds({0}, {3});
    p.set_var_types({VarType::INTEGER});

    Solution lp_sol;
    lp_sol.primal_values = {3.7};
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(nullptr);
    auto result = h.simple_rounding(p, lp_sol);
    if (result.has_value())
        EXPECT_LE(result->primal_values[0], 3.0 + 1e-8);
}

// ── Fractional diving ─────────────────────────────────────────────────────────
TEST(Heuristics, FractionalDivingWithMockSolver) {
    Problem p = make_simple_mip();

    // Mock solver that returns integer-feasible solution after one dive
    auto mock = std::make_shared<MockLPSolver>();
    mock->fixed_sol.primal_values = {2.0, 1.0};  // integer after first fix
    mock->fixed_sol.status = SolverStatus::OPTIMAL;

    Solution lp_sol;
    lp_sol.primal_values = {1.7, 1.4};
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(mock.get());
    auto result = h.fractional_diving(p, lp_sol, 10);

    // Should find integer solution {2,1} since mock always returns it
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(result->primal_values[0], 2.0, 0.5);
}

TEST(Heuristics, FractionalDivingInfeasibleMock) {
    Problem p = make_simple_mip();

    auto mock = std::make_shared<MockLPSolver>();
    mock->fixed_status = SolverStatus::INFEASIBLE;  // always infeasible

    Solution lp_sol;
    lp_sol.primal_values = {1.7, 1.4};
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(mock.get());
    auto result = h.fractional_diving(p, lp_sol, 10);
    EXPECT_FALSE(result.has_value());
}

// ── Feasibility pump ──────────────────────────────────────────────────────────
TEST(Heuristics, FeasibilityPumpAlreadyInteger) {
    Problem p = make_simple_mip();

    // LP solution is already integer-feasible: {2,2}
    auto mock = std::make_shared<MockLPSolver>();
    mock->fixed_sol.primal_values = {2.0, 2.0};
    mock->fixed_sol.status = SolverStatus::OPTIMAL;

    Solution lp_sol;
    lp_sol.primal_values = {2.0, 2.0};  // already integer
    lp_sol.status = SolverStatus::OPTIMAL;

    Heuristics h(mock.get());
    auto result = h.feasibility_pump(p, lp_sol, 10);
    // Should immediately return the integer-feasible solution
    ASSERT_TRUE(result.has_value());
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
