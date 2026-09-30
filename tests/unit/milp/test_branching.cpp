/// Unit tests for the Branching module.
#include <gtest/gtest.h>
#include "milp/branching.h"
#include "core/problem.h"
#include "core/solution.h"
#include "core/types.h"

using namespace suplex;

// ── Helpers ───────────────────────────────────────────────────────────────────
static Solution make_lp_sol(const std::vector<Real>& vals) {
    Solution s;
    s.status = SolverStatus::OPTIMAL;
    s.primal_values = vals;
    s.objective_value = 0.0;
    return s;
}

// ── PseudoCosts ───────────────────────────────────────────────────────────────
TEST(PseudoCosts, UpdateDown) {
    PseudoCosts pc(3);
    pc.update_down(0, 2.6, 5.0);   // frac = 0.6
    EXPECT_EQ(pc.down_count[0], 1);
    // delta / frac = 5.0 / 0.6 = 8.33...
    EXPECT_NEAR(pc.down_cost[0], 8.333, 0.01);
}

TEST(PseudoCosts, UpdateUp) {
    PseudoCosts pc(3);
    pc.update_up(1, 2.3, 7.0);    // frac_up = ceil - val = 0.7
    EXPECT_EQ(pc.up_count[1], 1);
    EXPECT_NEAR(pc.up_cost[1], 10.0, 0.01);  // 7.0 / 0.7
}

TEST(PseudoCosts, ReliabilityThreshold) {
    PseudoCosts pc(2);
    EXPECT_FALSE(pc.is_reliable(0, 8));

    for (int i = 0; i < 8; ++i) {
        pc.update_down(0, 1.5, 1.0);
        pc.update_up(0, 1.5, 1.0);
    }
    EXPECT_TRUE(pc.is_reliable(0, 8));
}

TEST(PseudoCosts, ScorePositive) {
    PseudoCosts pc(2);
    pc.update_down(0, 1.5, 2.0);  // down_cost = 4.0
    pc.update_up  (0, 1.5, 3.0);  // up_cost   = 6.0
    Real score = pc.score(0, 1.5);
    EXPECT_GT(score, 0.0);
}

// ── Branching variable selection ──────────────────────────────────────────────
TEST(Branching, MostFractionalSelected) {
    // Create a simple problem with 3 integer variables
    Problem p;
    p.set_dimensions(1, 3);
    std::vector<Index> r = {0,0,0}, c = {0,1,2};
    std::vector<Real>  v = {1,1,1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 3, r, c, v));
    p.set_objective({1,1,1});
    p.set_row_bounds({-INF}, {10});
    p.set_col_bounds({0,0,0}, {10,10,10});
    p.set_var_types({VarType::INTEGER, VarType::INTEGER, VarType::INTEGER});

    // x0=1.1 (frac=0.1), x1=2.5 (frac=0.5 → most fractional), x2=3.8 (frac=0.2)
    Solution sol = make_lp_sol({1.1, 2.5, 3.8});

    Branching br(3);
    Index selected = br.select_variable(p, sol, nullptr, 0.0, 0);
    EXPECT_EQ(selected, 1);  // x1 is most fractional
}

TEST(Branching, IntegerVariableSkipped) {
    Problem p;
    p.set_dimensions(1, 2);
    std::vector<Index> r = {0,0}, c = {0,1};
    std::vector<Real>  v = {1,1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 2, r, c, v));
    p.set_objective({1,1});
    p.set_row_bounds({-INF}, {10});
    p.set_col_bounds({0,0}, {10,10});
    p.set_var_types({VarType::INTEGER, VarType::INTEGER});

    // Both already integer
    Solution sol = make_lp_sol({3.0, 4.0});
    Branching br(2);
    Index selected = br.select_variable(p, sol, nullptr, 0.0, 0);
    EXPECT_EQ(selected, -1);  // no fractional variables
}

TEST(Branching, ContinuousSkipped) {
    Problem p;
    p.set_dimensions(1, 2);
    std::vector<Index> r = {0,0}, c = {0,1};
    std::vector<Real>  v = {1,1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 2, r, c, v));
    p.set_objective({1,1});
    p.set_row_bounds({-INF}, {10});
    p.set_col_bounds({0,0}, {10,10});
    p.set_var_types({VarType::CONTINUOUS, VarType::INTEGER});

    // x0 continuous (fractional ok), x1 integer and fractional → select x1
    Solution sol = make_lp_sol({1.7, 2.4});
    Branching br(2);
    Index selected = br.select_variable(p, sol, nullptr, 0.0, 0);
    EXPECT_EQ(selected, 1);
}

TEST(Branching, PseudoCostUpdateRecorded) {
    Branching br(3);
    // After observing pseudo-costs, reliability should increase
    for (int i = 0; i < 10; ++i)
        br.record_observation(0, BranchDirection::DOWN, 2.4, 100.0, 105.0);
    for (int i = 0; i < 10; ++i)
        br.record_observation(0, BranchDirection::UP, 2.4, 100.0, 103.0);

    // Now variable 0 should be considered reliable
    PseudoCosts pc(3);
    for (int i = 0; i < 10; ++i) {
        pc.update_down(0, 2.4, 5.0);
        pc.update_up  (0, 2.4, 3.0);
    }
    EXPECT_TRUE(pc.is_reliable(0, 8));
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
