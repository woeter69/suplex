#include "test_framework.h"
#include "src/core/problem.h"
#include "src/core/solution.h"
#include "src/simplex/primal_simplex.h"
#include "src/simplex/dual_simplex.h"
#include <cmath>

namespace suplex::test {

void test_primal_simplex_tiny() {
    // Problem:
    // minimize -3 x_1 - 2 x_2
    // subject to:
    //   x_1 + x_2 <= 4    =>  -INF <= x_1 + x_2 <= 4
    //   x_1 - x_2 <= 2    =>  -INF <= x_1 - x_2 <= 2
    //   0 <= x_1 <= INF
    //   0 <= x_2 <= INF
    // Expected optimal: x_1 = 3, x_2 = 1, obj = -11.0

    Problem prob(2, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -3.0);
    prob.set_objective_coeff(1, -2.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);

    prob.set_row_bounds(0, -INF, 4.0);
    prob.set_row_bounds(1, -INF, 2.0);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, -1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);

    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT(sol.is_optimal());
    TEST_ASSERT_NEAR(sol.objective_value, -11.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[0], 3.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[1], 1.0, 1e-6);
}

void test_dual_simplex_tiny() {
    Problem prob(2, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -3.0);
    prob.set_objective_coeff(1, -2.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);

    prob.set_row_bounds(0, -INF, 4.0);
    prob.set_row_bounds(1, -INF, 2.0);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, -1.0);
    prob.set_constraint_matrix(trip.to_csc());

    DualSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT(sol.is_optimal());
    TEST_ASSERT_NEAR(sol.objective_value, -11.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[0], 3.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[1], 1.0, 1e-6);
}

void test_unbounded() {
    // minimize -x_1 - x_2
    // subject to:
    //   x_1 - x_2 <= 1
    //   x_1, x_2 >= 0
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -1.0);
    prob.set_objective_coeff(1, -1.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);
    prob.set_row_bounds(0, -INF, 1.0);

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, -1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::UNBOUNDED);
}

void test_infeasible() {
    // minimize x_1
    // subject to:
    //   x_1 <= 2
    //   x_1 >= 4
    Problem prob(2, 1);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, 1.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_row_bounds(0, -INF, 2.0);
    prob.set_row_bounds(1, 4.0, INF);

    TripletMatrix trip(2, 1);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 0, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::INFEASIBLE);
}

void test_warm_start() {
    // Base problem:
    // minimize -3 x_1 - 2 x_2
    // s.t. x_1 + x_2 <= 4
    //      x_1 - x_2 <= 2
    // x_1, x_2 >= 0
    // Solved at x_1 = 3, x_2 = 1.
    //
    // Branch on x_1: add bound x_1 <= 2.0
    // Expected new solution: x_1 = 2, x_2 = 2, obj = -3*2 - 2*2 = -10.0

    Problem prob(2, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -3.0);
    prob.set_objective_coeff(1, -2.0);
    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);
    prob.set_row_bounds(0, -INF, 4.0);
    prob.set_row_bounds(1, -INF, 2.0);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, -1.0);
    prob.set_constraint_matrix(trip.to_csc());

    DualSimplexSolver solver;
    Solution sol;
    solver.solve(prob, sol);
    TEST_ASSERT(sol.is_optimal());
    std::vector<Index> parent_basis = solver.get_basis();

    // Now branch on x_1 <= 2.0
    prob.set_col_bounds(0, 0.0, 2.0);

    Solution child_sol;
    SolverStatus child_stat = solver.solve_from_basis(prob, parent_basis, child_sol);

    TEST_ASSERT(child_stat == SolverStatus::OPTIMAL);
    TEST_ASSERT(child_sol.is_optimal());
    TEST_ASSERT_NEAR(child_sol.objective_value, -10.0, 1e-6);
    TEST_ASSERT_NEAR(child_sol.primal_values[0], 2.0, 1e-6);
    TEST_ASSERT_NEAR(child_sol.primal_values[1], 2.0, 1e-6);
}

REGISTER_TEST(test_primal_simplex_tiny);
REGISTER_TEST(test_dual_simplex_tiny);
REGISTER_TEST(test_unbounded);
REGISTER_TEST(test_infeasible);
REGISTER_TEST(test_warm_start);

} // namespace suplex::test
