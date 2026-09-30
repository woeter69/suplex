#include "test_framework.h"
#include "src/core/sparse_matrix.h"
#include "src/core/problem.h"
#include "src/core/solution.h"
#include "src/simplex/lu_factor.h"
#include "src/simplex/lu_update.h"
#include "src/simplex/primal_simplex.h"
#include "src/simplex/dual_simplex.h"
#include "src/simplex/basis.h"
#include "src/simplex/pricing.h"
#include <cmath>
#include <vector>
#include <numeric>

namespace suplex::test {

// ============================================================
// Sparse Matrix Boundary Conditions
// ============================================================

void test_empty_csc_from_triplets() {
    // Zero rows, zero cols, empty triplet vectors
    std::vector<Index> rows, cols;
    std::vector<Real> vals;
    SparseMatrixCSC csc = SparseMatrixCSC::from_triplets(0, 0, rows, cols, vals);
    TEST_ASSERT(csc.num_rows == 0);
    TEST_ASSERT(csc.num_cols == 0);
    TEST_ASSERT(csc.nnz == 0);
}

void test_csc_from_triplets_no_entries() {
    // Matrix has dimensions but no non-zeros
    std::vector<Index> rows, cols;
    std::vector<Real> vals;
    SparseMatrixCSC csc = SparseMatrixCSC::from_triplets(5, 5, rows, cols, vals);
    TEST_ASSERT(csc.num_rows == 5);
    TEST_ASSERT(csc.num_cols == 5);
    TEST_ASSERT(csc.nnz == 0);
    TEST_ASSERT_NEAR(csc.get(0, 0), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csc.get(4, 4), 0.0, 1e-15);
}

void test_csc_1x1_matrix() {
    std::vector<Index> rows = {0};
    std::vector<Index> cols = {0};
    std::vector<Real> vals = {42.0};
    SparseMatrixCSC csc = SparseMatrixCSC::from_triplets(1, 1, rows, cols, vals);
    TEST_ASSERT(csc.nnz == 1);
    TEST_ASSERT_NEAR(csc.get(0, 0), 42.0, 1e-15);

    Real x = 3.0;
    Real y = 0.0;
    csc.multiply_vec(&x, &y);
    TEST_ASSERT_NEAR(y, 126.0, 1e-12);
}

void test_csc_identity_matrix() {
    TripletMatrix trip(4, 4);
    for (Index i = 0; i < 4; ++i) {
        trip.add_entry(i, i, 1.0);
    }
    SparseMatrixCSC csc = trip.to_csc();
    TEST_ASSERT(csc.nnz == 4);
    for (Index i = 0; i < 4; ++i) {
        for (Index j = 0; j < 4; ++j) {
            Real expected = (i == j) ? 1.0 : 0.0;
            TEST_ASSERT_NEAR(csc.get(i, j), expected, 1e-15);
        }
    }
}

void test_csc_out_of_bounds_get() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 5.0);
    SparseMatrixCSC csc = trip.to_csc();
    // Out-of-range accesses should return 0.0 without crashing
    TEST_ASSERT_NEAR(csc.get(-1, 0), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csc.get(0, -1), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csc.get(2, 0), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csc.get(0, 2), 0.0, 1e-15);
}

void test_csc_spmv_zero_vector() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 2.0);
    trip.add_entry(1, 1, 3.0);
    trip.add_entry(2, 2, 4.0);
    SparseMatrixCSC csc = trip.to_csc();

    Real x[3] = {0.0, 0.0, 0.0};
    Real y[3] = {999.0, 999.0, 999.0};
    csc.multiply_vec(x, y);
    for (int i = 0; i < 3; ++i) {
        TEST_ASSERT_NEAR(y[i], 0.0, 1e-15);
    }
}

void test_csc_dot_col_out_of_range() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    SparseMatrixCSC csc = trip.to_csc();

    Real x[2] = {1.0, 1.0};
    TEST_ASSERT_NEAR(csc.dot_col(-1, x), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csc.dot_col(2, x), 0.0, 1e-15);
}

void test_csc_multiply_col_out_of_range() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    SparseMatrixCSC csc = trip.to_csc();

    Real y[2] = {0.0, 0.0};
    // Should not crash for out-of-range column
    csc.multiply_col(-1, 1.0, y);
    csc.multiply_col(2, 1.0, y);
    TEST_ASSERT_NEAR(y[0], 0.0, 1e-15);
    TEST_ASSERT_NEAR(y[1], 0.0, 1e-15);
}

void test_csc_multiply_col_zero_scalar() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 5.0);
    trip.add_entry(1, 0, 10.0);
    SparseMatrixCSC csc = trip.to_csc();

    Real y[2] = {1.0, 2.0};
    csc.multiply_col(0, 0.0, y);
    // Should remain unchanged because scalar is zero
    TEST_ASSERT_NEAR(y[0], 1.0, 1e-15);
    TEST_ASSERT_NEAR(y[1], 2.0, 1e-15);
}

void test_csc_multiple_duplicates_same_entry() {
    // Many duplicates at the same position should sum
    std::vector<Index> rows = {0, 0, 0, 0, 0};
    std::vector<Index> cols = {0, 0, 0, 0, 0};
    std::vector<Real> vals = {1.0, 2.0, 3.0, 4.0, 5.0};
    SparseMatrixCSC csc = SparseMatrixCSC::from_triplets(1, 1, rows, cols, vals);
    TEST_ASSERT(csc.nnz == 1);
    TEST_ASSERT_NEAR(csc.get(0, 0), 15.0, 1e-12);
}

void test_csr_empty_multiply() {
    SparseMatrixCSR csr(3, 3, 0);
    Real x[3] = {1.0, 2.0, 3.0};
    Real y[3] = {999.0, 999.0, 999.0};
    csr.multiply_vec(x, y);
    for (int i = 0; i < 3; ++i) {
        TEST_ASSERT_NEAR(y[i], 0.0, 1e-15);
    }
}

void test_csr_out_of_bounds_get() {
    SparseMatrixCSR csr(2, 2, 0);
    TEST_ASSERT_NEAR(csr.get(-1, 0), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csr.get(0, -1), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csr.get(2, 0), 0.0, 1e-15);
    TEST_ASSERT_NEAR(csr.get(0, 2), 0.0, 1e-15);
}

void test_triplet_clear() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 2.0);
    trip.clear();
    TEST_ASSERT(trip.num_rows == 0);
    TEST_ASSERT(trip.num_cols == 0);
    TEST_ASSERT(trip.rows.empty());
    TEST_ASSERT(trip.cols.empty());
    TEST_ASSERT(trip.values.empty());
}

void test_csc_clear() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 1.0);
    SparseMatrixCSC csc = trip.to_csc();
    csc.clear();
    TEST_ASSERT(csc.num_rows == 0);
    TEST_ASSERT(csc.num_cols == 0);
    TEST_ASSERT(csc.nnz == 0);
}

// ============================================================
// Problem Boundary Conditions
// ============================================================

void test_problem_zero_dimensions() {
    Problem prob(0, 0);
    TEST_ASSERT(prob.num_rows() == 0);
    TEST_ASSERT(prob.num_cols() == 0);
    TEST_ASSERT(prob.validate());
}

void test_problem_validate_mismatched() {
    Problem prob(2, 3);
    // Mess up the objective size
    std::vector<Real> bad_obj = {1.0, 2.0};
    prob.set_objective(std::move(bad_obj));
    TEST_ASSERT(!prob.validate());
}

void test_problem_validate_inverted_bounds() {
    Problem prob(1, 1);
    prob.set_col_bounds(0, 10.0, 5.0); // lower > upper
    TEST_ASSERT(!prob.validate());
}

void test_problem_set_bounds_out_of_range() {
    Problem prob(2, 2);
    // Out of range set should be silently ignored
    prob.set_col_bounds(-1, 0.0, 1.0);
    prob.set_col_bounds(2, 0.0, 1.0);
    prob.set_row_bounds(-1, 0.0, 1.0);
    prob.set_row_bounds(2, 0.0, 1.0);
    prob.set_objective_coeff(-1, 5.0);
    prob.set_objective_coeff(2, 5.0);
    // Should still validate cleanly
    TEST_ASSERT(prob.validate());
}

void test_problem_is_mip() {
    Problem prob(1, 3);
    TEST_ASSERT(!prob.is_mip());
    std::vector<VarType> types = {VarType::CONTINUOUS, VarType::INTEGER, VarType::BINARY};
    prob.set_var_types(types);
    TEST_ASSERT(prob.is_mip());
    TEST_ASSERT(prob.num_integers() == 2);
}

void test_problem_quadratic() {
    Problem prob(1, 2);
    TEST_ASSERT(!prob.is_qp());
    TEST_ASSERT(prob.quadratic_objective() == nullptr);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 1.0);
    prob.set_quadratic_objective(trip.to_csc());

    TEST_ASSERT(prob.is_qp());
    TEST_ASSERT(prob.quadratic_objective() != nullptr);
    TEST_ASSERT(prob.quadratic_objective()->num_rows == 2);
}

// ============================================================
// Solution Boundary Conditions
// ============================================================

void test_solution_reset() {
    Solution sol;
    sol.status = SolverStatus::OPTIMAL;
    sol.objective_value = 42.0;
    sol.primal_values = {1.0, 2.0};
    sol.num_iterations = 100;
    sol.reset();
    TEST_ASSERT(sol.status == SolverStatus::NOT_STARTED);
    TEST_ASSERT_NEAR(sol.objective_value, 0.0, 1e-15);
    TEST_ASSERT(sol.primal_values.empty());
    TEST_ASSERT(sol.num_iterations == 0);
}

void test_solution_is_feasible() {
    Solution sol;
    sol.status = SolverStatus::OPTIMAL;
    TEST_ASSERT(sol.is_feasible());

    sol.status = SolverStatus::ITERATION_LIMIT;
    TEST_ASSERT(sol.is_feasible());

    sol.status = SolverStatus::TIME_LIMIT;
    TEST_ASSERT(sol.is_feasible());

    sol.status = SolverStatus::INFEASIBLE;
    TEST_ASSERT(!sol.is_feasible());

    sol.status = SolverStatus::UNBOUNDED;
    TEST_ASSERT(!sol.is_feasible());

    sol.status = SolverStatus::NUMERICAL_ERROR;
    TEST_ASSERT(!sol.is_feasible());
}

// ============================================================
// LU Factorization Boundary Conditions
// ============================================================

void test_lu_1x1_factorize() {
    TripletMatrix trip(1, 1);
    trip.add_entry(0, 0, 7.0);
    SparseMatrixCSC B = trip.to_csc();

    SparseLU lu;
    TEST_ASSERT(lu.factorize_matrix(B) == SolverStatus::OPTIMAL);
    TEST_ASSERT(lu.is_valid());

    Real b = 21.0;
    Real x = 0.0;
    lu.ftran(&b, &x);
    TEST_ASSERT_NEAR(x, 3.0, 1e-12);

    Real c = 14.0;
    Real y = 0.0;
    lu.btran(&c, &y);
    TEST_ASSERT_NEAR(y, 2.0, 1e-12);
}

void test_lu_identity_factorize() {
    Index m = 5;
    TripletMatrix trip(m, m);
    for (Index i = 0; i < m; ++i) {
        trip.add_entry(i, i, 1.0);
    }
    SparseMatrixCSC B = trip.to_csc();

    SparseLU lu;
    TEST_ASSERT(lu.factorize_matrix(B) == SolverStatus::OPTIMAL);

    std::vector<Real> b = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<Real> x(m, 0.0);
    lu.ftran(b.data(), x.data());
    for (Index i = 0; i < m; ++i) {
        TEST_ASSERT_NEAR(x[i], b[i], 1e-12);
    }
}

void test_lu_near_singular() {
    // Matrix with a very small but non-zero entry
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 0.0);
    trip.add_entry(1, 0, 0.0);
    trip.add_entry(1, 1, 1e-15); // extremely small

    SparseMatrixCSC B = trip.to_csc();
    SparseLU lu;
    SolverStatus stat = lu.factorize_matrix(B);
    // Should detect near-singularity
    TEST_ASSERT(stat == SolverStatus::NUMERICAL_ERROR);
}

void test_lu_zero_dimension() {
    SparseMatrixCSC B;
    B.num_rows = 0;
    B.num_cols = 0;
    B.col_start = {};
    SparseLU lu;
    // 0-dim is handled specially in factorize_matrix: m_=0, num_cols != m_ => error
    // But factorize with empty basis should return OPTIMAL
    std::vector<Index> empty_basis;
    SolverStatus stat = lu.factorize(B, empty_basis);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
}

void test_lu_dimension_mismatch() {
    TripletMatrix trip(2, 3);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 1.0);
    SparseMatrixCSC B = trip.to_csc();
    SparseLU lu;
    SolverStatus stat = lu.factorize_matrix(B);
    TEST_ASSERT(stat == SolverStatus::NUMERICAL_ERROR);
}

void test_lu_btran_after_update() {
    // Verify BTRAN works correctly after an eta update
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 2.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 1, 3.0);
    trip.add_entry(2, 0, 1.0); trip.add_entry(2, 2, 4.0);

    SparseMatrixCSC B = trip.to_csc();
    std::vector<Index> basis = {0, 1, 2};
    LUManager mgr;
    TEST_ASSERT(mgr.refactorize(B, basis) == SolverStatus::OPTIMAL);

    // New column entering at position 1: [2, 5, 0]^T
    Real a_ent[3] = {2.0, 5.0, 0.0};
    std::vector<Real> alpha(3);
    mgr.ftran(a_ent, alpha.data());
    TEST_ASSERT(mgr.update(1, alpha) == SolverStatus::OPTIMAL);

    // Now the basis is:
    // [ 2  2  0 ]
    // [ 0  5  0 ]
    // [ 1  0  4 ]
    // Test BTRAN: solve B_new^T * y = c
    Real c[3] = {1.0, 0.0, 0.0};
    Real y[3] = {0.0, 0.0, 0.0};
    mgr.btran(c, y);

    // Verify: B_new^T * y should equal c
    // B_new^T = [ 2  0  1 ]
    //           [ 2  5  0 ]
    //           [ 0  0  4 ]
    // Check: 2*y[0] + 0*y[1] + 1*y[2] = 1
    //        2*y[0] + 5*y[1] + 0*y[2] = 0
    //        0*y[0] + 0*y[1] + 4*y[2] = 0
    // From row 3: y[2] = 0
    // From row 1: 2*y[0] = 1 => y[0] = 0.5
    // From row 2: 2*0.5 + 5*y[1] = 0 => y[1] = -0.2
    TEST_ASSERT_NEAR(y[0], 0.5, 1e-8);
    TEST_ASSERT_NEAR(y[1], -0.2, 1e-8);
    TEST_ASSERT_NEAR(y[2], 0.0, 1e-8);
}

void test_lu_update_near_zero_pivot() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 1.0);
    SparseMatrixCSC B = trip.to_csc();
    std::vector<Index> basis = {0, 1};

    LUManager mgr;
    mgr.refactorize(B, basis);

    // Update with alpha that has near-zero pivot at leaving position
    std::vector<Real> alpha = {1e-15, 1.0};
    SolverStatus stat = mgr.update(0, alpha);
    TEST_ASSERT(stat == SolverStatus::NUMERICAL_ERROR);
}

void test_lu_update_out_of_range() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 1.0);
    SparseMatrixCSC B = trip.to_csc();
    std::vector<Index> basis = {0, 1};

    LUManager mgr;
    mgr.refactorize(B, basis);

    std::vector<Real> alpha = {1.0, 1.0};
    TEST_ASSERT(mgr.update(-1, alpha) == SolverStatus::NUMERICAL_ERROR);
    TEST_ASSERT(mgr.update(2, alpha) == SolverStatus::NUMERICAL_ERROR);
}

void test_lu_should_refactorize() {
    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(1, 1, 1.0);
    SparseMatrixCSC B = trip.to_csc();
    std::vector<Index> basis = {0, 1};

    LUManager mgr;
    mgr.refactorize(B, basis);
    mgr.set_max_updates(3);

    TEST_ASSERT(!mgr.should_refactorize());

    std::vector<Real> alpha1 = {2.0, 0.0};
    mgr.update(0, alpha1);
    TEST_ASSERT(!mgr.should_refactorize());

    std::vector<Real> alpha2 = {0.0, 3.0};
    mgr.update(1, alpha2);
    TEST_ASSERT(!mgr.should_refactorize());

    std::vector<Real> alpha3 = {4.0, 0.0};
    mgr.update(0, alpha3);
    TEST_ASSERT(mgr.should_refactorize());
}

// ============================================================
// Basis Boundary Conditions
// ============================================================

void test_basis_pivot_out_of_range() {
    Basis basis(3, 6);
    basis.set_logical_basis(3);
    // Out-of-range pivot should be silently ignored
    basis.pivot(-1, 0, BasisStatus::AT_LOWER);
    basis.pivot(3, 0, BasisStatus::AT_LOWER);
    basis.pivot(0, -1, BasisStatus::AT_LOWER);
    basis.pivot(0, 6, BasisStatus::AT_LOWER);
    // Original basis should be intact
    TEST_ASSERT(basis.basic_var_at(0) == 3);
    TEST_ASSERT(basis.basic_var_at(1) == 4);
    TEST_ASSERT(basis.basic_var_at(2) == 5);
}

void test_basis_set_status_out_of_range() {
    Basis basis(2, 4);
    // Out-of-range should not crash
    basis.set_status(-1, BasisStatus::FIXED);
    basis.set_status(4, BasisStatus::FIXED);
}

// ============================================================
// Pricing Boundary Conditions
// ============================================================

void test_pricing_bland_rule() {
    Pricing pricing;
    pricing.init(5, 2, PricingStrategy::BLAND);

    std::vector<Real> rc = {0.0, -0.5, -0.1, 0.0, -0.3};
    std::vector<BasisStatus> status = {
        BasisStatus::BASIC, BasisStatus::AT_LOWER, BasisStatus::AT_LOWER,
        BasisStatus::BASIC, BasisStatus::AT_LOWER
    };
    Index entering = pricing.select_entering(rc, status);
    // Bland's rule: first eligible (index 1)
    TEST_ASSERT(entering == 1);
}

void test_pricing_no_eligible() {
    Pricing pricing;
    pricing.init(4, 2, PricingStrategy::DANTZIG);

    // All reduced costs are non-negative for AT_LOWER
    std::vector<Real> rc = {0.0, 0.5, 0.1, 0.0};
    std::vector<BasisStatus> status = {
        BasisStatus::BASIC, BasisStatus::AT_LOWER, BasisStatus::AT_LOWER,
        BasisStatus::BASIC
    };
    Index entering = pricing.select_entering(rc, status);
    TEST_ASSERT(entering == -1);
}

void test_pricing_free_variable() {
    Pricing pricing;
    pricing.init(3, 1, PricingStrategy::DANTZIG);

    std::vector<Real> rc = {0.0, 0.5, 0.0};
    std::vector<BasisStatus> status = {
        BasisStatus::BASIC, BasisStatus::FREE_ZERO, BasisStatus::FIXED
    };
    // Free variable with rc > EPS_OPTIMALITY should be eligible
    Index entering = pricing.select_entering(rc, status);
    TEST_ASSERT(entering == 1);
}

void test_pricing_fixed_variable_skipped() {
    Pricing pricing;
    pricing.init(3, 1, PricingStrategy::DANTZIG);

    std::vector<Real> rc = {0.0, -1.0, 0.0};
    std::vector<BasisStatus> status = {
        BasisStatus::BASIC, BasisStatus::FIXED, BasisStatus::BASIC
    };
    // Fixed variables should never be selected
    Index entering = pricing.select_entering(rc, status);
    TEST_ASSERT(entering == -1);
}

// ============================================================
// Simplex Solver Boundary Conditions
// ============================================================

void test_simplex_maximization() {
    // maximize 2 x_1 + 3 x_2
    // subject to x_1 + x_2 <= 4
    //            x_1 + 3 x_2 <= 6
    //            x_1, x_2 >= 0
    // Optimal: x_1 = 3, x_2 = 1, obj = 2*3 + 3*1 = 9
    Problem prob(2, 2);
    prob.set_objective_sense(ObjectiveSense::MAXIMIZE);
    prob.set_objective_coeff(0, 2.0);
    prob.set_objective_coeff(1, 3.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);

    prob.set_row_bounds(0, -INF, 4.0);
    prob.set_row_bounds(1, -INF, 6.0);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, 3.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.objective_value, 9.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[0], 3.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[1], 1.0, 1e-6);
}

void test_simplex_equality_constraints() {
    // minimize x_1 + x_2
    // subject to x_1 + x_2 = 5  (row_lower = row_upper = 5)
    //            x_1, x_2 >= 0
    // Optimal: any x_1 + x_2 = 5 with obj = 5
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, 1.0);
    prob.set_objective_coeff(1, 1.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);

    prob.set_row_bounds(0, 5.0, 5.0); // equality

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.objective_value, 5.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[0] + sol.primal_values[1], 5.0, 1e-6);
}

void test_simplex_fixed_variable() {
    // minimize x_1 + x_2
    // subject to x_1 + x_2 <= 10
    //            x_1 = 3 (fixed)
    //            x_2 >= 0
    // Optimal: x_1 = 3, x_2 = 0, obj = 3
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, 1.0);
    prob.set_objective_coeff(1, 1.0);

    prob.set_col_bounds(0, 3.0, 3.0); // fixed
    prob.set_col_bounds(1, 0.0, INF);

    prob.set_row_bounds(0, -INF, 10.0);

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.primal_values[0], 3.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[1], 0.0, 1e-6);
    TEST_ASSERT_NEAR(sol.objective_value, 3.0, 1e-6);
}

void test_simplex_single_variable() {
    // minimize 5 * x_1
    // subject to x_1 >= 2
    // x_1 unbounded above
    // Optimal: x_1 = 2, obj = 10
    Problem prob(1, 1);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, 5.0);
    prob.set_col_bounds(0, 0.0, INF);
    prob.set_row_bounds(0, 2.0, INF);

    TripletMatrix trip(1, 1);
    trip.add_entry(0, 0, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.primal_values[0], 2.0, 1e-6);
    TEST_ASSERT_NEAR(sol.objective_value, 10.0, 1e-6);
}

void test_simplex_upper_bounded_variables() {
    // minimize -x_1 - x_2
    // subject to x_1 + x_2 <= 100 (slack)
    //            0 <= x_1 <= 3
    //            0 <= x_2 <= 4
    // Optimal: x_1 = 3, x_2 = 4, obj = -7
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -1.0);
    prob.set_objective_coeff(1, -1.0);

    prob.set_col_bounds(0, 0.0, 3.0);
    prob.set_col_bounds(1, 0.0, 4.0);

    prob.set_row_bounds(0, -INF, 100.0);

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.primal_values[0], 3.0, 1e-6);
    TEST_ASSERT_NEAR(sol.primal_values[1], 4.0, 1e-6);
    TEST_ASSERT_NEAR(sol.objective_value, -7.0, 1e-6);
}

void test_simplex_dual_infeasible_detected() {
    // Same infeasible problem but via dual simplex
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

    DualSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::INFEASIBLE);
}

void test_simplex_iteration_limit() {
    // A valid LP but force iteration limit to 0
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -1.0);
    prob.set_objective_coeff(1, -1.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);
    prob.set_row_bounds(0, -INF, 10.0);

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    solver.set_iteration_limit(0);
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::ITERATION_LIMIT);
}

void test_simplex_larger_problem() {
    // 5-variable, 3-constraint LP
    // minimize -2x1 - 4x2 - 3x3 - x4 - x5
    // subject to:
    //   3x1 + x2 + x3 + x4       <= 10
    //   x1 + 2x2 + x3      + x5  <= 8
    //   x1 + x2 + 3x3             <= 6
    //   all vars >= 0
    Problem prob(3, 5);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -2.0);
    prob.set_objective_coeff(1, -4.0);
    prob.set_objective_coeff(2, -3.0);
    prob.set_objective_coeff(3, -1.0);
    prob.set_objective_coeff(4, -1.0);

    for (Index j = 0; j < 5; ++j) {
        prob.set_col_bounds(j, 0.0, INF);
    }
    prob.set_row_bounds(0, -INF, 10.0);
    prob.set_row_bounds(1, -INF, 8.0);
    prob.set_row_bounds(2, -INF, 6.0);

    TripletMatrix trip(3, 5);
    trip.add_entry(0, 0, 3.0); trip.add_entry(0, 1, 1.0); trip.add_entry(0, 2, 1.0); trip.add_entry(0, 3, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, 2.0); trip.add_entry(1, 2, 1.0); trip.add_entry(1, 4, 1.0);
    trip.add_entry(2, 0, 1.0); trip.add_entry(2, 1, 1.0); trip.add_entry(2, 2, 3.0);
    prob.set_constraint_matrix(trip.to_csc());

    {
        PrimalSimplexSolver solver;
        Solution sol;
        SolverStatus stat = solver.solve(prob, sol);
        TEST_ASSERT(stat == SolverStatus::OPTIMAL);
        // Verify feasibility: Ax <= b
        for (Index i = 0; i < 3; ++i) {
            TEST_ASSERT(sol.row_activities[i] <= prob.row_upper()[i] + 1e-6);
        }
        // All primal values non-negative
        for (Index j = 0; j < 5; ++j) {
            TEST_ASSERT(sol.primal_values[j] >= -1e-6);
        }
    }

    {
        DualSimplexSolver solver;
        Solution sol;
        SolverStatus stat = solver.solve(prob, sol);
        TEST_ASSERT(stat == SolverStatus::OPTIMAL);
        for (Index i = 0; i < 3; ++i) {
            TEST_ASSERT(sol.row_activities[i] <= prob.row_upper()[i] + 1e-6);
        }
        for (Index j = 0; j < 5; ++j) {
            TEST_ASSERT(sol.primal_values[j] >= -1e-6);
        }
    }
}

void test_simplex_primal_dual_agree() {
    // Solve the same problem with both solvers and verify they agree
    Problem prob(2, 3);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -5.0);
    prob.set_objective_coeff(1, -4.0);
    prob.set_objective_coeff(2, -3.0);

    for (Index j = 0; j < 3; ++j) {
        prob.set_col_bounds(j, 0.0, INF);
    }
    prob.set_row_bounds(0, -INF, 10.0);
    prob.set_row_bounds(1, -INF, 8.0);

    TripletMatrix trip(2, 3);
    trip.add_entry(0, 0, 6.0); trip.add_entry(0, 1, 4.0); trip.add_entry(0, 2, 2.0);
    trip.add_entry(1, 0, 3.0); trip.add_entry(1, 1, 5.0); trip.add_entry(1, 2, 4.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver primal;
    Solution psol;
    SolverStatus pstat = primal.solve(prob, psol);

    DualSimplexSolver dual;
    Solution dsol;
    SolverStatus dstat = dual.solve(prob, dsol);

    TEST_ASSERT(pstat == SolverStatus::OPTIMAL);
    TEST_ASSERT(dstat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(psol.objective_value, dsol.objective_value, 1e-6);
}

void test_simplex_zero_objective() {
    // minimize 0 (feasibility only)
    // subject to x_1 + x_2 <= 5
    //            x_1, x_2 >= 0
    Problem prob(1, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, 0.0);
    prob.set_objective_coeff(1, 0.0);

    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);
    prob.set_row_bounds(0, -INF, 5.0);

    TripletMatrix trip(1, 2);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 1.0);
    prob.set_constraint_matrix(trip.to_csc());

    PrimalSimplexSolver solver;
    Solution sol;
    SolverStatus stat = solver.solve(prob, sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT_NEAR(sol.objective_value, 0.0, 1e-6);
}

void test_dual_simplex_warm_start_tighten_both() {
    // Warm start: tighten both bounds on a variable
    Problem prob(2, 2);
    prob.set_objective_sense(ObjectiveSense::MINIMIZE);
    prob.set_objective_coeff(0, -3.0);
    prob.set_objective_coeff(1, -2.0);
    prob.set_col_bounds(0, 0.0, INF);
    prob.set_col_bounds(1, 0.0, INF);
    prob.set_row_bounds(0, -INF, 10.0);
    prob.set_row_bounds(1, -INF, 8.0);

    TripletMatrix trip(2, 2);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, 2.0);
    prob.set_constraint_matrix(trip.to_csc());

    DualSimplexSolver solver;
    Solution sol;
    solver.solve(prob, sol);
    TEST_ASSERT(sol.is_optimal());
    std::vector<Index> parent_basis = solver.get_basis();

    // Tighten bounds
    prob.set_col_bounds(0, 1.0, 5.0);
    prob.set_col_bounds(1, 0.0, 3.0);

    Solution child_sol;
    SolverStatus stat = solver.solve_from_basis(prob, parent_basis, child_sol);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    // Verify feasibility
    TEST_ASSERT(child_sol.primal_values[0] >= 1.0 - 1e-6);
    TEST_ASSERT(child_sol.primal_values[0] <= 5.0 + 1e-6);
    TEST_ASSERT(child_sol.primal_values[1] >= -1e-6);
    TEST_ASSERT(child_sol.primal_values[1] <= 3.0 + 1e-6);
}

// ============================================================
// to_string coverage
// ============================================================

void test_to_string_coverage() {
    // SolverStatus
    TEST_ASSERT(to_string(SolverStatus::OPTIMAL) == std::string("OPTIMAL"));
    TEST_ASSERT(to_string(SolverStatus::INFEASIBLE) == std::string("INFEASIBLE"));
    TEST_ASSERT(to_string(SolverStatus::UNBOUNDED) == std::string("UNBOUNDED"));
    TEST_ASSERT(to_string(SolverStatus::INF_OR_UNBD) == std::string("INF_OR_UNBD"));
    TEST_ASSERT(to_string(SolverStatus::NOT_STARTED) == std::string("NOT_STARTED"));
    TEST_ASSERT(to_string(SolverStatus::ITERATION_LIMIT) == std::string("ITERATION_LIMIT"));
    TEST_ASSERT(to_string(SolverStatus::TIME_LIMIT) == std::string("TIME_LIMIT"));
    TEST_ASSERT(to_string(SolverStatus::NUMERICAL_ERROR) == std::string("NUMERICAL_ERROR"));
    TEST_ASSERT(to_string(SolverStatus::USER_INTERRUPT) == std::string("USER_INTERRUPT"));

    // BasisStatus
    TEST_ASSERT(std::string(to_string(BasisStatus::BASIC)) == "BASIC");
    TEST_ASSERT(std::string(to_string(BasisStatus::AT_LOWER)) == "AT_LOWER");
    TEST_ASSERT(std::string(to_string(BasisStatus::AT_UPPER)) == "AT_UPPER");
    TEST_ASSERT(std::string(to_string(BasisStatus::FIXED)) == "FIXED");
    TEST_ASSERT(std::string(to_string(BasisStatus::FREE_ZERO)) == "FREE_ZERO");

    // PricingStrategy
    TEST_ASSERT(std::string(to_string(PricingStrategy::DANTZIG)) == "DANTZIG");
    TEST_ASSERT(std::string(to_string(PricingStrategy::DEVEX)) == "DEVEX");
    TEST_ASSERT(std::string(to_string(PricingStrategy::STEEPEST_EDGE)) == "STEEPEST_EDGE");
    TEST_ASSERT(std::string(to_string(PricingStrategy::BLAND)) == "BLAND");
}

// ============================================================
// Registration
// ============================================================

REGISTER_TEST(test_empty_csc_from_triplets);
REGISTER_TEST(test_csc_from_triplets_no_entries);
REGISTER_TEST(test_csc_1x1_matrix);
REGISTER_TEST(test_csc_identity_matrix);
REGISTER_TEST(test_csc_out_of_bounds_get);
REGISTER_TEST(test_csc_spmv_zero_vector);
REGISTER_TEST(test_csc_dot_col_out_of_range);
REGISTER_TEST(test_csc_multiply_col_out_of_range);
REGISTER_TEST(test_csc_multiply_col_zero_scalar);
REGISTER_TEST(test_csc_multiple_duplicates_same_entry);
REGISTER_TEST(test_csr_empty_multiply);
REGISTER_TEST(test_csr_out_of_bounds_get);
REGISTER_TEST(test_triplet_clear);
REGISTER_TEST(test_csc_clear);
REGISTER_TEST(test_problem_zero_dimensions);
REGISTER_TEST(test_problem_validate_mismatched);
REGISTER_TEST(test_problem_validate_inverted_bounds);
REGISTER_TEST(test_problem_set_bounds_out_of_range);
REGISTER_TEST(test_problem_is_mip);
REGISTER_TEST(test_problem_quadratic);
REGISTER_TEST(test_solution_reset);
REGISTER_TEST(test_solution_is_feasible);
REGISTER_TEST(test_lu_1x1_factorize);
REGISTER_TEST(test_lu_identity_factorize);
REGISTER_TEST(test_lu_near_singular);
REGISTER_TEST(test_lu_zero_dimension);
REGISTER_TEST(test_lu_dimension_mismatch);
REGISTER_TEST(test_lu_btran_after_update);
REGISTER_TEST(test_lu_update_near_zero_pivot);
REGISTER_TEST(test_lu_update_out_of_range);
REGISTER_TEST(test_lu_should_refactorize);
REGISTER_TEST(test_basis_pivot_out_of_range);
REGISTER_TEST(test_basis_set_status_out_of_range);
REGISTER_TEST(test_pricing_bland_rule);
REGISTER_TEST(test_pricing_no_eligible);
REGISTER_TEST(test_pricing_free_variable);
REGISTER_TEST(test_pricing_fixed_variable_skipped);
REGISTER_TEST(test_simplex_maximization);
REGISTER_TEST(test_simplex_equality_constraints);
REGISTER_TEST(test_simplex_fixed_variable);
REGISTER_TEST(test_simplex_single_variable);
REGISTER_TEST(test_simplex_upper_bounded_variables);
REGISTER_TEST(test_simplex_dual_infeasible_detected);
REGISTER_TEST(test_simplex_iteration_limit);
REGISTER_TEST(test_simplex_larger_problem);
REGISTER_TEST(test_simplex_primal_dual_agree);
REGISTER_TEST(test_simplex_zero_objective);
REGISTER_TEST(test_dual_simplex_warm_start_tighten_both);
REGISTER_TEST(test_to_string_coverage);

} // namespace suplex::test
