#include "test_framework.h"
#include "src/simplex/lu_factor.h"
#include "src/simplex/lu_update.h"
#include <cmath>

namespace suplex::test {

void test_lu_small() {
    // 4x4 non-singular matrix
    // [ 2  1  0  0 ]
    // [ 1  3  1  0 ]
    // [ 0  2  4  1 ]
    // [ 0  0  1  5 ]
    TripletMatrix trip(4, 4);
    trip.add_entry(0, 0, 2.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, 3.0); trip.add_entry(1, 2, 1.0);
    trip.add_entry(2, 1, 2.0); trip.add_entry(2, 2, 4.0); trip.add_entry(2, 3, 1.0);
    trip.add_entry(3, 2, 1.0); trip.add_entry(3, 3, 5.0);

    SparseMatrixCSC B = trip.to_csc();

    SparseLU lu;
    SolverStatus stat = lu.factorize_matrix(B);
    TEST_ASSERT(stat == SolverStatus::OPTIMAL);
    TEST_ASSERT(lu.is_valid());

    // Expected x = [1, 2, 3, 4]
    // row 0: 2*1 + 1*2 = 4
    // row 1: 1*1 + 3*2 + 1*3 = 10
    // row 2: 2*2 + 4*3 + 1*4 = 20 (wait: 2*2 + 4*3 + 1*4 = 20, let's verify row 2: 2*x_1 + 4*x_2 + 1*x_3 = 2*2 + 4*3 + 1*4 = 20!)
    // row 3: 1*3 + 5*4 = 23
    Real b[4] = {4.0, 10.0, 20.0, 23.0};
    Real x[4] = {0.0, 0.0, 0.0, 0.0};

    lu.ftran(b, x);

    TEST_ASSERT_NEAR(x[0], 1.0, 1e-8);
    TEST_ASSERT_NEAR(x[1], 2.0, 1e-8);
    TEST_ASSERT_NEAR(x[2], 3.0, 1e-8);
    TEST_ASSERT_NEAR(x[3], 4.0, 1e-8);
}

void test_lu_btran() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 3.0); trip.add_entry(0, 1, 2.0);
    trip.add_entry(1, 0, 1.0); trip.add_entry(1, 1, 4.0); trip.add_entry(1, 2, 1.0);
    trip.add_entry(2, 1, 1.0); trip.add_entry(2, 2, 5.0);

    SparseMatrixCSC B = trip.to_csc();

    SparseLU lu;
    TEST_ASSERT(lu.factorize_matrix(B) == SolverStatus::OPTIMAL);

    Real c[3] = {7.0, 16.0, 17.0};
    // B^T y = c
    Real y[3] = {0.0, 0.0, 0.0};
    lu.btran(c, y);

    // Verify B^T * y == c
    Real b_t_y[3] = {0.0, 0.0, 0.0};
    B.multiply_transpose_vec(y, b_t_y);

    TEST_ASSERT_NEAR(b_t_y[0], c[0], 1e-8);
    TEST_ASSERT_NEAR(b_t_y[1], c[1], 1e-8);
    TEST_ASSERT_NEAR(b_t_y[2], c[2], 1e-8);
}

void test_lu_singular() {
    // 3x3 singular matrix (col 0 == col 1)
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 1.0); trip.add_entry(0, 1, 1.0); trip.add_entry(0, 2, 2.0);
    trip.add_entry(1, 0, 2.0); trip.add_entry(1, 1, 2.0); trip.add_entry(1, 2, 4.0);
    trip.add_entry(2, 0, 3.0); trip.add_entry(2, 1, 3.0); trip.add_entry(2, 2, 6.0);

    SparseMatrixCSC B = trip.to_csc();

    SparseLU lu;
    SolverStatus stat = lu.factorize_matrix(B);
    TEST_ASSERT(stat == SolverStatus::NUMERICAL_ERROR);
}

void test_lu_update() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 2.0); trip.add_entry(0, 1, 1.0);
    trip.add_entry(1, 1, 3.0);
    trip.add_entry(2, 0, 1.0); trip.add_entry(2, 2, 4.0);

    SparseMatrixCSC B = trip.to_csc();
    std::vector<Index> basis = {0, 1, 2};

    LUManager mgr;
    TEST_ASSERT(mgr.refactorize(B, basis) == SolverStatus::OPTIMAL);

    // New column entering at position 1: [2, 5, 0]^T
    // Solve alpha = B^{-1} * a_entering
    Real a_ent[3] = {2.0, 5.0, 0.0};
    std::vector<Real> alpha(3);
    mgr.ftran(a_ent, alpha.data());

    TEST_ASSERT(mgr.update(1, alpha) == SolverStatus::OPTIMAL);

    // New basis matrix B_new has column 1 replaced with [2, 5, 0]^T:
    // [ 2  2  0 ]
    // [ 0  5  0 ]
    // [ 1  0  4 ]
    // Solve B_new * x = b for b = [6, 15, 11]
    // x_1 = 3, 2*x_0 + 2*3 = 6 => x_0 = 0, 1*0 + 4*x_2 = 11 => wait: 4*x_2 = 12 => b_2 = 12
    Real b[3] = {6.0, 15.0, 12.0};
    Real x[3] = {0.0, 0.0, 0.0};
    mgr.ftran(b, x);

    // Expected x = [0, 3, 3]
    TEST_ASSERT_NEAR(x[0], 0.0, 1e-8);
    TEST_ASSERT_NEAR(x[1], 3.0, 1e-8);
    TEST_ASSERT_NEAR(x[2], 3.0, 1e-8);
}

REGISTER_TEST(test_lu_small);
REGISTER_TEST(test_lu_btran);
REGISTER_TEST(test_lu_singular);
REGISTER_TEST(test_lu_update);

} // namespace suplex::test
