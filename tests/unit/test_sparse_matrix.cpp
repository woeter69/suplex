#include "test_framework.h"
#include "src/core/sparse_matrix.h"
#include <vector>
#include <cmath>

namespace suplex::test {

void test_csc_from_triplets() {
    // 3x3 matrix with duplicate at (0, 0)
    std::vector<Index> rows = {0, 1, 2, 0, 1};
    std::vector<Index> cols = {0, 0, 1, 0, 2};
    std::vector<Real> vals  = {2.0, 3.0, 4.0, 5.0, 6.0};

    SparseMatrixCSC csc = SparseMatrixCSC::from_triplets(3, 3, rows, cols, vals);

    TEST_ASSERT(csc.num_rows == 3);
    TEST_ASSERT(csc.num_cols == 3);
    TEST_ASSERT(csc.nnz == 4);

    // Entry at (0, 0) should be 2.0 + 5.0 = 7.0
    TEST_ASSERT_NEAR(csc.get(0, 0), 7.0, 1e-12);
    TEST_ASSERT_NEAR(csc.get(1, 0), 3.0, 1e-12);
    TEST_ASSERT_NEAR(csc.get(2, 1), 4.0, 1e-12);
    TEST_ASSERT_NEAR(csc.get(1, 2), 6.0, 1e-12);
    TEST_ASSERT_NEAR(csc.get(0, 1), 0.0, 1e-12);
}

void test_csc_spmv() {
    TripletMatrix trip(3, 4);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 2, 2.0);
    trip.add_entry(1, 1, 3.0);
    trip.add_entry(1, 3, 4.0);
    trip.add_entry(2, 0, 5.0);
    trip.add_entry(2, 2, 6.0);

    SparseMatrixCSC csc = trip.to_csc();

    Real x[4] = {1.0, 2.0, 3.0, 4.0};
    Real y[3] = {0.0, 0.0, 0.0};

    csc.multiply_vec(x, y);

    // Row 0: 1*1 + 2*3 = 7
    // Row 1: 3*2 + 4*4 = 22
    // Row 2: 5*1 + 6*3 = 23
    TEST_ASSERT_NEAR(y[0], 7.0, 1e-12);
    TEST_ASSERT_NEAR(y[1], 22.0, 1e-12);
    TEST_ASSERT_NEAR(y[2], 23.0, 1e-12);
}

void test_csc_transpose_spmv() {
    TripletMatrix trip(3, 3);
    trip.add_entry(0, 0, 1.0);
    trip.add_entry(0, 1, 2.0);
    trip.add_entry(1, 0, 3.0);
    trip.add_entry(2, 2, 4.0);

    SparseMatrixCSC csc = trip.to_csc();

    Real x[3] = {2.0, 1.0, 3.0};
    Real y[3] = {0.0, 0.0, 0.0};

    csc.multiply_transpose_vec(x, y);

    // A^T has columns as original rows:
    // col 0: A[0,0]*2 + A[1,0]*1 = 1*2 + 3*1 = 5
    // col 1: A[0,1]*2 = 2*2 = 4
    // col 2: A[2,2]*3 = 4*3 = 12
    TEST_ASSERT_NEAR(y[0], 5.0, 1e-12);
    TEST_ASSERT_NEAR(y[1], 4.0, 1e-12);
    TEST_ASSERT_NEAR(y[2], 12.0, 1e-12);
}

void test_csc_to_csr_roundtrip() {
    TripletMatrix trip(4, 4);
    trip.add_entry(0, 1, 10.0);
    trip.add_entry(1, 0, 20.0);
    trip.add_entry(2, 2, 30.0);
    trip.add_entry(3, 3, 40.0);
    trip.add_entry(0, 3, 50.0);

    SparseMatrixCSC csc1 = trip.to_csc();
    SparseMatrixCSR csr = csc1.to_csr();
    SparseMatrixCSC csc2 = csr.to_csc();

    TEST_ASSERT(csc1.nnz == csc2.nnz);
    TEST_ASSERT(csc1.num_rows == csc2.num_rows);
    TEST_ASSERT(csc1.num_cols == csc2.num_cols);

    for (Index i = 0; i < 4; ++i) {
        for (Index j = 0; j < 4; ++j) {
            TEST_ASSERT_NEAR(csc1.get(i, j), csc2.get(i, j), 1e-12);
            TEST_ASSERT_NEAR(csc1.get(i, j), csr.get(i, j), 1e-12);
        }
    }
}

REGISTER_TEST(test_csc_from_triplets);
REGISTER_TEST(test_csc_spmv);
REGISTER_TEST(test_csc_transpose_spmv);
REGISTER_TEST(test_csc_to_csr_roundtrip);

} // namespace suplex::test
