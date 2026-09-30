#pragma once

#include "types.h"
#include <vector>
#include <cassert>

namespace suplex {

struct SparseMatrixCSR;

/**
 * @brief Compressed Sparse Column (CSC) matrix format.
 * Optimal for column-oriented operations standard in revised simplex.
 */
struct SparseMatrixCSC {
    Index num_rows = 0;
    Index num_cols = 0;
    Index nnz = 0;
    std::vector<Index> col_start;
    std::vector<Index> row_index;
    std::vector<Real> values;

    SparseMatrixCSC() = default;
    SparseMatrixCSC(Index nrows, Index ncols, Index nonzeros = 0);

    Index col_begin(Index j) const {
        assert(j >= 0 && j < num_cols);
        return col_start[j];
    }

    Index col_end(Index j) const {
        assert(j >= 0 && j < num_cols);
        return col_start[j + 1];
    }

    Index col_length(Index j) const {
        assert(j >= 0 && j < num_cols);
        return col_start[j + 1] - col_start[j];
    }

    Real get(Index row, Index col) const;

    void multiply_vec(const Real* x, Real* y) const;
    void multiply_transpose_vec(const Real* x, Real* y) const;
    void multiply_col(Index j, Real scalar, Real* y) const;
    Real dot_col(Index j, const Real* x) const;

    static SparseMatrixCSC from_triplets(
        Index nrows,
        Index ncols,
        const std::vector<Index>& rows,
        const std::vector<Index>& cols,
        const std::vector<Real>& vals
    );

    SparseMatrixCSR to_csr() const;
    void clear();
};

/**
 * @brief Compressed Sparse Row (CSR) matrix format.
 * Useful for row-oriented operations, presolve, and GPU routines.
 */
struct SparseMatrixCSR {
    Index num_rows = 0;
    Index num_cols = 0;
    Index nnz = 0;
    std::vector<Index> row_start;
    std::vector<Index> col_index;
    std::vector<Real> values;

    SparseMatrixCSR() = default;
    SparseMatrixCSR(Index nrows, Index ncols, Index nonzeros = 0);

    Index row_begin(Index i) const {
        assert(i >= 0 && i < num_rows);
        return row_start[i];
    }

    Index row_end(Index i) const {
        assert(i >= 0 && i < num_rows);
        return row_start[i + 1];
    }

    Index row_length(Index i) const {
        assert(i >= 0 && i < num_rows);
        return row_start[i + 1] - row_start[i];
    }

    Real get(Index row, Index col) const;

    void multiply_vec(const Real* x, Real* y) const;
    SparseMatrixCSC to_csc() const;
    void clear();
};

/**
 * @brief Intermediate triplet (COO) matrix for incremental assembly.
 */
struct TripletMatrix {
    Index num_rows = 0;
    Index num_cols = 0;
    std::vector<Index> rows;
    std::vector<Index> cols;
    std::vector<Real> values;

    TripletMatrix() = default;
    TripletMatrix(Index nrows, Index ncols) : num_rows(nrows), num_cols(ncols) {}

    void add_entry(Index row, Index col, Real val);
    void reserve(Index nnz_est);
    void clear();
    SparseMatrixCSC to_csc() const;
};

} // namespace suplex
