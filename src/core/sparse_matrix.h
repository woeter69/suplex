#pragma once
#include "types.h"
#include <vector>
#include <optional>

namespace suplex {

// ── Forward declarations ──────────────────────────────────────────────────────
struct SparseMatrixCSC;
struct SparseMatrixCSR;

// ── Compressed Sparse Column ─────────────────────────────────────────────────
/// Column-major sparse storage. Column j spans row_index[col_start[j]..col_start[j+1]).
struct SparseMatrixCSC {
    Index num_rows = 0;
    Index num_cols = 0;
    Index nnz      = 0;
    std::vector<Index> col_start;  ///< size num_cols+1
    std::vector<Index> row_index;  ///< size nnz
    std::vector<Real>  values;     ///< size nnz

    Index col_begin(Index j) const { return col_start[j]; }
    Index col_end  (Index j) const { return col_start[j + 1]; }
    Index col_length(Index j) const { return col_start[j + 1] - col_start[j]; }

    // Implemented by Person 1 — stubs so Person 3 can compile
    void multiply_vec(const Real* x, Real* y) const;
    void multiply_transpose_vec(const Real* x, Real* y) const;
    void multiply_col(Index j, Real scalar, Real* y) const;
    Real dot_col(Index j, const Real* x) const;

    SparseMatrixCSR to_csr() const;

    static SparseMatrixCSC from_triplets(
        Index nrows, Index ncols,
        const std::vector<Index>& rows,
        const std::vector<Index>& cols,
        const std::vector<Real>&  vals);
};

// ── Compressed Sparse Row ────────────────────────────────────────────────────
/// Row-major sparse storage. Row i spans col_index[row_start[i]..row_start[i+1]).
struct SparseMatrixCSR {
    Index num_rows = 0;
    Index num_cols = 0;
    Index nnz      = 0;
    std::vector<Index> row_start;  ///< size num_rows+1
    std::vector<Index> col_index;  ///< size nnz
    std::vector<Real>  values;     ///< size nnz

    Index row_begin(Index i) const { return row_start[i]; }
    Index row_end  (Index i) const { return row_start[i + 1]; }
    Index row_length(Index i) const { return row_start[i + 1] - row_start[i]; }

    void multiply_vec(const Real* x, Real* y) const;
    SparseMatrixCSC to_csc() const;
};

// ── Triplet (COO) format ─────────────────────────────────────────────────────
/// Coordinate format for easy construction. Convert to CSC before solving.
struct TripletMatrix {
    Index num_rows = 0;
    Index num_cols = 0;
    std::vector<Index> rows;
    std::vector<Index> cols;
    std::vector<Real>  values;

    void add_entry(Index row, Index col, Real val) {
        rows.push_back(row);
        cols.push_back(col);
        values.push_back(val);
    }
    void reserve(Index nnz_estimate) {
        rows.reserve(nnz_estimate);
        cols.reserve(nnz_estimate);
        values.reserve(nnz_estimate);
    }
    SparseMatrixCSC to_csc() const;
};

} // namespace suplex
