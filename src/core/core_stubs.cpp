// Stub implementations for Person 1's sparse matrix operations.
// Replace with Person 1's actual src/core/sparse_matrix.cpp when integrating.
#include "sparse_matrix.h"
#include "problem.h"
#include <algorithm>
#include <numeric>
#include <cassert>

namespace suplex {

// ── Problem implementation stubs ──────────────────────────────────────────────

void Problem::set_dimensions(Index rows, Index cols) {
    A_.num_rows = rows;
    A_.num_cols = cols;
    c_.assign(cols, 0.0);
    col_lower_.assign(cols, 0.0);
    col_upper_.assign(cols, INF);
    row_lower_.assign(rows, -INF);
    row_upper_.assign(rows,  INF);
    var_types_.assign(cols, VarType::CONTINUOUS);
}

void Problem::set_col_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper) {
    col_lower_ = std::move(lower);
    col_upper_ = std::move(upper);
}

void Problem::set_row_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper) {
    row_lower_ = std::move(lower);
    row_upper_ = std::move(upper);
}

Index Problem::num_integers() const {
    Index n = 0;
    for (auto t : var_types_)
        if (t == VarType::INTEGER || t == VarType::BINARY) ++n;
    return n;
}

bool Problem::is_mip() const { return num_integers() > 0; }

bool Problem::validate() const {
    if (A_.num_rows != static_cast<Index>(row_lower_.size())) return false;
    if (A_.num_rows != static_cast<Index>(row_upper_.size())) return false;
    if (A_.num_cols != static_cast<Index>(col_lower_.size())) return false;
    if (A_.num_cols != static_cast<Index>(col_upper_.size())) return false;
    if (A_.num_cols != static_cast<Index>(c_.size()))         return false;
    if (A_.num_cols != static_cast<Index>(var_types_.size())) return false;
    return true;
}

// ── Sparse matrix stub ops (to be replaced by Person 1) ──────────────────────

void SparseMatrixCSC::multiply_vec(const Real* x, Real* y) const {
    std::fill(y, y + num_rows, 0.0);
    for (Index j = 0; j < num_cols; ++j) {
        if (x[j] == 0.0) continue;
        for (Index k = col_start[j]; k < col_start[j + 1]; ++k)
            y[row_index[k]] += values[k] * x[j];
    }
}

void SparseMatrixCSC::multiply_transpose_vec(const Real* x, Real* y) const {
    for (Index j = 0; j < num_cols; ++j) {
        Real s = 0.0;
        for (Index k = col_start[j]; k < col_start[j + 1]; ++k)
            s += values[k] * x[row_index[k]];
        y[j] = s;
    }
}

void SparseMatrixCSC::multiply_col(Index j, Real scalar, Real* y) const {
    for (Index k = col_start[j]; k < col_start[j + 1]; ++k)
        y[row_index[k]] += scalar * values[k];
}

Real SparseMatrixCSC::dot_col(Index j, const Real* x) const {
    Real s = 0.0;
    for (Index k = col_start[j]; k < col_start[j + 1]; ++k)
        s += values[k] * x[row_index[k]];
    return s;
}

SparseMatrixCSR SparseMatrixCSC::to_csr() const {
    SparseMatrixCSR csr;
    csr.num_rows = num_rows;
    csr.num_cols = num_cols;
    csr.nnz      = nnz;
    csr.row_start.assign(num_rows + 1, 0);
    csr.col_index.resize(nnz);
    csr.values.resize(nnz);

    for (Index k = 0; k < nnz; ++k)
        csr.row_start[row_index[k] + 1]++;

    for (Index i = 0; i < num_rows; ++i)
        csr.row_start[i + 1] += csr.row_start[i];

    std::vector<Index> cursor(csr.row_start.begin(), csr.row_start.begin() + num_rows);
    for (Index j = 0; j < num_cols; ++j)
        for (Index k = col_start[j]; k < col_start[j + 1]; ++k) {
            Index i   = row_index[k];
            Index pos = cursor[i]++;
            csr.col_index[pos] = j;
            csr.values[pos]    = values[k];
        }
    return csr;
}

SparseMatrixCSC SparseMatrixCSC::from_triplets(
        Index nrows, Index ncols,
        const std::vector<Index>& rows,
        const std::vector<Index>& cols,
        const std::vector<Real>&  vals) {

    SparseMatrixCSC csc;
    csc.num_rows = nrows;
    csc.num_cols = ncols;
    csc.col_start.assign(ncols + 1, 0);

    Index ne = static_cast<Index>(vals.size());
    for (Index k = 0; k < ne; ++k)
        csc.col_start[cols[k] + 1]++;
    for (Index j = 0; j < ncols; ++j)
        csc.col_start[j + 1] += csc.col_start[j];
    csc.nnz = csc.col_start[ncols];
    csc.row_index.resize(csc.nnz);
    csc.values.resize(csc.nnz);

    std::vector<Index> cursor(csc.col_start.begin(), csc.col_start.begin() + ncols);
    for (Index k = 0; k < ne; ++k) {
        Index pos = cursor[cols[k]]++;
        csc.row_index[pos] = rows[k];
        csc.values[pos]    = vals[k];
    }
    // Sort each column by row index
    for (Index j = 0; j < ncols; ++j) {
        Index s = csc.col_start[j], e = csc.col_start[j + 1];
        // simple insertion sort (columns usually small)
        for (Index a = s + 1; a < e; ++a)
            for (Index b = a; b > s && csc.row_index[b] < csc.row_index[b - 1]; --b) {
                std::swap(csc.row_index[b], csc.row_index[b - 1]);
                std::swap(csc.values[b],    csc.values[b - 1]);
            }
    }
    return csc;
}

SparseMatrixCSC TripletMatrix::to_csc() const {
    return SparseMatrixCSC::from_triplets(num_rows, num_cols, rows, cols, values);
}

void SparseMatrixCSR::multiply_vec(const Real* x, Real* y) const {
    for (Index i = 0; i < num_rows; ++i) {
        Real s = 0.0;
        for (Index k = row_start[i]; k < row_start[i + 1]; ++k)
            s += values[k] * x[col_index[k]];
        y[i] = s;
    }
}

SparseMatrixCSC SparseMatrixCSR::to_csc() const {
    // Transpose of CSR → CSC
    SparseMatrixCSC csc;
    csc.num_rows = num_rows;
    csc.num_cols = num_cols;
    csc.nnz      = nnz;
    csc.col_start.assign(num_cols + 1, 0);
    csc.row_index.resize(nnz);
    csc.values.resize(nnz);

    for (Index k = 0; k < nnz; ++k)
        csc.col_start[col_index[k] + 1]++;
    for (Index j = 0; j < num_cols; ++j)
        csc.col_start[j + 1] += csc.col_start[j];

    std::vector<Index> cursor(csc.col_start.begin(), csc.col_start.begin() + num_cols);
    for (Index i = 0; i < num_rows; ++i)
        for (Index k = row_start[i]; k < row_start[i + 1]; ++k) {
            Index j   = col_index[k];
            Index pos = cursor[j]++;
            csc.row_index[pos] = i;
            csc.values[pos]    = values[k];
        }
    return csc;
}

} // namespace suplex
