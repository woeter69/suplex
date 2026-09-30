#include "sparse_matrix.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace suplex {

SparseMatrixCSC::SparseMatrixCSC(Index nrows, Index ncols, Index nonzeros)
    : num_rows(nrows), num_cols(ncols), nnz(nonzeros) {
    col_start.resize(num_cols + 1, 0);
    if (nnz > 0) {
        row_index.resize(nnz, 0);
        values.resize(nnz, 0.0);
    }
}

void SparseMatrixCSC::clear() {
    num_rows = 0;
    num_cols = 0;
    nnz = 0;
    col_start.clear();
    row_index.clear();
    values.clear();
}

Real SparseMatrixCSC::get(Index row, Index col) const {
    if (col < 0 || col >= num_cols || row < 0 || row >= num_rows) {
        return 0.0;
    }
    const Index begin = col_start[col];
    const Index end = col_start[col + 1];
    auto it = std::lower_bound(row_index.begin() + begin, row_index.begin() + end, row);
    if (it != row_index.begin() + end && *it == row) {
        Index idx = static_cast<Index>(std::distance(row_index.begin(), it));
        return values[idx];
    }
    return 0.0;
}

void SparseMatrixCSC::multiply_vec(const Real* x, Real* y) const {
    std::fill(y, y + num_rows, 0.0);
    for (Index j = 0; j < num_cols; ++j) {
        const Real xj = x[j];
        if (std::abs(xj) <= EPS_ZERO) {
            continue;
        }
        const Index begin = col_start[j];
        const Index end = col_start[j + 1];
        for (Index k = begin; k < end; ++k) {
            y[row_index[k]] += values[k] * xj;
        }
    }
}

void SparseMatrixCSC::multiply_transpose_vec(const Real* x, Real* y) const {
    for (Index j = 0; j < num_cols; ++j) {
        Real sum = 0.0;
        const Index begin = col_start[j];
        const Index end = col_start[j + 1];
        for (Index k = begin; k < end; ++k) {
            sum += values[k] * x[row_index[k]];
        }
        y[j] = sum;
    }
}

void SparseMatrixCSC::multiply_col(Index j, Real scalar, Real* y) const {
    if (j < 0 || j >= num_cols || std::abs(scalar) <= EPS_ZERO) {
        return;
    }
    const Index begin = col_start[j];
    const Index end = col_start[j + 1];
    for (Index k = begin; k < end; ++k) {
        y[row_index[k]] += scalar * values[k];
    }
}

Real SparseMatrixCSC::dot_col(Index j, const Real* x) const {
    if (j < 0 || j >= num_cols) {
        return 0.0;
    }
    Real sum = 0.0;
    const Index begin = col_start[j];
    const Index end = col_start[j + 1];
    for (Index k = begin; k < end; ++k) {
        sum += values[k] * x[row_index[k]];
    }
    return sum;
}

SparseMatrixCSC SparseMatrixCSC::from_triplets(
    Index nrows,
    Index ncols,
    const std::vector<Index>& rows,
    const std::vector<Index>& cols,
    const std::vector<Real>& vals
) {
    SparseMatrixCSC csc;
    csc.num_rows = nrows;
    csc.num_cols = ncols;
    csc.col_start.assign(ncols + 1, 0);

    const size_t num_entries = rows.size();
    if (num_entries == 0 || ncols == 0) {
        csc.nnz = 0;
        return csc;
    }

    // Step 1: Count nonzeros per column
    std::vector<Index> col_counts(ncols, 0);
    for (size_t k = 0; k < num_entries; ++k) {
        Index c = cols[k];
        if (c >= 0 && c < ncols && rows[k] >= 0 && rows[k] < nrows) {
            col_counts[c]++;
        }
    }

    // Step 2: Prefix sum for col_start
    csc.col_start[0] = 0;
    for (Index j = 0; j < ncols; ++j) {
        csc.col_start[j + 1] = csc.col_start[j] + col_counts[j];
    }

    const Index total_candidates = csc.col_start[ncols];
    std::vector<Index> cursor = csc.col_start;
    csc.row_index.resize(total_candidates);
    csc.values.resize(total_candidates);

    // Step 3: Populate arrays
    for (size_t k = 0; k < num_entries; ++k) {
        Index r = rows[k];
        Index c = cols[k];
        if (c >= 0 && c < ncols && r >= 0 && r < nrows) {
            Index pos = cursor[c]++;
            csc.row_index[pos] = r;
            csc.values[pos] = vals[k];
        }
    }

    // Step 4: Sort and sum duplicates per column
    std::vector<std::pair<Index, Real>> col_entries;
    std::vector<Index> new_col_start(ncols + 1, 0);
    std::vector<Index> final_rows;
    std::vector<Real> final_vals;
    final_rows.reserve(total_candidates);
    final_vals.reserve(total_candidates);

    for (Index j = 0; j < ncols; ++j) {
        Index start = csc.col_start[j];
        Index end = csc.col_start[j + 1];
        col_entries.clear();
        for (Index k = start; k < end; ++k) {
            col_entries.emplace_back(csc.row_index[k], csc.values[k]);
        }
        std::sort(col_entries.begin(), col_entries.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        new_col_start[j] = static_cast<Index>(final_rows.size());
        for (size_t k = 0; k < col_entries.size(); ++k) {
            if (!final_rows.empty() &&
                static_cast<Index>(final_rows.size()) > new_col_start[j] &&
                final_rows.back() == col_entries[k].first) {
                final_vals.back() += col_entries[k].second;
            } else {
                final_rows.push_back(col_entries[k].first);
                final_vals.push_back(col_entries[k].second);
            }
        }
    }
    new_col_start[ncols] = static_cast<Index>(final_rows.size());

    csc.col_start = std::move(new_col_start);
    csc.row_index = std::move(final_rows);
    csc.values = std::move(final_vals);
    csc.nnz = static_cast<Index>(csc.values.size());

    return csc;
}

SparseMatrixCSR SparseMatrixCSC::to_csr() const {
    SparseMatrixCSR csr;
    csr.num_rows = num_rows;
    csr.num_cols = num_cols;
    csr.nnz = nnz;
    csr.row_start.assign(num_rows + 1, 0);
    csr.col_index.resize(nnz);
    csr.values.resize(nnz);

    if (nnz == 0 || num_rows == 0) {
        return csr;
    }

    // Step 1: Count row frequencies
    std::vector<Index> row_count(num_rows, 0);
    for (Index k = 0; k < nnz; ++k) {
        row_count[row_index[k]]++;
    }

    // Step 2: Prefix sum to row_start
    csr.row_start[0] = 0;
    for (Index i = 0; i < num_rows; ++i) {
        csr.row_start[i + 1] = csr.row_start[i] + row_count[i];
    }

    // Step 3: Populate CSR
    std::vector<Index> cursor = csr.row_start;
    for (Index j = 0; j < num_cols; ++j) {
        Index begin = col_start[j];
        Index end = col_start[j + 1];
        for (Index k = begin; k < end; ++k) {
            Index r = row_index[k];
            Index pos = cursor[r]++;
            csr.col_index[pos] = j;
            csr.values[pos] = values[k];
        }
    }

    return csr;
}

SparseMatrixCSR::SparseMatrixCSR(Index nrows, Index ncols, Index nonzeros)
    : num_rows(nrows), num_cols(ncols), nnz(nonzeros) {
    row_start.resize(num_rows + 1, 0);
    if (nnz > 0) {
        col_index.resize(nnz, 0);
        values.resize(nnz, 0.0);
    }
}

void SparseMatrixCSR::clear() {
    num_rows = 0;
    num_cols = 0;
    nnz = 0;
    row_start.clear();
    col_index.clear();
    values.clear();
}

Real SparseMatrixCSR::get(Index row, Index col) const {
    if (row < 0 || row >= num_rows || col < 0 || col >= num_cols) {
        return 0.0;
    }
    const Index begin = row_start[row];
    const Index end = row_start[row + 1];
    auto it = std::lower_bound(col_index.begin() + begin, col_index.begin() + end, col);
    if (it != col_index.begin() + end && *it == col) {
        Index idx = static_cast<Index>(std::distance(col_index.begin(), it));
        return values[idx];
    }
    return 0.0;
}

void SparseMatrixCSR::multiply_vec(const Real* x, Real* y) const {
    for (Index i = 0; i < num_rows; ++i) {
        Real sum = 0.0;
        const Index begin = row_start[i];
        const Index end = row_start[i + 1];
        for (Index k = begin; k < end; ++k) {
            sum += values[k] * x[col_index[k]];
        }
        y[i] = sum;
    }
}

SparseMatrixCSC SparseMatrixCSR::to_csc() const {
    SparseMatrixCSC csc;
    csc.num_rows = num_rows;
    csc.num_cols = num_cols;
    csc.nnz = nnz;
    csc.col_start.assign(num_cols + 1, 0);
    csc.row_index.resize(nnz);
    csc.values.resize(nnz);

    if (nnz == 0 || num_cols == 0) {
        return csc;
    }

    std::vector<Index> col_count(num_cols, 0);
    for (Index k = 0; k < nnz; ++k) {
        col_count[col_index[k]]++;
    }

    csc.col_start[0] = 0;
    for (Index j = 0; j < num_cols; ++j) {
        csc.col_start[j + 1] = csc.col_start[j] + col_count[j];
    }

    std::vector<Index> cursor = csc.col_start;
    for (Index i = 0; i < num_rows; ++i) {
        Index begin = row_start[i];
        Index end = row_start[i + 1];
        for (Index k = begin; k < end; ++k) {
            Index c = col_index[k];
            Index pos = cursor[c]++;
            csc.row_index[pos] = i;
            csc.values[pos] = values[k];
        }
    }

    return csc;
}

void TripletMatrix::add_entry(Index row, Index col, Real val) {
    rows.push_back(row);
    cols.push_back(col);
    values.push_back(val);
}

void TripletMatrix::reserve(Index nnz_est) {
    rows.reserve(nnz_est);
    cols.reserve(nnz_est);
    values.reserve(nnz_est);
}

void TripletMatrix::clear() {
    num_rows = 0;
    num_cols = 0;
    rows.clear();
    cols.clear();
    values.clear();
}

SparseMatrixCSC TripletMatrix::to_csc() const {
    return SparseMatrixCSC::from_triplets(num_rows, num_cols, rows, cols, values);
}

} // namespace suplex
