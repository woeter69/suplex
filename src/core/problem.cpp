#include "problem.h"
#include <cassert>

namespace suplex {

Problem::Problem(Index nrows, Index ncols) {
    set_dimensions(nrows, ncols);
}

void Problem::set_dimensions(Index rows, Index cols) {
    num_rows_ = rows;
    num_cols_ = cols;
    A_.num_rows = rows;
    A_.num_cols = cols;
    A_.col_start.assign(cols + 1, 0);
    c_.assign(cols, 0.0);
    col_lower_.assign(cols, 0.0);
    col_upper_.assign(cols, INF);
    row_lower_.assign(rows, -INF);
    row_upper_.assign(rows, INF);
    var_types_.assign(cols, VarType::CONTINUOUS);
}

Index Problem::num_integers() const {
    Index count = 0;
    for (VarType vt : var_types_) {
        if (vt == VarType::INTEGER || vt == VarType::BINARY) {
            count++;
        }
    }
    return count;
}

bool Problem::is_mip() const {
    return num_integers() > 0;
}

void Problem::set_constraint_matrix(SparseMatrixCSC&& A) {
    A_ = std::move(A);
    num_rows_ = A_.num_rows;
    num_cols_ = A_.num_cols;
}

void Problem::set_constraint_matrix(const SparseMatrixCSC& A) {
    A_ = A;
    num_rows_ = A_.num_rows;
    num_cols_ = A_.num_cols;
}

void Problem::set_objective(std::vector<Real>&& c) {
    c_ = std::move(c);
}

void Problem::set_objective(const std::vector<Real>& c) {
    c_ = c;
}

void Problem::set_col_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper) {
    col_lower_ = std::move(lower);
    col_upper_ = std::move(upper);
}

void Problem::set_col_bounds(const std::vector<Real>& lower, const std::vector<Real>& upper) {
    col_lower_ = lower;
    col_upper_ = upper;
}

void Problem::set_row_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper) {
    row_lower_ = std::move(lower);
    row_upper_ = std::move(upper);
}

void Problem::set_row_bounds(const std::vector<Real>& lower, const std::vector<Real>& upper) {
    row_lower_ = lower;
    row_upper_ = upper;
}

void Problem::set_var_types(std::vector<VarType>&& types) {
    var_types_ = std::move(types);
}

void Problem::set_var_types(const std::vector<VarType>& types) {
    var_types_ = types;
}

void Problem::set_quadratic_objective(SparseMatrixCSC&& Q) {
    Q_ = std::move(Q);
}

void Problem::set_quadratic_objective(const SparseMatrixCSC& Q) {
    Q_ = Q;
}

void Problem::set_col_bounds(Index col, Real lower, Real upper) {
    if (col >= 0 && col < num_cols_) {
        col_lower_[col] = lower;
        col_upper_[col] = upper;
    }
}

void Problem::set_row_bounds(Index row, Real lower, Real upper) {
    if (row >= 0 && row < num_rows_) {
        row_lower_[row] = lower;
        row_upper_[row] = upper;
    }
}

void Problem::set_objective_coeff(Index col, Real coeff) {
    if (col >= 0 && col < num_cols_) {
        c_[col] = coeff;
    }
}

bool Problem::validate() const {
    if (num_rows_ < 0 || num_cols_ < 0) return false;
    if (A_.num_rows != num_rows_ || A_.num_cols != num_cols_) return false;
    if (static_cast<Index>(A_.col_start.size()) != num_cols_ + 1) return false;
    if (static_cast<Index>(c_.size()) != num_cols_) return false;
    if (static_cast<Index>(col_lower_.size()) != num_cols_) return false;
    if (static_cast<Index>(col_upper_.size()) != num_cols_) return false;
    if (static_cast<Index>(row_lower_.size()) != num_rows_) return false;
    if (static_cast<Index>(row_upper_.size()) != num_rows_) return false;
    if (static_cast<Index>(var_types_.size()) != num_cols_) return false;

    for (Index j = 0; j < num_cols_; ++j) {
        if (col_lower_[j] > col_upper_[j]) return false;
    }
    for (Index i = 0; i < num_rows_; ++i) {
        if (row_lower_[i] > row_upper_[i]) return false;
    }
    if (Q_.has_value()) {
        if (Q_->num_rows != num_cols_ || Q_->num_cols != num_cols_) return false;
    }
    return true;
}

} // namespace suplex
