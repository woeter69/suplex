#include "basis.h"
#include <cmath>
#include <algorithm>

namespace suplex {

Basis::Basis(Index num_rows, Index total_vars) {
    init(num_rows, total_vars);
}

void Basis::init(Index num_rows, Index total_vars) {
    num_rows_ = num_rows;
    total_vars_ = total_vars;
    basic_vars_.assign(num_rows, -1);
    var_status_.assign(total_vars, BasisStatus::AT_LOWER);
    basic_row_.assign(total_vars, -1);
}

void Basis::set_status(Index var, BasisStatus status) {
    if (var >= 0 && var < total_vars_) {
        var_status_[var] = status;
        if (status != BasisStatus::BASIC) {
            basic_row_[var] = -1;
        }
    }
}

void Basis::pivot(Index leaving_row, Index entering_var, BasisStatus leaving_status) {
    if (leaving_row < 0 || leaving_row >= num_rows_ || entering_var < 0 || entering_var >= total_vars_) {
        return;
    }
    Index leaving_var = basic_vars_[leaving_row];
    if (leaving_var >= 0 && leaving_var < total_vars_) {
        var_status_[leaving_var] = leaving_status;
        basic_row_[leaving_var] = -1;
    }

    basic_vars_[leaving_row] = entering_var;
    var_status_[entering_var] = BasisStatus::BASIC;
    basic_row_[entering_var] = leaving_row;
}

void Basis::set_logical_basis(Index num_cols) {
    for (Index i = 0; i < num_rows_; ++i) {
        Index slack_var = num_cols + i;
        basic_vars_[i] = slack_var;
        var_status_[slack_var] = BasisStatus::BASIC;
        basic_row_[slack_var] = i;
    }
    for (Index j = 0; j < num_cols; ++j) {
        var_status_[j] = BasisStatus::AT_LOWER;
        basic_row_[j] = -1;
    }
}

void Basis::crash(
    const SparseMatrixCSC& A_augmented,
    const std::vector<Real>& col_lower,
    const std::vector<Real>& col_upper
) {
    Index m = num_rows_;
    Index n_aug = total_vars_;
    Index num_structural = n_aug - m;

    std::vector<bool> row_assigned(m, false);
    std::vector<bool> col_assigned(n_aug, false);

    // Pass 1: Singleton columns in structural part
    for (Index j = 0; j < num_structural; ++j) {
        if (col_lower[j] == col_upper[j]) continue; // skip fixed variables

        Index count = 0;
        Index last_row = -1;
        Real last_val = 0.0;
        Index begin = A_augmented.col_begin(j);
        Index end = A_augmented.col_end(j);

        for (Index k = begin; k < end; ++k) {
            Index r = A_augmented.row_index[k];
            if (r < m && !row_assigned[r]) {
                count++;
                last_row = r;
                last_val = A_augmented.values[k];
            }
        }

        if (count == 1 && std::abs(last_val) > EPS_PIVOT) {
            basic_vars_[last_row] = j;
            var_status_[j] = BasisStatus::BASIC;
            basic_row_[j] = last_row;
            row_assigned[last_row] = true;
            col_assigned[j] = true;
        }
    }

    // Pass 2: Two-element columns
    for (Index j = 0; j < num_structural; ++j) {
        if (col_assigned[j] || col_lower[j] == col_upper[j]) continue;

        Index count = 0;
        Index best_row = -1;
        Real best_val = 0.0;
        Index begin = A_augmented.col_begin(j);
        Index end = A_augmented.col_end(j);

        for (Index k = begin; k < end; ++k) {
            Index r = A_augmented.row_index[k];
            if (r < m && !row_assigned[r]) {
                count++;
                Real v = A_augmented.values[k];
                if (std::abs(v) > std::abs(best_val)) {
                    best_row = r;
                    best_val = v;
                }
            }
        }

        if (count == 2 && best_row != -1 && std::abs(best_val) > EPS_PIVOT) {
            basic_vars_[best_row] = j;
            var_status_[j] = BasisStatus::BASIC;
            basic_row_[j] = best_row;
            row_assigned[best_row] = true;
            col_assigned[j] = true;
        }
    }

    // Fill all remaining unassigned rows with logical slack variables
    for (Index i = 0; i < m; ++i) {
        if (!row_assigned[i]) {
            Index slack_col = num_structural + i;
            basic_vars_[i] = slack_col;
            var_status_[slack_col] = BasisStatus::BASIC;
            basic_row_[slack_col] = i;
            row_assigned[i] = true;
            col_assigned[slack_col] = true;
        }
    }

    // Set non-basic variables to appropriate initial status
    for (Index j = 0; j < n_aug; ++j) {
        if (var_status_[j] != BasisStatus::BASIC) {
            if (col_lower[j] == col_upper[j]) {
                var_status_[j] = BasisStatus::FIXED;
            } else if (col_lower[j] > -INF) {
                var_status_[j] = BasisStatus::AT_LOWER;
            } else if (col_upper[j] < INF) {
                var_status_[j] = BasisStatus::AT_UPPER;
            } else {
                var_status_[j] = BasisStatus::FREE_ZERO;
            }
        }
    }
}

} // namespace suplex
