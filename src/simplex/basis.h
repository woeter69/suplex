#pragma once

#include "src/core/types.h"
#include "src/core/sparse_matrix.h"
#include <vector>

namespace suplex {

enum class BasisStatus : uint8_t {
    BASIC = 0,
    AT_LOWER = 1,
    AT_UPPER = 2,
    FIXED = 3,
    FREE_ZERO = 4
};

inline const char* to_string(BasisStatus s) {
    switch (s) {
        case BasisStatus::BASIC: return "BASIC";
        case BasisStatus::AT_LOWER: return "AT_LOWER";
        case BasisStatus::AT_UPPER: return "AT_UPPER";
        case BasisStatus::FIXED: return "FIXED";
        case BasisStatus::FREE_ZERO: return "FREE_ZERO";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Manages basis state, variable partitioning, and initial basis crash procedures.
 */
class Basis {
public:
    Basis() = default;
    Basis(Index num_rows, Index total_vars);

    void init(Index num_rows, Index total_vars);

    Index num_rows() const { return num_rows_; }
    Index total_vars() const { return total_vars_; }

    const std::vector<Index>& basic_vars() const { return basic_vars_; }
    const std::vector<BasisStatus>& var_status() const { return var_status_; }
    const std::vector<Index>& basic_row() const { return basic_row_; }

    BasisStatus status(Index var) const { return var_status_[var]; }
    Index basic_var_at(Index row) const { return basic_vars_[row]; }
    Index row_of_basic_var(Index var) const { return basic_row_[var]; }
    bool is_basic(Index var) const { return var_status_[var] == BasisStatus::BASIC; }

    void set_status(Index var, BasisStatus status);
    void pivot(Index leaving_row, Index entering_var, BasisStatus leaving_status);

    /**
     * @brief Performs triangular basis crash on the augmented matrix [A, I].
     */
    void crash(const SparseMatrixCSC& A_augmented,
               const std::vector<Real>& col_lower,
               const std::vector<Real>& col_upper);

    /**
     * @brief Sets logical slack variables (columns num_cols to num_cols + num_rows - 1) as basic.
     */
    void set_logical_basis(Index num_cols);

private:
    Index num_rows_ = 0;
    Index total_vars_ = 0;
    std::vector<Index> basic_vars_;       // size m: basic_vars_[i] = var in row i
    std::vector<BasisStatus> var_status_; // size total_vars
    std::vector<Index> basic_row_;        // size total_vars: -1 if non-basic
};

} // namespace suplex
