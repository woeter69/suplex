#pragma once

#include "types.h"
#include "sparse_matrix.h"
#include <vector>
#include <optional>
#include <string>

namespace suplex {

/**
 * @brief Primary mathematical optimization problem container.
 * Supports linear constraints with ranged bounds: row_lower <= Ax <= row_upper,
 * variable bounds: col_lower <= x <= col_upper, and optional quadratic objective.
 */
class Problem {
public:
    Problem() = default;
    Problem(Index nrows, Index ncols);

    // Dimension inspectors
    Index num_rows() const { return num_rows_; }
    Index num_cols() const { return num_cols_; }
    Index num_integers() const;
    bool is_mip() const;
    bool is_qp() const { return Q_.has_value(); }

    // Const accessors
    const SparseMatrixCSC& constraint_matrix() const { return A_; }
    const std::vector<Real>& objective() const { return c_; }
    const std::vector<Real>& col_lower() const { return col_lower_; }
    const std::vector<Real>& col_upper() const { return col_upper_; }
    const std::vector<Real>& row_lower() const { return row_lower_; }
    const std::vector<Real>& row_upper() const { return row_upper_; }
    const std::vector<VarType>& var_types() const { return var_types_; }
    ObjectiveSense sense() const { return sense_; }
    const SparseMatrixCSC* quadratic_objective() const {
        return Q_.has_value() ? &Q_.value() : nullptr;
    }

    // Mutable accessors for algorithmic transforms & scaling
    SparseMatrixCSC& mutable_constraint_matrix() { return A_; }
    std::vector<Real>& mutable_objective() { return c_; }
    std::vector<Real>& mutable_col_lower() { return col_lower_; }
    std::vector<Real>& mutable_col_upper() { return col_upper_; }
    std::vector<Real>& mutable_row_lower() { return row_lower_; }
    std::vector<Real>& mutable_row_upper() { return row_upper_; }
    std::vector<VarType>& mutable_var_types() { return var_types_; }

    // Setters / construction
    void set_dimensions(Index rows, Index cols);
    void set_objective_sense(ObjectiveSense s) { sense_ = s; }
    void set_constraint_matrix(SparseMatrixCSC&& A);
    void set_constraint_matrix(const SparseMatrixCSC& A);
    void set_objective(std::vector<Real>&& c);
    void set_objective(const std::vector<Real>& c);
    void set_col_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper);
    void set_col_bounds(const std::vector<Real>& lower, const std::vector<Real>& upper);
    void set_row_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper);
    void set_row_bounds(const std::vector<Real>& lower, const std::vector<Real>& upper);
    void set_var_types(std::vector<VarType>&& types);
    void set_var_types(const std::vector<VarType>& types);
    void set_quadratic_objective(SparseMatrixCSC&& Q);
    void set_quadratic_objective(const SparseMatrixCSC& Q);

    // Helpers
    void set_col_bounds(Index col, Real lower, Real upper);
    void set_row_bounds(Index row, Real lower, Real upper);
    void set_objective_coeff(Index col, Real coeff);

    // Validation
    bool validate() const;

private:
    Index num_rows_ = 0;
    Index num_cols_ = 0;
    SparseMatrixCSC A_;
    std::vector<Real> c_;
    std::vector<Real> col_lower_;
    std::vector<Real> col_upper_;
    std::vector<Real> row_lower_;
    std::vector<Real> row_upper_;
    std::vector<VarType> var_types_;
    ObjectiveSense sense_ = ObjectiveSense::MINIMIZE;
    std::optional<SparseMatrixCSC> Q_;
};

} // namespace suplex
