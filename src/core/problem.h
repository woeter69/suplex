#pragma once
#include "types.h"
#include "sparse_matrix.h"
#include <vector>
#include <optional>
#include <string>

namespace suplex {

// ── Problem representation ────────────────────────────────────────────────────
///
/// Represents a general mathematical program:
///
///   min/max  c^T x  (+  0.5 x^T Q x  for QP)
///   s.t.     row_lower <= A x <= row_upper
///            col_lower <= x <= col_upper
///            x_j ∈ Z  for integer variables
///
/// Ranged rows (row_lower, row_upper) are used instead of Ax = b to avoid
/// an explosion of slack variables and to represent all constraint types
/// uniformly:
///   equality   : row_lower[i] == row_upper[i]
///   <=  (leq)  : row_lower[i] == -INF
///   >=  (geq)  : row_upper[i] ==  INF
///   ranged     : both finite and different
///
class Problem {
public:
    // ── Construction ────────────────────────────────────────────────────────
    void set_dimensions(Index rows, Index cols);
    void set_objective_sense(ObjectiveSense s) { sense_ = s; }
    void set_constraint_matrix(SparseMatrixCSC&& A) { A_ = std::move(A); }
    void set_objective(std::vector<Real>&& c)        { c_ = std::move(c); }
    void set_col_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper);
    void set_row_bounds(std::vector<Real>&& lower, std::vector<Real>&& upper);
    void set_var_types(std::vector<VarType>&& types) { var_types_ = std::move(types); }
    void set_quadratic_objective(SparseMatrixCSC&& Q) { Q_ = std::move(Q); }

    // ── Dimension queries ────────────────────────────────────────────────────
    Index num_rows()     const { return A_.num_rows; }
    Index num_cols()     const { return A_.num_cols; }
    Index num_integers() const;   ///< Count of INTEGER + BINARY variables
    bool  is_mip()       const;   ///< Has any integer variables?
    bool  is_qp()        const { return Q_.has_value(); }

    // ── Data accessors (const refs) ──────────────────────────────────────────
    const SparseMatrixCSC&   constraint_matrix()      const { return A_; }
    const std::vector<Real>& objective()              const { return c_; }
    const std::vector<Real>& col_lower()              const { return col_lower_; }
    const std::vector<Real>& col_upper()              const { return col_upper_; }
    const std::vector<Real>& row_lower()              const { return row_lower_; }
    const std::vector<Real>& row_upper()              const { return row_upper_; }
    const std::vector<VarType>& var_types()           const { return var_types_; }
    ObjectiveSense            sense()                 const { return sense_; }
    const SparseMatrixCSC*   quadratic_objective()    const {
        return Q_.has_value() ? &Q_.value() : nullptr;
    }

    // ── Validation ──────────────────────────────────────────────────────────
    bool validate() const;   ///< Check dimensional consistency

private:
    SparseMatrixCSC          A_;
    std::vector<Real>        c_;
    std::vector<Real>        col_lower_;
    std::vector<Real>        col_upper_;
    std::vector<Real>        row_lower_;
    std::vector<Real>        row_upper_;
    std::vector<VarType>     var_types_;
    ObjectiveSense           sense_ = ObjectiveSense::MINIMIZE;
    std::optional<SparseMatrixCSC> Q_;
};

} // namespace suplex
