#pragma once

#include "src/core/types.h"
#include "lu_factor.h"
#include <vector>
#include <memory>

namespace suplex {

/**
 * @brief Representation of an elementary column update (Eta matrix).
 * E = I + (eta - e_p) * e_p^T
 */
struct EtaMatrix {
    Index leaving_pos;
    Real pivot_val;
    std::vector<Index> indices; // nonzeros of eta
    std::vector<Real> values;
};

/**
 * @brief Manages basis factorization and sequential Forrest-Tomlin / Product Form Eta updates.
 */
class LUManager {
public:
    LUManager() = default;

    /**
     * @brief Performs full refactorization of basis columns.
     */
    SolverStatus refactorize(const SparseMatrixCSC& A, const std::vector<Index>& basis);

    /**
     * @brief Adds an update transformation after a pivot.
     * @param leaving_pos Position in the basis being replaced (0 to m-1).
     * @param alpha Pivot column: alpha = B^{-1} * a_entering (from FTRAN).
     * @return SolverStatus::OPTIMAL on success, NUMERICAL_ERROR on near-zero pivot.
     */
    SolverStatus update(Index leaving_pos, const std::vector<Real>& alpha);

    /**
     * @brief Forward solve: Solves B_current * x = b.
     */
    void ftran(const Real* b, Real* x) const;
    void ftran(std::vector<Real>& b_in_x_out) const;

    /**
     * @brief Backward solve: Solves B_current^T * y = c.
     */
    void btran(const Real* c, Real* y) const;
    void btran(std::vector<Real>& c_in_y_out) const;

    /**
     * @brief Evaluates whether a complete refactorization should be triggered.
     */
    bool should_refactorize() const;

    Index num_updates() const { return static_cast<Index>(etas_.size()); }
    Index dimension() const { return base_lu_.dimension(); }
    bool is_valid() const { return base_lu_.is_valid(); }
    const SparseLU& base_lu() const { return base_lu_; }

    void set_max_updates(Index max_up) { max_updates_ = max_up; }

private:
    SparseLU base_lu_;
    std::vector<EtaMatrix> etas_;
    Index max_updates_ = 50;
    Real max_eta_growth_ = 1.0;
};

} // namespace suplex
