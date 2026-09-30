#pragma once

#include "src/core/types.h"
#include "src/core/sparse_matrix.h"
#include <vector>

namespace suplex {

/**
 * @brief Sparse LU Factorization using Markowitz ordering and threshold partial pivoting.
 * Factors P * B * Q = L * U, where P and Q are permutation matrices,
 * L is unit lower triangular, and U is upper triangular.
 */
class SparseLU {
public:
    SparseLU() = default;

    /**
     * @brief Factorizes basis matrix B defined by selected columns of A.
     * @param A Full constraint matrix.
     * @param basis List of m column indices of A forming B.
     * @param threshold Markowitz threshold (default 0.1).
     * @return SolverStatus::OPTIMAL on success, SolverStatus::NUMERICAL_ERROR on singularity.
     */
    SolverStatus factorize(const SparseMatrixCSC& A, const std::vector<Index>& basis, Real threshold = 0.1);

    /**
     * @brief Direct factorization of an m x m CSC matrix.
     */
    SolverStatus factorize_matrix(const SparseMatrixCSC& B, Real threshold = 0.1);

    /**
     * @brief Forward solve: Solves B * x = b.
     */
    void ftran(const Real* b, Real* x) const;
    void ftran(std::vector<Real>& b_in_x_out) const;

    /**
     * @brief Backward solve: Solves B^T * y = c.
     */
    void btran(const Real* c, Real* y) const;
    void btran(std::vector<Real>& c_in_y_out) const;

    /**
     * @brief Performs iterative refinement to polish solution to B * x = b.
     */
    void iterative_refine(const SparseMatrixCSC& A, const std::vector<Index>& basis,
                          const Real* b, Real* x, int max_rounds = 3) const;

    Index dimension() const { return m_; }
    bool is_valid() const { return is_valid_; }
    Real condition_estimate() const { return condition_est_; }

    const std::vector<Index>& row_perm() const { return p_; }
    const std::vector<Index>& col_perm() const { return q_; }
    const std::vector<Index>& row_perm_inv() const { return p_inv_; }
    const std::vector<Index>& col_perm_inv() const { return q_inv_; }

private:
    Index m_ = 0;
    bool is_valid_ = false;
    Real condition_est_ = 1.0;

    std::vector<Index> p_;     // p_[k] = original row of pivot step k
    std::vector<Index> q_;     // q_[k] = original col of pivot step k
    std::vector<Index> p_inv_; // p_inv_[row] = step
    std::vector<Index> q_inv_; // q_inv_[col] = step

    // L stored by rows in step coordinates: row k contains pairs (l, L_kl) for l < k
    std::vector<std::vector<std::pair<Index, Real>>> L_rows_;
    // U stored by rows in step coordinates: row k contains pairs (l, U_kl) for l > k
    std::vector<std::vector<std::pair<Index, Real>>> U_rows_;
    // U diagonal in step coordinates
    std::vector<Real> U_diag_;
};

} // namespace suplex
