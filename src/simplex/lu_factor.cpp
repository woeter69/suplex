#include "lu_factor.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace suplex {

SolverStatus SparseLU::factorize(const SparseMatrixCSC& A, const std::vector<Index>& basis, Real threshold) {
    m_ = static_cast<Index>(basis.size());
    is_valid_ = false;
    if (m_ == 0) {
        return SolverStatus::OPTIMAL;
    }

    TripletMatrix trip(m_, m_);
    trip.reserve(m_ * 10);
    for (Index j = 0; j < m_; ++j) {
        Index col_idx = basis[j];
        Index begin = A.col_begin(col_idx);
        Index end = A.col_end(col_idx);
        for (Index k = begin; k < end; ++k) {
            Index r = A.row_index[k];
            if (r < m_) {
                trip.add_entry(r, j, A.values[k]);
            }
        }
    }
    SparseMatrixCSC B = trip.to_csc();
    return factorize_matrix(B, threshold);
}

SolverStatus SparseLU::factorize_matrix(const SparseMatrixCSC& B, Real threshold) {
    m_ = B.num_rows;
    is_valid_ = false;
    if (m_ == 0 || B.num_cols != m_) {
        return SolverStatus::NUMERICAL_ERROR;
    }

    p_.assign(m_, -1);
    q_.assign(m_, -1);
    p_inv_.assign(m_, -1);
    q_inv_.assign(m_, -1);
    U_diag_.assign(m_, 0.0);
    L_rows_.assign(m_, {});
    U_rows_.assign(m_, {});

    // Working dense matrix M initialized from B
    std::vector<std::vector<Real>> M(m_, std::vector<Real>(m_, 0.0));
    Real max_b_elem = 1.0;
    for (Index j = 0; j < m_; ++j) {
        Index begin = B.col_begin(j);
        Index end = B.col_end(j);
        for (Index k = begin; k < end; ++k) {
            Index r = B.row_index[k];
            Real v = B.values[k];
            M[r][j] = v;
            max_b_elem = std::max(max_b_elem, std::abs(v));
        }
    }

    std::vector<bool> row_active(m_, true);
    std::vector<bool> col_active(m_, true);

    std::vector<Index> row_deg(m_, 0);
    std::vector<Index> col_deg(m_, 0);
    for (Index i = 0; i < m_; ++i) {
        for (Index j = 0; j < m_; ++j) {
            if (std::abs(M[i][j]) > EPS_ZERO) {
                row_deg[i]++;
                col_deg[j]++;
            }
        }
    }

    Real max_u_elem = 0.0;

    for (Index step = 0; step < m_; ++step) {
        // Step 1: Markowitz search
        Index best_row = -1;
        Index best_col = -1;
        int64_t best_cost = std::numeric_limits<int64_t>::max();
        Real best_pivot_val = 0.0;

        for (Index j = 0; j < m_; ++j) {
            if (!col_active[j] || col_deg[j] == 0) continue;

            Real max_abs_in_col = 0.0;
            for (Index i = 0; i < m_; ++i) {
                if (row_active[i]) {
                    max_abs_in_col = std::max(max_abs_in_col, std::abs(M[i][j]));
                }
            }

            if (max_abs_in_col < EPS_PIVOT) continue;

            const Real thresh_val = threshold * max_abs_in_col;

            for (Index i = 0; i < m_; ++i) {
                if (!row_active[i]) continue;
                Real abs_v = std::abs(M[i][j]);
                if (abs_v >= thresh_val && abs_v >= EPS_PIVOT) {
                    int64_t cost = static_cast<int64_t>(row_deg[i] - 1) * static_cast<int64_t>(col_deg[j] - 1);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_row = i;
                        best_col = j;
                        best_pivot_val = M[i][j];
                        if (best_cost == 0) break;
                    }
                }
            }
            if (best_cost == 0) break;
        }

        // Secondary search if threshold was too strict
        if (best_row == -1) {
            for (Index j = 0; j < m_; ++j) {
                if (!col_active[j]) continue;
                for (Index i = 0; i < m_; ++i) {
                    if (!row_active[i]) continue;
                    if (std::abs(M[i][j]) >= EPS_PIVOT) {
                        best_row = i;
                        best_col = j;
                        best_pivot_val = M[i][j];
                        break;
                    }
                }
                if (best_row != -1) break;
            }
        }

        if (best_row == -1 || std::abs(best_pivot_val) < EPS_PIVOT) {
            return SolverStatus::NUMERICAL_ERROR;
        }

        p_[step] = best_row;
        q_[step] = best_col;
        p_inv_[best_row] = step;
        q_inv_[best_col] = step;
        U_diag_[step] = best_pivot_val;
        max_u_elem = std::max(max_u_elem, std::abs(best_pivot_val));

        row_active[best_row] = false;
        col_active[best_col] = false;

        // Eliminate in remaining active rows
        for (Index r = 0; r < m_; ++r) {
            if (!row_active[r]) continue;
            Real val_rj = M[r][best_col];
            if (std::abs(val_rj) <= EPS_ZERO) continue;

            Real mult = val_rj / best_pivot_val;
            M[r][best_col] = mult; // store multiplier in-place

            for (Index c = 0; c < m_; ++c) {
                if (!col_active[c]) continue;
                Real val_pc = M[best_row][c];
                if (std::abs(val_pc) > EPS_ZERO) {
                    M[r][c] -= mult * val_pc;
                    if (std::abs(M[r][c]) <= EPS_ZERO) {
                        M[r][c] = 0.0;
                    } else {
                        max_u_elem = std::max(max_u_elem, std::abs(M[r][c]));
                    }
                }
            }
        }

        // Update active degrees
        for (Index i = 0; i < m_; ++i) {
            if (!row_active[i]) {
                row_deg[i] = 0;
            } else {
                Index deg = 0;
                for (Index j = 0; j < m_; ++j) {
                    if (col_active[j] && std::abs(M[i][j]) > EPS_ZERO) deg++;
                }
                row_deg[i] = deg;
            }
        }
        for (Index j = 0; j < m_; ++j) {
            if (!col_active[j]) {
                col_deg[j] = 0;
            } else {
                Index deg = 0;
                for (Index i = 0; i < m_; ++i) {
                    if (row_active[i] && std::abs(M[i][j]) > EPS_ZERO) deg++;
                }
                col_deg[j] = deg;
            }
        }
    }

    // Now convert M into row-wise L and U structures in step coordinates
    for (Index k = 0; k < m_; ++k) {
        Index r = p_[k];
        U_diag_[k] = M[r][q_[k]];

        // Entries for L: l < k
        for (Index l = 0; l < k; ++l) {
            Index c = q_[l];
            Real val = M[r][c];
            if (std::abs(val) > EPS_ZERO) {
                L_rows_[k].emplace_back(l, val);
            }
        }

        // Entries for U: l > k
        for (Index l = k + 1; l < m_; ++l) {
            Index c = q_[l];
            Real val = M[r][c];
            if (std::abs(val) > EPS_ZERO) {
                U_rows_[k].emplace_back(l, val);
            }
        }
    }

    condition_est_ = max_u_elem / std::max(1.0, max_b_elem);
    is_valid_ = true;
    return SolverStatus::OPTIMAL;
}

void SparseLU::ftran(const Real* b, Real* x) const {
    if (!is_valid_ || m_ == 0) return;

    // 1. Forward solve L * z = P * b
    // z_k = b[p_[k]] - sum_{l < k} L_{kl} * z_l
    std::vector<Real> z(m_);
    for (Index k = 0; k < m_; ++k) {
        Real sum = b[p_[k]];
        for (const auto& [l, val] : L_rows_[k]) {
            sum -= val * z[l];
        }
        z[k] = sum;
    }

    // 2. Backward solve U * w = z
    // w_k = (z_k - sum_{l > k} U_{kl} * w_l) / U_{diag}[k]
    std::vector<Real> w(m_, 0.0);
    for (Index k = m_ - 1; k >= 0; --k) {
        Real sum = z[k];
        for (const auto& [l, val] : U_rows_[k]) {
            sum -= val * w[l];
        }
        w[k] = sum / U_diag_[k];
    }

    // 3. x = Q * w => x[q_[k]] = w[k]
    for (Index k = 0; k < m_; ++k) {
        x[q_[k]] = w[k];
    }
}

void SparseLU::ftran(std::vector<Real>& b_in_x_out) const {
    std::vector<Real> x(m_);
    ftran(b_in_x_out.data(), x.data());
    b_in_x_out = std::move(x);
}

void SparseLU::btran(const Real* c, Real* y) const {
    if (!is_valid_ || m_ == 0) return;

    // 1. Forward solve U^T * v = Q^T * c
    // v starts as c permuted by Q: v_k = c[q_[k]]
    // For each k from 0 to m-1:
    // v_k /= U_diag_[k]
    // for each (l, U_kl) with l > k:
    // v_l -= U_kl * v_k
    std::vector<Real> v(m_);
    for (Index k = 0; k < m_; ++k) {
        v[k] = c[q_[k]];
    }

    for (Index k = 0; k < m_; ++k) {
        v[k] /= U_diag_[k];
        Real vk = v[k];
        if (std::abs(vk) <= EPS_ZERO) continue;
        for (const auto& [l, val] : U_rows_[k]) {
            v[l] -= val * vk;
        }
    }

    // 2. Backward solve L^T * u = v
    // u starts as v
    // For each k from m-1 down to 0:
    // for each (l, L_kl) with l < k:
    // u_l -= L_kl * u_k
    std::vector<Real> u = v;
    for (Index k = m_ - 1; k >= 0; --k) {
        Real uk = u[k];
        if (std::abs(uk) <= EPS_ZERO) continue;
        for (const auto& [l, val] : L_rows_[k]) {
            u[l] -= val * uk;
        }
    }

    // 3. y = P^T * u => y[p_[k]] = u[k]
    for (Index k = 0; k < m_; ++k) {
        y[p_[k]] = u[k];
    }
}

void SparseLU::btran(std::vector<Real>& c_in_y_out) const {
    std::vector<Real> y(m_);
    btran(c_in_y_out.data(), y.data());
    c_in_y_out = std::move(y);
}

void SparseLU::iterative_refine(const SparseMatrixCSC& A, const std::vector<Index>& basis,
                                const Real* b, Real* x, int max_rounds) const {
    if (!is_valid_ || m_ == 0) return;

    std::vector<Real> r(m_);
    std::vector<Real> d(m_);

    for (int round = 0; round < max_rounds; ++round) {
        for (Index i = 0; i < m_; ++i) {
            r[i] = b[i];
        }
        for (Index j = 0; j < m_; ++j) {
            const Real xj = x[j];
            if (std::abs(xj) <= EPS_ZERO) continue;
            Index col_idx = basis[j];
            Index begin = A.col_begin(col_idx);
            Index end = A.col_end(col_idx);
            for (Index k = begin; k < end; ++k) {
                Index row = A.row_index[k];
                if (row < m_) {
                    r[row] -= A.values[k] * xj;
                }
            }
        }

        Real max_r = 0.0;
        Real max_b = 1.0;
        for (Index i = 0; i < m_; ++i) {
            max_r = std::max(max_r, std::abs(r[i]));
            max_b = std::max(max_b, std::abs(b[i]));
        }

        if (max_r <= EPS_FEASIBILITY * max_b) {
            break;
        }

        ftran(r.data(), d.data());
        for (Index j = 0; j < m_; ++j) {
            x[j] += d[j];
        }
    }
}

} // namespace suplex
