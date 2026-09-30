#include "lu_update.h"
#include <cmath>
#include <algorithm>

namespace suplex {

SolverStatus LUManager::refactorize(const SparseMatrixCSC& A, const std::vector<Index>& basis) {
    etas_.clear();
    max_eta_growth_ = 1.0;
    SolverStatus status = base_lu_.factorize(A, basis);
    Index m = static_cast<Index>(basis.size());
    max_updates_ = std::max<Index>(50, static_cast<Index>(std::sqrt(m * 2)));
    return status;
}

SolverStatus LUManager::update(Index leaving_pos, const std::vector<Real>& alpha) {
    Index m = base_lu_.dimension();
    if (leaving_pos < 0 || leaving_pos >= m) {
        return SolverStatus::NUMERICAL_ERROR;
    }

    Real pivot = alpha[leaving_pos];
    if (std::abs(pivot) < EPS_PIVOT) {
        return SolverStatus::NUMERICAL_ERROR;
    }

    EtaMatrix eta;
    eta.leaving_pos = leaving_pos;
    eta.pivot_val = pivot;

    for (Index i = 0; i < m; ++i) {
        Real val = alpha[i];
        if (std::abs(val) > EPS_ZERO || i == leaving_pos) {
            eta.indices.push_back(i);
            eta.values.push_back(val);
            max_eta_growth_ = std::max(max_eta_growth_, std::abs(val));
        }
    }

    etas_.push_back(std::move(eta));
    return SolverStatus::OPTIMAL;
}

void LUManager::ftran(const Real* b, Real* x) const {
    Index m = base_lu_.dimension();
    if (m == 0) return;

    // 1. Solve base factorization: x = B_0^{-1} * b
    base_lu_.ftran(b, x);

    // 2. Apply Eta transformations in forward order: x = E_{K-1}^{-1} ... E_0^{-1} x
    for (const auto& eta : etas_) {
        Index p = eta.leaving_pos;
        Real xp = x[p];
        if (std::abs(xp) <= EPS_ZERO) continue;

        Real mult = xp / eta.pivot_val;
        // x[p] is updated to mult
        // for i != p: x[i] -= mult * eta[i]
        const size_t nnz = eta.indices.size();
        for (size_t k = 0; k < nnz; ++k) {
            Index i = eta.indices[k];
            if (i != p) {
                x[i] -= mult * eta.values[k];
            }
        }
        x[p] = mult;
    }
}

void LUManager::ftran(std::vector<Real>& b_in_x_out) const {
    std::vector<Real> x(base_lu_.dimension());
    ftran(b_in_x_out.data(), x.data());
    b_in_x_out = std::move(x);
}

void LUManager::btran(const Real* c, Real* y) const {
    Index m = base_lu_.dimension();
    if (m == 0) return;

    // To solve B_K^T y = c where B_K = B_0 E_0 ... E_{K-1}:
    // (B_0 E_0 ... E_{K-1})^T y = c
    // E_{K-1}^T ... E_0^T (B_0^T y) = c
    // Let w_K = c.
    // For k = K-1 down to 0: w_k = (E_k^{-1})^T w_{k+1}.
    // Then B_0^T y = w_0 => y = B_0^{-T} w_0.

    std::vector<Real> w(c, c + m);

    // Apply (E_k^{-1})^T in reverse order
    for (int64_t k = static_cast<int64_t>(etas_.size()) - 1; k >= 0; --k) {
        const auto& eta = etas_[k];
        Index p = eta.leaving_pos;

        // (E^{-1})^T w:
        // for i != p: w'[i] = w[i]
        // w'[p] = (w[p] - sum_{i != p} eta[i] * w[i]) / eta[p]
        Real sum = w[p];
        const size_t nnz = eta.indices.size();
        for (size_t idx = 0; idx < nnz; ++idx) {
            Index i = eta.indices[idx];
            if (i != p) {
                sum -= eta.values[idx] * w[i];
            }
        }
        w[p] = sum / eta.pivot_val;
    }

    base_lu_.btran(w.data(), y);
}

void LUManager::btran(std::vector<Real>& c_in_y_out) const {
    std::vector<Real> y(base_lu_.dimension());
    btran(c_in_y_out.data(), y.data());
    c_in_y_out = std::move(y);
}

bool LUManager::should_refactorize() const {
    if (static_cast<Index>(etas_.size()) >= max_updates_) {
        return true;
    }
    if (max_eta_growth_ > 1e10) {
        return true;
    }
    return false;
}

} // namespace suplex
