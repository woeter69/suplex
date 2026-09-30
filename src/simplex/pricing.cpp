#include "pricing.h"
#include <cmath>
#include <algorithm>

namespace suplex {

Pricing::Pricing(PricingStrategy strategy) : strategy_(strategy) {}

void Pricing::init(Index total_vars, Index num_rows, PricingStrategy strategy) {
    total_vars_ = total_vars;
    num_rows_ = num_rows;
    strategy_ = strategy;
    reset_weights();
}

void Pricing::reset_weights() {
    weights_.assign(total_vars_, 1.0);
    iterations_since_reset_ = 0;
}

Index Pricing::select_entering(
    const std::vector<Real>& reduced_costs,
    const std::vector<BasisStatus>& var_status
) {
    Index best_j = -1;
    Real best_score = -1.0;

    for (Index j = 0; j < total_vars_; ++j) {
        BasisStatus st = var_status[j];
        if (st == BasisStatus::BASIC || st == BasisStatus::FIXED) {
            continue;
        }

        Real rc = reduced_costs[j];
        bool eligible = false;
        Real violation = 0.0;

        if (st == BasisStatus::AT_LOWER) {
            if (rc < -EPS_OPTIMALITY) {
                eligible = true;
                violation = -rc;
            }
        } else if (st == BasisStatus::AT_UPPER) {
            if (rc > EPS_OPTIMALITY) {
                eligible = true;
                violation = rc;
            }
        } else if (st == BasisStatus::FREE_ZERO) {
            if (std::abs(rc) > EPS_OPTIMALITY) {
                eligible = true;
                violation = std::abs(rc);
            }
        }

        if (!eligible) continue;

        if (strategy_ == PricingStrategy::BLAND) {
            return j; // Bland's rule: take first eligible index
        }

        Real score = 0.0;
        if (strategy_ == PricingStrategy::DANTZIG) {
            score = violation;
        } else {
            // DEVEX or STEEPEST_EDGE
            Real w = (j < static_cast<Index>(weights_.size())) ? weights_[j] : 1.0;
            w = std::max(w, 1e-4);
            score = (violation * violation) / w;
        }

        if (score > best_score) {
            best_score = score;
            best_j = j;
        }
    }

    return best_j;
}

void Pricing::update(
    Index entering_var,
    Index leaving_row,
    const std::vector<Real>& alpha
) {
    iterations_since_reset_++;
    Index reset_thresh = std::max<Index>(100, total_vars_ / 5);
    if (iterations_since_reset_ > reset_thresh) {
        reset_weights();
        return;
    }

    if (strategy_ != PricingStrategy::DEVEX && strategy_ != PricingStrategy::STEEPEST_EDGE) {
        return;
    }

    if (leaving_row < 0 || leaving_row >= static_cast<Index>(alpha.size())) {
        return;
    }

    Real pivot_val = alpha[leaving_row];
    if (std::abs(pivot_val) < EPS_PIVOT) return;

    // Devex weight update
    Real pivot_sq = pivot_val * pivot_val;
    for (Index j = 0; j < total_vars_; ++j) {
        Real alpha_j = (j < static_cast<Index>(alpha.size())) ? alpha[j] : 0.0;
        Real candidate = 1.0 + (alpha_j * alpha_j) / pivot_sq;
        weights_[j] = std::max(0.01, std::max(candidate, 0.5 * weights_[j]));
    }

    if (entering_var >= 0 && entering_var < total_vars_) {
        weights_[entering_var] = 1.0;
    }
}

} // namespace suplex
