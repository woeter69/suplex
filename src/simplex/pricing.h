#pragma once

#include "src/core/types.h"
#include "basis.h"
#include <vector>

namespace suplex {

enum class PricingStrategy : uint8_t {
    DANTZIG = 0,
    DEVEX = 1,
    STEEPEST_EDGE = 2,
    BLAND = 3
};

inline const char* to_string(PricingStrategy s) {
    switch (s) {
        case PricingStrategy::DANTZIG: return "DANTZIG";
        case PricingStrategy::DEVEX: return "DEVEX";
        case PricingStrategy::STEEPEST_EDGE: return "STEEPEST_EDGE";
        case PricingStrategy::BLAND: return "BLAND";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Pricing engine supporting Dantzig, Devex, Steepest-Edge, and Bland's anti-cycling rule.
 */
class Pricing {
public:
    Pricing() = default;
    explicit Pricing(PricingStrategy strategy);

    void init(Index total_vars, Index num_rows, PricingStrategy strategy = PricingStrategy::DEVEX);

    /**
     * @brief Selects the entering variable index based on reduced costs and variable statuses.
     * @return Index of entering variable, or -1 if current solution is optimal.
     */
    Index select_entering(
        const std::vector<Real>& reduced_costs,
        const std::vector<BasisStatus>& var_status
    );

    /**
     * @brief Updates pricing weights after a pivot step.
     */
    void update(
        Index entering_var,
        Index leaving_row,
        const std::vector<Real>& alpha
    );

    void reset_weights();
    PricingStrategy strategy() const { return strategy_; }
    void set_strategy(PricingStrategy strat) { strategy_ = strat; }

private:
    Index total_vars_ = 0;
    Index num_rows_ = 0;
    PricingStrategy strategy_ = PricingStrategy::DEVEX;
    std::vector<Real> weights_;
    Index iterations_since_reset_ = 0;
};

} // namespace suplex
