#pragma once
#include "../core/types.h"
#include "../core/problem.h"
#include "../core/solution.h"
#include <vector>
#include <string>

namespace suplex {

// ── Cut types ─────────────────────────────────────────────────────────────────
enum class CutType : uint8_t {
    GOMORY,    ///< Gomory mixed-integer cut
    MIR,       ///< Mixed-integer rounding cut
    CLIQUE,    ///< Clique cut from conflict graph
    COVER,     ///< Knapsack cover cut
    FLOW_COVER ///< Flow cover cut
};

// ── A single cutting plane ────────────────────────────────────────────────────
/// Represents: sum(coefficients[k] * x[indices[k]]) >= rhs
struct Cut {
    std::vector<Index> indices;
    std::vector<Real>  coefficients;
    Real               rhs       = 0.0;
    CutType            type      = CutType::GOMORY;
    Real               efficacy  = 0.0;  ///< Normalised distance from LP solution
    int                age       = 0;    ///< Rounds since last active in basis
    std::string        name;
};

// ── Cut pool ──────────────────────────────────────────────────────────────────
class CutPool {
public:
    void   add(Cut cut);
    void   age_all();                      ///< Increment age of all cuts
    void   remove_old(int max_age = 10);   ///< Remove stale cuts
    Index  size() const { return static_cast<Index>(cuts_.size()); }
    const  std::vector<Cut>& cuts() const { return cuts_; }

private:
    std::vector<Cut> cuts_;
};

// ── Cut generators ────────────────────────────────────────────────────────────
class CutGenerator {
public:
    explicit CutGenerator(Index num_cols);

    /// Generate Gomory cuts from optimal simplex basis.
    /// Requires the basis (from LPSolver::get_basis()) and tableau access.
    std::vector<Cut> gomory_cuts(const Problem&  p,
                                 const Solution& lp_sol,
                                 const std::vector<Index>& basis,
                                 int max_cuts = 50);

    /// Generate Mixed-Integer Rounding cuts.
    std::vector<Cut> mir_cuts   (const Problem&  p,
                                 const Solution& lp_sol,
                                 int max_cuts = 50);

    /// Generate knapsack cover cuts.
    std::vector<Cut> cover_cuts (const Problem&  p,
                                 const Solution& lp_sol,
                                 int max_cuts = 30);

    /// Generate clique cuts from binary conflict graph.
    std::vector<Cut> clique_cuts(const Problem&  p,
                                 const Solution& lp_sol,
                                 int max_cuts = 20);

    /// Select the best subset from a candidate pool.
    /// Filters by efficacy, orthogonality, and sparsity.
    std::vector<Cut> select_cuts(std::vector<Cut>& candidates,
                                 int max_to_add = 100);

    void set_min_efficacy(Real e)    { min_efficacy_ = e; }

private:
    Index num_cols_;
    Real  min_efficacy_   = 0.01;

    // ── Helpers ───────────────────────────────────────────────────────────
    Real  compute_efficacy(const Cut& cut, const Solution& lp_sol) const;
    Real  compute_orthogonality(const Cut& a, const Cut& b) const;

    // Gomory internals
    bool  extract_tableau_row(Index basic_pos,
                              const std::vector<Index>& basis,
                              const Problem& p,
                              std::vector<Real>& row_coeffs,
                              Real& rhs) const;

    // Cover cut internals
    bool  find_minimal_cover(const std::vector<Index>& knapsack_vars,
                             const std::vector<Real>&  knapsack_coeffs,
                             Real rhs,
                             std::vector<Index>& cover) const;
    Real  lifting_coefficient(const std::vector<Index>& cover,
                              const std::vector<Real>&  coeffs,
                              Index new_var_idx, Real new_var_coeff,
                              Real knapsack_rhs) const;

    // MIR internals
    Cut   mir_from_row(const std::vector<Real>& row_coeffs,
                       const std::vector<Index>& col_map,
                       Real rhs,
                       const Problem& p,
                       const Solution& lp_sol) const;
};

} // namespace suplex
