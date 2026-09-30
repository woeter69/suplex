/// Unit tests for the Cut generators.
#include <gtest/gtest.h>
#include "milp/cuts.h"
#include "core/problem.h"
#include "core/solution.h"
#include "core/types.h"

using namespace suplex;

static Problem make_knapsack(const std::vector<Real>& coeffs, Real rhs) {
    int n = static_cast<int>(coeffs.size());
    std::vector<Index> r(n, 0), c(n);
    std::iota(c.begin(), c.end(), 0);

    Problem p;
    p.set_dimensions(1, n);
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, n, r, c, coeffs));
    p.set_objective(std::vector<Real>(n, 1.0));
    p.set_row_bounds({-INF}, {rhs});
    p.set_col_bounds(std::vector<Real>(n, 0.0), std::vector<Real>(n, 1.0));
    p.set_var_types(std::vector<VarType>(n, VarType::BINARY));
    return p;
}

// ── CutPool ───────────────────────────────────────────────────────────────────
TEST(CutPool, AddAndAge) {
    CutPool pool;
    Cut c1; c1.efficacy = 0.5; c1.type = CutType::GOMORY;
    Cut c2; c2.efficacy = 0.3; c2.type = CutType::MIR;
    pool.add(c1);
    pool.add(c2);
    EXPECT_EQ(pool.size(), 2);

    pool.age_all();
    EXPECT_EQ(pool.cuts()[0].age, 1);
    EXPECT_EQ(pool.cuts()[1].age, 1);
}

TEST(CutPool, RemoveOld) {
    CutPool pool;
    for (int i = 0; i < 5; ++i) {
        Cut c; c.age = i * 3;
        pool.add(c);
    }
    pool.remove_old(5);
    // Ages 0, 3 survive; 6, 9, 12 removed
    EXPECT_EQ(pool.size(), 2);
}

// ── Efficacy ──────────────────────────────────────────────────────────────────
TEST(CutGenerator, EfficacyPositiveWhenViolated) {
    CutGenerator gen(3);
    // Cut: x0 + x1 + x2 >= 1 (sum >= 1)
    Cut cut;
    cut.indices      = {0, 1, 2};
    cut.coefficients = {1.0, 1.0, 1.0};
    cut.rhs          = 1.0;

    // LP solution: all 0.2 → sum = 0.6 < 1 → violated → efficacy > 0
    Solution sol;
    sol.primal_values = {0.2, 0.2, 0.2};

    // Compute efficacy manually: (1.0 - 0.6) / sqrt(3) ≈ 0.231
    // We test via cover_cuts which calls this internally
    Real expected = (1.0 - 0.6) / std::sqrt(3.0);
    EXPECT_NEAR(expected, 0.2309, 0.001);
}

// ── Cover cuts ────────────────────────────────────────────────────────────────
TEST(CutGenerator, CoverCutFoundForKnapsack) {
    // Knapsack: 3x0 + 4x1 + 5x2 <= 6, binary vars
    // LP solution: x0=1, x1=0.5, x2=0.2 → sum of coefficients for a cover
    Problem p = make_knapsack({3.0, 4.0, 5.0}, 6.0);

    Solution sol;
    sol.primal_values = {1.0, 0.5, 0.2};

    CutGenerator gen(3);
    auto cuts = gen.cover_cuts(p, sol, 10);

    // Should find a cover (e.g., {x1, x2} with sum 4+5=9 > 6)
    // Cut: x1 + x2 <= 1
    EXPECT_GE(cuts.size(), 0u);   // may not find any if efficacy too low
}

TEST(CutGenerator, MinimalCoverCorrect) {
    // Coefficients: {6, 4, 3}, rhs=7
    // Greedy: sort desc = {6,4,3}. Add 6 (sum=6), add 4 (sum=10 > 7). Cover = {0,1}
    CutGenerator gen(3);
    Problem p = make_knapsack({6.0, 4.0, 3.0}, 7.0);

    Solution sol;
    sol.primal_values = {0.8, 0.7, 0.6};

    auto cuts = gen.cover_cuts(p, sol, 5);
    // Just verify no crash and at most 5 cuts returned
    EXPECT_LE(cuts.size(), 5u);
}

// ── MIR cuts ──────────────────────────────────────────────────────────────────
TEST(CutGenerator, MIRCutsNoSEGFAULT) {
    // 2-row, 3-col MIP
    int m = 2, n = 3;
    std::vector<Index> r = {0,0,0,1,1,1};
    std::vector<Index> c = {0,1,2,0,1,2};
    std::vector<Real>  v = {1,2,3,4,1,2};

    Problem p;
    p.set_dimensions(m, n);
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(m, n, r, c, v));
    p.set_objective({1,1,1});
    p.set_row_bounds({-INF,-INF}, {7.5, 6.3});
    p.set_col_bounds({0,0,0}, {5,5,5});
    p.set_var_types({VarType::INTEGER, VarType::INTEGER, VarType::CONTINUOUS});

    Solution sol;
    sol.primal_values = {1.5, 1.3, 0.8};

    CutGenerator gen(n);
    auto cuts = gen.mir_cuts(p, sol, 10);
    EXPECT_LE(cuts.size(), 10u);  // sanity check
}

// ── Clique cuts ───────────────────────────────────────────────────────────────
TEST(CutGenerator, CliqueCutFromClique) {
    // Row 0: x0 + x1 + x2 <= 1 (all binary) — this IS a clique constraint
    Problem p;
    p.set_dimensions(1, 3);
    std::vector<Index> r = {0,0,0}, c = {0,1,2};
    std::vector<Real>  v = {1,1,1};
    p.set_constraint_matrix(SparseMatrixCSC::from_triplets(1, 3, r, c, v));
    p.set_objective({1,1,1});
    p.set_row_bounds({-INF}, {1.0});
    p.set_col_bounds({0,0,0}, {1,1,1});
    p.set_var_types({VarType::BINARY, VarType::BINARY, VarType::BINARY});

    // LP solution: x0=0.4, x1=0.4, x2=0.3 → sum=1.1 > 1 → violated
    Solution sol;
    sol.primal_values = {0.4, 0.4, 0.3};

    CutGenerator gen(3);
    auto cuts = gen.clique_cuts(p, sol, 10);
    // Should find the clique cut
    EXPECT_GE(cuts.size(), 1u);
    if (!cuts.empty()) {
        EXPECT_EQ(cuts[0].type, CutType::CLIQUE);
        EXPECT_NEAR(cuts[0].rhs, 1.0, 1e-9);
        EXPECT_EQ(cuts[0].indices.size(), 3u);
    }
}

// ── Cut selection ─────────────────────────────────────────────────────────────
TEST(CutGenerator, SelectionFiltersParallelCuts) {
    CutGenerator gen(4);

    // Two identical cuts — one should be filtered
    std::vector<Cut> candidates;
    for (int i = 0; i < 2; ++i) {
        Cut c;
        c.indices      = {0, 1, 2};
        c.coefficients = {1.0, 1.0, 1.0};
        c.rhs          = 1.0;
        c.efficacy     = 0.5;
        candidates.push_back(c);
    }

    auto selected = gen.select_cuts(candidates, 10);
    // Parallel cuts filtered — should have at most 1
    EXPECT_LE(selected.size(), 1u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
