#pragma once
#include <cstdint>
#include <limits>

namespace suplex {

/// Floating point type used throughout the solver (64-bit for numerical stability)
using Real  = double;
/// Index type for matrix dimensions and variable/constraint indices
using Index = int32_t;

// ── Numerical constants ──────────────────────────────────────────────────────
constexpr Real  INF              = 1e30;   ///< Numerical infinity sentinel
constexpr Real  EPS_ZERO         = 1e-12;  ///< General zero tolerance
constexpr Real  EPS_PIVOT        = 1e-10;  ///< Minimum acceptable pivot element
constexpr Real  EPS_FEASIBILITY  = 1e-8;   ///< Primal/dual feasibility tolerance
constexpr Real  EPS_OPTIMALITY   = 1e-8;   ///< Reduced-cost optimality tolerance
constexpr Real  EPS_INTEGER      = 1e-6;   ///< Integrality tolerance
constexpr Index MAX_ITER         = 1000000; ///< Default iteration limit
constexpr Real  TIME_LIMIT_DEFAULT = 3600.0; ///< Default time limit (seconds)

// ── Enumerations ─────────────────────────────────────────────────────────────

enum class ObjectiveSense : uint8_t {
    MINIMIZE,
    MAXIMIZE
};

enum class VarType : uint8_t {
    CONTINUOUS,
    INTEGER,
    BINARY
};

enum class BoundType : uint8_t {
    FREE,   ///< No bounds
    LOWER,  ///< Only lower bound
    UPPER,  ///< Only upper bound
    BOTH,   ///< Both bounds finite
    FIXED   ///< Lower == Upper
};

enum class SolverStatus : uint8_t {
    NOT_STARTED,
    OPTIMAL,
    INFEASIBLE,
    UNBOUNDED,
    INF_OR_UNBD,
    ITERATION_LIMIT,
    TIME_LIMIT,
    NUMERICAL_ERROR,
    USER_INTERRUPT
};

enum class LogLevel : uint8_t {
    TRACE,
    DEBUG,
    INFO,
    WARN,
    ERROR,
    OFF
};

} // namespace suplex
