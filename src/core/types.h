#pragma once

#include <cstdint>
#include <limits>
#include <string>

namespace suplex {

// Floating point type used throughout the solver
using Real = double;
using Index = int32_t;

// Numerical Constants
constexpr Real INF = 1e30;
constexpr Real EPS_ZERO = 1e-12;
constexpr Real EPS_PIVOT = 1e-10;
constexpr Real EPS_FEASIBILITY = 1e-8;
constexpr Real EPS_OPTIMALITY = 1e-8;
constexpr Real EPS_INTEGER = 1e-6;
constexpr Index MAX_ITER = 1000000;
constexpr Real TIME_LIMIT = 3600.0;

enum class ObjectiveSense : uint8_t {
    MINIMIZE = 0,
    MAXIMIZE = 1
};

enum class VarType : uint8_t {
    CONTINUOUS = 0,
    INTEGER = 1,
    BINARY = 2
};

enum class BoundType : uint8_t {
    FREE = 0,
    LOWER = 1,
    UPPER = 2,
    BOTH = 3,
    FIXED = 4
};

enum class SolverStatus : uint8_t {
    NOT_STARTED = 0,
    OPTIMAL = 1,
    INFEASIBLE = 2,
    UNBOUNDED = 3,
    INF_OR_UNBD = 4,
    ITERATION_LIMIT = 5,
    TIME_LIMIT = 6,
    NUMERICAL_ERROR = 7,
    USER_INTERRUPT = 8
};

enum class LogLevel : uint8_t {
    TRACE = 0,
    DEBUG = 1,
    INFO = 2,
    WARN = 3,
    ERROR = 4,
    OFF = 5
};

inline std::string to_string(SolverStatus status) {
    switch (status) {
        case SolverStatus::NOT_STARTED: return "NOT_STARTED";
        case SolverStatus::OPTIMAL: return "OPTIMAL";
        case SolverStatus::INFEASIBLE: return "INFEASIBLE";
        case SolverStatus::UNBOUNDED: return "UNBOUNDED";
        case SolverStatus::INF_OR_UNBD: return "INF_OR_UNBD";
        case SolverStatus::ITERATION_LIMIT: return "ITERATION_LIMIT";
        case SolverStatus::TIME_LIMIT: return "TIME_LIMIT";
        case SolverStatus::NUMERICAL_ERROR: return "NUMERICAL_ERROR";
        case SolverStatus::USER_INTERRUPT: return "USER_INTERRUPT";
        default: return "UNKNOWN";
    }
}

} // namespace suplex
