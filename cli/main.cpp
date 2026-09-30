#include "src/api/suplex_solver.h"
#include "src/io/sol_writer.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <thread>

namespace {

struct Options {
    suplex::Algorithm algorithm = suplex::Algorithm::AUTO;
    suplex::LogLevel log_level = suplex::LogLevel::INFO;
    bool presolve = true;
    bool gpu = false;
    bool stats = false;
    double time_limit = suplex::TIME_LIMIT;
    double mip_gap = 1e-4;
    int threads = 0;
    std::string solution_file;
    std::string problem_file;
};

void usage(std::ostream& out) {
    out << "Usage: suplex [options] <problem_file>\n\n"
        << "Options:\n"
        << "  -h, --help              Show this help message\n"
        << "  -v, --version           Show version\n"
        << "  --algorithm <algo>      primal, dual, ipm, or auto\n"
        << "  --presolve <on|off>     Enable or disable presolve\n"
        << "  --time-limit <seconds>  Positive solve time limit\n"
        << "  --threads <n>           Worker count (0 means automatic)\n"
        << "  --gpu <on|off>          Request GPU acceleration\n"
        << "  --log-level <level>     trace, debug, info, warn, error, off\n"
        << "  --mip-gap <value>       Nonnegative relative MIP gap\n"
        << "  --solution <file>       Write a solution file\n"
        << "  --stats                 Print detailed statistics\n";
}

bool parse_switch(const std::string& value, bool& result) {
    if (value == "on") { result = true; return true; }
    if (value == "off") { result = false; return true; }
    return false;
}

bool parse_options(int argc, char** argv, Options& options, std::string& error, bool& early_exit) {
    early_exit = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") { usage(std::cout); early_exit = true; return true; }
        if (arg == "-v" || arg == "--version") { std::cout << "Suplex 0.1.0\n"; early_exit = true; return true; }
        if (arg == "--stats") { options.stats = true; continue; }
        if (!arg.empty() && arg[0] != '-') {
            if (!options.problem_file.empty()) { error = "more than one problem file was supplied"; return false; }
            options.problem_file = arg; continue;
        }
        if (i + 1 >= argc) { error = "missing value for " + arg; return false; }
        const std::string value = argv[++i];
        try {
            if (arg == "--algorithm") {
                if (value == "auto") options.algorithm = suplex::Algorithm::AUTO;
                else if (value == "primal") options.algorithm = suplex::Algorithm::PRIMAL_SIMPLEX;
                else if (value == "dual") options.algorithm = suplex::Algorithm::DUAL_SIMPLEX;
                else if (value == "ipm") options.algorithm = suplex::Algorithm::IPM;
                else { error = "unknown algorithm: " + value; return false; }
            } else if (arg == "--presolve") {
                if (!parse_switch(value, options.presolve)) { error = "--presolve expects on or off"; return false; }
            } else if (arg == "--gpu") {
                if (!parse_switch(value, options.gpu)) { error = "--gpu expects on or off"; return false; }
            } else if (arg == "--time-limit") {
                options.time_limit = std::stod(value);
                if (options.time_limit <= 0.0) { error = "time limit must be positive"; return false; }
            } else if (arg == "--threads") {
                options.threads = std::stoi(value);
                if (options.threads < 0) { error = "thread count must be nonnegative"; return false; }
            } else if (arg == "--mip-gap") {
                options.mip_gap = std::stod(value);
                if (options.mip_gap < 0.0) { error = "MIP gap must be nonnegative"; return false; }
            } else if (arg == "--solution") options.solution_file = value;
            else if (arg == "--log-level") {
                if (value == "trace") options.log_level = suplex::LogLevel::TRACE;
                else if (value == "debug") options.log_level = suplex::LogLevel::DEBUG;
                else if (value == "info") options.log_level = suplex::LogLevel::INFO;
                else if (value == "warn") options.log_level = suplex::LogLevel::WARN;
                else if (value == "error") options.log_level = suplex::LogLevel::ERROR;
                else if (value == "off") options.log_level = suplex::LogLevel::OFF;
                else { error = "unknown log level: " + value; return false; }
            } else { error = "unknown option: " + arg; return false; }
        } catch (const std::exception&) { error = "invalid value for " + arg + ": " + value; return false; }
    }
    if (options.problem_file.empty()) { error = "a .mps or .lp problem file is required"; return false; }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Options options;
    std::string error;
    bool early_exit = false;
    if (!parse_options(argc, argv, options, error, early_exit)) {
        std::cerr << "suplex: " << error << "\nTry 'suplex --help' for usage.\n";
        return 2;
    }
    if (early_exit) return 0;

    suplex::Suplex solver;
    std::string extension = std::filesystem::path(options.problem_file).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    bool loaded = false;
    if (extension == ".mps" || extension == ".qps") loaded = solver.read_mps(options.problem_file);
    else if (extension == ".lp") loaded = solver.read_lp(options.problem_file);
    else { std::cerr << "suplex: unsupported file extension: " << extension << '\n'; return 2; }
    if (!loaded) { std::cerr << "suplex: " << solver.last_error() << '\n'; return 2; }

    const suplex::Problem& problem = solver.problem();
    std::cout << "Problem: " << problem.num_rows() << " rows, " << problem.num_cols()
              << " columns, " << problem.constraint_matrix().nnz << " nonzeros, "
              << problem.num_integers() << " integer variables\n";
    solver.set_algorithm(options.algorithm);
    solver.set_presolve(options.presolve);
    solver.set_time_limit(options.time_limit);
    solver.set_threads(options.threads);
    solver.set_gpu(options.gpu);
    solver.set_log_level(options.log_level);
    solver.set_mip_gap(options.mip_gap);
    const suplex::SolverStatus status = solver.solve();
    const suplex::Solution& solution = solver.solution();
    std::cout << "Status: " << suplex::to_string(status) << '\n';
    if (solution.is_feasible()) std::cout << "Objective: " << solution.objective_value << '\n';
    if (!solver.last_error().empty()) std::cerr << "suplex: " << solver.last_error() << '\n';
    if (options.stats) {
        std::cout << "Iterations: " << solution.num_iterations << "\nSolve time: "
                  << solution.solve_time_seconds << " seconds\n";
        if (problem.is_mip()) std::cout << "Nodes: " << solution.num_nodes << "\nMIP gap: " << solution.mip_gap << '\n';
    }
    if (!options.solution_file.empty()) {
        suplex::SolutionWriter writer;
        if (!writer.write(options.solution_file, solution, solver.column_names(), solver.row_names())) {
            std::cerr << "suplex: " << writer.last_error() << '\n'; return 2;
        }
    }
    return status == suplex::SolverStatus::OPTIMAL ? 0 : 1;
}
