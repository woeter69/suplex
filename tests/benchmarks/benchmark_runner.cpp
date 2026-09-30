#include "src/api/suplex_solver.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct BenchmarkResult {
    std::string name;
    suplex::Index rows = 0, columns = 0, nonzeros = 0;
    suplex::SolverStatus status = suplex::SolverStatus::NOT_STARTED;
    double objective = 0.0;
    double reference = std::numeric_limits<double>::quiet_NaN();
    double gap = std::numeric_limits<double>::infinity();
    double seconds = 0.0;
    int64_t iterations = 0, nodes = 0;
    bool passed = false;
};

std::unordered_map<std::string, double> read_references(const std::filesystem::path& filename) {
    std::unordered_map<std::string, double> references;
    std::ifstream input(filename);
    for (std::string line; std::getline(input, line);) {
        if (line.empty() || line[0] == '#') continue;
        const auto comma = line.find(',');
        if (comma == std::string::npos || line.substr(0, comma) == "problem") continue;
        try { references[line.substr(0, comma)] = std::stod(line.substr(comma + 1)); }
        catch (const std::exception&) { std::cerr << "Skipping invalid reference row: " << line << '\n'; }
    }
    return references;
}

void write_csv(const std::filesystem::path& filename, const std::vector<BenchmarkResult>& results) {
    std::ofstream out(filename);
    out << "problem,rows,columns,nonzeros,status,objective,known_optimal,relative_gap,seconds,iterations,nodes,passed\n";
    out << std::setprecision(17);
    for (const auto& r : results) out << r.name << ',' << r.rows << ',' << r.columns << ',' << r.nonzeros
        << ',' << suplex::to_string(r.status) << ',' << r.objective << ',' << r.reference << ',' << r.gap
        << ',' << r.seconds << ',' << r.iterations << ',' << r.nodes << ',' << (r.passed ? 1 : 0) << '\n';
}

void write_json(const std::filesystem::path& filename, const std::vector<BenchmarkResult>& results) {
    std::ofstream out(filename);
    out << "[\n";
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        out << "  {\"problem\":\"" << r.name << "\",\"status\":\"" << suplex::to_string(r.status)
            << "\",\"rows\":" << r.rows << ",\"columns\":" << r.columns << ",\"nonzeros\":" << r.nonzeros
            << ",\"objective\":" << r.objective << ",\"relative_gap\":";
        if (std::isfinite(r.gap)) out << r.gap; else out << "null";
        out
            << ",\"seconds\":" << r.seconds << ",\"iterations\":" << r.iterations
            << ",\"nodes\":" << r.nodes << ",\"passed\":" << (r.passed ? "true" : "false") << '}';
        out << (i + 1 == results.size() ? "\n" : ",\n");
    }
    out << "]\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: suplex_benchmark <problem-directory> <reference.csv> [output-prefix]\n";
        return 2;
    }
    const std::filesystem::path directory = argv[1];
    const auto references = read_references(argv[2]);
    const std::filesystem::path prefix = argc == 4 ? argv[3] : "suplex-benchmark";
    std::vector<BenchmarkResult> results;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        std::string extension = entry.path().extension().string();
        if (extension != ".mps" && extension != ".lp" && extension != ".qps") continue;
        suplex::Suplex solver;
        const bool loaded = extension == ".lp" ? solver.read_lp(entry.path().string()) : solver.read_mps(entry.path().string());
        BenchmarkResult result;
        result.name = entry.path().stem().string();
        const auto reference = references.find(result.name);
        if (reference != references.end()) result.reference = reference->second;
        if (loaded) {
            result.rows = solver.problem().num_rows(); result.columns = solver.problem().num_cols();
            result.nonzeros = solver.problem().constraint_matrix().nnz;
            solver.set_time_limit(3600.0); solver.set_log_level(suplex::LogLevel::WARN);
            result.status = solver.solve();
            const auto& solution = solver.solution();
            result.objective = solution.objective_value; result.seconds = solution.solve_time_seconds;
            result.iterations = solution.num_iterations; result.nodes = solution.num_nodes;
            if (std::isfinite(result.reference)) {
                result.gap = std::abs(result.objective - result.reference) / std::max(1.0, std::abs(result.reference));
                result.passed = result.status == suplex::SolverStatus::OPTIMAL && result.gap <= 1e-6;
            }
        } else std::cerr << entry.path() << ": " << solver.last_error() << '\n';
        std::cout << std::left << std::setw(24) << result.name << std::setw(18) << suplex::to_string(result.status)
                  << std::right << std::setw(12) << result.seconds << " s\n";
        results.push_back(result);
    }
    write_csv(prefix.string() + ".csv", results);
    write_json(prefix.string() + ".json", results);
    const std::size_t passed = static_cast<std::size_t>(std::count_if(results.begin(), results.end(),
        [](const BenchmarkResult& result) { return result.passed; }));
    std::cout << "Passed " << passed << '/' << results.size() << " referenced problems\n";
    return 0;
}
