#include "sol_writer.h"

#include "src/core/types.h"

#include <fstream>
#include <iomanip>

namespace suplex {

bool SolutionWriter::write(const std::string& filename, const Solution& solution,
                           const std::vector<std::string>& column_names,
                           const std::vector<std::string>& row_names) {
    error_.clear();
    std::ofstream output(filename);
    if (!output) { error_ = "cannot open solution file for writing: " + filename; return false; }
    output << "Suplex Solution File\n"
           << "Status: " << to_string(solution.status) << '\n'
           << "Objective: " << std::setprecision(17) << solution.objective_value << "\n\n"
           << "Variables:\n" << std::scientific << std::setprecision(8);
    for (std::size_t j = 0; j < solution.primal_values.size(); ++j) {
        const std::string name = j < column_names.size() ? column_names[j] : "x" + std::to_string(j + 1);
        output << "  " << name << "    " << solution.primal_values[j] << '\n';
    }
    output << "\nConstraints:\n";
    for (std::size_t i = 0; i < solution.row_activities.size(); ++i) {
        const std::string name = i < row_names.size() ? row_names[i] : "c" + std::to_string(i + 1);
        const Real dual = i < solution.dual_values.size() ? solution.dual_values[i] : 0.0;
        output << "  " << name << "    Activity: " << solution.row_activities[i]
               << "  Dual: " << dual << '\n';
    }
    output << "\nStatistics:\n"
           << "  Iterations: " << solution.num_iterations << '\n'
           << "  Solve Time: " << std::fixed << std::setprecision(6)
           << solution.solve_time_seconds << " seconds\n";
    if (solution.num_nodes > 0) {
        output << "  Nodes: " << solution.num_nodes << '\n'
               << "  Gap: " << std::setprecision(4) << solution.mip_gap * 100.0 << "%\n";
    }
    if (!output) { error_ = "failed while writing solution file: " + filename; return false; }
    return true;
}

} // namespace suplex
