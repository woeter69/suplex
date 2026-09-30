#pragma once

#include "src/core/solution.h"

#include <string>
#include <vector>

namespace suplex {

/** Writes a stable, human-readable Suplex solution file. */
class SolutionWriter {
public:
    bool write(const std::string& filename, const Solution& solution,
               const std::vector<std::string>& column_names = {},
               const std::vector<std::string>& row_names = {});
    std::string last_error() const { return error_; }

private:
    std::string error_;
};

} // namespace suplex
