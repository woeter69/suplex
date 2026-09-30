#pragma once

#include "src/core/problem.h"

#include <string>
#include <vector>

namespace suplex {

/** Common interface implemented by optimization problem file readers. */
class ProblemReader {
public:
    virtual ~ProblemReader() = default;
    /** Parse filename into problem, replacing its previous contents. */
    virtual bool read(const std::string& filename, Problem& problem) = 0;
    /** Return the most recent diagnostic, including a line number where possible. */
    virtual std::string last_error() const = 0;
    /** Names in the same order as the parsed problem columns. */
    virtual const std::vector<std::string>& column_names() const = 0;
    /** Names in the same order as the parsed problem rows. */
    virtual const std::vector<std::string>& row_names() const = 0;
};

} // namespace suplex
