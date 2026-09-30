#pragma once

#include "problem_reader.h"

namespace suplex {

/** Reader for fixed-column and whitespace-separated MPS files. */
class MPSReader final : public ProblemReader {
public:
    bool read(const std::string& filename, Problem& problem) override;
    std::string last_error() const override { return error_; }
    const std::vector<std::string>& column_names() const override { return col_names_; }
    const std::vector<std::string>& row_names() const override { return row_names_; }

private:
    std::string error_;
    std::vector<std::string> col_names_;
    std::vector<std::string> row_names_;
};

} // namespace suplex
