#include "mps_parser.h"

#include "src/core/sparse_matrix.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>

namespace suplex {
namespace {

enum class Section { NONE, ROWS, COLUMNS, RHS, RANGES, BOUNDS, OBJSENSE };

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

std::vector<std::string> fields(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> result;
    for (std::string token; input >> token;) result.push_back(std::move(token));
    return result;
}

bool number(const std::string& text, Real& value) {
    std::string normalized = text;
    std::replace(normalized.begin(), normalized.end(), 'D', 'E');
    std::replace(normalized.begin(), normalized.end(), 'd', 'e');
    char* end = nullptr;
    errno = 0;
    value = std::strtod(normalized.c_str(), &end);
    return errno != ERANGE && end != normalized.c_str() && *end == '\0';
}

std::string unquote(std::string text) {
    if (text.size() >= 2 && ((text.front() == '\'' && text.back() == '\'') ||
                            (text.front() == '"' && text.back() == '"'))) {
        return text.substr(1, text.size() - 2);
    }
    return text;
}

} // namespace

bool MPSReader::read(const std::string& filename, Problem& problem) {
    error_.clear();
    col_names_.clear();
    row_names_.clear();

    std::ifstream input(filename);
    if (!input) {
        error_ = "cannot open MPS file: " + filename;
        return false;
    }

    Section section = Section::NONE;
    bool saw_rows = false, saw_columns = false, saw_end = false;
    bool integer_block = false;
    ObjectiveSense sense = ObjectiveSense::MINIMIZE;
    std::string objective_row;
    std::unordered_map<std::string, Index> row_index;
    std::unordered_map<std::string, Index> col_index;
    std::vector<char> row_type;
    std::vector<Real> objective;
    std::vector<Real> col_lower, col_upper;
    std::vector<Real> row_lower, row_upper;
    std::vector<VarType> var_types;
    std::vector<Index> triplet_rows, triplet_cols;
    std::vector<Real> triplet_values;
    std::string selected_rhs, selected_range, selected_bounds;

    auto fail = [&](std::size_t line_number, const std::string& line,
                    const std::string& message) {
        error_ = "MPS line " + std::to_string(line_number) + ": " + message;
        if (!line.empty()) error_ += " [" + line + "]";
        return false;
    };
    auto add_column = [&](const std::string& name) {
        const Index index = static_cast<Index>(col_names_.size());
        col_index.emplace(name, index);
        col_names_.push_back(name);
        objective.push_back(0.0);
        col_lower.push_back(0.0);
        col_upper.push_back(INF);
        var_types.push_back(integer_block ? VarType::INTEGER : VarType::CONTINUOUS);
        return index;
    };

    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '*') continue;
        const std::vector<std::string> token = fields(line);
        if (token.empty()) continue;
        const std::string keyword = upper(token.front());

        if (keyword == "NAME") { section = Section::NONE; continue; }
        if (keyword == "OBJSENSE") {
            if (token.size() > 1) {
                const std::string value = upper(token[1]);
                if (value == "MAX" || value == "MAXIMIZE") sense = ObjectiveSense::MAXIMIZE;
                else if (value == "MIN" || value == "MINIMIZE") sense = ObjectiveSense::MINIMIZE;
                else return fail(line_number, line, "expected MIN or MAX after OBJSENSE");
                section = Section::NONE;
            } else section = Section::OBJSENSE;
            continue;
        }
        if (keyword == "ROWS") { section = Section::ROWS; saw_rows = true; continue; }
        if (keyword == "COLUMNS") { section = Section::COLUMNS; saw_columns = true; continue; }
        if (keyword == "RHS") { section = Section::RHS; continue; }
        if (keyword == "RANGES") { section = Section::RANGES; continue; }
        if (keyword == "BOUNDS") { section = Section::BOUNDS; continue; }
        if (keyword == "ENDATA") { saw_end = true; break; }

        if (section == Section::OBJSENSE) {
            if (keyword == "MAX" || keyword == "MAXIMIZE") sense = ObjectiveSense::MAXIMIZE;
            else if (keyword == "MIN" || keyword == "MINIMIZE") sense = ObjectiveSense::MINIMIZE;
            else return fail(line_number, line, "expected MIN or MAX after OBJSENSE");
            section = Section::NONE;
            continue;
        }

        if (section == Section::ROWS) {
            if (token.size() != 2 || token[0].size() != 1 ||
                std::string("NLGE").find(static_cast<char>(std::toupper(token[0][0]))) == std::string::npos) {
                return fail(line_number, line, "invalid ROWS entry");
            }
            const char type = static_cast<char>(std::toupper(token[0][0]));
            const std::string& name = token[1];
            if (type == 'N') {
                if (objective_row.empty()) objective_row = name;
                continue;
            }
            if (row_index.contains(name)) return fail(line_number, line, "duplicate row name " + name);
            const Index index = static_cast<Index>(row_names_.size());
            row_index.emplace(name, index);
            row_names_.push_back(name);
            row_type.push_back(type);
            row_lower.push_back(type == 'G' || type == 'E' ? 0.0 : -INF);
            row_upper.push_back(type == 'L' || type == 'E' ? 0.0 : INF);
            continue;
        }

        if (section == Section::COLUMNS) {
            if (token.size() >= 3 && upper(unquote(token[1])) == "MARKER") {
                const std::string marker = upper(unquote(token[2]));
                if (marker == "INTORG") integer_block = true;
                else if (marker == "INTEND") integer_block = false;
                else return fail(line_number, line, "unknown integer marker " + token[2]);
                continue;
            }
            if (token.size() != 3 && token.size() != 5) {
                return fail(line_number, line, "COLUMNS entry needs one or two row/value pairs");
            }
            Index col;
            const auto found_col = col_index.find(token[0]);
            if (found_col == col_index.end()) col = add_column(token[0]);
            else {
                col = found_col->second;
                if (integer_block) var_types[static_cast<std::size_t>(col)] = VarType::INTEGER;
            }
            for (std::size_t pos = 1; pos < token.size(); pos += 2) {
                Real value = 0.0;
                if (!number(token[pos + 1], value)) return fail(line_number, line, "invalid numeric value");
                if (token[pos] == objective_row) objective[static_cast<std::size_t>(col)] += value;
                else {
                    const auto row = row_index.find(token[pos]);
                    if (row == row_index.end()) return fail(line_number, line, "unknown row " + token[pos]);
                    triplet_rows.push_back(row->second);
                    triplet_cols.push_back(col);
                    triplet_values.push_back(value);
                }
            }
            continue;
        }

        if (section == Section::RHS || section == Section::RANGES) {
            if (token.size() != 3 && token.size() != 5) {
                return fail(line_number, line, "entry needs one or two row/value pairs");
            }
            std::string& selected = section == Section::RHS ? selected_rhs : selected_range;
            if (selected.empty()) selected = token[0];
            if (token[0] != selected) continue;
            for (std::size_t pos = 1; pos < token.size(); pos += 2) {
                Real value = 0.0;
                if (!number(token[pos + 1], value)) return fail(line_number, line, "invalid numeric value");
                if (token[pos] == objective_row) continue;
                const auto row = row_index.find(token[pos]);
                if (row == row_index.end()) return fail(line_number, line, "unknown row " + token[pos]);
                const std::size_t i = static_cast<std::size_t>(row->second);
                const char type = row_type[i];
                if (section == Section::RHS) {
                    if (type == 'L') row_upper[i] = value;
                    else if (type == 'G') row_lower[i] = value;
                    else row_lower[i] = row_upper[i] = value;
                } else {
                    const Real magnitude = std::abs(value);
                    if (type == 'L') row_lower[i] = row_upper[i] - magnitude;
                    else if (type == 'G') row_upper[i] = row_lower[i] + magnitude;
                    else if (value >= 0.0) row_upper[i] = row_lower[i] + magnitude;
                    else row_lower[i] = row_upper[i] - magnitude;
                }
            }
            continue;
        }

        if (section == Section::BOUNDS) {
            if (token.size() < 3 || token.size() > 4) return fail(line_number, line, "invalid BOUNDS entry");
            if (selected_bounds.empty()) selected_bounds = token[1];
            if (token[1] != selected_bounds) continue;
            const auto found = col_index.find(token[2]);
            if (found == col_index.end()) return fail(line_number, line, "unknown column " + token[2]);
            const std::size_t col = static_cast<std::size_t>(found->second);
            const std::string type = upper(token[0]);
            Real value = 0.0;
            const bool needs_value = type == "UP" || type == "LO" || type == "FX" || type == "LI" || type == "UI";
            if (needs_value && (token.size() != 4 || !number(token[3], value))) {
                return fail(line_number, line, "bound type " + type + " requires a numeric value");
            }
            if (type == "UP") col_upper[col] = value;
            else if (type == "LO") col_lower[col] = value;
            else if (type == "FX") col_lower[col] = col_upper[col] = value;
            else if (type == "FR") { col_lower[col] = -INF; col_upper[col] = INF; }
            else if (type == "MI") col_lower[col] = -INF;
            else if (type == "PL") col_upper[col] = INF;
            else if (type == "BV") { col_lower[col] = 0.0; col_upper[col] = 1.0; var_types[col] = VarType::BINARY; }
            else if (type == "LI") { col_lower[col] = value; var_types[col] = VarType::INTEGER; }
            else if (type == "UI") { col_upper[col] = value; var_types[col] = VarType::INTEGER; }
            else return fail(line_number, line, "unknown bound type " + type);
            continue;
        }

        return fail(line_number, line, "data outside a recognized section");
    }

    if (!saw_rows) return fail(line_number, "", "missing ROWS section");
    if (!saw_columns) return fail(line_number, "", "missing COLUMNS section");
    if (!saw_end) return fail(line_number, "", "missing ENDATA");
    for (std::size_t j = 0; j < col_lower.size(); ++j) {
        if (col_lower[j] > col_upper[j]) return fail(line_number, "", "inconsistent bounds for " + col_names_[j]);
    }

    const Index rows = static_cast<Index>(row_names_.size());
    const Index cols = static_cast<Index>(col_names_.size());
    Problem parsed(rows, cols);
    parsed.set_objective_sense(sense);
    parsed.set_constraint_matrix(SparseMatrixCSC::from_triplets(
        rows, cols, triplet_rows, triplet_cols, triplet_values));
    parsed.set_objective(std::move(objective));
    parsed.set_col_bounds(std::move(col_lower), std::move(col_upper));
    parsed.set_row_bounds(std::move(row_lower), std::move(row_upper));
    parsed.set_var_types(std::move(var_types));
    if (!parsed.validate()) return fail(line_number, "", "parsed problem is dimensionally invalid");
    problem = std::move(parsed);
    return true;
}

} // namespace suplex
