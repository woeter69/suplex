#include "lp_parser.h"

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

enum class LPSection { NONE, OBJECTIVE, CONSTRAINTS, BOUNDS, GENERALS, BINARIES };

std::string trim(const std::string& value) {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

bool parse_number(const std::string& text, Real& value) {
    char* end = nullptr;
    errno = 0;
    value = std::strtod(text.c_str(), &end);
    return errno != ERANGE && end != text.c_str() && *end == '\0';
}

std::vector<std::string> words(const std::string& text) {
    std::istringstream input(text);
    std::vector<std::string> result;
    for (std::string word; input >> word;) result.push_back(std::move(word));
    return result;
}

std::vector<std::string> bound_words(const std::string& text) {
    std::string spaced;
    for (std::size_t i = 0; i < text.size(); ++i) {
        if ((text[i] == '<' || text[i] == '>') && i + 1 < text.size() && text[i + 1] == '=') {
            spaced += ' '; spaced += text[i]; spaced += "= "; ++i;
        } else if (text[i] == '=') {
            spaced += " = ";
        } else spaced += text[i];
    }
    return words(spaced);
}

bool parse_linear_expression(const std::string& text,
                             std::vector<std::pair<std::string, Real>>& terms,
                             std::string& error) {
    std::size_t pos = 0;
    Real sign = 1.0;
    while (pos < text.size()) {
        while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
        if (pos == text.size()) break;
        if (text[pos] == '+') { sign = 1.0; ++pos; continue; }
        if (text[pos] == '-') { sign = -1.0; ++pos; continue; }

        Real coefficient = 1.0;
        bool has_number = false;
        if (std::isdigit(static_cast<unsigned char>(text[pos])) || text[pos] == '.') {
            char* end = nullptr;
            coefficient = std::strtod(text.c_str() + pos, &end);
            if (end == text.c_str() + pos) { error = "invalid coefficient"; return false; }
            pos = static_cast<std::size_t>(end - text.c_str());
            has_number = true;
            while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
            if (pos < text.size() && text[pos] == '*') {
                ++pos;
                while (pos < text.size() && std::isspace(static_cast<unsigned char>(text[pos]))) ++pos;
            }
        }
        if (pos >= text.size() || !(std::isalpha(static_cast<unsigned char>(text[pos])) || text[pos] == '_')) {
            if (has_number) { // LP objective constants do not affect the optimizer.
                sign = 1.0;
                continue;
            }
            error = "expected a variable name";
            return false;
        }
        const std::size_t begin = pos++;
        while (pos < text.size()) {
            const unsigned char c = static_cast<unsigned char>(text[pos]);
            if (!(std::isalnum(c) || text[pos] == '_' || text[pos] == '.' || text[pos] == '$')) break;
            ++pos;
        }
        terms.emplace_back(text.substr(begin, pos - begin), sign * coefficient);
        sign = 1.0;
    }
    return true;
}

struct ParsedRow {
    std::string name;
    std::vector<std::pair<std::string, Real>> terms;
    char relation = 'E';
    Real rhs = 0.0;
};

bool split_constraint(const std::string& text, ParsedRow& row, std::string& error) {
    std::string body = trim(text);
    const auto colon = body.find(':');
    if (colon != std::string::npos) {
        row.name = trim(body.substr(0, colon));
        body = trim(body.substr(colon + 1));
    }
    std::size_t relation_pos = body.find("<=");
    std::size_t relation_len = 2;
    if (relation_pos != std::string::npos) row.relation = 'L';
    else {
        relation_pos = body.find(">=");
        if (relation_pos != std::string::npos) row.relation = 'G';
        else {
            relation_pos = body.find('=');
            relation_len = 1;
            row.relation = 'E';
        }
    }
    if (relation_pos == std::string::npos) { error = "constraint has no <=, >=, or = relation"; return false; }
    if (!parse_number(trim(body.substr(relation_pos + relation_len)), row.rhs)) {
        error = "constraint right-hand side is not numeric";
        return false;
    }
    return parse_linear_expression(body.substr(0, relation_pos), row.terms, error);
}

} // namespace

bool LPReader::read(const std::string& filename, Problem& problem) {
    error_.clear();
    col_names_.clear();
    row_names_.clear();
    std::ifstream input(filename);
    if (!input) { error_ = "cannot open LP file: " + filename; return false; }

    LPSection section = LPSection::NONE;
    ObjectiveSense sense = ObjectiveSense::MINIMIZE;
    bool saw_objective = false, saw_end = false;
    std::string objective_text;
    std::vector<ParsedRow> rows;
    std::vector<std::string> bound_lines, general_names, binary_names;
    std::string line;
    std::size_t line_number = 0;
    auto fail = [&](const std::string& message) {
        error_ = "LP line " + std::to_string(line_number) + ": " + message;
        return false;
    };

    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::string clean = trim(line);
        if (clean.empty() || clean.front() == '\\') continue;
        const std::string heading = lower(clean);
        if (heading == "minimize" || heading == "minimum" || heading == "min") {
            section = LPSection::OBJECTIVE; sense = ObjectiveSense::MINIMIZE; saw_objective = true; continue;
        }
        if (heading == "maximize" || heading == "maximum" || heading == "max") {
            section = LPSection::OBJECTIVE; sense = ObjectiveSense::MAXIMIZE; saw_objective = true; continue;
        }
        if (heading == "subject to" || heading == "such that" || heading == "st" || heading == "s.t.") {
            section = LPSection::CONSTRAINTS; continue;
        }
        if (heading == "bounds") { section = LPSection::BOUNDS; continue; }
        if (heading == "generals" || heading == "general" || heading == "integers" || heading == "integer") {
            section = LPSection::GENERALS; continue;
        }
        if (heading == "binaries" || heading == "binary") { section = LPSection::BINARIES; continue; }
        if (heading == "end") { saw_end = true; break; }

        if (section == LPSection::OBJECTIVE) {
            const auto colon = clean.find(':');
            if (objective_text.empty() && colon != std::string::npos) clean = trim(clean.substr(colon + 1));
            objective_text += " " + clean;
        } else if (section == LPSection::CONSTRAINTS) {
            ParsedRow row;
            std::string diagnostic;
            if (!split_constraint(clean, row, diagnostic)) return fail(diagnostic);
            if (row.name.empty()) row.name = "c" + std::to_string(rows.size() + 1);
            rows.push_back(std::move(row));
        } else if (section == LPSection::BOUNDS) bound_lines.push_back(clean);
        else if (section == LPSection::GENERALS) {
            const auto names = words(clean); general_names.insert(general_names.end(), names.begin(), names.end());
        } else if (section == LPSection::BINARIES) {
            const auto names = words(clean); binary_names.insert(binary_names.end(), names.begin(), names.end());
        } else return fail("content appears before a recognized section");
    }
    if (!saw_objective) return fail("missing Minimize or Maximize section");
    if (!saw_end) return fail("missing End marker");

    std::vector<std::pair<std::string, Real>> objective_terms;
    std::string diagnostic;
    if (!parse_linear_expression(objective_text, objective_terms, diagnostic)) return fail(diagnostic);
    std::unordered_map<std::string, Index> col_index;
    auto ensure_col = [&](const std::string& name) {
        const auto found = col_index.find(name);
        if (found != col_index.end()) return found->second;
        const Index index = static_cast<Index>(col_names_.size());
        col_index.emplace(name, index);
        col_names_.push_back(name);
        return index;
    };
    for (const auto& [name, value] : objective_terms) { (void)value; ensure_col(name); }
    for (const ParsedRow& row : rows) for (const auto& [name, value] : row.terms) { (void)value; ensure_col(name); }
    for (const std::string& name : general_names) ensure_col(name);
    for (const std::string& name : binary_names) ensure_col(name);

    // Bounds may introduce columns, so identify their variable token before sizing arrays.
    for (const std::string& bound : bound_lines) {
        for (const std::string& token : bound_words(bound)) {
            if (!token.empty() && (std::isalpha(static_cast<unsigned char>(token[0])) || token[0] == '_') && lower(token) != "free") {
                ensure_col(token);
            }
        }
    }

    const Index ncols = static_cast<Index>(col_names_.size());
    const Index nrows = static_cast<Index>(rows.size());
    std::vector<Real> objective(static_cast<std::size_t>(ncols), 0.0);
    std::vector<Real> col_lower(static_cast<std::size_t>(ncols), 0.0);
    std::vector<Real> col_upper(static_cast<std::size_t>(ncols), INF);
    std::vector<VarType> var_types(static_cast<std::size_t>(ncols), VarType::CONTINUOUS);
    std::vector<Real> row_lower(static_cast<std::size_t>(nrows), -INF);
    std::vector<Real> row_upper(static_cast<std::size_t>(nrows), INF);
    std::vector<Index> triplet_rows, triplet_cols;
    std::vector<Real> triplet_values;
    for (const auto& [name, value] : objective_terms) objective[static_cast<std::size_t>(col_index.at(name))] += value;
    for (Index i = 0; i < nrows; ++i) {
        const ParsedRow& row = rows[static_cast<std::size_t>(i)];
        row_names_.push_back(row.name);
        if (row.relation == 'L') row_upper[static_cast<std::size_t>(i)] = row.rhs;
        else if (row.relation == 'G') row_lower[static_cast<std::size_t>(i)] = row.rhs;
        else row_lower[static_cast<std::size_t>(i)] = row_upper[static_cast<std::size_t>(i)] = row.rhs;
        for (const auto& [name, value] : row.terms) {
            triplet_rows.push_back(i); triplet_cols.push_back(col_index.at(name)); triplet_values.push_back(value);
        }
    }

    for (const std::string& bound : bound_lines) {
        const auto token = bound_words(bound);
        if (token.size() == 2 && lower(token[1]) == "free") {
            const auto j = static_cast<std::size_t>(col_index.at(token[0]));
            col_lower[j] = -INF; col_upper[j] = INF; continue;
        }
        Real a = 0.0, b = 0.0;
        if (token.size() == 5 && token[1] == "<=" && token[3] == "<=" &&
            parse_number(token[0], a) && parse_number(token[4], b)) {
            const auto j = static_cast<std::size_t>(col_index.at(token[2])); col_lower[j] = a; col_upper[j] = b; continue;
        }
        if (token.size() == 3 && token[1] == "=" && parse_number(token[2], a)) {
            const auto j = static_cast<std::size_t>(col_index.at(token[0])); col_lower[j] = col_upper[j] = a; continue;
        }
        if (token.size() == 3 && (token[1] == "<=" || token[1] == ">=") && parse_number(token[2], a)) {
            const auto j = static_cast<std::size_t>(col_index.at(token[0]));
            if (token[1] == "<=") col_upper[j] = a; else col_lower[j] = a; continue;
        }
        if (token.size() == 3 && (token[1] == "<=" || token[1] == ">=") && parse_number(token[0], a)) {
            const auto j = static_cast<std::size_t>(col_index.at(token[2]));
            if (token[1] == "<=") col_lower[j] = a; else col_upper[j] = a; continue;
        }
        return fail("unsupported bound: " + bound);
    }
    for (const std::string& name : general_names) var_types[static_cast<std::size_t>(col_index.at(name))] = VarType::INTEGER;
    for (const std::string& name : binary_names) {
        const auto j = static_cast<std::size_t>(col_index.at(name));
        var_types[j] = VarType::BINARY; col_lower[j] = 0.0; col_upper[j] = 1.0;
    }
    for (std::size_t j = 0; j < col_lower.size(); ++j) if (col_lower[j] > col_upper[j]) return fail("inconsistent bounds for " + col_names_[j]);

    Problem parsed(nrows, ncols);
    parsed.set_objective_sense(sense);
    parsed.set_constraint_matrix(SparseMatrixCSC::from_triplets(nrows, ncols, triplet_rows, triplet_cols, triplet_values));
    parsed.set_objective(std::move(objective));
    parsed.set_col_bounds(std::move(col_lower), std::move(col_upper));
    parsed.set_row_bounds(std::move(row_lower), std::move(row_upper));
    parsed.set_var_types(std::move(var_types));
    if (!parsed.validate()) return fail("parsed problem is dimensionally invalid");
    problem = std::move(parsed);
    return true;
}

} // namespace suplex
