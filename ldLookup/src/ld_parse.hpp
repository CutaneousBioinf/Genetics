#ifndef _LDLOOKUP_LD_PARSE_HPP_
#define _LDLOOKUP_LD_PARSE_HPP_

#include <algorithm> // std::max
#include <iostream> // std::istream
#include <string>
#include <string_view>
#include <vector>

struct LDPair {
    std::string index_snp;
    std::string ld_snp;
    double index_maf;
    double r_squared;
};

struct LDPairRequirements {
    char column_separator;
    size_t index_snp_column;
    size_t ld_snp_column;
    size_t index_maf_column;
    size_t r_squared_column;
    size_t max_index_snp_size;
    double min_r_squared;
};

inline bool parse_ld_pair(
    const std::string &line,
    LDPair &ld_pair,
    const LDPairRequirements &reqs
) {
    size_t last_column = std::max({
        reqs.index_snp_column,
        reqs.ld_snp_column,
        reqs.index_maf_column,
        reqs.r_squared_column
    });

    size_t start(line.find_first_not_of(reqs.column_separator));
    size_t column(0);
    size_t end;
    while (start != std::string::npos && column <= last_column) {
        end = line.find(reqs.column_separator, start);
        size_t column_size(end - start);

        if (column == reqs.index_snp_column) {
            ld_pair.index_snp = line.substr(start, column_size);
        } else if (column == reqs.ld_snp_column) {
            ld_pair.ld_snp = line.substr(start, column_size);
        } else if (column == reqs.index_maf_column) {
            try {
                ld_pair.index_maf = std::stod(line.substr(start, column_size));
            } catch (std::invalid_argument &ignore) {
                return false;
            }
        } else if (column == reqs.r_squared_column) {
            try {
                ld_pair.r_squared = std::stod(line.substr(start, column_size));
            } catch (std::invalid_argument &ignore) {
                return false;
            }
        }

        start = line.find_first_not_of(reqs.column_separator, end);
        column++;
    }

    return column > last_column;
}

inline bool validate_ld_pair(
    const LDPair &ld_pair,
    const LDPairRequirements &reqs
) {
    return ld_pair.index_snp.size() <= reqs.max_index_snp_size
        && reqs.min_r_squared <= ld_pair.r_squared
        && ld_pair.r_squared <= 1;
}

#endif