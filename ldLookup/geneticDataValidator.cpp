#include <stdexcept>

#include "geneticDataValidator.hpp"

bool GeneticDataValidator::validate(const std::string& line) {
    size_t end = 0;
    size_t start = line.find_first_not_of(delimiter, end);
    size_t col = 0;
    // Being lazy is fast, so process as few columns as possible.
    while (start != std::string::npos && col <= last_important_column) {
        end = line.find(delimiter, start);

        if (col == snp_a_index) {
            data.snp_a = line.substr(start, end - start);
            // Overlong snp_a strings are invalid.
            if (data.snp_a.size() > max_key_size) {
                return false;
            }
        } else if (col == snp_b_index) {
            data.snp_b = line.substr(start, end - start);
        } else if (col == maf_index) {
            // Non-double maf values are invalid.
            try {
                data.maf = std::stod(line.substr(start, end - start));
            } catch (std::invalid_argument& e) {
                return false;
            }
        } else if (col == r2_index) {
            // Non-double r2 values or r2 values below min_r2
            // are invalid.
            try {
                data.r2 = std::stod(line.substr(start, end - start));
                if (data.r2 < min_r2) {
                    return false;
                }
            } catch (std::invalid_argument& e) {
                return false;
            }
        }

        start = line.find_first_not_of(delimiter, end);
        col++;
    }

    // If `line` does not contain columns for every required
    // value, it is invalid.
    return col > last_important_column;
}
