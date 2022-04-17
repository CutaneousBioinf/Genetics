#ifndef _LDLOOKUP_PARSE_VARIANTS_HPP_
#define _LDLOOKUP_PARSE_VARIANTS_HPP_

#include <string>
#include <vector>

/**
 *  TODO: Document!
 */
struct IndexVariantSummary {
	std::string variant_id;
	double maf;
	size_t n_surrogates;
};

/**
 *  TODO: Document!
 */
struct LDPair {
    std::string index_variant_id;
    std::string ld_variant_id;
    double index_variant_maf;
	double r2;

	bool is_in_ld(double r2_threshold) const;
};

/**
 *  TODO: Document!
 */
struct LDPairParser {
	/**
	 *  TODO: Document!
	 */
	char delimiter;
	size_t index_variant_id_column;
	size_t ld_variant_id_column;
	size_t index_variant_maf_column;
	size_t r2_column;
	double r2_threshold_for_ld;

	/**
	 *  TODO: Document!
	 */
	size_t get_max_column() const;

	/**
	 *  TODO: Document!
	 */
	bool parse_pair(const std::string& s, LDPair& pair) const;
};

/**
 *  TODO: Document!
 */
template <typename F1, typename F2, typename F3>
void iterate_ld_data(
    const std::string& ld_data_path,
    const LDPairParser& parser,
    F1 on_ld_pair,
    F2 on_new_index_variant,
    F3 on_invalid_line);

/**
 *  TODO: Document!
 */
template <typename F1>
void iterate_variants(
    const std::string& variants_file,
    const std::vector<std::string>& additional_variants,
    F1 on_variant);

/*************************************************/
/*************************************************/
/****             Implementations             ****/
/*************************************************/
/*************************************************/

#include <algorithm>  // std::max
#include <fstream>    // std::ifstream

#include "string_ops.hpp"

inline bool LDPair::is_in_ld(double r2_threshold) const {
	return r2 >= r2_threshold;
}

inline size_t LDPairParser::get_max_column() const {
	return std::max({
		index_variant_id_column,
		ld_variant_id_column,
		index_variant_maf_column,
		r2_column
    });
}

inline bool LDPairParser::parse_pair(
	const std::string& s,
	LDPair& pair) const {
	std::string_view sv(s);
    std::vector<std::string_view> vec(split(sv, delimiter));

	size_t max_col = get_max_column();
    if (vec.size() < max_col) {
        return false;
    }

	// Set MAF field.
    // MAF must be a double.
	try {
		std::string col_str(vec.at(index_variant_maf_column-1));
		pair.index_variant_maf = std::stod(col_str);
	} catch (std::invalid_argument& ignore) {
		return false;
	}

	// MAF must be in [0, 0.5].
	if (pair.index_variant_maf > 0.5 || pair.index_variant_maf < 0) {
		return false;
	}

	// Set r-squared field.
	// R2 must be a double.
	try {
		std::string col_str(vec.at(r2_column-1));
		pair.r2 = std::stod(col_str);
	} catch (std::invalid_argument& ignore) {
		return false;
	}

	// R2 must be in [0, 1].
	if (pair.r2 > 1 || pair.r2 < 0) {
		return false;
	}

	// Set index_variant_id and ld_variant_id. No validation is performed here.
    pair.index_variant_id = std::string(vec.at(index_variant_id_column-1));
	pair.ld_variant_id = std::string(vec.at(ld_variant_id_column-1));
    return true;
}

template <typename F1, typename F2, typename F3>
inline void iterate_ld_data(
    const std::string& ld_data_file,
    const LDPairParser& parser,
    F1 on_ld_pair,
    F2 on_index_variant_summary,
    F3 on_invalid_line) {
	// Open the LD data.
	std::ifstream ld_data(ld_data_file);
	if (!ld_data) {
		throw std::runtime_error("Failed to open '" + ld_data_file + "'");
	}
	
	// Initialize curr_summary with placeholder/junk values.
	IndexVariantSummary curr_summary{ "", 0.0, 0 };
	bool found_variant = false;
	std::string line;
	LDPair parsed_line;

	// Loop over lines of ld_data.
	while (std::getline(ld_data, line)) {
		// Parse each line into an LDPair.
		if (!parser.parse_pair(line, parsed_line)) {
			on_invalid_line(line);
			continue;
		}

		// Check if we found a new index variant.
		if (parsed_line.index_variant_id != curr_summary.variant_id) {
			// Report summary statistics on the previous index variant.
			// found_variant ensures such a variant exists: we don't
			// want to report our placeholder values when we encounter
			// our first valid variant.
			if (found_variant) {
				on_index_variant_summary(curr_summary);
			} else {
				found_variant = true;
			}

			// Set the current index variant to the new index variant.
			curr_summary.variant_id = parsed_line.index_variant_id;
			curr_summary.maf = parsed_line.index_variant_maf;
			curr_summary.n_surrogates = 0;
		}

		// Report pairs in LD.
		if (parsed_line.is_in_ld(parser.r2_threshold_for_ld)) {
			curr_summary.n_surrogates++;
			on_ld_pair(parsed_line);
		}
	}

	// Report summary statistics for the final index variant.
	if (found_variant) {
		on_index_variant_summary(curr_summary);
	}
}

template <typename F1>
inline void iterate_variants(
    const std::string& variants_file,
    const std::vector<std::string>& additional_variants,
    F1 on_variant) {
	if (variants_file.size()) {
		// Open variants_file.
		std::ifstream file(variants_file, std::ios_base::in);
		if (!file.is_open()) {
			throw std::runtime_error("Failed to open '" + variants_file + "'");
		}

		// For each line in variants_file, call on_variant.
		std::string line;
		while (std::getline(file, line)) {
			on_variant(line);
		}
	}

    // For each line in additional_variants, call on_variant.
    for (const std::string& var : additional_variants) {
        on_variant(var);
    }
}

#endif