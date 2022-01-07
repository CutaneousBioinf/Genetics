#ifndef LDLOOKUP_GENETICDATA_HPP
#define LDLOOKUP_GENETICDATA_HPP

#include <algorithm> // std::lower_bound, std::max
#include <iterator> // std::next, std::prev
#include <map> // std::map
#include <stdexcept> // std::invalid_argument
#include <string> // std::stod, std::string
#include <utility> // std::make_pair
#include <vector> // std::vector

#include "global.hpp"

/** 
 * A simple container class for linkage disequilibrium data.
 * 
 * Data Members:
 *  index_snp_id - Index SNP ID
 *  ld_snp_id - ID of SNP in LD with index SNP
 *  ld_score - LD score between SNP_A and SNP_B
 *  maf - Minor allele frequency of index SNP
**/
struct GeneticData {
    std::string index_snp_id;
    std::string ld_snp_id;
    double ld_score;
    double maf;
};

/* Parses valid lines of linkage disequilibrium (LD) data to GeneticData. */
class GeneticDataValidator {
	public:
        // Separates fields within lines of LD data.
		char delimiter;
        // Indices of columns containing fields with which to populate GeneticData.
        // Column indices are zero-indexed and increment when one or more delimiter
        // characters are encountered in a line of LD data.
		size_t index_snp_id_col;
		size_t ld_snp_id_col;
		size_t ld_score_col;
		size_t maf_col;
        // Lines of LD data with R2 values below this threshold are skipped.
        double threshold_ld_score;
        // Lines with index_snp_id fields exceeding this length are skipped.
        size_t index_snp_id_max_length;
        // Parsed, validated LD data
		GeneticData data;

		GeneticDataValidator(const char delimiter,
                             const size_t index_snp_id_col,
                             const size_t ld_snp_id_col,
                             const size_t ld_score_col,
                             const size_t maf_col,
                             const double threshold_ld_score,
                             const size_t index_snp_id_max_length) :
                        delimiter(delimiter),
                        index_snp_id_col(index_snp_id_col),
                        ld_snp_id_col(ld_snp_id_col),
                        ld_score_col(ld_score_col),
                        maf_col(maf_col),
                        threshold_ld_score(threshold_ld_score),
                        index_snp_id_max_length(index_snp_id_max_length),
                        last_important_column(std::max({ index_snp_id_col,
                                                         ld_snp_id_col,
                                                         ld_score_col,
                                                         maf_col }))
                        {}

        /** 
         * Attempts to parse one line of LD data. Returns true if
         * the line is valid.
         * 
         * If the line is valid, the data member will contain values parsed
         * from the line. If the line is invalid, the data member may contain
         * junk.
         **/
		bool validate(const std::string& line) {
            size_t end = 0;
            size_t start = line.find_first_not_of(delimiter, end);
            size_t col = 0;

            // We process columns as little as possible.
            while (start != std::string::npos && col <= last_important_column) {
                end = line.find(delimiter, start);

                if (col == index_snp_id_col) {
                    // index_snp_id_col must be below maximum length
                    data.index_snp_id = line.substr(start, end - start);
                    if (data.index_snp_id.size() > index_snp_id_max_length) {
                        return false;
                    }
                } else if (col == ld_snp_id_col) {
                    data.ld_snp_id = line.substr(start, end - start);
                } else {
                    try {
                        if (col == maf_col) {
                            // maf_col must be a double
                            data.maf = std::stod(line.substr(start, end - start));
                        } else if (col == ld_score_col) {
                            // ld_score_col must be a double greater than or
                            // equal to the threshold ld_score
                            data.ld_score = std::stod(line.substr(start, end - start));
                            if (data.ld_score < threshold_ld_score) {
                                return false;
                            }
                        }
                    } catch (std::invalid_argument& ignore) {
                        // Catches invalid maf and ld_score fields
                        return false;
                    }
                }

                start = line.find_first_not_of(delimiter, end);
                col++;
            }

            // The line must contain columns for every data field.
            return col > last_important_column;
        }

	private:
		size_t last_important_column;
};

/**
 * Chunks a one-dimensional frequency distribution into bins of roughly
 * uniform frequencies.
 * 
 * Parameters:
 *  hist: Map representing a histogram. Keys of the map represent the lower
 *      bound value of one bin of the histogram. Values of the map represent
 *      the number of observations in the range [key, next greater key).
 *      ex. The data 0,1,1,2,3,5,8 could be represented by the histograms
 *      { 0: 1, 1: 2, 2: 1, 3: 1, 5: 1, 8: 1 } or {0: 4, 3: 2, 6: 1}.
 *  n_bins: Number of uniform bins to create.
 * 
 * Returns:
 *  A histogram with n_bins bins and a roughly even number of items in each
 *  bin. This histogram represents the same frequency distribution as hist.
 * 
 * Example:
 *  bin_histogram({ 0: 1, 1: 2, 2: 1, 3: 1, 5: 1, 8: 1 }, 2) -> {0: 4, 3: 3}
 **/
template <class Key>
std::map<Key, size_t> bin_histogram(
    const std::map<Key, size_t>& hist,
    const size_t n_bins
) {
    if (!n_bins || hist.empty()) {
        return std::map<Key, size_t>();
    }

    // Find the total number of items in the histogram.
    size_t total = 0;
    for (const auto& [k, v] : hist) {
        total += v;
    }

    std::map<Key, size_t> bins;
    // Minimum number of items per uniform bin
    size_t bin_spacing = (total + (n_bins / 2)) / n_bins;
    // Number of items in the current uniform bin
    size_t bin_size = hist.begin()->second;
    // Lower bound of the current uniform bin
    Key last_cutpoint = hist.begin()->first;
    for (auto it = std::next(hist.begin()); it != hist.end(); it++) {
        if (bin_size >= bin_spacing) {
            bins.insert(std::make_pair(last_cutpoint, bin_size));
            bin_size = 0;
            last_cutpoint = it->first;
        }

        bin_size += it->second;
    }

    if (bin_size) {
        bins.insert(std::make_pair(last_cutpoint, bin_size));
    }

    return bins;
}

/* Extracts the keys of a map into a vector */
template <class Key, class Value>
std::vector<Key> get_map_keys(const std::map<Key, Value>& map) {
    std::vector<Key> ret;
    for (const auto& [k, v] : map) {
        ret.push_back(k);
    }
    return ret;
}

/* Gets last item of a vector sorted in ascending order that is <=`t`. */
template <class T>
typename std::vector<T>::const_iterator const
get_last_lte(T t, const std::vector<T>& vec) {
    auto first_gt = std::upper_bound(vec.begin(), vec.end(), t);
    if (first_gt == vec.begin()) {
        return vec.end();
    } else {
        return std::prev(first_gt);
    }
}

#endif
