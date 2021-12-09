#ifndef LDLOOKUP_GENETICDATA_HPP
#define LDLOOKUP_GENETICDATA_HPP

#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "global.hpp"

struct GeneticData {
    std::string snp_a;
    std::string snp_b;
    double r2;
    double maf;
};

/** Provides fast line-oriented parsing of genetic records.*/
class GeneticDataValidator {
	public:
		char delimiter;
		size_t snp_a_index;
		size_t snp_b_index;
		size_t r2_index;
		size_t maf_index;
		double min_r2;
        size_t max_key_size;
		GeneticData data;

		GeneticDataValidator(const char delimiter,
                             const size_t snp_a_index,
                             const size_t snp_b_index,
                             const size_t r2_index,
                             const size_t maf_index,
                             const double min_r2,
                             const size_t max_key_size) :
                        delimiter(delimiter),
                        snp_a_index(snp_a_index),
                        snp_b_index(snp_b_index),
                        r2_index(r2_index),
                        maf_index(maf_index),
                        min_r2(min_r2),
                        max_key_size(max_key_size),
                        last_important_column(std::max({ snp_a_index,
                                                         snp_b_index,
                                                         maf_index,
                                                         r2_index }))
                        {}

        /** Parses one line of genetic data into C++ types.
         * 
         * If `line` is valid, the `data` member will contain
         * the parsed data from `line` after this call. If
         * `line` is invalid, the `data` member may contain
         * junk.
         */
		bool validate(const std::string& line);

	private:
		size_t last_important_column;
};

template <class Key>
std::map<Key, size_t> bin_histogram(
    const std::map<Key, size_t>& hist,
    const size_t n_bins
) {
    if (!n_bins || hist.empty()) {
        return std::map<Key, size_t>();
    }

    size_t total = 0;
    for (const auto& [k, v] : hist) {
        total += v;
    }

    std::map<Key, size_t> bins;
    size_t bin_spacing = (total + (n_bins / 2)) / n_bins;
    size_t bin_size = hist.begin()->second;
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

template <class Key, class Value>
std::vector<Key> get_map_keys(const std::map<Key, Value>& map) {
    std::vector<Key> ret;
    for (const auto& [k, v] : map) {
        ret.push_back(k);
    }
    return ret;
}

template <class T>
typename std::vector<T>::const_iterator const
get_first_lte(T t, const std::vector<T> vec) {
    auto first_gt = std::lower_bound(vec.begin(), vec.end(), t);
    if (first_gt == vec.begin()) {
        return vec.end();
    } else {
        return first_gt--;
    }
}

#endif
