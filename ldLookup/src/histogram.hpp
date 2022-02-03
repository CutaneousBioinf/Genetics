#ifndef _LDLOOKUP_HISTOGRAM_HPP
#define _LDLOOKUP_HISTOGRAM_HPP

#include <algorithm> // std::upper_bound
#include <map>
#include <vector>

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