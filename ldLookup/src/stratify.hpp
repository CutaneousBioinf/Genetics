#ifndef _LDLOOKUP_STRATIFY_HPP_
#define _LDLOOKUP_STRATIFY_HPP_

#include <stddef.h>  // size_t

#include <map>        // std::map
#include <stdexcept>  // std::runtime_error
#include <vector>

/* Custom Exception for Histogram */
struct histogram_error : std::runtime_error {
	histogram_error(const std::string &msg="")
	    : std::runtime_error("Histogram Error: " + msg) {}
};

/**
 * TODO: Document!
 */
template <typename K>
class Histogram {
    public:
        /**
         * TODO: Document!
         */
        size_t get_count(K key) const;

        /**
         * TODO: Document!
         */
        K get_stratum(K key) const;

        /**
         * TODO: Document!
         */
        void increase_count(K key, size_t increase_by=1);

        /**
         * TODO: Document!
         */
        Histogram<K> stratify(size_t bins) const;

        /**
         * TODO: Document!
         */
        std::vector<K> strata() const;

        /**
         * TODO: Document!
         */
        size_t total_count() const;
    
    private:
        std::map<K, size_t> histogram;
};

/*************************************************/
/*************************************************/
/****             Implementations             ****/
/*************************************************/
/*************************************************/

template <typename K>
inline size_t Histogram<K>::get_count(K key) const {
    auto it = histogram.find(key);
    if (it == histogram.end()) {
        throw histogram_error("get_count() Called on Nonexistent Key");
    }
    return it->second;
}

template <typename K>
inline K Histogram<K>::get_stratum(K key) const {
    auto it = histogram.lower_bound(key);
    if (it == histogram.begin()) {
        throw histogram_error("get_stratum() Called on Out-of-Range Key");
    }
    return (--it)->first;
}

template <typename K>
inline void Histogram<K>::increase_count(K key, size_t increase_by) {
    auto emplace_pair = histogram.emplace(key, 0);
    emplace_pair.first->second += increase_by;
}

template <typename K>
inline Histogram<K> Histogram<K>::stratify(size_t bins) const {
    if (bins == 0 || total_count() == 0) {
        throw histogram_error("stratify() Called on Empty Histogram");
    }

    Histogram<K> restratified;
    size_t bin_size = total_count() / bins;
    size_t running_count = 0;
    for (auto it = histogram.rbegin(); it != histogram.rend(); it++) {
        running_count += it->second;
        if (running_count >= bin_size) {
            restratified.histogram.emplace(it->first, running_count);
            running_count = 0;
        }
    }
    
    if (running_count != 0) {
        restratified.histogram.emplace(
            histogram.begin()->first,
            running_count);
    }

    return restratified;
}

template <typename K>
inline std::vector<K> Histogram<K>::strata() const {
    std::vector<K> strata;
    for (auto it = histogram.begin(); it != histogram.end(); it++) {
        strata.emplace_back(it->first);
    }
    return strata;
}

template <typename K>
inline size_t Histogram<K>::total_count() const {
    size_t count = 0;
    for (auto it = histogram.begin(); it != histogram.end(); it++) {
        count += it->second;
    }
    return count;
}

#endif