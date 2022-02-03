#ifndef _LDLOOKUP_SPLIT_HPP_
#define _LDLOOKUP_SPLIT_HPP_

#include <random> // std::random_device, std::mt19337, std::uniform_int_dist...
#include <string>
#include <string_view>
#include <vector>

inline std::vector<std::string> split_to_strings(
    const std::string &s,
    const char delimiter
) {
    std::vector<std::string> ret;
    ret.reserve(s.length() / 2);
    size_t start = s.find_first_not_of(delimiter);
    size_t end;
    while (start != std::string::npos) {
	    end = s.find(delimiter, start);
	    ret.emplace_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
    }
    return ret;
}

inline std::vector<std::string_view> split_to_string_views(
    const std::string_view &s,
    const char delimiter
) {
    std::vector<std::string_view> ret;
    ret.reserve(s.length() / 2);
    size_t start = s.find_first_not_of(delimiter);
    size_t end;
    while (start != std::string::npos) {
	    end = s.find(delimiter, start);
	    ret.emplace_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
    }
    return ret;
}

template <typename T> inline
std::vector<T> sample_with_replacement(const std::vector<T>& vec, size_t k) {
    std::vector<T> ret;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, vec.size() - 1);

    for (size_t i = 0; i < k; i++) {
        ret.push_back(vec.at(dist(gen)));
    }
    
    return ret;
}

#endif