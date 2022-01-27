#include <random>

#include "utils.hpp"

template std::vector<std::string_view> sample_with_replacement(
    const std::vector<std::string_view>& vec, size_t k
);
template std::vector<std::string> sample_with_replacement(
    const std::vector<std::string>& vec, size_t k
);

std::vector<std::string> split_str(const std::string& s, char delimiter) {
    std::vector<std::string> ret;
    ret.reserve(s.length() / 2); // for performance reasons
    size_t end = 0;
    size_t start = s.find_first_not_of(delimiter, end);
    while (start != std::string::npos) {
	    end = s.find(delimiter, start);
	    ret.emplace_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
    }
    return ret;
}

std::vector<std::string_view> split_str(
    const std::string_view& s, 
    char delimiter
) {
    std::vector<std::string_view> ret;
    ret.reserve(s.length() / 2); // for performance reasons
    size_t end = 0;
    size_t start = s.find_first_not_of(delimiter, end);
    while (start != std::string::npos) {
	    end = s.find(delimiter, start);
	    ret.emplace_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
    }
    return ret;
}

template <typename T>
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