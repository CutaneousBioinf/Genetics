#ifndef LDLOOKUP_UTIL_HPP
#define LDLOOKUP_UTIL_HPP

#include <string>
#include <string_view>
#include <vector>
#include <stdexcept>

#define CHECK_FAIL(file, msg) if (!(file)) { throw std::runtime_error((msg)); }

/* Splits `s` at a `delimiter`, ignoring consecutive/leading/trailing delimiters. */
std::vector<std::string> split_str(const std::string& s, const char delimiter);

/* Splits `s` at a `delimiter`, ignoring consecutive/leading/trailing delimiters. */
std::vector<std::string_view> split_str(const std::string_view& s, const char delimiter);

/* Randomly selects k items from vec with replacement. */
template <typename T>
std::vector<T> sample_with_replacement(const std::vector<T>& vec, size_t k);

#endif
