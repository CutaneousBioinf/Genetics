#ifndef _LDLOOKUP_STRING_OPS_HPP_
#define _LDLOOKUP_STRING_OPS_HPP_

#include <stddef.h>  // size_t

#include <string>
#include <string_view>
#include <vector>

/**
 *  TODO: Document!
 */
std::string join(
    const std::vector<std::string> &vec,
    const char delimiter,
    const bool trailing = true);

/**
 *  TODO: Document!
 */
std::vector<std::string> sample(
    const std::string &s,
	const char delimiter,
    size_t k);

/**
 *  TODO: Document!
 */
template <typename StringLike>
std::vector<StringLike> split(
    const StringLike &s,
    const char delimiter);

/*************************************************/
/*************************************************/
/****         Template Implementation         ****/
/*************************************************/
/*************************************************/

template <typename StringLike>
inline std::vector<StringLike> split(
    const StringLike &s,
    const char delimiter) {
	// Container for split strings
	std::vector<StringLike> ret;

	// Reserve space for performance reasons
	ret.reserve(s.size() / 2);

	// Repeatedly find substrings between occurences of delimiter
	size_t start = s.find_first_not_of(delimiter);
	size_t end;
	while (start != std::string::npos) {
		end = s.find(delimiter, start);
		ret.emplace_back(s.substr(start, end - start));
		start = s.find_first_not_of(delimiter, end);
	}

	// Deallocate unused space
	ret.shrink_to_fit();
	return ret;
}

#endif