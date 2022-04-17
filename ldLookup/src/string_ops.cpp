#include "string_ops.hpp"

#include <random>  // std::mt19937, random_device, uniform_int_distribution

using std::string;
using std::vector;

/**
 *  TODO: Document!
 */
template <typename T>
inline vector<T> sample_with_replacement(const vector<T> &vec, size_t k) {
	// Container for sampled items
	vector<T> sampled;

	// Code to generate random indices into vec
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, vec.size() - 1);

	// Repeatedly generate random indices into vec.
	// Add corresponding items to the sample.
	for (size_t i = 0; i < k; i++) {
		sampled.emplace_back(vec.at(dist(gen)));
	}

	return sampled;
}

string join(
    const vector<string> &vec,
    const char delimiter,
    const bool trailing) {
	// Handle empty vector case.
	if (!vec.size()) {
		return std::string();
	}

	// This will be the joined string.
	std::string ret = vec.at(0);

	// Do the joining.
	for (auto it = ++vec.begin(); it != vec.end(); it++) {
		ret.push_back(delimiter);
		ret += *it;
	}

	// Add a trailing delimiter if needed.
	if (trailing) {
		ret.push_back(delimiter);
	}

	return ret;
}

vector<string> sample(const string &s, const char delimiter, size_t k) {
	// Sample random entries from s without constructing lots of strings.
	std::string_view line_view = s;
	vector<std::string_view> split_lv = split(line_view, delimiter);
	auto samples = sample_with_replacement(split_lv, k);

	vector<string> ret;

	// Reserve space for performance reasons.
	ret.reserve(k);

	// Convert the random entries to strings.
	ret.insert(ret.begin(), samples.begin(), samples.end());
	return ret;
}