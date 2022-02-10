#include <cassert>
#include <iostream>
#include <map>
#include <fstream>
#include "histogram.hpp"

#define LOG_TEST(success, name) std::cout << name << ": " << ((success) ? "Passed" : "Failed") << '\n'

template <typename Map>
bool maps_equal(const Map& m1, const Map& m2) {
    return m1.size() == m2.size()
           && std::equal(m1.begin(), m1.end(), m2.begin());
}

template <typename Map>
void pretty_print_map(const Map& m) {
    std::cout << "Map Items:\n\t";
    for (const auto& [k, v] : m) {
        std::cout << k << " : " << v << "\n\t";
    }
    std::cout << "\n";
}

void test_bin_histogram(const std::map<int, size_t>& histogram) {
    auto binned = bin_histogram(histogram, 10);
    std::map<int, size_t> expected {
        std::make_pair(0, 30),
        std::make_pair(11, 31),
        std::make_pair(48, 30),
        std::make_pair(119, 30),
        std::make_pair(200, 30),
        std::make_pair(295, 30),
        std::make_pair(401, 30),
        std::make_pair(533, 30),
        std::make_pair(706, 30),
        std::make_pair(865, 29)
    };
    
    assert(maps_equal(binned, expected));
}

void test_bin_histogram_large_bin(const std::map<int, size_t>& histogram) {
    auto binned = bin_histogram(histogram, 1);
    std::map<int, size_t> expected { std::make_pair(0, 300) };
    assert(maps_equal(binned, expected));
}

void test_bin_histogram_small_bins(const std::map<int, size_t>& histogram) {
    auto binned = bin_histogram(histogram, histogram.size()*3);
    assert(maps_equal(binned, histogram));
}

void test_get_map_keys() {
    std::map<long, int> test_map = { {1, 0}, {2, 1} };
    std::vector<long> keys = { 1, 2 };
    for (int i = 0; i < 50; i++) {
        auto fib = keys.at(i) + keys.at(i+1);
        keys.push_back(fib);
        test_map.emplace(fib, i);
    }

    assert(keys == get_map_keys(test_map));
}

void test_get_last_lte() {
    std::vector<int> lte_vec = { 0, 1, 2, 3 };
    assert(get_last_lte(4, lte_vec) == lte_vec.begin()+3);
    assert(get_last_lte(3, lte_vec) == lte_vec.begin()+3);
    assert(get_last_lte(2, lte_vec) == lte_vec.begin()+2);
    assert(get_last_lte(1, lte_vec) == lte_vec.begin()+1);
    assert(get_last_lte(0, lte_vec) == lte_vec.begin()+0);
    assert(get_last_lte(-1, lte_vec) == lte_vec.end());
}

int main(int argc, char** argv) {
    std::map<int, size_t> histogram;

    std::fstream fin("tst/hist_test_data");
    if (!fin) {
        std::cout << "Failed to open hist_test_data";
        return 1;
    }

    int n;
    while (fin >> n) {
        histogram.emplace(n, 0).first->second++;
    }

    std::cout << "Test bin_histogram (Simulated Data): " << std::flush;
    test_bin_histogram(histogram);
    std::cout << "PASS\nTest bin_histogram (Large Bin): " << std::flush;
    test_bin_histogram_large_bin(histogram);
    std::cout << "PASS\nTest bin_histogram (Small Bins): " << std::flush;
    test_bin_histogram_small_bins(histogram);
    std::cout << "PASS\nTest get_map_keys: " << std::flush;
    test_get_map_keys();
    std::cout << "PASS\nTest get_last_lte: " << std::flush;
    test_get_last_lte();
    std::cout << "PASS\n";
    return 0;
}