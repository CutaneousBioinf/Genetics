#include <cmath>
#include <map>
#include <iterator>

#include "ldLookup.hpp"

LDLookup::LDLookup(const std::string& name) {
    ld_pairs.reset(new VectorDiskHash(name + LD_PAIR_EXTENSION));
    bins.reset(new VectorDiskHash(name + BINS_EXTENSION));
}

LDLookup::LDLookup(const std::string& name,
                   const std::string& source_path,
                   GeneticDataValidator validator) {
    auto ld_pairs_path = name + LD_PAIR_EXTENSION;
    ld_pairs.reset(new VectorDiskHash(ld_pairs_path,
                                      validator.max_key_size));

    auto bins_path = name + BINS_EXTENSION;
    bins.reset(new VectorDiskHash(bins_path, MAX_BIN_KEY_SIZE));

    std::fstream data(source_path, std::ios_base::in);
    CHECK_FAIL(data, "Error opening file '" + source_path + "'");

    // On our first pass through the data, populate the ld_pair table.
    // Track the number of genetic markers associated with particular
    // bins. A bin is identified by a number of LD surrogates and a
    // range of MAF values.
    std::map<std::string, size_t> bin_sizes;
    GeneticData last_key;
    std::string line;
    size_t surrogate_count;
    bool key_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        // A new key has been encountered, so store bin information
        // for the last key.
        if (last_key.snp_a.compare(validator.data.snp_a)) {
            if (key_found) {
                auto key = to_bin_key(surrogate_count, last_key.maf);
                auto emplace_pair(bin_sizes.emplace(key, 0));
                emplace_pair.first->second += last_key.snp_a.size() + 1;
            } else {
                key_found = true;
            }

            last_key = validator.data;
            surrogate_count = 0;
        }

        // Populate the ld_pair table.
        ld_pairs->insert(validator.data.snp_a, validator.data.snp_b);
        // Update bin information.
        surrogate_count++;
    }

    // Reserve bin space.
    for (auto it = bin_sizes.begin(); it != bin_sizes.end(); it++) {
        bins->reserve(it->first, it->second);
    }

    // On our second pass through the data, populate the bin table.
    data.clear();
    data.seekg(0, std::ios::beg);
    last_key = {};
    surrogate_count = 0;
    key_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        // A new key has been encountered, so write it to the bin table.
        if (last_key.snp_a.compare(validator.data.snp_a)) {
            if (key_found) {
                auto key = to_bin_key(surrogate_count, last_key.maf);
                bins->insert(key, last_key.snp_a);
            } else {
                key_found = true;
            }

            last_key = validator.data;
            surrogate_count = 0;
        }

        // Update bin information.
        surrogate_count++;
    }
}

std::vector<std::string> LDLookup::find_ld(const std::string& key) {
    return ld_pairs->get(key);
}

std::vector<std::string> LDLookup::find_similar(const int surrogate_count,
                                                const double maf) {
    try {
        return bins->get(to_bin_key(surrogate_count, maf));
    } catch (std::runtime_error& e) {
        return std::vector<std::string>();
    }
}

std::string LDLookup::to_bin_key(const int surrogate_count,
                                 const double maf) {
    const int maf_bin = floor(256 * log2(maf));
    return std::to_string(surrogate_count) + " " + std::to_string(maf_bin);
}
