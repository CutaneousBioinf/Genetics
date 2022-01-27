#include <algorithm>  // std::random_shuffle
#include <random>
#include <stdexcept>  // std::invalid_argument,  std::out_of_range, std::runtime_error

#include "tables.hpp"

std::string serialize_surrogates(size_t surrogate_count) {
    return std::to_string(surrogate_count);
}

std::string serialize_maf(double maf) {
    return std::to_string(maf);
}

size_t deserialize_surrogates(const std::string& surrogate_count) {
    return std::stoul(surrogate_count);
}

double deserialize_maf(const std::string& maf) {
    return std::stod(maf);
}

BinsTable::BinsTable(const std::string& name) {
    // Open the table.
    table.reset(new VectorDiskHash(name));

    try {
        // Load LD pairs bins from the table.
        for (auto surrogates : table->get(LD_BINS_KEY)) {
            if (surrogates.find("LOW") == std::string::npos) {
                surrogate_quantiles.emplace_back(
                    deserialize_surrogates(surrogates)
                );
            }
        }

        // Load MAF bins from the table.
        for (auto maf : table->get(MAF_BINS_KEY)) {
            if (maf.find("LOW") == std::string::npos) {
                maf_quantiles.emplace_back(deserialize_maf(maf));
            }
        }
    } catch (std::invalid_argument& e) {
        throw std::runtime_error("BinsTable Corrupted: Invalid bins");
    } catch (std::out_of_range& e) {
        throw std::runtime_error("BinsTable Corrupted: Invalid bins");
    } catch (std::runtime_error& e) {
        throw std::runtime_error("BinsTable Corrupted: Invalid bins");
    }
}

BinsTable::BinsTable(
    const std::string& name,
    const std::vector<size_t>& surrogate_quantiles,
    const std::vector<double>& maf_quantiles
) {
    // Prepare the table.
    table.reset(new VectorDiskHash(name, MAX_KEY_SIZE));

    this->surrogate_quantiles = surrogate_quantiles;
    this->maf_quantiles = maf_quantiles;

    // Save LD pairs bins in the table.
    for (auto surrogates : surrogate_quantiles) {
        table->insert(LD_BINS_KEY, serialize_surrogates(surrogates));
    }

    // Save MAF bins in the table.
    for (auto maf : maf_quantiles) {
        table->insert(MAF_BINS_KEY, serialize_maf(maf));
    }
}

void BinsTable::reserve(
    const std::string& bin,
    size_t space
) {
    table->reserve(bin, space);
}

void BinsTable::insert(
    const std::string& bin,
    const std::string& snp
) {
    table->insert(bin, snp);
}

std::vector<std::string> BinsTable::get(const std::string& bin) {
    if (bin.find("LOW") != std::string::npos) {
        return std::vector<std::string>();
    }
    return table->get(bin);
}

std::vector<std::string> BinsTable::get_random(
    const std::string& bin,
    size_t n_random
) {
    auto snp_vec(get(bin));
    if (snp_vec.size() == 0) {
        return std::vector<std::string>();
    }

    std::vector<std::string> random_snps;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, snp_vec.size() - 1);
    for (size_t i = 0; i < n_random; i++) {
        random_snps.push_back(snp_vec.at(dist(gen)));
    }

    return random_snps;
}

std::string BinsTable::bin(
    size_t surrogate_count,
    double maf
) {
    auto surrogate_bin(get_last_lte(surrogate_count, surrogate_quantiles));
    if (surrogate_bin == surrogate_quantiles.end()) {
        return "LOW";
    }

    auto maf_bin(get_last_lte(maf, maf_quantiles));
    if (maf_bin == maf_quantiles.end()) {
        return "LOW";
    }

    auto key(
        serialize_surrogates(*surrogate_bin)
        + FIELD_SEPARATOR
        + serialize_maf(*maf_bin)
    );
    return key;
}

SNPTable::SNPTable(const std::string& name) {
    table.reset(new VectorDiskHash(name));
}

SNPTable::SNPTable(const std::string& name, const size_t max_key_size) {
    table.reset(new VectorDiskHash(name, max_key_size));
}

void SNPTable::insert(
    const std::string& snp,
    const size_t surrogate_count,
    const double maf
) {
    auto val(
        serialize_surrogates(surrogate_count)
        + FIELD_SEPARATOR
        + serialize_maf(maf)
    );
    table->insert(snp, val);
}

std::pair<size_t, double> SNPTable::get(const std::string& snp) {
    try {
        auto val(split_str(table->get(snp).at(0), FIELD_SEPARATOR));
        auto surrogate_count(deserialize_surrogates(val.at(0)));
        auto maf(deserialize_maf(val.at(1)));
        return std::pair<size_t, double>(surrogate_count, maf);
    } catch (std::invalid_argument& e) {
        auto msg("SNPTable Corrupted: Malformed Key '" + snp + "'");
        throw std::runtime_error(msg);
    } catch (std::out_of_range& e) {
        auto msg("SNPTable Corrupted: Malformed Key '" + snp + "'");
        throw std::runtime_error(msg);
    }
}

LDTable::LDTable(const std::string& name) {
    table.reset(new VectorDiskHash(name));
}

LDTable::LDTable(const std::string& name, const size_t max_key_size) {
    table.reset(new VectorDiskHash(name, max_key_size));
}

void LDTable::insert(const std::string& snp, const std::string& ld_snp) {
    table->insert(snp, ld_snp);
}

std::vector<std::string> LDTable::get(const std::string& snp) {
    return table->get(snp);
}
