#include <algorithm>  // std::random_shuffle
#include <stdexcept>  // std::invalid_argument, out_of_range, and runtime_error

#include "tables.hpp"

std::string serialize_ld(size_t n_ld_pairs) {
    return std::to_string(n_ld_pairs);
}

std::string serialize_maf(double maf) {
    return std::to_string(maf);
}

size_t deserialize_ld(const std::string& n_ld_pairs) {
    return std::stoul(n_ld_pairs);
}

double deserialize_maf(const std::string& maf) {
    return std::stod(maf);
}

BinsTable::BinsTable(const std::string& name) {
    // Open the table.
    table.reset(new VectorDiskHash(name));

    try {
        // Load LD pairs bins from the table.
        for (auto ld : table->get(LD_BINS_KEY)) {
            if (ld.find("LOW") == std::string::npos) {
                ld_quantiles.emplace_back(deserialize_ld(ld));
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
    const std::vector<size_t>& ld_quantiles,
    const std::vector<double>& maf_quantiles
) {
    // Prepare the table.
    table.reset(new VectorDiskHash(name, MAX_KEY_SIZE));

    this->ld_quantiles = ld_quantiles;
    this->maf_quantiles = maf_quantiles;

    // Save LD pairs bins in the table.
    for (auto ld : ld_quantiles) {
        table->insert(LD_BINS_KEY, serialize_ld(ld));
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
    std::random_shuffle(snp_vec.begin(), snp_vec.end());
    if (n_random >= snp_vec.size()) {
        return snp_vec;
    } else {
        return std::vector<std::string>(
            snp_vec.begin(),
            snp_vec.begin()+n_random
        );
    }
}

std::string BinsTable::bin(
    size_t n_ld_pairs,
    double maf
) {
    auto ld_bin(get_last_lte(n_ld_pairs, ld_quantiles));
    if (ld_bin == ld_quantiles.end()) {
        return "LOW";
    }

    auto maf_bin(get_last_lte(maf, maf_quantiles));
    if (maf_bin == maf_quantiles.end()) {
        return "LOW";
    }

    auto key(
        serialize_ld(*ld_bin) + FIELD_SEPARATOR + serialize_maf(*maf_bin)
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
    const size_t n_ld_pairs,
    const double maf
) {
    auto val(
        serialize_ld(n_ld_pairs) + FIELD_SEPARATOR + serialize_maf(maf)
    );
    table->insert(snp, val);
}

std::pair<size_t, double> SNPTable::get(const std::string& snp) {
    try {
        auto val(split_str(table->get(snp).at(0), FIELD_SEPARATOR));
        auto n_ld_pairs(deserialize_ld(val.at(0)));
        auto maf(deserialize_maf(val.at(1)));
        return std::pair<size_t, double>(n_ld_pairs, maf);
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

void LDTable::insert(const std::string& snp, const std::string& ld_surrogate) {
    table->insert(snp, ld_surrogate);
}

std::vector<std::string> LDTable::get(const std::string& snp) {
    return table->get(snp);
}
