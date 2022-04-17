#include "tables.hpp"

using std::string;
using std::vector;

LDTable::LDTable(const Options& opts)
: VectorDiskHash(opts) {}

StrataTable::StrataTable(
    const string& file_path,
    const string& table_path,
    Histogram<size_t> n_surrogates_strata_in,
    Histogram<double> maf_strata_in)
    : n_surrogates_strata(n_surrogates_strata_in),
      maf_strata(maf_strata_in) {
	table.reset(new VectorDiskHash({
        file_path, table_path, MAX_KEY_SIZE, true
    }));

    n_surrogates_strata.increase_count(0, 0);
    for (size_t n_surrogates : n_surrogates_strata.strata()) {
        table->append(N_SURROGATES_KEY, std::to_string(n_surrogates));
    }

    maf_strata.increase_count(0.0, 0);
    for (double maf : maf_strata.strata()) {
        table->append(MAF_KEY, std::to_string(maf));
    }
}

StrataTable::StrataTable(
    const string& file_path,
    const string& table_path) {
	table.reset(new VectorDiskHash({
        file_path, table_path, MAX_KEY_SIZE, false
    }));

    try {
        for (string n_surrogates_str : table->lookup(N_SURROGATES_KEY)) {
            auto n_surrogates = static_cast<size_t>(std::stoull(n_surrogates_str));
            n_surrogates_strata.increase_count(n_surrogates, 0);
        }

        for (string maf_str : table->lookup(MAF_KEY)) {
            auto maf = std::stod(maf_str);
            maf_strata.increase_count(maf, 0);
        }
    } catch (std::invalid_argument& e) {
        throw std::runtime_error("StrataTable Constructor: Unreadable Strata");
    } catch (std::out_of_range& e) {
        throw std::runtime_error("StrataTable Constructor: Unreadable Strata");
    }
}

void StrataTable::append(const IndexVariantSummary& summary) {
    table->append(get_stratum(summary), summary.variant_id);
}

vector<string> StrataTable::lookup(const IndexVariantSummary& summary) {
    return table->lookup(get_stratum(summary));
}

vector<string> StrataTable::lookup_sample(
    const IndexVariantSummary& summary,
    const size_t k) {
    return table->lookup_sample(get_stratum(summary), k);
}

void StrataTable::reserve(
    const Histogram<StrataTable::Stratum>& strata_sizes) {
    for (auto &stratum : strata_sizes.strata()) {
        size_t bytes_to_reserve = strata_sizes.get_count(stratum);
        table->reserve(stratum, bytes_to_reserve);
    }
}

StrataTable::Stratum
StrataTable::get_stratum(const IndexVariantSummary& summary) {
    size_t surr_stratum = n_surrogates_strata.get_stratum(summary.n_surrogates);
    double maf_stratum = maf_strata.get_stratum(summary.maf);
    return std::to_string(surr_stratum) + " " + std::to_string(maf_stratum);
}

SummaryTable::SummaryTable(const Options& opts) {
    table.reset(new VectorDiskHash(opts));
}

void SummaryTable::append(const IndexVariantSummary& summary) {
    vector<string> values = {
        std::to_string(summary.maf),
        std::to_string(summary.n_surrogates)
    };
    table->append(summary.variant_id, values);
}

IndexVariantSummary SummaryTable::lookup(const string &index_variant_id) {
    IndexVariantSummary ret{ index_variant_id, 0.0, 0 };
    vector<string> lookup_values = table->lookup(index_variant_id);
    try {
        ret.maf = std::stod(lookup_values.at(0));
        auto n_surrogates = std::stoull(lookup_values.at(1));
        ret.n_surrogates = static_cast<size_t>(n_surrogates);
        return ret;
    } catch (std::invalid_argument& e) {
        auto msg = "SummaryTable lookup(): Corrupted Value for Key ";
        throw std::runtime_error(msg + index_variant_id);
    } catch (std::out_of_range& e) {
        auto msg = "SummaryTable lookup(): Corrupted Value for Key ";
        throw std::runtime_error(msg + index_variant_id);
    }
}
