#include "CLI11.hpp"
#include "geneticData.hpp"
#include "tables.hpp"

const std::string BINSTABLE_EXT = "_binstable";
const std::string LDTABLE_EXT = "_ldtable";
const std::string SNPTABLE_EXT = "_snptable";

void pretty_print_vector(const std::vector<std::string>& values) {
    std::string sep = "";
    for (auto v : values) {
        std::cout << sep << v;
        sep = " ";
    }

    std::cout << "\n";
}

void pretty_print_key_and_values(
    const std::string key,
    const std::vector<std::string>& values
) {
    std::cout << "Key: " << key << "\nValues:\n";
    pretty_print_vector(values);
    std::cout << "\n";
}

/**
 * Extracts valid lines of linkage disequilibrium data from a file.
 * 
 * Parameters:
 *  source_path - File to read from
 *  validator - Line parser that ignores malformed lines and genetic
 *      markers that are not in LD.
 *  F1 - Function to be called when new LD data is parsed. F1 is
 *      passed a GeneticData object representing the new data.
 *  F2 - Function to be called for each unique index SNP. F2 is passed
 *      the index SNP ID, the index SNP's MAF, and the number of SNPs
 *      in LD with the index SNP.
**/
template <typename F1, typename F2>
void iterate_genetic_data(
    const std::string& source_path,
    GeneticDataValidator validator,
    F1 on_data,
    F2 on_new_snp
) {
    std::fstream data(source_path, std::ios_base::in);
    CHECK_FAIL(data, "Error opening file '" + source_path + "'");

    GeneticData last_record;
    std::string line;
    size_t ld_pairs = 0;
    bool record_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        if (last_record.index_snp_id.compare(validator.data.index_snp_id)
            && record_found) {
            on_new_snp(
                last_record.index_snp_id,
                last_record.maf,
                ld_pairs
            );
            last_record = validator.data;
            ld_pairs = 0;
        }

        on_data(validator.data);
        ld_pairs++;
        record_found = true;
    }

    if (record_found) {
        on_new_snp(
            last_record.index_snp_id,
            last_record.maf,
            ld_pairs
        );
    }
}

/**
 * Creates all of ldLookup's data tables.
 * 
 * Parameters:
 *  source_path - File containing linkage disequilibrium data.
 *  name - A string used to identify and reopen tables created by
 *      this call.
 *  validator - Line parser that ignores malformed lines and genetic
 *      markers that are not in LD.
 *  n_ld_bins - Approximate number of BinsTable divisions by number
 *      of LD surrogates.
 *  n_maf_bins - Approximate number of BinsTable divisions by MAF.
**/  
void create_tables(
    const std::string& source_path,
    const std::string& table,
    GeneticDataValidator validator,
    size_t n_ld_bins,
    size_t n_maf_bins
) {
    // First Pass: Populate LDTable and SNPTable.
    // Also, gather distribution data for MAF/LD.
    LDTable ldt(table + LDTABLE_EXT, validator.index_snp_id_max_length);
    SNPTable snpt(table + SNPTABLE_EXT, validator.index_snp_id_max_length);

    std::map<size_t, size_t> ld_counts_hist;
    std::map<double, size_t> maf_counts_hist;

    auto first_pass_on_data = [&ldt] (const GeneticData& data) {
        ldt.insert(data.index_snp_id, data.ld_snp_id);
    };

    auto first_pass_on_new_snp = [&snpt, &ld_counts_hist, &maf_counts_hist] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        snpt.insert(index_snp_id, n_ld_pairs, maf);
        ld_counts_hist.emplace(n_ld_pairs, 0).first->second++;
        maf_counts_hist.emplace(maf, 0).first->second++;
    };

    iterate_genetic_data(
        source_path,
        validator,
        first_pass_on_data,
        first_pass_on_new_snp
    );

    // Determine bin cutpoints.
    auto ld_quantiles(get_map_keys(
        bin_histogram(ld_counts_hist, n_ld_bins)
    ));

    auto maf_quantiles(get_map_keys(
        bin_histogram(maf_counts_hist, n_maf_bins)
    ));

    // Initialize BinsTable.
    BinsTable bt(table + BINSTABLE_EXT, ld_quantiles, maf_quantiles);

    // Second Pass: Determine exact bin sizes.
    std::map<std::string, size_t> bin_sizes;

    auto second_pass_on_data = [] (const GeneticData& data) {};

    auto second_pass_on_new_snp = [&bt, &bin_sizes] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        auto key(bt.bin(n_ld_pairs, maf));
        auto size_added(index_snp_id.size() + 1);
        bin_sizes.emplace(key, 0).first->second += size_added;
    };

    iterate_genetic_data(
        source_path,
        validator,
        second_pass_on_data,
        second_pass_on_new_snp
    );

    // Reserve space in BinsTable.
    for (const auto& [bin, size] : bin_sizes) {
        bt.reserve(bin, size);
    }

    // Third Pass: Populate BinsTable.
    auto third_pass_on_data = [] (const GeneticData& data) {};

    auto third_pass_on_new_snp = [&bt, &bin_sizes] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        bt.insert(bt.bin(n_ld_pairs, maf), index_snp_id);
    };

    iterate_genetic_data(
        source_path,
        validator,
        third_pass_on_data,
        third_pass_on_new_snp
    );
}


void do_retrieval(
    LDTable& ldt,
    const std::string& keys_path,
    const std::vector<std::string>& keys
) {
    if (keys_path.size()) {
        std::fstream f(keys_path, std::ios_base::in);
        CHECK_FAIL(f, "Failed to open file " + keys_path);

        std::string line;
        while (f) {
            getline(f, line);
            std::vector<std::string> values = ldt.get(line);
            pretty_print_key_and_values(line, values);
        }
    }

    for (auto key : keys) {
        std::vector<std::string> values = ldt.get(key);
        pretty_print_key_and_values(key, values);
    }
}


void do_bin(
    BinsTable& bt,
    const double maf,
    const size_t surrogate_count
) {
    auto values = bt.get(bt.bin(surrogate_count, maf));
    std::cout << "Key: " << maf << " " << surrogate_count<< "\nValues:\n";
    pretty_print_vector(values);
}

void do_bin_snp(
    BinsTable& bt,
    SNPTable& snpt,
    const std::string& keys_path,
    const std::vector<std::string>& keys
) {
    if (keys_path.size()) {
        std::fstream f(keys_path, std::ios_base::in);
        CHECK_FAIL(f, "Failed to open file " + keys_path);

        std::string line;
        while (f) {
            getline(f, line);
            auto maf_and_surrogates = snpt.get(line);
            auto values = bt.get(bt.bin(
                maf_and_surrogates.first,
                maf_and_surrogates.second
            ));
            pretty_print_key_and_values(line, values);
        }
    }

    for (auto key : keys) {
        auto maf_and_surrogates = snpt.get(key);
        auto values = bt.get(bt.bin(
            maf_and_surrogates.first,
            maf_and_surrogates.second
        ));
        pretty_print_key_and_values(key, values);
    }
}


int main(int argc, char** argv) {
    std::string description;

    description = "ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants";
    CLI::App app{description};
    app.option_defaults()->always_capture_default();
    app.require_subcommand(1);

    std::string name;
    description = "Identifier for tables operated on by this command";
    app.add_option("name,-n,--name", name, description)->required();

    // 'create' Subcommand
    auto create = app.add_subcommand("create", "Create ldLookup tables from LD data");

    std::string create_path;
    description = "File containing source LD data";
    create->add_option("path,-p,--path", create_path, description)->required();

    float threshold_ld_score{ 0 };
    description = "Minimum LD score, as determined by r-squared value, for genetic markers to be processed";
    create->add_option("-t,--threshold", threshold_ld_score, description);

    size_t index_snp_id_col{ 2 };
    description = "Zero-based index to column of data containing index SNP IDs";
    create->add_option("-I,--index-snp-id-column", index_snp_id_col, description);

    size_t maf_col{ 3 };
    description = "Zero-based index to column of data containing MAFs for index SNPs";
    create->add_option("-M,--maf-col", maf_col, description);

    size_t ld_snp_id_col{ 6 };
    description = "Zero-based index to column of data containing IDs of SNPs in LD with index SNP";
    create->add_option("-L,--ld-snp-id-column", ld_snp_id_col, description);

    size_t ld_score_col{ 8 };
    description = "Zero-based index to column of data containing r-squared values";
    create->add_option("-R,--r2-index", ld_score_col, description);

    char delimiter{ ' ' };
    description = "Character used to separate columns of data";
    create->add_option("-d,--delimiter", delimiter, description);

    size_t max_key_size{ 200 };
    description = "Maximum index SNP ID length in bytes";
    create->add_option("-k,--keys", max_key_size, description);

    size_t n_ld_bins{ 15 };
    description = "Approximate number of groupings for index SNPs by number of LD surrogates";
    create->add_option("-l,--ld-bin-count", n_ld_bins, description);

    size_t n_maf_bins{ 15 };
    description = "Approximate number of groupings for index SNPs by MAF ";
    create->add_option("-m,--maf-bin-count", n_maf_bins, description);

    // 'Retrieve' Subcommand
    description = "Retrieve SNPs in LD with a particular index SNP";
    auto retrieve = app.add_subcommand("retrieve", description);

    std::string retrieve_path;
    description = "File containing newline-separated index SNP IDs to retrieve";
    retrieve->add_option("path,-p,--path", retrieve_path, description);

    std::vector<std::string> retrieve_keys;
    description = "Space-separated SNP IDs to retrieve";
    retrieve->add_option("keys,-k,--keys", retrieve_keys, description);

    // 'Bin' Subcommand
    description = "Retrieve markers with similar MAF and number of LD surrogates";
    auto bin = app.add_subcommand("bin", description);
    bin->option_defaults()->always_capture_default(false);

    double maf;
    bin->add_option("maf,-m,--maf", maf, description)->required();

    size_t surrogate_count;
    bin->add_option("surrogates,-s,--surrogates", surrogate_count, description)->required();

    // 'Bin_SNP' Subcommand
    description = "Retrieve markers with similar MAF and number of LD surrogates to an index marker";
    auto bin_snp = app.add_subcommand("bin_snp", description);
    bin_snp->option_defaults()->always_capture_default(false);

    std::string bin_snp_path;
    description = "File containing newline-separated index SNP IDs to bin";
    bin_snp->add_option("path,-p,--path", bin_snp_path, description);

    std::vector<std::string> bin_snp_keys;
    description = "Space-separated SNP IDs to bin";
    bin_snp->add_option("keys,-k,--keys", bin_snp_keys, description);

    CLI11_PARSE(app, argc, argv);


    // Application Logic
    try {
        if (*create) {
            GeneticDataValidator validator(
                delimiter,
                index_snp_id_col,
                ld_snp_id_col,
                ld_score_col,
                maf_col,
                threshold_ld_score,
                max_key_size
            );
            create_tables(create_path, name, validator, n_ld_bins, n_maf_bins);
        } else {
            LDTable ldt(name + LDTABLE_EXT);
            BinsTable bt(name + BINSTABLE_EXT);
            SNPTable snpt(name + SNPTABLE_EXT);

            if (*retrieve) {
                do_retrieval(ldt, retrieve_path, retrieve_keys);
            } else if (*bin) {
                do_bin(bt, maf, surrogate_count);
            } else if (*bin_snp) {
                do_bin_snp(bt, snpt, bin_snp_path, bin_snp_keys);
            } else {
                throw std::runtime_error("Unknown subcommand passed to the CLI");
            }
        }
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
