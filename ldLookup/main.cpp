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

template <typename F1, typename F2>
void iterate_genetic_data(
    const std::string& source_path,
    GeneticDataValidator validator,
    F1 on_data,
    F2 on_new_snp
) {
    std::fstream data(source_path, std::ios_base::in);
    CHECK_FAIL(data, "Error opening file '" + source_path + "'");

    GeneticData last_key;
    std::string line;
    size_t ld_pairs = 0;
    bool key_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        if (last_key.snp_a.compare(validator.data.snp_a) && key_found) {
            on_new_snp(last_key, ld_pairs);
            last_key = validator.data;
            ld_pairs = 0;
        }

        on_data(validator.data);
        ld_pairs++;
        key_found = true;
    }

    if (key_found) {
        on_new_snp(last_key, ld_pairs);
    }
}

void create_tables(
    const std::string& source_path,
    const std::string& table,
    GeneticDataValidator validator,
    size_t n_ld_bins,
    size_t n_maf_bins
) {
    // First Pass: Populate LDTable and SNPTable.
    // Also, gather distribution data for MAF/LD.
    LDTable ldt(table + LDTABLE_EXT, validator.max_key_size);
    SNPTable snpt(table + SNPTABLE_EXT, validator.max_key_size);

    std::map<size_t, size_t> ld_counts_hist;
    std::map<double, size_t> maf_counts_hist;

    auto first_pass_on_data = [&ldt] (const GeneticData& data) {
        ldt.insert(data.snp_a, data.snp_b);
    };

    auto first_pass_on_new_snp = [&snpt, &ld_counts_hist, &maf_counts_hist] (
        const GeneticData& data, size_t n_ld_pairs
    ) {
        snpt.insert(data.snp_a, n_ld_pairs, data.maf);
        ld_counts_hist.emplace(n_ld_pairs, 0).first->second++;
        maf_counts_hist.emplace(data.maf, 0).first->second++;
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
        const GeneticData& data, size_t n_ld_pairs
    ) {
        auto key(bt.bin(n_ld_pairs, data.maf));
        auto size_added(data.snp_a.size() + 1);
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
        const GeneticData& data, size_t n_ld_pairs
    ) {
        bt.insert(bt.bin(n_ld_pairs, data.maf), data.snp_a);
    };

    iterate_genetic_data(
        source_path,
        validator,
        third_pass_on_data,
        third_pass_on_new_snp
    );
}

int main(int argc, char** argv) {
    std::string description;

    description = "ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants";
    CLI::App app{description};
    app.option_defaults()->always_capture_default();
    app.require_subcommand(1);

    // Global Options
    std::string table;
    description = "Name of lookup table to search or create";
    app.add_option("table,-t,--t", table, description)->required();

    // 'Create' Subcommand
    auto create = app.add_subcommand("create", "Create new lookup table");

    std::string source_path;
    description = "File containing source data";
    create->add_option("source_path", source_path, description)->required();

    float min_r2{ 0 };
    description = "Minimum r-squared value to include a key-value pair in the table";
    create->add_option("-r,--r2-threshold", min_r2, description);

    size_t key_index{ 2 };
    description = "Zero-based index to column of data containing lookup table keys";
    create->add_option("-K,--key-index", key_index, description);

    size_t maf_index{ 3 };
    description = "Zero-based index to column of data containing MAF values";
    create->add_option("-M,--maf-index", maf_index, description);

    size_t value_index{ 6 };
    description = "Zero-based index to column of data containing lookup table values";
    create->add_option("-V,--value-index", value_index, description);

    size_t r2_index{ 8 };
    description = "Zero-based index to column of data containing r-squared values";
    create->add_option("-R,--r2-index", r2_index, description);

    char delimiter{ ' ' };
    description = "Character used to separate columns of data";
    create->add_option("-d,--delimiter", delimiter, description);

    size_t max_key_size{ 200 };
    description = "Maximum key size in bytes";
    create->add_option("-k,--keys", max_key_size, description);

    // 'Retrieve' Subcommand
    description = "Retrieve values from existing lookup table";
    auto retrieve = app.add_subcommand("retrieve", description);

    std::string file;
    description = "File containing alleles to look up (one per line)";
    retrieve->add_option("file,-f,--file", file, description);

    std::vector<std::string> keys;
    description = "Alleles to look up";
    retrieve->add_option("keys,-k,--keys", keys, description);

    // 'Bin' Subcommand
    description = "Retrieve markers with similar MAF and number of LD surrogates";
    auto bin = app.add_subcommand("bin", description);
    bin->option_defaults()->always_capture_default(false);

    double maf;
    description = "MAF value to bin";
    bin->add_option("maf,-m,--maf", maf, description)->required();

    size_t surrogate_count;
    description = "Number of LD surrogates";
    bin->add_option("surrogates,-s,--surrogates", surrogate_count, description)->required();

    CLI11_PARSE(app, argc, argv);


    // Application Logic
    try {
        if (*create) {
            GeneticDataValidator validator {
                delimiter, key_index, value_index, r2_index, 
                maf_index, min_r2, max_key_size
            };
            create_tables(source_path, table, validator, 15, 15);
        } else if (*retrieve) {
            LDTable ldt(table + LDTABLE_EXT);

            if (file.size()) {
                std::fstream f(file, std::ios_base::in);
                CHECK_FAIL(f, "Failed to open file " + file);

                std::string line;
                std::vector<std::string> values;
                while (f) {
                    getline(f, line);
                    values = ldt.get(line);
                    std::cout << "Key: " << line << "\nValues: ";
                    pretty_print_vector(values);
                }
            }

            std::vector<std::string> values;
            for (auto key : keys) {
                values = ldt.get(key);
                std::cout << "Key: " << key << "\nValues: ";
                pretty_print_vector(values);
            }
        } else if (*bin) {
            BinsTable bt(table + BINSTABLE_EXT);
            auto values = bt.get(bt.bin(surrogate_count, maf));
            std::cout << "Values: ";
            pretty_print_vector(values);
        } else {
            throw std::runtime_error("Unknown subcommand passed to the CLI");
        }
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
