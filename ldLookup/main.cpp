#include "CLI11.hpp"
#include "ldLookup.hpp"

void pretty_print_vector(const std::vector<std::string>& values) {
    int i = 0;
    for (auto v : values) {
        if (++i%5 == 0) {
            std::cout << "\n";
        }

        std::cout << v << " ";
    }

    std::cout << "\n";
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
            GeneticDataValidator v{
                delimiter, key_index, value_index, r2_index, maf_index,
                min_r2, max_key_size
            };
            LDLookup(table, source_path, v);
        } else if (*retrieve) {
            LDLookup ldl(table);

            if (file.size()) {
                std::fstream f(file, std::ios_base::in);
                CHECK_FAIL(f, "Failed to open file " + file);

                std::string line;
                std::vector<std::string> values;
                while (f) {
                    getline(f, line);
                    values = ldl.find_ld(line);
                    std::cout << "Key: " << line << "\nValues: ";
                    pretty_print_vector(values);
                }
            }

            std::vector<std::string> values;
            for (auto key : keys) {
                values = ldl.find_ld(key);
                std::cout << "Key: " << key << "\nValues: ";
                pretty_print_vector(values);
            }
        } else if (*bin) {
            LDLookup ldl(table);
            auto values = ldl.find_similar(surrogate_count, maf);
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
