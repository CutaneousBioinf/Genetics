#include "CLI11.hpp"
#include "ldLookup.hpp"

int main (int argc, char** argv) {
    std::string description;

    description = "ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants";
    CLI::App app{description};
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
    description = "Minimum r-squared value for a key-value pair to be included in the table";
    create->add_option("-r,--r2-threshold", min_r2, description);

    size_t key_index{ 2 };
    description = "Zero-based index to column of data containing lookup table keys";
    create->add_option("-K,--key-index", key_index, description);

    size_t value_index{ 6 };
    description = "Zero-based index to column of data containing lookup table values";
    create->add_option("-V,--value-index", value_index, description);

    size_t r2_index{ 8 };
    description = "Zero-based index to column of data containing r-squared values";
    create->add_option("-R,--r2-index", r2_index, description);

    char delimiter{ ' ' };
    description = "Character used to separate columns of data";
    create->add_option("-d,--delimiter", delimiter, description);


    // 'Retrieve' Subcommand
    description = "Retrieve values existing from lookup table";
    auto retrieve = app.add_subcommand("retrieve", description);

    std::string key;
    description = "Allele to look up";
    retrieve->add_option("key,-k,--key", key, description)->required();

    CLI11_PARSE(app, argc, argv);


    // Application Logic
    try {
        if (*create) {
            RecordParser rp{key_index, value_index, r2_index, delimiter, min_r2};
            LDTable::create_table(table, source_path, rp);
        } else if (*retrieve) {
            auto opened_table = LDTable(table);
            for (auto v : opened_table.get(key)) {
                std::cout << "Value: " << v << "\n";
            }
        } else {
            throw std::runtime_error("Unknown subcommand passed to the CLI");
        }
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}