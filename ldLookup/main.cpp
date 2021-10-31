#include "CLI11.hpp"
#include "ldLookup.hpp"

void pretty_print_table_entry(const std::string& key,
                              const std::vector<std::string>& values) {
    std::cout << "Key: " << key << "\n" << "Values: ";
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

    size_t value_index{ 6 };
    description = "Zero-based index to column of data containing lookup table values";
    create->add_option("-V,--value-index", value_index, description);

    size_t r2_index{ 8 };
    description = "Zero-based index to column of data containing r-squared values";
    create->add_option("-R,--r2-index", r2_index, description);

    char delimiter{ ' ' };
    description = "Character used to separate columns of data";
    create->add_option("-d,--delimiter", delimiter, description);

    size_t max_key_length{ 200 };
    description = "Maximum key length in bytes";
    create->add_option("-k,--keys", max_key_length, description);

    // 'Retrieve' Subcommand
    description = "Retrieve values existing from lookup table";
    auto retrieve = app.add_subcommand("retrieve", description);

    std::string file;
    description = "File containing alleles to look up (one per line)";
    retrieve->add_option("file,-f,--file", file, description);

    std::vector<std::string> keys;
    description = "Alleles to look up";
    retrieve->add_option("keys,-k,--keys", keys, description);

    CLI11_PARSE(app, argc, argv);


    // Application Logic
    try {
        if (*create) {
            RecordParser rp{key_index, value_index, r2_index, delimiter, min_r2};
            LDTable::create_table(table, source_path, rp, max_key_length);
        } else if (*retrieve) {
            auto opened_table = LDTable(table);

            if (file.size()) {
                std::fstream f(file, std::ios_base::in);
                CHECK_FAIL(f, "Failed to open file " + file);

                std::string line;
                std::vector<std::string> values;
                while (f) {
                    getline(f, line);
                    values = opened_table.get(line);
                    pretty_print_table_entry(line, values);
                }
            }

            std::vector<std::string> values;
            for (auto key : keys) {
                values = opened_table.get(key);
                pretty_print_table_entry(key, values);
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
