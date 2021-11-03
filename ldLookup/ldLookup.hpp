#ifndef LDLOOKUP_LDLOOKUP_HPP
#define LDLOOKUP_LDLOOKUP_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"

#define CHECK_FAIL(file, msg) if ((file).fail()) { throw std::runtime_error((msg)); }

/** Standard function to split a string at a delimiter. */
std::vector<std::string> split_str(const std::string& s, const char delimiter);

/** Provides fast line-oriented parsing of genetic records.*/
class RecordParser {
    public:
        std::string snp_a;
        std::string snp_b;
        float r2;

        /** Constructs a RecordParser.
         * 
         * snp_a_index - Index to column of data containing SNP A values.
         * snp_b_index - Index to column of data containing SNP B values.
         * r2_index - Index to column of data containing r-squared values
         * delimiter - Character used to separate columns of data
         * min_r2 - Minimum r-squared value to include a key-value pair in the table
         */
        RecordParser(const size_t snp_a_index = 2,
                     const size_t snp_b_index = 6,
                     const size_t r2_index = 8,
                     const char delimiter = ' ',
                     const float min_r2 = 0) :
                     snp_a_index(snp_a_index),
                     snp_b_index(snp_b_index),
                     r2_index(r2_index),
                     delimiter(delimiter),
                     min_r2(min_r2) {
            last_col_index = std::max({ snp_a_index, snp_b_index, r2_index });
        }

        /** Parses one genetic record.
         * 
         * line - One genetic record in space-delimited format.
         * 
         * When a correctly formed line is passed, parse_line returns true. The
         * class attributes snp_a, snp_b, and r2 are updated to the values
         * parsed from `line`. When a malformed line is passed, parse_line
         * returns true and no guarantee about snp_a, snp_b, or r2 is made.
         */
        bool parse_line(const std::string& line);

    private:
        int snp_a_index;
        int snp_b_index;
        int r2_index;
        char delimiter;
        float min_r2;
        int last_col_index;
};

/** Persistently maps genetic markers in linkage disequilibrium. */
class LDTable {
    public:
        LDTable(const std::string& name) {
            open(name);
        }

        LDTable(const std::string& name,
                const std::string& source_path,
                RecordParser parser,
                const size_t max_key_length) {
            create(name, source_path, parser, max_key_length);
        }

        /** Opens an existing LDTable by name. */
        void open(const std::string& name);

        /** Creates a new LDTable.
         * 
         * name - The name of the new LDTable.
         * source_path - Path to data used to populate the table.
         * parser - RecordParser to transform text data into C++ types.
         * max_key_length - Maximum table key length in bytes.
         */
        void create(const std::string& name,
                    const std::string& source_path,
                    RecordParser parser,
                    const size_t max_key_length);

        /** Fetches values associated with a key */
        std::vector<std::string> get(const std::string& key);

    private:
        static const char KEY_DELIMITER = '\n';
        static const char VALUE_DELIMITER = '\t';

        std::shared_ptr<dht::DiskHash<std::streampos>> hashtable;
        std::fstream file;

        void write_key(const std::string& key);
        void write_value(const std::string& value);
        void write_max_key_length(const size_t max_key_length);
        const size_t read_max_key_length();
};

#endif
