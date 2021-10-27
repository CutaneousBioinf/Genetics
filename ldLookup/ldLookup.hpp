#ifndef __ldLookup_hpp__
#define __ldLookup_hpp__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"

#define CHECK_FAIL(file, msg) if ((file).fail()) { throw std::runtime_error((msg)); }

/** Standard function to split a string at a delimiter. */
std::vector<std::string> split_str(const std::string& s, char delimiter);

/** Provides fast line-oriented parsing of genetic records.*/
class SNPRecord {
    public:
        std::string snp_a;
        std::string snp_b;
        float r2;

        /** Parses one genetic record, updating an SNPRecord.
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
        static const int SNP_A_INDEX = 2;
        static const int SNP_B_INDEX = 6;
        static const int R2_INDEX = 8;
        static const char DELIMITER = ' ';
};

/** Persistently maps genetic markers in linkage disequilibrium. */
class LDTable {
    public:
        /** Opens an existing LDTable by its name */
        LDTable(const std::string& name);

        /** Fetches values associated with a key */
        std::vector<std::string> get(const std::string& key);

        /** Creates a new LDTable.
         * 
         * name - The name of the new LDTable.
         * source_path - Path to data used to populate the table.
         * min_r2 - Lowest R-squared value to include in the table.
         */
        static void create_table(const std::string& name,
                                 const std::string& source_path,
                                 const float min_r2);

    private:
        static const char KEY_DELIMITER = '\n';
        static const char VALUE_DELIMITER = '\t';
        static const size_t MAX_KEY_LENGTH = 200;

        dht::DiskHash<std::streampos> hashtable;
        std::fstream file;
};

#endif