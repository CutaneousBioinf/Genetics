#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"

#define CHECK_FAIL(file, msg) if ((file).fail()) { throw std::runtime_error((msg)); }

/** Standard function to split a string at a delimiter. */
std::vector<std::string> split_str(const std::string& s, char delimiter) {
    std::vector<std::string> ret;
    size_t prev = 0;
    size_t curr = s.find(delimiter);
    while (curr != std::string::npos) {
        ret.push_back(s.substr(prev, curr - prev));
        prev = curr + 1;
        curr = s.find(delimiter, prev);
    }

    return ret;
}

/** Provides fast line-oriented parsing of genetic records.*/
struct SNPRecordParser {

    static const int SNP_A_INDEX = 2;
    static const int SNP_B_INDEX = 6;
    static const int R2_INDEX = 8;
    static const char DELIMITER = ' ';

    std::string snp_a;
    std::string snp_b;
    float r2;

    /** Parses one genetic record and updates the state of the SNPRecordParser.
     * 
     * line - One genetic record in space-delimited format.
     * 
     * When a correctly formed line is passed, parse_line returns true. The
     * class attributes snp_a, snp_b, and r2 are updated to the values
     * parsed from `line`. When a malformed line is passed, parse_line
     * returns true and no guarantee about snp_a, snp_b, or r2 is made.
     */
    bool parse_line(const std::string& line) {
        size_t prev = 0;
        size_t curr = line.find(DELIMITER);
        int col = 0;
        while (curr != std::string::npos) {
            if (col == SNP_A_INDEX) {
                snp_a = line.substr(prev, curr - prev);
            } else if (col == SNP_B_INDEX) {
                snp_b = line.substr(prev, curr - prev);
            } else if (col == R2_INDEX) {
                try {
                    r2 = std::stof(line.substr(prev, curr - prev));
                } catch (std::invalid_argument& e) {
                    return false;
                }
            }

            prev = curr + 1;
            curr = line.find(DELIMITER, prev);
            col++;
        }
        // If there weren't values for all of snp_a, snp_b, and r2,
        // return false.
        return col > R2_INDEX;
    }
};

/** Persistently maps genetic markers in linkage disequilibrium. */
class LDTable {
    
    private:
        inline static const char KEY_DELIMITER = '\n';
        inline static const char VALUE_DELIMITER = '\t';
        inline static const size_t MAX_KEY_LENGTH = 200;
        // Minimum R2 value for a record to be included in the table.
        // This will be parameterized/removed later in development.
        inline static constexpr float THRESHOLD = 0.8;

        dht::DiskHash<std::streampos> hashtable;
        std::fstream file;

    public:
        /** Opens an existing LDTable by its name */
        LDTable(const std::string& name)
        : hashtable(dht::DiskHash<std::streampos>((name + ".dht").c_str(),
                                                  MAX_KEY_LENGTH,
                                                  dht::DHOpenRO)) {
            auto path = name + ".dat";
            file.open(path, std::ios_base::in);
            CHECK_FAIL(file, "Error opening file '" + path + "'");
        }

        /** Fetches values associated with a key */
        std::vector<std::string> get(const std::string& key) {
            std::streampos* loc = hashtable.lookup(key.c_str());
            if (loc == NULL) {
                throw std::runtime_error("Nonexistent key '"+ key + "'");
            }

            file.seekg(*loc);
            std::string line;
            getline(file, line);
            CHECK_FAIL(file, "Failed to get key '"+ key + "'");
            return split_str(line, VALUE_DELIMITER);
        }

        /** Creates a new LDTable.
         * 
         * name - The name of the new LDTable.
         * source_path - Path to data used to populate the table.
         * 
         * Only correctly formed genetic records with R2 values >= 0.8
         * will be included in the table.
         */
        static void create_table(const std::string& name,
                                 const std::string& source_path) {
            std::fstream data(source_path, std::ios_base::in);
            CHECK_FAIL(data, "Error opening file '" + source_path + "'");

            // Create the LDTable's files on disk.
            auto table_path = name + ".dat";
            auto mode = std::ios_base::out | std::ios_base::app;
            std::fstream table(table_path, mode);
            CHECK_FAIL(table, "Error opening file '" + table_path + "'");

            auto dht_path = (name + ".dht").c_str();
            auto hashtable(dht::DiskHash<std::streampos>(dht_path,
                                                         MAX_KEY_LENGTH,
                                                         dht::DHOpenRW));
            // Populate the LDTable.
            std::string line, last_snp_a;
            SNPRecordParser parser;
            while (data) {
                getline(data, line);
                if (!parser.parse_line(line) || parser.r2 < THRESHOLD) {
                    continue;
                } else if (parser.snp_a.size() > MAX_KEY_LENGTH) {
                    std::string msg = "Key '" + parser.snp_a + "' too long";
                    throw std::runtime_error(msg);
                }
                
                // Consider adding error-checking for cases where key contains
                // whitespace (\n, \t, \r\n, etc).

                if (last_snp_a.compare(parser.snp_a)) {
                    table.write(&KEY_DELIMITER, 1);
                    std::streampos loc = table.tellp();
                    if (!hashtable.insert(parser.snp_a.c_str(), loc)) {
                        std::string msg = "Duplicate key '" + parser.snp_a;
                        throw std::runtime_error(msg);
                    }
                }

                last_snp_a = parser.snp_a;

                table.write(parser.snp_b.c_str(), parser.snp_b.size());
                table.write(&VALUE_DELIMITER, 1);
                CHECK_FAIL(table, "Possible corruption: table not writable");
            }
        }
};

int main (int argc, char** argv) {
    const std::string usage = "USAGE: ldlookupx create [table_name] [source_file]\n       ldlookupx get [table_name] [key]\n";
    if (argc != 4) {
        std::cout << usage;
        return 1;
    }

    try {
        const std::string command = std::string(argv[1]);
        const std::string table = std::string(argv[2]);
        const std::string arg = std::string(argv[3]);
        if (!command.compare("create")) {
            LDTable::create_table(table, arg);
        } else if (!command.compare("get")) {
            LDTable ldt(table);
            for (std::string s : ldt.get(arg)) {
                std::cout << "Value: " << s << "\n";
            }
        } else {
            std::cout << usage;
            return 1;
        }
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }
}