#include "ldLookup.hpp"

std::vector<std::string> split_str(const std::string& s, char delimiter) {
    std::vector<std::string> ret;
	size_t end = 0;
    size_t start = s.find_first_not_of(delimiter, end);
	while (start != std::string::npos) {
		end = s.find(delimiter, start);
		ret.push_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
	}

    return ret;
}

bool SNPRecord::parse_line(const std::string& line) {
    size_t end = 0;
    size_t start = line.find_first_not_of(DELIMITER, end);
    int col = 0;
    while (start != std::string::npos) {
        end = line.find(DELIMITER, start);

        if (col == SNP_A_INDEX) {
            snp_a = line.substr(start, end - start);
        } else if (col == SNP_B_INDEX) {
            snp_b = line.substr(start, end - start);
        } else if (col == R2_INDEX) {
            try {
                r2 = std::stof(line.substr(start, end - start));
            } catch (std::invalid_argument& e) {
                return false;
            }
        }

        start = line.find_first_not_of(DELIMITER, end);
        col++;
    }
    // If there weren't values for all of snp_a, snp_b, and r2,
    // return false.
    return col > R2_INDEX;
}

LDTable::LDTable(const std::string& name)
: hashtable(dht::DiskHash<std::streampos>((name + ".dht").c_str(),
                                          MAX_KEY_LENGTH,
                                          dht::DHOpenRO)) {
    auto path = name + ".dat";
    file.open(path, std::ios_base::in);
    CHECK_FAIL(file, "Error opening file '" + path + "'");
}

std::vector<std::string> LDTable::get(const std::string& key) {
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

void LDTable::create_table(const std::string& name,
                           const std::string& source_path,
                           const float min_r2) {
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
    SNPRecord parser;
    while (data) {
        getline(data, line);
        if (!parser.parse_line(line) || parser.r2 < min_r2) {
            continue;
        } else if (parser.snp_a.size() > MAX_KEY_LENGTH) {
            std::string msg = "Key '" + parser.snp_a + "' too long";
            throw std::runtime_error(msg);
        }
        
        // Consider adding error-checking for cases where key contains
        // whitespace (\n, \t, \r\n, etc).

        if (last_snp_a.compare(parser.snp_a)) {
            table.put(KEY_DELIMITER);
            std::streampos loc = table.tellp();
            if (!hashtable.insert(parser.snp_a.c_str(), loc)) {
                std::string msg = "Duplicate key '" + parser.snp_a;
                throw std::runtime_error(msg);
            }
        }

        last_snp_a = parser.snp_a;

        table.write(parser.snp_b.c_str(), parser.snp_b.size());
        table.put(VALUE_DELIMITER);
        CHECK_FAIL(table, "Possible corruption: table not writable");
    }
}