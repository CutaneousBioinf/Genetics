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

RecordParser::RecordParser(size_t key_index, size_t value_index,
                           size_t r2_index, char delimiter,
                           float min_r2) {
    this->key_index = key_index;
    this->value_index = value_index;
    this->r2_index = r2_index;
    this->delimiter = delimiter;
    this->min_r2 = min_r2;
    this->last_col_index = std::max({ key_index, value_index, r2_index });
}

bool RecordParser::parse_line(const std::string& line) {
    size_t end = 0;
    size_t start = line.find_first_not_of(delimiter, end);
    int col = 0;
    while (start != std::string::npos && col <= last_col_index) {
        end = line.find(delimiter, start);

        if (col == key_index) {
            snp_a = line.substr(start, end - start);
        } else if (col == value_index) {
            snp_b = line.substr(start, end - start);
        } else if (col == r2_index) {
            try {
                r2 = std::stof(line.substr(start, end - start));
                if (r2 < min_r2) {
                    return false;
                }
            } catch (std::invalid_argument& e) {
                return false;
            }
        }

        start = line.find_first_not_of(delimiter, end);
        col++;
    }
    
    return col > last_col_index;
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
                           RecordParser parser) {
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
    while (data) {
        getline(data, line);
        if (!parser.parse_line(line)) {
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