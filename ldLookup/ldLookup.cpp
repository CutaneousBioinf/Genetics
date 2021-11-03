#include <map>

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

bool RecordParser::parse_line(const std::string& line) {
    size_t end = 0;
    size_t start = line.find_first_not_of(delimiter, end);
    int col = 0;
    while (start != std::string::npos && col <= last_col_index) {
        end = line.find(delimiter, start);

        if (col == snp_a_index) {
            snp_a = line.substr(start, end - start);
        } else if (col == snp_b_index) {
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

void LDTable::open(const std::string& name) {
    auto path = name + ".dat";
    file.open(path, std::ios_base::in);
    CHECK_FAIL(file, "Error opening file '" + path + "'");

    // Max key length is the first line of the data file.
    path = name + ".dht";
    hashtable.reset(new dht::DiskHash<std::streampos>(path.c_str(),
                                                      read_max_key_length(),
                                                      dht::DHOpenRO));
}

void LDTable::create(const std::string& name, const std::string& source_path,
                     RecordParser parser, const size_t max_key_length) {
    std::fstream data(source_path, std::ios_base::in);
    CHECK_FAIL(data, "Error opening file '" + source_path + "'");

    // Create the LDTable's files on disk.
    auto path = name + ".dat";
    file.open(path, std::ios_base::out | std::ios_base::app);
    CHECK_FAIL(file, "Error opening file '" + path + "'");

    path = name + ".dht";
    hashtable.reset(new dht::DiskHash<std::streampos>(path.c_str(),
                                                      max_key_length,
                                                      dht::DHOpenRW));

    // Max key length must be persisted so the hashtable can be reopened.
    write_max_key_length(max_key_length);

    // Track the number of keys associated with some number of values.
    // (e.g. # of keys with one value, # of keys with two values...)
    std::map<int, int> num_values_map;
    int num_values;

    // Track the number of keys associated with some range of MAFs.
    // (e.g. # of keys with MF 0-0.001, # of keys with MF 0.001-0.002...).
    // Due to the hazards of floating point keys, keys are multiplied by 10K
    // and rounded down.
    std::map<int, int> maf_values_map;

    // Populate the LDTable.
    std::string line, last_snp_a;
    while (data) {
        getline(data, line);
        if (!parser.parse_line(line)) {
            continue;
        }

        if (last_snp_a.compare(parser.snp_a)) {
            write_key(parser.snp_a);
            last_snp_a = parser.snp_a;

            auto emplace_pair(num_values_map.emplace(num_values, 0));
            emplace_pair.first->second++;
            num_values = 0;
        }

        // MAF Parsing to be Implemented
        // map_safe_maf = floor(parser.maf*10_000)
        // auto emplace_pair(maf_values_map.emplace(parser.maf, 0))
        // emplace_pair.first->second++;

        write_value(parser.snp_b);
        num_values++;
        CHECK_FAIL(file, "Possible corruption: table not writable");
    }
}

std::vector<std::string> LDTable::get(const std::string& key) {
    std::streampos* loc = hashtable->lookup(key.c_str());
    if (loc == NULL) {
        throw std::runtime_error("Nonexistent key '"+ key + "'");
    }

    file.seekg(*loc);
    std::string line;
    getline(file, line, KEY_DELIMITER);
    CHECK_FAIL(file, "Failed to get key '"+ key + "'");
    return split_str(line, VALUE_DELIMITER);
}

void LDTable::write_key(const std::string& key) {
    file.put(KEY_DELIMITER);
    std::streampos loc = file.tellp();
    if (!hashtable->insert(key.c_str(), loc)) {
        std::string msg = "Duplicate key '" + key;
        throw std::runtime_error("Duplicate key '" + key);
    }
}

void LDTable::write_value(const std::string& value) {
    file.write(value.c_str(), value.size());
    file.put(VALUE_DELIMITER);
}

void LDTable::write_max_key_length(const size_t max_key_length) {
    auto s = std::to_string(max_key_length);
    file.write(s.c_str(), s.size());
    file.put(KEY_DELIMITER);
}

const size_t LDTable::read_max_key_length() {
    std::string line;
    getline(file, line, KEY_DELIMITER);
    try {
        return std::stoi(line);
    } catch (std::invalid_argument& e) {
        std::string msg = "Possible corruption: incorrect table format";
        throw std::runtime_error(msg);
    }
}
