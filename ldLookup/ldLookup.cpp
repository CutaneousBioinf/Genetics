#include "ldLookup.hpp"

#include <cmath>
#include <map>
#include <iterator>

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

VectorDiskHash::VectorDiskHash(const std::string& name) {
    auto path = name + DATA_EXTENSION;
    file.open(path, std::ios_base::binary | std::ios_base::in);
    CHECK_FAIL(file, "Error opening file '" + path + "'");

    path = name + DHT_EXTENSION;
    hashtable.reset(new dht::DiskHash<Location>(path.c_str(),
                                                read_max_key_size(),
                                                dht::DHOpenRO));
}

VectorDiskHash::VectorDiskHash(const std::string& name,
                               const size_t max_key_size) {
    auto path = name + DATA_EXTENSION;
    auto mode = std::ios_base::binary | std::ios_base::trunc | 
                std::ios_base::in | std::ios_base::out;
    file.open(path, mode);
    CHECK_FAIL(file, "Error opening file '" + path + "'");

    path = name + DHT_EXTENSION;
    hashtable.reset(new dht::DiskHash<Location>(path.c_str(),
                                                max_key_size,
                                                dht::DHOpenRW));
    // `max_key_size` must be persisted to reopen `hashtable`.
    write_max_key_size(max_key_size);
}

void VectorDiskHash::reserve(const std::string& key, const size_t space) {
    if (hashtable->lookup(key.c_str()) != NULL) {
        std::string msg = "Cannot reserve space for existing key '" + key + "'";
        throw std::runtime_error(msg);
    }

    // Insert the new key at the end of the file.
    file.seekp(0, std::ios_base::end);
    file.put(KEY_DELIMITER);
    // Add the key to the hashtable.
    // loc.start should point to right after the KEY_DELIMITER.
    Location loc{file.tellp(), file.tellp(), space};
    hashtable->insert(key.c_str(), loc);
    // Unused reserved space is represented with KEY_DELIMITERs.
    std::fill_n(std::ostream_iterator<char>(file), space, KEY_DELIMITER);
}

void VectorDiskHash::insert(const std::string& key, const std::string& value) {
    const char* value_str = value.c_str();
    const size_t value_size = value.size();

    if (key.compare(end_of_file_key)) {
        // Keys that aren't at the end of the file must a) be new keys
        // or b) be reserved keys.
        Location* loc = hashtable->lookup(key.c_str());

        if (loc == NULL) {
            // `key` is a new key
            // Insert it at the end of the file.
            file.seekp(0, std::ios_base::end);
            file.put(KEY_DELIMITER);
            // Add the key to the hashtable .
            // loc should point to right after KEY_DELIMITER.
            Location loc{file.tellp(), 0, 0};
            hashtable->insert(key.c_str(), loc);
            // Add the new value.
            file.write(value_str, value_size);
            end_of_file_key = key;
            return;
        } else {
            // `key` is a previously reserved key
            // +1 accounts for the length of VALUE_DELIMITER
            const size_t write_size = value_size + 1;
            if (write_size > loc->reserve_remaining) {
                std::string msg = "Insufficient space reserved for key '" + key + "'";
                throw std::runtime_error(msg);
            }
            // Write the new value to the file. Sometimes, we add an
            // extra leading delimiter, but it doesn't matter.
            file.seekp(loc->write_location);
            file.put(VALUE_DELIMITER);
            file.write(value_str, value_size);
            // Update reservation information for this key.
            loc->write_location = file.tellp();
            loc->reserve_remaining -= write_size;
        }
    } else {
        // The last key in the file can always be easily updated.
        // It's already in the hashtable, so no need to add it.
        // Just write the new value to the file.
        file.seekp(0, std::ios_base::end);
        file.put(VALUE_DELIMITER);
        file.write(value_str, value_size);
    }
}

void VectorDiskHash::write_max_key_size(const size_t max_key_size) {
    auto s = std::to_string(max_key_size);
    file.write(s.c_str(), s.size());
    file.put(KEY_DELIMITER);
}

size_t VectorDiskHash::read_max_key_size() {
    std::string line;
    getline(file, line, KEY_DELIMITER);
    try {
        return std::stoi(line);
    } catch (std::invalid_argument& e) {
        std::string msg = "Possible corruption: incorrect table format";
        throw std::runtime_error(msg);
    }
}

std::vector<std::string> VectorDiskHash::get(const std::string& key) {
    Location* loc = hashtable->lookup(key.c_str());
    if (loc == NULL) {
        throw std::runtime_error("Nonexistent key '"+ key + "'");
    }

    file.seekg(loc->start);
    std::string line;
    getline(file, line, KEY_DELIMITER);
    CHECK_FAIL(file, "Failed to get key '"+ key + "'");
    // Since insert() sometimes adds extra leading delimiters, it is
    // important that split_str ignore them (which it does).
    return split_str(line, VALUE_DELIMITER);
}

bool GeneticDataValidator::validate(const std::string& line) {
    size_t end = 0;
    size_t start = line.find_first_not_of(delimiter, end);
    size_t col = 0;
    // Being lazy is fast, so process as few columns as possible.
    while (start != std::string::npos && col <= last_important_column) {
        end = line.find(delimiter, start);

        if (col == snp_a_index) {
            data.snp_a = line.substr(start, end - start);
            // Overlong snp_a strings are invalid.
            if (data.snp_a.size() > max_key_size) {
                return false;
            }
        } else if (col == snp_b_index) {
            data.snp_b = line.substr(start, end - start);
        } else if (col == maf_index) {
            // Non-double maf values are invalid.
            try {
                data.maf = std::stod(line.substr(start, end - start));
            } catch (std::invalid_argument& e) {
                return false;
            }
        } else if (col == r2_index) {
            // Non-double r2 values or r2 values below min_r2
            // are invalid.
            try {
                data.r2 = std::stod(line.substr(start, end - start));
                if (data.r2 < min_r2) {
                    return false;
                }
            } catch (std::invalid_argument& e) {
                return false;
            }
        }

        start = line.find_first_not_of(delimiter, end);
        col++;
    }

    // If `line` does not contain columns for every required
    // value, it is invalid.
    return col > last_important_column;
}

LDLookup::LDLookup(const std::string& name) {
    ld_pairs.reset(new VectorDiskHash(name + LD_PAIR_EXTENSION));
    bins.reset(new VectorDiskHash(name + BINS_EXTENSION));
}

LDLookup::LDLookup(const std::string& name,
                   const std::string& source_path,
                   GeneticDataValidator validator) {
    auto ld_pairs_path = name + LD_PAIR_EXTENSION;
    ld_pairs.reset(new VectorDiskHash(ld_pairs_path,
                                      validator.max_key_size));

    auto bins_path = name + BINS_EXTENSION;
    bins.reset(new VectorDiskHash(bins_path, MAX_BIN_KEY_SIZE));

    std::fstream data(source_path, std::ios_base::in);
    CHECK_FAIL(data, "Error opening file '" + source_path + "'");

    // On our first pass through the data, populate the ld_pair table.
    // Track the number of genetic markers associated with particular
    // bins. A bin is identified by a number of LD surrogates and a
    // range of MAF values.
    std::map<std::string, size_t> bin_sizes;
    GeneticDataValidator::ProcessedData last_key;
    std::string line;
    size_t surrogate_count;
    bool key_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        // A new key has been encountered, so store bin information
        // for the last key.
        if (last_key.snp_a.compare(validator.data.snp_a)) {
            if (key_found) {
                auto key = to_bin_key(surrogate_count, last_key.maf);
                auto emplace_pair(bin_sizes.emplace(key, 0));
                emplace_pair.first->second += last_key.snp_a.size() + 1;
            } else {
                key_found = true;
            }

            last_key = validator.data;
            surrogate_count = 0;
        }

        // Populate the ld_pair table.
        ld_pairs->insert(validator.data.snp_a, validator.data.snp_b);
        // Update bin information.
        surrogate_count++;
    }

    // Reserve bin space.
    for (auto it = bin_sizes.begin(); it != bin_sizes.end(); it++) {
        bins->reserve(it->first, it->second);
    }

    // On our second pass through the data, populate the bin table.
    data.clear();
    data.seekg(0, std::ios::beg);
    last_key = {};
    surrogate_count = 0;
    key_found = false;

    while (data) {
        getline(data, line);
        if (!validator.validate(line)) {
            continue;
        }

        // A new key has been encountered, so write it to the bin table.
        if (last_key.snp_a.compare(validator.data.snp_a)) {
            if (key_found) {
                auto key = to_bin_key(surrogate_count, last_key.maf);
                bins->insert(key, last_key.snp_a);
            } else {
                key_found = true;
            }

            last_key = validator.data;
            surrogate_count = 0;
        }

        // Update bin information.
        surrogate_count++;
    }
}

std::vector<std::string> LDLookup::find_ld(const std::string& key) {
    return ld_pairs->get(key);
}

std::vector<std::string> LDLookup::find_similar(const int surrogate_count,
                                                const double maf) {
    try {
        return bins->get(to_bin_key(surrogate_count, maf));
    } catch (std::runtime_error& e) {
        return std::vector<std::string>();
    }
}

std::string LDLookup::to_bin_key(const int surrogate_count,
                                 const double maf) {
    const int maf_bin = floor(256 * log2(maf));
    return std::to_string(surrogate_count) + " " + std::to_string(maf_bin);
}
