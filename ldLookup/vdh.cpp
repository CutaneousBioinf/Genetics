#include <iterator>

#include "vdh.hpp"

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
