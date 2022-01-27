#include <iterator> // std::fill_n

#include "vdh.hpp"

std::vector<std::string> split_str(const std::string& s, char delimiter) {
    std::vector<std::string> ret;
    ret.reserve(s.length() / 2); // for performance reasons
    size_t end = 0;
    size_t start = s.find_first_not_of(delimiter, end);
    while (start != std::string::npos) {
	    end = s.find(delimiter, start);
	    ret.emplace_back(s.substr(start, end - start));
        start = s.find_first_not_of(delimiter, end);
    }
    return ret;
}

VectorDiskHash::VectorDiskHash(const std::string& name) {
    auto data_path = name + DATA_EXTENSION;
    auto mode = std::ios_base::binary | std::ios_base::in;
    file.open(data_path, mode);
    CHECK_FAIL(file, "Error Opening File '" + data_path + "'");

    this->name = name;
    this->max_key_size = read_max_key_size();

    auto dht_path = name + DHT_EXTENSION;
    hashtable.reset(new dht::DiskHash<Location>(
        dht_path.c_str(),
        max_key_size,
        dht::DHOpenRO
    ));
}

VectorDiskHash::VectorDiskHash(
    const std::string& name, const size_t max_key_size
) {
    auto data_path = name + DATA_EXTENSION;
    auto mode = std::ios_base::binary
              | std::ios_base::trunc 
              | std::ios_base::in 
              | std::ios_base::out;
    file.open(data_path, mode);
    CHECK_FAIL(
        file, 
        "VectorDiskHash '" + name + "': Error Opening File '" + data_path + "'"
    );

    this->name = name;
    this->max_key_size = max_key_size;

    write_max_key_size(max_key_size);

    auto dht_path = name + DHT_EXTENSION;
    hashtable.reset(new dht::DiskHash<Location>(
        dht_path.c_str(),
        max_key_size,
        dht::DHOpenRW
    ));
}

void VectorDiskHash::reserve(const std::string& key, const size_t space) {
    if (hashtable->lookup(key.c_str()) != NULL) {
        throw std::runtime_error(
            "VectorDiskHash '" + name 
            + "': Cannot reserve() Existing Key '" + key + "'"
        );
    }

    // Reserve space at the end of the file for the new key.
    // Write KEY_DELIMITER to indicate a new key.
    file.seekp(0, std::ios_base::end);
    file.put(KEY_DELIMITER);
    // Create the hashtable entry for key.
    // Do this here because loc.start must point to the start of the
    // reserved space.
    Location loc{file.tellp(), file.tellp(), space};
    hashtable->insert(key.c_str(), loc);
    // Write KEY_DELIMITER characters representing unused reserved space.
    std::fill_n(std::ostream_iterator<char>(file), space, KEY_DELIMITER);
}

void VectorDiskHash::insert(const std::string& key, const std::string& value) {
    const char* value_str = value.c_str();
    const size_t value_size = value.size();

    if (key.compare(end_of_file_key)) {
        // Key must be reserved or totally new.
        Location* loc = hashtable->lookup(key.c_str());

        if (loc == NULL) {
            // Key is new.
            // Insert it at the end of the file.
            file.seekp(0, std::ios_base::end);
            file.put(KEY_DELIMITER);
            // Create the hashtable entry for key.
            // Do this here because loc.start must point to the start
            // of the key's values.
            Location start_loc{file.tellp(), 0, 0};
            hashtable->insert(key.c_str(), start_loc);
            // Insert the new value.
            file.write(value_str, value_size);
            end_of_file_key = key;
        } else {
            // Key is a previously reserved key.
            const size_t write_size = value_size + sizeof(VALUE_DELIMITER);
            if (write_size > loc->reserve_remaining) {
                throw std::runtime_error(
                    "VectorDiskHash '" + name +
                    "': Out of reserve()d Space for Key '" + key + "'"
                );
            }
            // Insert the new value.
            // The first value associated with a key will have an
            // extra leading delimiter, but it doesn't matter.
            file.seekp(loc->write_location);
            file.put(VALUE_DELIMITER);
            file.write(value_str, value_size);
            // Update reservation information for this key.
            loc->write_location = file.tellp();
            loc->reserve_remaining -= write_size;
        }
    } else {
        // The key is the last key in the file, so it's already
        // in the hashtable.
        // Inser the new value.
        file.seekp(0, std::ios_base::end);
        file.put(VALUE_DELIMITER);
        file.write(value_str, value_size);
    }

    CHECK_FAIL(
        file,
        "VectorDiskHash '" + name + "': Write Failed"
    );
}

std::vector<std::string> VectorDiskHash::get(const std::string& key) {
    Location* loc = hashtable->lookup(key.c_str());
    if (loc == NULL) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Nonexistent Key '" + key + "'"
        );
    }

    file.seekg(loc->start);
    std::string line;
    read_until_key_delimiter(line);
    // split_str ignores leading delimiters.
    // This is important because insert() can produce extra leading
    // delimiters.
    return split_str(line, VALUE_DELIMITER);
}

size_t VectorDiskHash::read_max_key_size() {
    std::string line;
    read_until_key_delimiter(line);
    try {
        return std::stoi(line);
    } catch (std::invalid_argument& e) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Invalid Max Key Size"
        );
    }
}

void VectorDiskHash::write_max_key_size(const size_t max_key_size) {
    auto mks = std::to_string(max_key_size);
    file.write(mks.c_str(), mks.size());
    file.put(KEY_DELIMITER);
}

void VectorDiskHash::read_until_key_delimiter(std::string& buff) {
    if (!getline(file, buff, KEY_DELIMITER)) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Read Failed"
        );
    }
}

std::string VectorDiskHash::get_name() {
    return name;
}

size_t VectorDiskHash::get_max_key_size() {
    return max_key_size;
}
