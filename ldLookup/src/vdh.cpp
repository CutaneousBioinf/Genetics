#include <iterator> // std::fill_n
#include <random> // std::mt19337, std::uniform_int_distribution, std::random_device
#include <string_view> // std::string_view

#include "split.hpp"
#include "vdh.hpp"

VectorDiskHash::VectorDiskHash(const std::string& name)
: name(name) {
    // Open the data file.
    auto data_path(name + DATA_EXTENSION);
    auto mode(std::ios_base::binary | std::ios_base::in);
    file.open(data_path, mode);
    if (!file) {
        throw std::runtime_error(
            "VectorDiskHash '" + name +
            "': Failed to Open '" + data_path + "'"
        );
    }
    
    // Read max_key_size from the data file.
    this->max_key_size = read_max_key_size();

    // Open the hashtable.
    auto dht_path(name + DHT_EXTENSION);
    try {
        hashtable.reset(new dht::DiskHash<Location>(
            dht_path.c_str(), max_key_size, dht::DHOpenRO
        ));
    } catch (std::exception &e) {
        throw std::runtime_error(
            "VectorDiskHash '" + name +
            "': Failed to Open '" + dht_path + "'"
        );
    }
}

VectorDiskHash::VectorDiskHash(
    const std::string& name, const size_t max_key_size
) : name(name), max_key_size(max_key_size) {
    // Open the data file.
    auto data_path = name + DATA_EXTENSION;
    auto mode = std::ios_base::binary
              | std::ios_base::trunc 
              | std::ios_base::in 
              | std::ios_base::out;
    file.open(data_path, mode);
    if (!file) {
        throw std::runtime_error(
            "VectorDiskHash '" + name +
            "': Failed to Open '" + data_path + "'"
        );
    }

    // Save max_key_size to the data file.
    write_max_key_size(max_key_size);

    // Open the hashtable.
    auto dht_path = name + DHT_EXTENSION;
    try {
        hashtable.reset(new dht::DiskHash<Location>(
            dht_path.c_str(), max_key_size, dht::DHOpenRW
        ));
    } catch (std::exception &e) {
        throw std::runtime_error(
            "VectorDiskHash '" + name +
            "': Failed to Open '" + dht_path + "'"
        );
    }
}

void VectorDiskHash::reserve(const std::string& key, const size_t space) {
    // Prevent reserve() calls on existing keys.
    if (hashtable->lookup(key.c_str()) != NULL) {
        throw std::runtime_error(
            "VectorDiskHash '" + name 
            + "': Cannot reserve() Existing Key '" + key + "'"
        );
    }

    // Reserve space at the end of the file for the new key's values.
    // First, write KEY_DELIMITER to begin a new set of values.
    file.seekp(0, std::ios_base::end);
    file.put(KEY_DELIMITER);
    // We need to store the location of the reserved space in the hashtable.
    // Since loc.start must be the start of the reserved space, we construct
    // loc before writing characters to represent the unused reserved space.
    Location loc{ file.tellp(), file.tellp(), space} ;
    hashtable->insert(key.c_str(), loc);
    // Write KEY_DELIMITER characters representing the unused reserved space.
    std::fill_n(std::ostream_iterator<char>(file), space, KEY_DELIMITER);
}

void VectorDiskHash::insert(const std::string& key, const std::string& value) {
    const char* value_str(value.c_str());
    const size_t value_size(value.size());

    if (key.compare(end_of_file_key)) {
        // key != end_of_file_key. Thus, key is either reserve()d or completely
        // new. We can check whether is reserve()d or new with this lookup.
        Location* loc = hashtable->lookup(key.c_str());

        if (loc == NULL) {
            // Key is completely new.
            // First, mark the end of the end_of_file_key with a KEY_DELIMITER.
            file.seekp(0, std::ios_base::end);
            file.put(KEY_DELIMITER);
            // Store the location of the new key's values in the hashtable.
            // Since loc.start must be the start of the values, we construct
            // loc before writing the values.
            Location start_loc{file.tellp(), 0, 0};
            hashtable->insert(key.c_str(), start_loc);
            // Write the values.
            file.write(value_str, value_size);
            // key is the new end_of_file key.
            end_of_file_key = key;
        } else {
            // Key is a previously reserved key.
            // Assert that the reserved space is large enough to hold the
            // new value.
            const size_t write_size(value_size + sizeof(VALUE_DELIMITER));
            if (write_size > loc->reserve_remaining) {
                throw std::runtime_error(
                    "VectorDiskHash '" + name +
                    "': Out of reserve()d Space for Key '" + key + "'"
                );
            }
            // Insert the new value.
            file.seekp(loc->write_location);
            file.put(VALUE_DELIMITER);
            file.write(value_str, value_size);
            // Update reservation information for this key.
            loc->write_location = file.tellp();
            loc->reserve_remaining -= write_size;
        }
    } else {
        // key == end_of_file_key, so it's already in the hashtable.
        // Insert the new value.
        file.seekp(0, std::ios_base::end);
        file.put(VALUE_DELIMITER);
        file.write(value_str, value_size);
    }

    if (!file) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Write Failed"
        );
    }
}

std::vector<std::string> VectorDiskHash::get(const std::string& key) {
    const Location* loc(hashtable->lookup(key.c_str()));
    // lookup returns NULL if the key doesn't exist.
    if (loc == NULL) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Nonexistent Key '" + key + "'"
        );
    }
    // key exists, so read its values from the file.
    file.seekg(loc->start);
    std::string line;
    read_until_key_delimiter(line);
    return split_to_strings(line, VALUE_DELIMITER);
}

std::vector<std::string> VectorDiskHash::get_random(
    const std::string& key, 
    size_t n_random
) {
    const Location* loc(hashtable->lookup(key.c_str()));
    // lookup returns NULL if the key doesn't exist.
    if (loc == NULL) {
        throw std::runtime_error(
            "VectorDiskHash '" + name + "': Nonexistent Key '" + key + "'"
        );
    }
    // key exists, so read its values from the file.
    file.seekg(loc->start);
    std::string line;
    read_until_key_delimiter(line);
    // Sample random entries from the line without constructing a bunch
    // of strings.
    std::string_view line_view(line);
    auto split_views(split_to_string_views(line_view, VALUE_DELIMITER));
    auto random_choices(sample_with_replacement(split_views, n_random));
    // Convert the randomly-selected entries from string_views to strings.
    std::vector<std::string> ret;
    ret.reserve(n_random);
    ret.insert(ret.begin(), random_choices.begin(), random_choices.end());
    return ret;
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
