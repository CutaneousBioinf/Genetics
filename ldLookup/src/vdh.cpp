#include <iterator>  // std::ostream_iterator
#include <string_view>

#include "string_ops.hpp"
#include "vdh.hpp"

using std::string;
using std::vector;

VectorDiskHash::VectorDiskHash(const Options &opts)
    : options(opts), eof_key("") {
	if (options.max_key_size) {
		options.max_key_size++;
	}
	
	// Set up 'file'.
	auto file_mask = std::ios_base::binary | std::ios_base::in;
	if (options.create) {
		// Ensure 'file_path' doesn't already exist.
		if (open_file(options.file_path, file_mask)) {
			throw vdh_mode_error(options, "File Already Exists");
		}

		// Update file_mask to create 'file'.
		file_mask |= std::ios_base::out | std::ios_base::trunc;
	}

	// Create/open 'file'.
	if (!open_file(options.file_path, file_mask)) {
		throw vdh_internal_error(options, "Failed to Open File");
	}

	// 'file' should throw on failure to simplify error checking.
	file.exceptions(std::fstream::badbit | std::fstream::failbit);

	if (options.create) {
		// Persist max_key_size.
		string mks = std::to_string(options.max_key_size);
		file.write(mks.c_str(), mks.size());
		file.put(KEY_DELIMITER);
	} else if (!options.max_key_size) {
		// Load the persisted max_key_size. We need it to open 'table'.
		string mks;
        std::getline(file, mks, KEY_DELIMITER);
		try {
			options.max_key_size = static_cast<size_t>(std::stoull(mks));
		} catch (std::invalid_argument &e) {
			throw vdh_internal_error(options, "Failed to Read Max Key Size");
		}
	}

	// Set up 'table'.
	auto dht_mask = dht::DHOpenRO;
	if (options.create) {
		// Ensure 'table_path' doesn't exist.
		if (open_table(options.table_path, options.max_key_size, dht_mask)) {
			throw vdh_mode_error(options, "Table Already Exists");
		}

		// Update dht_mask to create 'table'.
		dht_mask = dht::DHOpenRW;
	}

	// Create/open 'table'.
	if (!open_table(options.table_path, options.max_key_size, dht_mask)) {
		throw vdh_internal_error(options, "Failed to Open Table");
	}
}

void VectorDiskHash::append(
    const string &key,
    const vector<string> &values) {
	// Begin Function Body
	if (!options.create) {
		string msg = "append(): VectorDiskHash is Read-Only";
		throw vdh_mode_error(options, msg);
	} else if (key.size() > options.max_key_size) {
		string msg = "append(): Key Too Long - " + key;
		throw vdh_key_error(options, msg);
	}

	// Serialize the values vector.
	string serialized = join(values, VALUE_DELIMITER, true);
	size_t serialized_size = serialized.size();

	if (key == eof_key) {
		file.seekp(0, std::ios_base::end);
		file.put(VALUE_DELIMITER);
		file.write(serialized.data(), serialized_size);
	} else {
		Location *loc = table->lookup(key.c_str());
		if (loc == nullptr) {
			file.seekp(0, std::ios_base::end);
			file.put(KEY_DELIMITER);

			Location new_loc;
			new_loc.start = file.tellp();
			new_loc.bytes_reserved = 0;

			file.write(serialized.data(), serialized_size);
			new_loc.write_location = file.tellp();

			table->insert(key.c_str(), new_loc);

			eof_key = key;
		} else if (serialized_size > loc->bytes_reserved) {
			string msg = "append(): Key Out of Reserved Space - " + key;
			throw vdh_value_error(options, msg);
		} else {
			file.seekp(loc->write_location);
			file.write(serialized.data(), serialized_size);
			loc->bytes_reserved -= serialized_size;
			loc->write_location = file.tellp();
		}
	}
}

void VectorDiskHash::append(const string &key, const string &values) {
	append(key, vector<string>({values}));
}

Options VectorDiskHash::get_options() const {
	return options;
}

bool VectorDiskHash::is_member(const string &key) const {
	return table->is_member(key.c_str());
}

vector<string> VectorDiskHash::lookup(const string &key) {
	string serialized(read_serialized_vector(key));
	return split(serialized, VALUE_DELIMITER);
}

vector<string> VectorDiskHash::lookup_sample(const string &key, size_t k) {
	string serialized(read_serialized_vector(key));
	return sample(serialized, VALUE_DELIMITER, k);
}

void VectorDiskHash::reserve(const string &key, size_t bytes_to_reserve) {
	if (!options.create) {
		string msg = "reserve(): VectorDiskHash is Read-Only";
		throw vdh_mode_error(options, msg);
	} else if (key.size() > options.max_key_size) {
		string msg = "reserve(): Key Too Long - " + key;
		throw vdh_key_error(options, msg);
	} else if (is_member(key)) {
		string msg = "reserve(): Key is_member - " + key;
		throw vdh_key_error(options, msg);
	}

	file.seekp(0, std::ios_base::end);
	file.put(KEY_DELIMITER);

	Location loc;
	loc.start = file.tellp();
	loc.write_location = file.tellp();
	loc.bytes_reserved = bytes_to_reserve;

	table->insert(key.c_str(), loc);
	std::fill_n(
	    std::ostream_iterator<char>(file),
	    loc.bytes_reserved,
	    KEY_DELIMITER);
}

bool VectorDiskHash::open_file(
    const string &file_path,
    std::_Ios_Openmode mask) {
	// Open file with the specified mask.
	file.open(file_path, mask);

	// Return whether there were any errors.
	return file.good();
}

bool VectorDiskHash::open_table(
    const string &table_path,
    const size_t max_key_size,
    dht::OpenMode mask) {
	// Open table with the specified settings and return
	// whether there were any errors.
	try {
		table.reset(new dht::DiskHash<Location>(
		    table_path.c_str(), max_key_size, mask));
		return true;
	} catch (std::exception &e) {
		return false;
	}
}

string VectorDiskHash::read_serialized_vector(const string &key) {
	Location *loc = table->lookup(key.c_str());
	// Check that key is in the VectorDiskHash.
	if (loc == nullptr) {
		string msg = "lookup(): Nonexistent Key - " + key;
		throw vdh_key_error(options, msg);
	}

	file.seekg(loc->start);
	string serialized;
	std::getline(file, serialized, KEY_DELIMITER);
	return serialized;
}