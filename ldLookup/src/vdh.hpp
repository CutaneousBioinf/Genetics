#ifndef _LDLOOKUP_VDH_HPP_
#define _LDLOOKUP_VDH_HPP_

#include <stddef.h>  // size_t

#include <fstream>    // std::fstream, streampos, _Ios_openmode
#include <memory>     // std::shared_ptr
#include <stdexcept>  // std::runtime_error
#include <string>
#include <vector>

#include "diskhash/src/diskhash.hpp"

/**
 * VectorDiskHashes are made of two files, the paths to which must be
 * specified to create or open the VectorDiskHash. file_path and
 * table_path are NOT interchangeable.
 * 
 * If 'create' is true:
 * - Both files will be created.
 * - If either file already exists, vdh_mode_error will be thrown.
 * If 'create' is false:
 * - If either file does not exist, vdh_mode_error will be thrown.
 */
struct Options {
	std::string file_path;
	std::string table_path;
    unsigned long max_key_size;
	bool create;
};

/* Custom Exceptions for VectorDiskHash */
struct vdh_error : std::runtime_error {
	vdh_error(const Options &opts, const std::string &msg)
	    : std::runtime_error(
	          "VectorDiskHash Error\nFile Path: " +
			  opts.file_path + 
			  "\nTable Path: " + 
			  opts.table_path + 
			  "\nInfo: " +
			  msg) {}
};

struct vdh_internal_error : vdh_error {
	vdh_internal_error(
	    const Options &opts,
	    const std::string &msg = "Internal Error")
	    : vdh_error(opts, msg) {}
};

struct vdh_mode_error : vdh_error {
	vdh_mode_error(
	    const Options &opts,
	    const std::string &msg = "Mode Error")
	    : vdh_error(opts, msg) {}
};

struct vdh_key_error : vdh_error {
	vdh_key_error(
	    const Options &opts,
	    const std::string &msg = "Key Error")
	    : vdh_error(opts, msg) {}
};

struct vdh_value_error : vdh_error {
	vdh_value_error(
	    const Options &opts,
	    const std::string &msg = "Value Error")
	    : vdh_error(opts, msg) {}
};

/* Persistent, write-once map from strings to string vectors. */
class VectorDiskHash {
   public:
	/**
     * EFFECTS: Creates or opens a VectorDiskHash.
     * THROWS: vdh_mode_error (see Options class).
     *         exception on operation failure.
     * NOTE: Also see documentation for the Options class.
     */
	VectorDiskHash(const Options &opts);

    /**
     * EFFECTS: Associates 'key' with 'values'.
     * THROWS: vdh_mode_error if VectorDiskHash is read-only.
     *         vdh_value_error if 'values' overruns available space for 'key'
     *                      OR if a string in 'values' contains RESERVED_CHAR.
     *         exception on operation failure.
     */
	void append(
	    const std::string &key,
	    const std::vector<std::string> &values);

	/**
     * EFFECTS: Equivalent to append(key, std::vector<std::string>({values}))
     */
	void append(const std::string &key, const std::string &values);

	/**
     * EFFECTS: Returns options used to create the VectorDiskHash.
     */
	Options get_options() const;

    /**
     * EFFECTS: Returns whether reserve(key...) or append(key...) were called.
     */
    bool is_member(const std::string& key) const;

    /**
     * EFFECTS: Retrieves values associated with 'key'.
     * THROWS: vdh_key_error if !is_member(key).
     *         exception on operation failure.
     */
	std::vector<std::string> lookup(const std::string &key);

    /**
     * EFFECTS: Randomly samples 'k' values associated with 'key'
     *          (with replacement).
     * THROWS: vdh_key_error if !is_member(key).
     *         exception on operation failure.
     */
	std::vector<std::string> lookup_sample(const std::string &key, size_t k);

	/**
     * EFFECTS: Allocates fixed amount of storage for values associated with
     *          'key'.
     * THROWS: vdh_mode_error if VectorDiskHash is read-only.
     *         vdh_key_error if key.size() > get_options().max_key_size
     *                    OR if is_member(key).
     *         exception on operation failure.
     */
	void reserve(const std::string &key, size_t bytes_to_reserve);

   private:
	const char KEY_DELIMITER = '\n';
	const char VALUE_DELIMITER = '\t';

    /* Stores location of a serialized vector of strings in 'file'. */
	struct Location {
		std::streampos start;
		std::streampos write_location;
		size_t bytes_reserved;
	};

	/* Stores options used to create the SDH. */
	Options options;

	/* Stores serialized values. */
	std::fstream file;

	/* Maps keys to the locations of serialized values in 'file'. */
	std::shared_ptr<dht::DiskHash<Location>> table;
	
	/* Stores the newest non-reserve()d key. */
	std::string eof_key;

	bool open_file(const std::string &file_path, std::_Ios_Openmode mask);
	bool open_table(
	    const std::string &open_table,
	    const size_t max_key_size,
	    dht::OpenMode mask);
	std::string read_serialized_vector(const std::string& key);
};

#endif