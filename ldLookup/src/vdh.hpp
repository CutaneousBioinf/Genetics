#ifndef LDLOOKUP_VDH_HPP
#define LDLOOKUP_VDH_HPP

#include <fstream>
#include <memory> // std::shared_ptr
#include <string>
#include <stdexcept> // std::runtime_error
#include <vector>

#include "diskhash/diskhash.hpp"

struct VectorDiskHashError : public std::runtime_error {
    VectorDiskHashError(
        const std::string& name,
        const std::string& what = "No Details"
    ) : std::runtime_error("VectorDiskHash Error in '" + name + "': " + what) {}
};

struct InternalError : public VectorDiskHashError {
    InternalError(
        const std::string& name,
        const std::string& what = "No Details"
    ) : VectorDiskHashError(name, "Internal Error - " + what) {}
};

struct KeyError : public VectorDiskHashError {
    KeyError(
        const std::string& name,
        const std::string& what = "No Details"
    ) : VectorDiskHashError(name, "Key Error - " + what) {}
};

struct ModeError : public VectorDiskHashError {
    ModeError(
        const std::string& name,
        const std::string& what = "No Details"
    ) : VectorDiskHashError(name, "Mode Error - " + what) {}
};

/* Persistent, write-once map of strings to arrays of strings. */
class VectorDiskHash {
    public:
        /**
         * EFFECTS: Opens an existing VectorDiskHash.
         * THROWS: internal_error on failure.
         */
        VectorDiskHash(const std::string& name, const std::string& dir="");

        /**
         * EFFECTS: Creates a new VectorDiskHash.
         * THROWS: mode_error if the VectorDiskHash already exists.
         *         internal_error on failure.
         */
        VectorDiskHash(
            const std::string& name,
            const size_t max_key_size,
            const std::string& dir=""
        );

        /**
         * EFFECTS: Pre-allocates space for values associated with a key.
         * THROWS: key_error if key.size() > get_max_key_size().
         *         key_error if key was already insert()ed or reserve()d.
         *         mode_error if the VectorDiskHash was opened from existing
         *           files.
         *         internal_error on failure.
         * NOTE:
         * All values associated with a particular key must be
         * insert()ed consecutively. Alternatively, space can
         * be reserve()d for those values, and they can be inserted
         * as long as space remains.
         * 
         * `space` should be as large as the sum of .size()s of all values
         * associated with `key` plus one byte for each value (storage
         * overhead).
         */
        void reserve(const std::string& key, const size_t space);

        /**
         * EFFECTS: Associates a value with a key.
         * THROWS: key_error if key.size() exceeds get_max_key_size().
         *         key_error if key is reserved and value.size()+1 exceeds the
         *           amount of reserved space.
         *         mode_error if the VectorDiskHash was opened from existing
         *           files.
         *         internal_error on failure.
         */
		void insert(const std::string& key, const std::string& value);

        /**
         * EFFECTS: Retrieves all values associated with a key. 
         * THROWS: key_error if key does not exist.
         *         internal_error on failure.
        */
		std::vector<std::string> get(const std::string& key);

        /**
         * EFFECTS: Randomly samples n_random values associated with a key
         * (with replacement).
         * THROWS: key_error if key does not exist.
         *         internal_error on failure.
        */
        std::vector<std::string> get_random(const std::string& key, size_t n_random);

        std::string get_name();

        size_t get_max_key_size();

    private:
        inline static const char KEY_DELIMITER = '\n';
        inline static const char VALUE_DELIMITER = '\t';
        inline static const std::string DATA_EXTENSION = ".vdhdat";
        inline static const std::string DHT_EXTENSION = ".vdhdht";

        struct Location {
            std::streampos start;
            std::streampos write_location;
            size_t reserve_remaining;
        };

        std::shared_ptr<dht::DiskHash<Location>> hashtable;
        std::fstream file;
        std::string end_of_file_key;
        std::string name;
        const bool read_only;
        size_t max_key_size;

        size_t read_max_key_size();
        void write_max_key_size(const size_t max_key_size);
        void read_until_key_delimiter(std::string& buff);
};

#endif
