#ifndef LDLOOKUP_VDH_HPP
#define LDLOOKUP_VDH_HPP

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"
#include "global.hpp"

/* Splits `s` at a `delimiter`, ignoring consecutive/leading/trailing delimiters. */
std::vector<std::string> split_str(const std::string& s, const char delimiter);

/* Persistent, write-once map of strings to arrays of strings. */
class VectorDiskHash {
    public:
        /** Opens an existing VectorDiskHash.
         * 
         * Throws std::runtime_error on failure or if the VectorDiskHash
         * does not exist.
        */
        VectorDiskHash(const std::string& name);

        /* Creates a new VectorDiskHash. */
        VectorDiskHash(const std::string& name, const size_t max_key_size);

        /** Pre-allocates space for values associated with a key.
         * 
         * All values associated with a particular key must be
         * insert()ed consecutively. Alternatively, space can
         * be reserve()d for those values, and they can be inserted
         * as long as space remains.
         * 
         * `space` should be as large as the sum of size()s of all values
         * associated with `key` plus one byte for each value to account
         * for overhead.
         * 
         * Throws std::runtime_error if key already exists or if the
         * VectorDiskHash was opened from existing files rather than
         * created.
         */
        void reserve(const std::string& key, const size_t space);

        /** Associates a value with a key.
         * 
         * Throws std::runtime_error if the VectorDiskHash was opened
         * from existing files rather than created.
         */
		void insert(const std::string& key, const std::string& value);

        /** Retrieves values associated with a key. 
         * 
         * Throws std::runtime_error if key does not exist or
         * cannot be retrieved.
        */
		std::vector<std::string> get(const std::string& key);

    private:
        inline static const char KEY_DELIMITER = '\n';
        inline static const char VALUE_DELIMITER = '\t';
        inline static const std::string DATA_EXTENSION = ".mvdhdat";
        inline static const std::string DHT_EXTENSION = ".mvdhdht";

        struct Location {
            std::streampos start;
            std::streampos write_location;
            size_t reserve_remaining;
        };

        std::shared_ptr<dht::DiskHash<Location>> hashtable;
        std::fstream file;
        std::string end_of_file_key;


        void write_max_key_size(const size_t max_key_size);
        size_t read_max_key_size();
};

#endif
