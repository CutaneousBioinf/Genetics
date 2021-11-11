#ifndef LDLOOKUP_LDLOOKUP_HPP
#define LDLOOKUP_LDLOOKUP_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"

#define CHECK_FAIL(file, msg) if ((file).fail()) { throw std::runtime_error((msg)); }

// A general note: when operations fail, these calsses/functions usually throw
// std::runtime_error.

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

/** Provides fast line-oriented parsing of genetic records.*/
class GeneticDataValidator {
	public:
		struct ProcessedData {
			std::string snp_a;
			std::string snp_b;
			double r2;
			double maf;
		};

		char delimiter;
		size_t snp_a_index;
		size_t snp_b_index;
		size_t r2_index;
		size_t maf_index;
		double min_r2;
        size_t max_key_size;
		ProcessedData data;

		GeneticDataValidator(const char delimiter,
                             const size_t snp_a_index,
                             const size_t snp_b_index,
                             const size_t r2_index,
                             const size_t maf_index,
                             const double min_r2,
                             const size_t max_key_size) :
                        delimiter(delimiter),
                        snp_a_index(snp_a_index),
                        snp_b_index(snp_b_index),
                        r2_index(r2_index),
                        maf_index(maf_index),
                        min_r2(min_r2),
                        max_key_size(max_key_size),
                        last_important_column(std::max({ snp_a_index,
                                                         snp_b_index,
                                                         maf_index,
                                                         r2_index }))
                        {}

        /** Parses one line of genetic data into C++ types.
         * 
         * If `line` is valid, the `data` member will contain
         * the parsed data from `line` after this call. If
         * `line` is invalid, the `data` member may contain
         * junk.
         */
		bool validate(const std::string& line);

	private:
		size_t last_important_column;
};


/** Provides interface to genetic data. */
class LDLookup {
	public:
        /* Opens an existing LDLookup. */
		LDLookup(const std::string& name);

        /* Creates a new LDLookup. */
		LDLookup(const std::string& name,
                 const std::string& source_path,
                 GeneticDataValidator validator);

        /* Finds markers in linkage disequilibrium with a key marker. */
		std::vector<std::string> find_ld(const std::string& key);

        /* Find markers with similar MAF and number of LD surrogates. */
		std::vector<std::string> find_similar(const int surrogate_count, const double maf);
	
	private:
		inline static const std::string LD_PAIR_EXTENSION = "_ld_pairs";
        inline static const std::string BINS_EXTENSION = "_bins";
        static const size_t MAX_BIN_KEY_SIZE = 32;

        std::shared_ptr<VectorDiskHash> ld_pairs;
        std::shared_ptr<VectorDiskHash> bins;
		
        /* Converts an MAF/LD surrogates pair into a key for `bins`. */
		std::string to_bin_key(const int surrogate_count, const double maf);
};

#endif
