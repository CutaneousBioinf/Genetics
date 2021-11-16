#ifndef LDLOOKUP_LDLOOKUP_HPP
#define LDLOOKUP_LDLOOKUP_HPP

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "diskhash/diskhash.hpp"
#include "geneticDataValidator.hpp"
#include "global.hpp"
#include "vdh.hpp"

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
