#ifndef LDLOOKUP_GENETICDATAVALIDATOR_HPP
#define LDLOOKUP_GENETICDATAVALIDATOR_HPP

#include <algorithm>
#include <string>

#include "global.hpp"

struct GeneticData {
    std::string snp_a;
    std::string snp_b;
    double r2;
    double maf;
};

/** Provides fast line-oriented parsing of genetic records.*/
class GeneticDataValidator {
	public:
		char delimiter;
		size_t snp_a_index;
		size_t snp_b_index;
		size_t r2_index;
		size_t maf_index;
		double min_r2;
        size_t max_key_size;
		GeneticData data;

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

#endif
