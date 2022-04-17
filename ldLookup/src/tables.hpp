#ifndef _LDLOOKUP_TABLES_HPP_
#define _LDLOOKUP_TABLES_HPP_

#include <stddef.h>  // size_t

#include <memory>  // std::shared_ptr

#include "parse_variants.hpp"
#include "stratify.hpp"
#include "vdh.hpp"

/**
 * TODO: Document!
 */
class LDTable : public VectorDiskHash {
   public:
	/**
	 * TODO: Document!
	 */
	LDTable(const Options& opts);
};

/**
 * TODO: Document!
 */
class StrataTable {
   public:
	using Stratum = std::string;

	/**
	 * TODO: Document!
	 */
	StrataTable(
	    const std::string& file_path,
        const std::string& table_path,
	    Histogram<size_t> n_surrogates_strata_in,
	    Histogram<double> maf_strata_in);

    /**
	 * TODO: Document!
	 */
	StrataTable(
	    const std::string& file_path,
        const std::string& table_path);

	/**
	 * TODO: Document!
	 */
	void append(const IndexVariantSummary& summary);

	/**
	 * TODO: Document!
	 */
	Stratum get_stratum(const IndexVariantSummary& summary);

	/**
	 * TODO: Document!
	 */
	std::vector<std::string> lookup(const IndexVariantSummary& summary);

	/**
	 * TODO: Document!
	 */
	std::vector<std::string> lookup_sample(
		const IndexVariantSummary& summary,
		const size_t k);

	/**
	 * TODO: Document!
	 */
	void reserve(const Histogram<Stratum>& strata_sizes);

   private:
	const size_t MAX_KEY_SIZE = 64;
	const std::string N_SURROGATES_KEY = "__N_SURROGATES_KEY__";
	const std::string MAF_KEY = "__MAF_KEY__";

	Histogram<size_t> n_surrogates_strata;
	Histogram<double> maf_strata;
	std::shared_ptr<VectorDiskHash> table;
};

/**
 * TODO: Document!
 */
class SummaryTable {
   public:
	/**
	 * TODO: Document!
	 */
	SummaryTable(const Options& opts);

	/**
	 * TODO: Document!
	 */
	void append(const IndexVariantSummary& summary);

	/**
	 * TODO: Document!
	 */
	IndexVariantSummary lookup(const std::string& index_variant_id);

   private:
	std::shared_ptr<VectorDiskHash> table;
};

#endif