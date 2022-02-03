#ifndef LDLOOKUP_TABLES_HPP
#define LDLOOKUP_TABLES_HPP

#include <map> // std::map
#include <memory> // std::shared_ptr
#include <string> // std::string
#include <vector> // std::vector

#include "vdh.hpp"

// Character used to separate fields of composite keys
const char FIELD_SEPARATOR = ' ';

/* Convert a number of LD surrogates to a string */
std::string serialize_surrogates(size_t surrogate_count);

/* Convert a minor allele frequency to a string */
std::string serialize_maf(double maf);

/* Recover a number of LD surrogates from a serialize_ld string */
size_t deserialize_surrogates(const std::string& surrogate_count);

/* Recover a minor allele frequency from a serialize_maf string */
double deserialize_maf(const std::string& maf);

/**
 * Allows fuzzy lookup of SNPs by MAF and number of LD surrogates.
 * 
 * BinsTable is initialized with two sets of bin boundaries. One
 * is used to cluster SNPs by their MAF. The other is used to cluster
 * SNPs by the number of SNPs they are in LD with. Then, by querying
 * the table for an MAF or number of LD surrogates, all SNPs with
 * both a similar MAF and a similar number of LD surrogates can be found.
 **/
class BinsTable {
    public:
        /**
         * Opens an existing BinsTable.
         * 
         * Parameters:
         *  name - Name of table to open
        **/
        BinsTable(const std::string& name);

        /**
         * Creates a new BinsTable.
         * 
         * Parameters:
         *  name - Name of table to create
         *  surrogate_quantiles - Lower bounds for LD surrogate bins
         *  maf_quantiles - Lower bounds for MAF bins
        **/
        BinsTable(
            const std::string& name,
            const std::vector<size_t>& surrogate_quantiles,
            const std::vector<double>& maf_quantiles
        );

        /**
         * Pre-allocates space for values in a bin.
         * 
         * Reserve slightly more space than is necessary--one byte per value,
         * to be exact.
         * 
         * Parameters:
         *  bin - String from BinsTable.bin()
         *  space - Space to reserve, in bytes, for bin
        **/
        void reserve(
            const std::string& bin,
            size_t space
        );

        /**
         * Associates an SNP with a bin.
         * 
         * Parameters:
         *  bin - String from BinsTable.bin()
         *  snp - SNP to place in bin
        **/
        void insert(
            const std::string& bin,
            const std::string& snp
        );

        /**
         * Retrives all SNPs associated with a bin.
         * 
         * Parameters:
         *  bin - String from BinsTable.bin()
        **/
        std::vector<std::string> get(
            const std::string& bin
        );

        /**
         * Samples a random set of SNPs associated with a bin (with replacement).
         * 
         * Parameters:
         *  bin - String from BinsTable.bin()
         *  n_random - Number of SNPs to randomly choose from the bin.
        **/
        std::vector<std::string> get_random(
            const std::string& bin,
            size_t n_random=1
        );

        /* Translates a number of LD surrogates and an MAF into a valid key. */ 
        std::string bin(size_t surrogate_count, double maf);

    private:
        // Maximum length of serialize_surrogate_count + FIELD_SEPARATOR
        // + serialize_maf
        static const size_t MAX_KEY_SIZE = 64;
        // Special keys used to store bin cutpoints
        inline static const std::string LD_BINS_KEY = "__RESERVED_KEY_LD_BINS__";
        inline static const std::string MAF_BINS_KEY = "__RESERVED_KEY_MAF_BINS__";

        inline static const std::string BINSTABLE_EXT = "_bt";

        std::shared_ptr<VectorDiskHash> table;
        std::vector<size_t> surrogate_quantiles;
        std::vector<double> maf_quantiles;
};

/* Allows lookup of MAF and number of LD surrogates by SNP */
class SNPTable {
    public:
        /**
         * Opens an existing SNPTable.
         * 
         * Parameters:
         *  name - Name of table to open
        **/
        SNPTable(const std::string& name);

        /**
         * Creates a new SNPTable.
         * 
         * Parameters:
         *  name - Name of table to create
         *  max_key_size - Maximum length of an SNP key
        **/
        SNPTable(const std::string& name, const size_t max_key_size);

        /* Associates an SNP with an MAF and number of LD surrogates. */
        void insert(
            const std::string& snp,
            const size_t surrogate_count,
            const double maf
        );

        /* Retrieves number of LD surrogates and MAF using an SNP */
        std::pair<size_t, double> get(const std::string& snp);

    private:
        std::shared_ptr<VectorDiskHash> table;
        inline static const std::string SNPTABLE_EXT = "_snpt";
};

/* Allows lookup of markers in LD with a particular SNP */
class LDTable {
    public:
        /**
         * Opens an existing LDTable.
         * 
         * Parameters:
         *  name - Name of table to open
        **/
        LDTable(const std::string& name);

        /**
         * Creates a new LDTable.
         * 
         * Parameters:
         *  name - Name of table to create
         *  max_key_size - Maximum length of an SNP key
        **/
        LDTable(const std::string& name, const size_t max_key_size);

        /**
         * Associate an SNP with a genetic marker.
         * 
         * ld_snps associated with a particular snp must be inserted
         * consecutively. That is, do this:
         *  insert(A, B)
         *  insert(A, C)
         *  insert(D, E)
         * and not this:
         *  insert(A, B)
         *  insert(D, E)
         *  insert(A, C)
         * 
         * Parameters:
         *  snp - Index SNP ID
         *  ld_snp - SNP ID of marker in LD with index SNP
        **/
        void insert(const std::string& snp, const std::string& ld_snp);

        /* Retrieve all SNPs previously associated with an index SNP */
        std::vector<std::string> get(const std::string& snp);

    private:
        std::shared_ptr<VectorDiskHash> table;
        inline static const std::string LDTABLE_EXT = "_ldt";
};

#endif