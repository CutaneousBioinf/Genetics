#ifndef LDLOOKUP_TABLES_HPP
#define LDLOOKUP_TABLES_HPP

#include <map> // std::map
#include <memory>
#include <string> // std::string
#include <vector> // std::vector

#include "geneticData.hpp"
#include "vdh.hpp"

const char FIELD_SEPARATOR = ' ';

std::string serialize_ld(size_t n_ld_pairs);

std::string serialize_maf(double maf);

size_t deserialize_ld(const std::string& n_ld_pairs);

double deserialize_maf(const std::string& maf);

/* Allows fuzzy lookup of SNPs by MAF and number of markers in LD */
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
         *  ld_quantiles - Lower bounds for LD bins.
         *  maf_quantiles - Lower bounds for MAF bins.
        **/
        BinsTable(
            const std::string& name,
            const std::vector<size_t>& ld_quantiles,
            const std::vector<double>& maf_quantiles
        );

        void reserve(
            const std::string& bin,
            size_t space
        );

        void insert(
            const std::string& bin,
            const std::string& snp
        );

        std::vector<std::string> get(
            const std::string& bin
        );

        std::vector<std::string> get_random(
            const std::string& bin,
            size_t n_random=1
        );

        std::string bin(size_t n_ld_pairs, double maf);

    private:
        static const size_t MAX_KEY_SIZE = 64;
        inline static const std::string LD_BINS_KEY = "__RESERVED_KEY_LD_BINS__";
        inline static const std::string MAF_BINS_KEY = "__RESERVED_KEY_MAF_BINS__";

        std::shared_ptr<VectorDiskHash> table;
        std::vector<size_t> ld_quantiles;
        std::vector<double> maf_quantiles;
};

class SNPTable {
    public:
        SNPTable(const std::string& name);

        SNPTable(const std::string& name, const size_t max_key_size);

        void insert(
            const std::string& snp,
            const size_t n_ld_pairs,
            const double maf
        );

        std::pair<size_t, double> get(const std::string& snp);

    private:
        std::shared_ptr<VectorDiskHash> table;
};

class LDTable {
    public:
        LDTable(const std::string& name);

        LDTable(const std::string& name, const size_t max_key_size);

        void insert(const std::string& snp, const std::string& ld_surrogate);

        std::vector<std::string> get(const std::string& snp);

    private:
        std::shared_ptr<VectorDiskHash> table;
};

#endif