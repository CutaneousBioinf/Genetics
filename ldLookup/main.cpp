#include <algorithm>   // std::find
#include <filesystem>  // std::filesystem::create_directory
#include <iostream>    // std::cout
#include <fstream>     // std::ifstream
#include <memory>      // std::shared_ptr
#include <utility>     // std::pair
#include <string>
#include <vector>

#include <stddef.h>    // size_t

#include "CLI11.hpp"
#include "parse_variants.hpp"
#include "stratify.hpp"
#include "tables.hpp"

using std::string;
using std::vector;

/*************************************************/
/*************************************************/
/***                 Constants                 ***/
/*************************************************/
/*************************************************/

const string LD_TABLE_FILE_PATH = "ld.vdhdat";
const string LD_TABLE_TABLE_PATH = "ld.vdhdht";
const string STRATA_TABLE_FILE_PATH = "strata.vdhdat";
const string STRATA_TABLE_TABLE_PATH = "strata.vdhdht";
const string SUMMARY_TABLE_FILE_PATH = "summary.vdhdat";
const string SUMMARY_TABLE_TABLE_PATH = "summary.vdhdht";

/*************************************************/
/*************************************************/
/***                  Options                  ***/
/*************************************************/
/*************************************************/

struct SubcommandOptsSetup {
    string dir;
    string src;
    char delimiter = ' ';
    std::string index_variant_id_column = "SNP_A";
    std::string ld_variant_id_column = "SNP_B";
    std::string index_variant_maf_column = "MAF_A";
    std::string r2_column = "R2";
    double r2_threshold_for_ld = 0.0;
    size_t index_variants_per_ld_bin = 0;
    size_t n_ld_bins = 0;
    size_t index_variants_per_maf_bin = 0;
    size_t n_maf_bins = 0;
};

struct SubcommandOptsGetVariantsInLDWith {
    string dir;
    string key_variants_file = "";
    vector<string> key_variants = vector<string>();
};

struct SubcommandOptsGetVariantsSimilarTo {
    string dir;
    string key_variants_file = "";
    vector<string> key_variants = vector<string>();
};

struct SubcommandOptsGetVariantsWithStatsLike {
    string dir;
    double target_maf;
    size_t target_surrogate_count;
};

struct SubcommandOptsGetVariantStatistics {
    string dir;
    string key_variants_file = "";
    vector<string> key_variants = vector<string>();
};

struct SubcommandOptsSample {
    string dir;
    string key_variants_file = "";
    vector<string> key_variants = vector<string>();
    size_t n_samples = 1;
};

/*************************************************/
/*************************************************/
/***                   Types                   ***/
/*************************************************/
/*************************************************/

struct SetupFirstIterationResults {
    Histogram<size_t> n_surrogates_hist;
    Histogram<double> maf_hist;
    size_t max_index_variant_size;
};

struct Tables {
    std::shared_ptr<LDTable> ld_t;
    std::shared_ptr<StrataTable> strata_t;
    std::shared_ptr<SummaryTable> summary_t;
};

/*************************************************/
/*************************************************/
/***                  Helpers                  ***/
/*************************************************/
/*************************************************/

void on_invalid_cb(const string& line) {
    std::cout << "Ignoring Invalid Line: " << line << std::endl;
}

Tables open_tables(const std::string& dir_str) {
    std::filesystem::path dir(dir_str);

	Tables ret;
	ret.ld_t.reset(new LDTable(
	    {dir / LD_TABLE_FILE_PATH,
	     dir / LD_TABLE_TABLE_PATH,
	     0,
	     false}));
	ret.strata_t.reset(new StrataTable(
	    dir / STRATA_TABLE_FILE_PATH,
	    dir / STRATA_TABLE_TABLE_PATH));
	ret.summary_t.reset(new SummaryTable(
	    {dir / SUMMARY_TABLE_FILE_PATH,
	     dir / SUMMARY_TABLE_TABLE_PATH,
	     0,
	     false}));

	return ret;
}

void print_to_columns(
    const string& first_col,
    const vector<string>& second_col) {
    for (const string& s : second_col) {
        std::cout << first_col << '\t' << s << '\n';
    }
}

string read_first_line(string path) {
    // Open the source input file.
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open '" + path + "'");
    }

    // Read the header row of file.
    string header;
    std::getline(file, header);

    return header;
}

/*************************************************/
/*************************************************/
/***                   Logic                   ***/
/*************************************************/
/*************************************************/

LDPairParser create_parser(std::shared_ptr<SubcommandOptsSetup> opts) {
    // Read column names from the header row.
    string header = read_first_line(opts->src);
    vector<string> col_names = split(header, opts->delimiter);

    // This is the parser we will return.
    LDPairParser parser;
    parser.delimiter = opts->delimiter;
    parser.r2_threshold_for_ld = opts->r2_threshold_for_ld;

    // Convert column strings to column indices.
    vector<string> str_args = {
        opts->index_variant_id_column,
        opts->index_variant_maf_column,
        opts->ld_variant_id_column,
        opts->r2_column
    };

    vector<size_t*> col_idx_fields = {
        &parser.index_variant_id_column,
        &parser.index_variant_maf_column,
        &parser.ld_variant_id_column,
        &parser.r2_column
    };

    for (size_t i = 0; i < str_args.size(); i++) {
        string arg = str_args[i];
        try {
            // Treat the user argument as a number representing
            // a column index.
            *col_idx_fields[i] = std::stoi(arg);
        } catch (std::invalid_argument& ignore) {
            // On failure, treat the user argument as the name
            // of a column, and determine its index.
            auto it = std::find(col_names.begin(), col_names.end(), arg);
            if (it == col_names.end()) {
                throw std::invalid_argument("Invalid Column: " + arg);
            }
            *col_idx_fields[i] = it - col_names.begin() + 1;
        }
    }

    return parser;  
}

SetupFirstIterationResults do_setup_first_iteration(
    const LDPairParser& parser,
    std::shared_ptr<SubcommandOptsSetup> opts
) {
    SetupFirstIterationResults results;
    results.max_index_variant_size = 0;

	auto it1_on_ld_pair_cb = [](const LDPair& pair) {
        (void)pair;
    };

	auto it1_on_new_index_variant_cb = [&](const IndexVariantSummary& summary) {
        // Update maximum variant size.
        if (summary.variant_id.size() > results.max_index_variant_size) {
            results.max_index_variant_size = summary.variant_id.size();
        }

        // Update distribution information for StrataTable.
        results.maf_hist.increase_count(summary.maf);
        results.n_surrogates_hist.increase_count(summary.n_surrogates);
    };

	iterate_ld_data(
	    opts->src,
	    parser,
	    it1_on_ld_pair_cb,
	    it1_on_new_index_variant_cb,
	    on_invalid_cb);

    // Determine strata for MAF.
    size_t n_maf_bins = opts->n_maf_bins;
    if (opts->index_variants_per_maf_bin != 0) {
        auto total = results.maf_hist.total_count();
		n_maf_bins = total / opts->index_variants_per_maf_bin;
	}
    results.maf_hist = results.maf_hist.stratify(n_maf_bins);

    // Determine strata for number of LD surrogates.
    size_t n_ld_bins = opts->n_ld_bins;
    if (opts->index_variants_per_ld_bin != 0) {
        auto total = results.n_surrogates_hist.total_count();
        n_ld_bins = total / opts->index_variants_per_ld_bin;
    }
    results.n_surrogates_hist = results.n_surrogates_hist.stratify(n_ld_bins);

    return results;
}

void do_setup(std::shared_ptr<SubcommandOptsSetup> opts) {
    // Create an string-to-LDPair parser based on opts.
    LDPairParser parser = create_parser(opts);

    // Do first iteration.
    SetupFirstIterationResults results = do_setup_first_iteration(parser, opts);

    // Create the output directory.
    std::filesystem::path dir(opts->dir);
    if (std::filesystem::exists(dir)) {
        throw std::runtime_error("Directory Already Exists: " + opts->dir);
    }
    std::filesystem::create_directory(dir);

    // Initialize LDTable, StrataTable, and SummaryTable.
    Options ld_table_options{
        dir / LD_TABLE_FILE_PATH,
        dir / LD_TABLE_TABLE_PATH,
        results.max_index_variant_size,
        true };
    
    Options summary_table_options{
        dir / SUMMARY_TABLE_FILE_PATH,
        dir / SUMMARY_TABLE_TABLE_PATH,
        results.max_index_variant_size,
        true };
    
    LDTable ld_t(ld_table_options);
    StrataTable strata_t(
        dir / STRATA_TABLE_FILE_PATH,
        dir / STRATA_TABLE_TABLE_PATH,
        results.n_surrogates_hist,
        results.maf_hist);
    SummaryTable summary_t(summary_table_options);

    // Second Iteration Over Data:
    // - Determine space needed for strata.
    Histogram<StrataTable::Stratum> strata_sizes;
	auto it2_on_ld_pair_cb = [](const LDPair& pair) {
        (void)pair;
    };

	auto it2_on_new_index_variant_cb = [&](const IndexVariantSummary& summary) {
        StrataTable::Stratum stratum = strata_t.get_stratum(summary);
        size_t variant_size = summary.variant_id.size()+1;
        strata_sizes.increase_count(stratum, variant_size);
    };

    iterate_ld_data(
	    opts->src,
	    parser,
	    it2_on_ld_pair_cb,
	    it2_on_new_index_variant_cb,
	    on_invalid_cb);
    
    // Reserve strata on-disk.
    strata_t.reserve(strata_sizes);

    // Third Iteration Over Data:
    // - Populate LDTable, StatsTable, and StrataTable.
    auto it3_on_ld_pair_cb = [&](const LDPair& pair) {
        ld_t.append(pair.index_variant_id, pair.ld_variant_id);
	};

	auto it3_on_new_index_variant_cb = [&](const IndexVariantSummary& summary) {
        strata_t.append(summary);
        summary_t.append(summary);
    };

    iterate_ld_data(
	    opts->src,
	    parser,
	    it3_on_ld_pair_cb,
	    it3_on_new_index_variant_cb,
	    on_invalid_cb);
}

void do_get_variants_in_ld_with(
    std::shared_ptr<SubcommandOptsGetVariantsInLDWith> opts) {
    
    Tables tables = open_tables(opts->dir);
    std::cout << "Variant ID\tVariant ID of LD Surrogate\n";
    auto on_variant = [&](string variant) {
        print_to_columns(variant, tables.ld_t->lookup(variant));
    };

    iterate_variants(
        opts->key_variants_file,
        opts->key_variants,
        on_variant);
}

void do_get_variants_similar_to(
    std::shared_ptr<SubcommandOptsGetVariantsSimilarTo> opts) {
    
    Tables tables = open_tables(opts->dir);
    std::cout << "Variant ID\tVariant ID of Similar Variant\n";
    auto on_variant = [&](string variant) {
        IndexVariantSummary stats = tables.summary_t->lookup(variant);
        print_to_columns(variant, tables.strata_t->lookup(stats));
    };

    iterate_variants(
        opts->key_variants_file,
        opts->key_variants,
        on_variant);
}

void do_get_variants_with_stats_like(
    std::shared_ptr<SubcommandOptsGetVariantsWithStatsLike> opts) {
    
    Tables tables = open_tables(opts->dir);
    IndexVariantSummary stats;
    stats.n_surrogates = opts->target_surrogate_count;
    stats.maf = opts->target_maf;
    
    std::cout << "Target MAF\tTarget # LD Surrogates\tVariant ID\n";
    for (auto &s : tables.strata_t->lookup(stats)) {
        std::cout << opts->target_maf << '\t';
        std::cout << opts->target_surrogate_count << '\t';
        std::cout << s << '\n';
    }
}

void do_get_variant_statistics(
    std::shared_ptr<SubcommandOptsGetVariantStatistics> opts) {
    
    Tables tables = open_tables(opts->dir);
    std::cout << "Variant ID\t# LD Surrogates\tMAF\n";
    auto on_variant = [&](string variant) {
        IndexVariantSummary stats = tables.summary_t->lookup(variant);
        std::cout << variant << '\t' << stats.n_surrogates << '\t';
        std::cout << stats.maf << '\n';
    };

    iterate_variants(
        opts->key_variants_file,
        opts->key_variants,
        on_variant);
}

void do_sample(std::shared_ptr<SubcommandOptsSample> opts) {
    Tables tables = open_tables(opts->dir);
    std::cout << "Sample #\tVariant ID\tVariant ID of Similar Variant\n";
    auto on_variant = [&](string variant) {
        IndexVariantSummary stats = tables.summary_t->lookup(variant);
        auto sampled = tables.strata_t->lookup_sample(stats, opts->n_samples);
        for (size_t i = 0; i < sampled.size(); i++) {
            std::cout << i+1 << '\t' << variant << '\t' << sampled.at(i) << '\n';
        }
    };

    iterate_variants(
        opts->key_variants_file,
        opts->key_variants,
        on_variant);
}

/*************************************************/
/*************************************************/
/***                    CLI                    ***/
/*************************************************/
/*************************************************/

void subcommand_setup(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsSetup>());
    auto cmd(app.add_subcommand("setup", "Create a new lookup table"));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory in which to store the lookup table"
    )->check(CLI::NonexistentPath)->required();

    cmd->add_option(
        "src,-s,--src",
        opts->src,
        "File from which to read LD data"
    )->check(CLI::ExistingFile)->required();

    cmd->add_option(
        "-d,--delimiter",
        opts->delimiter,
        "Character that separates columns of LD data"
    );

    cmd->add_option(
        "-I,--index-id-column",
        opts->index_variant_id_column,
        "Column of LD data containing index variant IDs"
    );

    cmd->add_option(
        "-L,--ld-id-column",
        opts->ld_variant_id_column,
        "Column of LD data containing IDs for variants in LD with the index variant"
    );

    cmd->add_option(
        "-M,--index-maf-column",
        opts->index_variant_maf_column,
        "Column of LD data containing MAFs of index variants"
    );

    cmd->add_option(
        "-R,--r2-column",
        opts->r2_column,
        "Column of LD data containing r-squared values"
    );

    cmd->add_option(
        "-t,--r2-threshold-for-ld",
        opts->r2_threshold_for_ld,
        "Minimum r-squared value for a variant pair to be considered 'in LD'"
    )->check(CLI::Range(0.0, 1.0));

    // LD Bin Group
    auto ld_group = cmd->add_option_group("ld_bins");
    ld_group->add_option(
        "--index-variants-per-ld-bin",
        opts->index_variants_per_ld_bin,
        "Approximate size of each strata when index variants are stratified by number of LD surrogates"
    )->check(CLI::PositiveNumber);

    ld_group->add_option(
        "--n-ld-bins",
        opts->n_ld_bins,
        "Approximate number of strata when index variants are stratified by number of LD surrogates"
    )->check(CLI::PositiveNumber);

    ld_group->require_option(1);
    // End LD Bin Group

    // MAF Bin Group
    auto maf_group = cmd->add_option_group("maf_bins");
    maf_group->add_option(
        "--index-variants-per-maf-bin",
        opts->index_variants_per_maf_bin,
        "Approximate size of each strata when index variants are stratified by MAF"
    )->check(CLI::PositiveNumber);

    maf_group->add_option(
        "--n-maf-bins",
        opts->n_maf_bins,
        "Approximate number of strata when index variants are stratified by MAF"
    )->check(CLI::PositiveNumber);

    maf_group->require_option(1);
    // End MAF Bin Group

    cmd->callback([opts]() {
        do_setup(opts);
    });
}

void subcommand_get_variants_in_ld_with(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetVariantsInLDWith>());
    auto cmd(app.add_subcommand(
        "get_variants_in_ld_with", 
        "Get variants in LD with specified key variants"
    ));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory where lookup table is stored"
    )->check(CLI::ExistingDirectory)->required();

    cmd->add_option(
        "key_variants_file,-f,--key-variants-file",
        opts->key_variants_file,
        "File containing newline-separated index variant IDs"
    )->check(CLI::ExistingFile);

    cmd->add_option(
        "key_variants,-k,--key-variants",
        opts->key_variants,
        "Space-separated index variant IDs"
    );

    cmd->callback([opts]() {
        do_get_variants_in_ld_with(opts);
    });
}

void subcommand_get_variants_similar_to(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetVariantsSimilarTo>());
    auto cmd(app.add_subcommand(
        "get_variants_similar_to",
        "Get variants with MAF and number of LD surrogates similar to those of specified key variants"
    ));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory where lookup table is stored"
    )->check(CLI::ExistingDirectory)->required();

    cmd->add_option(
        "key_variants_file,-f,--key-variants-file",
        opts->key_variants_file,
        "File containing newline-separated index variant IDs"
    )->check(CLI::ExistingFile);

    cmd->add_option(
        "key_variants,-k,--key-variants",
        opts->key_variants,
        "Space-separated index variant IDs"
    );

    cmd->callback([opts]() {
        do_get_variants_similar_to(opts);
    });
}

void subcommand_get_variants_with_stats_like(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetVariantsWithStatsLike>());
    auto cmd(app.add_subcommand(
        "get_variants_with_stats_like",
        "Get variants with MAF and number of LD surrogates near specified targets"
    ));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory where lookup table is stored"
    )->check(CLI::ExistingDirectory)->required();

    cmd->add_option(
        "target_maf,-m,--target-maf",
        opts->target_maf, 
        "Target MAF value"
    )->required();

    cmd->add_option(
        "target_n_ld_surrogates,-n,--target-n-ld-surrogates",
        opts->target_surrogate_count,
        "Target number of LD surrogates"
    )->required();

    cmd->callback([opts]() {
        do_get_variants_with_stats_like(opts);
    });
}

void subcommand_get_variant_statistics(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetVariantStatistics>());
    auto cmd(app.add_subcommand(
        "get_variant_statistics",
        "Get MAF and number of LD surrogates of specified key variants"
    ));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory where lookup table is stored"
    );

    cmd->add_option(
        "key_variants_file,-f,--key-variants-file",
        opts->key_variants_file,
        "File containing newline-separated index variant IDs"
    )->check(CLI::ExistingFile);

    cmd->add_option(
        "key_variants,-k,--key-variants",
        opts->key_variants,
        "Space-separated index variant IDs"
    );

    cmd->callback([opts]() {
        do_get_variant_statistics(opts);
    });
}

void subcommand_sample(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsSample>());
    auto cmd(app.add_subcommand(
        "sample",
        "Randomly sample variants with MAF and number of LD surrogates similar to those of specified key variants"
    ));

    cmd->add_option(
        "dir,--dir",
        opts->dir,
        "Directory where lookup table is stored"
    )->check(CLI::ExistingDirectory)->required();

    cmd->add_option(
        "key_variants_file,-f,--key-variants-file",
        opts->key_variants_file,
        "File containing newline-separated index variant IDs"
    )->check(CLI::ExistingFile);

    cmd->add_option(
        "key_variants,-k,--key-variants",
        opts->key_variants,
        "Space-separated index variant IDs"
    );

    cmd->add_option(
        "-n,--n-samples",
        opts->n_samples,
        "Number of samples to take for each variant"
    );

    cmd->callback([opts]() {
        do_sample(opts);
    });
}

/************************************************/
/************************************************/
/***                   MAIN                   ***/
/************************************************/
/************************************************/

int main(int argc, char** argv) {
    CLI::App app("ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants");
    app.option_defaults()->always_capture_default();
    app.require_subcommand(1);

    subcommand_setup(app);
    subcommand_get_variants_in_ld_with(app);
    subcommand_get_variants_similar_to(app);
    subcommand_get_variants_with_stats_like(app);
    subcommand_get_variant_statistics(app);
    subcommand_sample(app);

    // Application Logic
    try {
        CLI11_PARSE(app, argc, argv);
    } catch (std::exception &e) {
        std::cout << "ERROR:\n" << e.what() << "\n";
        return 1;
    }

    return 0;
}