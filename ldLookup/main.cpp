#include "CLI11.hpp"
#include "histogram.hpp"
#include "ld_parse.hpp"
#include "tables.hpp"

/*************************************************/
/***             Utility Functions             ***/
/*************************************************/

template <typename Key>
void print_key_and_values(
    const Key key,
    const std::vector<std::string>& values
) {
    for (const auto& val : values) {
        std:: cout << key << '\t' << val << '\n';
    }
}

/**
 * Extracts valid lines of linkage disequilibrium data from a file.
 * 
 * Parameters:
 *  source_path - File to read from
 *  reqs - Requirements for parsing and validating LD data.
 *  F1 - Function to be called when new LD data is parsed. F1 is
 *      passed an LDPair object representing the new data.
 *  F2 - Function to be called for each unique index SNP. F2 is passed
 *      the index SNP ID, the index SNP's MAF, and the number of SNPs
 *      in LD with the index SNP.
**/
template <typename F1, typename F2>
void iterate_genetic_data(
    const std::string& source_path,
    LDPairRequirements reqs,
    F1 on_data,
    F2 on_new_snp
) {
    std::fstream data(source_path, std::ios_base::in);
    if (!data) {
        throw std::runtime_error("Failed to open '" + source_path + "'");
    }

    std::string line;
    size_t ld_pairs(0);
    LDPair last_new_snp{ "", "", 0, 0 };
    LDPair curr{ "", "", 0 ,0 };

    while (getline(data, line)) {
        if (!parse_ld_pair(line, curr, reqs)) {
            std::cout << "Did not parse malformed line " << line << '\n';
            continue;
        } else if (!validate_ld_pair(curr, reqs)) {
            continue;
        }

        if (last_new_snp.index_snp.compare(curr.index_snp)) {
            if (ld_pairs) {
                on_new_snp(
                    last_new_snp.index_snp,
                    last_new_snp.index_maf,
                    ld_pairs
                );
            }
            last_new_snp = curr;
            ld_pairs = 0;
        }

        on_data(curr);
        ld_pairs++;
    }

    if (ld_pairs) {
        on_new_snp(
            last_new_snp.index_snp,
            last_new_snp.index_maf,
            ld_pairs
        );
    }
}

/**
 * Creates all of ldLookup's data tables. Requires that the tables
 * do not already exist.
 * 
 * Parameters:
 *  source_path - File containing linkage disequilibrium data.
 *  dir - Path to directory in which to create tables
 *  reqs - Requirements for parsing and validating LD data.
 *  n_ld_bins - Approximate number of BinsTable divisions by number
 *      of LD surrogates.
 *  n_maf_bins - Approximate number of BinsTable divisions by MAF.
**/  
void create_tables(
    const std::string& source_path,
    const std::string& dir,
    LDPairRequirements reqs,
    size_t n_ld_bins,
    size_t n_maf_bins
) {
    // First Pass: Populate LDTable and SNPTable.
    // Also, gather distribution data for MAF/LD.
    LDTable ldt(reqs.max_index_snp_size, dir);
    SNPTable snpt(reqs.max_index_snp_size, dir);

    std::map<size_t, size_t> ld_counts_hist;
    std::map<double, size_t> maf_counts_hist;

    auto first_pass_on_data([&ldt] (const LDPair& ld_pair) {
        ldt.insert(ld_pair.index_snp, ld_pair.ld_snp);
    });

    auto first_pass_on_new_snp(
        [&snpt, &ld_counts_hist, &maf_counts_hist] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        snpt.insert(index_snp_id, n_ld_pairs, maf);
        ld_counts_hist.emplace(n_ld_pairs, 0).first->second++;
        maf_counts_hist.emplace(maf, 0).first->second++;
    });

    iterate_genetic_data(
        source_path,
        reqs,
        first_pass_on_data,
        first_pass_on_new_snp
    );

    // Determine bin cutpoints.
    auto ld_quantiles(get_map_keys(
        bin_histogram(ld_counts_hist, n_ld_bins)
    ));

    auto maf_quantiles(get_map_keys(
        bin_histogram(maf_counts_hist, n_maf_bins)
    ));

    // Initialize BinsTable.
    BinsTable bt(ld_quantiles, maf_quantiles, dir);

    // Second Pass: Determine exact bin sizes.
    std::map<std::string, size_t> bin_sizes;

    auto second_pass_on_data([] (const LDPair& data) {});

    auto second_pass_on_new_snp([&bt, &bin_sizes] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        auto key(bt.bin(n_ld_pairs, maf));
        auto size_added(index_snp_id.size() + 1);
        bin_sizes.emplace(key, 0).first->second += size_added;
    });

    iterate_genetic_data(
        source_path,
        reqs,
        second_pass_on_data,
        second_pass_on_new_snp
    );

    // Reserve space in BinsTable.
    for (const auto& bin_and_size : bin_sizes) {
        bt.reserve(bin_and_size.first, bin_and_size.second);
    }

    // Third Pass: Populate BinsTable.
    auto third_pass_on_data([] (const LDPair& data) {});

    auto third_pass_on_new_snp([&bt, &bin_sizes] (
        const std::string& index_snp_id, const double maf, size_t n_ld_pairs
    ) {
        bt.insert(bt.bin(n_ld_pairs, maf), index_snp_id);
    });

    iterate_genetic_data(
        source_path,
        reqs,
        third_pass_on_data,
        third_pass_on_new_snp
    );
}

/** Extracts lines from a file with a callback for each line.
 * 
 * path - File to iterate through
 * extra_lines - Additional lines to call the callback on.
 * on_line - Callback function that takes a string parameter.
 * 
 */
template<typename F>
void iterate_lines(
    const std::string& path,
    const std::vector<std::string>& extra_lines,
    F on_line
) {
    if (path.size()) {
        std::fstream f(path, std::ios_base::in);
        if (!f) {
        throw std::runtime_error("Failed to open '" + path + "'");
    }

        std::string line;
        while (getline(f, line)) {
            on_line(line);
        }
    }

    for (auto key : extra_lines) {
        on_line(key);
    }
}


/************************************************/
/***          Command Line Interface          ***/
/************************************************/

struct SubcommandOptsCreate {
    std::string dir;
    std::string data_path;
    LDPairRequirements reqs{ ' ', 3, 7, 4, 9, 200, 0.0 };
    size_t n_ld_bins = 15;
    size_t n_maf_bins = 15;
};

struct SubcommandOptsGetLD {
    std::string dir;
    std::string path_to_markers = "";
    std::vector<std::string> cli_markers = std::vector<std::string>();
};

struct SubcommandOptsGetSimilarByValue {
    std::string dir;
    double target_maf;
    size_t target_surrogate_count;
};

struct SubcommandOptsGetSimilarBySNP {
    std::string dir;
    std::string path_to_markers = "";
    std::vector<std::string> cli_markers = std::vector<std::string>();
};

struct SubcommandOptsDistribute {
    std::string dir;
    std::string path_to_markers = "";
    std::vector<std::string> cli_markers = std::vector<std::string>();
    size_t n_distributions = 1;
};

void subcommand_create(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsCreate>());
    auto create(app.add_subcommand("create", "Create a new dataset from LD data"));

    create->add_option(
        "dir,--dir",
        opts->dir,
        "Specifies a directory that will contain the created dataset"
    )->required();

    create->add_option(
        "path_to_data,-p,--path-to-data",
        opts->data_path, 
    "File containing LD data"
    )->required();

    create->add_option(
        "-t,--threshold",
        opts->reqs.min_r_squared,
        "LD score (as measured by R^2) below which a marker will be ignored"
    );

    create->add_option(
        "-I,--index-snp-id-column",
        opts->reqs.index_snp_column,
        "Index to column of data containing index SNP IDs"
    );

    create->add_option(
        "-M,--maf-col",
        opts->reqs.index_maf_column,
        "Index to column of data containing MAFs for index SNPs"
    );

    create->add_option(
        "-L,--ld-snp-id-column",
        opts->reqs.ld_snp_column,
        "Index to column of data containing IDs for SNPs in LD with the index SNP"
    );

    create->add_option(
        "-R,--r2-index",
        opts->reqs.r_squared_column,
        "Index to column of data containing r-squared values"
    );

    create->add_option(
        "-d,--delimiter",
        opts->reqs.column_separator,
        "Character used to separate columns of data"
    );

    create->add_option(
        "-k,--key-size-limit", 
        opts->reqs.max_index_snp_size, 
        "Maximum index SNP ID length in bytes"
    );

    create->add_option(
        "-l,--ld-bin-count",
        opts->n_ld_bins,
        "Approximate number of groupings for index SNPs by number of LD surrogates"
    );

    create->add_option(
        "-m,--maf-bin-count",
        opts->n_maf_bins,
        "Approximate number of groupings for index SNPs by MAF"
    );

    create->callback([opts]() {
        // These parameters are passed as/default to one-indexed columns.
        // They need to be zero-indexed.
        opts->reqs.index_snp_column--;
        opts->reqs.ld_snp_column--;
        opts->reqs.index_maf_column--;
        opts->reqs.r_squared_column--;
        create_tables(
            opts->data_path,
            opts->dir,
            opts->reqs,
            opts->n_ld_bins,
            opts->n_maf_bins
        );
    });
}

void subcommand_get_ld(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetLD>());
    auto get_ld(app.add_subcommand(
        "get_ld", "Get markers in LD with a particular index marker"
    ));

    get_ld->add_option(
        "dir,--dir",
        opts->dir,
        "Specifies the location of the dataset to operate on (see 'ldLookup create')"
    )->required();

    get_ld->add_option(
        "path_to_markers,-p,--path-to-markers",
        opts->path_to_markers,
        "File containing newline-separated index SNP IDs"
    );

    get_ld->add_option(
        "markers,-m,--markers",
        opts->cli_markers,
        "Additional space-separated SNP IDs"
    );

    get_ld->callback([opts]() {
        LDTable ldt(opts->dir);
        auto do_get_ld([&ldt] (const std::string& snp) {
            std::vector<std::string> values = ldt.get(snp);
            print_key_and_values(snp, values);
        });
        iterate_lines(opts->path_to_markers, opts->cli_markers, do_get_ld);
    });
}

void subcommand_get_similar_by_value(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetSimilarByValue>());
    auto similar_by_value(app.add_subcommand(
        "similar_by_value",
        "Retrieve markers with MAF and number of LD surrogates near target values"
    ));

    similar_by_value->add_option(
        "dir,--dir",
        opts->dir,
        "Specifies the location of the dataset to operate on (see 'ldLookup create')"
    )->required();

    similar_by_value->add_option(
        "target_maf,-m,--target-maf",
        opts->target_maf, 
        "Target MAF value"
    )->required();

    similar_by_value->add_option(
        "target_surrogates,-s,--target-surrogates",
        opts->target_surrogate_count,
        "Target number of LD surrogates"
    )->required();

    similar_by_value->callback([opts]() {
        BinsTable bt(opts->dir);
        SNPTable snpt(opts->dir);
        auto values(bt.get(
            bt.bin(opts->target_surrogate_count, opts->target_maf)
        ));

        std::string key(serialize_maf(opts->target_maf) 
                        + " " 
                        + serialize_surrogates(opts->target_surrogate_count));
        print_key_and_values(key, values);
    });
}

void subcommand_get_similar_by_snp(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsGetSimilarBySNP>());
    auto similar_by_snp(app.add_subcommand(
        "similar_by_snp",
        "Retrieve markers with MAF and number of LD surrogates similar to a key marker"
    ));

    similar_by_snp->add_option(
        "dir,--dir",
        opts->dir,
        "Specifies the location of the dataset to operate on (see 'ldLookup create')"
    )->required();

    similar_by_snp->add_option(
        "path_to_markers,-p,--path-to-markers",
        opts->path_to_markers,
        "File containing newline-separated index SNP IDs"
    );

    similar_by_snp->add_option(
        "markers,-m,--markers",
        opts->cli_markers,
        "Additional space-separated SNP IDs"
    );

    similar_by_snp->callback([opts]() {
        BinsTable bt(opts->dir);
        SNPTable snpt(opts->dir);

        auto do_get_similar([&bt, &snpt] (const std::string& snp) {
            auto maf_and_surrogates(snpt.get(snp));
            auto values(bt.get(bt.bin(
                maf_and_surrogates.first,
                maf_and_surrogates.second
            )));
            print_key_and_values(snp, values);
        });

        iterate_lines(
            opts->path_to_markers,
            opts->cli_markers,
            do_get_similar
        );
    });
}

void subcommand_distribute(CLI::App& app) {
    auto opts(std::make_shared<SubcommandOptsDistribute>());
    auto distribute(app.add_subcommand(
        "distribute",
        "Generate distributions for hypothesis testing"
    ));

    distribute->add_option(
        "dir,--dir",
        opts->dir,
        "Specifies the location of the dataset to operate on (see 'ldLookup create')"
    )->required();

    distribute->add_option(
        "path_to_markers,-p,--path-to-markers",
        opts->path_to_markers,
        "File containing newline-separated index SNP IDs"
    );

    distribute->add_option(
        "markers,-m,--markers",
        opts->cli_markers,
        "Additional space-separated SNP IDs"
    );

    distribute->add_option(
        "-d,--distributions",
        opts->n_distributions,
        "Number of distributions to generate"
    );

    distribute->callback([opts]() {
        BinsTable bt(opts->dir);
        SNPTable snpt(opts->dir);
        std::map<size_t, std::vector<std::string>> dists;

        for (size_t i = 0; i < opts->n_distributions; i++) {
            dists.emplace(i, std::vector<std::string>());
        }

        auto do_distribute([&] (const std::string& snp) {
            auto maf_and_surrogates = snpt.get(snp);
            auto values(bt.get_random(
                bt.bin(
                    maf_and_surrogates.first,
                    maf_and_surrogates.second
                ),
                opts->n_distributions
            ));

            for (size_t i = 0; i < values.size(); i++) {
                dists[i].push_back(values.at(i));
            }
        });

        iterate_lines(
            opts->path_to_markers,
            opts->cli_markers,
            do_distribute
        );

        for (auto const& number_and_dist : dists) {
            print_key_and_values(
                "Distribution #" + std::to_string(number_and_dist.first+1),
                number_and_dist.second
            );
        }
    });
}

int main(int argc, char** argv) {
    std::string description("ldLookup - lookup and analysis of linkage disequilibrium (LD) between genetic variants");
    CLI::App app{description};
    app.option_defaults()->always_capture_default();
    app.require_subcommand(1);

    subcommand_create(app);
    subcommand_get_ld(app);
    subcommand_get_similar_by_value(app);
    subcommand_get_similar_by_snp(app);
    subcommand_distribute(app);

    // Application Logic
    try {
        CLI11_PARSE(app, argc, argv);
    } catch (std::exception &e) {
        std::cout << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
