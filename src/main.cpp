#include <cli.hpp>
#include <fmindex.hpp>
#include <input.hpp>
#include <output.hpp>
#include <pex.hpp>
#include <search.hpp>
#include <verification.hpp>

#include <atomic>
#include <exception>
#include <filesystem>
#include <ranges>
#include <span>
#include <vector>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <ivsigma/ivsigma.h>

int main(int argc, char** argv) {
    fmt::println("/\\ \\/ /\\ \\/ /\\ welcome to floxer /\\ \\/ /\\ \\/ /\\");
    
    cli::options opt;
    try {
        opt = cli::parse_and_validate_options(argc, argv);
    } catch (std::exception const & e) {
        fmt::print(stderr, "[CLI PARSER ERROR]\n{}\n", e.what());
        return -1;
    }
    
    fmt::println("--> reading reference sequences from {} ... ", opt.reference_sequence.c_str());

    std::vector<input::reference_record> references;
    try {
        references = input::read_references(opt.reference_sequence);
    } catch (std::exception const& e) {
        fmt::print(
            stderr,
            "[INPUT ERROR]\nAn error occured while trying to read the reference from "
            "the file {}.\n{}\n",
            opt.reference_sequence.c_str(),
            e.what()
        );
        return -1;
    }

    fmt::println("     ... done.");

    fmindex index;
    if (!opt.index_path.empty() && std::filesystem::exists(opt.index_path)) {
        try {
            fmt::println(" --> loading index from {} ... ", opt.index_path.c_str());
            
            index = input::load_index(opt.index_path);

            fmt::println("      ... done.");
        } catch (std::exception const& e) {
            fmt::print(
                stderr,
                "[INPUT ERROR]\nAn error occured while trying to load the index from "
                "the file {}.\n{}\n",
                opt.index_path.c_str(),
                e.what()
            );
            return -1;
        }
    } else {
        fmt::println(
            " --> building index with {} thread{} ... ",
            opt.num_threads,
            opt.num_threads == 1 ? "" : "s"
        );

        size_t constexpr suffix_array_sampling_rate = 16; 
        index = fmindex(
            references | std::views::transform(&input::reference_record::rank_sequence),
            suffix_array_sampling_rate,
            opt.num_threads
        );

        fmt::println("      ... done.");

        if (!opt.index_path.empty()) {
            try {
                fmt::println(" --> saving index to {} ... ", opt.index_path.c_str());

                output::save_index(index, opt.index_path);

                fmt::println("     ... done.");
            } catch (std::exception const& e) {
                fmt::print(
                    stderr,
                    "[OUTPUT WARNING]\nAn error occured while trying to write the index to "
                    "the file {}.\nContinueing without saving the index.\n{}\n",
                    opt.index_path.c_str(),
                    e.what()
                );
            }
        }
    }

    fmt::println("  --> reading queries from {} ... ", opt.queries.c_str());

    std::vector<input::query_record> fastq_queries;
    try {
        fastq_queries = input::read_queries(opt.queries);
    } catch (std::exception const& e) {
        fmt::print(
            stderr,
            "[INPUT ERROR]\nAn error occured while trying to read the queries from "
            "the file {}.\n{}\n",
            opt.queries.c_str(),
            e.what()
        );
        return -1;
    }

    fmt::println("       ... done.");

    auto sam_output = output::sam_output(opt.output_path, references);

    search::search_scheme_cache scheme_cache;
    pex_tree_cache tree_cache;

    fmt::println(
        "   --> aligning queries with {} thread{} and writing output file to {} ... ",
        opt.num_threads,
        opt.num_threads == 1 ? "" : "s",
        opt.output_path.c_str()
    );

    // some workarounds for handling errors in threads
    std::atomic_bool encountered_error = false;
    std::vector<std::exception_ptr> exceptions{};

    #pragma omp parallel for \
        num_threads(opt.num_threads) \
        default(none) \
        private(tree_cache, scheme_cache) \
        shared(fastq_queries, opt, references, index, sam_output, exceptions, encountered_error) \
        schedule(static)
    for (size_t i = 0; i < fastq_queries.size(); ++i) {
        if (encountered_error) {
            continue;
        }

        try {
            auto const& fastq_query = fastq_queries[i];
            size_t const query_num_errors = fastq_query.num_errors_from_user_config(opt);

            auto const tree_config = pex_tree_config {
                .total_query_length = fastq_query.sequence_length,
                .query_num_errors = query_num_errors,
                .leaf_max_num_errors = opt.pex_leaf_num_errors
            };
            auto const& tree = tree_cache.get(tree_config);

            auto alignments = verification::fastq_query_alignments(references.size());

            bool is_reverse_complement = false;
            tree.search(
                references,
                fastq_query.rank_sequence,
                alignments,
                is_reverse_complement,
                scheme_cache,
                index
            );

            auto const reverse_complement_fastq_query_rank_sequence = 
                ivs::reverse_complement_rank<ivs::d_dna4>(fastq_query.rank_sequence);
            is_reverse_complement = true;

            tree.search(
                references,
                reverse_complement_fastq_query_rank_sequence,
                alignments,
                is_reverse_complement,
                scheme_cache,
                index
            );

            #pragma omp critical
            sam_output.output_for_query(
                fastq_query,
                references,
                alignments
            );
        } catch (...) {
            #pragma omp critical
            exceptions.emplace_back(std::current_exception());

            encountered_error = true;
        }
    }

    for (auto const& e : exceptions) {
        try {
            std::rethrow_exception(e);
        } catch (std::exception const& e) {
            fmt::print(
                stderr,
                "[ERROR]\nAn error occured while a thread was aligning reads or writing output to "
                "the file {}.\nThe output file is likely incomplete and invalid.\n{}\n",
                opt.output_path.c_str(),
                e.what()
            );
        }
    } 

    if (!exceptions.empty()) {
        return -1;
    }

    fmt::println("        ... done.");

    return 0;
}
