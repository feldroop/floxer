#include <about_floxer.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <filesystem>
#include <iterator>
#include <optional>
#include <random>
#include <ranges>
#include <unordered_map>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include <sharg/all.hpp>

#include <spdlog/fmt/ranges.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

struct chromosome_t {
    std::string name;
    seqan3::dna4_vector sequence;
};

using genome_t = std::vector<chromosome_t>;

genome_t create_genome(
    size_t const chromosome_length,
    size_t const num_chromosomes,
    std::mt19937& random_generator
) {
    genome_t genome(num_chromosomes);
    std::uniform_int_distribution<uint8_t> base_rank_distribution(0, 3);

    for (size_t i = 0; i < num_chromosomes; ++i) {
        genome[i].name = fmt::format("chromosome_{}", i);
        auto& chromosome_sequence = genome[i].sequence;
        chromosome_sequence.resize(chromosome_length);

        for (auto& base : chromosome_sequence) {
            base.assign_rank(base_rank_distribution(random_generator));
        }
    }

    return genome;
}

void write_genome_to_file(genome_t const& genome, std::filesystem::path const& genome_path) {
    seqan3::sequence_file_output genome_output{genome_path};

    for (auto const& chromosome : genome) {
        genome_output.emplace_back(chromosome.sequence, chromosome.name);
    }
}

enum class mutation_kind_t {
    mismatch, insertion, deletion
};

static const std::array<mutation_kind_t, 3> possible_mutation_kinds = {
    mutation_kind_t::mismatch,
    mutation_kind_t::insertion,
    mutation_kind_t::deletion
};

struct mutation_t {
    size_t index;
    mutation_kind_t kind;
    seqan3::dna4 new_base;
};

uint8_t choose_distinct_rank(uint8_t const generated_rank, uint8_t const origin_rank) {
    return generated_rank >= origin_rank ?
        generated_rank + 1 :
        generated_rank;
}

void create_and_write_reads(
    size_t const base_read_length,
    size_t const num_reads,
    double const error_rate,
    genome_t const& genome,
    size_t const chromosome_length,
    std::mt19937& random_generator,
    std::filesystem::path const& read_path
) {
    size_t const num_errors = error_rate * base_read_length;

    seqan3::dna4_vector read_sequence{};

    std::uniform_int_distribution<size_t> chromosome_index_distribution(0, genome.size() - 1);
    std::uniform_int_distribution<size_t> read_start_index_distribution(0, chromosome_length - base_read_length - 1);
    std::uniform_int_distribution<uint8_t> index_up_to_2_distribution(0, 2);
    std::uniform_int_distribution<uint8_t> index_up_to_3_distribution(0, 3);

    std::vector<size_t> read_indices{};
    std::vector<size_t> mutation_indices{};
    std::vector<mutation_t> mutations{};

    using namespace seqan3::literals;
    std::vector<seqan3::phred42> read_quality(base_read_length + num_errors, 'I'_phred42);

    seqan3::sequence_file_output read_output{read_path};

    for (size_t read_id = 0; read_id < num_reads; ++read_id) {
        size_t const chromosome_index = chromosome_index_distribution(random_generator);
        auto const& chromosome_sequence = genome[chromosome_index].sequence;

        size_t const read_start_index = read_start_index_distribution(random_generator);

        auto read_indices_v = std::ranges::iota_view{read_start_index, read_start_index + base_read_length};
        // this sadly is necessary, because iota_view not not seem to
        // be a forward random access iterator or something like that
        std::copy(
            read_indices_v.begin(),
            read_indices_v.end(),
            std::back_inserter(read_indices)
        );

        // the mutations work as follows: exactly num_errors distinct indices are chosen for a mutation
        // for deletions, the origin base of that index is deleted
        // for mismatches, the origin base of that index is definitely changed (never stays the same)
        // for insertions, the origin base of that index is kept and a new random base is inserted after it
        // known limitation: neighboring insertions can undo deletions and vice versa
        std::sample(
            read_indices.begin(),
            read_indices.end(),
            std::back_inserter(mutation_indices),
            num_errors,
            random_generator
        );

        for (auto const mutation_index : mutation_indices) {
            auto const mutation_kind = possible_mutation_kinds[index_up_to_2_distribution(random_generator)];

            seqan3::dna4 new_base{};
            auto const origin_rank = chromosome_sequence[mutation_index].to_rank();
            switch (mutation_kind) {
                case mutation_kind_t::mismatch:
                    new_base.assign_rank(
                        choose_distinct_rank(index_up_to_2_distribution(random_generator), origin_rank)
                    );
                    break;

                case mutation_kind_t::insertion:
                    new_base.assign_rank(index_up_to_3_distribution(random_generator));
                    break;

                case mutation_kind_t::deletion:
                default:
                    break;
            }

            mutations.push_back(mutation_t {
                .index = mutation_index,
                .kind = mutation_kind,
                .new_base = new_base
            });
        }

        assert(mutations.size() == num_errors);

        std::ranges::sort(mutations, {}, [] (mutation_t const& mutation) { return mutation.index; });

        size_t mutation_index_in_mutations = 0;
        size_t origin_index_offset = 0;

        while (origin_index_offset < base_read_length) {
            size_t const curr_origin_index = read_start_index + origin_index_offset;
            auto const origin_base = chromosome_sequence[curr_origin_index];

            if (mutation_index_in_mutations == mutations.size()) {
                read_sequence.push_back(origin_base);
                ++origin_index_offset;
                continue;
            }

            auto const& mutation = mutations[mutation_index_in_mutations];
            if (mutation.index != curr_origin_index) {
                read_sequence.push_back(origin_base);
                ++origin_index_offset;
                continue;
            }

            switch (mutation.kind) {
                case mutation_kind_t::mismatch:
                    read_sequence.push_back(mutation.new_base);
                    break;

                case mutation_kind_t::insertion:
                    read_sequence.push_back(origin_base);
                    read_sequence.push_back(mutation.new_base);
                    break;

                case mutation_kind_t::deletion:
                default:
                    break;
            }

            ++origin_index_offset;
            ++mutation_index_in_mutations;
        }

        std::string read_name = fmt::format(
            "id_{}_chromosome_{}_position_{}_max_errors_{}",
            read_id,
            chromosome_index,
            read_start_index,
            num_errors
        );

        auto const quality_span = std::span(read_quality).subspan(0, read_sequence.size());
        read_output.emplace_back(read_sequence, read_name, quality_span);

        read_sequence.clear();
        read_indices.clear();
        mutation_indices.clear();
        mutations.clear();
    }
}

int create_data_set(sharg::parser& parser) {
    parser.info.description = {
        "Simulates a genome (uniform base distribution) "
        "and a set of long reads with the configured amount of edit distance errors."
    };

    std::filesystem::path genome_path;
    std::filesystem::path read_path;

    size_t chromosome_length = 50'000'000;
    size_t num_chromosomes = 10;
    size_t base_read_length = 20'000;
    size_t num_reads = 8000;

    double error_rate = 0.07;

    size_t random_seed = 7267281;

    parser.add_option(genome_path, sharg::config{
        .short_id = 'g',
        .long_id = "genomes",
        .description = "The path for the genome.",
        .required = true,
        .validator = sharg::output_file_validator{{
            "fa", "fasta", "fna", "ffn", "fas", "faa", "mpfa", "frn"
        }}
    });

    parser.add_option(read_path, sharg::config{
        .short_id = 'r',
        .long_id = "reads",
        .description = "The path for the reads.",
        .required = true,
        .validator = sharg::output_file_validator{{ "fq", "fastq" }}
    });

    parser.add_option(chromosome_length, sharg::config{
        .short_id = 'c',
        .long_id = "chromosome-length",
        .description = "The length of each chromosome in the genome.",
        .validator = sharg::arithmetic_range_validator{0, 1'000'000'000}
    });

    parser.add_option(num_chromosomes, sharg::config{
        .short_id = 'n',
        .long_id = "num-chromosomes",
        .description = "The number of chromosomes in the genome.",
        .validator = sharg::arithmetic_range_validator{1, 50}
    });

    parser.add_option(base_read_length, sharg::config{
        .short_id = 'l',
        .long_id = "read-length",
        .description = "The length of each simulated read (might differ in the end due to indels).",
        .validator = sharg::arithmetic_range_validator{0, 1'000'000}
    });

    parser.add_option(num_reads, sharg::config{
        .short_id = 'm',
        .long_id = "num-reads",
        .description = "The number of reads in the output.",
        .validator = sharg::arithmetic_range_validator{1, 100'000}
    });

    parser.add_option(error_rate, sharg::config{
        .short_id = 'e',
        .long_id = "error-rate",
        .description =
            "The error rate of the read. Reads will usually have "
            "floor(error_rate * read_length) edit distance errors. "
            "Sometimes, there are less errors when neighboring insertions "
            "undo deletions and vice versa.",
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.add_option(random_seed, sharg::config{
        .short_id = 's',
        .long_id = "random-seed",
        .description = "The random seed to use for all randomization."
    });

    parser.parse();

    if (chromosome_length <= base_read_length) {
        spdlog::error(
            "Chromomsome length {} must be larger than read length {}",
            chromosome_length,
            base_read_length
        );

        return -1;
    }

    std::mt19937 random_generator(random_seed);

    auto const genome = create_genome(chromosome_length, num_chromosomes, random_generator);

    write_genome_to_file(genome, genome_path);

    create_and_write_reads(
        base_read_length,
        num_reads,
        error_rate,
        genome,
        chromosome_length,
        random_generator,
        read_path
    );

    return 0;
}

struct alignment_origin {
    size_t const chromosome_id;
    size_t const position;
    size_t const max_num_errors;
};

alignment_origin parse_query_id(std::string_view id) {
    std::vector<std::string> parts{};
    while (!id.empty()) {
        size_t const next_underscore_index = id.find("_");
        size_t const next_split_index = next_underscore_index == std::string_view::npos ?
            id.size() :
            next_underscore_index;

        parts.emplace_back(id.substr(0, next_split_index));
        id.remove_prefix(std::min(next_split_index + 1, id.size()));
    }

    assert(parts[0] == "id");
    assert(parts[2] == "chromosome");
    assert(parts[4] == "position");
    assert(parts[6] == "max");
    assert(parts[7] == "errors");

    return alignment_origin {
        .chromosome_id = static_cast<size_t>(std::stoi(parts[3])),
        .position = static_cast<size_t>(std::stoi(parts[5])),
        .max_num_errors = static_cast<size_t>(std::stoi(parts[8]))
    };
}

size_t parse_chromosome_id(std::string_view chromosome_id) {
    size_t const next_underscore_index = chromosome_id.find("_");
    assert(next_underscore_index != std::string_view::npos);

    chromosome_id.remove_prefix(next_underscore_index + 1);

    return static_cast<size_t>(std::stoi(std::string(chromosome_id)));
}

struct alignment_data {
    size_t const chromosome_id;
    size_t const position;
    size_t const num_errors;
};

int verify_alignments(sharg::parser& parser) {
    parser.info.description = {
        "For a previously simluated data set, verify whether an aligner mapped the reads correctly."
    };

    std::filesystem::path input_path{};
    size_t allowed_pos_diff{};

    parser.add_option(input_path, sharg::config{
        .short_id = 'a',
        .long_id = "alignments",
        .description = "The sam file that should be analyzed.",
        .required = true,
        .validator = sharg::input_file_validator{}
    });

    parser.add_option(allowed_pos_diff, sharg::config{
        .short_id = 'p',
        .long_id = "allowed-pos_diff",
        .description = "If an alignment position is shifted by this amount, it counts as found optimal."
    });

    parser.parse();

    seqan3::sam_file_input input{input_path};
    std::unordered_map<std::string, std::vector<alignment_data>> alignments_by_query_id{};

    using namespace seqan3::literals;
    for (auto const& record : input) {
        alignments_by_query_id[record.id()].push_back(alignment_data {
            .chromosome_id = parse_chromosome_id(input.header().ref_ids()[record.reference_id().value()]),
            .position = static_cast<size_t>(record.reference_position().value()),
            .num_errors = static_cast<size_t>(record.tags().get<"NM"_tag>())
        });
    }

    fmt::print("queries = [\n");

    for (auto const& [query_id, alignments] : alignments_by_query_id) {
        auto const origin = parse_query_id(query_id);

        size_t pos_diff = std::numeric_limits<uint32_t>::max();
        size_t pos_diff_higher_num_errors = std::numeric_limits<uint32_t>::max();

        for (auto const& alignment : alignments) {
            if (origin.chromosome_id != alignment.chromosome_id) {
                continue;
            }

            if (alignment.num_errors > origin.max_num_errors) {
                pos_diff_higher_num_errors = std::min(
                    std::max(origin.position, alignment.position) - std::min(origin.position, alignment.position),
                    pos_diff_higher_num_errors
                );
            } else {
                pos_diff = std::min(
                    std::max(origin.position, alignment.position) - std::min(origin.position, alignment.position),
                    pos_diff
                );
            }

            if (pos_diff == 0) {
                break;
            }
        }

        fmt::print("    {{ id = \"{}\", status = {{ ", query_id);

        if (pos_diff <= allowed_pos_diff) {
            fmt::print("FoundOptimal = {{}}");
        } else if (pos_diff == std::numeric_limits<size_t>::max() && pos_diff_higher_num_errors == std::numeric_limits<size_t>::max()) {
            fmt::print("NotFound = {{}}");
        } else {
            fmt::print(
                "FoundSuboptimal = {{ pos_diff_expected_num_errors = {}, pos_diff_higher_num_errors = {} }}",
                pos_diff,
                pos_diff_higher_num_errors
            );
        }

        fmt::print(" }} }},\n");
    }

    fmt::print("]\n");
    return 0;
}

int main(int argc, char** argv) {
    sharg::parser top_level_parser{ "simulated_dataset", argc, argv, sharg::update_notifications::off, { "create", "verify" } };

    top_level_parser.info.author = about_floxer::author;
    top_level_parser.info.email = about_floxer::email;
    top_level_parser.info.url = about_floxer::url;
    top_level_parser.info.short_description = "Simulate a genome and long reads, then verify whether an aligner mapped the reads correctly.";
    top_level_parser.info.synopsis = {
        "./simulate_simple_dataset create --genome genome.fasta --reads reads.fastq",
        "./simulate_simple_dataset verify --alignments alignments.sam",
    };
    top_level_parser.info.version = "1.0.0";
    top_level_parser.info.date = "31.05.2024";

    top_level_parser.parse();
    sharg::parser& sub_parser = top_level_parser.get_sub_parser();

    if (sub_parser.info.app_name == std::string_view{"simulated_dataset-create"}) {
        return create_data_set(sub_parser);
    } else if (sub_parser.info.app_name == std::string_view{"simulated_dataset-verify"}) {
        return verify_alignments(sub_parser);
    } else {
        spdlog::info("unknown subcommand: {}", sub_parser.info.app_name);
        return -1;
    }
}
