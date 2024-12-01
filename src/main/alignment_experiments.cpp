#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>


#include <string>
#include <string_view>
#include <array>   // std::array
#include <utility> // std::index_sequence

template <std::size_t...Idxs>
constexpr auto substring_as_array(std::string_view str, std::index_sequence<Idxs...>)
{
  return std::array{str[Idxs]..., '\n'};
}

template <typename T>
constexpr auto type_name_array()
{
  constexpr auto prefix   = std::string_view{"with T = "};
  constexpr auto suffix   = std::string_view{"]"};
  constexpr auto function = std::string_view{__PRETTY_FUNCTION__};

  constexpr auto start = function.find(prefix) + prefix.size();
  constexpr auto end = function.rfind(suffix);

  static_assert(start < end);

  constexpr auto name = function.substr(start, (end - start));
  return substring_as_array(name, std::make_index_sequence<name.size()>{});
}

template <typename T>
struct type_name_holder {
  static inline constexpr auto value = type_name_array<T>();
};

template <typename T>
constexpr auto type_name() -> std::string_view
{
  constexpr auto& value = type_name_holder<T>::value;
  return std::string_view{value.data(), value.size()};
}

template <
    bool compute_score_matrix_,
    typename database_t,
    typename query_t,
    typename align_cfg_t,
    typename word_t,
    typename is_semi_global_t,
    typename traits_t =
        seqan3::detail::default_edit_distance_trait_type<database_t, query_t, align_cfg_t, is_semi_global_t, word_t>>
struct edit_traits_type : traits_t
{
    static constexpr bool compute_score_matrix = compute_score_matrix_;
    static constexpr bool compute_matrix = compute_score_matrix || traits_t::compute_trace_matrix;
};

class aligner {
public:
    aligner();

    void align(std::span<const uint8_t> const seq1, std::span<const uint8_t> const seq2);

    static const config_base = auto config = seqan3::align_cfg::method_global{
        seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
        seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
        seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}
    }
    | seqan3::align_cfg::edit_scheme
    | seqan3::align_cfg::min_score{-3}; // problem: this changes for each invocation

    static auto const config_only_score = config_base | seqan3::align_cfg::output_score{};
    static auto const config_full_output = config_base
        | seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_begin_position{}
        | seqan3::align_cfg::output_alignment{};

    using only_score_config_with_score_type_t = decltype(config_only_score);
    using full_output_config_with_score_type_t = decltype(config_full_output);

    using sequence_t = std::span<const uint8_t>;

    using only_score_alignment_result_value_t = typename seqan3::detail::
        align_result_selector<sequence_t, sequence_t, std::remove_reference_t<only_score_config_with_score_type_t>>::type;

    using only_score_alignment_result_t = seqan3::alignment_result<only_score_alignment_result_value_t>;

    auto config_only_score_with_result_type = config_only_score | seqan3::align_cfg::detail::result_type<only_score_alignment_result_t>{};
    auto config_full_output_with_result_type = config_full_output | seqan3::align_cfg::detail::result_type<alignment_result_t>{};

    using only_score_align_config_with_result_type_t = decltype(config_only_score_with_result_type);

    using only_score_edit_traits = edit_traits_type<true, sequence_t, sequence_t, only_score_align_config_with_result_type_t>;

    using only_score_algorithm_t =
        seqan3::detail::edit_distance_unbanded<sequence_t, sequence_t, only_score_align_config_with_result_type_t, only_score_edit_traits>;
};

aligner::aligner() {

}

void align(std::span<const uint8_t> const seq1, std::span<const uint8_t> const seq2) {
    auto sequences = std::tie(seq1, seq2);
    auto seq_view = std::views::single(sequences) | seqan3::detail::all;
}

int main() {
    std::vector<uint8_t> const seq1 {
        0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,1,
        0,0,0,0,0,0,0,0,0,0
    };
    std::vector<uint8_t> const seq2 {
        1,1,1,1,1,1
    };
    // std::vector<uint8_t> seq3 {
    //     1,1,1,1,1,2,1,1,1
    // };

    auto config = seqan3::align_cfg::method_global{
        seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
        seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
        seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}
    }
    | seqan3::align_cfg::edit_scheme
    | seqan3::align_cfg::min_score{-3}
    | seqan3::align_cfg::output_score{}
    | seqan3::align_cfg::output_begin_position{}
    | seqan3::align_cfg::output_alignment{};

    auto const seq1_span = std::span<const uint8_t>(seq1);
    auto const seq2_span = std::span<const uint8_t>(seq2);

    // auto alignment_results = seqan3::align_pairwise(sequences, config);
    // auto alignment = *alignment_results.begin();

    // seqan3::debug_stream << alignment << '\n';
    // seqan3::debug_stream << alignment.alignment() << '\n';
    // seqan3::debug_stream << seqan3::cigar_from_alignment(alignment.alignment(), {}, true) << '\n';

    seqan3::debug_stream << test_res << '\n';
    seqan3::debug_stream << test_res.alignment() << '\n';
    seqan3::debug_stream << seqan3::cigar_from_alignment(test_res.alignment(), {}, true) << '\n';

    return 0;
}
