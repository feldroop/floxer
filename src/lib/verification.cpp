#include <math.hpp>
#include <verification.hpp>

#include <stdexcept>

namespace verification {

void query_verifier::verify() {
    switch (kind) {
        case pex::verification_kind_t::direct_full:
            direct_full_verification();
            break;

        case pex::verification_kind_t::hierarchical:
            hierarchical_verification();
            break;

        default:
            throw std::runtime_error("Internal bug in verification kind (should not happen)");
    }
}

void query_verifier::direct_full_verification() {
    if (root_was_already_verified()) {
        return;
    }

    auto const root_reference_span_config = compute_root_reference_span_config();
    internal::try_to_align_pex_node_query_with_reference_span(
        pex_tree.root(),
        reference,
        root_reference_span_config,
        query,
        orientation,
        without_cigar,
        alignments,
        stats
    );

    auto && [lock, verified_intervals] = already_verified_intervals.lock_unique();
    verified_intervals.insert(root_reference_span_config.as_half_open_interval());
}

void query_verifier::hierarchical_verification() {
    if (root_was_already_verified()) {
        return;
    }

    auto const root_reference_span_config = compute_root_reference_span_config();

    // case for when the whole PEX tree is just a single root
    if (pex_leaf_node.is_root()) {
        [[maybe_unused]] auto const outcome = internal::try_to_align_pex_node_query_with_reference_span(
            pex_leaf_node,
            reference,
            root_reference_span_config,
            query,
            orientation,
            without_cigar,
            alignments,
            stats
        );
        assert(outcome == alignment::alignment_outcome::alignment_exists);

        {
            auto && [lock, verified_intervals] = already_verified_intervals.lock_unique();
            verified_intervals.insert(root_reference_span_config.as_half_open_interval());
        }

        return;
    }

    size_t const seed_query_index_from = pex_leaf_node.query_index_from;
    auto curr_pex_node = pex_tree.get_parent_of_child(pex_leaf_node);

    while (true) {
        auto const reference_span_config = internal::compute_reference_span_start_and_length(
            anchor,
            curr_pex_node,
            seed_query_index_from,
            reference.rank_sequence.size(),
            curr_pex_node.is_root() ? extra_verification_ratio : 0.0
        );

        static constexpr size_t MAX_REF_SPAN_LENGTH_WITHOUT_CHECKING_INTERVALS = 512;
        // we ask again, because another thread might have done it
        // this is only done when the reference span is not super small
        if (
            reference_span_config.length > MAX_REF_SPAN_LENGTH_WITHOUT_CHECKING_INTERVALS
            && root_was_already_verified()
        ) {
            return;
        }

        auto const outcome = internal::try_to_align_pex_node_query_with_reference_span(
            curr_pex_node,
            reference,
            reference_span_config,
            query,
            orientation,
            without_cigar,
            alignments,
            stats
        );

        if (curr_pex_node.is_root()) {
            auto && [lock, verified_intervals] = already_verified_intervals.lock_unique();
            verified_intervals.insert(reference_span_config.as_half_open_interval());
        }

        if (outcome == alignment::alignment_outcome::no_adequate_alignment_exists || curr_pex_node.is_root()) {
            break;
        }

        curr_pex_node = pex_tree.get_parent_of_child(curr_pex_node);
    }
}

bool query_verifier::root_was_already_verified() const {
    auto const root_reference_span_config = compute_root_reference_span_config();

    auto const root_interval_to_verify_without_extra_length = root_reference_span_config
        .as_half_open_interval()
        .trim_from_both_sides(root_reference_span_config.applied_extra_verification_length_per_side);

    auto && [lock, verified_intervals] = already_verified_intervals.lock_shared();

    if (verified_intervals.contains(root_interval_to_verify_without_extra_length)) {
        // we have already verified the interval where the whole query could be found according to this anchor
        stats.add_reference_span_size_avoided_root(root_reference_span_config.length);

        return true;
    }

    return false;
}

internal::span_config query_verifier::compute_root_reference_span_config() const {
    return internal::compute_reference_span_start_and_length(
        anchor,
        pex_tree.root(),
        pex_leaf_node.query_index_from,
        reference.rank_sequence.size(),
        extra_verification_ratio
    );
}

namespace internal {

intervals::half_open_interval span_config::as_half_open_interval() const {
    return intervals::half_open_interval{
        .start = offset,
        .end = offset + length
    };
}

span_config compute_reference_span_start_and_length(
    search::anchor_t const& anchor,
    pex::pex_tree::node const& pex_node,
    size_t const leaf_query_index_from,
    size_t const full_reference_length,
    double const extra_verification_ratio
) {
    size_t const verification_interval_base_length = pex_node.length_of_query_span() + 2 * pex_node.num_errors + 1;
    size_t const extra_verification_length =
        math::floating_point_error_aware_ceil(verification_interval_base_length * extra_verification_ratio);

    int64_t const start_signed = static_cast<int64_t>(anchor.reference_position) -
        static_cast<int64_t>(leaf_query_index_from - pex_node.query_index_from) -
        static_cast<int64_t>(pex_node.num_errors) -
        static_cast<int64_t>(extra_verification_length);

    size_t const reference_span_start = start_signed >= 0 ? start_signed : 0;
    size_t const reference_span_length = std::min(
        verification_interval_base_length + 2 * extra_verification_length,
        full_reference_length - reference_span_start
    );

    return span_config {
        .offset = reference_span_start,
        .length = reference_span_length,
        .applied_extra_verification_length_per_side = extra_verification_length
    };
}

alignment::alignment_outcome try_to_align_pex_node_query_with_reference_span(
    pex::pex_tree::node const& pex_node,
    input::reference_record const& reference,
    span_config const reference_span_config,
    std::span<const uint8_t> const query,
    alignment::query_orientation const orientation,
    bool const without_cigar,
    alignment::query_alignments& alignments,
    statistics::search_and_alignment_statistics& stats
) {
    auto const this_node_query_span = query.subspan(
        pex_node.query_index_from,
        pex_node.length_of_query_span()
    );

    auto const reference_subspan = std::span<const uint8_t>(reference.rank_sequence).subspan(
        reference_span_config.offset,
        reference_span_config.length
    );

    auto mode = alignment::alignment_mode::only_verify_existance;
    if (pex_node.is_root()) {
        if (without_cigar) {
            mode = alignment::alignment_mode::verify_and_return_alignment_without_cigar;
        } else {
            mode = alignment::alignment_mode::verify_and_return_alignment_with_cigar;
        }
    }

    auto const config = alignment::alignment_config {
        .reference_span_offset = reference_span_config.offset,
        .num_allowed_errors = pex_node.num_errors,
        .orientation = orientation,
        .mode = mode
    };

    auto const alignment_result = alignment::align(
        reference_subspan,
        this_node_query_span,
        config
    );

    if (alignment_result.alignment.has_value()) {
        assert(pex_node.is_root());
        assert(alignment_result.outcome == alignment::alignment_outcome::alignment_exists);

        alignments.insert(
            std::move(alignment_result.alignment.value()),
            reference.internal_id
        );
    }

    if (pex_node.is_root()) {
        stats.add_reference_span_size_aligned_root(reference_span_config.length);
    } else {
        stats.add_reference_span_size_aligned_inner_node(reference_span_config.length);
    }

    return alignment_result.outcome;
}

} // namespace internal

} // verification
