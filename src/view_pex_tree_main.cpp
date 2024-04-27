#include <about_floxer.hpp>
#include <pex.hpp>

#include <cmath>
#include <limits>

#include <sharg/all.hpp>
#include <spdlog/fmt/fmt.h>

int main(int argc, char** argv) {
    sharg::parser parser{ "view-pex-tree", argc, argv, sharg::update_notifications::off };

    parser.info.author = about_floxer::author;
    parser.info.description = { "Given PEX tree config, print resulting tree to stdout in DOT format." };
    parser.info.email = about_floxer::email;
    parser.info.url = about_floxer::url;
    parser.info.short_description = "View PEX tree in DOT format";
    parser.info.synopsis = {
        "./view-pex-tree --query-length 100 --query-errors 7 --seed-errors 2",
        "./view-pex-tree --query-length 100 --query-error-probability 0.07 --seed-errors 2",
    };
    parser.info.version = "1.0.0";
    parser.info.date = "27.04.2024";

    size_t total_query_length;
    size_t const given_query_num_errors_default = std::numeric_limits<size_t>::max();
    size_t given_query_num_errors = given_query_num_errors_default;
    double const query_error_probability_default = NAN;
    double query_error_probability = query_error_probability_default;
    size_t leaf_max_num_errors = 2;

    parser.add_option(total_query_length, sharg::config{
        .short_id = 'q', 
        .long_id = "query-length", 
        .description = "The length of the query for which the PEX tree should be built.",
        .required = true
    });

    parser.add_option(given_query_num_errors, sharg::config{
        .short_id = 'e', 
        .long_id = "query-errors", 
        .description = "The number of errors allowed in each query. This is only used if no error "
            "probability is given. Either this or an error probability must be given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0, 4096}
    });

    parser.add_option(query_error_probability, sharg::config{
        .short_id = 'p', 
        .long_id = "query-error-probability",
        .description = "The error probability in the queries, per base. If this is given, it is used "
            "rather than the fixed number of errors. Either this or a fixed number of errors must be "
            "given.",
        .default_message = "no default",
        .validator = sharg::arithmetic_range_validator{0.00001, 0.99999}
    });

    parser.add_option(leaf_max_num_errors, sharg::config{
        .short_id = 's', 
        .long_id = "seed-errors", 
        .description = "The number of errors in the leaves of the PEX tree that are used as seeds. "
            "The sequences will be searched with this parameter using the FM-index.",
        .validator = sharg::arithmetic_range_validator{0, 10}
    });

    parser.parse();

    if (
        given_query_num_errors == given_query_num_errors_default &&
        query_error_probability == query_error_probability_default
    ) {
        throw std::runtime_error(
            "Either a fixed number of errors in the query or an error probability must be given."
        );
    }

    auto const query_num_errors = !std::isnan(query_error_probability) ?
        static_cast<size_t>(
            std::ceil(total_query_length * query_error_probability)
        ) :
        given_query_num_errors;


    auto const config = pex::pex_tree_config{
        .total_query_length = total_query_length,
        .query_num_errors = query_num_errors,
        .leaf_max_num_errors = leaf_max_num_errors
    };

    auto const tree = pex::pex_tree(config);

    fmt::print("{}", tree.dot_statement());
    
    return 0;
}