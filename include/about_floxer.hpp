#pragma once

#include <string>

namespace about_floxer {

static const std::string program_name = "floxer";
static const std::string version = "0.1.0";
static const std::string author = "Felix Leander Droop";
static const std::string email = "felix.droop@fu-berlin.de";
static const std::string version_date = "09.03.2024";
static const std::string short_description = "A longread aligner using pidgeonhole exact (PEX) seeds and the FM-index.";
static const std::string long_description = 
    "floxer is an exact longread aligner using FM-index search with optimal search schemes, "
    "the PEX (pidgeonhole exact) hierarchical verification scheme and "
    "a highly parallel alignment implementation.";
static const std::string url = "https://github.com/feldroop/floxer";

} // namespace floxer
