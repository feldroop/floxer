#pragma once

#include <fmindex-collection/fmindex/BiFMIndex.h>
#include <fmindex-collection/occtable/EPR.h>

size_t constexpr Sigma = 6; // DNA + N + Sentinel
using Table = fmindex_collection::occtable::EprV2_16<Sigma>;
using fmindex = fmindex_collection::BiFMIndex<Table>;
