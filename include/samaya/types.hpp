#pragma once

#include <cstdint>
#include <limits>

namespace samaya {

// Row / column index. 32 bits covers ~2 billion rows or columns.
using Index = std::int32_t;
// Nonzero offset. 64 bits so matrices with more than 2^31 nonzeros stay addressable.
using NnzIndex = std::int64_t;

inline constexpr double kInf = std::numeric_limits<double>::infinity();

// Values at or beyond this magnitude in input files are treated as infinite (MPS convention).
inline constexpr double kInfiniteBound = 1e30;

}  // namespace samaya
