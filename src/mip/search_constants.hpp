#pragma once

// Constants shared by the sequential and the parallel tree search.

namespace samaya::search {

// Without an incumbent a plunge dives this deep; with one it continues while the child's bound
// stays within this fraction of the gap above the best open bound.
constexpr int kMaxPlungeDepth = 1000;
constexpr double kPlungeGapFraction = 0.3;
// Iteration budget of the LP over the continuous columns after the integers are rounded and
// fixed (round_and_solve at nodes and in the heuristics).
constexpr long long kRoundAndSolveIterations = 2000;

}  // namespace samaya::search
