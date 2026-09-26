#pragma once

#include <cstdint>
#include <vector>

#include "samaya/model.hpp"

namespace samaya {

struct FeasibilityJumpResult {
  bool found = false;
  std::vector<double> x;  // The feasible point when found.
  long long steps = 0;
  long long work = 0;
};

// Feasibility Jump (Luteberget and Sartor, Math. Prog. Comp. 2023): an LP-free local search for
// a point satisfying the first `rows` rows of `model` within the column bounds, integral in the
// integer columns. Each step takes a random violated row and moves one of its columns to the
// value that minimizes the weighted violation of the rows it appears in (a convex piecewise
// linear function of that value, minimized at a breakpoint); at a local minimum the weights of
// the violated rows increase. Stops after `max_work` column/row visits or `seconds`. The search
// starts from `start` where it is given and finite (clamped to the bounds, integers rounded),
// else from the bound nearest 0.
FeasibilityJumpResult feasibility_jump(const Model& model, Index rows,
                                       const std::vector<double>& lower,
                                       const std::vector<double>& upper, double seconds,
                                       long long max_work, std::uint64_t seed,
                                       const std::vector<double>* start = nullptr);

}  // namespace samaya
