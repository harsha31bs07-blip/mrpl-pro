#pragma once

// Reference MILP solver for tests: depth-first branch-and-bound on the dense reference simplex,
// branching on the first fractional integer column. Exponential, only for small bounded models;
// it shares no code with the solver's search.

#include <cmath>
#include <vector>

#include "reference_lp.hpp"
#include "samaya/model.hpp"

namespace samaya::test {

struct ReferenceMilpResult {
  enum class Status { kOptimal, kInfeasible, kUnbounded, kNodeLimit } status = Status::kInfeasible;
  double objective = 0.0;
  long long nodes = 0;
};

inline ReferenceMilpResult reference_milp(const Model& model, long long node_limit = 200000) {
  ReferenceMilpResult result;
  const double sense = model.sense == ObjSense::kMaximize ? -1.0 : 1.0;
  double best = kInf;  // Minimization.
  struct Box {
    std::vector<double> lower;
    std::vector<double> upper;
  };
  std::vector<Box> stack{{model.col_lower, model.col_upper}};
  Model node = model;
  while (!stack.empty()) {
    if (++result.nodes > node_limit) {
      result.status = ReferenceMilpResult::Status::kNodeLimit;
      return result;
    }
    Box box = std::move(stack.back());
    stack.pop_back();
    node.col_lower = box.lower;
    node.col_upper = box.upper;
    const ReferenceResult lp = ReferenceLp(node).solve();
    if (lp.status == ReferenceResult::Status::kInfeasible) continue;
    if (lp.status == ReferenceResult::Status::kUnbounded) {
      result.status = ReferenceMilpResult::Status::kUnbounded;
      return result;
    }
    const double value = sense * lp.objective;
    if (value >= best - 1e-9 * (1.0 + std::fabs(best))) continue;
    Index branch = -1;
    for (Index j = 0; j < model.num_cols(); ++j) {
      if (model.col_type[j] == VarType::kInteger &&
          std::fabs(lp.x[j] - std::round(lp.x[j])) > 1e-7) {
        branch = j;
        break;
      }
    }
    if (branch < 0) {
      best = value;
      continue;
    }
    const double v = lp.x[branch];
    Box up = box;
    up.lower[branch] = std::ceil(v);
    box.upper[branch] = std::floor(v);
    stack.push_back(std::move(up));
    stack.push_back(std::move(box));
  }
  if (best < kInf) {
    result.status = ReferenceMilpResult::Status::kOptimal;
    result.objective = sense * best;
  }
  return result;
}

}  // namespace samaya::test
