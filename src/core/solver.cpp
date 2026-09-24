#include "samaya/solver.hpp"

#include <string>

#include "core/log.hpp"

#ifndef SAMAYA_VERSION_STRING
#define SAMAYA_VERSION_STRING "0.0.0"
#endif

namespace samaya {

const char* version() { return SAMAYA_VERSION_STRING; }

Result Solver::solve(const Model& model) const {
  const Timer timer;
  const Logger log(params_.log_level);
  Result result;

  if (std::string error = model.validate(); !error.empty()) {
    result.status = Status::kInvalidModel;
    result.message = std::move(error);
    log.log(1, "Invalid model: %s", result.message.c_str());
    result.solve_seconds = timer.seconds();
    return result;
  }

  log.log(1, "samaya %s: %s with %d rows, %d cols, %lld nonzeros", version(),
          to_string(model.problem_class()), model.num_rows(), model.num_cols(),
          static_cast<long long>(model.A.nnz()));

  // Dispatch point for presolve -> {simplex, barrier, PDLP, branch-and-cut} -> postsolve -> verify.
  // The algorithms land in PLAN.md phases 1-4; until then report that honestly.
  result.status = Status::kNotImplemented;
  result.message = std::string(to_string(model.problem_class())) + " solver not implemented yet";
  log.log(1, "%s", result.message.c_str());
  result.solve_seconds = timer.seconds();
  return result;
}

}  // namespace samaya
