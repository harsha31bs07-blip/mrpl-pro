#include "samaya/solver.hpp"

#include <algorithm>
#include <string>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "samaya/verify.hpp"

#ifndef SAMAYA_VERSION_STRING
#define SAMAYA_VERSION_STRING "0.0.0"
#endif

namespace samaya {

const char* version() { return SAMAYA_VERSION_STRING; }

namespace {

Status to_status(SimplexStatus s) {
  switch (s) {
    case SimplexStatus::kOptimal: return Status::kOptimal;
    case SimplexStatus::kInfeasible: return Status::kInfeasible;
    case SimplexStatus::kUnbounded: return Status::kUnbounded;
    case SimplexStatus::kIterationLimit: return Status::kIterationLimit;
    case SimplexStatus::kTimeLimit: return Status::kTimeLimit;
    case SimplexStatus::kNumericalError: return Status::kNumericalError;
  }
  return Status::kNumericalError;
}

// Runs the independent verifier on an LP outcome. Returns true when the claim is proven (or is a
// limit status, which claims nothing).
bool verify_lp(const Model& model, const LpResult& lp, Result& result) {
  VerifyReport report;
  switch (lp.status) {
    case SimplexStatus::kOptimal:
      report = verify_lp_optimality(model, lp.col_value, lp.row_dual);
      break;
    case SimplexStatus::kInfeasible:
      report = verify_infeasibility(model, lp.dual_ray);
      break;
    case SimplexStatus::kUnbounded:
      report = verify_unbounded_ray(model, lp.primal_ray);
      break;
    default:
      return true;
  }
  result.max_primal_violation = std::max(report.max_bound_violation, report.max_row_violation);
  result.max_dual_violation = report.max_dual_violation;
  result.message = report.message;
  return report.ok;
}

// Solves an LP with the dual simplex. If the outcome fails verification the solve is repeated
// with tighter tolerances and then without scaling; an outcome that never verifies is reported
// as a numerical error rather than as a solution.
void solve_lp_model(const Model& model, const Params& params, const Logger& log, Result& result) {
  if (params.lp_method == LpMethod::kBarrier || params.lp_method == LpMethod::kPdlp) {
    log.log(1, "Requested LP method not available yet; using the dual simplex");
  }
  LpSolveOptions attempts[3];
  attempts[0] = lp_options_from_params(params);
  attempts[1] = attempts[0];
  attempts[1].simplex.primal_tol = std::min(attempts[1].simplex.primal_tol, 1e-9);
  attempts[1].simplex.dual_tol = std::min(attempts[1].simplex.dual_tol, 1e-9);
  attempts[2] = attempts[1];
  attempts[2].scale = false;
  const int num_attempts = params.verify ? 3 : 1;

  LpResult lp;
  bool verified = false;
  long long iterations = 0;
  for (int a = 0; a < num_attempts; ++a) {
    lp = solve_lp(model, attempts[a], log);
    iterations += lp.iterations;
    if (!params.verify) break;
    verified = verify_lp(model, lp, result);
    if (verified) break;
    log.log(1, "Verification of %s result failed (%s); retrying with stricter settings",
            to_string(lp.status), result.message.c_str());
  }

  result.status = to_status(lp.status);
  result.simplex_iterations = iterations;
  result.verified = verified && (lp.status == SimplexStatus::kOptimal ||
                                 lp.status == SimplexStatus::kInfeasible ||
                                 lp.status == SimplexStatus::kUnbounded);
  if (params.verify && !verified) {
    result.message = std::string("result failed verification: ") + result.message;
    result.status = Status::kNumericalError;
  } else if (verified) {
    result.message.clear();
  }

  result.col_value = std::move(lp.col_value);
  result.row_activity = std::move(lp.row_activity);
  if (lp.status == SimplexStatus::kOptimal) {
    result.objective = lp.objective;
    result.dual_bound = lp.objective;
    result.row_dual = std::move(lp.row_dual);
    result.col_dual = std::move(lp.col_dual);
  }
  result.infeasibility_certificate = std::move(lp.dual_ray);
  result.unbounded_ray = std::move(lp.primal_ray);
}

}  // namespace

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

  if (std::string crossed = model.crossed_bounds(); !crossed.empty()) {
    // Infeasible by inspection; the crossed bound itself is the proof.
    result.status = Status::kInfeasible;
    result.verified = true;
    result.message = std::move(crossed);
    log.log(1, "Infeasible: %s", result.message.c_str());
    result.solve_seconds = timer.seconds();
    return result;
  }

  switch (model.problem_class()) {
    case ProblemClass::kLP:
      solve_lp_model(model, params_, log, result);
      break;
    case ProblemClass::kQP:
    case ProblemClass::kMILP:
    case ProblemClass::kMIQP:
      // Branch-and-cut and the QP solvers land in later phases (PLAN.md).
      result.status = Status::kNotImplemented;
      result.message = std::string(to_string(model.problem_class())) + " solver not implemented yet";
      break;
  }

  result.solve_seconds = timer.seconds();
  if (result.status == Status::kOptimal) {
    log.log(1, "Optimal objective %.12g after %lld simplex iterations, %.3f s%s", result.objective,
            result.simplex_iterations, result.solve_seconds,
            result.verified ? " (verified)" : "");
  } else {
    log.log(1, "Status %s after %lld simplex iterations, %.3f s%s%s", to_string(result.status),
            result.simplex_iterations, result.solve_seconds, result.message.empty() ? "" : ": ",
            result.message.c_str());
  }
  return result;
}

}  // namespace samaya
