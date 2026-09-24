#pragma once

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "samaya/model.hpp"
#include "samaya/params.hpp"
#include "samaya/status.hpp"

namespace samaya {

struct Result {
  Status status = Status::kNotSolved;
  std::string message;

  // Objective value of the returned solution and the best proven bound (equal for LP/QP optima).
  double objective = std::numeric_limits<double>::quiet_NaN();
  double dual_bound = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> col_value;     // x
  std::vector<double> row_activity;  // A x
  std::vector<double> row_dual;      // y (LP/QP only)
  std::vector<double> col_dual;      // reduced costs (LP/QP only)

  double solve_seconds = 0.0;
  long long simplex_iterations = 0;
  long long barrier_iterations = 0;
  long long nodes = 0;
};

class Solver {
 public:
  explicit Solver(Params params = {}) : params_(std::move(params)) {}

  const Params& params() const { return params_; }
  Params& params() { return params_; }

  Result solve(const Model& model) const;

 private:
  Params params_;
};

// Semantic version of the library, e.g. "0.1.0".
const char* version();

}  // namespace samaya
