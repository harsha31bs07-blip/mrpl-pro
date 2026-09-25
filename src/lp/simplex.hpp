#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "core/log.hpp"
#include "linalg/basis_factor.hpp"
#include "samaya/sparse_matrix.hpp"
#include "samaya/types.hpp"

namespace samaya {

// LP in computational form:  min cost'v  s.t.  [A  -I] v = 0,  lower <= v <= upper,
// where v = (x, r) holds the n structural variables followed by the m row activities
// ("logicals"). Logical n+i carries the bounds of row i.
struct LpProblem {
  Index m = 0;
  Index n = 0;
  SparseMatrix A;   // m x n.
  SparseMatrix At;  // Transpose of A, i.e. A in row-wise form.
  std::vector<double> cost;   // n + m.
  std::vector<double> lower;  // n + m.
  std::vector<double> upper;  // n + m.
};

enum class SimplexStatus : std::uint8_t {
  kOptimal,
  kInfeasible,
  kUnbounded,
  kIterationLimit,
  kTimeLimit,
  kNumericalError,
};

const char* to_string(SimplexStatus status);

enum class VarStatus : std::uint8_t { kBasic, kAtLower, kAtUpper, kAtZero };

// Work counters, mainly for tests and logs: a healthy dual simplex does nearly all of its work in
// dual iterations, with few cost shifts and little primal cleanup.
struct SimplexStats {
  long long dual_iterations = 0;
  long long primal_iterations = 0;
  long long bound_flips = 0;
  long long cost_shifts = 0;
  long long refactorizations = 0;
  long long rebuilds_on_mismatch = 0;
};

struct SimplexOptions {
  double primal_tol = 1e-7;
  double dual_tol = 1e-7;
  double pivot_tol = 1e-7;
  long long max_iterations = -1;  // -1: automatic, proportional to the problem size.
  double time_limit = kInf;
  bool perturb = true;
  int refactor_interval = 100;
  std::uint64_t seed = 1;
};

// Bounded dual simplex with dual steepest-edge pricing, a Harris two-pass ratio test, cost
// perturbation and shifting against degeneracy, and a primal simplex for cleanup.
//
// Phase 1 uses the auxiliary-problem approach: bounds are replaced by boxes ([0,1] for
// lower-bounded, [-1,0] for upper-bounded, [-1000,1000] for free variables, [0,0] otherwise) and
// the dual simplex minimizes the sum of dual infeasibilities. If the original problem is dual
// infeasible, a zero-cost dual simplex decides primal feasibility and the primal simplex then
// proves unboundedness.
class Simplex {
 public:
  Simplex(const LpProblem& lp, const SimplexOptions& options, const Logger& log);

  SimplexStatus solve();
  // Starts from the given variable statuses (n + m, exactly m basic) instead of the slack basis:
  // bounds are restored, dual infeasibilities are removed by bound flips and cost shifts, and the
  // dual simplex continues. Falls back to solve() if the statuses do not form a usable basis.
  SimplexStatus solve(const std::vector<VarStatus>& start);

  // All vectors are in the (scaled) space of the LpProblem.
  const std::vector<double>& values() const { return x_; }         // n + m.
  const std::vector<double>& duals() const { return y_; }          // m.
  const std::vector<double>& reduced_costs() const { return d_; }  // n + m.
  const std::vector<VarStatus>& status() const { return status_; }
  // Farkas ray y (length m) after kInfeasible: 0 is outside the range of
  // y'A x - y'r over the variable bounds.
  const std::vector<double>& dual_ray() const { return dual_ray_; }
  // Improving direction (length n + m) after kUnbounded.
  const std::vector<double>& primal_ray() const { return primal_ray_; }
  long long iterations() const { return iterations_; }
  const SimplexStats& stats() const { return stats_; }

  // Dual steepest-edge weights as maintained by the updates, and recomputed from scratch as
  // ||e_r' B^-1||^2 (m btran solves; for tests).
  const std::vector<double>& dse_weights() const { return dse_weight_; }
  std::vector<double> exact_dse_weights() const;

 private:
  enum class LoopResult : std::uint8_t { kDone, kLostFeasibility };

  void reset();
  bool rebuild();
  void compute_primal();
  void compute_dual();
  void compute_pivot_row();
  void set_value_from_status(Index j);
  void place_nonbasic_dual_feasible();
  bool correct_dual_infeasibilities();
  void perturb_costs(double base);
  Index count_dual_infeasibilities_unboxed() const;
  double max_primal_infeasibility() const;
  double max_dual_infeasibility() const;
  bool limit_reached(SimplexStatus& status) const;

  SimplexStatus dual_phase1();
  SimplexStatus phase2();
  SimplexStatus solve_dual_infeasible();
  SimplexStatus dual_loop();
  SimplexStatus primal_loop(LoopResult& result);

  Index choose_leaving_row() const;
  Index dual_ratio_test(double direction, double slope);
  void apply_flips();
  Index choose_entering_column() const;

  template <typename F>
  void for_column(Index j, F&& f) const;
  void load_column(Index j, std::vector<double>& out) const;

  const LpProblem& lp_;
  SimplexOptions options_;
  const Logger& log_;
  Index m_;
  Index n_;
  Index nt_;
  long long max_iterations_;
  Timer timer_;

  BasisFactor factor_;
  std::vector<Index> basic_;
  std::vector<Index> position_;
  std::vector<VarStatus> status_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> d_;
  std::vector<double> cost_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<double> dse_weight_;
  std::vector<double> devex_weight_;  // Primal Devex reference weights (n + m).

  std::vector<double> rho_;
  std::vector<double> tau_;
  std::vector<double> alpha_row_;
  std::vector<Index> row_nz_;       // Indices of the nonzeros of alpha_row_.
  std::vector<char> in_row_nz_;
  std::vector<Index> candidates_;   // Ratio test workspace.
  std::vector<Index> flips_;        // Boxed variables flipped by the long-step ratio test.
  std::vector<double> alpha_col_;
  std::vector<double> spike_;
  std::vector<double> dual_ray_;
  std::vector<double> primal_ray_;

  long long iterations_ = 0;
  SimplexStats stats_;
  std::mt19937_64 rng_;
};

}  // namespace samaya
