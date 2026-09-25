#pragma once

#include <map>
#include <memory>
#include <vector>

#include "core/log.hpp"
#include "linalg/scaling.hpp"
#include "lp/simplex.hpp"
#include "samaya/model.hpp"
#include "samaya/status.hpp"

namespace samaya {

struct MipOptions {
  double time_limit = kInf;
  long long node_limit = -1;  // -1: unlimited.
  double rel_gap = 1e-4;
  double abs_gap = 1e-6;
  double integrality_tol = 1e-6;
  double feasibility_tol = 1e-6;  // Row and bound tolerance for accepting a solution.
  // Pseudocosts count as reliable after this many observations in each direction; before that a
  // candidate is strong-branched.
  int reliability = 4;
  int max_strong_branching = 8;       // Strong-branching candidates per node.
  long long strong_iterations = 200;  // Dual simplex iterations per strong-branching child.
};

struct MipOutcome {
  // kOptimal (gap closed), kInfeasible, kInfeasibleOrUnbounded, kUnbounded, kTimeLimit,
  // kNodeLimit or kNumericalError (a node LP could not be solved, so the search is incomplete).
  Status status = Status::kNotSolved;
  double objective = kInf;  // Of `x`, in the model's sense; +-inf without a solution.
  double bound = -kInf;     // Best proven bound, in the model's sense.
  std::vector<double> x;    // Best integer solution found (empty if none).
  long long nodes = 0;
  long long lp_iterations = 0;
  long long strong_branching_iterations = 0;
  int heuristic_solutions = 0;
};

// LP-based branch-and-bound for mixed-integer linear programs.
//
// Node relaxations are solved by the dual simplex, warm-started from the parent's optimal basis,
// on one scaled copy of the LP whose column bounds change from node to node. Each node first runs
// bound propagation on the integer columns. Branching uses pseudocosts with reliability
// initialization by strong branching (the product score); strong branching that proves a child
// infeasible or cut off tightens the node instead. Nodes are selected by best bound, with plunging
// (a depth-first dive into a child) while the dive stays promising. Heuristics: lock-based simple
// rounding at every node, and rounding followed by an LP over the continuous columns at the root
// and periodically. A candidate solution is accepted only if it satisfies the rows, bounds and
// integrality within the tolerances.
class BranchAndBound {
 public:
  BranchAndBound(const Model& model, const MipOptions& options, const Logger& log);
  ~BranchAndBound();

  MipOutcome solve();

 private:
  struct BoundChange {
    Index col;
    double lower;
    double upper;
  };
  struct Node {
    double bound = -kInf;  // Lower bound (minimization) inherited from the parent.
    int depth = 0;
    std::vector<BoundChange> path;  // Applied in order on top of the root bounds.
    std::shared_ptr<const std::vector<VarStatus>> basis;
    Index branch_col = -1;  // The branching that created this node, for pseudocosts.
    bool branch_up = false;
    double branch_distance = 0.0;
  };
  enum class NodeResult : std::uint8_t { kPruned, kBranched, kStopped, kFailed, kUnbounded };

  // Relaxation.
  void set_bound(Index j, double lower, double upper);
  void apply_node_bounds(const Node& node);
  SimplexStatus solve_relaxation(const std::vector<VarStatus>* start, long long iteration_limit);
  double relaxation_objective() const;
  std::vector<VarStatus> current_basis() const;

  // Search.
  NodeResult process_node(Node& node, std::vector<Node>& children);
  bool propagate(std::vector<Index> changed, std::vector<BoundChange>* record);
  Index select_branching(const std::vector<Index>& fractional, const std::vector<double>& x,
                         double objective, const std::vector<VarStatus>& basis,
                         bool& node_infeasible, BoundChange& tighten, bool& has_tighten);
  double pseudocost(Index j, bool up) const;
  void record_pseudocost(Index j, bool up, double gain, double distance);
  double cutoff() const;
  double effective_bound(double objective) const;
  double best_bound(double extra) const;
  bool time_up() const;

  // Solutions.
  bool is_fractional(double v) const;
  bool try_solution(std::vector<double> x);
  void simple_rounding(const std::vector<double>& x);
  void round_and_solve(const std::vector<double>& x, const std::vector<VarStatus>& basis);

  const Model& model_;
  MipOptions options_;
  const Logger& log_;
  Timer timer_;
  Index m_;
  Index n_;
  double sense_;
  std::vector<double> cost_;  // Minimization costs.
  double offset_;
  bool integral_objective_ = false;
  std::vector<Index> integers_;
  SparseMatrix At_;
  std::vector<int> down_locks_;
  std::vector<int> up_locks_;

  Logger lp_log_{0};
  Scaling scaling_;
  LpProblem lp_;
  std::unique_ptr<Simplex> simplex_;
  std::vector<double> x_;  // Relaxation solution in the original space.
  SimplexStatus lp_status_ = SimplexStatus::kNumericalError;

  std::vector<double> root_lower_;
  std::vector<double> root_upper_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::vector<Index> touched_;
  std::vector<char> is_touched_;
  std::vector<VarStatus> root_basis_;

  std::vector<char> row_mark_;  // Propagation queue membership.

  std::vector<double> pc_sum_[2];  // [0] down, [1] up: summed gain per unit change.
  std::vector<int> pc_count_[2];
  double pc_total_sum_[2] = {0.0, 0.0};
  long long pc_total_count_[2] = {0, 0};

  std::multimap<double, Node> open_;
  double incumbent_value_ = kInf;
  std::vector<double> incumbent_;
  double pruned_bound_ = kInf;  // Smallest bound among nodes pruned by the cutoff.
  MipOutcome outcome_;
  bool incomplete_ = false;
};

}  // namespace samaya
