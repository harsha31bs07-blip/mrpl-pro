// Primal heuristics of the branch-and-bound: the feasibility pump, LP diving and the RENS and
// RINS sub-MIPs. They share the search's scaled LP and simplex; every candidate goes through
// try_solution, so a heuristic can only ever add verified solutions.

#include <algorithm>
#include <cmath>
#include <functional>

#include "mip/branch_and_bound.hpp"
#include "mip/feasibility_jump.hpp"
#include "presolve/presolve.hpp"

namespace samaya {

namespace {

// Budgets. Root dives and the pump may use this multiple of the root LP's iterations (at least
// the minimum). Tree dives start every kDiveFrequency nodes, backing off to kMaxDiveInterval
// while they fail, as long as all heuristic LP iterations stay below kTreeIterationFraction of
// the node LP iterations plus kTreeIterationBase.
constexpr double kRootIterationFactor = 1.0;
constexpr long long kMinHeuristicIterations = 1000;
constexpr long long kDiveFrequency = 10;
constexpr long long kMaxDiveInterval = 1000;
constexpr double kTreeIterationFraction = 0.1;
constexpr long long kTreeIterationBase = 5000;
// Feasibility pump: iterations, decay of the original objective's weight, and the number of
// columns flipped when the rounding repeats (between the two values, at random).
constexpr int kPumpMaxIterations = 100;
constexpr double kPumpAlphaDecay = 0.9;
constexpr int kPumpMinFlips = 10;
constexpr int kPumpMaxFlips = 30;
constexpr std::size_t kPumpHistory = 100;
// Sub-MIPs run only if they fix at least this fraction of the integer columns, with this node
// limit, at most this fraction of the remaining time and at most this many seconds.
constexpr double kRensMinFixed = 0.2;
constexpr double kRinsMinFixed = 0.3;
constexpr long long kSubMipNodes = 500;
constexpr double kSubMipTimeFraction = 0.1;
constexpr double kSubMipMaxSeconds = 30.0;
constexpr long long kFirstRinsNode = 100;
// Feasibility Jump: at most this many seconds (and this share of the time limit), and this many
// row/column visits per matrix nonzero plus a minimum.
constexpr double kJumpMaxSeconds = 5.0;
constexpr double kJumpTimeFraction = 0.05;
constexpr double kJumpWorkPerNonzero = 20000.0;
constexpr long long kJumpMinWork = 10000000;
// RINS runs only while all sub-MIPs together took at most this share of the elapsed time plus
// this many seconds (refsched_8c_6p_3u_26t spent 26 of 60 s in sub-MIPs without it; 46 s with).
constexpr double kSubMipTimeShare = 0.1;
constexpr double kSubMipFreeSeconds = 1.0;
// Keeps the pseudocost ratio finite when a direction costs nothing.
constexpr double kScoreFloorDive = 1e-6;

}  // namespace

// Feasibility Jump before the root LP (as in HiGHS and Xpress): LP-free, so it can find a first
// solution where the LP-based heuristics do not (general integers, equality-heavy models).
void BranchAndBound::run_feasibility_jump() {
  double seconds = kJumpMaxSeconds;
  if (std::isfinite(options_.time_limit)) {
    seconds = std::min(seconds, kJumpTimeFraction * options_.time_limit);
  }
  const auto work = static_cast<long long>(kJumpWorkPerNonzero *
                                           static_cast<double>(model_.A.nnz())) +
                    kJumpMinWork;
  const Timer timer;
  const FeasibilityJumpResult jump =
      feasibility_jump(model_, original_rows_, lower_, upper_, seconds, work, 12345u);
  const bool improved = jump.found && try_solution(jump.x);
  if (improved) {
    ++outcome_.heuristic_solutions;
    jump_incumbent_value_ = incumbent_value_;
  }
  log_.log(1, "MIP feasibility jump: %s after %lld steps, %.2f s%s",
           jump.found ? "found a solution" : "no solution", jump.steps, timer.seconds(),
           jump.found && !improved ? " (not accepted)" : "");
}

void BranchAndBound::undo_bounds(std::size_t mark) {
  const bool logging = logging_undo_;
  logging_undo_ = false;
  while (undo_log_.size() > mark) {
    const BoundChange c = undo_log_.back();
    undo_log_.pop_back();
    set_bound(c.col, c.lower, c.upper);
  }
  logging_undo_ = logging;
}

void BranchAndBound::run_heuristics(const Node& node, const std::vector<double>& x,
                                    const std::vector<VarStatus>& basis) {
  if (integers_.empty()) return;
  const long long used = outcome_.heuristic_lp_iterations;
  if (node.depth == 0) {
    const long long budget = std::max(
        kMinHeuristicIterations,
        static_cast<long long>(kRootIterationFactor * static_cast<double>(outcome_.lp_iterations)));
    if (incumbent_.empty() || incumbent_value_ == jump_incumbent_value_) {
      feasibility_pump(x, basis, budget);
    }
    for (const DiveRule rule : {DiveRule::kCoefficient, DiveRule::kFractional,
                                DiveRule::kPseudocost, DiveRule::kGuided}) {
      if (time_up()) return;
      if (rule == DiveRule::kGuided && incumbent_.empty()) continue;
      dive(rule, x, basis, budget);
    }
    if (options_.sub_mip_heuristics && !time_up()) rens(x);
    log_.log(1, "MIP root heuristics: %d solutions, incumbent %.10g, %lld LP iterations, %.1f s",
             outcome_.heuristic_solutions, sense_ * incumbent_value_,
             outcome_.heuristic_lp_iterations - used, timer_.seconds());
    return;
  }
  if (outcome_.nodes >= next_dive_node_) {
    const double node_iterations =
        static_cast<double>(outcome_.lp_iterations - outcome_.heuristic_lp_iterations);
    const long long allowed = kTreeIterationBase +
                              static_cast<long long>(kTreeIterationFraction * node_iterations) -
                              outcome_.heuristic_lp_iterations;
    if (allowed > 0) {
      DiveRule rule = static_cast<DiveRule>(next_dive_rule_++ % 4);
      if (rule == DiveRule::kGuided && incumbent_.empty()) rule = DiveRule::kFractional;
      const int before = outcome_.heuristic_solutions;
      dive(rule, x, basis, allowed);
      // Back off after dives that find nothing: the interval doubles up to the maximum and a
      // success resets it.
      dive_interval_ = outcome_.heuristic_solutions > before
                           ? kDiveFrequency
                           : std::min(kMaxDiveInterval, 2 * dive_interval_);
    }
    next_dive_node_ = outcome_.nodes + dive_interval_;
  }
  if (options_.sub_mip_heuristics && !incumbent_.empty() && outcome_.nodes >= next_rins_node_ &&
      incumbent_value_ < rins_incumbent_ &&
      sub_mip_seconds_ <= kSubMipTimeShare * timer_.seconds() + kSubMipFreeSeconds) {
    rins_incumbent_ = incumbent_value_;
    next_rins_node_ = std::max(kFirstRinsNode, 2 * outcome_.nodes);
    rins(x);
  }
}

// LP diving: fix one fractional column per LP in the rule's direction, propagate and re-solve;
// on an infeasible or cut-off LP try the other direction once, then give up.
void BranchAndBound::dive(DiveRule rule, const std::vector<double>& x0,
                          const std::vector<VarStatus>& basis0, long long budget) {
  const long long start = outcome_.lp_iterations;
  undo_log_.clear();
  logging_undo_ = true;
  std::vector<double> x = x0;
  std::vector<VarStatus> basis = basis0;
  long long steps = 0;
  for (;;) {
    Index best = -1;
    bool best_up = false;
    double best_score = kInf;
    bool all_roundable = true;
    for (const Index j : integers_) {
      const double v = x[j];
      if (lower_[j] == upper_[j] || !is_fractional(v)) continue;
      if (std::ceil(v) <= lower_[j] || std::floor(v) >= upper_[j]) {
        // Within the primal tolerance of an integral bound: rounding changes no bound, so a dive
        // step on it would repeat forever. Snap it instead.
        x[j] = std::clamp(std::round(v), lower_[j], upper_[j]);
        continue;
      }
      const double f = v - std::floor(v);
      const int locks_down = down_locks_[j];
      const int locks_up = up_locks_[j];
      if (locks_down > 0 && locks_up > 0) all_roundable = false;
      bool up = f > 0.5;
      double score = std::min(f, 1.0 - f);
      switch (rule) {
        case DiveRule::kFractional: break;
        case DiveRule::kCoefficient:
          // Fewest locks first: the direction that endangers the fewest rows.
          if (locks_down != locks_up) up = locks_up < locks_down;
          score += std::min(locks_down, locks_up);
          break;
        case DiveRule::kPseudocost: {
          const double down_cost = pseudocost(j, false) * f;
          const double up_cost = pseudocost(j, true) * (1.0 - f);
          up = f > 0.7 || (f >= 0.3 && up_cost < down_cost);
          const double chosen = up ? up_cost : down_cost;
          const double other = up ? down_cost : up_cost;
          // The most decided column first: the other direction costs much more.
          score = -(other + kScoreFloorDive) / (chosen + kScoreFloorDive);
          break;
        }
        case DiveRule::kGuided:
          up = incumbent_[j] > v;
          score = std::fabs(v - incumbent_[j]);
          break;
      }
      if (score < best_score) {
        best_score = score;
        best = j;
        best_up = up;
      }
    }
    if (best < 0) {
      if (try_solution(x)) {
        ++outcome_.heuristic_solutions;
      } else {
        round_and_solve(x, basis);  // Integral within the tolerance but not exactly.
      }
      break;
    }
    if (all_roundable) {
      simple_rounding(x);
      break;
    }
    // Each step counts at least one unit of effort: a step that moves a column by one without an
    // LP iteration (a bound flip) could otherwise walk across a huge integer domain for free.
    ++steps;
    if (time_up() || outcome_.lp_iterations - start + steps >= budget) break;

    bool ok = false;
    bool abort = false;
    for (int attempt = 0; attempt < 2 && !ok && !abort; ++attempt) {
      const bool up = attempt == 0 ? best_up : !best_up;
      const std::size_t mark = undo_log_.size();
      const double v = x[best];
      if (up) {
        set_bound(best, std::ceil(v), upper_[best]);
      } else {
        set_bound(best, lower_[best], std::floor(v));
      }
      if (propagate({best}, nullptr)) {
        const long long left = std::max(1LL, budget - (outcome_.lp_iterations - start + steps));
        const SimplexStatus status = solve_relaxation(&basis, left);
        if (status == SimplexStatus::kOptimal) {
          ok = effective_bound(relaxation_objective()) < cutoff();
        } else if (status != SimplexStatus::kInfeasible) {
          abort = true;  // A limit or numerical trouble: stop the dive.
        }
      }
      if (!ok) undo_bounds(mark);
    }
    if (!ok) break;
    basis = simplex_->status();
    x = x_;
  }
  undo_bounds(0);
  logging_undo_ = false;
  outcome_.heuristic_lp_iterations += outcome_.lp_iterations - start;
}

// Objective feasibility pump (Fischetti, Glover and Lodi; Achterberg and Berthold): alternate
// between rounding the LP point and the LP point closest to the rounding in the L1 distance over
// integer columns at a bound, blended with the original objective whose weight decays.
// General integers strictly inside their bounds do not enter the distance.
void BranchAndBound::feasibility_pump(const std::vector<double>& x0,
                                      const std::vector<VarStatus>& basis0, long long budget) {
  const long long start = outcome_.lp_iterations;
  const std::vector<double> saved_cost(lp_.cost.begin(), lp_.cost.begin() + n_);
  double norm = 0.0;
  for (Index j = 0; j < n_; ++j) norm += cost_[j] * cost_[j];
  norm = std::sqrt(norm);
  const double objective_scale =
      norm > 0.0 ? std::sqrt(static_cast<double>(integers_.size())) / norm : 0.0;

  std::vector<double> x = x0;
  std::vector<VarStatus> basis = basis0;
  std::vector<double> rounded(static_cast<std::size_t>(n_), 0.0);
  std::vector<double> previous;
  std::vector<std::size_t> history;
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const auto flip = [&](Index j) {
    const double r = rounded[j];
    rounded[j] = std::clamp(x[j] > r ? r + 1.0 : r - 1.0, lower_[j], upper_[j]);
  };
  bool found = false;
  double alpha = 1.0;
  for (int iteration = 0; iteration < kPumpMaxIterations; ++iteration) {
    for (const Index j : integers_) rounded[j] = std::clamp(std::round(x[j]), lower_[j], upper_[j]);
    bool same = !previous.empty();
    for (std::size_t k = 0; same && k < integers_.size(); ++k) {
      same = rounded[integers_[k]] == previous[k];
    }
    if (same) {
      // A 1-cycle: flip the columns furthest from their rounding.
      std::vector<std::pair<double, Index>> far;
      for (const Index j : integers_) {
        const double d = std::fabs(x[j] - rounded[j]);
        if (d > 0.0) far.emplace_back(-d, j);
      }
      const int flips = std::uniform_int_distribution<int>(kPumpMinFlips, kPumpMaxFlips)(rng_);
      const std::size_t count = std::min(far.size(), static_cast<std::size_t>(flips));
      std::partial_sort(far.begin(), far.begin() + static_cast<std::ptrdiff_t>(count), far.end());
      for (std::size_t k = 0; k < count; ++k) flip(far[k].second);
    }
    std::size_t hash = 0;
    for (const Index j : integers_) {
      hash = hash * 1000003u ^ std::hash<double>{}(rounded[j]);
    }
    if (std::find(history.begin(), history.end(), hash) != history.end()) {
      // A longer cycle: perturb at random, more strongly near 0.5.
      for (const Index j : integers_) {
        const double t = std::max(0.0, unit(rng_) - 0.3);
        if (std::fabs(x[j] - rounded[j]) + t > 0.5) flip(j);
      }
    }
    history.push_back(hash);
    if (history.size() > kPumpHistory) history.erase(history.begin());
    previous.clear();
    for (const Index j : integers_) previous.push_back(rounded[j]);

    for (Index j = 0; j < n_; ++j) {
      double c = alpha * objective_scale * cost_[j];
      if (model_.col_type[j] == VarType::kInteger && lower_[j] < upper_[j]) {
        if (rounded[j] == lower_[j]) {
          c += 1.0;
        } else if (rounded[j] == upper_[j]) {
          c -= 1.0;
        }
      }
      lp_.cost[j] = c * scaling_.col[j];
    }
    const long long left = budget - (outcome_.lp_iterations - start);
    if (left <= 0 || time_up()) break;
    if (solve_relaxation(&basis, left) != SimplexStatus::kOptimal) break;
    basis = simplex_->status();
    x = x_;
    if (std::none_of(integers_.begin(), integers_.end(),
                     [&](Index j) { return is_fractional(x[j]); })) {
      found = try_solution(x);
      if (found) ++outcome_.heuristic_solutions;
      break;
    }
    alpha *= kPumpAlphaDecay;
  }
  for (Index j = 0; j < n_; ++j) lp_.cost[j] = saved_cost[j];
  outcome_.heuristic_lp_iterations += outcome_.lp_iterations - start;
  log_.log(2, "mip: feasibility pump %s after %lld LP iterations", found ? "succeeded" : "failed",
           outcome_.lp_iterations - start);
  // The last point's rounding with the continuous columns re-optimized.
  if (!found && !time_up()) round_and_solve(x, basis);
}

void BranchAndBound::rens(const std::vector<double>& x) {
  std::vector<double> lower = lower_;
  std::vector<double> upper = upper_;
  std::size_t fixed = 0;
  for (const Index j : integers_) {
    if (!is_fractional(x[j])) {
      lower[j] = upper[j] = std::clamp(std::round(x[j]), lower_[j], upper_[j]);
      ++fixed;
    } else {
      lower[j] = std::max(lower_[j], std::floor(x[j]));
      upper[j] = std::min(upper_[j], std::ceil(x[j]));
    }
  }
  if (static_cast<double>(fixed) < kRensMinFixed * static_cast<double>(integers_.size())) return;
  sub_mip(lower, upper, "RENS");
}

void BranchAndBound::rins(const std::vector<double>& x) {
  std::vector<double> lower = root_lower_;
  std::vector<double> upper = root_upper_;
  std::size_t fixed = 0;
  for (const Index j : integers_) {
    if (std::fabs(x[j] - incumbent_[j]) <= options_.integrality_tol) {
      lower[j] = upper[j] = incumbent_[j];
      ++fixed;
    }
  }
  if (static_cast<double>(fixed) < kRinsMinFixed * static_cast<double>(integers_.size())) return;
  sub_mip(lower, upper, "RINS");
}

// Solves the model (with the root cuts, which are valid everywhere) under the given column
// bounds by a nested, presolved search with small limits and hands its solution to try_solution.
void BranchAndBound::sub_mip(const std::vector<double>& lower, const std::vector<double>& upper,
                             const char* name) {
  const double remaining = remaining_time();
  const double seconds = std::min(kSubMipMaxSeconds, kSubMipTimeFraction * remaining);
  if (!(seconds > 0.0)) return;
  const Timer timer;
  struct Account {  // Adds this sub-MIP's time to the total however it returns.
    const Timer& timer;
    double& total;
    ~Account() { total += timer.seconds(); }
  } account{timer, sub_mip_seconds_};
  Model sub = model_;
  sub.col_lower = lower;
  sub.col_upper = upper;
  PresolveOptions presolve_options;
  presolve_options.mip = true;
  presolve_options.integrality_tol = options_.integrality_tol;
  Presolve presolve(sub, presolve_options);
  if (presolve.run() == PresolveStatus::kInfeasible) return;
  const Model& reduced = presolve.reduced();

  MipOutcome out;
  if (reduced.num_cols() == 0) {
    out.x.clear();
  } else {
    MipOptions options = options_;
    options.time_limit = seconds - timer.seconds();
    options.node_limit = kSubMipNodes;
    options.sub_mip_heuristics = false;
    options.threads = 1;
    options.debug_solution.clear();
    options.objective_cutoff.reset();
    if (!incumbent_.empty()) options.objective_cutoff = sense_ * cutoff();
    BranchAndBound search(reduced, options, lp_log_);
    out = search.solve();
    if (out.x.empty()) {
      log_.log(2, "mip: %s found no solution (%s, %lld nodes, %.2f s)", name,
               to_string(out.status), out.nodes, timer.seconds());
      return;
    }
  }
  std::vector<double> x;
  std::vector<double> unused_y;
  presolve.postsolve(out.x, {}, x, unused_y);
  const bool improved = try_solution(std::move(x));
  if (improved) {
    ++outcome_.heuristic_solutions;
    jump_incumbent_value_ = incumbent_value_;
  }
  log_.log(2, "mip: %s %s after %lld nodes, %.2f s", name,
           improved ? "improved the incumbent" : "found nothing better", out.nodes,
           timer.seconds());
}

}  // namespace samaya
