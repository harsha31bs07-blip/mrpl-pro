// Root reductions of the branch-and-bound that need the integrality of the columns: coefficient
// tightening on rows with binaries, and probing on binaries (fixings and implied bounds). Both
// keep every integer-feasible point of the model, so solutions are unaffected; the verifier
// still checks the final solution against the original model.

#include <algorithm>
#include <cmath>
#include <utility>

#include "mip/branch_and_bound.hpp"

namespace samaya {

namespace {

// Probing stops after this many probes or this share of the time limit (at most the cap).
constexpr int kMaxProbes = 5000;
constexpr double kProbingTimeFraction = 0.05;
constexpr double kProbingMaxSeconds = 10.0;
// Coefficients change only by more than this (relative), so round-off never rewrites a row.
constexpr double kMinCoefficientChange = 1e-9;

bool is_binary(double lower, double upper) { return lower == 0.0 && upper == 1.0; }

}  // namespace

// For a one-sided row sum a_j x_j <= b with maximum activity M > b and a binary x_k:
// - a_k > 0 and M - a_k < b: at x_k = 0 the row is redundant, so a_k and b both shrink by
//   d = b - (M - a_k); the row is unchanged at x_k = 1 and still redundant at x_k = 0.
// - a_k < 0 and M + a_k < b: at x_k = 1 the row is redundant, so a_k grows by d = b - (M + a_k)
//   toward zero; the row is unchanged at x_k = 0.
// A >= row is handled as its negation. Returns the number of changed coefficients.
int BranchAndBound::tighten_coefficients() {
  const auto rstart = At_.col_start();
  const auto rindex = At_.row_index();
  const auto rval = At_.values();
  std::vector<Triplet> t;
  t.reserve(static_cast<std::size_t>(model_.A.nnz()));
  int changed = 0;
  for (Index i = 0; i < m_; ++i) {
    const double rl = model_.row_lower[i];
    const double ru = model_.row_upper[i];
    const bool upper_row = ru < kInf && rl == -kInf;
    const bool lower_row = rl > -kInf && ru == kInf;
    const NnzIndex begin = rstart[i];
    const NnzIndex end = rstart[i + 1];
    if (!upper_row && !lower_row) {
      for (NnzIndex p = begin; p < end; ++p) t.push_back({i, rindex[p], rval[p]});
      continue;
    }
    // Work with sigma * row <= b.
    const double sigma = upper_row ? 1.0 : -1.0;
    double b = sigma * (upper_row ? ru : rl);
    double max_activity = 0.0;
    for (NnzIndex p = begin; p < end; ++p) {
      const double a = sigma * rval[p];
      max_activity += a * (a > 0.0 ? upper_[rindex[p]] : lower_[rindex[p]]);
    }
    std::vector<double> coefficient(static_cast<std::size_t>(end - begin));
    for (NnzIndex p = begin; p < end; ++p) coefficient[p - begin] = sigma * rval[p];
    if (std::isfinite(max_activity) && max_activity > b) {
      for (NnzIndex p = begin; p < end; ++p) {
        const Index j = rindex[p];
        if (model_.col_type[j] != VarType::kInteger || !is_binary(lower_[j], upper_[j])) continue;
        double& a = coefficient[p - begin];
        const double scale = kMinCoefficientChange * (1.0 + std::fabs(a) + std::fabs(b));
        if (a > 0.0 && max_activity - a < b - scale) {
          const double d = b - (max_activity - a);
          a -= d;
          b -= d;
          max_activity -= d;
          ++changed;
        } else if (a < 0.0 && max_activity + a < b - scale) {
          a += b - (max_activity + a);
          ++changed;
        }
      }
    }
    for (NnzIndex p = begin; p < end; ++p) {
      const double a = sigma * coefficient[p - begin];
      if (a != 0.0) t.push_back({i, rindex[p], a});
    }
    if (upper_row) {
      model_.row_upper[i] = b;
    } else {
      model_.row_lower[i] = -b;
    }
  }
  if (changed > 0) {
    model_.A = SparseMatrix::from_triplets(m_, n_, std::move(t));
    build_relaxation();
  }
  return changed;
}

// Probing: each binary is fixed to 0 and to 1 in turn and the fixing propagated. A side that is
// infeasible fixes the binary to the other; bounds implied on both sides hold at the root.
bool BranchAndBound::probe() {
  double seconds = kProbingMaxSeconds;
  if (std::isfinite(options_.time_limit)) {
    seconds = std::min(seconds, kProbingTimeFraction * options_.time_limit);
  }
  const Timer timer;
  std::vector<Index> candidates;
  const auto cstart = model_.A.col_start();
  for (const Index j : integers_) {
    if (is_binary(lower_[j], upper_[j])) candidates.push_back(j);
  }
  // Longest columns first: their fixings reach the most rows.
  std::stable_sort(candidates.begin(), candidates.end(), [&](Index a, Index b) {
    return cstart[a + 1] - cstart[a] > cstart[b + 1] - cstart[b];
  });

  std::vector<double> side_lower[2];
  std::vector<double> side_upper[2];
  std::vector<Index> side_cols[2];
  std::vector<char> seen(static_cast<std::size_t>(n_), 0);
  std::vector<Index> position(static_cast<std::size_t>(n_), -1);
  int probes = 0;
  undo_log_.clear();
  for (const Index j : candidates) {
    if (probes >= kMaxProbes || timer.seconds() > seconds || time_up()) break;
    if (!is_binary(lower_[j], upper_[j])) continue;  // Fixed by an earlier probe.
    ++probes;
    bool feasible[2] = {false, false};
    for (int v = 0; v < 2; ++v) {
      logging_undo_ = true;
      set_bound(j, v, v);
      feasible[v] = propagate({j}, nullptr);
      side_cols[v].clear();
      side_lower[v].clear();
      side_upper[v].clear();
      if (feasible[v]) {
        // The undo log names every changed column; its current bounds are this side's.
        for (const BoundChange& c : undo_log_) {
          if (c.col == j || seen[c.col]) continue;
          seen[c.col] = 1;
          side_cols[v].push_back(c.col);
          side_lower[v].push_back(lower_[c.col]);
          side_upper[v].push_back(upper_[c.col]);
        }
        for (const Index k : side_cols[v]) seen[k] = 0;
      }
      undo_bounds(0);
      logging_undo_ = false;
    }
    if (!feasible[0] && !feasible[1]) return false;
    if (!feasible[0] || !feasible[1]) {
      const double v = feasible[0] ? 0.0 : 1.0;
      set_bound(j, v, v);
      ++outcome_.probing_fixed;
      if (!propagate({j}, nullptr)) return false;
      continue;
    }
    // Bounds implied by both sides: the union of the two sides' boxes.
    std::vector<std::pair<Index, std::pair<double, double>>> implied;
    for (std::size_t a = 0; a < side_cols[0].size(); ++a) {
      position[side_cols[0][a]] = static_cast<Index>(a);
    }
    for (std::size_t b = 0; b < side_cols[1].size(); ++b) {
      const Index k = side_cols[1][b];
      if (position[k] < 0) continue;
      const auto a = static_cast<std::size_t>(position[k]);
      const double lo = std::min(side_lower[0][a], side_lower[1][b]);
      const double up = std::max(side_upper[0][a], side_upper[1][b]);
      if (lo > lower_[k] || up < upper_[k]) {
        implied.push_back({k, {std::max(lo, lower_[k]), std::min(up, upper_[k])}});
      }
    }
    for (const Index k : side_cols[0]) position[k] = -1;
    std::vector<Index> changed;
    for (const auto& [k, box] : implied) {
      set_bound(k, box.first, box.second);
      changed.push_back(k);
      ++outcome_.probing_tightened;
    }
    if (!changed.empty() && !propagate(changed, nullptr)) return false;
  }
  log_.log(1, "MIP probing: %d probes, %d columns fixed, %d bounds tightened, %.2f s", probes,
           outcome_.probing_fixed, outcome_.probing_tightened, timer.seconds());
  return true;
}

}  // namespace samaya
