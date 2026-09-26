#include "mip/feasibility_jump.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <utility>

#include "core/log.hpp"

namespace samaya {

namespace {

// A row counts as satisfied within this relative tolerance (tighter than the solver's acceptance
// tolerance, so a point found here passes try_solution).
constexpr double kRowTol = 1e-9;
// Columns sampled from a long violated row per step.
constexpr std::size_t kMaxCandidates = 64;
constexpr double kWeightIncrement = 1.0;
// Time is checked every this many steps.
constexpr long long kClockInterval = 256;

}  // namespace

FeasibilityJumpResult feasibility_jump(const Model& model, Index rows,
                                       const std::vector<double>& lower,
                                       const std::vector<double>& upper, double seconds,
                                       long long max_work, std::uint64_t seed,
                                       const std::vector<double>* start) {
  FeasibilityJumpResult result;
  const Index n = model.num_cols();
  const SparseMatrix at = model.A.transpose();
  const auto cstart = model.A.col_start();
  const auto cindex = model.A.row_index();
  const auto cvalue = model.A.values();
  const auto rstart = at.col_start();
  const auto rindex = at.row_index();
  const Timer timer;
  std::mt19937_64 rng(seed);

  std::vector<double> x(static_cast<std::size_t>(n));
  const bool has_start = start != nullptr && start->size() == static_cast<std::size_t>(n);
  for (Index j = 0; j < n; ++j) {
    if (has_start && std::isfinite((*start)[j])) {
      double v = (*start)[j];
      if (model.col_type[j] == VarType::kInteger) v = std::round(v);
      x[j] = std::clamp(v, lower[j], upper[j]);
      continue;
    }
    double v = std::clamp(0.0, lower[j], upper[j]);
    if (model.col_type[j] == VarType::kInteger) v = std::ceil(v - 1e-9);
    x[j] = std::clamp(v, lower[j], upper[j]);
  }
  std::vector<double> activity(static_cast<std::size_t>(rows), 0.0);
  for (Index j = 0; j < n; ++j) {
    for (NnzIndex p = cstart[j]; p < cstart[j + 1]; ++p) {
      if (cindex[p] < rows) activity[cindex[p]] += cvalue[p] * x[j];
    }
  }
  std::vector<double> weight(static_cast<std::size_t>(rows), 1.0);
  const auto violation = [&](Index i, double act) {
    const double rl = model.row_lower[i];
    const double ru = model.row_upper[i];
    if (act < rl - kRowTol * (1.0 + std::fabs(rl))) return rl - act;
    if (act > ru + kRowTol * (1.0 + std::fabs(ru))) return act - ru;
    return 0.0;
  };
  std::vector<Index> violated;
  std::vector<Index> position(static_cast<std::size_t>(rows), -1);
  const auto update_row = [&](Index i) {
    const bool bad = violation(i, activity[i]) > 0.0;
    if (bad && position[i] < 0) {
      position[i] = static_cast<Index>(violated.size());
      violated.push_back(i);
    } else if (!bad && position[i] >= 0) {
      const Index last = violated.back();
      violated[position[i]] = last;
      position[last] = position[i];
      violated.pop_back();
      position[i] = -1;
    }
  };
  for (Index i = 0; i < rows; ++i) update_row(i);

  // Weighted violation of column j's rows if x_j were t.
  const auto cost = [&](Index j, double t) {
    double c = 0.0;
    for (NnzIndex p = cstart[j]; p < cstart[j + 1]; ++p) {
      const Index i = cindex[p];
      if (i >= rows) continue;
      c += weight[i] * violation(i, activity[i] + cvalue[p] * (t - x[j]));
    }
    return c;
  };
  // The value of column j minimizing that cost: each row contributes w |a| times the distance of
  // t to the interval where it is satisfied, so the minimum is where the slope turns >= 0.
  std::vector<std::pair<double, double>> events;
  const auto jump = [&](Index j) {
    events.clear();
    double slope = 0.0;
    for (NnzIndex p = cstart[j]; p < cstart[j + 1]; ++p) {
      const Index i = cindex[p];
      if (i >= rows) continue;
      const double a = cvalue[p];
      const double r = activity[i] - a * x[j];
      const double w = weight[i] * std::fabs(a);
      double lo = a > 0.0 ? (model.row_lower[i] - r) / a : (model.row_upper[i] - r) / a;
      double hi = a > 0.0 ? (model.row_upper[i] - r) / a : (model.row_lower[i] - r) / a;
      if (std::isnan(lo)) lo = -kInf;
      if (std::isnan(hi)) hi = kInf;
      if (lo > -kInf) {
        slope -= w;
        events.push_back({lo, w});
      }
      if (hi < kInf) events.push_back({hi, w});
    }
    std::sort(events.begin(), events.end());
    double t = x[j];
    if (slope < 0.0) {
      for (std::size_t k = 0; k < events.size(); ++k) {
        slope += events[k].second;
        if (slope >= 0.0) {
          // Flat stretch up to the next event when the slope is exactly 0: stay close to x_j.
          const double left = events[k].first;
          const double right = slope == 0.0 && k + 1 < events.size() ? events[k + 1].first : left;
          t = std::clamp(x[j], left, right);
          break;
        }
      }
    } else if (!events.empty()) {
      t = std::min(x[j], events.front().first);  // Already at a minimum when x_j is left of it.
    }
    t = std::clamp(t, lower[j], upper[j]);
    if (model.col_type[j] == VarType::kInteger) {
      const double down = std::clamp(std::floor(t + 1e-9), lower[j], upper[j]);
      const double up = std::clamp(std::ceil(t - 1e-9), lower[j], upper[j]);
      t = cost(j, down) <= cost(j, up) ? down : up;
    }
    return t;
  };

  std::vector<Index> candidates;
  while (!violated.empty()) {
    if (result.work >= max_work) break;
    if (result.steps % kClockInterval == 0 && timer.seconds() > seconds) break;
    ++result.steps;
    const Index i = violated[std::uniform_int_distribution<std::size_t>(0, violated.size() - 1)(rng)];
    candidates.clear();
    for (NnzIndex p = rstart[i]; p < rstart[i + 1]; ++p) {
      const Index j = rindex[p];
      if (lower[j] < upper[j]) candidates.push_back(j);
    }
    if (candidates.size() > kMaxCandidates) {
      std::shuffle(candidates.begin(), candidates.end(), rng);
      candidates.resize(kMaxCandidates);
    }
    Index best = -1;
    double best_value = 0.0;
    double best_score = -kInf;
    for (const Index j : candidates) {
      const double t = jump(j);
      result.work += 3 * (cstart[j + 1] - cstart[j]);
      if (t == x[j]) continue;
      const double score = cost(j, x[j]) - cost(j, t);
      if (score > best_score) {
        best_score = score;
        best = j;
        best_value = t;
      }
    }
    if (!(best_score > 1e-12)) {
      // Local minimum: the violated rows weigh more from now on.
      for (const Index v : violated) weight[v] += kWeightIncrement;
      result.work += static_cast<long long>(violated.size());
      if (best < 0) continue;
    }
    const double delta = best_value - x[best];
    x[best] = best_value;
    for (NnzIndex p = cstart[best]; p < cstart[best + 1]; ++p) {
      const Index r = cindex[p];
      if (r >= rows) continue;
      activity[r] += cvalue[p] * delta;
      update_row(r);
    }
    result.work += cstart[best + 1] - cstart[best];
    if (violated.empty()) {
      // Recompute the activities exactly before trusting the incremental ones.
      std::fill(activity.begin(), activity.end(), 0.0);
      for (Index j = 0; j < n; ++j) {
        for (NnzIndex p = cstart[j]; p < cstart[j + 1]; ++p) {
          if (cindex[p] < rows) activity[cindex[p]] += cvalue[p] * x[j];
        }
      }
      for (Index r = 0; r < rows; ++r) update_row(r);
    }
  }
  if (violated.empty()) {
    result.found = true;
    result.x = std::move(x);
  }
  return result;
}

}  // namespace samaya
