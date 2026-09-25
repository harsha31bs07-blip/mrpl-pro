#include "mip/cuts.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace samaya {

namespace {

constexpr double kMinFractionality = 0.01;  // Gomory: basic value at least this far from integer.
constexpr double kTinyCoefficient = 1e-12;
constexpr double kMaxDynamism = 1e6;    // Largest ratio of cut coefficient magnitudes.
constexpr double kRhsRelaxation = 1e-9;  // Relative safety margin on the right-hand side.
constexpr double kMinViolation = 1e-6;
constexpr int kMaxMirDivisors = 8;
constexpr NnzIndex kMaxRowLength = 1000;  // MIR and cover rows longer than this are skipped.

bool is_integer_col(const Model& model, Index j) {
  return model.col_type[j] == VarType::kInteger;
}

// Accumulates a cut as dense column coefficients plus a right-hand side.
class CutBuilder {
 public:
  explicit CutBuilder(const CutContext& ctx)
      : ctx_(ctx),
        dense_(static_cast<std::size_t>(ctx.model.num_cols()), 0.0),
        mark_(static_cast<std::size_t>(ctx.model.num_cols()), 0) {}

  void add(Index j, double v) {
    if (!mark_[j]) {
      mark_[j] = 1;
      touched_.push_back(j);
    }
    dense_[j] += v;
  }
  // Adds v times row i's activity.
  void add_row(Index i, double v) {
    const auto start = ctx_.At.col_start();
    const auto index = ctx_.At.row_index();
    const auto value = ctx_.At.values();
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) add(index[p], v * value[p]);
  }
  // Returns  sum dense_j x_j  (sign) rhs  as a >= cut: sign +1 keeps it, -1 negates a <= cut.
  Cut take(double rhs, double sign) {
    Cut cut;
    for (const Index j : touched_) {
      if (dense_[j] != 0.0) {
        cut.index.push_back(j);
        cut.value.push_back(sign * dense_[j]);
      }
      dense_[j] = 0.0;
      mark_[j] = 0;
    }
    touched_.clear();
    cut.lower = sign * rhs;
    return cut;
  }
  void clear() { take(0.0, 1.0); }

 private:
  const CutContext& ctx_;
  std::vector<double> dense_;
  std::vector<char> mark_;
  std::vector<Index> touched_;
};

double violation_of(const CutContext& ctx, const Cut& cut) {
  double activity = 0.0;
  for (std::size_t k = 0; k < cut.index.size(); ++k) activity += cut.value[k] * ctx.x[cut.index[k]];
  return cut.lower - activity;
}

double norm_of(const Cut& cut) {
  double s = 0.0;
  for (const double v : cut.value) s += v * v;
  return std::sqrt(s);
}

}  // namespace

bool finalize_cut(const CutContext& ctx, Cut& cut) {
  if (cut.index.empty() || !std::isfinite(cut.lower)) return false;
  double max_abs = 0.0;
  for (const double v : cut.value) {
    if (!std::isfinite(v)) return false;
    max_abs = std::max(max_abs, std::fabs(v));
  }
  if (max_abs == 0.0) return false;
  // Drop tiny coefficients: pi_j x_j <= max(pi_j l_j, pi_j u_j), so the rest must cover the
  // right-hand side minus that amount.
  std::size_t kept = 0;
  double min_abs = kInf;
  for (std::size_t k = 0; k < cut.index.size(); ++k) {
    const Index j = cut.index[k];
    const double v = cut.value[k];
    if (std::fabs(v) < kTinyCoefficient * max_abs) {
      const double largest = std::max(v * ctx.lower[j], v * ctx.upper[j]);
      if (!std::isfinite(largest)) return false;
      cut.lower -= largest;
      continue;
    }
    min_abs = std::min(min_abs, std::fabs(v));
    cut.index[kept] = j;
    cut.value[kept] = v;
    ++kept;
  }
  cut.index.resize(kept);
  cut.value.resize(kept);
  if (kept == 0 || max_abs / min_abs > kMaxDynamism) return false;
  cut.lower -= kRhsRelaxation * std::max(1.0, std::fabs(cut.lower));
  const double violation = violation_of(ctx, cut);
  if (violation <= kMinViolation * (1.0 + std::fabs(cut.lower))) return false;
  cut.efficacy = violation / norm_of(cut);
  return true;
}

bool gomory_mixed_integer_cut(const CutContext& ctx, Index k, const std::vector<double>& row,
                              const std::vector<VarStatus>& status, Cut& cut) {
  const Model& model = ctx.model;
  const Index n = model.num_cols();
  const Index m = model.num_rows();
  // Tableau row: x_k + sum_N a_j v_j = 0. With v_j = bound_j +- t_j (t_j >= 0) it becomes
  // x_k + sum a'_j t_j = b'.
  struct Term {
    Index var;
    double a;      // a'_j
    double bound;  // The bound t_j is measured from.
    bool at_upper;
  };
  std::vector<Term> terms;
  double b = 0.0;
  double row_max = 0.0;
  for (Index j = 0; j < n + m; ++j) {
    if (j != k) row_max = std::max(row_max, std::fabs(row[j]));
  }
  for (Index j = 0; j < n + m; ++j) {
    if (j == k || status[j] == VarStatus::kBasic) continue;
    const double a = row[j];
    if (std::fabs(a) <= kTinyCoefficient * std::max(1.0, row_max)) continue;
    const double lo = j < n ? ctx.lower[j] : model.row_lower[j - n];
    const double up = j < n ? ctx.upper[j] : model.row_upper[j - n];
    if (status[j] == VarStatus::kAtZero) return false;  // Free nonbasic: no t_j >= 0 form.
    const bool at_upper = status[j] == VarStatus::kAtUpper;
    const double bound = at_upper ? up : lo;
    if (!std::isfinite(bound)) return false;
    b -= a * bound;
    if (lo == up) continue;  // Fixed: t_j = 0.
    terms.push_back({j, at_upper ? -a : a, bound, at_upper});
  }
  const double xk = ctx.x[k];
  if (std::fabs(b - xk) > 1e-6 * (1.0 + std::fabs(xk))) return false;  // Inaccurate row.
  const double f0 = b - std::floor(b);
  if (f0 < kMinFractionality || f0 > 1.0 - kMinFractionality) return false;

  // GMI:  sum_int min(f_j/f0, (1-f_j)/(1-f0)) t_j
  //     + sum_cont (a_j > 0 ? a_j/f0 : -a_j/(1-f0)) t_j >= 1.
  CutBuilder builder(ctx);
  double rhs = 1.0;
  for (const Term& t : terms) {
    const bool integral = t.var < n && is_integer_col(model, t.var) &&
                          t.bound == std::round(t.bound);
    double coef;
    if (integral) {
      const double fj = t.a - std::floor(t.a);
      coef = fj <= f0 ? fj / f0 : (1.0 - fj) / (1.0 - f0);
    } else {
      coef = t.a >= 0.0 ? t.a / f0 : -t.a / (1.0 - f0);
    }
    if (coef == 0.0) continue;
    // t = v - bound (at lower) or bound - v (at upper).
    const double v_coef = t.at_upper ? -coef : coef;
    rhs += v_coef * t.bound;
    if (t.var < n) {
      builder.add(t.var, v_coef);
    } else {
      builder.add_row(t.var - n, v_coef);
    }
  }
  cut = builder.take(rhs, 1.0);
  return true;
}

void separate_mir(const CutContext& ctx, std::vector<Cut>& cuts) {
  const Model& model = ctx.model;
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  CutBuilder builder(ctx);
  struct IntTerm {
    Index col;
    double c;  // Coefficient of t_j >= 0.
    double t;  // t_j at the LP point.
    double bound;
    bool at_upper;
  };
  struct ContTerm {
    Index col;
    double h;
    double bound;
    bool at_upper;
  };
  std::vector<IntTerm> ints;
  std::vector<ContTerm> conts;
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] > kMaxRowLength || start[i + 1] == start[i]) continue;
    for (int side = 0; side < 2; ++side) {
      // side 0: a x <= ru;  side 1: -a x <= -rl.
      const double bound = side == 0 ? model.row_upper[i] : model.row_lower[i];
      if (!std::isfinite(bound)) continue;
      const double sign = side == 0 ? 1.0 : -1.0;
      double beta = sign * bound;
      ints.clear();
      conts.clear();
      bool ok = true;
      for (NnzIndex p = start[i]; p < start[i + 1] && ok; ++p) {
        const Index j = index[p];
        const double c = sign * value[p];
        const double lo = ctx.lower[j];
        const double up = ctx.upper[j];
        const double xj = ctx.x[j];
        const bool use_lower =
            std::isfinite(lo) && (!std::isfinite(up) || xj - lo <= up - xj);
        if (!use_lower && !std::isfinite(up)) {
          ok = false;
          break;
        }
        const double b = use_lower ? lo : up;
        beta -= c * b;
        if (lo == up) continue;
        const double coef = use_lower ? c : -c;
        if (is_integer_col(model, j)) {
          ints.push_back({j, coef, use_lower ? xj - lo : up - xj, b, !use_lower});
        } else {
          conts.push_back({j, coef, b, !use_lower});
        }
      }
      if (!ok || ints.empty()) continue;

      // Divisors: coefficients of integer columns away from their bound.
      std::vector<double> divisors;
      for (const IntTerm& t : ints) {
        if (t.t <= 1e-6) continue;
        const double d = std::fabs(t.c);
        if (d < 1e-6) continue;
        if (std::none_of(divisors.begin(), divisors.end(),
                         [&](double e) { return std::fabs(e - d) <= 1e-9 * d; })) {
          divisors.push_back(d);
        }
        if (static_cast<int>(divisors.size()) >= kMaxMirDivisors) break;
      }
      const auto try_divisor = [&](double delta, Cut& out) {
        const double scaled = beta / delta;
        const double f = scaled - std::floor(scaled);
        if (f < 0.05 || f > 0.95) return false;
        // sum F(c_j/delta) t_j + sum_{h<0} h_j s_j / (delta (1-f)) <= floor(beta/delta).
        double rhs = std::floor(scaled);
        for (const IntTerm& t : ints) {
          const double a = t.c / delta;
          const double fa = a - std::floor(a);
          const double coef = std::floor(a) + std::max(0.0, fa - f) / (1.0 - f);
          if (coef == 0.0) continue;
          const double x_coef = t.at_upper ? -coef : coef;
          rhs += x_coef * t.bound;  // coef * t with t = x - l or u - x.
          builder.add(t.col, x_coef);
        }
        for (const ContTerm& s : conts) {
          if (s.h >= 0.0) continue;
          const double coef = s.h / (delta * (1.0 - f));
          const double x_coef = s.at_upper ? -coef : coef;
          rhs += x_coef * s.bound;
          builder.add(s.col, x_coef);
        }
        out = builder.take(rhs, -1.0);
        return finalize_cut(ctx, out);
      };
      Cut best;
      double best_delta = 0.0;
      for (const double delta : divisors) {
        Cut c;
        if (try_divisor(delta, c) && c.efficacy > best.efficacy) {
          best = std::move(c);
          best_delta = delta;
        }
      }
      if (best_delta > 0.0) {
        for (const double factor : {0.5, 0.25, 0.125}) {
          Cut c;
          if (try_divisor(best_delta * factor, c) && c.efficacy > best.efficacy) {
            best = std::move(c);
          }
        }
        cuts.push_back(std::move(best));
      }
    }
  }
}

void separate_knapsack_covers(const CutContext& ctx, std::vector<Cut>& cuts) {
  const Model& model = ctx.model;
  const auto start = ctx.At.col_start();
  const auto index = ctx.At.row_index();
  const auto value = ctx.At.values();
  CutBuilder builder(ctx);
  struct Item {
    Index col;
    double w;
    double y;  // LP value of the (possibly complemented) binary.
    bool complemented;
  };
  std::vector<Item> items;
  for (Index i = 0; i < ctx.original_rows; ++i) {
    if (start[i + 1] - start[i] > kMaxRowLength || start[i + 1] - start[i] < 2) continue;
    for (int side = 0; side < 2; ++side) {
      const double bound = side == 0 ? model.row_upper[i] : model.row_lower[i];
      if (!std::isfinite(bound)) continue;
      const double sign = side == 0 ? 1.0 : -1.0;
      double capacity = sign * bound;
      items.clear();
      bool ok = true;
      for (NnzIndex p = start[i]; p < start[i + 1]; ++p) {
        const Index j = index[p];
        const double w = sign * value[p];
        const bool binary = is_integer_col(model, j) && ctx.lower[j] == 0.0 && ctx.upper[j] == 1.0;
        if (!binary) {
          // Relax to the smallest contribution.
          const double least = w > 0.0 ? w * ctx.lower[j] : w * ctx.upper[j];
          if (!std::isfinite(least)) {
            ok = false;
            break;
          }
          capacity -= least;
          continue;
        }
        if (w > 0.0) {
          items.push_back({j, w, ctx.x[j], false});
        } else {
          capacity -= w;  // w x = w - w (1 - x).
          items.push_back({j, -w, 1.0 - ctx.x[j], true});
        }
      }
      if (!ok || items.size() < 2 || capacity < 0.0) continue;
      const double total = std::accumulate(items.begin(), items.end(), 0.0,
                                           [](double s, const Item& it) { return s + it.w; });
      const double eps = 1e-9 * (1.0 + std::fabs(capacity));
      if (total <= capacity + eps) continue;

      // Greedy cover: prefer items with large LP value per unit weight.
      std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return (1.0 - a.y) * b.w < (1.0 - b.y) * a.w;
      });
      std::vector<std::size_t> cover;
      double weight = 0.0;
      for (std::size_t t = 0; t < items.size() && weight <= capacity + eps; ++t) {
        cover.push_back(t);
        weight += items[t].w;
      }
      if (weight <= capacity + eps) continue;
      // Make it minimal, removing items with the smallest LP value first.
      std::vector<std::size_t> order = cover;
      std::sort(order.begin(), order.end(),
                [&](std::size_t a, std::size_t b) { return items[a].y < items[b].y; });
      for (const std::size_t t : order) {
        if (weight - items[t].w > capacity + eps) {
          weight -= items[t].w;
          cover.erase(std::find(cover.begin(), cover.end(), t));
        }
      }
      double lhs = 0.0;
      double max_w = 0.0;
      for (const std::size_t t : cover) {
        lhs += items[t].y;
        max_w = std::max(max_w, items[t].w);
      }
      const double size = static_cast<double>(cover.size());
      if (lhs <= size - 1.0 + 1e-6) continue;
      // Extended cover: every item at least as heavy as the heaviest cover item joins.
      std::vector<char> in_cover(items.size(), 0);
      for (const std::size_t t : cover) in_cover[t] = 1;
      double rhs = size - 1.0;  // sum_E y <= |C| - 1, with y = x or 1 - x.
      for (std::size_t t = 0; t < items.size(); ++t) {
        if (!in_cover[t] && items[t].w < max_w) continue;
        if (items[t].complemented) {
          builder.add(items[t].col, -1.0);
          rhs -= 1.0;
        } else {
          builder.add(items[t].col, 1.0);
        }
      }
      Cut cut = builder.take(rhs, -1.0);
      if (finalize_cut(ctx, cut)) cuts.push_back(std::move(cut));
    }
  }
}

}  // namespace samaya
