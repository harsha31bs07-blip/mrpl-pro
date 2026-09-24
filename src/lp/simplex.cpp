#include "lp/simplex.hpp"

#include <algorithm>
#include <cmath>

namespace samaya {

const char* to_string(SimplexStatus status) {
  switch (status) {
    case SimplexStatus::kOptimal: return "optimal";
    case SimplexStatus::kInfeasible: return "infeasible";
    case SimplexStatus::kUnbounded: return "unbounded";
    case SimplexStatus::kIterationLimit: return "iteration_limit";
    case SimplexStatus::kTimeLimit: return "time_limit";
    case SimplexStatus::kNumericalError: return "numerical_error";
  }
  return "unknown";
}

namespace {

constexpr double kMinDseWeight = 1e-4;
constexpr double kPhase1FreeBound = 1000.0;
constexpr double kPerturbationBase = 5e-7;
// Relative disagreement between the pivot element computed from the column (ftran) and the row
// (btran) that triggers a refactorization.
constexpr double kPivotConsistencyTol = 1e-7;

BasisFactor::Options factor_options(const SimplexOptions& options) {
  BasisFactor::Options f;
  f.max_updates = options.refactor_interval;
  return f;
}

}  // namespace

Simplex::Simplex(const LpProblem& lp, const SimplexOptions& options, const Logger& log)
    : lp_(lp),
      options_(options),
      log_(log),
      m_(lp.m),
      n_(lp.n),
      nt_(lp.n + lp.m),
      factor_(lp.A, factor_options(options)),
      rng_(options.seed) {
  max_iterations_ = options.max_iterations >= 0
                        ? options.max_iterations
                        : std::max<long long>(100000, 100LL * static_cast<long long>(nt_));
  const auto m = static_cast<std::size_t>(m_);
  const auto nt = static_cast<std::size_t>(nt_);
  basic_.resize(m);
  position_.assign(nt, -1);
  status_.assign(nt, VarStatus::kAtLower);
  x_.assign(nt, 0.0);
  y_.assign(m, 0.0);
  d_.assign(nt, 0.0);
  dse_weight_.assign(m, 1.0);
  alpha_row_.assign(nt, 0.0);
  alpha_col_.assign(m, 0.0);
}

template <typename F>
void Simplex::for_column(Index j, F&& f) const {
  if (j < n_) {
    const auto start = lp_.A.col_start();
    const auto index = lp_.A.row_index();
    const auto value = lp_.A.values();
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) f(index[p], value[p]);
  } else {
    f(j - n_, -1.0);
  }
}

void Simplex::load_column(Index j, std::vector<double>& out) const {
  out.assign(static_cast<std::size_t>(m_), 0.0);
  for_column(j, [&](Index i, double a) { out[i] = a; });
}

void Simplex::set_value_from_status(Index j) {
  switch (status_[j]) {
    case VarStatus::kAtLower: x_[j] = lower_[j]; break;
    case VarStatus::kAtUpper: x_[j] = upper_[j]; break;
    case VarStatus::kAtZero: x_[j] = 0.0; break;
    case VarStatus::kBasic: break;
  }
}

bool Simplex::rebuild() {
  for (int attempt = 0; attempt < 4; ++attempt) {
    ++stats_.refactorizations;
    if (factor_.factorize(basic_) == 0) {
      compute_primal();
      compute_dual();
      return true;
    }
    // Replace dependent basis columns with the logicals of the rows left unpivoted.
    const auto& positions = factor_.singular_positions();
    const auto& rows = factor_.unpivoted_rows();
    log_.log(3, "simplex: basis rank deficiency %zu, repairing", positions.size());
    for (std::size_t k = 0; k < positions.size(); ++k) {
      const Index pos = positions[k];
      const Index entering = n_ + rows[k];
      const Index leaving = basic_[pos];
      if (position_[entering] != -1) return false;
      const bool has_lower = lower_[leaving] > -kInf;
      const bool has_upper = upper_[leaving] < kInf;
      if (has_lower && (!has_upper || std::fabs(x_[leaving] - lower_[leaving]) <=
                                          std::fabs(x_[leaving] - upper_[leaving]))) {
        status_[leaving] = VarStatus::kAtLower;
      } else if (has_upper) {
        status_[leaving] = VarStatus::kAtUpper;
      } else {
        status_[leaving] = VarStatus::kAtZero;
      }
      set_value_from_status(leaving);
      position_[leaving] = -1;
      basic_[pos] = entering;
      position_[entering] = pos;
      status_[entering] = VarStatus::kBasic;
      dse_weight_[pos] = 1.0;
    }
  }
  return false;
}

void Simplex::compute_primal() {
  std::vector<double> rhs(static_cast<std::size_t>(m_), 0.0);
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const double v = x_[j];
    if (v == 0.0) continue;
    for_column(j, [&](Index i, double a) { rhs[i] -= a * v; });
  }
  factor_.ftran(rhs);
  for (Index k = 0; k < m_; ++k) x_[basic_[k]] = rhs[k];
}

void Simplex::compute_dual() {
  std::vector<double> cb(static_cast<std::size_t>(m_));
  for (Index k = 0; k < m_; ++k) cb[k] = cost_[basic_[k]];
  factor_.btran(cb);
  y_ = std::move(cb);
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) {
      d_[j] = 0.0;
      continue;
    }
    double s = cost_[j];
    for_column(j, [&](Index i, double a) { s -= a * y_[i]; });
    d_[j] = s;
  }
}

void Simplex::compute_pivot_row() {
  std::fill(alpha_row_.begin(), alpha_row_.end(), 0.0);
  const auto start = lp_.At.col_start();
  const auto index = lp_.At.row_index();
  const auto value = lp_.At.values();
  for (Index i = 0; i < m_; ++i) {
    const double r = rho_[i];
    if (r == 0.0) continue;
    for (NnzIndex p = start[i]; p < start[i + 1]; ++p) alpha_row_[index[p]] += r * value[p];
    alpha_row_[n_ + i] = -r;
  }
}

void Simplex::place_nonbasic_dual_feasible() {
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (has_lower && has_upper) {
      status_[j] = (lower_[j] == upper_[j] || d_[j] >= 0.0) ? VarStatus::kAtLower
                                                             : VarStatus::kAtUpper;
    } else if (has_lower) {
      status_[j] = VarStatus::kAtLower;
    } else if (has_upper) {
      status_[j] = VarStatus::kAtUpper;
    } else {
      status_[j] = VarStatus::kAtZero;
    }
    set_value_from_status(j);
  }
}

bool Simplex::correct_dual_infeasibilities() {
  bool flipped = false;
  bool changed = false;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    const VarStatus st = status_[j];
    if (st == VarStatus::kBasic || lower_[j] == upper_[j]) continue;
    const double dj = d_[j];
    const bool infeasible = (st == VarStatus::kAtLower && dj < -tol) ||
                            (st == VarStatus::kAtUpper && dj > tol) ||
                            (st == VarStatus::kAtZero && std::fabs(dj) > tol);
    if (!infeasible) continue;
    changed = true;
    if (lower_[j] > -kInf && upper_[j] < kInf) {
      status_[j] = dj >= 0.0 ? VarStatus::kAtLower : VarStatus::kAtUpper;
      set_value_from_status(j);
      flipped = true;
    } else {
      // Shift the cost so that d_j = 0; shifts are removed before optimality is declared.
      ++stats_.cost_shifts;
      cost_[j] -= dj;
      d_[j] = 0.0;
    }
  }
  if (flipped) compute_primal();
  return changed;
}

void Simplex::perturb_costs() {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (Index j = 0; j < nt_; ++j) {
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (lower_[j] == upper_[j] || (!has_lower && !has_upper)) continue;
    const double xi = kPerturbationBase * (1.0 + std::fabs(cost_[j])) * (1.0 + unit(rng_));
    switch (status_[j]) {
      case VarStatus::kAtLower: cost_[j] += xi; break;
      case VarStatus::kAtUpper: cost_[j] -= xi; break;
      case VarStatus::kAtZero: break;
      case VarStatus::kBasic:
        if (has_lower && !has_upper) {
          cost_[j] += xi;
        } else if (has_upper && !has_lower) {
          cost_[j] -= xi;
        } else {
          cost_[j] += cost_[j] >= 0.0 ? xi : -xi;
        }
        break;
    }
  }
}

Index Simplex::count_dual_infeasibilities_unboxed() const {
  Index count = 0;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    if (status_[j] == VarStatus::kBasic) continue;
    const bool has_lower = lower_[j] > -kInf;
    const bool has_upper = upper_[j] < kInf;
    if (has_lower && has_upper) continue;
    const double dj = d_[j];
    if ((has_lower && dj < -tol) || (has_upper && dj > tol) ||
        (!has_lower && !has_upper && std::fabs(dj) > tol)) {
      ++count;
    }
  }
  return count;
}

double Simplex::max_primal_infeasibility() const {
  double worst = 0.0;
  for (Index k = 0; k < m_; ++k) {
    const Index j = basic_[k];
    worst = std::max({worst, lower_[j] - x_[j], x_[j] - upper_[j]});
  }
  return worst;
}

double Simplex::max_dual_infeasibility() const {
  double worst = 0.0;
  for (Index j = 0; j < nt_; ++j) {
    if (lower_[j] == upper_[j]) continue;
    switch (status_[j]) {
      case VarStatus::kAtLower: worst = std::max(worst, -d_[j]); break;
      case VarStatus::kAtUpper: worst = std::max(worst, d_[j]); break;
      case VarStatus::kAtZero: worst = std::max(worst, std::fabs(d_[j])); break;
      case VarStatus::kBasic: break;
    }
  }
  return worst;
}

bool Simplex::limit_reached(SimplexStatus& status) const {
  if (iterations_ >= max_iterations_) {
    status = SimplexStatus::kIterationLimit;
    return true;
  }
  if ((iterations_ & 31) == 0 && timer_.seconds() > options_.time_limit) {
    status = SimplexStatus::kTimeLimit;
    return true;
  }
  return false;
}

Index Simplex::choose_leaving_row() const {
  Index best = -1;
  double best_score = 0.0;
  const double tol = options_.primal_tol;
  for (Index k = 0; k < m_; ++k) {
    const Index j = basic_[k];
    double infeasibility = 0.0;
    if (x_[j] < lower_[j] - tol) {
      infeasibility = lower_[j] - x_[j];
    } else if (x_[j] > upper_[j] + tol) {
      infeasibility = x_[j] - upper_[j];
    } else {
      continue;
    }
    const double score = infeasibility * infeasibility / dse_weight_[k];
    if (score > best_score) {
      best_score = score;
      best = k;
    }
  }
  return best;
}

Index Simplex::dual_ratio_test(double direction) const {
  // Along the dual step t >= 0, d_j(t) = d_j - t * a_j with a_j = direction * alpha_rj.
  const double tol = options_.dual_tol;
  const double pivot_tol = options_.pivot_tol;
  const auto is_candidate = [&](Index j, double a) {
    const VarStatus st = status_[j];
    if (st == VarStatus::kBasic || lower_[j] == upper_[j] || std::fabs(a) < pivot_tol) return false;
    if (a > 0.0) return st == VarStatus::kAtLower || st == VarStatus::kAtZero;
    return st == VarStatus::kAtUpper || st == VarStatus::kAtZero;
  };

  // Pass 1: largest step keeping every d_j feasible within the tolerance.
  double max_step = kInf;
  for (Index j = 0; j < nt_; ++j) {
    const double a = direction * alpha_row_[j];
    if (!is_candidate(j, a)) continue;
    const double bound = a > 0.0 ? (d_[j] + tol) / a : (d_[j] - tol) / a;
    max_step = std::min(max_step, bound);
  }
  if (max_step == kInf) return -1;

  // Pass 2: among the candidates within that step, the largest pivot.
  Index best = -1;
  double best_abs = 0.0;
  for (Index j = 0; j < nt_; ++j) {
    const double a = direction * alpha_row_[j];
    if (!is_candidate(j, a)) continue;
    if (d_[j] / a <= max_step && std::fabs(a) > best_abs) {
      best_abs = std::fabs(a);
      best = j;
    }
  }
  return best;
}

Index Simplex::choose_entering_column() const {
  Index best = -1;
  double best_score = 0.0;
  const double tol = options_.dual_tol;
  for (Index j = 0; j < nt_; ++j) {
    if (lower_[j] == upper_[j]) continue;
    double score = 0.0;
    switch (status_[j]) {
      case VarStatus::kAtLower: score = -d_[j]; break;
      case VarStatus::kAtUpper: score = d_[j]; break;
      case VarStatus::kAtZero: score = std::fabs(d_[j]); break;
      case VarStatus::kBasic: continue;
    }
    if (score > tol && score > best_score) {
      best_score = score;
      best = j;
    }
  }
  return best;
}

SimplexStatus Simplex::dual_loop() {
  if (!rebuild()) return SimplexStatus::kNumericalError;
  correct_dual_infeasibilities();
  bool fresh = true;
  const auto refresh = [&] {
    if (!rebuild()) return false;
    correct_dual_infeasibilities();
    fresh = true;
    return true;
  };

  for (;;) {
    SimplexStatus limit;
    if (limit_reached(limit)) return limit;
    if (factor_.should_refactor() && !refresh()) return SimplexStatus::kNumericalError;

    const Index r = choose_leaving_row();
    if (r < 0) {
      if (fresh) return SimplexStatus::kOptimal;
      if (!refresh()) return SimplexStatus::kNumericalError;
      continue;
    }
    const Index p = basic_[r];
    const bool to_lower = x_[p] < lower_[p];
    const double bound = to_lower ? lower_[p] : upper_[p];
    const double direction = to_lower ? -1.0 : 1.0;
    const double delta = x_[p] - bound;

    rho_.assign(static_cast<std::size_t>(m_), 0.0);
    rho_[r] = 1.0;
    factor_.btran(rho_);
    compute_pivot_row();

    const Index q = dual_ratio_test(direction);
    if (q < 0) {
      if (!fresh) {
        if (!refresh()) return SimplexStatus::kNumericalError;
        continue;
      }
      dual_ray_ = rho_;
      return SimplexStatus::kInfeasible;
    }

    load_column(q, alpha_col_);
    factor_.ftran(alpha_col_, &spike_);
    const double a_col = alpha_col_[r];
    const double a_row = alpha_row_[q];
    if (std::fabs(a_col - a_row) > kPivotConsistencyTol * (1.0 + std::fabs(a_col))) {
      if (!fresh) {
        ++stats_.rebuilds_on_mismatch;
        if (!refresh()) return SimplexStatus::kNumericalError;
        continue;
      }
      log_.log(3, "simplex: pivot mismatch %.3e vs %.3e after refactor", a_col, a_row);
    }
    if (std::fabs(a_col) < 1e-11) return SimplexStatus::kNumericalError;

    tau_ = rho_;
    factor_.ftran(tau_);

    // Dual step. A slightly infeasible d_q would make the step go the wrong way; shift its cost
    // instead so the step is zero.
    double theta_d = d_[q] / a_row;
    if (direction * theta_d < 0.0) {
      ++stats_.cost_shifts;
      cost_[q] -= d_[q];
      theta_d = 0.0;
    }
    if (theta_d != 0.0) {
      for (Index j = 0; j < nt_; ++j) {
        if (status_[j] != VarStatus::kBasic && alpha_row_[j] != 0.0) d_[j] -= theta_d * alpha_row_[j];
      }
    }
    d_[q] = 0.0;
    d_[p] = -theta_d;

    // Primal step: p moves exactly to the violated bound.
    const double theta_p = delta / a_col;
    for (Index k = 0; k < m_; ++k) x_[basic_[k]] -= theta_p * alpha_col_[k];
    x_[q] += theta_p;
    x_[p] = bound;

    // Dual steepest-edge weights (Forrest–Goldfarb update).
    const double w_r = dse_weight_[r];
    for (Index k = 0; k < m_; ++k) {
      if (k == r || alpha_col_[k] == 0.0) continue;
      const double ratio = alpha_col_[k] / a_col;
      dse_weight_[k] =
          std::max(dse_weight_[k] + ratio * (ratio * w_r - 2.0 * tau_[k]), kMinDseWeight);
    }
    dse_weight_[r] = std::max(w_r / (a_col * a_col), kMinDseWeight);

    basic_[r] = q;
    position_[q] = r;
    status_[q] = VarStatus::kBasic;
    position_[p] = -1;
    status_[p] = (to_lower || lower_[p] == upper_[p]) ? VarStatus::kAtLower : VarStatus::kAtUpper;
    ++iterations_;
    ++stats_.dual_iterations;
    if (factor_.update(r, spike_)) {
      fresh = false;
    } else if (!refresh()) {
      return SimplexStatus::kNumericalError;
    }
  }
}

SimplexStatus Simplex::primal_loop(LoopResult& result) {
  result = LoopResult::kDone;
  if (!rebuild()) return SimplexStatus::kNumericalError;
  bool fresh = true;
  const auto refresh = [&] {
    if (!rebuild()) return false;
    fresh = true;
    return true;
  };
  const double ptol = options_.primal_tol;
  if (max_primal_infeasibility() > ptol) {
    result = LoopResult::kLostFeasibility;
    return SimplexStatus::kOptimal;
  }

  for (;;) {
    SimplexStatus limit;
    if (limit_reached(limit)) return limit;
    if (factor_.should_refactor()) {
      if (!refresh()) return SimplexStatus::kNumericalError;
      if (max_primal_infeasibility() > ptol) {
        result = LoopResult::kLostFeasibility;
        return SimplexStatus::kOptimal;
      }
    }

    const Index q = choose_entering_column();
    if (q < 0) {
      if (fresh) return SimplexStatus::kOptimal;
      if (!refresh()) return SimplexStatus::kNumericalError;
      if (max_primal_infeasibility() > ptol) {
        result = LoopResult::kLostFeasibility;
        return SimplexStatus::kOptimal;
      }
      continue;
    }
    const double dir = d_[q] < 0.0 ? 1.0 : -1.0;
    load_column(q, alpha_col_);
    factor_.ftran(alpha_col_, &spike_);

    // Harris two-pass ratio test on x_B(t) = x_B - t * dir * alpha.
    double max_step = kInf;
    for (Index k = 0; k < m_; ++k) {
      const double a = dir * alpha_col_[k];
      if (std::fabs(a) < options_.pivot_tol) continue;
      const Index j = basic_[k];
      if (a > 0.0 && lower_[j] > -kInf) {
        max_step = std::min(max_step, (x_[j] - lower_[j] + ptol) / a);
      } else if (a < 0.0 && upper_[j] < kInf) {
        max_step = std::min(max_step, (x_[j] - upper_[j] - ptol) / a);
      }
    }
    Index r = -1;
    double step = kInf;
    double best_abs = 0.0;
    if (max_step < kInf) {
      for (Index k = 0; k < m_; ++k) {
        const double a = dir * alpha_col_[k];
        if (std::fabs(a) < options_.pivot_tol) continue;
        const Index j = basic_[k];
        double ratio;
        if (a > 0.0 && lower_[j] > -kInf) {
          ratio = (x_[j] - lower_[j]) / a;
        } else if (a < 0.0 && upper_[j] < kInf) {
          ratio = (x_[j] - upper_[j]) / a;
        } else {
          continue;
        }
        if (ratio <= max_step && std::fabs(a) > best_abs) {
          best_abs = std::fabs(a);
          r = k;
          step = std::max(ratio, 0.0);
        }
      }
    }
    const double range = upper_[q] - lower_[q];

    if (r < 0 && !(range < kInf)) {
      if (!fresh) {
        if (!refresh()) return SimplexStatus::kNumericalError;
        continue;
      }
      primal_ray_.assign(static_cast<std::size_t>(nt_), 0.0);
      primal_ray_[q] = dir;
      for (Index k = 0; k < m_; ++k) primal_ray_[basic_[k]] = -dir * alpha_col_[k];
      return SimplexStatus::kUnbounded;
    }

    if (range <= step) {
      // Bound flip of the entering variable; the basis is unchanged.
      for (Index k = 0; k < m_; ++k) x_[basic_[k]] -= dir * range * alpha_col_[k];
      status_[q] = status_[q] == VarStatus::kAtLower ? VarStatus::kAtUpper : VarStatus::kAtLower;
      set_value_from_status(q);
      ++iterations_;
      ++stats_.primal_iterations;
      ++stats_.bound_flips;
      continue;
    }

    const Index p = basic_[r];
    const double a_col = alpha_col_[r];
    const bool p_to_lower = dir * a_col > 0.0;

    rho_.assign(static_cast<std::size_t>(m_), 0.0);
    rho_[r] = 1.0;
    factor_.btran(rho_);
    compute_pivot_row();
    const double a_row = alpha_row_[q];
    if (std::fabs(a_col - a_row) > kPivotConsistencyTol * (1.0 + std::fabs(a_col)) && !fresh) {
      ++stats_.rebuilds_on_mismatch;
      if (!refresh()) return SimplexStatus::kNumericalError;
      continue;
    }

    const double theta_d = d_[q] / a_row;
    for (Index j = 0; j < nt_; ++j) {
      if (status_[j] != VarStatus::kBasic && alpha_row_[j] != 0.0) d_[j] -= theta_d * alpha_row_[j];
    }
    d_[q] = 0.0;
    d_[p] = -theta_d;

    for (Index k = 0; k < m_; ++k) x_[basic_[k]] -= step * dir * alpha_col_[k];
    x_[q] += step * dir;
    x_[p] = p_to_lower ? lower_[p] : upper_[p];

    basic_[r] = q;
    position_[q] = r;
    status_[q] = VarStatus::kBasic;
    position_[p] = -1;
    status_[p] = (p_to_lower || lower_[p] == upper_[p]) ? VarStatus::kAtLower : VarStatus::kAtUpper;
    dse_weight_[r] = 1.0;
    ++iterations_;
    ++stats_.primal_iterations;
    if (factor_.update(r, spike_)) {
      fresh = false;
    } else if (!refresh()) {
      return SimplexStatus::kNumericalError;
    }
  }
}

std::vector<double> Simplex::exact_dse_weights() const {
  std::vector<double> weights(static_cast<std::size_t>(m_));
  std::vector<double> row;
  for (Index r = 0; r < m_; ++r) {
    row.assign(static_cast<std::size_t>(m_), 0.0);
    row[r] = 1.0;
    factor_.btran(row);
    double sum = 0.0;
    for (const double v : row) sum += v * v;
    weights[r] = sum;
  }
  return weights;
}

SimplexStatus Simplex::dual_phase1() {
  for (Index j = 0; j < nt_; ++j) {
    const bool has_lower = lp_.lower[j] > -kInf;
    const bool has_upper = lp_.upper[j] < kInf;
    if (has_lower && has_upper) {
      lower_[j] = upper_[j] = 0.0;
    } else if (has_lower) {
      lower_[j] = 0.0;
      upper_[j] = 1.0;
    } else if (has_upper) {
      lower_[j] = -1.0;
      upper_[j] = 0.0;
    } else {
      lower_[j] = -kPhase1FreeBound;
      upper_[j] = kPhase1FreeBound;
    }
  }
  place_nonbasic_dual_feasible();
  SimplexStatus status = dual_loop();
  lower_ = lp_.lower;
  upper_ = lp_.upper;
  cost_ = lp_.cost;
  log_.log(2, "simplex: dual phase 1 %s after %lld iterations", to_string(status), iterations_);
  // The auxiliary problem is always feasible (v = 0), so "infeasible" means numerical trouble.
  if (status == SimplexStatus::kInfeasible) return SimplexStatus::kNumericalError;
  if (status != SimplexStatus::kOptimal) return status;
  compute_dual();
  return SimplexStatus::kOptimal;
}

SimplexStatus Simplex::phase2() {
  for (int round = 0; round < 8; ++round) {
    SimplexStatus status = dual_loop();
    log_.log(2, "simplex: dual phase 2 %s after %lld iterations", to_string(status), iterations_);
    if (status != SimplexStatus::kOptimal) return status;

    // Remove perturbations and shifts; the basis stays primal feasible.
    if (cost_ != lp_.cost) {
      cost_ = lp_.cost;
      compute_dual();
    }
    if (max_dual_infeasibility() <= options_.dual_tol) return SimplexStatus::kOptimal;

    LoopResult result;
    status = primal_loop(result);
    log_.log(2, "simplex: primal cleanup %s after %lld iterations", to_string(status), iterations_);
    if (result == LoopResult::kLostFeasibility) {
      std::fill(dse_weight_.begin(), dse_weight_.end(), 1.0);
      continue;
    }
    return status;
  }
  return SimplexStatus::kNumericalError;
}

SimplexStatus Simplex::solve_dual_infeasible() {
  // Decide primal feasibility with a zero objective (trivially dual feasible) ...
  std::fill(cost_.begin(), cost_.end(), 0.0);
  compute_dual();
  place_nonbasic_dual_feasible();
  compute_primal();
  SimplexStatus status = dual_loop();
  log_.log(2, "simplex: feasibility search %s after %lld iterations", to_string(status),
           iterations_);
  if (status != SimplexStatus::kOptimal) return status;

  // ... then the primal simplex from the feasible basis finds the unbounded ray.
  cost_ = lp_.cost;
  compute_dual();
  LoopResult result;
  status = primal_loop(result);
  if (result == LoopResult::kLostFeasibility) return phase2();
  return status;
}

SimplexStatus Simplex::solve() {
  timer_ = Timer();
  iterations_ = 0;
  cost_ = lp_.cost;
  lower_ = lp_.lower;
  upper_ = lp_.upper;

  // Slack basis.
  for (Index i = 0; i < m_; ++i) {
    basic_[i] = n_ + i;
    position_[n_ + i] = i;
    status_[n_ + i] = VarStatus::kBasic;
  }
  for (Index j = 0; j < n_; ++j) {
    position_[j] = -1;
    if (lower_[j] > -kInf) {
      status_[j] = VarStatus::kAtLower;
    } else if (upper_[j] < kInf) {
      status_[j] = VarStatus::kAtUpper;
    } else {
      status_[j] = VarStatus::kAtZero;
    }
    set_value_from_status(j);
  }
  if (!rebuild()) return SimplexStatus::kNumericalError;

  if (count_dual_infeasibilities_unboxed() > 0) {
    const SimplexStatus status = dual_phase1();
    if (status != SimplexStatus::kOptimal) return status;
    if (count_dual_infeasibilities_unboxed() > 0) return solve_dual_infeasible();
  }
  place_nonbasic_dual_feasible();
  compute_primal();
  if (options_.perturb) perturb_costs();
  return phase2();
}

}  // namespace samaya
