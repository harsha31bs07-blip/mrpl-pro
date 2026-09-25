#pragma once

// Textbook dense two-phase tableau simplex with Bland's rule, in long double. It shares no code
// with the production solver and is used only to cross-check it on small LPs.

#include <cmath>
#include <utility>
#include <vector>

#include "samaya/model.hpp"

namespace samaya::test {

struct ReferenceResult {
  enum class Status { kOptimal, kInfeasible, kUnbounded } status = Status::kInfeasible;
  double objective = 0.0;
  std::vector<double> x;  // An optimal point (optimal status only).
};

class ReferenceLp {
 public:
  explicit ReferenceLp(const Model& model) : model_(model) {}

  ReferenceResult solve() {
    build_standard_form();
    return run();
  }

 private:
  using Real = long double;
  static constexpr Real kEps = 1e-9L;

  struct Row {
    std::vector<Real> coef;  // Over standard-form variables.
    Real rhs;
    bool equality;
  };

  // x_j = shift_j + sum_k coef * p_k with p >= 0.
  struct Substitution {
    Real shift = 0;
    std::vector<std::pair<int, Real>> terms;
  };

  void build_standard_form() {
    const Index n = model_.num_cols();
    subs_.resize(static_cast<std::size_t>(n));
    std::vector<std::pair<int, Real>> upper_limits;  // p_k <= value
    for (Index j = 0; j < n; ++j) {
      const double lo = model_.col_lower[j];
      const double up = model_.col_upper[j];
      Substitution& s = subs_[j];
      if (lo > -kInf) {
        s.shift = lo;
        s.terms.push_back({num_vars_, 1});
        if (up < kInf) upper_limits.push_back({num_vars_, static_cast<Real>(up) - lo});
        ++num_vars_;
      } else if (up < kInf) {
        s.shift = up;
        s.terms.push_back({num_vars_++, -1});
      } else {
        s.terms.push_back({num_vars_++, 1});
        s.terms.push_back({num_vars_++, -1});
      }
    }

    const Real sense = model_.sense == ObjSense::kMaximize ? -1 : 1;
    cost_.assign(static_cast<std::size_t>(num_vars_), 0);
    constant_ = sense * static_cast<Real>(model_.obj_offset);
    for (Index j = 0; j < n; ++j) {
      const Real c = sense * static_cast<Real>(model_.obj[j]);
      constant_ += c * subs_[j].shift;
      for (const auto& [k, f] : subs_[j].terms) cost_[k] += c * f;
    }

    // Dense rows of A.
    const Index m = model_.num_rows();
    std::vector<std::vector<Real>> dense(static_cast<std::size_t>(m), std::vector<Real>(n, 0));
    for (Index j = 0; j < n; ++j) {
      for (auto p = model_.A.col_start()[j]; p < model_.A.col_start()[j + 1]; ++p) {
        dense[model_.A.row_index()[p]][j] = model_.A.values()[p];
      }
    }
    for (Index i = 0; i < m; ++i) {
      std::vector<Real> coef(static_cast<std::size_t>(num_vars_), 0);
      Real shift = 0;
      for (Index j = 0; j < n; ++j) {
        if (dense[i][j] == 0) continue;
        shift += dense[i][j] * subs_[j].shift;
        for (const auto& [k, f] : subs_[j].terms) coef[k] += dense[i][j] * f;
      }
      const double lo = model_.row_lower[i];
      const double up = model_.row_upper[i];
      if (lo == up) {
        rows_.push_back({coef, static_cast<Real>(up) - shift, true});
        continue;
      }
      if (up < kInf) rows_.push_back({coef, static_cast<Real>(up) - shift, false});
      if (lo > -kInf) {
        std::vector<Real> neg = coef;
        for (Real& v : neg) v = -v;
        rows_.push_back({neg, -(static_cast<Real>(lo) - shift), false});
      }
    }
    for (const auto& [k, limit] : upper_limits) {
      std::vector<Real> coef(static_cast<std::size_t>(num_vars_), 0);
      coef[k] = 1;
      rows_.push_back({coef, limit, false});
    }
  }

  // Pivots the tableau on (r, c).
  void pivot(int r, int c) {
    const Real p = tab_[r][c];
    for (Real& v : tab_[r]) v /= p;
    for (std::size_t i = 0; i < tab_.size(); ++i) {
      if (static_cast<int>(i) == r) continue;
      const Real f = tab_[i][c];
      if (f == 0) continue;
      for (std::size_t k = 0; k < tab_[i].size(); ++k) tab_[i][k] -= f * tab_[r][k];
    }
    basis_[r] = c;
  }

  // Minimizes cost over the current tableau with Bland's rule. Returns false if unbounded.
  bool optimize(const std::vector<Real>& cost, const std::vector<bool>& allowed) {
    const int cols = static_cast<int>(cost.size());
    for (;;) {
      int enter = -1;
      for (int j = 0; j < cols && enter < 0; ++j) {
        if (!allowed[j]) continue;
        Real rc = cost[j];
        for (std::size_t r = 0; r < tab_.size(); ++r) rc -= cost[basis_[r]] * tab_[r][j];
        if (rc < -kEps) enter = j;
      }
      if (enter < 0) return true;
      int leave = -1;
      Real best = 0;
      for (std::size_t r = 0; r < tab_.size(); ++r) {
        if (tab_[r][enter] <= kEps) continue;
        const Real ratio = tab_[r][cols] / tab_[r][enter];
        if (leave < 0 || ratio < best - kEps ||
            (ratio <= best + kEps && basis_[r] < basis_[leave])) {
          best = ratio;
          leave = static_cast<int>(r);
        }
      }
      if (leave < 0) return false;
      pivot(leave, enter);
    }
  }

  ReferenceResult run() {
    // Columns: standard variables, one slack/surplus per inequality, one artificial per row
    // that needs it; the last entry of each tableau row is the right-hand side.
    const int nv = num_vars_;
    int num_slack = 0;
    for (const Row& r : rows_) num_slack += r.equality ? 0 : 1;
    const int rows = static_cast<int>(rows_.size());
    const int first_art = nv + num_slack;
    const int cols = first_art + rows;
    tab_.assign(static_cast<std::size_t>(rows), std::vector<Real>(cols + 1, 0));
    basis_.assign(static_cast<std::size_t>(rows), -1);
    int slack = nv;
    for (int r = 0; r < rows; ++r) {
      const Row& row = rows_[r];
      const Real sign = row.rhs < 0 ? -1 : 1;
      for (int k = 0; k < nv; ++k) tab_[r][k] = sign * row.coef[k];
      tab_[r][cols] = sign * row.rhs;
      if (!row.equality) {
        tab_[r][slack] = sign;
        if (sign > 0) basis_[r] = slack;
        ++slack;
      }
      if (basis_[r] < 0) {
        tab_[r][first_art + r] = 1;
        basis_[r] = first_art + r;
      }
    }

    std::vector<bool> allowed(static_cast<std::size_t>(cols), true);
    std::vector<Real> phase1(static_cast<std::size_t>(cols), 0);
    for (int k = first_art; k < cols; ++k) phase1[k] = 1;
    optimize(phase1, allowed);
    Real infeasibility = 0;
    for (int r = 0; r < rows; ++r) {
      if (basis_[r] >= first_art) infeasibility += tab_[r][cols];
    }
    if (infeasibility > 1e-7L) return {ReferenceResult::Status::kInfeasible, 0.0, {}};

    // Drive artificials out of the basis; drop rows that turn out redundant.
    for (int r = rows - 1; r >= 0; --r) {
      if (basis_[r] < first_art) continue;
      int col = -1;
      for (int k = 0; k < first_art && col < 0; ++k) {
        if (std::fabs(tab_[r][k]) > kEps) col = k;
      }
      if (col >= 0) {
        pivot(r, col);
      } else {
        tab_.erase(tab_.begin() + r);
        basis_.erase(basis_.begin() + r);
      }
    }

    for (int k = first_art; k < cols; ++k) allowed[k] = false;
    std::vector<Real> phase2(static_cast<std::size_t>(cols), 0);
    for (int k = 0; k < nv; ++k) phase2[k] = cost_[k];
    if (!optimize(phase2, allowed)) return {ReferenceResult::Status::kUnbounded, 0.0, {}};

    Real objective = constant_;
    for (std::size_t r = 0; r < tab_.size(); ++r) objective += phase2[basis_[r]] * tab_[r][cols];
    if (model_.sense == ObjSense::kMaximize) objective = -objective;
    std::vector<Real> p(static_cast<std::size_t>(cols), 0);
    for (std::size_t r = 0; r < tab_.size(); ++r) p[basis_[r]] = tab_[r][cols];
    std::vector<double> x(subs_.size());
    for (std::size_t j = 0; j < subs_.size(); ++j) {
      Real v = subs_[j].shift;
      for (const auto& [k, f] : subs_[j].terms) v += f * p[k];
      x[j] = static_cast<double>(v);
    }
    return {ReferenceResult::Status::kOptimal, static_cast<double>(objective), std::move(x)};
  }

  const Model& model_;
  std::vector<Substitution> subs_;
  int num_vars_ = 0;
  std::vector<Real> cost_;
  Real constant_ = 0;
  std::vector<Row> rows_;
  std::vector<std::vector<Real>> tab_;
  std::vector<int> basis_;
};

}  // namespace samaya::test
