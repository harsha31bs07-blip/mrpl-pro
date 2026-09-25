#include "linalg/scaling.hpp"

#include <algorithm>
#include <cmath>

namespace samaya {
namespace {

// Nearest power of two to a positive finite value.
double round_to_power_of_two(double v) {
  if (!(v > 0.0) || !std::isfinite(v)) return 1.0;
  return std::ldexp(1.0, static_cast<int>(std::lround(std::log2(v))));
}

// Clamp scale factors so that pathological input cannot overflow the scaled model.
double clamp_factor(double v) { return std::clamp(v, 0x1p-40, 0x1p40); }

}  // namespace

Scaling Scaling::identity(Index rows, Index cols) {
  Scaling s;
  s.row.assign(static_cast<std::size_t>(rows), 1.0);
  s.col.assign(static_cast<std::size_t>(cols), 1.0);
  return s;
}

SparseMatrix Scaling::apply(const SparseMatrix& A) const {
  std::vector<Triplet> t;
  t.reserve(static_cast<std::size_t>(A.nnz()));
  const auto start = A.col_start();
  const auto index = A.row_index();
  const auto value = A.values();
  for (Index j = 0; j < A.cols(); ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      t.push_back({index[p], j, value[p] * row[index[p]] * col[j]});
    }
  }
  return SparseMatrix::from_triplets(A.rows(), A.cols(), std::move(t));
}

void Scaling::scale_cols(std::span<double> v) const {
  for (std::size_t j = 0; j < v.size(); ++j) v[j] /= col[j];
}
void Scaling::unscale_cols(std::span<double> v) const {
  for (std::size_t j = 0; j < v.size(); ++j) v[j] *= col[j];
}
void Scaling::scale_rows(std::span<double> v) const {
  for (std::size_t i = 0; i < v.size(); ++i) v[i] *= row[i];
}
void Scaling::unscale_rows(std::span<double> v) const {
  for (std::size_t i = 0; i < v.size(); ++i) v[i] /= row[i];
}
void Scaling::scale_costs(std::span<double> v) const {
  for (std::size_t j = 0; j < v.size(); ++j) v[j] *= col[j];
}
void Scaling::unscale_row_duals(std::span<double> v) const {
  for (std::size_t i = 0; i < v.size(); ++i) v[i] *= row[i];
}
void Scaling::unscale_col_duals(std::span<double> v) const {
  for (std::size_t j = 0; j < v.size(); ++j) v[j] /= col[j];
}

double scaled_magnitude_ratio(const SparseMatrix& A, const Scaling& s) {
  double lo = kInf;
  double hi = 0.0;
  const auto start = A.col_start();
  const auto index = A.row_index();
  const auto value = A.values();
  for (Index j = 0; j < A.cols(); ++j) {
    for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
      const double a = std::fabs(value[p]) * s.row[index[p]] * s.col[j];
      lo = std::min(lo, a);
      hi = std::max(hi, a);
    }
  }
  return hi > 0.0 ? hi / lo : 1.0;
}

Scaling compute_scaling(const SparseMatrix& A, const ScalingOptions& options) {
  const Index m = A.rows();
  const Index n = A.cols();
  Scaling s = Scaling::identity(m, n);
  if (A.nnz() == 0) return s;

  const auto start = A.col_start();
  const auto index = A.row_index();
  const auto value = A.values();
  std::vector<double> row_min(static_cast<std::size_t>(m));
  std::vector<double> row_max(static_cast<std::size_t>(m));

  const auto scan_rows = [&] {
    std::fill(row_min.begin(), row_min.end(), kInf);
    std::fill(row_max.begin(), row_max.end(), 0.0);
    for (Index j = 0; j < n; ++j) {
      for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
        const double a = std::fabs(value[p]) * s.col[j];
        row_min[index[p]] = std::min(row_min[index[p]], a);
        row_max[index[p]] = std::max(row_max[index[p]], a);
      }
    }
  };

  // Iterated geometric-mean scaling: each row, then each column, is divided by
  // sqrt(min * max) of its current magnitudes.
  double ratio = scaled_magnitude_ratio(A, s);
  for (int pass = 0; pass < options.max_geometric_passes && ratio > 1.0; ++pass) {
    const Scaling previous = s;
    scan_rows();
    for (Index i = 0; i < m; ++i) {
      if (row_max[i] > 0.0) s.row[i] = clamp_factor(1.0 / std::sqrt(row_min[i] * row_max[i]));
    }
    for (Index j = 0; j < n; ++j) {
      double lo = kInf;
      double hi = 0.0;
      for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
        const double a = std::fabs(value[p]) * s.row[index[p]];
        lo = std::min(lo, a);
        hi = std::max(hi, a);
      }
      if (hi > 0.0) s.col[j] = clamp_factor(1.0 / std::sqrt(lo * hi));
    }
    const double new_ratio = scaled_magnitude_ratio(A, s);
    if (new_ratio > ratio) {
      s = previous;  // Never make things worse.
      break;
    }
    const bool small_gain = new_ratio > options.min_improvement * ratio;
    ratio = new_ratio;
    if (small_gain) break;
  }

  if (options.equilibrate) {
    scan_rows();
    for (Index i = 0; i < m; ++i) {
      if (row_max[i] > 0.0) s.row[i] = clamp_factor(1.0 / row_max[i]);
    }
    for (Index j = 0; j < n; ++j) {
      double hi = 0.0;
      for (NnzIndex p = start[j]; p < start[j + 1]; ++p) {
        hi = std::max(hi, std::fabs(value[p]) * s.row[index[p]]);
      }
      if (hi > 0.0) s.col[j] = clamp_factor(1.0 / hi);
    }
  }

  for (double& f : s.row) f = round_to_power_of_two(f);
  for (double& f : s.col) f = round_to_power_of_two(f);
  return s;
}

}  // namespace samaya
