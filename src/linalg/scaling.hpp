#pragma once

#include <span>
#include <vector>

#include "samaya/sparse_matrix.hpp"

namespace samaya {

// Row and column scale factors for a constraint matrix. The scaled matrix is R A C with
// R = diag(row) and C = diag(col). Every factor is an exact power of two, so scaling and unscaling
// introduce no rounding error.
//
// For a model  min c'x  s.t.  L <= A x <= U,  l <= x <= u  the scaled model in x~ = C^-1 x is
//   min (C c)'x~  s.t.  R L <= (R A C) x~ <= R U,  C^-1 l <= x~ <= C^-1 u.
// Its solution maps back as  x = C x~,  row activity = R^-1 activity~,  row duals y = R y~,
// reduced costs d = C^-1 d~.
struct Scaling {
  std::vector<double> col;
  std::vector<double> row;

  static Scaling identity(Index rows, Index cols);

  // Scaled matrix R A C.
  SparseMatrix apply(const SparseMatrix& A) const;

  // In-place conversions between original and scaled quantities (infinities are preserved).
  void scale_cols(std::span<double> col_values) const;    // bounds, primal x: divide by col
  void unscale_cols(std::span<double> col_values) const;  // multiply by col
  void scale_rows(std::span<double> row_values) const;    // row bounds, activities: multiply by row
  void unscale_rows(std::span<double> row_values) const;  // divide by row
  void scale_costs(std::span<double> costs) const;        // multiply by col
  void unscale_row_duals(std::span<double> duals) const;  // multiply by row
  void unscale_col_duals(std::span<double> duals) const;  // divide by col
};

struct ScalingOptions {
  int max_geometric_passes = 10;
  // Stop geometric passes once a pass improves the max/min magnitude ratio by less than this.
  double min_improvement = 0.9;
  // Finish with one row + column equilibration pass so every row and column max is near 1.
  bool equilibrate = true;
};

// Iterated geometric-mean scaling followed by equilibration, rounded to powers of two.
Scaling compute_scaling(const SparseMatrix& A, const ScalingOptions& options = {});

// max |a_ij| / min |a_ij| over the nonzeros of R A C (1 for an empty matrix).
double scaled_magnitude_ratio(const SparseMatrix& A, const Scaling& scaling);

}  // namespace samaya
