#pragma once

#include <iosfwd>

#include "samaya/model.hpp"

namespace samaya {

// Magnitude range of the finite nonzero values in a set; min_abs > max_abs when the set is empty.
struct ValueRange {
  double min_abs = kInf;
  double max_abs = 0.0;

  void add(double v);
  bool empty() const { return min_abs > max_abs; }
};

// Size and numerics summary of a model. The coefficient ranges are the first thing to look at
// when diagnosing an ill-conditioned model.
struct ModelStats {
  Index rows = 0;
  Index cols = 0;
  NnzIndex nnz = 0;
  NnzIndex q_nnz = 0;
  Index integers = 0;
  Index binaries = 0;
  Index free_cols = 0;
  Index equality_rows = 0;
  Index ranged_rows = 0;
  ValueRange matrix;
  ValueRange objective;
  ValueRange bounds;
  ValueRange rhs;
};

ModelStats compute_stats(const Model& model);
void print_stats(std::ostream& out, const Model& model, const ModelStats& stats);

}  // namespace samaya
