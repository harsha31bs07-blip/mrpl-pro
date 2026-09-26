#pragma once

#include <vector>

#include "lp/simplex.hpp"
#include "samaya/model.hpp"

namespace samaya {

// A cutting plane  sum_k value[k] * x[index[k]] >= lower  in the original (unscaled) columns.
struct Cut {
  std::vector<Index> index;
  std::vector<double> value;
  double lower = 0.0;
  double efficacy = 0.0;  // Violation by the LP point divided by the coefficient norm.
  int pool_index = -1;    // Position in the search's cut pool when taken from there.
};

// What the separators see: the model's rows (including earlier cuts), the column bounds in
// force at the root, and the LP solution, all in original units.
struct CutContext {
  const Model& model;
  const SparseMatrix& At;  // model.A row-wise.
  const std::vector<double>& lower;
  const std::vector<double>& upper;
  const std::vector<double>& x;
  Index original_rows;  // Rows before the first cut; only these seed MIR and cover cuts.
};

// Gomory mixed-integer cut from the tableau row of an integer column k that is basic at a
// fractional value. `row` is the tableau row in original units over all n + m variables (row
// activities for the logicals), `status` the nonbasic statuses. Returns false if no safe cut.
bool gomory_mixed_integer_cut(const CutContext& ctx, Index k, const std::vector<double>& row,
                              const std::vector<VarStatus>& status, Cut& cut);

// Complemented mixed-integer rounding (c-MIR) cuts from single rows, trying several divisors.
void separate_mir(const CutContext& ctx, std::vector<Cut>& cuts);

// Aggregated c-MIR (Marchand and Wolsey): rows are combined to eliminate continuous columns far
// from their bounds, and each aggregation is tried with bound substitution, including variable
// bounds x <= c y and x >= c y on binaries (flow-cover strength on fixed-charge models).
void separate_aggregated_mir(const CutContext& ctx, std::vector<Cut>& cuts);

// Lifted (extended) knapsack cover cuts from rows over binary columns; other columns are
// relaxed to their bounds.
void separate_knapsack_covers(const CutContext& ctx, std::vector<Cut>& cuts);

// Makes a cut numerically safe or rejects it: drops tiny coefficients (relaxing the right-hand
// side through the column bounds), rejects large coefficient ranges, relaxes the right-hand side
// slightly and computes the efficacy. Returns false if the cut should not be used.
bool finalize_cut(const CutContext& ctx, Cut& cut);

}  // namespace samaya
