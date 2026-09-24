#include "samaya/model_stats.hpp"

#include <cmath>
#include <cstdio>
#include <ostream>

namespace samaya {

void ValueRange::add(double v) {
  const double a = std::fabs(v);
  if (a == 0.0 || !std::isfinite(a)) return;
  if (a < min_abs) min_abs = a;
  if (a > max_abs) max_abs = a;
}

ModelStats compute_stats(const Model& model) {
  ModelStats s;
  s.rows = model.num_rows();
  s.cols = model.num_cols();
  s.nnz = model.A.nnz();
  s.q_nnz = model.Q.nnz();

  for (Index j = 0; j < s.cols; ++j) {
    const double lo = model.col_lower[j];
    const double up = model.col_upper[j];
    if (model.col_type[j] == VarType::kInteger) {
      ++s.integers;
      if (lo >= 0.0 && up <= 1.0) ++s.binaries;
    }
    if (lo == -kInf && up == kInf) ++s.free_cols;
    s.objective.add(model.obj[j]);
    s.bounds.add(lo);
    s.bounds.add(up);
  }
  for (Index i = 0; i < s.rows; ++i) {
    const double lo = model.row_lower[i];
    const double up = model.row_upper[i];
    if (lo == up) {
      ++s.equality_rows;
    } else if (lo > -kInf && up < kInf) {
      ++s.ranged_rows;
    }
    s.rhs.add(lo);
    s.rhs.add(up);
  }
  for (double v : model.A.values()) s.matrix.add(v);
  return s;
}

namespace {

void print_range(std::ostream& out, const char* label, const ValueRange& r) {
  char buf[96];
  if (r.empty()) {
    std::snprintf(buf, sizeof buf, "  %-10s [ - ]\n", label);
  } else {
    std::snprintf(buf, sizeof buf, "  %-10s [%.0e, %.0e]\n", label, r.min_abs, r.max_abs);
  }
  out << buf;
}

}  // namespace

void print_stats(std::ostream& out, const Model& model, const ModelStats& s) {
  out << "Model " << (model.name.empty() ? "<unnamed>" : model.name) << " ("
      << to_string(model.problem_class()) << ", "
      << (model.sense == ObjSense::kMinimize ? "minimize" : "maximize") << ")\n"
      << "  rows " << s.rows << " (" << s.equality_rows << " equality, " << s.ranged_rows
      << " ranged), cols " << s.cols << " (" << s.integers << " integer, " << s.binaries
      << " binary, " << s.free_cols << " free), nonzeros " << s.nnz;
  if (s.q_nnz > 0) out << ", Q nonzeros " << s.q_nnz;
  out << "\nCoefficient ranges\n";
  print_range(out, "matrix", s.matrix);
  print_range(out, "objective", s.objective);
  print_range(out, "bounds", s.bounds);
  print_range(out, "rhs", s.rhs);
}

}  // namespace samaya
