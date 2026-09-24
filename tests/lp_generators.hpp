#pragma once

// Random LP generators for differential testing against the reference simplex.

#include <algorithm>
#include <random>
#include <vector>

#include "samaya/model.hpp"

namespace samaya::test {

enum class LpFamily {
  kFeasible,    // Row bounds built around a point inside the column bounds.
  kRandom,      // Independent random bounds: often infeasible.
  kDegenerate,  // Many rows tight at one point, small integer data, many zero costs.
  kLoose,       // Few finite bounds: often unbounded.
};

inline const char* to_string(LpFamily f) {
  switch (f) {
    case LpFamily::kFeasible: return "feasible";
    case LpFamily::kRandom: return "random";
    case LpFamily::kDegenerate: return "degenerate";
    case LpFamily::kLoose: return "loose";
  }
  return "?";
}

inline Model random_lp(LpFamily family, int max_rows, int max_cols, std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto chance = [&](double p) { return std::bernoulli_distribution(p)(rng); };

  Model model;
  const int m = uniform_int(family == LpFamily::kDegenerate ? 2 : 1, max_rows);
  const int n = uniform_int(1, max_cols);
  model.sense = chance(0.3) ? ObjSense::kMaximize : ObjSense::kMinimize;
  model.obj_offset = chance(0.2) ? uniform_int(-10, 10) : 0.0;

  const double density = std::uniform_real_distribution<double>(0.25, 0.9)(rng);
  std::vector<Triplet> t;
  std::vector<std::vector<double>> dense(m, std::vector<double>(n, 0.0));
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (!chance(density)) continue;
      int v = uniform_int(-4, 4);
      if (v == 0) v = 1;
      const double value = family == LpFamily::kDegenerate || chance(0.6)
                               ? v
                               : v * std::uniform_real_distribution<double>(0.5, 1.5)(rng);
      dense[i][j] = value;
      t.push_back({i, j, value});
    }
  }
  model.A = SparseMatrix::from_triplets(m, n, std::move(t));

  const double free_chance = family == LpFamily::kLoose ? 0.5 : 0.1;
  std::vector<double> point(n);
  for (int j = 0; j < n; ++j) {
    const double c = family == LpFamily::kDegenerate && chance(0.4) ? 0 : uniform_int(-5, 5);
    model.obj.push_back(c);
    model.col_type.push_back(VarType::kContinuous);
    double lo = -kInf;
    double up = kInf;
    const double r = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
    const int a = uniform_int(-5, 5);
    const int b = a + uniform_int(0, 6);
    if (r < free_chance) {
      // Free.
    } else if (r < 0.45) {
      lo = family == LpFamily::kLoose ? a : 0;
    } else if (r < 0.65) {
      lo = a;
      up = b;
    } else if (r < 0.75) {
      up = b;
    } else if (r < 0.8) {
      lo = up = a;
    } else {
      lo = a;
      if (family != LpFamily::kLoose) up = a + 10;
    }
    model.col_lower.push_back(lo);
    model.col_upper.push_back(up);
    const double from = lo > -kInf ? lo : (up < kInf ? up - 5 : -5);
    const double to = up < kInf ? up : from + 5;
    point[j] = family == LpFamily::kDegenerate
                   ? std::round(std::uniform_real_distribution<double>(from, to)(rng))
                   : std::uniform_real_distribution<double>(from, to)(rng);
    point[j] = std::clamp(point[j], from, to);
  }

  for (int i = 0; i < m; ++i) {
    double activity = 0.0;
    for (int j = 0; j < n; ++j) activity += dense[i][j] * point[j];
    const double r = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
    double lo = -kInf;
    double up = kInf;
    if (family == LpFamily::kRandom) {
      const int a = uniform_int(-10, 10);
      if (r < 0.3) {
        up = a;
      } else if (r < 0.6) {
        lo = a;
      } else if (r < 0.8) {
        lo = up = a;
      } else {
        lo = a;
        up = a + uniform_int(0, 8);
      }
    } else {
      const bool tight = family == LpFamily::kDegenerate ? chance(0.8) : chance(0.3);
      const double slack = tight ? 0.0 : std::uniform_real_distribution<double>(0.0, 5.0)(rng);
      if (r < 0.35) {
        up = activity + slack;
      } else if (r < 0.7) {
        lo = activity - slack;
      } else if (r < 0.85) {
        lo = up = activity;
      } else if (r < 0.95) {
        lo = activity - slack;
        up = activity + slack + 1.0;
      }
      // Otherwise a free row.
    }
    model.row_lower.push_back(lo);
    model.row_upper.push_back(up);
  }
  return model;
}

}  // namespace samaya::test
