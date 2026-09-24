#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "linalg/scaling.hpp"
#include "test_framework.hpp"

using samaya::Scaling;
using samaya::SparseMatrix;
using samaya::Triplet;

namespace {

bool is_power_of_two(double v) {
  int exponent = 0;
  return v > 0.0 && std::frexp(v, &exponent) == 0.5;
}

// Random sparse matrix whose rows and columns carry wildly different magnitudes.
SparseMatrix badly_scaled_matrix(int m, int n, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unit(0.5, 2.0);
  std::uniform_int_distribution<int> exponent(-6, 6);
  std::vector<double> rs(m), cs(n);
  for (double& r : rs) r = std::pow(10.0, exponent(rng));
  for (double& c : cs) c = std::pow(10.0, exponent(rng));
  std::vector<Triplet> t;
  std::bernoulli_distribution keep(0.3);
  for (int j = 0; j < n; ++j) {
    t.push_back({j % m, j, rs[j % m] * cs[j] * unit(rng)});
    for (int i = 0; i < m; ++i) {
      if (keep(rng)) t.push_back({i, j, rs[i] * cs[j] * unit(rng)});
    }
  }
  return SparseMatrix::from_triplets(m, n, std::move(t));
}

}  // namespace

TEST(scaling_factors_are_powers_of_two_and_improve_ratio) {
  for (unsigned seed = 1; seed <= 20; ++seed) {
    const SparseMatrix a = badly_scaled_matrix(15, 20, seed);
    const Scaling s = samaya::compute_scaling(a);
    for (double f : s.row) CHECK(is_power_of_two(f));
    for (double f : s.col) CHECK(is_power_of_two(f));
    const double before = samaya::scaled_magnitude_ratio(a, Scaling::identity(15, 20));
    const double after = samaya::scaled_magnitude_ratio(a, s);
    CHECK(before > 1e6);
    CHECK(after < 1e3);
  }
}

TEST(scaling_equilibrates_row_and_column_max) {
  const SparseMatrix a = badly_scaled_matrix(10, 12, 7);
  const SparseMatrix scaled = samaya::compute_scaling(a).apply(a);
  std::vector<double> col_max(12, 0.0), row_max(10, 0.0);
  for (int j = 0; j < 12; ++j) {
    for (auto p = scaled.col_start()[j]; p < scaled.col_start()[j + 1]; ++p) {
      const double v = std::fabs(scaled.values()[p]);
      col_max[j] = std::max(col_max[j], v);
      row_max[scaled.row_index()[p]] = std::max(row_max[scaled.row_index()[p]], v);
    }
  }
  // Power-of-two rounding keeps each max within a factor of two of 1.
  for (double v : col_max) CHECK(v > 0.49 && v < 2.01);
  for (double v : row_max) CHECK(v > 0.2 && v < 4.01);
}

TEST(scaling_round_trip_is_exact) {
  const SparseMatrix a = badly_scaled_matrix(8, 9, 3);
  const Scaling s = samaya::compute_scaling(a);
  std::vector<double> x{1.1, -2.3, 0.0, 5e7, -1e-9, samaya::kInf, -samaya::kInf, 3.0, 4.0};
  const std::vector<double> original = x;
  s.scale_cols(x);
  s.unscale_cols(x);
  CHECK(x == original);
  std::vector<double> r(original.begin(), original.begin() + 8);
  const std::vector<double> r0 = r;
  s.scale_rows(r);
  s.unscale_rows(r);
  CHECK(r == r0);
}

TEST(scaling_handles_empty_and_zero_rows) {
  const SparseMatrix empty(3, 4);
  const Scaling s = samaya::compute_scaling(empty);
  CHECK_EQ(s.row.size(), 3u);
  CHECK_EQ(s.col[3], 1.0);
  const SparseMatrix a = SparseMatrix::from_triplets(3, 2, {{0, 0, 1e4}, {0, 1, 1e-4}});
  const Scaling t = samaya::compute_scaling(a);
  CHECK_EQ(t.row[1], 1.0);
  CHECK(samaya::scaled_magnitude_ratio(a, t) <= 2.0);
}
