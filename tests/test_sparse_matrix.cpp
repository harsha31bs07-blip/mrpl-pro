#include <stdexcept>
#include <vector>

#include "samaya/sparse_matrix.hpp"
#include "test_framework.hpp"

using samaya::SparseMatrix;
using samaya::Triplet;

namespace {

// [ 1 0 2 ]
// [ 0 3 0 ]
SparseMatrix small_matrix() {
  return SparseMatrix::from_triplets(2, 3, {{1, 1, 3.0}, {0, 2, 2.0}, {0, 0, 1.0}});
}

}  // namespace

TEST(sparse_from_triplets_sorts_into_csc) {
  const SparseMatrix a = small_matrix();
  CHECK_EQ(a.rows(), 2);
  CHECK_EQ(a.cols(), 3);
  CHECK_EQ(a.nnz(), 3);
  const std::vector<samaya::NnzIndex> start(a.col_start().begin(), a.col_start().end());
  CHECK(start == (std::vector<samaya::NnzIndex>{0, 1, 2, 3}));
  CHECK_EQ(a.row_index()[2], 0);
  CHECK_EQ(a.values()[1], 3.0);
}

TEST(sparse_from_triplets_sums_duplicates_and_drops_zeros) {
  const SparseMatrix a = SparseMatrix::from_triplets(
      3, 2, {{2, 0, 1.0}, {0, 0, 5.0}, {2, 0, 2.0}, {1, 1, 4.0}, {1, 1, -4.0}, {0, 1, 0.0}});
  CHECK_EQ(a.nnz(), 2);
  CHECK_EQ(a.row_index()[0], 0);
  CHECK_EQ(a.row_index()[1], 2);
  CHECK_EQ(a.values()[1], 3.0);
  CHECK_EQ(a.col_start()[2], 2);
}

TEST(sparse_from_triplets_rejects_out_of_range) {
  CHECK_THROWS(SparseMatrix::from_triplets(2, 2, {{2, 0, 1.0}}), std::invalid_argument);
  CHECK_THROWS(SparseMatrix::from_triplets(2, 2, {{0, -1, 1.0}}), std::invalid_argument);
}

TEST(sparse_transpose) {
  const SparseMatrix t = small_matrix().transpose();
  CHECK_EQ(t.rows(), 3);
  CHECK_EQ(t.cols(), 2);
  CHECK_EQ(t.nnz(), 3);
  // Column 0 of the transpose is row 0 of the original: entries at rows 0 and 2.
  CHECK_EQ(t.col_start()[1], 2);
  CHECK_EQ(t.row_index()[0], 0);
  CHECK_EQ(t.row_index()[1], 2);
  CHECK_EQ(t.values()[1], 2.0);
  CHECK_EQ(t.row_index()[2], 1);
}

TEST(sparse_multiply) {
  const SparseMatrix a = small_matrix();
  const std::vector<double> x{1.0, 2.0, 3.0};
  std::vector<double> y(2, -1.0);
  a.multiply(x, y);
  CHECK_EQ(y[0], 7.0);
  CHECK_EQ(y[1], 6.0);

  const std::vector<double> u{1.0, -1.0};
  std::vector<double> v(3);
  a.multiply_transpose(u, v);
  CHECK_EQ(v[0], 1.0);
  CHECK_EQ(v[1], -3.0);
  CHECK_EQ(v[2], 2.0);
}

TEST(sparse_empty_matrix) {
  const SparseMatrix a(0, 0);
  CHECK(a.empty());
  CHECK_EQ(a.transpose().nnz(), 0);
}
