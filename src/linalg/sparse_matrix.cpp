#include "samaya/sparse_matrix.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>

namespace samaya {

SparseMatrix::SparseMatrix(Index rows, Index cols)
    : rows_(rows), cols_(cols), start_(static_cast<std::size_t>(cols) + 1, 0) {
  if (rows < 0 || cols < 0) throw std::invalid_argument("SparseMatrix: negative dimension");
}

SparseMatrix SparseMatrix::from_triplets(Index rows, Index cols, std::vector<Triplet> triplets) {
  SparseMatrix m(rows, cols);
  for (const Triplet& t : triplets) {
    if (t.row < 0 || t.row >= rows || t.col < 0 || t.col >= cols) {
      throw std::invalid_argument("SparseMatrix: triplet (" + std::to_string(t.row) + ", " +
                                  std::to_string(t.col) + ") out of range");
    }
  }
  std::sort(triplets.begin(), triplets.end(), [](const Triplet& a, const Triplet& b) {
    return a.col != b.col ? a.col < b.col : a.row < b.row;
  });

  m.index_.reserve(triplets.size());
  m.value_.reserve(triplets.size());
  std::size_t k = 0;
  for (Index j = 0; j < cols; ++j) {
    while (k < triplets.size() && triplets[k].col == j) {
      const Index i = triplets[k].row;
      double v = 0.0;
      while (k < triplets.size() && triplets[k].col == j && triplets[k].row == i) {
        v += triplets[k].value;
        ++k;
      }
      if (v != 0.0) {
        m.index_.push_back(i);
        m.value_.push_back(v);
      }
    }
    m.start_[static_cast<std::size_t>(j) + 1] = static_cast<NnzIndex>(m.index_.size());
  }
  return m;
}

SparseMatrix SparseMatrix::transpose() const {
  SparseMatrix t(cols_, rows_);
  t.index_.resize(index_.size());
  t.value_.resize(value_.size());

  // Count entries per row, then prefix-sum into column starts of the transpose.
  for (Index i : index_) ++t.start_[static_cast<std::size_t>(i) + 1];
  for (std::size_t i = 0; i < static_cast<std::size_t>(rows_); ++i) t.start_[i + 1] += t.start_[i];

  // Scattering columns in increasing order keeps row indices of the transpose sorted.
  std::vector<NnzIndex> next(t.start_.begin(), t.start_.end() - 1);
  for (Index j = 0; j < cols_; ++j) {
    for (NnzIndex p = start_[j]; p < start_[j + 1]; ++p) {
      const NnzIndex q = next[index_[p]]++;
      t.index_[q] = j;
      t.value_[q] = value_[p];
    }
  }
  return t;
}

void SparseMatrix::multiply(std::span<const double> x, std::span<double> y) const {
  assert(x.size() == static_cast<std::size_t>(cols_));
  assert(y.size() == static_cast<std::size_t>(rows_));
  std::fill(y.begin(), y.end(), 0.0);
  for (Index j = 0; j < cols_; ++j) {
    const double xj = x[j];
    if (xj == 0.0) continue;
    for (NnzIndex p = start_[j]; p < start_[j + 1]; ++p) y[index_[p]] += value_[p] * xj;
  }
}

void SparseMatrix::multiply_transpose(std::span<const double> x, std::span<double> y) const {
  assert(x.size() == static_cast<std::size_t>(rows_));
  assert(y.size() == static_cast<std::size_t>(cols_));
  for (Index j = 0; j < cols_; ++j) {
    double sum = 0.0;
    for (NnzIndex p = start_[j]; p < start_[j + 1]; ++p) sum += value_[p] * x[index_[p]];
    y[j] = sum;
  }
}

}  // namespace samaya
