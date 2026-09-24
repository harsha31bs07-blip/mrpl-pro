#pragma once

#include <span>
#include <vector>

#include "samaya/types.hpp"

namespace samaya {

struct Triplet {
  Index row;
  Index col;
  double value;
};

// Compressed sparse column (CSC) matrix. Row indices within each column are sorted and unique,
// and explicit zeros are never stored.
class SparseMatrix {
 public:
  SparseMatrix() = default;
  SparseMatrix(Index rows, Index cols);

  // Builds a CSC matrix from unordered triplets. Duplicate entries are summed; entries that are
  // (or sum to) exactly zero are dropped. Throws std::invalid_argument on out-of-range indices.
  static SparseMatrix from_triplets(Index rows, Index cols, std::vector<Triplet> triplets);

  Index rows() const { return rows_; }
  Index cols() const { return cols_; }
  NnzIndex nnz() const { return static_cast<NnzIndex>(value_.size()); }
  bool empty() const { return value_.empty(); }

  // col_start()[j] .. col_start()[j+1] is the range of column j in row_index() / values().
  std::span<const NnzIndex> col_start() const { return start_; }
  std::span<const Index> row_index() const { return index_; }
  std::span<const double> values() const { return value_; }

  // Returns the transpose, i.e. the same matrix in compressed sparse row form.
  SparseMatrix transpose() const;

  // y = A x. x has cols() entries, y has rows() entries.
  void multiply(std::span<const double> x, std::span<double> y) const;
  // y = A^T x. x has rows() entries, y has cols() entries.
  void multiply_transpose(std::span<const double> x, std::span<double> y) const;

 private:
  Index rows_ = 0;
  Index cols_ = 0;
  std::vector<NnzIndex> start_{0};
  std::vector<Index> index_;
  std::vector<double> value_;
};

}  // namespace samaya
