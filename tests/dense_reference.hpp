#pragma once

// Dense reference implementations used to cross-check the sparse production code. They favour
// obviously-correct over fast and are only meant for small problems.

#include <cmath>
#include <utility>
#include <vector>

namespace samaya::test {

using DenseMatrix = std::vector<std::vector<double>>;  // Row-major.

// Solves M z = b by Gaussian elimination with partial pivoting. Returns false if M is singular.
inline bool dense_solve(DenseMatrix M, std::vector<double> b, std::vector<double>& z) {
  const std::size_t n = b.size();
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; ++r) {
      if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
    }
    if (std::fabs(M[piv][c]) < 1e-13) return false;
    std::swap(M[piv], M[c]);
    std::swap(b[piv], b[c]);
    for (std::size_t r = c + 1; r < n; ++r) {
      const double f = M[r][c] / M[c][c];
      if (f == 0.0) continue;
      for (std::size_t k = c; k < n; ++k) M[r][k] -= f * M[c][k];
      b[r] -= f * b[c];
    }
  }
  z.assign(n, 0.0);
  for (std::size_t r = n; r-- > 0;) {
    double s = b[r];
    for (std::size_t k = r + 1; k < n; ++k) s -= M[r][k] * z[k];
    z[r] = s / M[r][r];
  }
  return true;
}

inline DenseMatrix transpose(const DenseMatrix& M) {
  DenseMatrix t(M.empty() ? 0 : M[0].size(), std::vector<double>(M.size()));
  for (std::size_t r = 0; r < M.size(); ++r) {
    for (std::size_t c = 0; c < M[r].size(); ++c) t[c][r] = M[r][c];
  }
  return t;
}

inline std::vector<double> multiply(const DenseMatrix& M, const std::vector<double>& x) {
  std::vector<double> y(M.size(), 0.0);
  for (std::size_t r = 0; r < M.size(); ++r) {
    for (std::size_t c = 0; c < x.size(); ++c) y[r] += M[r][c] * x[c];
  }
  return y;
}

}  // namespace samaya::test
