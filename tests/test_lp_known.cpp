// Small LPs with known answers, including classic degenerate / cycling examples.

#include <cmath>
#include <string>
#include <vector>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "samaya/io.hpp"
#include "test_framework.hpp"

using samaya::kInf;
using samaya::LpResult;
using samaya::Model;
using samaya::SimplexStatus;

namespace {

LpResult solve(const Model& model) {
  const samaya::Logger quiet(0);
  return samaya::solve_lp(model, samaya::LpSolveOptions{}, quiet);
}

Model make_model(int m, int n, std::vector<samaya::Triplet> a, std::vector<double> obj,
                 std::vector<double> col_lower, std::vector<double> col_upper,
                 std::vector<double> row_lower, std::vector<double> row_upper) {
  Model model;
  model.A = samaya::SparseMatrix::from_triplets(m, n, std::move(a));
  model.obj = std::move(obj);
  model.col_lower = std::move(col_lower);
  model.col_upper = std::move(col_upper);
  model.col_type.assign(static_cast<std::size_t>(n), samaya::VarType::kContinuous);
  model.row_lower = std::move(row_lower);
  model.row_upper = std::move(row_upper);
  return model;
}

}  // namespace

TEST(lp_known_tiny_lp_file) {
  const Model model = samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/tiny_lp.mps");
  const LpResult r = solve(model);
  REQUIRE(r.status == SimplexStatus::kOptimal);
  CHECK_NEAR(r.objective, 11.0, 1e-9);
  CHECK_NEAR(r.col_value[0], 3.0, 1e-9);
  CHECK_NEAR(r.col_value[1], 1.0, 1e-9);
}

TEST(lp_known_beale_cycling_example) {
  // Beale (1955): cycles under Dantzig's rule without anti-cycling measures.
  //   min -3/4 x1 + 20 x2 - 1/2 x3 + 6 x4
  //   s.t. 1/4 x1 -  8 x2 -     x3 + 9 x4 <= 0
  //        1/2 x1 - 12 x2 - 1/2 x3 + 3 x4 <= 0
  //                             x3        <= 1,   x >= 0.   Optimum -5/4.
  const Model model = make_model(
      3, 4,
      {{0, 0, 0.25}, {0, 1, -8}, {0, 2, -1}, {0, 3, 9}, {1, 0, 0.5}, {1, 1, -12}, {1, 2, -0.5},
       {1, 3, 3}, {2, 2, 1}},
      {-0.75, 20, -0.5, 6}, {0, 0, 0, 0}, {kInf, kInf, kInf, kInf}, {-kInf, -kInf, -kInf},
      {0, 0, 1});
  const LpResult r = solve(model);
  REQUIRE(r.status == SimplexStatus::kOptimal);
  CHECK_NEAR(r.objective, -1.25, 1e-9);
}

TEST(lp_known_infeasible_with_certificate) {
  // x + y >= 4, x + y <= 2.
  const Model model = make_model(2, 2, {{0, 0, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}}, {1, 1},
                                 {0, 0}, {kInf, kInf}, {4, -kInf}, {kInf, 2});
  const LpResult r = solve(model);
  REQUIRE(r.status == SimplexStatus::kInfeasible);
  REQUIRE(r.dual_ray.size() == 2);
  // Any Farkas multiplier must weight both rows equally with opposite signs.
  CHECK(std::fabs(r.dual_ray[0] + r.dual_ray[1]) <= 1e-9 * std::fabs(r.dual_ray[0]));
  CHECK(r.dual_ray[0] != 0.0);
}

TEST(lp_known_unbounded_with_ray) {
  // min -x - y  s.t.  x - y <= 1,  x, y >= 0.
  const Model model =
      make_model(1, 2, {{0, 0, 1}, {0, 1, -1}}, {-1, -1}, {0, 0}, {kInf, kInf}, {-kInf}, {1});
  const LpResult r = solve(model);
  REQUIRE(r.status == SimplexStatus::kUnbounded);
  REQUIRE(r.primal_ray.size() == 2);
  CHECK(-r.primal_ray[0] - r.primal_ray[1] < 0.0);  // Improving.
  CHECK(r.primal_ray[0] >= 0.0 && r.primal_ray[1] >= 0.0);
  CHECK(r.primal_ray[0] - r.primal_ray[1] <= 1e-12);  // Stays feasible for the row.
}

TEST(lp_known_models_without_rows) {
  Model model = make_model(0, 2, {}, {1, -1}, {2, -3}, {5, 4}, {}, {});
  LpResult r = solve(model);
  REQUIRE(r.status == SimplexStatus::kOptimal);
  CHECK_NEAR(r.objective, 2.0 - 4.0, 1e-12);

  model.sense = samaya::ObjSense::kMaximize;
  r = solve(model);
  REQUIRE(r.status == SimplexStatus::kOptimal);
  CHECK_NEAR(r.objective, 5.0 + 3.0, 1e-12);

  model.col_upper[0] = kInf;
  r = solve(model);
  CHECK(r.status == SimplexStatus::kUnbounded);
}

TEST(lp_known_free_variables_and_equalities) {
  // min x + 2y + 3z  s.t.  x + y + z = 6,  x - y = 1 (ranged: 1 <= x - y <= 1),  y - z >= -1,
  // x, y, z free.  Solution follows from minimizing along the free directions.
  const Model model = make_model(
      3, 3, {{0, 0, 1}, {0, 1, 1}, {0, 2, 1}, {1, 0, 1}, {1, 1, -1}, {2, 1, 1}, {2, 2, -1}},
      {1, 2, 3}, {-kInf, -kInf, -kInf}, {kInf, kInf, kInf}, {6, 1, -1}, {6, 1, kInf});
  const LpResult r = solve(model);
  // Substituting x = y + 1, z = 5 - 2y: objective = y + 1 + 2y + 15 - 6y = 16 - 3y, and
  // y - z = 3y - 5 >= -1 gives y >= 4/3; so the LP is unbounded as y grows.
  CHECK(r.status == SimplexStatus::kUnbounded);

  Model bounded = model;
  bounded.col_upper[1] = 3.0;  // y <= 3: optimum at y = 3, objective 16 - 9 = 7.
  const LpResult rb = solve(bounded);
  REQUIRE(rb.status == SimplexStatus::kOptimal);
  CHECK_NEAR(rb.objective, 7.0, 1e-9);
  CHECK_NEAR(rb.col_value[0], 4.0, 1e-9);
  CHECK_NEAR(rb.col_value[2], -1.0, 1e-9);
}

TEST(lp_known_badly_scaled) {
  // The same LP as tiny_lp with rows and columns scaled by powers of ten up to 1e6.
  //   max 3e-3 x' + 2e4 y'  with  x = 1e-3 x', y = 1e4 y'
  //   s.t. 1e5 (x + y) <= 4e5,  1e-4 (x + 3y) <= 6e-4,  x' <= 3e3.
  const Model model = make_model(
      2, 2, {{0, 0, 1e2}, {0, 1, 1e9}, {1, 0, 1e-7}, {1, 1, 3.0}}, {3e-3, 2e4}, {0, 0},
      {3e3, kInf}, {-kInf, -kInf}, {4e5, 6e-4});
  Model max_model = model;
  max_model.sense = samaya::ObjSense::kMaximize;
  const LpResult r = solve(max_model);
  REQUIRE(r.status == SimplexStatus::kOptimal);
  CHECK_NEAR(r.objective, 11.0, 1e-8);
  CHECK_NEAR(r.col_value[0], 3e3, 1e-6);
  CHECK_NEAR(r.col_value[1], 1e-4, 1e-12);
}
