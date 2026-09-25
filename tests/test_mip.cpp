// Branch-and-bound against the reference MILP solver (depth-first branch-and-bound on the dense
// reference simplex) on random bounded MILPs, plus limits and known models.

#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "core/log.hpp"
#include "mip/branch_and_bound.hpp"
#include "reference_milp.hpp"
#include "samaya/io.hpp"
#include "samaya/solver.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::Index;
using samaya::kInf;
using samaya::Model;
using samaya::Status;
using samaya::test::ReferenceMilpResult;

namespace {

enum class MilpFamily { kMixed, kPureInteger, kKnapsack, kEquality };

// Bounded random MILP. Rows are built around a point that is integral in the integer columns, so
// most models are feasible; kEquality adds equality rows over integer columns (often infeasible).
Model random_milp(MilpFamily family, std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const auto chance = [&](double p) { return std::bernoulli_distribution(p)(rng); };
  const auto real = [&](double lo, double hi) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
  };
  Model model;
  model.sense = chance(0.4) ? samaya::ObjSense::kMaximize : samaya::ObjSense::kMinimize;
  const int n = family == MilpFamily::kKnapsack ? uniform_int(4, 14) : uniform_int(2, 10);
  std::vector<double> point(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    const bool integer = family != MilpFamily::kMixed || chance(0.6);
    const int lo = family == MilpFamily::kKnapsack ? 0 : uniform_int(-3, 2);
    const int span = family == MilpFamily::kKnapsack ? uniform_int(1, 3) : uniform_int(1, 6);
    model.col_lower.push_back(lo);
    model.col_upper.push_back(lo + span);
    model.col_type.push_back(integer ? samaya::VarType::kInteger : samaya::VarType::kContinuous);
    point[j] = integer ? uniform_int(lo, lo + span) : real(lo, lo + span);
    model.obj.push_back(family == MilpFamily::kKnapsack ? uniform_int(1, 20)
                                                        : uniform_int(-9, 9) * real(0.5, 1.5));
  }
  if (family == MilpFamily::kKnapsack) model.sense = samaya::ObjSense::kMaximize;

  std::vector<samaya::Triplet> t;
  const int m = family == MilpFamily::kKnapsack ? uniform_int(1, 3) : uniform_int(1, 8);
  for (int i = 0; i < m; ++i) {
    double act = 0.0;
    double weight_sum = 0.0;
    for (int j = 0; j < n; ++j) {
      if (family != MilpFamily::kKnapsack && !chance(0.6)) continue;
      double a = family == MilpFamily::kKnapsack ? uniform_int(1, 25) : uniform_int(-6, 6);
      if (a == 0.0) a = 3.0;
      if (family == MilpFamily::kMixed && chance(0.3)) a *= real(0.5, 1.5);
      t.push_back({i, j, a});
      act += a * point[j];
      weight_sum += a * model.col_upper[j];
    }
    double lo = -kInf;
    double up = kInf;
    if (family == MilpFamily::kKnapsack) {
      up = std::floor(weight_sum * real(0.3, 0.7));
    } else if (family == MilpFamily::kEquality && chance(0.5)) {
      // Shifted right-hand sides may leave no integer solution even when the LP is feasible.
      const double r = real(0, 1);
      lo = up = act + (r < 0.3 ? uniform_int(-2, 2) : r < 0.5 ? real(0.1, 0.9) : 0.0);
    } else {
      const double slack = chance(0.5) ? 0.0 : real(0, 4);
      const double r = real(0, 1);
      if (r < 0.4) {
        up = act + slack;
      } else if (r < 0.8) {
        lo = act - slack;
      } else {
        lo = act - slack;
        up = act + slack + 1;
      }
      // Fractional right-hand sides make the LP bound weaker than the integer optimum.
      if (chance(0.5) && up < kInf) up += real(0, 0.9);
    }
    model.row_lower.push_back(lo);
    model.row_upper.push_back(up);
  }
  model.A = samaya::SparseMatrix::from_triplets(m, n, std::move(t));
  return model;
}

const char* name(MilpFamily f) {
  switch (f) {
    case MilpFamily::kMixed: return "mixed";
    case MilpFamily::kPureInteger: return "pure";
    case MilpFamily::kKnapsack: return "knapsack";
    case MilpFamily::kEquality: return "equality";
  }
  return "?";
}

struct MilpTally {
  int optimal = 0;
  int infeasible = 0;
  int failures = 0;
  long long nodes = 0;
};

MilpTally cross_check(MilpFamily family, int count, unsigned seed, bool presolve) {
  std::mt19937 rng(seed);
  MilpTally tally;
  samaya::Params params;
  params.log_level = 0;
  params.presolve = presolve;
  params.mip_rel_gap = 0.0;
  params.mip_abs_gap = 1e-9;
  for (int k = 0; k < count; ++k) {
    const Model model = random_milp(family, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status == ReferenceMilpResult::Status::kNodeLimit) continue;
    const samaya::Result got = samaya::Solver(params).solve(model);
    bool ok = false;
    if (ref.status == ReferenceMilpResult::Status::kOptimal) {
      ++tally.optimal;
      ok = got.status == Status::kOptimal && got.verified &&
           std::fabs(got.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)) &&
           std::fabs(got.dual_bound - got.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective));
    } else if (ref.status == ReferenceMilpResult::Status::kInfeasible) {
      ++tally.infeasible;
      ok = got.status == Status::kInfeasible;
    }
    tally.nodes += got.nodes;
    if (!ok) {
      ++tally.failures;
      std::fprintf(stderr, "  %s #%d (%dx%d): reference %d %.9g, solver %s %.9g bound %.9g\n",
                   name(family), k, model.num_rows(), model.num_cols(),
                   static_cast<int>(ref.status), ref.objective, samaya::to_string(got.status),
                   got.objective, got.dual_bound);
    }
  }
  std::printf("  %-9s presolve %d: %d optimal, %d infeasible, %d failures, %lld nodes\n",
              name(family), presolve, tally.optimal, tally.infeasible, tally.failures,
              tally.nodes);
  return tally;
}

// Set partitioning with a planted partition: every row is covered exactly once. Rounding the LP
// point rarely gives a partition, so these need the pump, diving or a sub-MIP.
Model set_partitioning(std::mt19937& rng) {
  const auto uniform_int = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  const int m = uniform_int(15, 30);
  std::vector<std::vector<int>> columns;
  // The planted partition: consecutive blocks of 1-3 rows.
  for (int i = 0; i < m;) {
    const int size = std::min(m - i, uniform_int(1, 3));
    std::vector<int> col;
    for (int k = 0; k < size; ++k) col.push_back(i + k);
    columns.push_back(col);
    i += size;
  }
  const int extra = uniform_int(2 * m, 4 * m);
  for (int k = 0; k < extra; ++k) {
    std::vector<int> col;
    const int size = uniform_int(2, 4);
    while (static_cast<int>(col.size()) < size) {
      const int r = uniform_int(0, m - 1);
      if (std::find(col.begin(), col.end(), r) == col.end()) col.push_back(r);
    }
    columns.push_back(col);
  }
  std::shuffle(columns.begin(), columns.end(), rng);
  Model model;
  std::vector<samaya::Triplet> t;
  for (std::size_t j = 0; j < columns.size(); ++j) {
    // Costs slightly below the size favour the random (non-planted) columns in the LP.
    model.obj.push_back(static_cast<double>(columns[j].size()) *
                        (0.8 + 0.4 * (uniform_int(0, 99) / 100.0)));
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
    for (const int r : columns[j]) t.push_back({r, static_cast<Index>(j), 1.0});
  }
  model.row_lower.assign(static_cast<std::size_t>(m), 1.0);
  model.row_upper.assign(static_cast<std::size_t>(m), 1.0);
  model.A = samaya::SparseMatrix::from_triplets(m, static_cast<Index>(columns.size()),
                                                std::move(t));
  return model;
}

}  // namespace

TEST(mip_heuristics_find_solutions_at_the_root) {
  // The root alone (one node, no cuts): simple rounding and round-and-solve against the pump,
  // diving and RENS. Every solution must pass the verifier on the model.
  const samaya::Logger quiet(0);
  int found[2] = {0, 0};
  int models = 0;
  std::mt19937 rng(31);
  for (int k = 0; k < 40; ++k) {
    const Model model = set_partitioning(rng);
    ++models;
    for (int h = 0; h < 2; ++h) {
      samaya::MipOptions options;
      options.cuts = false;
      options.node_limit = 1;
      options.heuristics = h == 1;
      const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
      if (out.x.empty()) continue;
      ++found[h];
      const samaya::VerifyReport report = samaya::verify_primal(model, out.x, {});
      CHECK(report.ok);
      CHECK(out.bound <= out.objective + 1e-9);
    }
  }
  std::printf("  %d set-partitioning models: solution at the root %d without heuristics, %d "
              "with\n", models, found[0], found[1]);
  CHECK(found[1] >= 36);
  CHECK(found[1] > found[0]);
}

TEST(mip_objective_cutoff_accepts_only_better_solutions) {
  const samaya::Logger quiet(0);
  std::mt19937 rng(41);
  int checked = 0;
  for (int k = 0; k < 120; ++k) {
    const Model model = random_milp(k % 2 ? MilpFamily::kMixed : MilpFamily::kKnapsack, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    ++checked;
    const double worse = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
    samaya::MipOptions options;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    // Nothing is better than the optimum (the margin covers round-off between the reference's
    // objective and ours).
    options.objective_cutoff = ref.objective - worse * 1e-6 * (1.0 + std::fabs(ref.objective));
    CHECK(samaya::BranchAndBound(model, options, quiet).solve().x.empty());
    // A looser cutoff still finds the optimum.
    options.objective_cutoff = ref.objective + worse * (1.0 + std::fabs(ref.objective));
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    CHECK(out.status == Status::kOptimal);
    CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)));
  }
  CHECK(checked > 60);
}

TEST(mip_random_models_match_reference) {
  for (const bool presolve : {true, false}) {
    const MilpTally mixed = cross_check(MilpFamily::kMixed, 250, 11, presolve);
    const MilpTally pure = cross_check(MilpFamily::kPureInteger, 250, 12, presolve);
    const MilpTally knapsack = cross_check(MilpFamily::kKnapsack, 200, 13, presolve);
    const MilpTally equality = cross_check(MilpFamily::kEquality, 250, 14, presolve);
    CHECK_EQ(mixed.failures + pure.failures + knapsack.failures + equality.failures, 0);
    CHECK(mixed.optimal > 150);
    CHECK(pure.optimal > 150);
    CHECK(knapsack.optimal > 150);
    CHECK(equality.infeasible > 20);
    CHECK(knapsack.nodes > 200);  // The search actually branches.
  }
}

TEST(mip_cuts_never_separate_an_optimal_solution) {
  // Root cuts are checked against a known optimal solution (the reference's); a valid cut can
  // never separate it. The search runs without presolve so the solution applies directly.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long cuts = 0;
  int models = 0;
  int tightened = 0;
  for (const MilpFamily family :
       {MilpFamily::kMixed, MilpFamily::kPureInteger, MilpFamily::kKnapsack}) {
    std::mt19937 rng(100 + static_cast<unsigned>(family));
    for (int k = 0; k < 200; ++k) {
      const Model model = random_milp(family, rng);
      const ReferenceMilpResult ref = samaya::test::reference_milp(model);
      if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
      samaya::MipOptions options;
      options.rel_gap = 0.0;
      options.abs_gap = 1e-9;
      options.debug_solution = ref.x;
      samaya::BranchAndBound search(model, options, quiet);
      const samaya::MipOutcome out = search.solve();
      ++models;
      violations += out.debug_cut_violations;
      cuts += out.cut_rounds > 0 ? out.cuts_added : 0;
      const double sense = model.sense == samaya::ObjSense::kMaximize ? -1.0 : 1.0;
      if (sense * out.root_bound_cuts > sense * out.root_bound + 1e-9) ++tightened;
      CHECK(out.status == Status::kOptimal);
      CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    }
  }
  std::printf("  %d models, %lld cuts kept, %d root bounds tightened, %lld violations\n",
              models, cuts, tightened, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 400);
  CHECK(tightened > 100);
}

TEST(mip_root_reductions_keep_an_optimal_solution) {
  // Coefficient tightening and probing change the root bounds and rows; a known optimal
  // solution must stay feasible for them, and the search must still find the optimum.
  const samaya::Logger quiet(0);
  long long violations = 0;
  long long reductions = 0;
  int models = 0;
  for (const MilpFamily family :
       {MilpFamily::kMixed, MilpFamily::kPureInteger, MilpFamily::kKnapsack}) {
    std::mt19937 rng(200 + static_cast<unsigned>(family));
    for (int k = 0; k < 200; ++k) {
      Model model = random_milp(family, rng);
      // Binaries make the reductions apply more often.
      if (k % 2 == 0) {
        for (Index j = 0; j < model.num_cols(); ++j) {
          if (model.col_type[j] == samaya::VarType::kInteger) {
            model.col_lower[j] = 0;
            model.col_upper[j] = 1;
          }
        }
      }
      const ReferenceMilpResult ref = samaya::test::reference_milp(model);
      if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
      samaya::MipOptions options;
      options.rel_gap = 0.0;
      options.abs_gap = 1e-9;
      options.cuts = false;
      options.debug_solution = ref.x;
      const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
      ++models;
      violations += out.debug_reduction_violations;
      reductions += out.coefficients_tightened + out.probing_fixed + out.probing_tightened;
      CHECK(out.status == Status::kOptimal);
      CHECK(std::fabs(out.objective - ref.objective) <= 1e-6 * (1 + std::fabs(ref.objective)));
    }
  }
  // 5 x1 + x2 + x3 <= 6 over binaries is x1 + x2 + x3 <= 2 on integer points.
  Model knap;
  knap.sense = samaya::ObjSense::kMaximize;
  knap.obj = {5, 3, 3};  // LP bound 10 on the original row (x2 = x3 = 1, x1 = 0.8), 8 after.
  knap.col_lower = {0, 0, 0};
  knap.col_upper = {1, 1, 1};
  knap.col_type.assign(3, samaya::VarType::kInteger);
  knap.row_lower = {-kInf};
  knap.row_upper = {6};
  knap.A = samaya::SparseMatrix::from_triplets(1, 3, {{0, 0, 5}, {0, 1, 1}, {0, 2, 1}});
  samaya::MipOptions options;
  options.debug_solution = {1, 1, 0};
  const samaya::MipOutcome out = samaya::BranchAndBound(knap, options, quiet).solve();
  CHECK(out.coefficients_tightened >= 1);
  CHECK_EQ(out.debug_reduction_violations, 0);
  CHECK_NEAR(out.objective, 8.0, 1e-9);
  CHECK_NEAR(out.root_bound, 8.0, 1e-9);  // The tightened row makes the LP integral.

  std::printf("  %d models, %lld reductions, %lld violations\n", models, reductions, violations);
  CHECK_EQ(violations, 0);
  CHECK(models > 400);
  CHECK(reductions > 50);
}

TEST(mip_random_models_match_reference_without_cuts) {
  std::mt19937 rng(21);
  samaya::Params params;
  params.log_level = 0;
  int failures = 0;
  for (int k = 0; k < 150; ++k) {
    const Model model = random_milp(k % 2 ? MilpFamily::kMixed : MilpFamily::kKnapsack, rng);
    const ReferenceMilpResult ref = samaya::test::reference_milp(model);
    if (ref.status != ReferenceMilpResult::Status::kOptimal) continue;
    samaya::MipOptions options;
    options.cuts = false;
    options.rel_gap = 0.0;
    options.abs_gap = 1e-9;
    const samaya::Logger quiet(0);
    const samaya::MipOutcome out = samaya::BranchAndBound(model, options, quiet).solve();
    if (out.status != Status::kOptimal ||
        std::fabs(out.objective - ref.objective) > 1e-6 * (1 + std::fabs(ref.objective))) {
      ++failures;
    }
  }
  CHECK_EQ(failures, 0);
}

TEST(mip_node_limit_reports_limit_with_valid_bound) {
  // A knapsack whose LP bound is loose needs more than one node.
  Model model;
  model.sense = samaya::ObjSense::kMaximize;
  const double w[] = {12, 7, 11, 8, 9, 6, 5, 13, 10, 4};
  const double v[] = {24, 13, 23, 15, 16, 11, 10, 25, 19, 7};
  std::vector<samaya::Triplet> t;
  for (int j = 0; j < 10; ++j) {
    model.obj.push_back(v[j]);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
    t.push_back({0, j, w[j]});
  }
  model.row_lower.push_back(-kInf);
  model.row_upper.push_back(26);
  model.A = samaya::SparseMatrix::from_triplets(1, 10, std::move(t));
  const ReferenceMilpResult ref = samaya::test::reference_milp(model);
  REQUIRE(ref.status == ReferenceMilpResult::Status::kOptimal);

  samaya::Params params;
  params.log_level = 0;
  params.presolve = false;
  params.mip_rel_gap = 0.0;
  params.node_limit = 1;
  const samaya::Result limited = samaya::Solver(params).solve(model);
  CHECK(limited.status == Status::kNodeLimit || limited.status == Status::kOptimal);
  CHECK(limited.dual_bound >= ref.objective - 1e-9);  // Maximization: bound from above.

  // An open-node cap stops the search cleanly with a valid bound instead of exhausting memory;
  // the soft cap switches to depth-first search first.
  const samaya::Logger quiet(0);
  samaya::MipOptions capped;
  capped.rel_gap = 0.0;
  capped.cuts = false;
  capped.max_open_nodes_soft = 2;
  capped.max_open_nodes = 3;
  const samaya::MipOutcome out = samaya::BranchAndBound(model, capped, quiet).solve();
  CHECK(out.status == Status::kNodeLimit || out.status == Status::kOptimal);
  CHECK(out.bound >= ref.objective - 1e-9);
  capped.max_open_nodes = 1000;
  const samaya::MipOutcome dfs = samaya::BranchAndBound(model, capped, quiet).solve();
  CHECK(dfs.status == Status::kOptimal);
  CHECK_NEAR(dfs.objective, ref.objective, 1e-9);

  params.node_limit = -1;
  const samaya::Result full = samaya::Solver(params).solve(model);
  CHECK(full.status == Status::kOptimal);
  CHECK(full.verified);
  CHECK_NEAR(full.objective, ref.objective, 1e-9);
}

TEST(mip_time_limit_is_respected) {
  // A market-share model (equality rows over binaries with slack penalties) is far too hard to
  // finish in the limit; the search must stop within it and report a valid bound.
  std::mt19937 rng(7);
  constexpr int kRows = 4;
  constexpr int kBinaries = 40;
  Model model;
  std::vector<samaya::Triplet> t;
  for (int j = 0; j < kBinaries; ++j) {
    model.obj.push_back(0.0);
    model.col_lower.push_back(0);
    model.col_upper.push_back(1);
    model.col_type.push_back(samaya::VarType::kInteger);
  }
  for (int i = 0; i < kRows; ++i) {
    double sum = 0.0;
    for (int j = 0; j < kBinaries; ++j) {
      const double a = std::uniform_int_distribution<int>(0, 99)(rng);
      t.push_back({i, j, a});
      sum += a;
    }
    // sum_j a_ij x_j + s_i^- - s_i^+ = floor(sum / 2), minimizing the slacks.
    for (int k = 0; k < 2; ++k) {
      t.push_back({i, static_cast<Index>(model.obj.size()), k == 0 ? 1.0 : -1.0});
      model.obj.push_back(1.0);
      model.col_lower.push_back(0);
      model.col_upper.push_back(kInf);
      model.col_type.push_back(samaya::VarType::kContinuous);
    }
    model.row_lower.push_back(std::floor(sum / 2));
    model.row_upper.push_back(std::floor(sum / 2));
  }
  model.A = samaya::SparseMatrix::from_triplets(kRows, static_cast<Index>(model.obj.size()),
                                                std::move(t));
  samaya::Params params;
  params.log_level = 0;
  params.time_limit = 0.5;
  params.node_limit = 5000000;  // A backstop so a search that ignores the clock still ends.
  const auto start = std::chrono::steady_clock::now();
  const samaya::Result result = samaya::Solver(params).solve(model);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  CHECK(result.status == Status::kTimeLimit || result.status == Status::kOptimal);
  // Generous slack for sanitizer builds; the overshoot this guards against was seconds.
  CHECK(elapsed <= params.time_limit + 0.5);
  if (result.status == Status::kTimeLimit && std::isfinite(result.objective)) {
    CHECK(result.verified);
    CHECK(result.dual_bound <= result.objective + 1e-9);
  }
}

TEST(mip_nearly_integral_lp_solution_is_repaired_not_pruned) {
  // Regression (cases/mrpl.py crude, small): the root LP after cuts is integral within the
  // tolerance, but rounding a cargo binary (coefficient 130 in a tank balance) breaks the row by
  // more than the row tolerance. The node must be repaired by re-solving the continuous columns,
  // not pruned as infeasible. HiGHS: optimal 708749.0576.
  const Model model =
      samaya::read_mps(std::string(SAMAYA_TEST_DATA_DIR) + "/mrpl_crude_small.mps");
  samaya::Params params;
  params.log_level = 0;
  const samaya::Result result = samaya::Solver(params).solve(model);
  CHECK(result.status == Status::kOptimal);
  CHECK(result.verified);
  CHECK(std::fabs(result.objective - 708749.0576) <= 1e-4 * 708749.0576);
}

TEST(mip_infeasible_and_unbounded_models) {
  // 2x = 1 with x integer: LP feasible, integer infeasible.
  Model infeasible = samaya::read_mps_from_string(
      "NAME X\nROWS\n N obj\n E r\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n x obj 1 r 2\n"
      " MARKER 'MARKER' 'INTEND'\nRHS\n rhs r 1\nBOUNDS\n UP b x 10\nENDATA\n");
  samaya::Params params;
  params.log_level = 0;
  CHECK(samaya::Solver(params).solve(infeasible).status == Status::kInfeasible);
  params.presolve = false;
  CHECK(samaya::Solver(params).solve(infeasible).status == Status::kInfeasible);

  // max x + y with x integer free: unbounded relaxation.
  Model unbounded = samaya::read_mps_from_string(
      "NAME X\nOBJSENSE\n MAX\nROWS\n N obj\n L r\nCOLUMNS\n MARKER 'MARKER' 'INTORG'\n"
      " x obj 1 r 1\n MARKER 'MARKER' 'INTEND'\n y obj 1 r -1\nRHS\n rhs r 3\nBOUNDS\n"
      " PL b x\nENDATA\n");
  const samaya::Result r = samaya::Solver(params).solve(unbounded);
  CHECK(r.status == Status::kUnbounded || r.status == Status::kInfeasibleOrUnbounded);
}
