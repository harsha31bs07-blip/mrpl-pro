// Branch-and-bound against the reference MILP solver (depth-first branch-and-bound on the dense
// reference simplex) on random bounded MILPs, plus limits and known models.

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "core/log.hpp"
#include "mip/branch_and_bound.hpp"
#include "reference_milp.hpp"
#include "samaya/io.hpp"
#include "samaya/solver.hpp"
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

}  // namespace

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
