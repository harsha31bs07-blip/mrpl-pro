// Differential tests: the sparse dual simplex against the dense reference simplex on randomly
// generated LPs.

#include <cmath>
#include <cstdio>
#include <random>
#include <string>

#include "core/log.hpp"
#include "lp/lp_solver.hpp"
#include "lp_generators.hpp"
#include "reference_lp.hpp"
#include "samaya/verify.hpp"
#include "test_framework.hpp"

using samaya::LpResult;
using samaya::LpSolveOptions;
using samaya::Model;
using samaya::SimplexStatus;
using samaya::test::LpFamily;
using samaya::test::ReferenceLp;
using samaya::test::ReferenceResult;

namespace {

bool same_status(SimplexStatus a, ReferenceResult::Status b) {
  switch (b) {
    case ReferenceResult::Status::kOptimal: return a == SimplexStatus::kOptimal;
    case ReferenceResult::Status::kInfeasible: return a == SimplexStatus::kInfeasible;
    case ReferenceResult::Status::kUnbounded: return a == SimplexStatus::kUnbounded;
  }
  return false;
}

const char* to_string(ReferenceResult::Status s) {
  switch (s) {
    case ReferenceResult::Status::kOptimal: return "optimal";
    case ReferenceResult::Status::kInfeasible: return "infeasible";
    case ReferenceResult::Status::kUnbounded: return "unbounded";
  }
  return "?";
}

double max_primal_violation(const Model& m, const LpResult& r) {
  double worst = 0.0;
  for (samaya::Index j = 0; j < m.num_cols(); ++j) {
    worst = std::max({worst, m.col_lower[j] - r.col_value[j], r.col_value[j] - m.col_upper[j]});
  }
  for (samaya::Index i = 0; i < m.num_rows(); ++i) {
    worst = std::max(
        {worst, m.row_lower[i] - r.row_activity[i], r.row_activity[i] - m.row_upper[i]});
  }
  return worst;
}

struct Tally {
  int optimal = 0;
  int infeasible = 0;
  int unbounded = 0;
  int mismatches = 0;
  samaya::SimplexStats work;  // Summed over the optimal solves.
};

// Solves `count` random LPs of the family and compares with the reference.
Tally cross_check(LpFamily family, int count, int max_rows, int max_cols, unsigned seed,
                  const LpSolveOptions& options) {
  std::mt19937 rng(seed);
  const samaya::Logger quiet(0);
  Tally tally;
  for (int k = 0; k < count; ++k) {
    const Model model = samaya::test::random_lp(family, max_rows, max_cols, rng);
    const ReferenceResult ref = ReferenceLp(model).solve();
    const LpResult got = samaya::solve_lp(model, options, quiet);
    bool ok = same_status(got.status, ref.status);
    if (ok && ref.status == ReferenceResult::Status::kOptimal) {
      ok = std::fabs(got.objective - ref.objective) <= 1e-6 * (1.0 + std::fabs(ref.objective)) &&
           max_primal_violation(model, got) <= 1e-6;
    }
    // Every outcome must also pass the independent verifier.
    std::string verify_message;
    if (ok) {
      samaya::VerifyReport report;
      if (got.status == SimplexStatus::kOptimal) {
        report = samaya::verify_lp_optimality(model, got.col_value, got.row_dual);
      } else if (got.status == SimplexStatus::kInfeasible) {
        report = samaya::verify_infeasibility(model, got.dual_ray);
      } else {
        report = samaya::verify_unbounded_ray(model, got.primal_ray);
      }
      ok = report.ok;
      verify_message = report.message;
    }
    if (got.status == SimplexStatus::kOptimal) {
      tally.work.dual_iterations += got.stats.dual_iterations;
      tally.work.primal_iterations += got.stats.primal_iterations;
      tally.work.cost_shifts += got.stats.cost_shifts;
      tally.work.rebuilds_on_mismatch += got.stats.rebuilds_on_mismatch;
    }
    switch (ref.status) {
      case ReferenceResult::Status::kOptimal: ++tally.optimal; break;
      case ReferenceResult::Status::kInfeasible: ++tally.infeasible; break;
      case ReferenceResult::Status::kUnbounded: ++tally.unbounded; break;
    }
    if (!ok) {
      ++tally.mismatches;
      std::fprintf(stderr,
                   "  %s #%d (%dx%d): reference %s %.9g, solver %s %.9g (viol %.2e) %s\n",
                   samaya::test::to_string(family), k, model.num_rows(), model.num_cols(),
                   to_string(ref.status), ref.objective, samaya::to_string(got.status),
                   got.objective,
                   got.status == SimplexStatus::kOptimal ? max_primal_violation(model, got) : 0.0,
                   verify_message.c_str());
    }
  }
  std::printf("  %-10s %4d LPs: %d optimal, %d infeasible, %d unbounded, %d mismatches | "
              "optimal solves: %lld dual it, %lld primal it, %lld shifts, %lld mismatch rebuilds\n",
              samaya::test::to_string(family), count, tally.optimal, tally.infeasible,
              tally.unbounded, tally.mismatches, tally.work.dual_iterations,
              tally.work.primal_iterations, tally.work.cost_shifts,
              tally.work.rebuilds_on_mismatch);
  // The fallbacks (cost shifting, primal cleanup) can mask a broken dual simplex, so also check
  // that the dual simplex does the work.
  CHECK(tally.work.primal_iterations <= tally.work.dual_iterations / 10 + 5);
  CHECK(tally.work.cost_shifts <= tally.work.dual_iterations / 5 + 5);
  CHECK_EQ(tally.work.rebuilds_on_mismatch, 0);
  return tally;
}

}  // namespace

TEST(lp_random_default_options) {
  const LpSolveOptions options;
  int optimal = 0, infeasible = 0, unbounded = 0;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kRandom, LpFamily::kDegenerate,
                     LpFamily::kLoose}) {
    const Tally t = cross_check(f, 250, 8, 8, 1000 + static_cast<unsigned>(f), options);
    CHECK_EQ(t.mismatches, 0);
    optimal += t.optimal;
    infeasible += t.infeasible;
    unbounded += t.unbounded;
  }
  // The families must actually exercise every outcome.
  CHECK(optimal > 200);
  CHECK(infeasible > 50);
  CHECK(unbounded > 50);
}

TEST(lp_random_larger_models) {
  const LpSolveOptions options;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kDegenerate}) {
    CHECK_EQ(cross_check(f, 40, 30, 30, 2000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_random_without_scaling_or_perturbation) {
  LpSolveOptions options;
  options.scale = false;
  options.simplex.perturb = false;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kDegenerate, LpFamily::kLoose}) {
    CHECK_EQ(cross_check(f, 150, 8, 8, 3000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_random_frequent_refactorization) {
  // A tiny refactorization interval exercises factorize/update switching on every model.
  LpSolveOptions options;
  options.simplex.refactor_interval = 2;
  for (LpFamily f : {LpFamily::kFeasible, LpFamily::kRandom, LpFamily::kDegenerate}) {
    CHECK_EQ(cross_check(f, 100, 20, 20, 4000 + static_cast<unsigned>(f), options).mismatches, 0);
  }
}

TEST(lp_dual_steepest_edge_weights_stay_exact) {
  // Stop the dual simplex after a few iterations and compare the updated weights with weights
  // recomputed from the factorization.
  std::mt19937 rng(77);
  const samaya::Logger quiet(0);
  int compared = 0;
  int wrong = 0;
  for (int k = 0; k < 60; ++k) {
    const Model model = samaya::test::random_lp(LpFamily::kFeasible, 25, 25, rng);
    for (long long limit : {3LL, 8LL, 20LL}) {
      samaya::LpProblem lp;
      lp.m = model.num_rows();
      lp.n = model.num_cols();
      lp.A = model.A;
      lp.At = model.A.transpose();
      lp.cost = model.obj;
      lp.cost.resize(static_cast<std::size_t>(lp.n + lp.m), 0.0);
      lp.lower = model.col_lower;
      lp.upper = model.col_upper;
      lp.lower.insert(lp.lower.end(), model.row_lower.begin(), model.row_lower.end());
      lp.upper.insert(lp.upper.end(), model.row_upper.begin(), model.row_upper.end());
      samaya::SimplexOptions options;
      options.max_iterations = limit;
      samaya::Simplex simplex(lp, options, quiet);
      if (simplex.solve() != SimplexStatus::kIterationLimit) continue;
      if (simplex.stats().primal_iterations > 0) continue;  // Primal steps reset the weights.
      const std::vector<double> exact = simplex.exact_dse_weights();
      const std::vector<double>& kept = simplex.dse_weights();
      for (std::size_t r = 0; r < exact.size(); ++r) {
        if (exact[r] < 1e-3) continue;  // Updates clamp weights at 1e-4.
        ++compared;
        if (std::fabs(kept[r] - exact[r]) > 1e-6 * exact[r]) ++wrong;
      }
    }
  }
  CHECK(compared > 500);
  CHECK_EQ(wrong, 0);
}
