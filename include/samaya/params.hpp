#pragma once

#include <cstdint>
#include <string>

#include "samaya/types.hpp"

namespace samaya {

enum class LpMethod : std::uint8_t {
  kAuto,        // Pick per model; concurrent when threads allow.
  kDualSimplex,
  kPrimalSimplex,
  kBarrier,     // Interior point + crossover.
  kPdlp,        // First-order primal-dual method, GPU when available.
  kConcurrent,  // Race dual simplex, barrier and PDLP.
};

struct Params {
  // Limits.
  double time_limit = kInf;  // Seconds of wall-clock time.
  long long node_limit = -1;  // -1 = unlimited.
  int threads = 0;            // 0 = hardware concurrency.

  // Tolerances.
  double primal_feasibility_tol = 1e-7;
  double dual_feasibility_tol = 1e-7;
  double integrality_tol = 1e-6;
  double mip_rel_gap = 1e-4;
  double mip_abs_gap = 1e-6;

  // Algorithm selection.
  LpMethod lp_method = LpMethod::kAuto;
  bool presolve = true;
  bool use_gpu = false;

  // Independently re-check every solution before reporting it.
  bool verify = true;

  // 0 = silent, 1 = summary, 2 = iteration log, 3 = debug.
  int log_level = 1;
};

}  // namespace samaya
