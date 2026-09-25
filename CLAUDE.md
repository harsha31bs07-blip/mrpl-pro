# Instructions for coding agents working in this repository

samaya is a from-scratch LP/MILP solver for SIH problem statement 26119 (MRPL). Read
`docs/ARCHITECTURE.md` before changing code; it defines the contracts and the work packages.

## Hard rules

- C++20, no external solver or linear-algebra libraries (no BLAS/LAPACK, SuiteSparse, Eigen,
  HiGHS, cuSPARSE for factorizations). HiGHS is allowed only in `bench/` as a baseline.
- Every reported result goes through `src/verify/`. Never weaken a verifier tolerance silently.
- Never skip, disable or loosen a test to make it pass. Fix the code or explain why the test is
  wrong.
- Stay inside your work package's files (ownership map in `docs/ARCHITECTURE.md` §6). Edits to
  `src/core/solver.cpp` must be small dispatch hunks.
- Sign conventions: `d = c - Aᵀy`; at a lower bound `d ≥ 0`, at an upper bound `d ≤ 0`; row duals
  behave like reduced costs of row activities (§4.2). Check with `verify_lp_optimality`.

## Before every commit

```sh
cmake --preset debug   && cmake --build --preset debug   && ctest --preset debug    # -Werror
cmake --preset release && cmake --build --preset release && ctest --preset release
cmake --preset asan    && cmake --build --preset asan    && ctest --preset asan
```

All three must pass with zero warnings. Report the exact test output, not a summary of it.

## How to test new algorithms

- Differential tests against `tests/reference_lp.hpp` / `tests/reference_milp.hpp` on random
  families from `tests/lp_generators.hpp` (feasible, infeasible, unbounded, degenerate, boxed).
- Tests use `tests/test_framework.hpp` (`TEST`, `CHECK`, `CHECK_EQ`, `CHECK_NEAR`, `REQUIRE`);
  register new files in `tests/CMakeLists.txt`, sources in the root `CMakeLists.txt`.
- After a test passes, plant a bug (flip a sign, drop a term), confirm the test fails, then
  revert. Mention this check in the commit message.

## Style

Match the surrounding code: 2-space indent, `snake_case`, `member_` suffix, `kConstant` names,
named constants for tolerances, a short comment on every non-obvious numerical choice.
