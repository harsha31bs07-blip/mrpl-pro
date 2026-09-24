# samaya

A sovereign LP / MILP / QP optimization solver core, built from mathematical foundations for
SIH problem statement 26119 (MRPL). No external solver library is used. See [PLAN.md](PLAN.md)
for the architecture, algorithms, benchmarks and timeline.

## Status

| Component | State |
|---|---|
| Model, sparse matrix (CSC), model statistics | done |
| MPS / QPS reader (free + fixed, ranges, all common bound types, integer markers, QUADOBJ/QMATRIX) | done |
| C++ API, C API, `samaya` CLI (`--stats`, `--json`) | done |
| Unit tests, CI (GCC + Clang, Debug / Release / ASan+UBSan) | done |
| Benchmark harness + instance download script | done |
| Scaling, LU factorization, dual simplex | phase 1 |
| Presolve, barrier, PDLP (CPU + GPU) | phase 2–3 |
| Branch-and-cut | phase 4 |

Until the algorithms land, `solve()` validates the model and returns `not_implemented`.

## Build

Requires CMake 3.22+, Ninja and a C++20 compiler (GCC 11+ or Clang 14+).

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

Presets: `debug` (warnings as errors), `release`, `asan` (AddressSanitizer + UBSan), `cuda`
(reserved for the Phase 3 GPU kernels).

## Usage

```sh
build/release/apps/cli/samaya --stats model.mps     # sizes and coefficient ranges
build/release/apps/cli/samaya --json model.mps      # solve, JSON summary on the last line
build/release/apps/cli/samaya --help
```

C++:

```cpp
#include "samaya.hpp"

samaya::Model model = samaya::read_mps("model.mps");
samaya::Result result = samaya::Solver().solve(model);
```

C: see [`include/samaya_c.h`](include/samaya_c.h).

## Benchmarks

```sh
bench/fetch_instances.sh netlib miplib              # downloads into bench/instances/
bench/harness.py bench/instances/netlib --baseline highs --time-limit 300
```

Baseline solvers (for example HiGHS) are run as separate executables for comparison only.

## Layout

```
include/        public API (samaya.hpp, samaya/*.hpp, samaya_c.h)
src/core/       model, solver dispatch, status, logging, C API
src/io/         file readers
src/linalg/     sparse matrices (LU, Cholesky, AMD to follow)
apps/cli/       samaya command-line tool
tests/          unit tests (self-contained framework) and small instances
bench/          benchmark harness and instance download script
```

Later phases add `src/presolve`, `src/lp`, `src/qp`, `src/mip`, `src/gpu` and `src/verify`, as
described in PLAN.md §3.
