# MILP implementation plan (our part)

The teammate owns WP1 and WP2 (PDLP on the CPU, then on the GPU). This file covers the rest of
our work: making the MILP solver competitive (Track A), then the MRPL case studies (Track B).
It also says where each benchmark runs: in the cloud session or on the Dell G15.

## Where we start

The baseline is `main` at ba0046e. Full numbers are in [results/miplib.md](results/miplib.md).

| Measure (62 small "easy" MIPLIB 2017 instances) | Now |
|---|---|
| Solved in 10 min (first cut version) | 8/62 |
| Solved in 60 s | 3/62 |
| Feasible solution found in 60 s | 37/62 |
| Median gap to the known optimum at 60 s | 2.8% |
| Wrong answers | 0 |

**Why instances fail today:**
- **No feasible solution.** On 25/62 instances no feasible solution is found in 60 s, for
  example b1c1s1, cost266-UUE, glass4, lotsize, timtab1 and the rococo instances.
  - Cause: the only heuristics are simple rounding and round-and-solve.
- **Weak bounds.** Many instances have a good solution but a bound that doesn't move: mas74,
  p200x1188c, gen-ip*, pg5_34, neos5.
  - Cause: no flow covers, no probing or clique information, and cuts only at the root.
- **One core only.** All searches are single-threaded; the laptop has 16 hardware threads.

## Targets

| Measure | Now | Target |
|---|---|---|
| Solved in 10 min | 8/62 | **≥ 20/62** |
| Feasible solution in 60 s | 37/62 | **≥ 52/62** |
| Median gap at 60 s | 2.8% | **≤ 1%** |
| Wrong answers | 0 | **0**, always |

All targets are for a single thread. The parallel search (A6) is reported separately.

## Rules for every milestone

- **Tests.** Each milestone adds tests to `tests/test_mip.cpp`, or to a new
  `tests/test_mip_*.cpp`:
  - a differential check against `tests/reference_milp.hpp` on the random families;
  - a debug-solution check: no reduction, cut or heuristic may cut off the known optimum;
  - a planted-bug check, mentioned in the commit message.
- **Build.** The debug, release and asan presets must be green with zero warnings before every
  push.
- **Screening.** A 60 s run of the 62 instances, the new build against the previous one on the
  same machine. A change is kept only if:
  - there are no wrong answers;
  - solved and feasible counts do not drop;
  - the shifted geomean does not get worse by more than 3% (the noise level).
- **One PR per milestone.**
- **Full runs.** A 10-minute run at A0, after A3 and at the end (A7). Results go into
  `docs/results/miplib.md`, with the machine named.
- **Same machine.** Numbers from the laptop and from the cloud are never compared with each
  other. Every comparison reruns both builds on the same machine.

## Track A: MILP

Ordered by expected gain per unit of work.

### A0. Baseline with the tuned cuts (benchmark only)
- 10-minute run of the current `main` (it has not been run with the tuned cuts).
- If possible, run HiGHS single-threaded on the same machine so we know the gap to close.
- **Decides:** the reference numbers for everything below.

### A1. Small fixes in the search (`src/mip/branch_and_bound.cpp`)
- **Time limit:**
  - check the clock inside strong branching and the cut loop;
  - give the node LP an iteration or time budget so a run stops within 1 s of the limit
    (it can overshoot by up to 8 s today).
- **Root restart:** if the root cuts and propagation fix more than 20% of the integers,
  presolve again and restart once.
- **Unit tests** for the time-limit contract.

### A2. MIP presolve (`src/presolve/`, MIP mode)
- **Probing on binaries:** fix each to 0 and 1, then propagate.
  - A side that is infeasible fixes the column.
  - Collect implications and implied bounds.
  - Detect equal and complementary binaries and merge them.
  - A work limit keeps the cost down.
- **Coefficient tightening** on ≤ rows with binaries (a_j := min(a_j, rhs − minact) style).
- **Clique detection:** set-packing rows plus probing implications give a clique table, shared
  with A4 and propagation.
- **Implied integers:** a continuous column that is forced to be integral is marked as such.
- **Dual reductions for integer columns** and dominated-column fixing.
- **Tests:**
  - postsolve checks with the verifier;
  - debug-solution checks;
  - planted reductions on random MILPs (as in `test_presolve.cpp`).
- **Expected to help:** p200x1188c, gen-ip*, mas*, neos5, set-partitioning-style instances.

### A3. Primal heuristics (`src/mip/heuristics.{hpp,cpp}`, a new file)
- **Diving:** fractional, coefficient (lock-based), guided (toward the incumbent) and
  pseudocost diving.
  - Diving runs from the root and periodically in the tree, using the node LP with a node
    budget.
- **Feasibility pump** at the root, with perturbation and restarts.
  - This targets the instances where we find no solution at all.
- **Sub-MIP heuristics:**
  - RENS at the root: fix integers that are already integral, round the domains of the rest;
  - RINS in the tree: fix where the incumbent and the LP solution agree.
  - Both are solved with our own `BranchAndBound` on the reduced model, with a node and time
    limit. This needs A1's limits to be exact.
- **Shift-and-propagate** as a cheap LP-free start heuristic.
- **Scheduling:** each heuristic tracks its success rate and the time it has spent; the budget
  is about 10% of the solve time.
- **10-minute benchmark checkpoint.**

### A4. More cuts, and cuts in the tree (`src/mip/cuts.cpp`)
- **Flow covers** on variable-upper-bound structures (x ≤ u·y). These are the p200x1188c,
  fixed-charge and network-design instances.
- **Clique cuts** from the A2 clique table, and **implied-bound cuts** from probing.
- **Zero-half cuts** (simple version) for pure binary instances.
- **Cut pool:**
  - cuts are stored with their age;
  - they are separated again at nodes up to depth 5 or so, and every k nodes deeper;
  - local cuts stay valid only in their subtree.
- **Tests:** the debug-solution check on every separator, and the planted-bug check.

### A5. Tree search (`src/mip/branch_and_bound.cpp`)
- **Conflict analysis:** when a node is infeasible after propagation, learn a clause and
  propagate it.
- **Node LP speed:**
  - reuse the factorization when bounds change only by a few columns;
  - skip strong branching deep in the tree when the pseudocosts are reliable.
- **Symmetry** is out of scope unless A0–A4 leave time. It needs orbit detection, which is a
  large piece of work.

### A6. Parallel tree search (`src/mip/parallel*`, std::thread only)
- **Structure:**
  - each worker owns a copy of the scaled LP and a Simplex;
  - the open-node queue, incumbent and pseudocosts are shared, behind a mutex, in batches;
  - the root runs once, and workers start from its basis.
- **Deterministic mode** (fixed synchronization points) for reproducible tests; an
  opportunistic mode for speed.
- **Option:** `--threads N`, default 1, so single-thread results stay comparable.
- **Target:** at least 2.5× speed-up at 8 threads on the solved instances.
  - This is measured on the laptop, since the cloud machine has only 4 cores.

### A7. Final numbers and documentation
- **10-minute runs**, single-threaded and 8 threads, next to HiGHS on the same machine.
- **Update:**
  - `docs/results/miplib.md`;
  - the status table in `docs/ARCHITECTURE.md`;
  - a short "how the MILP solver works" section for the judges.

## Track B: MRPL case studies (WP5)

This starts after A3, when the solver finds solutions reliably, or earlier in parallel if
there is time.
- **Model families** (defined in ARCHITECTURE.md WP5):
  - multi-period planning (LP);
  - crude scheduling with cargo windows (MILP);
  - blending (LP);
  - a utility-system unit commitment (MILP).
- **Sizes:** three per family. Each size is cross-checked against HiGHS.
- **Demo:** `cases/demo.sh` generates, solves, verifies and prints a readable plan in under
  2 minutes.

Barrier (WP3) and QP (WP4) are not in this plan. We decide after A4 whether there is time for
them. They are large, and PDLP (the teammate's WP1) already gives us a second LP method.

## Where benchmarks run: cloud or laptop

**Before every benchmark I'll ask you "local or cloud?"**, with these estimates:

| Run | Cloud (4 cores, 4 at a time) | Dell G15 (5 at a time) |
|---|---|---|
| 60 s screening, one build | ≈ 17 min | ≈ 13 min, and faster per core |
| 60 s screening, A/B builds | ≈ 35 min | ≈ 26 min |
| 10 min full run, one build | ≈ 2.6 h | ≈ 2.1 h |
| HiGHS baseline, 10 min | ≈ 2.6 h | ≈ 2.1 h |
| Parallel search (A6) scaling | not possible (4 cores) | **laptop only** |

**When the laptop is worth it:**
- **Always for full 10-minute runs and A6.** While the laptop benchmarks, I keep coding in the
  cloud. In the cloud, a benchmark and my builds share 4 cores, which makes the timings noisy.
- **The cloud is fine for 60 s screening** between small steps, so you don't have to be at the
  laptop.
- **Memory:** 16 GB on the laptop, and one run can hold up to about 2 GB of open nodes. Keep it
  to 5 runs at a time.

### Running a benchmark on the laptop

On Linux or WSL2, from a clone of the repository:

```sh
git fetch origin && git checkout <branch I name> && git pull
cmake --preset release && cmake --build --preset release
bench/fetch_instances.sh miplib-list bench/miplib_small.test     # first time only
bench/run_miplib.sh 5 600 bench/miplib_small.test build/release/apps/cli/samaya \
  bench/results/<run-name>
```

Send the results back on a branch, and I'll fetch and write up the report:

```sh
git checkout -b results/<run-name>
git add -f bench/results/<run-name>/*.csv && git commit -m "Benchmark results: <run-name>"
git push -u origin results/<run-name>
```

**Keep the laptop in its performance power mode, plugged in, and otherwise idle during the
run.**

For the HiGHS baseline, first run `pip install highspy`. Then add `--baseline highspy` to a
`bench/harness.py` call; I'll give you the exact command when we get there.
