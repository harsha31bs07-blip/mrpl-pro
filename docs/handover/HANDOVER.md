# Handover: MILP work (from the cloud session, 2026-09-26)

Read `CLAUDE.md` first; its rules hold without exception. This note says where the MILP work
stands and how to continue it.

## Working rules agreed with the project owner

- **Stay on the fundamentals; never tune to the test set.**
  - Fixed instance sets, chosen before any run, and every result reported, bad ones included.
  - No growing or cherry-picking instances until the others look slow.
  - Parameters are chosen from the literature or from principle, and any tuning is stated.
- **Never weaken tests or verifier tolerances.** A new feature needs a differential test
  against `tests/reference_milp.hpp` and a planted bug that the test catches; say both in the
  commit message.
- **Before every commit:** the debug, release and asan presets, zero warnings, with the exact
  test output quoted.
- **Benchmarks:**
  - Never edit a script while a benchmark is running it.
  - Never touch `results/*` branches; the owner pushes laptop results there.
  - Don't build or test on a machine while it runs a timing benchmark.
- **No PR or merge without asking the owner.**
- **Kill any wait loops you leave behind.**

Branch: `claude/kind-hamilton-jr2b6u`. The goal is to be better than HiGHS, the leading
open-source solver, by honest means.

## What is in the branch

Recent MILP work, each item with tests and planted-bug checks (details in `docs/MILP_PLAN.md`):

- **Feasibility Jump** (`src/mip/feasibility_jump.cpp`): LP-free local search before the root
  cuts. The root still runs the objective pump while FJ's point is the incumbent.
- **Dive fix:** dives stop on columns already within the tolerance of an integral bound.
- **Cuts in the tree:** a pool of root cuts, plus fresh c-MIR/cover cuts every 10th depth
  derived with the root bounds. A node keeps its cuts only if they close 1% of its gap
  (they're removed exactly otherwise), and kept cuts may grow the LP by at most 25% of its rows.
  Bound A/B: 23 instances better, 9 worse. Largest gains neos17, enlight_hard, binkar10_1;
  largest loss neos-911970.
- **MIP start:**
  - `--mip-start FILE` in the CLI and `Params::mip_start` in the API.
  - Re-planning variants: `cases/mrpl.py generate --update SEED`.
  - Measured gain on our cases is small; see `cases/README.md`.
- **Verification write-up:** `docs/results/comparison.md`, last section.
- **Bound A/B tools:** `bench/ab_bound.sh` and `bench/ab_report.py`.

## Work in progress: conflict analysis (`docs/handover/conflict-analysis.patch`)

Apply it with `git apply docs/handover/conflict-analysis.patch`. It adds `src/mip/conflicts.cpp`
and a test.

**How it works:**
- An infeasible node LP gives a proof from its Farkas ray: (A'y)'x must lie within
  [min y's, max y's] over the row bounds.
- A node LP cut off by the incumbent gives a proof from its duals: d'x <= cutoff - offset -
  min y's, with d = c - A'y. The right-hand side follows the current cutoff.
- Both are valid for any multipliers, since they use only the rows and their global bounds.
- Proofs sparser than 15% of the columns + 10 are kept (at most 1000, the least recently
  used replaced) and propagated at the start of every node.

**Status:**
- **Test** `mip_conflicts_never_exclude_an_optimal_solution`:
  - 98 models, 556 proofs, 3746 nodes pruned by them, 0 violations, every optimum matches the
    reference.
  - A planted sign bug in the cutoff proof is caught (wrong optimum and a violated proof).
  - **Open:** a planted bug that drops the row scaling of the Farkas multipliers
    (`y[i] = ray[i]` instead of `ray[i] * scaling_.row[i]`) is **not** caught yet. The cause is
    known: in the test's families `analyze_infeasible_lp` is never called, because bound
    propagation prunes every infeasible node before its LP. All 556 proofs are cut-off proofs.
    Add a family whose node LPs become infeasible only through several rows together, e.g.
    equality rows that share continuous columns, keep the 1e-3..1e3 row factors, and require
    some infeasible-LP proofs in the test (count them in MipOutcome). Then the scaling planted
    bug must fail.
- **Quick trial** (60 s, against the tree-cuts build):
  - neos17: solved in 19.9 s. Before: unsolved; HiGHS and SCIP solve it.
  - timtab1: bound 499,978 -> 524,215, and a solution found.
  - neos-911970: bound 50.34 -> 50.89.
  - mas76 and pk1: unchanged.
- **Still to do:**
  1. Resolve the open test point.
  2. The three presets.
  3. A full bound A/B:
     `bench/ab_bound.sh OUT OLD NEW 60 && python3 bench/ab_report.py OUT`.
  4. Commit with the numbers.

## Next tasks, in the order agreed

1. **Finish conflict analysis** (above). Then consider proofs from infeasible strong-branching
   children too.
2. **Clique table and implied-bound cuts** for the weak big-M bounds (p200x1188c, mc11).
3. **Symmetry** (fhnw-binpack4-4, graph20-20-1rand): orbital fixing on a simple detector.
   Only if time allows.
4. **Structure cuts for the MRPL models**, e.g. (l,S) inequalities. The cases already solve in
   seconds, so this is low priority.
5. **Final results:**
   - Integrate the laptop runs (`results/compare-2`: netlib, cases, miplib60, miplib600,
     miplib600-t6) into `docs/results/comparison.md`.
   - Count SCIP's "gaplimit" status as optimal, including in older tables.
   - Keep numbers from different machines apart.

## How to measure

- **Build and test:** see `CLAUDE.md`.
- **Screening (MIPLIB, 62 instances at 60 s):** `bench/run_miplib.sh 4 60 bench/miplib_small.test
  BIN OUT`. Its CSV has no bound, so it measures only incumbents, which vary by several points
  run to run.
- **Judging cuts, conflicts or search changes:** `bench/ab_bound.sh`, both builds at the same
  time. Report better/worse counts and the mean and median bound distance, and name the
  largest losses too.
- **Comparing with the other solvers:** `bench/compare.sh SET` (header lists the sets); HiGHS
  gets the same thread count.
