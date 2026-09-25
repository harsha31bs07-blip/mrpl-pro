# samaya against the open-source solvers

samaya compared with HiGHS 1.15 (highspy), SCIP 10.0 (PySCIPOpt), CBC 2.10.11 and GLPK 5.0.
Everything was run with `bench/compare.sh` on the same machine, one instance at a time per shard,
every solver single-threaded with the same relative MIP gap (1e-4).
- Unsolved instances count at the time limit in the shifted geometric mean (shift 10 s).
- The test sets were fixed before any solver was run and are reported in full.

Machine: the 4-core cloud container (Intel Xeon, 2.8 GHz), 4 instances at a time. The laptop
runs (Dell G15) are added as they arrive; numbers from the two machines are never mixed.

## Netlib LP (93 instances, 300 s)

| Solver | Solved | Shifted geomean |
|---|---|---|
| HiGHS | 93/93 | 0.17 s |
| CBC (Clp) | 93/93 | 0.19 s |
| **samaya** | **93/93** | **0.32 s** |
| SCIP (SoPlex) | 93/93 | 0.41 s |
| GLPK | 92/93 | 0.76 s |

GLPK's one miss (e226) was a harness bug: an objective constant was lost when converting the
model to LP format for GLPK. It is fixed; after the fix all five agree on e226.

## MRPL case studies and generated MILPs (16 instances, 300 s)

These are `cases/` (planning LP, crude scheduling MILP, utility MILP; three sizes each) and the
seven generated refinery-scheduling, knapsack and facility models.

| Solver | Solved | Shifted geomean |
|---|---|---|
| SCIP | 16/16 | 3.18 s |
| CBC | 16/16 | 3.68 s |
| **samaya** | **16/16** | **3.71 s** |
| HiGHS | 16/16 | 4.63 s |
| GLPK | 12/16 | 19.6 s |

No solver contradicts another on any instance.

## MIPLIB 2017, 62 small "easy" instances, 60 s

| Solver | Solved | Shifted geomean |
|---|---|---|
| SCIP | 18/62 | 39.95 s |
| HiGHS | 16/62 | 42.49 s |
| CBC | 7/62 | 55.21 s |
| **samaya** (build c9a2a31) | **5/62** | **54.57 s** |

**Solved by samaya within 60 s but not by HiGHS or SCIP:**
- markshare_4_0: 12.1 s; no other solver.
- mas76: 53.6 s; CBC 34.7 s.
- pk1: 38.5 s; CBC 45.8 s.

**Solved by HiGHS or SCIP but not by samaya,** grouped by cause (found with solver logs and
profiles):

| Cause | Instances | What the logs show |
|---|---|---|
| Weak cuts on fixed-charge/flow models | p200x1188c, sp150x300d, exp-1-500-5-5, n5-3, mc11, beasleyC3 | Root bound far below the optimum (p200x1188c 7,323 vs 15,078; mc11 1,200 vs 11,689) |
| Strong branching dominates | beasleyC3, mc11, graph20-20-1rand | 77–88% of LP iterations in strong branching |
| Symmetry / feasibility | fhnw-binpack4-4, graph20-20-1rand, enlight_hard, neos-3381206-awhea | No solution found, or infeasibility not proven |
| Slow bound progress | binkar10_1, mik-250-20-75-4, nu25-pr12, pg, neos17, neos-911970 | Good solution, bound closing slowly |

The first two causes are being worked on (see docs/MILP_PLAN.md). Symmetry handling (orbital
fixing) is not planned before the deadline.

**Since this run (d560a26):**
- Aggregated c-MIR with variable bounds and a strong-branching cap solve exp-1-500-5-5 in 2.9 s.
- samaya-only 60 s screening on the same machine: 6/62. This is not yet a same-run comparison;
  the full comparison will be rerun on the final build.
