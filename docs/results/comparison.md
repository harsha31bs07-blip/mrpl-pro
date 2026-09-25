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

Two same-machine runs: before and after the overnight work of 25–26 September.

| Solver | Solved | Shifted geomean |
|---|---|---|
| SCIP | 18/62 | 39.91 s |
| HiGHS | 16/62 | 42.57 s |
| **samaya**, build cda171e | **7/62** | **50.99 s** |
| CBC | 7/62 | 55.26 s |
| samaya, build c9a2a31 (before the night) | 5/62 | 54.57 s |

The HiGHS, SCIP and CBC numbers above are from the second run; in the first they were 16, 18
and 7 with geomeans within 1%.

Solve times of every instance samaya solves (seconds; "–" = not solved in 60 s):

| Instance | samaya | HiGHS | SCIP | CBC |
|---|---|---|---|---|
| markshare_4_0 | **14.4** | – | – | – |
| mas76 | **25.5** | – | – | 35.0 |
| pk1 | **43.0** | – | – | 45.9 |
| neos859080 (infeasible) | **0.1** | 1.8 | 0.8 | – |
| exp-1-500-5-5 | 3.5 | 4.2 | **3.0** | – |
| app1-1 | 13.5 | 25.5 | 8.5 | **6.3** |
| sp150x300d | 0.8 | **0.1** | 0.4 | – |

- samaya solves three instances that neither HiGHS nor SCIP solves in 60 s (markshare_4_0,
  mas76, pk1), and proves neos859080 infeasible faster than any of them.
- HiGHS or SCIP solve 14 instances that samaya does not:

| Cause | Instances | What the logs show |
|---|---|---|
| Big-M network gaps | p200x1188c, mc11, beasleyC3, n5-3 | Root bound still far below the optimum (p200x1188c 6,020 vs 15,078; mc11 1,156 vs 11,689) |
| No solution for general-integer / symmetric models | neos-3381206-awhea, enlight_hard, fhnw-binpack4-4, graph20-20-1rand | No incumbent: the pump ignores interior general integers, no symmetry handling |
| Slow bound progress | binkar10_1, mik-250-20-75-4, nu25-pr12, pg, neos17, neos-911970 | Good solution, bound closing slowly (no cuts in the tree) |

**Fixed overnight (all in docs/MILP_PLAN.md):**
- Aggregated c-MIR with variable bounds: exp-1-500-5-5 in 3.5 s.
- Flow covers: sp150x300d in 0.8 s.
- Reduced-cost fixing and a root restart: mas76 from 53.6 s to 25.5 s.
- A cap on strong-branching effort.
- Three simplex speed-ups: a warm-started solve factorized three times, now once.
