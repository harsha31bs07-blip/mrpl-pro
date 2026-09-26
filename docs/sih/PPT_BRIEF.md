# Brief for the SIH 2026 presentation: samaya (PS 26119, MRPL)

This brief is for whoever builds the slides. It has the content for each slide, the numbers to
use, what to show visually and what to say. Keep the numbers exactly as written; the source
of each is in `docs/results/`. Items marked **[FILL]** are waiting on the final laptop runs or
the GPU teammate.

**Style:**
- 12–14 slides, about 7 minutes.
- One idea per slide, big numbers, little text.
- The name is **samaya** (Sovereign Accelerated Mathematical Algorithms for Yield &
  Allocation).

---

## Slide 1: Title
- **samaya:** an indigenous, GPU-accelerated optimization solver.
- SIH 2026, problem statement 26119 (Mangalore Refinery and Petrochemicals Ltd).
- Team name, members, institute. **[FILL]**
- Visual: a clean logo or wordmark, and a refinery silhouette in the background.

## Slide 2: The problem
- Refineries decide every day which crude to buy, how to run the units, how to blend
  products and how to schedule ships and tanks. These are **LP and MILP optimization
  problems**.
- Today they run on foreign commercial solvers (CPLEX, Gurobi) or open-source ones built
  abroad.
- **The need:** our own solver core, fast, trustworthy, with no external solver code inside.
- Visual: three icons (planning, scheduling, utilities) feeding a "solver" box.

## Slide 3: What we built
- A complete LP/MILP solver **written from scratch in C++20**.
- No external solver or linear-algebra library: our own sparse LU factorization, simplex,
  presolve and branch-and-cut.
- Interfaces: a command-line tool, a C and C++ API, and the standard MPS file format, so it
  works with existing models.
- **Every answer is independently verified** before it is reported.
- Visual: a layer diagram (the architecture on slide 4, simplified).

## Slide 4: Architecture
Draw this as a flow diagram:

```
 Model (MPS) -> Presolve -> LP: dual simplex (own sparse LU) --------------+
                          \-> MILP: branch-and-cut ------------------------+-> Verifier -> Result
                                 (cuts, heuristics, parallel tree search)   |
                 GPU: PDLP first-order LP method (large LPs) -------------+
```
- Say: every box is our code, and the verifier is independent of the solvers.

## Slide 5: How the LP engine works (one slide, not too deep)
- **Dual simplex** with dual steepest-edge pricing and a bound-flipping ratio test.
- **Our own sparse LU**: Markowitz ordering with Forrest–Tomlin updates.
- **Presolve**, which shrinks the model before solving, and **scaling**, which gives numerical
  stability.
- Visual: "model -> smaller model -> solve -> map back -> verify".

## Slide 6: How the MILP engine works
Branch-and-cut, the method the commercial solvers use, and every part of it is ours:
- **Cutting planes:** Gomory, mixed-integer rounding, knapsack covers and flow covers, now also
  inside the search tree.
- **Primal heuristics:** Feasibility Jump (2023 research), feasibility pump, diving, RENS/RINS.
- **Search:** reliability branching, reduced-cost fixing, restarts and a **multi-threaded tree
  search**; conflict analysis is being finished **[FILL: drop this if it isn't merged]**.
- Visual: a search tree with pruned branches in grey.

## Slide 7: GPU acceleration (teammate's part) **[FILL]**
- **PDLP**, a first-order LP method built from matrix-vector products that suits GPUs. It's
  the method Google (PDLP) and NVIDIA (cuOpt) use for huge LPs.
- Show: CPU vs GPU time on large LPs, as a bar chart **[FILL with measured numbers]**.
- Only numbers that were measured; say which GPU (RTX 3060).

## Slide 8: Trust: verified results
- Every result is checked by an **independent verifier** that shares no code with the solver
  and recomputes everything in extended precision:
  - LP optimum: primal and dual feasibility are checked.
  - "Infeasible": a mathematical certificate (Farkas proof) is checked.
  - "Unbounded": a ray is checked.
  - MILP solutions: every constraint and integrality is checked.
- Netlib: **28 of 29 infeasible models proven infeasible with a certificate.**
- Honest line for the speaker: "HiGHS checks its own answers; SCIP has an optional exact mode
  that is 7–10× slower. Ours is on by default and independent of the solver." The MILP bound
  itself isn't certified.
- Visual: a "verified" checkmark stamp on a result.

## Slide 9: Results on standard benchmarks (LP)
Same machine, every solver single-threaded, the same settings.

| Netlib LP (93 models) | Solved | Time (shifted geomean) |
|---|---|---|
| HiGHS | 93/93 | 0.17 s |
| CBC (Clp) | 93/93 | 0.19 s |
| **samaya** | **93/93** | **0.32 s** |
| SCIP (SoPlex) | 93/93 | 0.41 s |
| GLPK | 92/93 | 0.76 s |

- Say: samaya solves every Netlib model correctly, and is faster than SCIP and GLPK.
- Visual: a horizontal bar chart of the time column, with samaya highlighted.

## Slide 10: Results on refinery models (the MRPL-type cases)
Our case studies: refinery planning (LP), crude receipt scheduling (MILP) and captive
power/steam unit commitment (MILP), in three sizes each, plus 7 generated refinery and
facility models. 16 in total.

| Solver | Solved | Time (shifted geomean) |
|---|---|---|
| SCIP | 16/16 | 3.18 s |
| CBC | 16/16 | 3.68 s |
| **samaya** | **16/16** | **3.71 s** |
| HiGHS | 16/16 | 4.63 s |
| GLPK | 12/16 | 19.6 s |

- **[FILL]** Replace the table with the final laptop run (`results/compare-2`, "cases" set).
  The current build solves the largest crude-scheduling case in 2.9 s, down from 20 s when
  this table was measured.
- Say: faster than HiGHS on our refinery cases, and no solver contradicts another on any
  instance.
- Visual: a bar chart, plus one sample plan (a crude schedule Gantt chart from
  `cases/mrpl.py report`).

## Slide 11: Results on MIPLIB (hard, general MILP)
MIPLIB 2017 is the international MILP benchmark: 62 "easy" instances, 60 s each, same
machine.

| Solver | Solved in 60 s |
|---|---|
| SCIP | 18 |
| HiGHS | 16 |
| **samaya** | 7 |
| CBC | 7 |

- samaya solves **three instances that neither HiGHS nor SCIP solves in 60 s**
  (markshare_4_0, mas76, pk1), and proves neos859080 infeasible faster than all of them.
- It has given no wrong answer on any instance.
- It's improving fast. Since this run it gained Feasibility Jump (+5 instances with a
  solution), cuts in the tree, and conflict analysis in development (neos17 solved in 20 s; drop this if it isn't merged).
- **[FILL]** Replace with the final laptop run (60 s and 600 s, and 6 threads).
- Say it honestly: mature solvers took 10–15 years to build; we match CBC after weeks and beat
  HiGHS and SCIP on specific instances.

## Slide 12: Re-planning in a live refinery
- Refineries re-plan every day with revised data. samaya takes yesterday's plan as a start
  (`--mip-start`, matched by name) and keeps it if it's still feasible, so the plan doesn't
  change without reason. Otherwise it keeps the decisions and repairs the rest.
- Demo: `cases/mrpl.py generate --update 1`, then solve day 2 from day 1.
- Visual: "Day 1 plan -> new data -> Day 2 plan", with the unchanged decisions highlighted.

## Slide 13: Impact and what's next
- **Impact:** a sovereign solver core for the refining and petrochemical sector, with no
  licence dependency, verified answers and GPU acceleration for very large models.
- **Next:** symmetry handling, more cutting planes, GPU PDLP for large planning LPs, a convex
  QP (for blending and pooling), and Python bindings.

## Slide 14: Thank you / Q&A
- Repository link, team contacts. **[FILL]**

---

## Questions the jury may ask, and the answers
- **Is it really from scratch?** Yes. There's no external solver or factorization library;
  HiGHS, SCIP, CBC and GLPK are only used in the benchmark scripts as comparisons, never linked.
- **Why slower than HiGHS on MIPLIB?** MILP speed comes from decades of techniques: symmetry
  handling, many cut families, conflict analysis. We implemented the core ones in weeks, and
  each one we add closes part of the gap, which we measure on fixed test sets.
- **How do you know the answers are right?** The independent verifier, plus cross-checking
  every answer against the other solvers and against published optimal values. There's been
  no disagreement.
- **Were the tests picked to look good?** No. The test sets were fixed before any solver ran,
  and every result is reported, including the ones where we lose.

## Sources for every number
- `docs/results/comparison.md`: Netlib, cases, MIPLIB 60 s tables, and the verification
  section.
- `docs/results/netlib.md`: per-model Netlib results and the infeasibility certificates.
- `docs/MILP_PLAN.md`: progress of each MILP technique, with the measurements.
