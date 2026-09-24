# Implementation Plan: Indigenous GPU-Accelerated Optimization Solver

**SIH Problem Statement 26119, MRPL. Working name: `SAMAYA` (Sovereign Accelerated Mathematical Algorithms for Yield & Allocation)**

The rename is optional. The rest of the document uses `samaya` as the namespace.

---

## 1. Goals and non-goals

### What we will build
- A from-scratch solver core for **LP, MILP and convex QP**. It must not depend on HiGHS, CBC, GLPK, SCIP, CPLEX, Gurobi, OSQP or similar code.
- **Numerical robustness first.** That means scaling, presolve, stable LU and Cholesky factorizations, degeneracy handling, and verifying every reported solution.
- **Measurable GPU acceleration** where it pays off: a first-order LP method (PDLP) on very large LPs, plus GPU sparse kernels.
- **Multi-core parallelism:** a concurrent LP race (dual simplex, IPM and PDLP run together) and parallel branch-and-bound.
- **Interfaces:** a C API, a C++ API, Python bindings and a CLI. The solver reads MPS, free-MPS, LP and QPS files.
- **Reproducible benchmarks** on Netlib, MIPLIB 2017, Maros–Mészáros and QPLIB, plus refinery and supply-chain case studies. Results are compared against HiGHS and CBC, and against SCIP or GLPK if time allows.
- A **modular design** that can later take MIQP, NLP and MINLP without rewriting the core.

### What we will not build
- A GUI or a full modeling language. At most we add a thin Python modeling helper for the demo.
- A solver that beats CPLEX or Gurobi overall. The realistic target:
  - Solve 100% of Netlib correctly.
  - Solve the large majority of Maros–Mészáros.
  - Solve a clear, honestly reported subset of MIPLIB 2017 "benchmark/easy".
  - Show real GPU speedups on large LPs.

### What "from scratch" allows
- **Allowed:** C++ standard library, OpenMP, CUDA runtime, cuSPARSE and cuBLAS (as GPU kernel primitives only), pybind11 for bindings, GoogleTest, and CMake.
- **Not allowed:** any external LP, MIP or QP solver, and any external sparse direct factorization package (SuiteSparse, MUMPS, PARDISO, cuDSS). We write our own LU, Cholesky and AMD ordering. This is the part the problem statement cares about.
- **External solvers only for comparison:** HiGHS, CBC and SCIP. They are used only as external baselines in the benchmark harness and are never linked into the solver.

---

## 2. Technology choices

| Concern | Choice | Reason |
|---|---|---|
| Core language | **C++20** | Performance, control over memory layout, templates for index and value types |
| GPU | **CUDA 12** (optional build flag `SAMAYA_CUDA=ON`) | PDLP, SpMV, and batched vector operations |
| Parallelism | OpenMP plus a work-stealing task pool (`std::jthread`) | Concurrent LP, parallel tree search, parallel presolve and cuts |
| Build | CMake with presets; CI builds CPU-only and CUDA | |
| Bindings | C API (stable ABI), pybind11 Python module, CLI `samaya` | |
| Tests | GoogleTest, Python `pytest` for end-to-end tests | |
| Benchmarks | Python harness, results in CSV, performance profiles in matplotlib | |
| Precision | `double` by default, `long double` / double-double for refinement and verification, optional rational checker for small models | Proof of robustness |

---

## 3. Architecture

```
                    +------------------------------------------------+
  CLI / Python /    |                 samaya::Solver                  |
  C API  ---------> |  params, logging, time/limit control, callbacks |
                    +--------------------+---------------------------+
                                         |
          +------------------------------+-------------------------------+
          v                              v                               v
   +-------------+              +-----------------+              +---------------+
   |  io/        |  Model  -->  |  presolve/      |  reduced --> |  dispatcher   |
   |  MPS/LP/QPS |              |  + postsolve    |   model      |  (problem     |
   +-------------+              +-----------------+              |   class)      |
                                                                 +-------+-------+
                          +--------------------+---------------+---------+
                          v                    v               v
                  +---------------+    +---------------+ +-----------------+
                  | lp/           |    | qp/           | | mip/            |
                  |  dual simplex |    |  IPM (convex) | |  branch & cut   |
                  |  primal simplx|    |  active-set   | |  cuts, heur,    |
                  |  IPM + xover  |    |  (warm start) | |  propagation,   |
                  |  PDLP (GPU)   |    +---------------+ |  node queue     |
                  +-------+-------+                      +--------+--------+
                          |                                       |
                          +-------------------+-------------------+
                                              v
                         +-------------------------------------------+
                         | linalg/  sparse CSC/CSR, LU (Markowitz +  |
                         | Forrest–Tomlin), supernodal Cholesky, AMD,|
                         | hyper-sparse solves, GPU SpMV (cuSPARSE)  |
                         +-------------------------------------------+
```

### Repository layout

```
samaya/
  CMakeLists.txt  CMakePresets.json
  include/samaya/            public C++ API (model.hpp, solver.hpp, params.hpp, status.hpp)
  include/samaya_c.h         stable C API
  src/
    core/        Model, bounds, tolerances, logging, timer, params, thread pool
    io/          mps_reader, lp_reader, qps_reader, solution_writer
    linalg/      sparse_matrix, lu_factor, lu_update (FT), cholesky, amd, scaling,
                 dense_kernels, iterative_refinement
    presolve/    reductions/*, postsolve_stack
    lp/          dual_simplex, primal_simplex, pricing, ratio_test, crash,
                 ipm (hsd), crossover, pdlp
    qp/          qp_ipm, qp_active_set
    mip/         bnb_tree, node_queue, branching, propagation, cuts/*, heuristics/*,
                 cut_pool, conflict (later), symmetry (later)
    gpu/         cuda kernels, device vectors, pdlp_kernels.cu
    verify/      solution checker (double, extended, rational)
  apps/cli/      samaya CLI
  python/        pybind11 module + tiny modeling helper
  tests/         unit/, integration/, instances/ (small), regression golden files
  bench/         harness.py, instance lists, baselines, reports/
  cases/         refinery / blending / planning / supply-chain / unit-commitment models
  docs/          algorithm notes, API docs, benchmark report
```

---

## 4. Component design

### 4.1 Core model and I/O
- **Model:** `min c'x + ½x'Qx` subject to `L ≤ Ax ≤ U` and `l ≤ x ≤ u`, with `x_j ∈ ℤ` for integer columns.
  - A is stored column-wise (CSC) plus a lazily built row-wise copy (CSR).
  - Q is stored as its lower triangle in CSC.
- **Parsers:**
  - Fixed and free MPS, including RANGES, BOUNDS (MI/PL/BV/LI/UI/SC) and MARKER INTORG.
  - QPS (QUADOBJ/QMATRIX) and CPLEX-LP format.
  - All parsers stream input and stay robust on 10M+ nonzeros.
- **Outputs:** solution file (primal, dual, reduced costs, basis status) and a JSON run summary.

### 4.2 Numerical linear algebra (`linalg/`), the heart of robustness
1. **Scaling:** geometric-mean scaling (iterated), then equilibration. Keep the scale factors in powers of two so scaling adds no rounding error.
2. **Basis LU for simplex:**
   - Handle singletons first (row and column triangularization).
   - Factor the "nucleus" with Markowitz pivoting and threshold partial pivoting (τ = 0.1, adaptive up to 0.9 on numerical trouble).
   - **Forrest–Tomlin update**, with refactorization triggered by fill-in, an update count limit, or growth in the residual check.
   - **Hyper-sparse FTRAN/BTRAN** (Gilbert–Peierls DFS) for very sparse right-hand sides, which is critical for speed on large models.
   - **Rank-deficiency repair:** replace dependent columns with slacks and report them to the simplex layer.
3. **Cholesky for IPM:**
   - AMD ordering implemented ourselves.
   - Symbolic analysis with an elimination tree, then a **supernodal left-looking Cholesky** that uses dense BLAS-like micro-kernels we write.
   - Handle dense columns separately (Sherman–Morrison or product-form).
   - **Static regularization** (primal-dual regularized IPM) plus dynamic pivot boosting, followed by iterative refinement.
4. **Iterative refinement** in extended precision for final basic solutions and IPM steps.
5. **Condition estimation:** a Hager/Higham 1-norm estimator on the basis, logged and used to trigger stricter pivoting.

### 4.3 LP solvers (`lp/`)

**Bounded dual simplex.** This is the workhorse, and it is required for MIP warm starts.
- Crash basis: CPLEX-style triangular crash, or slack basis as a fallback.
- **Dual steepest-edge pricing** (exact initial weights, Forrest–Goldfarb update). Devex is the fallback.
- **Bound-flipping ratio test** (long-step) combined with a **Harris two-pass** tolerance ratio test.
- **Degeneracy handling:**
  - Cost perturbation, removed at the end, followed by a primal cleanup.
  - Bound shifting.
  - Randomized tie-breaking.
  - Anti-cycling detection with perturbation escalation.
- Phase 1 by the artificial-bounding ("big box") approach, or by composite pricing.
- Dual infeasibility cleanup with primal simplex. Status detection for infeasible and unbounded models, including Farkas certificates.

**Primal simplex.** Used for cleanup after perturbation, for crossover, and as the QP active-set foundation.

**Interior point method.**
- Homogeneous self-dual embedding with Mehrotra predictor-corrector and Gondzio multiple centrality correctors.
- Normal equations with our own Cholesky, or an augmented system with regularized LDLᵀ for models with dense columns.
- Robust detection of infeasibility and unboundedness through the HSD τ/κ variables.

**Crossover.** Takes the IPM solution to an optimal basis, as needed by MIP and for sensitivity analysis. Steps: primal push, dual push, then a cleanup simplex.

**PDLP (GPU and CPU).** A restarted primal-dual hybrid gradient method:
- Adaptive step size, primal weight updates and adaptive restarts.
- Diagonal preconditioning (Ruiz + Pock–Chambolle).
- Feasibility polishing.

On the GPU, the method is almost entirely SpMV (cuSPARSE) plus fused vector kernels. This is where the GPU shows **measurable** gains on LPs with millions of nonzeros. The PDLP result can seed crossover to reach a high-accuracy basic solution.

**Concurrent LP.** Dual simplex, IPM and PDLP (GPU) race on separate threads or devices. The first certified optimum wins and the others are cancelled.

### 4.4 Presolve and postsolve (`presolve/`)

Reductions are applied in rounds until they reach a fixed point:
- Empty and singleton rows and columns.
- Fixed columns.
- Forcing and redundant rows (activity bounds).
- Doubleton equations (substitution).
- Dominated and weakly dominated columns.
- Implied free column substitution.
- Duplicate rows and columns (hashing).
- Dual fixing.

MIP-specific reductions:
- Coefficient tightening.
- Bound tightening by propagation.
- Clique detection.
- Probing on binaries, which gives implications and fixings.
- Integer GCD tightening.

Every reduction pushes an entry onto a **postsolve stack**, which restores the primal, dual and basis (where valid). Presolve is followed by a numerical safety check: reductions that would create coefficients outside [1e-9, 1e9] are refused.

### 4.5 QP (`qp/`)
- **Convex QP IPM:** reuses the LP IPM with the regularized augmented system `[-(Q+X⁻¹Z+ρI) Aᵀ; A δI]`, factored by our own sparse LDLᵀ with quasidefinite ordering.
- **Primal active-set / QP simplex:** used for warm starts, which later serve MIQP branch-and-bound.
- Convexity check: attempt a Cholesky of Q with a small shift and report nonconvexity cleanly.

### 4.6 MILP branch-and-cut (`mip/`)

**Root node**
1. Presolve.
2. LP relaxation (concurrent LP), then cut rounds, with a stall criterion based on the relative gap improvement per round.
3. Root heuristics: simple and ZI rounding, shifting, **feasibility pump**, and fix-and-propagate with several variable orderings.
4. Reduced-cost fixing.

**Cutting planes**, all generated from scratch:
- Gomory mixed-integer cuts from the optimal tableau, with numerical safeguards: limits on dynamism, support size and violation, and rejection of near-parallel cuts.
- MIR / c-MIR with aggregation.
- Knapsack cover and lifted cover cuts.
- Flow cover cuts, which matter for fixed-charge supply-chain and refinery models.
- Clique cuts (from a conflict graph), implied bound cuts and zero-half cuts.
- A **cut pool** with efficacy, orthogonality, parallelism filtering and aging.

**Tree search**
- **Branching:** reliability pseudocost branching (strong branching initializes pseudocosts), with hybrid inference and pseudocost scoring.
- **Node selection:** best-estimate/best-bound hybrid with periodic plunging (depth-first dives) and a switch to best-bound when the gap is small.
- **Node LP:** warm-started dual simplex from the parent basis. Node bases are stored compactly (a diff from the parent).
- **Domain propagation** at every node: activity-based bound tightening and implications.

**Primal heuristics in the tree**
- Diving heuristics: fractional, coefficient, pseudocost, guided and vector-length.
- RINS, RENS and local branching, which run as sub-MIPs through a recursive call to the solver.

**Parallel tree search.** Deterministic, synchronized racing and work sharing: a shared node pool, the incumbent, the pseudocost table and the cut pool. A non-deterministic "opportunistic" mode is also available for maximum speed.

**Later extensions:** conflict analysis, orbital fixing and symmetry handling, and restarts after a root fixing threshold.

**Reported results:** primal bound, dual bound, gap, node count and a certificate check of the incumbent.

### 4.7 Solution verification (`verify/`)

Every returned solution is checked independently of the algorithm:
- Primal feasibility (absolute and relative).
- Integrality.
- Dual feasibility and complementarity for LP and QP.
- The objective, recomputed in extended precision.

A `--verify-exact` mode uses rational arithmetic to confirm optimal bases for small and medium Netlib instances. This gives a strong answer to "prove numerical robustness".

### 4.8 Extensibility hooks
- A `ProblemClass` enum and dispatcher. MIQP will plug in as B&B over `qp_active_set`.
- NLP/MINLP will later plug in as a new relaxation provider: interior-point NLP with our LDLᵀ, plus outer approximation.
- **Callbacks:**
  - Incumbent callback.
  - Lazy constraints.
  - User cuts.
  - Branching.
  - Progress logging.
- These let MRPL domain experts add refinery-specific logic without touching the core.

---

## 5. GPU strategy: measurable benefit only

| Component | GPU? | Justification |
|---|---|---|
| PDLP for large LPs | **Yes, primary GPU story** | Only SpMV and vector operations; published results show 10–100× speedups on huge LPs |
| Batched vector ops in IPM (residuals, step length) | Yes | Cheap to port |
| IPM Cholesky | Stretch goal | Supernodal dense blocks on cuBLAS; only on models with large supernodes |
| Simplex | No | Too sequential and too sparse; the GPU loses |
| MIP: parallel diving / heuristic LP batches with PDLP | Stretch | Many small LPs in parallel on one GPU |

Every GPU path has a CPU fallback. Benchmarks report CPU and GPU times side by side on the same instances.

---

## 6. Benchmarking and validation plan

| Set | Purpose | Target |
|---|---|---|
| **Netlib LP** (≈94 + Kennington) | Correctness, degeneracy (e.g. `pilot*`, `greenbea`, `degen3`, `perold`) | 100% solved, objective within 1e-8 relative to published optima |
| **Netlib infeasible set** | Infeasibility detection | 100% correctly detected |
| **Mittelmann LP feasibility benchmark** (large LPs) | Scalability, GPU PDLP | Solve most; show GPU speedup curve against nonzeros |
| **Maros–Mészáros QP** (138) | QP robustness, ill-conditioning | ≥ 90% solved to 1e-6 |
| **QPLIB** (convex subset) | QP | Report the solved fraction |
| **MIPLIB 2017 "easy" + benchmark subset** | MILP | Staged targets: 30 → 60 → 100+ instances optimal within 1 hr; report gap on the rest |
| **Industrial case studies** (`cases/`) | Relevance to MRPL | All solved; compare with HiGHS/CBC |

### Industrial case studies to model (from open literature)
1. **Crude oil blending and scheduling:** multi-period tanks, CDU feed quality specs. This is a linearized pooling (MILP) formulation.
2. **Refinery production planning:** multi-unit LP (CDU, VDU, FCC, reformer, hydrotreater) with product specs and demand.
3. **Gasoline blending:** LP/MILP with octane, RVP and sulfur specs, and a QP variant for deviation minimization.
4. **Supply chain / depot distribution:** multi-echelon fixed-charge network design, which exercises flow cover cuts.
5. **Power dispatch / unit commitment:** a MILP with min up/down times and a quadratic cost variant (MIQP-ready).
6. **Scaled generators:** instance generators for #1–#5 that produce 1k to 1M variables, to show scaling behavior.

### Metrics and reporting
- Shifted geometric mean of time (shift 10 s) and number solved.
- **Performance profiles** (Dolan–Moré).
- Final MIP gap and primal integral.
- Max primal and dual violation after verification.
- Identical hardware and thread counts across solvers.
- Baselines: **HiGHS and CBC** (required), SCIP (if time allows). Each is run through the same harness with the same time limits.

---

## 7. Robustness demonstrations

These are the evidence we need to show the jury:
1. **Degeneracy:** Netlib `degen2/3`, `pilot87` and highly degenerate refinery models. Show the iteration counts with perturbation and bound-flipping on and off.
2. **Ill-conditioning:** instances with coefficient ranges above 1e8. Show the scaling impact, residuals before and after iterative refinement, and condition estimates.
3. **Weak LP relaxations:** big-M refinery scheduling. Show the root gap closed by each cut family, and the result with and without presolve and propagation.
4. **Verified correctness:** exact rational verification of optimal bases on Netlib.
5. **Stress and fuzz testing:** random perturbations of instances, plus a differential comparison of objective values against HiGHS in CI.

---

## 8. Team structure (6 members, SIH format)

| Member | Ownership |
|---|---|
| **A: Linear algebra lead** | Sparse structures, LU + FT update, hyper-sparse solves, AMD, Cholesky, LDLᵀ |
| **B: Simplex lead** | Dual/primal simplex, pricing, ratio tests, degeneracy, crossover |
| **C: IPM / QP / GPU lead** | HSD IPM, QP IPM, PDLP (CUDA), concurrent LP |
| **D: Presolve + I/O** | Parsers, presolve/postsolve, scaling, verification module |
| **E: MIP lead** | B&B tree, branching, node selection, cuts, heuristics, parallel tree |
| **F: Benchmarks + cases + integration** | Harness, CI, industrial case models, Python API, CLI, report and demo |

---

## 9. Timeline (16 weeks to the grand finale, then a roadmap)

### Phase 0: Foundations (weeks 1–2)
- Repo, CMake, CI (GCC/Clang, sanitizers), coding standards, logging and parameters.
- Model data structure, and MPS/LP readers with tests on all Netlib files.
- Sparse matrix library and scaling.
- Benchmark harness skeleton that runs HiGHS and CBC baselines and stores CSV.

**Exit criterion:** all Netlib, MIPLIB and Maros–Mészáros files parse and their statistics match the published tables.

### Phase 1: A first correct LP (weeks 3–5)
- Dense-then-sparse LU with Markowitz pivoting; Forrest–Tomlin update.
- Bounded primal and dual simplex with Dantzig/Devex pricing and Harris ratio test.
- Solution verifier.

**Exit criterion:** 80% of Netlib solved correctly.

### Phase 2: A robust, fast LP (weeks 6–8)
- Dual steepest edge, bound-flipping ratio test, perturbation, and hyper-sparse solves.
- Crash basis.
- LP presolve and postsolve.
- HSD IPM with our own AMD + supernodal Cholesky, followed by crossover.
- CPU PDLP.

**Exit criterion:** 100% of Netlib, including the infeasible set; Kennington solved; LP times within a reasonable factor of HiGHS.

### Phase 3: GPU and QP (weeks 8–10, overlapping)
- CUDA PDLP and concurrent LP.
- QP IPM with LDLᵀ.
- Maros–Mészáros runs.

**Exit criteria:**
- GPU speedup demonstrated on at least 5 large Mittelmann LPs.
- At least 80% of Maros–Mészáros solved.

### Phase 4: MILP (weeks 9–13)
- B&B with dual simplex warm start, propagation, and pseudocost/reliability branching.
- Best-estimate plus plunging node selection.
- Gomory, MIR, cover, flow cover and clique cuts, with a cut pool.
- Heuristics: rounding, diving, feasibility pump and RINS.
- MIP presolve (probing, coefficient tightening).
- Parallel tree search.

**Exit criterion:** at least 30 MIPLIB 2017 "easy" instances solved to optimality in one hour or less; all case studies solved.

### Phase 5: Hardening and evidence (weeks 14–16)
- Parameter tuning across the whole test set.
- Numerical fixes for the failing instances.
- Rational verification mode.
- Fuzzing.
- Benchmark report with performance profiles.
- Python API polish.
- Demo script and slides.

**Exit criterion:** a frozen release tag, the report in `docs/`, and a reproducible `bench/run_all.sh`.

### After SIH: roadmap
- MIQP (B&B on the QP active-set method), conflict analysis, symmetry and restarts.
- A GPU Cholesky IPM.
- NLP via an interior-point method with filter line search; MINLP via outer approximation and spatial B&B.
- An exact/rational MIP mode for certification.
- A deployment pilot on MRPL planning models.

---

## 10. Hackathon demo storyline (grand finale)
1. **Transparency:** walk through the architecture and open the code for the dual simplex and the LU update.
2. **Correctness:** a live run of all Netlib instances, verified, with a rational checker spot-check.
3. **Robustness:** a degenerate / ill-conditioned instance, with features toggled on and off (perturbation, scaling, refinement).
4. **GPU:** a large LP with CPU IPM, CPU PDLP and GPU PDLP compared side by side, with a speedup chart.
5. **MILP:** a refinery crude-scheduling model with a live progress log (gap closing, cuts, heuristics). Compared against CBC and HiGHS.
6. **Report:** performance profiles and a solved-instance table against HiGHS and CBC.

---

## 11. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Our LU/Cholesky is slow or unstable | Build it first; test against dense references and random matrices; add threshold escalation and refactorization triggers |
| Simplex cycling or stalling on degenerate models | Perturbation, bound shifting, randomized pricing, IPM + crossover fallback via concurrent LP |
| MIP performance far behind commercial solvers | Focus on the case studies and the MIPLIB "easy" set; report honestly; put effort into presolve, propagation and cuts, which give the biggest gains |
| GPU work consumes schedule | GPU limited to PDLP plus vector kernels; CPU PDLP first, CUDA port second |
| Scope creep into NLP/MINLP | Architecture hooks only; implementation deferred to the roadmap |
| Numerical false claims (wrong "optimal") | A mandatory verifier on every solve; CI differential tests against HiGHS objectives |

---

## 12. Immediate next steps
1. Scaffold the repository layout from §3, with CMake, CI and a `samaya` CLI stub.
2. Implement `Model` + the MPS reader + the benchmark harness. Download the Netlib and MIPLIB subsets with a script, not committed.
3. Start the `linalg/` sparse LU in parallel with a dense-LU reference simplex, so the other members can start work immediately.
