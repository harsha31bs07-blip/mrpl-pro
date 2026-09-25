# MIPLIB 2017 results

The 62 MIPLIB 2017 benchmark instances with status "easy" and at most 20,000 nonzeros
(`bench/miplib_small.test`), each solved single-threaded on a 4-core Intel Xeon (2.8 GHz) cloud
machine, four instances at a time. Reproduce with:

```sh
bench/fetch_instances.sh miplib-list bench/miplib_small.test
bench/run_miplib.sh 4 600        # or 4 60 for the screening run
```

"Solved" means optimal within the default relative gap of 1e-4, or proven infeasible. Every
reported solution passed the verifier on the original model, and every solved instance agrees
with the published optimum (`miplib2017-v31.solu`): **no wrong answers in any run**.

## Summary

| Build | 10 min: solved | 60 s: solved | 60 s: feasible solution found | 60 s: median gap to optimum |
|---|---|---|---|---|
| Branch-and-bound, no cuts | 7/62 (1 crash: out of memory on mik-250-20-75-4) | 3/62 | 36/62 | 5.8% |
| With root cuts (first version) | **8/62** (pg newly solved), no crash | – | – | – |
| With root cuts, tuned | not rerun yet | 3/62 | **37/62** | **2.8%** |

* Root cuts (Gomory mixed-integer, c-MIR, knapsack covers) solve `pg` and give better
  solutions on more instances than they make worse (15 better, 12 worse at 60 s), and they
  halve the median gap between the best solution found and the known optimum.
* The first cut version slowed some solved instances (qap10 111 s → 276 s: cut rounds on a
  14 s root LP had no time budget; markshare_4_0 18 s → 38 s: cuts that did not move the bound
  were kept). The tuned version bounds the cut loop's time and drops cuts that do not help;
  markshare_4_0 is back to 18 s.
* The out-of-memory crash (6.7 GB of open nodes) is fixed: node paths are shared, and beyond
  2.5 million open nodes the search switches to depth first.
* This is plain branch-and-bound with root cuts: no cuts in the tree, no MIP presolve beyond
  the LP reductions, no parallel search. Most instances stop at the time limit. The next steps
  are in docs/ARCHITECTURE.md (flow covers, probing, better primal heuristics, parallel tree
  search).

## 10 minutes per instance

* cuts: solved 8/62, shifted geometric mean 470.22 s (unsolved instances count at the 600 s limit)
* no-cuts: solved 7/62, shifted geometric mean 457.92 s (unsolved instances count at the 600 s limit)

| instance | known optimum | cuts status | cuts objective | cuts s | no-cuts status | no-cuts objective | no-cuts s |
|---|---:|---|---:|---:|---|---:|---:|
| 50v-10 | 3311.179984 | time_limit | 3627.219989 | 600.22 | time_limit | 3794.219996 | 600.83 |
| app1-1 | -3 | optimal | -3 | 233.35 | optimal | -3 | 100.21 |
| assign1-5-8 | 212 | time_limit | 220 | 600.61 | time_limit | 214 | 600.56 |
| b1c1s1 | 24544.25 | time_limit |  | 600.05 | time_limit |  | 600.07 |
| beasleyC3 | 754 | time_limit | 1186 | 600.14 | time_limit | 855 | 600.14 |
| binkar10_1 | 6742.199884 | time_limit | 6751.550014 | 600.18 | time_limit | 6762.820018 | 600.19 |
| cost266-UUE | 25148940.56 | time_limit |  | 600.08 | time_limit | 27022909.7 | 600.11 |
| csched007 | 351 | time_limit | 472.9999999 | 600.28 | time_limit | 514 | 600.35 |
| csched008 | 173 | time_limit |  | 600.21 | time_limit | 186 | 600.29 |
| cvs16r128-89 | -97 | time_limit |  | 600.24 | time_limit | -20 | 600.04 |
| enlight_hard | 37 | time_limit | 37 | 601.02 | time_limit |  | 602.04 |
| exp-1-500-5-5 | 65887 | time_limit | 148832 | 600.28 | time_limit | 185493 | 600.44 |
| fastxgemm-n2r6s0t2 | 230 | time_limit | 236 | 600.05 | time_limit | 236 | 600.05 |
| fhnw-binpack4-4 | inf | time_limit |  | 600.28 | time_limit |  | 600.35 |
| fhnw-binpack4-48 | 0 | time_limit |  | 600.05 | time_limit |  | 600.05 |
| gen-ip002 | -4783.733392 | node_limit | -4776.045637 | 547.83 | time_limit | -4775.400852 | 607.32 |
| gen-ip054 | 6840.965642 | node_limit | 6848.339683 | 404.88 | time_limit | 6855.087531 | 608.03 |
| glass4 | 1200012600 | time_limit |  | 600.92 | time_limit |  | 600.86 |
| gmu-35-40 | -2406733.369 | time_limit | -2403082.723 | 600.78 | time_limit | -2399498.971 | 600.90 |
| gmu-35-50 | -2607958.33 | time_limit | -2601252.67 | 600.69 | time_limit | -2605348.8 | 600.68 |
| graph20-20-1rand | -9 | time_limit |  | 600.03 | time_limit |  | 600.03 |
| graphdraw-domain | 19685.99998 | time_limit | 22953 | 600.49 | time_limit | 37840 | 600.42 |
| ic97_potential | 3941.999931 | time_limit |  | 600.39 | time_limit |  | 600.30 |
| lotsize | 1480195 | time_limit |  | 600.10 | time_limit |  | 600.15 |
| mad | 0.0268 | time_limit | 0.1232 | 601.89 | time_limit | 0.212 | 602.46 |
| markshare_4_0 | 1 | optimal | 1 | 38.12 | optimal | 1 | 18.20 |
| mas74 | 11801.18572 | time_limit | 11801.18573 | 603.70 | time_limit | 11857.37126 | 601.92 |
| mas76 | 40005.05399 | optimal | 40005.05414 | 77.48 | optimal | 40005.05414 | 55.22 |
| mc11 | 11689 | time_limit | 18167 | 600.12 | time_limit | 13423 | 600.11 |
| mcsched | 211913 | time_limit | 215860 | 600.03 | time_limit | 211950 | 600.06 |
| mik-250-20-75-4 | -52301 | time_limit | -52301 | 600.20 | crash |  | 413.36 |
| milo-v12-6-r2-40-1 | 326481.1428 | time_limit |  | 600.08 | time_limit |  | 600.03 |
| n5-3 | 8105 | time_limit | 9998.161796 | 600.11 | time_limit | 9715 | 600.13 |
| neos-2657525-crna | 1.810748 | time_limit | 150.9285312 | 600.45 | time_limit | 46.0625012 | 600.73 |
| neos-3024952-loue | 26756 | time_limit |  | 600.04 | time_limit | 32266 | 600.23 |
| neos-3046615-murg | 1600 | time_limit | 1721 | 601.31 | time_limit |  | 600.85 |
| neos-3381206-awhea | 453 | time_limit |  | 600.18 | time_limit |  | 600.53 |
| neos-3627168-kasai | 988585.62 | time_limit |  | 600.15 | time_limit |  | 600.14 |
| neos-3754480-nidda | 12941.73839 | time_limit | 14457.60947 | 600.80 | time_limit | 14457.60947 | 600.54 |
| neos-4338804-snowy | 1471 | time_limit | 2836 | 600.21 | time_limit | 3994 | 600.29 |
| neos-4738912-atrato | 283627956.6 | time_limit | 320918644 | 600.03 | time_limit | 366313164.9 | 600.18 |
| neos-4954672-berkel | 2612710 | time_limit | 5265102 | 600.29 | time_limit | 5145677 | 600.29 |
| neos-5107597-kakapo | 3645 | time_limit |  | 600.03 | time_limit |  | 600.04 |
| neos-911970 | 54.76 | time_limit | 56.62999985 | 600.50 | time_limit | 60.82 | 601.02 |
| neos17 | 0.1500025774 | time_limit | 0.1516801798 | 600.15 | time_limit | 0.151943037 | 600.07 |
| neos5 | 15 | time_limit | 15 | 600.55 | time_limit | 15 | 600.99 |
| neos859080 | inf | infeasible |  | 0.17 | infeasible |  | 0.18 |
| nu25-pr12 | 53905 | optimal | 53905 | 93.84 | optimal | 53905 | 97.82 |
| p200x1188c | 15078 | time_limit | 15078 | 600.13 | time_limit | 15078 | 600.12 |
| pg | -8674.342607 | optimal | -8674.342708 | 239.17 | time_limit | -7088.799211 | 600.35 |
| pg5_34 | -14339.35345 | time_limit | -14324.80782 | 600.06 | time_limit | -14211.45426 | 600.52 |
| pk1 | 11 | optimal | 11 | 80.07 | optimal | 11 | 74.29 |
| qap10 | 340 | optimal | 340 | 276.06 | optimal | 340 | 111.46 |
| ran14x18-disj-8 | 3712 | time_limit | 3770 | 600.27 | time_limit | 3752.000003 | 600.25 |
| reblock115 | -36800603.23 | time_limit |  | 600.04 | time_limit |  | 600.15 |
| rococoB10-011000 | 19449 | time_limit |  | 600.03 | time_limit |  | 600.07 |
| rococoC10-001000 | 11460 | time_limit | 18130 | 600.03 | time_limit |  | 600.08 |
| sp150x300d | 69 | time_limit | 71 | 600.36 | time_limit | 69 | 600.37 |
| supportcase26 | 1745.123813 | time_limit | 1844.652692 | 600.79 | time_limit | 1925.373902 | 600.66 |
| timtab1 | 764772 | time_limit |  | 601.00 | time_limit |  | 600.90 |
| tr12-30 | 130596 | time_limit | 138460.9999 | 600.20 | time_limit |  | 600.37 |
| uct-subprob | 314 | time_limit | 367.9999999 | 600.07 | time_limit | 355 | 600.07 |

## 60 seconds per instance (screening)

* cuts-tuned: solved 3/62, shifted geometric mean 56.84 s (unsolved instances count at the 60 s limit)
* no-cuts: solved 3/62, shifted geometric mean 56.78 s (unsolved instances count at the 60 s limit)

| instance | known optimum | cuts-tuned status | cuts-tuned objective | cuts-tuned s | no-cuts status | no-cuts objective | no-cuts s |
|---|---:|---|---:|---:|---|---:|---:|
| 50v-10 | 3311.179984 | time_limit | 3925.369993 | 60.03 | time_limit | 3794.219996 | 60.17 |
| app1-1 | -3 | time_limit |  | 60.02 | time_limit |  | 60.03 |
| assign1-5-8 | 212 | time_limit | 220 | 60.03 | time_limit | 214 | 60.04 |
| b1c1s1 | 24544.25 | time_limit |  | 60.02 | time_limit |  | 60.02 |
| beasleyC3 | 754 | time_limit | 1254 | 60.01 | time_limit | 895 | 60.02 |
| binkar10_1 | 6742.199884 | time_limit | 6761.930013 | 60.04 | time_limit | 6762.820018 | 60.06 |
| cost266-UUE | 25148940.56 | time_limit |  | 60.02 | time_limit |  | 60.03 |
| csched007 | 351 | time_limit |  | 60.02 | time_limit | 514 | 60.04 |
| csched008 | 173 | time_limit | 186 | 60.02 | time_limit | 186 | 60.03 |
| cvs16r128-89 | -97 | time_limit |  | 60.31 | time_limit |  | 60.02 |
| enlight_hard | 37 | time_limit | 37 | 60.09 | time_limit |  | 60.29 |
| exp-1-500-5-5 | 65887 | time_limit | 159212 | 60.03 | time_limit | 212086 | 60.13 |
| fastxgemm-n2r6s0t2 | 230 | time_limit | 236 | 60.02 | time_limit | 527 | 60.02 |
| fhnw-binpack4-4 | inf | time_limit |  | 60.05 | time_limit |  | 60.04 |
| fhnw-binpack4-48 | 0 | time_limit |  | 60.01 | time_limit |  | 60.02 |
| gen-ip002 | -4783.733392 | time_limit | -4769.740589 | 60.74 | time_limit | -4772.326056 | 60.51 |
| gen-ip054 | 6840.965642 | time_limit | 6864.088667 | 61.10 | time_limit | 6869.871195 | 60.78 |
| glass4 | 1200012600 | time_limit |  | 60.27 | time_limit |  | 60.12 |
| gmu-35-40 | -2406733.369 | time_limit | -2399116.214 | 60.22 | time_limit | -2399116.214 | 60.14 |
| gmu-35-50 | -2607958.33 | time_limit | -2601559.66 | 60.17 | time_limit |  | 60.14 |
| graph20-20-1rand | -9 | time_limit |  | 60.03 | time_limit |  | 60.02 |
| graphdraw-domain | 19685.99998 | time_limit | 39996 | 60.05 | time_limit | 37840 | 60.05 |
| ic97_potential | 3941.999931 | time_limit |  | 60.05 | time_limit |  | 60.07 |
| lotsize | 1480195 | time_limit |  | 60.01 | time_limit |  | 60.04 |
| mad | 0.0268 | time_limit | 0.2168 | 60.37 | time_limit | 0.294 | 60.34 |
| markshare_4_0 | 1 | optimal | 1 | 18.06 | optimal | 1 | 17.97 |
| mas74 | 11801.18572 | time_limit | 11945.26305 | 60.59 | time_limit | 12015.08994 | 60.55 |
| mas76 | 40005.05399 | time_limit | 40005.05414 | 60.22 | optimal | 40005.05414 | 55.01 |
| mc11 | 11689 | time_limit | 20071 | 60.01 | time_limit | 13423 | 60.02 |
| mcsched | 211913 | time_limit | 215860 | 60.01 | time_limit | 218781 | 60.01 |
| mik-250-20-75-4 | -52301 | time_limit | -52301 | 60.04 | time_limit | -40314 | 60.67 |
| milo-v12-6-r2-40-1 | 326481.1428 | time_limit |  | 60.02 | time_limit |  | 60.03 |
| n5-3 | 8105 | time_limit | 10906.42462 | 60.01 | time_limit | 9725 | 60.02 |
| neos-2657525-crna | 1.810748 | time_limit | 160.4878077 | 60.13 | time_limit | 160.4878077 | 60.11 |
| neos-3024952-loue | 26756 | time_limit |  | 60.02 | time_limit |  | 60.03 |
| neos-3046615-murg | 1600 | time_limit |  | 60.28 | time_limit |  | 60.11 |
| neos-3381206-awhea | 453 | time_limit |  | 60.12 | time_limit |  | 60.21 |
| neos-3627168-kasai | 988585.62 | time_limit |  | 60.02 | time_limit |  | 60.02 |
| neos-3754480-nidda | 12941.73839 | time_limit | 14553.9693 | 60.12 | time_limit | 14553.9693 | 60.08 |
| neos-4338804-snowy | 1471 | time_limit | 3994 | 60.05 | time_limit | 3994 | 60.05 |
| neos-4738912-atrato | 283627956.6 | time_limit |  | 60.02 | time_limit | 366313164.9 | 60.07 |
| neos-4954672-berkel | 2612710 | time_limit | 5265102 | 60.03 | time_limit | 5327763.667 | 60.05 |
| neos-5107597-kakapo | 3645 | time_limit |  | 60.02 | time_limit |  | 60.02 |
| neos-911970 | 54.76 | time_limit | 96.81 | 60.40 | time_limit | 83.75 | 60.21 |
| neos17 | 0.1500025774 | time_limit | 0.1550515385 | 60.02 | time_limit | 0.1550515385 | 60.02 |
| neos5 | 15 | time_limit | 15 | 60.09 | time_limit | 15 | 60.09 |
| neos859080 | inf | infeasible |  | 0.19 | infeasible |  | 0.18 |
| nu25-pr12 | 53905 | time_limit | 53905 | 60.02 | time_limit | 53905 | 60.02 |
| p200x1188c | 15078 | time_limit | 15078 | 60.03 | time_limit | 15078 | 60.03 |
| pg | -8674.342607 | time_limit | -8652.219875 | 60.01 | time_limit | -7088.799211 | 60.15 |
| pg5_34 | -14339.35345 | time_limit | -14313.3255 | 60.01 | time_limit | -14211.45426 | 60.14 |
| pk1 | 11 | optimal | 11 | 58.30 | time_limit | 11 | 60.02 |
| qap10 | 340 | time_limit |  | 60.02 | time_limit |  | 60.02 |
| ran14x18-disj-8 | 3712 | time_limit | 3843.999999 | 60.02 | time_limit | 3867 | 60.02 |
| reblock115 | -36800603.23 | time_limit |  | 60.02 | time_limit |  | 60.02 |
| rococoB10-011000 | 19449 | time_limit |  | 60.02 | time_limit |  | 60.01 |
| rococoC10-001000 | 11460 | time_limit |  | 60.01 | time_limit |  | 60.02 |
| sp150x300d | 69 | time_limit | 71 | 60.04 | time_limit | 69 | 60.04 |
| supportcase26 | 1745.123813 | time_limit |  | 60.14 | time_limit |  | 60.10 |
| timtab1 | 764772 | time_limit |  | 60.07 | time_limit |  | 60.21 |
| tr12-30 | 130596 | time_limit | 138460.9999 | 60.02 | time_limit |  | 60.07 |
| uct-subprob | 314 | time_limit | 368 | 60.02 | time_limit | 355 | 60.02 |
