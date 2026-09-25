# Generated MILP instances

`bench/generate_lps.py --set mip` (refinery scheduling with cargo lots and unit on/off
decisions, multi-dimensional knapsack, capacitated facility location), solved single-threaded
with a relative gap of 1e-4 by samaya (with root cuts) and HiGHS 1.15.1 (highspy, comparison
only). Both solvers agree on every instance within the gap; every samaya solution passed the
verifier.

## Refinery scheduling, knapsack, facility location

* samaya: solved 7/7, shifted geometric mean 5.46 s (unsolved instances count at the 300 s limit)
* highspy: solved 7/7, shifted geometric mean 11.71 s (unsolved instances count at the 300 s limit)

| instance | samaya status | samaya objective | samaya s | highspy status | highspy objective | highspy s |
|---|---|---:|---:|---|---:|---:|
| facility_10x40 | optimal | 2909.771 | 0.01 | optimal | 2909.771 | 0.02 |
| facility_20x80 | optimal | 5709.306443 | 0.06 | optimal | 5709.306443 | 0.13 |
| knapsack_40x3 | optimal | 1426.75 | 0.10 | optimal | 1426.75 | 0.77 |
| knapsack_80x5 | optimal | 2394.59 | 25.94 | optimal | 2394.59 | 234.12 |
| refsched_6c_4p_2u_6t | optimal | 52496.87898 | 0.72 | optimal | 52496.87898 | 2.48 |
| refsched_6c_6p_3u_12t | optimal | 350726.8317 | 2.65 | optimal | 350725.3986 | 5.16 |
| refsched_8c_6p_3u_26t | optimal | 537893.4279 | 32.54 | optimal | 537909.2207 | 35.10 |
