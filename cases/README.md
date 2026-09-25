# MRPL case studies

Three refinery decision problems of the kind MRPL solves every day, generated from public
order-of-magnitude figures for a 15 MMTPA coastal refinery (about 45 kt of crude a day). They
are **not MRPL data**; the structure is what matters. Quantities are in kilotonnes (kt),
money in lakh rupees.

```sh
cmake --preset release && cmake --build --preset release
cases/demo.sh            # small sizes, under a second
cases/demo.sh large      # the largest sizes, about 25 s on a 4-core cloud machine
```

The demo generates the models (`cases/mrpl.py generate`), solves each with samaya (every
solution is checked by the independent verifier) and prints the plan (`cases/mrpl.py
report MODEL SOLUTION`).

| Model | Type | Decisions | Sizes (rows x columns, integers) |
|---|---|---|---|
| `plan` | LP | weekly crude slate, unit feeds, product routing, sales and exports | 4 weeks: 240 x 316; 13 weeks: 806 x 1,105; 52 weeks: 3,328 x 4,732 |
| `crude` | MILP | which day each crude cargo berths, spot cargoes, daily CDU charge by crude, tank stocks | 14 days (8 cargoes): 164 x 251, 83 integer; up to 30 days (22 cargoes): 354 x 652, 292 integer |
| `utility` | MILP | hourly on/off and load of boilers and a gas turbine, steam turbine, grid import and export | 24 h: 801 x 408, 192 integer; up to 96 h: 3,321 x 1,632, 768 integer |

## What each model captures

**Planning (LP).**
- Four to eight assayed crudes (API and sulfur); the CDU yields follow the API gravity.
- A 2.5% sulfur limit on the CDU blend (metallurgy).
- VDU, FCC, hydrocracker, naphtha reformer and diesel hydrotreater (DHDS), each with capacity
  and yields.
- Gasoline blending to BS-VI (RON ≥ 91, sulfur ≤ 10 ppm).
- All gasoil through the DHDS for 10 ppm diesel, and at most 20% LCO in the DHDS feed.
- Crude and product tanks.
- Domestic demand, plus exports at a discount.

The report includes the **marginal value of each binding unit capacity** (the LP duals). For
example, one more kt a week of reformer capacity is worth about 1,800 lakh rupees in the small
case. This is the kind of number that decides debottlenecking investments.

**Crude receipt scheduling (MILP).**
- Term cargoes (Aframax and Suezmax sizes) with arrival windows.
- One single-point mooring, so at most one cargo a day.
- Demurrage after the second day of a window, and a contract penalty for a cargo that is not
  lifted.
- Optional spot cargoes.
- Tank capacity per crude and for the tank farm.
- The CDU runs between its minimum and maximum throughput with a blend sulfur limit.
- At most one planned stop, and a closing stock.

**Captive power and steam (MILP).**
- Three boilers and a gas turbine with heat recovery, each with on/off decisions, minimum and
  maximum loads, startup costs, a minimum up time and ramp limits.
- A back-pressure steam turbine and an HP-to-MP letdown.
- Grid import with a peak and off-peak tariff, and export.
- Hourly steam and power demand with a day profile.

## Correctness

All nine models agree with HiGHS 1.15, used only as a reference: 9 of 9 status matches, with
objectives within the 1e-4 MIP gap.

```sh
python3 cases/mrpl.py generate --out cases/instances
python3 bench/harness.py cases/instances/*.mps --baseline highspy --time-limit 120
```

The crude scheduling case found a real bug in the MILP solver. A node whose LP solution was
integral only within the tolerance was pruned: a 0.9999995 cargo binary times 130 kt broke a
tank balance once rounded. The solver then reported "infeasible". It is fixed, and the model is
now a regression test (`tests/instances/mrpl_crude_small.mps`).
