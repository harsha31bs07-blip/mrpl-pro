#!/usr/bin/env python3
"""Generates structured LP and MILP instances as free-format MPS files.

LP families (feasible and bounded by construction):
  transport   balanced transportation problem, S sources x D sinks
  refinery    multi-period refinery planning: crude purchase, processing, product yields,
              inventories, CDU capacity and a sulfur blending constraint per period
  sparse      random sparse LP with row bounds built around a feasible point

MILP families:
  refsched    refinery scheduling: the refinery model with crude bought in whole cargo lots
              (integer), processing units switched on and off per period (binary) with minimum
              and maximum throughput, fixed running costs and startup costs
  knapsack    multi-dimensional 0-1 knapsack with correlated weights and profits
  facility    capacitated facility location: binary open decisions with fixed costs, continuous
              assignment of customer demand, capacity and linking constraints

Examples:
    bench/generate_lps.py --set lp  --out bench/instances/generated
    bench/generate_lps.py --set mip --out bench/instances/generated-mip
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path


class MpsWriter:
    """Accumulates an LP and writes it in free MPS format."""

    def __init__(self, name: str, maximize: bool = False) -> None:
        self.name = name
        self.maximize = maximize
        self.rows: list[tuple[str, str, float, float]] = []  # name, type, lo, up
        self.cols: dict[str, dict[str, float]] = {}
        self.obj: dict[str, float] = {}
        self.bounds: dict[str, tuple[float | None, float | None]] = {}
        self.integer: set[str] = set()

    def col(self, name: str, cost: float = 0.0, lo: float | None = 0.0, up: float | None = None,
            integer: bool = False):
        self.cols[name] = {}
        self.obj[name] = cost
        self.bounds[name] = (lo, up)
        if integer:
            self.integer.add(name)

    def row(self, name: str, terms: dict[str, float], lo: float | None, up: float | None):
        if lo is not None and up is not None and lo == up:
            kind = "E"
        elif up is not None and lo is None:
            kind = "L"
        elif lo is not None and up is None:
            kind = "G"
        else:
            kind = "R"  # Ranged: G with a range.
        self.rows.append((name, kind, lo, up))
        for c, v in terms.items():
            if v != 0.0:
                self.cols[c][name] = self.cols[c].get(name, 0.0) + v

    def write(self, path: Path) -> None:
        lines = [f"NAME {self.name}"]
        if self.maximize:
            lines += ["OBJSENSE", "    MAX"]
        lines.append("ROWS")
        lines.append(" N obj")
        for name, kind, _, _ in self.rows:
            lines.append(f" {'G' if kind == 'R' else kind} {name}")
        lines.append("COLUMNS")
        in_int = False
        for c, entries in self.cols.items():
            if (c in self.integer) != in_int:
                in_int = not in_int
                lines.append(f" M{len(lines)} 'MARKER' '{'INTORG' if in_int else 'INTEND'}'")
            # An empty column must still be declared, or a later BOUNDS entry names an unknown
            # column; an explicit zero objective entry does that.
            if self.obj[c] != 0.0 or not entries:
                lines.append(f" {c} obj {self.obj[c]:.12g}")
            for r, v in entries.items():
                lines.append(f" {c} {r} {v:.12g}")
        if in_int:
            lines.append(f" M{len(lines)} 'MARKER' 'INTEND'")
        lines.append("RHS")
        for name, kind, lo, up in self.rows:
            rhs = up if kind == "L" else lo
            if rhs:
                lines.append(f" rhs {name} {rhs:.12g}")
        ranged = [(n, lo, up) for n, k, lo, up in self.rows if k == "R"]
        if ranged:
            lines.append("RANGES")
            for name, lo, up in ranged:
                lines.append(f" rng {name} {up - lo:.12g}")
        lines.append("BOUNDS")
        for c, (lo, up) in self.bounds.items():
            if c in self.integer:
                # Readers disagree on default bounds of integer columns; always write both.
                lines.append(f" MI bnd {c}" if lo is None else f" LO bnd {c} {lo:.12g}")
                lines.append(f" PL bnd {c}" if up is None else f" UP bnd {c} {up:.12g}")
                continue
            if lo is None and up is None:
                lines.append(f" FR bnd {c}")
                continue
            if lo is None:
                lines.append(f" MI bnd {c}")
            elif lo != 0.0:
                lines.append(f" LO bnd {c} {lo:.12g}")
            if up is not None:
                lines.append(f" UP bnd {c} {up:.12g}")
        lines.append("ENDATA")
        path.write_text("\n".join(lines) + "\n")


def transport(sources: int, sinks: int, rng: random.Random) -> MpsWriter:
    lp = MpsWriter(f"TRANSPORT_{sources}x{sinks}")
    supply = [rng.randint(50, 150) for _ in range(sources)]
    demand = [0] * sinks
    remaining = sum(supply)
    for j in range(sinks - 1):
        demand[j] = remaining * rng.uniform(0.5, 1.5) / (sinks - j)
        demand[j] = round(min(demand[j], remaining), 3)
        remaining -= demand[j]
    demand[-1] = remaining
    for i in range(sources):
        for j in range(sinks):
            lp.col(f"x_{i}_{j}", cost=rng.randint(1, 100))
    for i in range(sources):
        lp.row(f"s_{i}", {f"x_{i}_{j}": 1.0 for j in range(sinks)}, None, supply[i])
    for j in range(sinks):
        lp.row(f"d_{j}", {f"x_{i}_{j}": 1.0 for i in range(sources)}, demand[j], None)
    return lp


def refinery(crudes: int, products: int, periods: int, rng: random.Random) -> MpsWriter:
    lp = MpsWriter(f"REFINERY_{crudes}c_{products}p_{periods}t", maximize=True)
    price = [rng.uniform(60, 120) for _ in range(products)]
    crude_cost = [rng.uniform(40, 80) for _ in range(crudes)]
    sulfur = [rng.uniform(0.1, 3.5) for _ in range(crudes)]
    max_sulfur = 1.8
    # Yields per crude sum to 0.97 (3% losses).
    yields = []
    for _ in range(crudes):
        w = [rng.random() for _ in range(products)]
        s = sum(w)
        yields.append([0.97 * v / s for v in w])
    cdu_cap = 1000.0
    for t in range(periods):
        for k in range(crudes):
            lp.col(f"buy_{k}_{t}", cost=-crude_cost[k] * rng.uniform(0.9, 1.1), up=rng.uniform(50, 200))
            lp.col(f"proc_{k}_{t}")
            lp.col(f"cinv_{k}_{t}", cost=-0.5, up=300.0)
        for p in range(products):
            lp.col(f"make_{p}_{t}")
            lp.col(f"sell_{p}_{t}", cost=price[p] * rng.uniform(0.9, 1.1), up=rng.uniform(40, 150))
            lp.col(f"pinv_{p}_{t}", cost=-0.8, up=200.0)
    for t in range(periods):
        for k in range(crudes):
            terms = {f"buy_{k}_{t}": 1.0, f"proc_{k}_{t}": -1.0, f"cinv_{k}_{t}": -1.0}
            if t > 0:
                terms[f"cinv_{k}_{t - 1}"] = 1.0
            lp.row(f"cbal_{k}_{t}", terms, 0.0, 0.0)
        for p in range(products):
            terms = {f"make_{p}_{t}": 1.0}
            for k in range(crudes):
                terms[f"proc_{k}_{t}"] = -yields[k][p]
            lp.row(f"yield_{p}_{t}", terms, 0.0, 0.0)
            terms = {f"make_{p}_{t}": 1.0, f"sell_{p}_{t}": -1.0, f"pinv_{p}_{t}": -1.0}
            if t > 0:
                terms[f"pinv_{p}_{t - 1}"] = 1.0
            lp.row(f"pbal_{p}_{t}", terms, 0.0, 0.0)
        lp.row(f"cdu_{t}", {f"proc_{k}_{t}": 1.0 for k in range(crudes)}, None, cdu_cap)
        lp.row(f"sulfur_{t}", {f"proc_{k}_{t}": sulfur[k] - max_sulfur for k in range(crudes)},
               None, 0.0)
    return lp


def sparse(rows: int, cols: int, rng: random.Random) -> MpsWriter:
    lp = MpsWriter(f"SPARSE_{rows}x{cols}")
    point = []
    for j in range(cols):
        lo, up = 0.0, rng.choice([10.0, 50.0, 100.0])
        lp.col(f"x{j}", cost=rng.randint(-10, 10), lo=lo, up=up)
        point.append(rng.uniform(0, 10.0))
    per_row = max(2, min(cols, 8))
    for i in range(rows):
        terms = {f"x{j}": float(rng.randint(-9, 9) or 1) for j in rng.sample(range(cols), per_row)}
        act = sum(v * point[int(c[1:])] for c, v in terms.items())
        slack = rng.choice([0.0, 0.0, rng.uniform(0, 20)])
        kind = rng.random()
        if kind < 0.4:
            lp.row(f"r{i}", terms, None, act + slack)
        elif kind < 0.8:
            lp.row(f"r{i}", terms, act - slack, None)
        else:
            lp.row(f"r{i}", terms, act, act)
    return lp


def refinery_scheduling(crudes: int, products: int, units: int, periods: int,
                        rng: random.Random) -> MpsWriter:
    """Refinery planning with cargo lots and unit commitment (on/off, startups)."""
    lp = MpsWriter(f"REFSCHED_{crudes}c_{products}p_{units}u_{periods}t", maximize=True)
    price = [rng.uniform(60, 120) for _ in range(products)]
    crude_cost = [rng.uniform(40, 80) for _ in range(crudes)]
    sulfur = [rng.uniform(0.1, 3.5) for _ in range(crudes)]
    max_sulfur = 1.8
    lot = [rng.choice([40.0, 50.0, 60.0]) for _ in range(crudes)]
    cap = [rng.uniform(250, 450) for _ in range(units)]
    min_frac = [rng.uniform(0.3, 0.5) for _ in range(units)]
    fixed = [rng.uniform(800, 2000) for _ in range(units)]
    startup = [rng.uniform(1500, 4000) for _ in range(units)]
    # Yields depend on the unit (deeper conversion units yield more of the light products).
    yields = []
    for u in range(units):
        per_crude = []
        for _ in range(crudes):
            w = [rng.random() * (1 + 0.3 * u * (p < products // 2)) for p in range(products)]
            total = sum(w)
            per_crude.append([0.96 * v / total for v in w])
        yields.append(per_crude)
    for t in range(periods):
        for k in range(crudes):
            lp.col(f"lots_{k}_{t}", cost=-crude_cost[k] * lot[k] * rng.uniform(0.9, 1.1), up=3,
                   integer=True)
            lp.col(f"cinv_{k}_{t}", cost=-0.5, up=300.0)
            for u in range(units):
                lp.col(f"proc_{k}_{u}_{t}")
        for u in range(units):
            lp.col(f"on_{u}_{t}", cost=-fixed[u], up=1, integer=True)
            lp.col(f"start_{u}_{t}", cost=-startup[u], up=1, integer=True)
        for p in range(products):
            lp.col(f"make_{p}_{t}")
            lp.col(f"sell_{p}_{t}", cost=price[p] * rng.uniform(0.9, 1.1), up=rng.uniform(60, 220))
            lp.col(f"pinv_{p}_{t}", cost=-0.8, up=200.0)
    for t in range(periods):
        for k in range(crudes):
            terms = {f"lots_{k}_{t}": lot[k], f"cinv_{k}_{t}": -1.0}
            for u in range(units):
                terms[f"proc_{k}_{u}_{t}"] = -1.0
            if t > 0:
                terms[f"cinv_{k}_{t - 1}"] = 1.0
            lp.row(f"cbal_{k}_{t}", terms, 0.0, 0.0)
        for u in range(units):
            throughput = {f"proc_{k}_{u}_{t}": 1.0 for k in range(crudes)}
            lp.row(f"cap_{u}_{t}", {**throughput, f"on_{u}_{t}": -cap[u]}, None, 0.0)
            lp.row(f"min_{u}_{t}", {**throughput, f"on_{u}_{t}": -min_frac[u] * cap[u]}, 0.0, None)
            terms = {f"start_{u}_{t}": 1.0, f"on_{u}_{t}": -1.0}
            if t > 0:
                terms[f"on_{u}_{t - 1}"] = 1.0
            lp.row(f"startup_{u}_{t}", terms, 0.0 if t > 0 else -1.0, None)
            lp.row(f"sulfur_{u}_{t}", {f"proc_{k}_{u}_{t}": sulfur[k] - max_sulfur
                                       for k in range(crudes)}, None, 0.0)
        for p in range(products):
            terms = {f"make_{p}_{t}": 1.0}
            for k in range(crudes):
                for u in range(units):
                    terms[f"proc_{k}_{u}_{t}"] = -yields[u][k][p]
            lp.row(f"yield_{p}_{t}", terms, 0.0, 0.0)
            terms = {f"make_{p}_{t}": 1.0, f"sell_{p}_{t}": -1.0, f"pinv_{p}_{t}": -1.0}
            if t > 0:
                terms[f"pinv_{p}_{t - 1}"] = 1.0
            lp.row(f"pbal_{p}_{t}", terms, 0.0, 0.0)
    return lp


def knapsack(items: int, dims: int, rng: random.Random) -> MpsWriter:
    """Multi-dimensional 0-1 knapsack; profits correlate with the weights."""
    lp = MpsWriter(f"KNAPSACK_{items}x{dims}", maximize=True)
    weights = [[rng.randint(1, 100) for _ in range(items)] for _ in range(dims)]
    for j in range(items):
        avg = sum(weights[d][j] for d in range(dims)) / dims
        lp.col(f"x{j}", cost=round(avg + rng.uniform(0, 20), 2), up=1, integer=True)
    for d in range(dims):
        capacity = round(0.5 * sum(weights[d]))
        lp.row(f"cap{d}", {f"x{j}": weights[d][j] for j in range(items)}, None, capacity)
    return lp


def facility(facilities: int, customers: int, rng: random.Random) -> MpsWriter:
    """Capacitated facility location with single-period demand."""
    lp = MpsWriter(f"FACILITY_{facilities}x{customers}")
    fpos = [(rng.random(), rng.random()) for _ in range(facilities)]
    cpos = [(rng.random(), rng.random()) for _ in range(customers)]
    demand = [rng.randint(5, 35) for _ in range(customers)]
    total = sum(demand)
    capacity = [rng.randint(int(3 * total / facilities), int(6 * total / facilities))
                for _ in range(facilities)]
    for f in range(facilities):
        lp.col(f"open{f}", cost=rng.randint(300, 900), up=1, integer=True)
    for f in range(facilities):
        for c in range(customers):
            dist = ((fpos[f][0] - cpos[c][0]) ** 2 + (fpos[f][1] - cpos[c][1]) ** 2) ** 0.5
            lp.col(f"x{f}_{c}", cost=round(10 * dist * demand[c], 3), up=1.0)
    for c in range(customers):
        lp.row(f"demand{c}", {f"x{f}_{c}": 1.0 for f in range(facilities)}, 1.0, 1.0)
    for f in range(facilities):
        terms = {f"x{f}_{c}": float(demand[c]) for c in range(customers)}
        terms[f"open{f}"] = -float(capacity[f])
        lp.row(f"cap{f}", terms, None, 0.0)
        for c in range(customers):
            lp.row(f"link{f}_{c}", {f"x{f}_{c}": 1.0, f"open{f}": -1.0}, None, 0.0)
    return lp


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", choices=["lp", "mip"], default="lp", help="instance set")
    parser.add_argument("--out", type=Path, help="output directory (default "
                        "bench/instances/generated or generated-mip)")
    parser.add_argument("--scale", type=int, default=1, help="size multiplier")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    if args.out is None:
        name = "generated" if args.set == "lp" else "generated-mip"
        args.out = Path(__file__).parent / "instances" / name
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    s = args.scale
    if args.set == "lp":
        instances = [
            transport(20 * s, 30 * s, rng),
            transport(40 * s, 60 * s, rng),
            refinery(10, 6, 12 * s, rng),
            refinery(20, 10, 52 * s, rng),
            sparse(500 * s, 800 * s, rng),
            sparse(1500 * s, 2500 * s, rng),
        ]
    else:
        instances = [
            refinery_scheduling(6, 4, 2, 6 * s, rng),
            refinery_scheduling(6, 6, 3, 12 * s, rng),
            refinery_scheduling(8, 6, 3, 26 * s, rng),
            knapsack(40 * s, 3, rng),
            knapsack(80 * s, 5, rng),
            facility(10 * s, 40 * s, rng),
            facility(20 * s, 80 * s, rng),
        ]
    for lp in instances:
        path = args.out / f"{lp.name.lower()}.mps"
        lp.write(path)
        print(f"wrote {path} ({len(lp.rows)} rows, {len(lp.cols)} cols)")


if __name__ == "__main__":
    main()
