#!/usr/bin/env python3
"""Generates structured LP instances as free-format MPS files.

Families (all feasible and bounded by construction):
  transport   balanced transportation problem, S sources x D sinks
  refinery    multi-period refinery planning: crude purchase, processing, product yields,
              inventories, CDU capacity and a sulfur blending constraint per period
  sparse      random sparse LP with row bounds built around a feasible point

Example:
    bench/generate_lps.py --out bench/instances/generated --scale 1
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

    def col(self, name: str, cost: float = 0.0, lo: float | None = 0.0, up: float | None = None):
        self.cols[name] = {}
        self.obj[name] = cost
        self.bounds[name] = (lo, up)

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
        for c, entries in self.cols.items():
            if self.obj[c] != 0.0:
                lines.append(f" {c} obj {self.obj[c]:.12g}")
            for r, v in entries.items():
                lines.append(f" {c} {r} {v:.12g}")
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "instances" / "generated")
    parser.add_argument("--scale", type=int, default=1, help="size multiplier")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    s = args.scale
    instances = [
        transport(20 * s, 30 * s, rng),
        transport(40 * s, 60 * s, rng),
        refinery(10, 6, 12 * s, rng),
        refinery(20, 10, 52 * s, rng),
        sparse(500 * s, 800 * s, rng),
        sparse(1500 * s, 2500 * s, rng),
    ]
    for lp in instances:
        path = args.out / f"{lp.name.lower()}.mps"
        lp.write(path)
        print(f"wrote {path} ({len(lp.rows)} rows, {len(lp.cols)} cols)")


if __name__ == "__main__":
    main()
