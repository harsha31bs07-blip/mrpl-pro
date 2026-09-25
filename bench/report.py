#!/usr/bin/env python3
"""Writes a Markdown table from one or more harness CSV files.

Example:
    bench/report.py bench/results/generated-mip.csv --title "Generated MILP set" \\
        --solu bench/instances/miplib_small/miplib2017.solu
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def shifted_geomean(values: list[float], shift: float = 10.0) -> float:
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def read_solu(path: Path) -> dict[str, tuple[str, float | None]]:
    known: dict[str, tuple[str, float | None]] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("=opt=", "=best=", "=inf=", "=unbd="):
            value = float(parts[2]) if len(parts) > 2 and parts[0] in ("=opt=", "=best=") else None
            known[parts[1]] = (parts[0].strip("="), value)
    return known


def fmt(value: str | None) -> str:
    if value in (None, "", "None"):
        return ""
    return f"{float(value):.10g}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("--title", default="Results")
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--solu", type=Path, help="MIPLIB .solu file with known optimal values")
    args = parser.parse_args()

    rows: dict[str, dict[str, dict]] = {}
    for path in args.csv:
        for r in csv.DictReader(path.open()):
            rows.setdefault(r["instance"], {})[r["solver"]] = r
    solvers = sorted({s for per in rows.values() for s in per}, key=lambda s: s != "samaya")
    known = read_solu(args.solu) if args.solu else {}

    out = [f"## {args.title}", ""]
    for solver in solvers:
        mine = [per[solver] for per in rows.values() if solver in per]
        solved = [r for r in mine if r["status"] in ("optimal", "infeasible", "unbounded")]
        times = [float(r["wall_seconds"]) if r in solved else args.time_limit for r in mine]
        out.append(f"* {solver}: solved {len(solved)}/{len(mine)}, shifted geometric mean "
                   f"{shifted_geomean(times):.2f} s (unsolved instances count at the "
                   f"{args.time_limit:g} s limit)")
    out.append("")
    header = "| instance |"
    rule = "|---|"
    if known:
        header += " known optimum |"
        rule += "---:|"
    for solver in solvers:
        header += f" {solver} status | {solver} objective | {solver} s |"
        rule += "---|---:|---:|"
    out += [header, rule]
    for name in sorted(rows):
        line = f"| {name} |"
        if known:
            kind, value = known.get(name, ("", None))
            line += f" {fmt(value) if kind == 'opt' else kind} |"
        for solver in solvers:
            r = rows[name].get(solver)
            if r is None:
                line += " | | |"
                continue
            line += (f" {r['status']} | {fmt(r['objective'])} |"
                     f" {float(r['wall_seconds']):.2f} |")
        out.append(line)
    print("\n".join(out))


if __name__ == "__main__":
    main()
