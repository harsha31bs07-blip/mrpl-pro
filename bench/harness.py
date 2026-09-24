#!/usr/bin/env python3
"""Runs samaya (and optionally baseline solvers) over a set of instances and writes a CSV.

Examples:
    bench/harness.py bench/instances/netlib --time-limit 60
    bench/harness.py bench/instances/netlib --baseline highs --out bench/results/netlib.csv

Baselines are external executables used only for comparison; they are never linked into samaya.
The summary reports solved counts and the shifted geometric mean of solve time (shift 10 s), the
standard metric in the Mittelmann benchmarks.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_SAMAYA = REPO / "build" / "release" / "apps" / "cli" / "samaya"
SOLVED = {"optimal", "infeasible", "unbounded"}


def run_samaya(exe: Path, instance: Path, time_limit: float, threads: int) -> dict:
    cmd = [str(exe), "--json", "--log-level", "0", "--time-limit", str(time_limit),
           "--threads", str(threads), str(instance)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "wall_seconds": time.perf_counter() - start}
    wall = time.perf_counter() - start
    lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if not lines:
        return {"status": "crash", "wall_seconds": wall, "message": proc.stderr.strip()[-200:]}
    result = json.loads(lines[-1])
    result["wall_seconds"] = wall
    return result


def run_highs(exe: str, instance: Path, time_limit: float, threads: int) -> dict:
    cmd = [exe, "--time_limit", str(time_limit), "--model_file", str(instance)]
    start = time.perf_counter()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=time_limit + 60)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "wall_seconds": time.perf_counter() - start}
    wall = time.perf_counter() - start
    out = proc.stdout
    status_match = re.search(r"Model status\s*:\s*(.+)", out)
    obj_match = re.search(r"Objective value\s*:\s*(\S+)", out)
    raw = status_match.group(1).strip().lower() if status_match else "unknown"
    status = {"optimal": "optimal", "infeasible": "infeasible", "unbounded": "unbounded",
              "time limit reached": "time_limit"}.get(raw, raw.replace(" ", "_"))
    return {"status": status, "objective": float(obj_match.group(1)) if obj_match else None,
            "solve_seconds": wall, "wall_seconds": wall}


def run_highspy(instance: Path, time_limit: float, threads: int) -> dict:
    import highspy  # Optional dependency, used only as a comparison baseline.

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(time_limit))
    h.setOptionValue("threads", int(threads))
    h.readModel(str(instance))
    start = time.perf_counter()
    h.run()
    wall = time.perf_counter() - start
    raw = h.modelStatusToString(h.getModelStatus()).lower()
    status = {"optimal": "optimal", "infeasible": "infeasible", "unbounded": "unbounded",
              "time limit reached": "time_limit"}.get(raw, raw.replace(" ", "_"))
    objective = h.getInfo().objective_function_value if status == "optimal" else None
    return {"status": status, "objective": objective, "solve_seconds": wall, "wall_seconds": wall}


def shifted_geomean(values: list[float], shift: float = 10.0) -> float:
    if not values:
        return float("nan")
    return math.exp(sum(math.log(v + shift) for v in values) / len(values)) - shift


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("instances", nargs="+", type=Path,
                        help="instance files or directories containing *.mps / *.qps")
    parser.add_argument("--samaya", type=Path, default=DEFAULT_SAMAYA)
    parser.add_argument("--baseline", action="append", default=[], choices=["highs", "highspy"],
                        help="also run a baseline solver: the highs executable on PATH, or the "
                             "highspy Python module")
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out", type=Path, default=REPO / "bench" / "results" / "results.csv")
    args = parser.parse_args()

    files: list[Path] = []
    for path in args.instances:
        if path.is_dir():
            files += sorted(p for p in path.iterdir() if p.suffix.lower() in {".mps", ".qps"})
        else:
            files.append(path)
    if not files:
        print("no instances found", file=sys.stderr)
        return 1
    if not args.samaya.exists():
        print(f"samaya binary not found at {args.samaya}; build with "
              "`cmake --preset release && cmake --build --preset release`", file=sys.stderr)
        return 1

    solvers = {"samaya": lambda f: run_samaya(args.samaya, f, args.time_limit, args.threads)}
    for name in args.baseline:
        if name == "highspy":
            solvers[name] = lambda f: run_highspy(f, args.time_limit, args.threads)
            continue
        exe = shutil.which(name)
        if exe is None:
            print(f"baseline '{name}' not found on PATH", file=sys.stderr)
            return 1
        solvers[name] = lambda f, exe=exe: run_highs(exe, f, args.time_limit, args.threads)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for instance in files:
        for solver, run in solvers.items():
            r = run(instance)
            row = {"instance": instance.stem, "solver": solver, "status": r.get("status"),
                   "objective": r.get("objective"), "solve_seconds": r.get("solve_seconds"),
                   "wall_seconds": round(r["wall_seconds"], 4), "nodes": r.get("nodes")}
            rows.append(row)
            print(f"{instance.stem:<20} {solver:<8} {row['status']:<16} obj={row['objective']}"
                  f" t={row['wall_seconds']:.3f}s", flush=True)

    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {args.out}")
    if len(solvers) > 1:
        by_instance: dict[str, dict[str, dict]] = {}
        for r in rows:
            by_instance.setdefault(r["instance"], {})[r["solver"]] = r
        disagreements = 0
        for inst, results in by_instance.items():
            base = [r for s, r in results.items() if s != "samaya"]
            mine = results.get("samaya")
            for other in base:
                same_status = mine["status"] == other["status"]
                same_obj = True
                if same_status and mine["status"] == "optimal":
                    a, b = mine["objective"], other["objective"]
                    same_obj = abs(a - b) <= 1e-6 * (1 + abs(b))
                if not (same_status and same_obj):
                    disagreements += 1
                    print(f"  DISAGREE {inst}: samaya {mine['status']} {mine['objective']} vs "
                          f"{other['solver']} {other['status']} {other['objective']}")
        print(f"Status/objective disagreements with baselines: {disagreements}")
    for solver in solvers:
        mine = [r for r in rows if r["solver"] == solver]
        solved = [r for r in mine if r["status"] in SOLVED]
        # Unsolved instances count at the time limit, as in the Mittelmann benchmarks.
        times = [r["wall_seconds"] if r["status"] in SOLVED else args.time_limit for r in mine]
        print(f"{solver:<8} solved {len(solved)}/{len(mine)}, "
              f"shifted geomean {shifted_geomean(times):.3f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
