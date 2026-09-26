#!/usr/bin/env python3
"""Report for bench/ab_bound.sh: solved counts, shifted geometric mean time (shift 10 s, unsolved
at the time limit) and, per instance, the distance of the proven bound from the known optimum
(|optimum - bound| / max(|optimum|, |bound|), 0 when solved, 100% without a finite bound).

    python3 bench/ab_report.py OUT [SOLU] [SECONDS]
"""
import json, glob, sys, math, os, statistics
d = sys.argv[1]
limit = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0  # The runs' time limit.
solu = {}
for l in open(sys.argv[2] if len(sys.argv) > 2 else 'bench/instances/miplib_small/miplib2017.solu'):
    p = l.split()
    if len(p) >= 3 and p[0] in ('=opt=', '=best='): solu[p[1]] = float(p[2])
names = sorted({os.path.basename(f).rsplit('.', 2)[0] for f in glob.glob(d + '/*.a.json')})
def load(n, s):
    try: return json.load(open(f'{d}/{n}.{s}.json'))
    except Exception: return None
def dgap(r, n):  # distance of the proven bound from the known optimum, relative
    if r is None or n not in solu: return None
    b = r.get('dual_bound')
    if r['status'] == 'optimal': return 0.0
    if b is None or not isinstance(b, (int, float)) or not math.isfinite(b): return 1.0
    o = solu[n]
    return min(1.0, abs(o - b) / max(abs(o), abs(b), 1e-9))
rows = []
solved = {'a': 0, 'b': 0}; times = {'a': [], 'b': []}
better = worse = 0; ga = []; gb = []
for n in names:
    a, b = load(n, 'a'), load(n, 'b')
    if a is None or b is None: continue
    for s, r in (('a', a), ('b', b)):
        ok = r['status'] in ('optimal', 'infeasible')
        solved[s] += ok
        times[s].append(r['solve_seconds'] if ok else limit)
    x, y = dgap(a, n), dgap(b, n)
    if x is None: continue
    ga.append(x); gb.append(y)
    if y < x - 1e-4: better += 1
    elif y > x + 1e-4: worse += 1
    if abs(x - y) > 1e-4: rows.append((n, x, y, a['nodes'], b['nodes']))
sgm = lambda t: math.exp(sum(math.log(v + 10) for v in t) / len(t)) - 10
print(f"instances {len(ga)}; solved A {solved['a']} B {solved['b']}; shifted geomean A {sgm(times['a']):.2f} B {sgm(times['b']):.2f}")
print(f"bound gap to the optimum: mean A {statistics.mean(ga):.2%} B {statistics.mean(gb):.2%}; median A {statistics.median(ga):.2%} B {statistics.median(gb):.2%}; B better {better}, worse {worse}")
for n, x, y, na, nb in sorted(rows, key=lambda r: r[2] - r[1]):
    print(f"  {n:24} {x:8.3%} -> {y:8.3%}   nodes {na} -> {nb}")
