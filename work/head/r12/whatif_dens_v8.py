"""Round 12: general clause density rho = m/n with the frozen v8 assembly (matchbal) and the head's round-12 search,
with the two corrections of GPT-6 Astra's review:
1. Small connected components are solved directly. With no isolated variables a connected component with c clauses
   of width <= 3 has at most 2c + 1 variables, so 2c - v >= -1 and the small components remove at most n_s/3 from
   2m - n: the large components carry at most 2m - 2n/3 + n/Lc wires. M_root = 2 rho - 2/3 + 3 rho/(B - 2) + 1/Lc
   (the last-but-one term for the chains of copies, as before).
2. Root balance by matchings with the exact unit degree for occurrence cap B: a unit has at most (2B + 1) Lc
   variables (variable-graph degree <= 2B), at most B (2B + 1) Lc designated clauses with two outside incidences
   each and B (2B + 1) Lc outside incidences of its own variables, so Delta <= 3 B (2B + 1) Lc and
   q = 2 Delta - 1 (B = 3 gives 126 Lc - 1). This replaces the old scaling of the two-stage constant by (B/3)^4.
Usage: python whatif_dens_v8.py <rho> <out> [targets]"""
import json
import math
import os
import sys
import time
import types

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
RHO = float(sys.argv[1])
OUT = sys.argv[2]
KEYS = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] else None
BCAP = 1024.0
XC = 3.0 * RHO / (BCAP - 2.0)
QEXPR = "(2.0 * 3.0 * %r * (2.0 * %r + 1.0) * Lc - 1.0)" % (BCAP, BCAP)


def load(path, name, repl):
    src = open(path, encoding="utf-8").read().replace("\r\n", "\n")
    for old, new in repl:
        assert src.count(old) == 1, (path, old, src.count(old))
        src = src.replace(old, new)
    mod = types.ModuleType(name)
    mod.__file__ = path
    sys.modules[name] = mod
    exec(compile(src, path, "exec"), mod.__dict__)
    return mod


m0 = "(2.0 * %r - 2.0 / 3.0 + %r + 1.0 / Lc) if fl[\"chunks\"] else (6.0 * %r - 3.0)" % (RHO, XC, RHO)
A = load(os.path.join(ROOT, "referee_appd", "appd_assembly.py"), "appd_assembly_public", [
    ("M = [(1.0 + 1.0 / Lc) if fl[\"chunks\"] else 3.0]", "M = [" + m0 + "]"),
    ("q = 126.0 * Lc - 1.0", "q = " + QEXPR),
    ("leaf = min(leaf, (1 + delta) * (1 + M[0] / 2) / r[0])", "leaf = min(leaf, (1 + delta) * (1 + %r + M[0] / 2) / r[0])" % XC),
])
m1 = "(2.0 * %r - 2.0 / 3.0 + %r + 1.0 / Lc) if chunks else (6.0 * %r - 3.0)" % (RHO, XC, RHO)
P = load(os.path.join(ROOT, "algo_appd", "appd_params.py"), "appd_params", [
    ("M = [(1.0 + 1.0 / Lc) if chunks else 3.0]", "M = [" + m1 + "]"),
    ("leaf = min(leaf, (1 + delta) * (1 + M[0] / 2) / r[0])", "leaf = min(leaf, (1 + delta) * (1 + %r + M[0] / 2) / r[0])" % XC),
    ("        q = 126.0 * Lc - 1.0\n", "        q = " + QEXPR + "\n"),
])


def ev(lam, key):
    p = P.choose_params(int(lam), key)
    s = A.exponent_of(p, int(lam))
    return (s["e"], p) if s["valid"] else (9.0, p)


res = {"rho": RHO, "bcap": BCAP, "M_root": "2 rho - 2/3 + 3 rho/(B-2) + 1/Lc", "q": QEXPR, "thresholds": {}, "params": {}}
t0 = time.time()
for key, target in A.TARGETS.items():
    if KEYS and key not in KEYS:
        continue
    lo, hi = 64, 10 ** 9
    e_hi, p_hi = ev(hi, key)
    if e_hi >= target:
        res["thresholds"][key] = None
        continue
    while hi - lo > 1:
        mid = int(math.sqrt(lo * hi)) if hi / lo > 1.01 else (lo + hi) // 2
        if mid in (lo, hi):
            mid = (lo + hi) // 2
        e, p = ev(mid, key)
        if e < target:
            hi, p_hi, e_hi = mid, p, e
        else:
            lo = mid
    res["thresholds"][key] = hi
    res["params"][key] = {k: p_hi.get(k) for k in ("levels", "r", "u", "a", "delta", "Lc", "kappa", "matchbal")}
    res["params"][key]["e"] = e_hi
    print(key, hi, res["params"][key], flush=True)
res["seconds"] = round(time.time() - t0, 1)
json.dump(res, open(OUT, "w"), indent=1)
print("DONE", json.dumps(res["thresholds"]), flush=True)
