"""Round 12 (head): the referee's threshold search (appd_gates.score, bisection from 64 to 1e9) on a candidate
directory with a given assembly, all code flags taken as verified. Usage: python scan12.py <cand_dir> <assembly.py> <out.json> [keys]"""
import importlib.util
import json
import math
import sys
import time

cand, asm, out = sys.argv[1], sys.argv[2], sys.argv[3]
keys = sys.argv[4].split(",") if len(sys.argv) > 4 else None
spec = importlib.util.spec_from_file_location("frozen_scan", asm)
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
sys.path.insert(0, cand)
import appd_params as P  # noqa: E402

calls, spent = 0, 0.0


def ev(lam, key):
    global calls, spent
    t0 = time.time()
    p = P.choose_params(int(lam), key)
    spent += time.time() - t0
    calls += 1
    s = A.exponent_of(p, int(lam))
    return (s["e"], p) if s["valid"] else (9.0, p)


res = {"thresholds": {}, "certificates": {}}
for key, target in A.TARGETS.items():
    if keys and key not in keys:
        continue
    lo, hi = 64, 10 ** 9
    e_hi, p_hi = ev(hi, key)
    if e_hi >= target:
        res["thresholds"][key] = None
        continue
    e_lo, _ = ev(lo, key)
    if e_lo < target:
        res["thresholds"][key] = lo
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
    res["certificates"][key] = {"lam": hi, "params": p_hi, "e": e_hi}
    print(key, hi, {k: p_hi.get(k) for k in ("r", "u", "a", "delta", "Lc", "kappa")}, e_hi, flush=True)
res["avg_s"] = spent / max(calls, 1)
res["calls"] = calls
json.dump(res, open(out, "w"), indent=1)
print("DONE", json.dumps(res["thresholds"]), "avg %.2f s" % res["avg_s"], flush=True)
