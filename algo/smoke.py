"""Quick public checks, about one minute.

Run from this folder with:  D:/ETH-loop/.venv/Scripts/python.exe smoke.py
It checks one construction, one exact count, and prints the cost exponent e at three sizes.
A lower e at lam = 2224173 means a better parameter choice.
"""
import sys
import time
import numpy as np
from construction import Aux, random_formula, build_design
from counting import Counter, brute_count
from assembly_public import exponent_of, TARGETS
from model import choose_params

t0 = time.time()
rep = Aux(2, 2, 2, 0).report()
keys = ("f_interp", "eigen", "orthogonal", "support", "P0_delta", "degree", "separation")
if not all(rep[k] for k in keys):
    print("SMOKE FAIL: construction checks", {k: rep[k] for k in keys})
    sys.exit(1)

rng = np.random.default_rng(11)
clauses = random_formula(5, 7, rng)
aux = Aux(2, 2, 2, 0)
des = build_design(aux, 5, clauses, rng)
C = Counter(aux, des)
data = {t: C.taylor_inclusion_exclusion(t) for t in range(1, aux.Q + 1)}
h0, bc = int(C.hermite_H0(data)), brute_count(5, clauses)
if h0 != bc:
    print(f"SMOKE FAIL: count {h0} != brute force {bc}")
    sys.exit(1)

for lam in (300000, 2224173, 5000000):
    t1 = time.time()
    p = choose_params(lam)
    s = exponent_of(p, lam)
    dt = time.time() - t1
    if not s["valid"]:
        print(f"SMOKE FAIL: invalid params at lam={lam}: {p} ({s['reason']})")
        sys.exit(1)
    if dt > 3.0:
        print(f"SMOKE FAIL: choose_params took {dt:.1f}s at lam={lam} (limit 3s)")
        sys.exit(1)
    print(f"lam={lam}: e={s['e']:.5f} params={p} ({dt:.2f}s)")
print(f"target e < {TARGETS['ppsz']:.5f};  smoke ok in {time.time() - t0:.1f}s")
