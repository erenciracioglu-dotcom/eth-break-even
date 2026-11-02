"""T25 self-check (round 10): with the bincross flag the search beats 2^(n/3) at lam = 652 (the search without it
needs 663)."""
import time
import appd_params as P
import appd_assembly_public as A

t0 = time.time()
p = P.choose_params(652, "third")
dt = time.time() - t0
s = A.exponent_of(p, 652)
print("valid", s["valid"], "e", s.get("e"), "bincross", p.get("bincross"), "r", p.get("r"), "seconds", round(dt, 2))
ok = bool(s["valid"]) and bool(p.get("bincross")) and s["e"] < 1.0 / 3.0 and dt < 20
print("T25 OK" if ok else "T25 FAIL")
