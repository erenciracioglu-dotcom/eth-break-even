"""T22 self-check: the improved search beats 2^(n/3) at lam = 670 (the current search needs 686)."""
import time
import appd_params as P
import appd_assembly_public as A

t0 = time.time()
p = P.choose_params(670, "third")
dt = time.time() - t0
s = A.exponent_of(p, 670)
print("valid", s["valid"], "e", s.get("e"), "seconds", round(dt, 2), "Lc", p.get("Lc"), "u", p.get("u"))
ok = bool(s["valid"]) and s["e"] < 1.0 / 3.0 and dt < 20
print("T22 OK" if ok else "T22 FAIL")
