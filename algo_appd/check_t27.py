"""T27 self-check (round 11): with the profile flag the search beats 2^(n/3) at lam = 600 (the search without it
needs 646)."""
import time
import appd_params as P
import appd_assembly_public as A

t0 = time.time()
p = P.choose_params(600, "third")
dt = time.time() - t0
s = A.exponent_of(p, 600)
print("valid", s["valid"], "e", s.get("e"), "profile", p.get("profile"), "r", p.get("r"), "seconds", round(dt, 2))
ok = bool(s["valid"]) and bool(p.get("profile")) and s["e"] < 1.0 / 3.0 and dt < 20
print("T27 OK" if ok else "T27 FAIL")
