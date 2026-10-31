"""Parameter choice for the Appendix D route (agent-editable). The referee calls choose_params(lam, target)
for every target in appd_assembly_public.TARGETS and re-scores each answer with its own frozen copy of the
assembly. USE_CODE_FLAGS switches on memo, modedp, chunks and leafvars, which only count once appd_core.py
passes the verify gate. Keep one call well under 3 s."""
import math
import appd_assembly_public as A

USE_CODE_FLAGS = True
R0 = (3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32, 48, 64)
R1 = (4, 6, 8, 12, 16)
RHO = (0.4, 0.6, 0.8)


def min_a(u, lp_in, M):
    """Smallest a with (2^lp_in + 4 M (a + 1)) (u - 1) / (4 M) <= 2^(a - u + 1) (the assembly's headroom test)."""
    a = max(u + 1, int(u - 1 + lp_in + math.log2(u - 1) - math.log2(4 * M)))
    while A.lse2(lp_in, math.log2(4 * M * (a + 1))) + math.log2(u - 1) - math.log2(4 * M) > a - u + 1:
        a += 1
    return a


def make(u0, r, rho, delta, Lc, flags):
    L = len(r)
    us = [u0]
    for d in range(1, L):
        us.append(max(2, int(rho * us[-1])))
    M = [(1.0 + 1.0 / Lc) if flags["chunks"] else 3.0]
    for d in range(L):
        M.append(M[-1] * 2 * (1 + delta) / r[d])
    a, lp_in = [], 0.0
    for d in range(L):
        a.append(min_a(us[d], lp_in, M[d]))
        h = M[d] / (us[d] - 1)
        lp_in = A.lse2(math.log2(6 * h) + a[d], lp_in, math.log2(4 * M[d] * (a[d] + 1)))
    return dict(levels=L, r=list(r), u=us, a=a, delta=delta, Lc=Lc, **flags)


def build(lam, r, rho, delta, Lc, flags):
    lo, hi = 2, max(3, int(lam))
    if not A.exponent_of(make(lo, r, rho, delta, Lc, flags), lam)["valid"]:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if A.exponent_of(make(mid, r, rho, delta, Lc, flags), lam)["valid"]:
            lo = mid
        else:
            hi = mid
    return make(lo, r, rho, delta, Lc, flags)


def choose_params(lam, target="third"):
    flags = dict(perclass=True, labelfree=True, memo=False, modedp=False, chunks=False, leafvars=False)
    if USE_CODE_FLAGS:
        flags.update(memo=True, modedp=True, chunks=True)
    Lc = 4096 if flags["chunks"] else 1
    best = None
    shapes = [(r0,) for r0 in R0] + [(r0, r1) for r0 in R0[:8] for r1 in R1]
    for shape in shapes:
        for rho in (RHO if len(shape) > 1 else (0.5,)):
            fl = dict(flags, leafvars=flags["chunks"] and len(shape) == 1)
            p = build(lam, shape, rho, 0.01, Lc, fl)
            if p is None:
                continue
            s = A.exponent_of(p, lam)
            if s["valid"] and (best is None or s["e"] < best[0]):
                best = (s["e"], p)
    return best[1] if best else dict(levels=1, r=[2], u=[2], a=[3], delta=0.01, Lc=1)
