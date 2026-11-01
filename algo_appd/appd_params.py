"""Parameter choice for the Appendix D route (agent-editable). The referee calls choose_params(lam, target)
for every target in appd_assembly_public.TARGETS and re-scores each answer with its own frozen copy of the
assembly. USE_CODE_FLAGS switches on memo, modedp, chunks and leafvars, which only count once appd_core.py
passes the verify gate. Keep one call well under 3 s.

Search (attempt 2): three levels allowed, free u_d per level (not rho * u_0), wider r grids, delta grid down
to 1e-4 and Lc up to 2^20. USE_CODE_FLAGS also switches on newton (appd-newton, the binomial mode basis), which
halves the mode-tuple constant cK in the replica (A.ck) and in the frozen assembly. Shapes are pre-ranked with a
cheap analytic replica of the assembly (closed-form greedy u from the padding budget + min_a), then the top
candidates are refined with the real model: bisection on u_0, and inside every trial a bisection that maximises
each remaining u_d. min_a is kept for a_d.
"""
import math
import appd_assembly_public as A

USE_CODE_FLAGS = True
LN2 = math.log(2.0)

# r grids: r_0 small (the root group count sets M_d / sqrt(u_d) everywhere), r_1 and r_2 wide and geometric
R0 = (4, 5, 6, 7, 8, 9, 10, 12, 16)
R1 = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 64, 96, 128, 192, 256, 384, 512)
R2 = (2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32, 64, 128)
R1_RANK = R1[:10]          # cheap pre-rank grid for level 1 of a 3-level shape
R2_RANK = R2[:10]          # cheap pre-rank grid for level 2
DELTAS = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)
LCS = (4096, 1 << 20)
TOPK = 10                  # candidates refined with the real assembly
UPMAX = 2 * 10 ** 6        # u_d < a_d <= 1e7


def min_a(u, lp_in, M):
    """Smallest a with (2^lp_in + 4 M (a + 1)) (u - 1) / (4 M) <= 2^(a - u + 1) (the assembly's headroom test)."""
    a = max(u + 1, int(u - 1 + lp_in + math.log2(u - 1) - math.log2(4 * M)))
    while A.lse2(lp_in, math.log2(4 * M * (a + 1))) + math.log2(u - 1) - math.log2(4 * M) > a - u + 1:
        a += 1
    return a


def _M(r, delta, Lc, chunks):
    M = [(1.0 + 1.0 / Lc) if chunks else 3.0]
    for d in range(len(r)):
        M.append(M[-1] * 2 * (1 + delta) / r[d])
    return M


def est_u(rd, Md, delta, lam, lp_in):
    """Closed-form greedy largest u_d allowed by the padding budget at depth d (a_d = min_a)."""
    P = math.log2(rd * (rd + 1) / 2.0)
    u = int((lam + math.log2(Md) + math.log2(delta)) / (P + 1.0))
    if u < 2:
        return 2
    for _ in range(6):
        a = min_a(u, lp_in, Md)
        nu = int((lam + math.log2(Md) + math.log2(delta) - a) / P)
        if nu >= u:
            break
        if nu < 2:
            return 2
        u = nu
    return max(2, u)


def rank_e(r, delta, Lc, lam, leafvars, cK=2.0):
    """Analytic replica: greedy u_d from the padding budget, a_d = min_a, cost exponent of the assembly."""
    L = len(r)
    M = _M(r, delta, Lc, True)
    dm = 2 * delta
    u, lp_in, cum, worst = [], 0.0, 0.0, 0.0
    for d in range(L):
        ud = est_u(r[d], M[d], delta, lam, lp_in)
        ad = min_a(ud, lp_in, M[d])
        u.append(ud)
        x = M[d] / math.sqrt(ud)
        ar = A.tup(ud, A.eps_dp(r[d]), dm, cK) * x
        ch = A.tup(ud, A.f_grp(r[d]), dm, cK) * x
        worst = max(worst, cum + ar)
        cum += ch
        h = M[d] / (ud - 1)
        lp_in = A.lse2(math.log2(6 * h) + ad, lp_in, math.log2(4 * M[d] * (ad + 1)))
    leaf = M[L]
    if leafvars:
        leaf = min(leaf, (1 + delta) * (1 + M[0] / 2) / r[0])
    return max(worst, cum + leaf), u


def make(u, r, delta, Lc, flags):
    """Params with a_d = min_a (same contract as the frozen headroom test)."""
    L = len(r)
    M = _M(r, delta, Lc, flags["chunks"])
    a, lp_in = [], 0.0
    for d in range(L):
        a.append(min_a(u[d], lp_in, M[d]))
        h = M[d] / (u[d] - 1)
        lp_in = A.lse2(math.log2(6 * h) + a[d], lp_in, math.log2(4 * M[d] * (a[d] + 1)))
    return dict(levels=L, r=list(r), u=list(u), a=a, delta=delta, Lc=Lc, **flags)


def build(r, delta, Lc, flags, lam, ustart=None):
    """Bisection on u_0; every trial maximises u_1..u_{L-1} (free per level) by bisection. min_a for a_d."""
    L = len(r)
    M = _M(r, delta, Lc, flags["chunks"])

    def score(u):
        return A.exponent_of(make(u, r, delta, Lc, flags), lam)

    def max_up(u, d, hi):
        lo = u[d]

        def ok(x):
            uu = list(u)
            uu[d] = x
            return score(uu)["valid"]

        if not ok(lo):
            return None
        if ok(hi):
            return hi
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ok(mid):
                lo = mid
            else:
                hi = mid
        return lo

    lp, u = 0.0, list(ustart) if ustart else [2] * L
    for d in range(L):
        u[d] = min(u[d], est_u(r[d], M[d], delta, lam, lp))
        lp = A.lse2(math.log2(6 * (M[d] / (u[d] - 1))) + min_a(u[d], lp, M[d]), lp,
                    math.log2(4 * M[d] * (min_a(u[d], lp, M[d]) + 1)))
    for d in range(L):
        m = max_up(u, d, UPMAX)
        if m is None:
            return None
        u[d] = m

    def full_ok(u0):
        uu = list(u)
        uu[0] = u0
        lp = A.lse2(math.log2(6 * (M[0] / (u0 - 1))) + min_a(u0, 0.0, M[0]), 0.0,
                    math.log2(4 * M[0] * (min_a(u0, 0.0, M[0]) + 1)))
        for d in range(1, L):
            uu[d] = min(max(2, est_u(r[d], M[d], delta, lam, lp)), UPMAX)
            ad = min_a(uu[d], lp, M[d])
            lp = A.lse2(math.log2(6 * (M[d] / (uu[d] - 1))) + ad, lp, math.log2(4 * M[d] * (ad + 1)))
        for d in range(1, L):
            m = max_up(uu, d, UPMAX)
            if m is None:
                return None
            uu[d] = m
        return uu if score(uu)["valid"] else None

    lo = u[0]
    hi = max(lo + 1, int(2 * u[0]))
    trial = full_ok(hi)
    guard = 0
    while trial is not None and guard < 40:
        lo = hi
        hi *= 2
        trial = full_ok(hi)
        guard += 1
    if full_ok(lo) is None:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if full_ok(mid) is not None:
            lo = mid
        else:
            hi = mid
    uu = full_ok(lo)
    if uu is None:
        return None
    for _ in range(2):                                   # coordinate ascent polish on all levels
        for d in range(L):
            m = max_up(uu, d, UPMAX)
            if m is None:
                break
            uu[d] = m
    s = score(uu)
    return (s["e"], make(uu, r, delta, Lc, flags)) if s["valid"] else None


def choose_params(lam, target="third"):
    flags = dict(perclass=True, labelfree=True, memo=False, modedp=False, chunks=False, leafvars=False,
                 newton=False, ldu=False)
    if USE_CODE_FLAGS:
        flags.update(memo=True, modedp=True, chunks=True, newton=True, ldu=True)
    lam = float(lam)
    shapes = []
    for r0 in R0:
        shapes.append((r0,))
    for r0 in R0:
        for r1 in R1:
            shapes.append((r0, r1))
    for r0 in R0[:6]:
        for r1 in R1_RANK:
            for r2 in R2_RANK:
                shapes.append((r0, r1, r2))
    cand = []
    cK = A.ck(flags)
    for shape in shapes:
        lv = flags["chunks"] and len(shape) == 1
        for delta in DELTAS:
            for Lc in LCS:
                e, _ = rank_e(list(shape), delta, Lc, lam, lv, cK)
                cand.append((e, shape, delta, Lc))
    cand.sort(key=lambda x: x[0])
    best = None
    seen = set()
    for _, shape, delta, Lc in cand[:TOPK]:
        key = (shape, delta, Lc)
        if key in seen:
            continue
        seen.add(key)
        fl = dict(flags, leafvars=flags["chunks"] and len(shape) == 1)
        p = build(list(shape), delta, Lc, fl, lam)
        if p is None:
            continue
        s = A.exponent_of(p[1], lam)
        if s["valid"] and (best is None or s["e"] < best[0]):
            best = (s["e"], p[1])
    if best:
        return best[1]
    return dict(levels=1, r=[2], u=[2], a=[3], delta=1e-4, Lc=4096, **flags)
