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

Search (attempt 3, round 9 / T22): est_u is now the smaller of the padding budget (est_u_pad) and the frozen
assembly's label-balance test (bal_u), which often binds first; the depth d, the previous level's a (0 at the
root) and Lc are threaded into every est_u call, and the grids are finer: 8 deltas down to 1e-4, 5 chunk sizes
from 64 to 2^20 and TOPK = 20 refined candidates.
"""
import math
import appd_assembly_public as A

USE_CODE_FLAGS = True
LN2 = math.log(2.0)

GUARD_KAPPA = 10           # appd-guarded: guard width H_d = 2^(a_d - GUARD_KAPPA) used by the frozen model
_KAPPA = None              # set by choose_params (None = the paper's Cauchy deformation)

# r grids: r_0 small (the root group count sets M_d / sqrt(u_d) everywhere), r_1 and r_2 wide and geometric
R0 = (4, 5, 6, 7, 8, 9, 10, 12, 16)
R1 = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 64, 96, 128, 192, 256, 384, 512)
R2 = (2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32, 64, 128)
R1_RANK = R1[:10]          # cheap pre-rank grid for level 1 of a 3-level shape
R2_RANK = R2[:10]          # cheap pre-rank grid for level 2
DELTAS = (1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 3e-2)
LCS = (64, 256, 1024, 4096, 1 << 20)
TOPK = 20                  # candidates refined with the real assembly
UPMAX = 2 * 10 ** 6        # u_d < a_d <= 1e7


def min_a(u, lp_in, M):
    """Smallest a with (2^lp_in + 4 M (a + 1)) (u - 1) / (4 M) <= 2^(a - u + 1) (the assembly's headroom test)."""
    a = max(u + 1, int(u - 1 + lp_in + math.log2(u - 1) - math.log2(4 * M)))
    a = max(a, u + (_KAPPA or 0))
    while A.need_bits(lp_in, M, a, u, _KAPPA) > a - u + 1:
        a += 1
    return a


def _M(r, delta, Lc, chunks):
    M = [(1.0 + 1.0 / Lc) if chunks else 3.0]
    for d in range(len(r)):
        M.append(M[-1] * 2 * (1 + delta) / r[d])
    return M


def est_u_pad(rd, Md, delta, lam, lp_in):
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


def bal_u(rd, Md, delta, lam, d, a_prev, Lc):
    """Largest integer u >= 2 that passes the frozen assembly's label-balance test (exponent_of, section
    'labels') at depth d: logM - lD >= slack with lD = 2 u log2(rd), logM = lam + log2(Md),
    lnD = max(lD ln2, 1); slack = log2(8192 Lc^2) + 2 log2(1/delta) + log2(ln(172 Lc) + lnD) at d = 0 and
    log2(8) + a_prev + 1 + 2 log2(1/delta) + log2(lnD) at d >= 1."""
    lr = math.log2(rd)
    two_lr = 2.0 * lr
    logM = lam + math.log2(Md)
    inv = 2.0 * math.log2(1.0 / delta)
    if d == 0:
        base = math.log2(8192.0 * Lc * Lc) + inv
        c1 = math.log(172.0 * Lc)
    else:
        base = math.log2(8.0) + a_prev + 1.0 + inv
        c1 = None

    def slack(u):
        lnD = max(2.0 * u * lr * LN2, 1.0)
        return base + (math.log2(c1 + lnD) if d == 0 else math.log2(lnD))

    def passes(u):
        lD = 2.0 * u * lr
        lnD = max(lD * LN2, 1.0)
        return logM - lD >= slack(u)

    u = max(2, int(logM / two_lr))
    for _ in range(6):
        nu = int((logM - slack(u)) / two_lr)
        if nu >= u:
            break
        if nu < 2:
            return 2
        u = nu
    while u > 2 and not passes(u):
        u -= 1
    return max(2, u)


def est_u(rd, Md, delta, lam, lp_in, d=0, a_prev=0, Lc=1 << 20):
    """Closed-form greedy largest u_d: the smaller of the padding budget and the label-balance test."""
    return min(est_u_pad(rd, Md, delta, lam, lp_in), bal_u(rd, Md, delta, lam, d, a_prev, Lc))


def rank_e(r, delta, Lc, lam, leafvars, cK=2.0, shellmm=False):
    """Analytic replica: greedy u_d from the padding budget, a_d = min_a, cost exponent of the assembly."""
    L = len(r)
    M = _M(r, delta, Lc, True)
    dm = 2 * delta
    u, lp_in, cum, worst, a_prev = [], 0.0, 0.0, 0.0, 0
    for d in range(L):
        ud = est_u(r[d], M[d], delta, lam, lp_in, d, a_prev, Lc)
        ad = min_a(ud, lp_in, M[d])
        a_prev = ad
        u.append(ud)
        cKd = cK * (A.gamma_guard(ad, _KAPPA) if _KAPPA is not None else 1.0)
        x = M[d] / math.sqrt(ud)
        ar = A.tup(ud, A.eps_mode(r[d], {"modedp": True, "shellmm": shellmm}), dm, cKd) * x
        ch = A.tup(ud, A.f_grp(r[d]), dm, cKd) * x
        worst = max(worst, cum + ar)
        cum += ch
        lp_in = A.next_lp(lp_in, M[d], ad, ud, _KAPPA)
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
        lp_in = A.next_lp(lp_in, M[d], a[d], u[d], _KAPPA)
    p = dict(levels=L, r=list(r), u=list(u), a=a, delta=delta, Lc=Lc, **flags)
    if _KAPPA is not None:
        p["kappa"] = [_KAPPA] * L
    return p


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
    a_prev = 0
    for d in range(L):
        u[d] = min(u[d], est_u(r[d], M[d], delta, lam, lp, d, a_prev, Lc))
        ad = min_a(u[d], lp, M[d])
        a_prev = ad
        lp = A.next_lp(lp, M[d], ad, u[d], _KAPPA)
    for d in range(L):
        m = max_up(u, d, UPMAX)
        if m is None:
            return None
        u[d] = m

    def full_ok(u0):
        uu = list(u)
        uu[0] = u0
        a_prev = min_a(u0, 0.0, M[0])
        lp = A.next_lp(0.0, M[0], a_prev, u0, _KAPPA)
        for d in range(1, L):
            uu[d] = min(max(2, est_u(r[d], M[d], delta, lam, lp, d, a_prev, Lc)), UPMAX)
            ad = min_a(uu[d], lp, M[d])
            a_prev = ad
            lp = A.next_lp(lp, M[d], ad, uu[d], _KAPPA)
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
    global _KAPPA
    flags = dict(perclass=True, labelfree=True, memo=False, modedp=False, chunks=False, leafvars=False,
                 newton=False, ldu=False, shellmm=False, guarded=False)
    if USE_CODE_FLAGS:
        flags.update(memo=True, modedp=True, chunks=True, newton=True, ldu=True, shellmm=True, guarded=True)
    _KAPPA = GUARD_KAPPA if flags["guarded"] else None
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
                e, _ = rank_e(list(shape), delta, Lc, lam, lv, cK, flags["shellmm"])
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
    return make([2], [2], 1e-4, 4096, flags)
