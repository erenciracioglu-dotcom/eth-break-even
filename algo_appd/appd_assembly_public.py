"""FROZEN Appendix D assembly, referee_appd version 1 (2026-11-01).

Cost exponent e (cost ~ 2^(e n), lam = log2 n, polynomial factors ignored) of the tensor-contraction algorithm
of Appendix D (Theorem D.1) on 3-CNF with at most three occurrences per variable, with the head-audited steps of
state/claims.md: appd-tight, appd-bernstein, appd-memo, appd-modedp, appd-chunks, appd-perclass, appd-leafvars,
appd-labelfree. All wire counts are per n.

params = {"levels": L in 1..3, "r": [r_0..r_{L-1}] groups per level, "u": [u_d] retained valuation levels
(integers >= 2), "a": [a_d] packet coordinates (L_d = 2^a_d, integers), "delta": balance tolerance,
"Lc": chunk size, and boolean flags memo, modedp, chunks, perclass, leafvars, labelfree}.
Flags memo, modedp, chunks and leafvars count only when the verify gate confirmed the agent's code for them
(see referee.py); perclass and labelfree are accounting facts (audited) and need no code.
"""
import math

LN2 = math.log(2.0)
TARGETS = {"third": 1.0 / 3.0, "occ3sat": math.log2(1.1199), "ppsz": math.log2(1.307031578),
           "sharp3sat": math.log2(1.637), "brute": 1.0}
PRIMARY = "third"            # fair counting competitor for this class: 2^(n/3) (pathwidth DP)
FLAGS = ("memo", "modedp", "chunks", "perclass", "leafvars", "labelfree")
CODE_FLAGS = ("memo", "modedp", "chunks", "leafvars")


def eps_dp(r):
    return max(2 * (k - 1) * (r - k) + 2 * (r - 1) for k in range(1, r + 1)) / float(r * r)


def f_all(r):
    return (r - 1.0) / r


def f_grp(r):
    return 2.0 * (r - 1) / (r * r)


def tup(u, f, dm, cK=2.0):
    """log2(#mode tuples) per unit M/sqrt(u): K = cK (1 + dm) f M / D per weight, p_t = 6 h D, h = M/(u - 1)."""
    return (2 * math.pi / LN2) * math.sqrt((u / (u - 1.0)) * cK * (1 + dm) * f)


def lse2(*xs):
    m = max(xs)
    return m + math.log2(sum(2.0 ** (x - m) for x in xs))


def _bad(why):
    return {"valid": False, "e": None, "why": why}


def exponent_of(params, lam):
    try:
        L = int(params["levels"])
        r = [int(x) for x in params["r"]]
        u = [int(x) for x in params["u"]]
        a = [int(x) for x in params["a"]]
        delta = float(params["delta"])
        Lc = int(params.get("Lc", 1))
        fl = {k: bool(params.get(k, False)) for k in FLAGS}
        lam = float(lam)
    except Exception as exc:  # noqa: BLE001
        return _bad(f"malformed params: {exc}")
    if not 1 <= L <= 3 or not (len(r) == len(u) == len(a) == L):
        return _bad("levels must be 1..3 with one r, u, a per level")
    if lam < 64:
        return _bad("lam below 64")
    if not 1e-4 <= delta <= 0.25:
        return _bad("delta outside [1e-4, 0.25]")
    if not 1 <= Lc <= 2 ** 20:
        return _bad("Lc outside [1, 2^20]")
    if fl["modedp"] and not fl["memo"]:
        return _bad("modedp needs memo (the elimination reads the child tables)")
    if fl["leafvars"] and not (fl["chunks"] and L == 1):
        return _bad("leafvars applies to the depth-1 leaves of a chunk network only")
    for d in range(L):
        if not 2 <= r[d] <= 4096:
            return _bad("r outside [2, 4096]")
        if not 2 <= u[d] < a[d] <= 10 ** 7:
            return _bad("need 2 <= u < a <= 1e7")
    M = [(1.0 + 1.0 / Lc) if fl["chunks"] else 3.0]
    for d in range(L):
        M.append(M[-1] * 2 * (1 + delta) / r[d])
    checks = {}
    # precision headroom (appd-tight item 2): (p_in + 4 M (a + 1)) / 2^(a - u + 1) <= 4 h, h = M/(u - 1)
    lp_in = 0.0                                   # root: P0 = n + 1, i.e. log2(p_in / n) = 0
    for d in range(L):
        need = lse2(lp_in, math.log2(4 * M[d] * (a[d] + 1))) + math.log2(u[d] - 1) - math.log2(4 * M[d])
        if need > a[d] - u[d] + 1:
            return _bad(f"precision headroom at depth {d}: need a - u + 1 >= {need:.2f}")
        h = M[d] / (u[d] - 1)
        # this call's worst-node precision (i = 0) is what its children inherit
        lp_in = lse2(math.log2(6 * h) + a[d], lp_in, math.log2(4 * M[d] * (a[d] + 1)))
    # labels (appd-bernstein, appd-perclass, appd-labelfree)
    for d in range(L):
        lab = u[d] if fl["labelfree"] else a[d]
        lr = math.log2(r[d])
        lD = 2 * lab * lr
        logM = lam + math.log2(M[d])
        lnD = max(lD * LN2, 1.0)
        if d == 0 and fl["chunks"]:
            slack = math.log2(8192.0 * Lc * Lc) + 2 * math.log2(1 / delta) + math.log2(math.log(172.0 * Lc) + lnD)
        else:
            ldeg = math.log2(3.0) if d == 0 else a[d - 1] + 1
            slack = math.log2(8.0) + ldeg + 2 * math.log2(1 / delta) + math.log2(lnD)
        if logM - lD < slack:
            return _bad(f"label balance at depth {d}: log2(M/D_arr) = {logM - lD:.1f} < {slack:.1f}")
        pad = a[d] + lab * (math.log2(r[d] * (r[d] + 1) / 2.0) if fl["perclass"] else 2 * lr)
        if pad > logM + math.log2(delta):
            return _bad(f"padding at depth {d}: {pad:.1f} > {logM + math.log2(delta):.1f}")
        checks[f"d{d}"] = dict(lD=round(lD, 1), pad_slack=round(logM + math.log2(delta) - pad, 2),
                               bal_slack=round(logM - lD - slack, 2))
    dm = 2 * delta if fl["perclass"] else delta
    worst, cum, terms = 0.0, 0.0, []
    for d in range(L):
        x = M[d] / math.sqrt(u[d])
        arith = tup(u[d], eps_dp(r[d]) if fl["modedp"] else f_all(r[d]), dm) * x
        if fl["memo"]:
            child = tup(u[d], f_grp(r[d]), dm) * x
            worst = max(worst, cum + arith)
            terms.append(dict(arith=cum + arith, child=child))
            cum += child
        else:
            terms.append(dict(arith=arith))
            cum += arith
    leaf = M[L]
    if fl["leafvars"]:
        leaf = min(leaf, (1 + delta) * (1 + M[0] / 2) / r[0])
    e = max(worst, cum + leaf)
    out = {"valid": True, "e": e, "terms": terms, "leaf": cum + leaf, "M_root": M[0], "checks": checks}
    out.update({k: fl[k] for k in FLAGS})
    return out
