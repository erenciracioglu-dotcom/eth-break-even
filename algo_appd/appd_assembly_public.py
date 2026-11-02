"""FROZEN Appendix D assembly, referee_appd version 7 (2026-11-02, round 11: profile flag).

Cost exponent e (cost ~ 2^(e n), lam = log2 n, polynomial factors ignored) of the tensor-contraction algorithm
of Appendix D (Theorem D.1) on 3-CNF with at most three occurrences per variable, with the head-audited steps of
state/claims.md: appd-tight, appd-bernstein, appd-memo, appd-modedp, appd-chunks, appd-perclass, appd-leafvars,
appd-labelfree, appd-newton, appd-ldu, appd-shellmm, appd-guarded, appd-bincross. All wire counts are per n.

params = {"levels": L in 1..3, "r": [r_0..r_{L-1}] groups per level, "u": [u_d] retained valuation levels
(integers >= 2), "a": [a_d] packet coordinates (L_d = 2^a_d, integers), "delta": balance tolerance,
"Lc": chunk size, and boolean flags memo, modedp, chunks, perclass, leafvars, labelfree, newton, ldu,
shellmm, guarded, bincross, profile} plus, with guarded, "kappa": [kappa_d] (guard width H_d = 2^(a_d - kappa_d), 2 <= kappa_d <= 30).
Flags memo, modedp, chunks, leafvars, newton, ldu, shellmm, guarded, bincross and profile count only when the verify gate confirmed the agent's code
for them (see referee.py); perclass and labelfree are accounting facts (audited) and need no code.
newton (appd-newton): crossing entries expanded in binomial moments, source forms sum_j binom(j, b) X_(y0 + c j);
the weight of mode b grows by v2(b!) >= b - log2(b + 1), so a class has about W/2 modes of weight <= W instead of W:
the mode-tuple constant cK = 2 of the power basis becomes cK = 1 (ck() below).
ldu (appd-ldu): the crossing block of each packet half is factored as a Cauchy matrix K = L D U over Z_(2)[[z]];
mode m = (U X)_m has weight a + 1 - i + 2b + 2 v2(b!) with b = floor(m/2^i), twice the Newton growth: cK = 1/2.
shellmm (appd-shellmm): every step of the group-by-group mode sum is a set of block matrix products, one per
triple of exact weights (shells) of the index sets x (processed-later), y (processed-current), z (current-later);
with Strassen on square blocks the mode fraction eps_dp(r) becomes eps_mm(r) (f_mm below, omega = log2 7).
guarded (appd-guarded, GPT-6 Astra): packets are virtual grids of Q = 2^a slots with empty guards of H = Q/2^kappa
slots at both ends and the deformation W~_xy(T) = l_y(x + T) (Lagrange basis on the grid, degree Q - 1 instead of
2Q - 1); nodes in [-H, H). The LDU modes and weights are unchanged, the Hermite loss halves and the packets hold
N = Q - 2H wires: K p scales by gamma_guard(a, kappa); headroom and inherited precision by need_bits, next_lp.
bincross (appd-bincross, GPT-6 Astra, head-verified): one binary crossing in the mode sum of every level with
r >= 5 groups. The two group-pair edges x = (g1, g3) and y = (g2, g4) are equality wires; on padded shell codes
of a and b bits the pair is SWAP^(m) (x) I, m = min(a, b), and SWAP^(m) = sum_i A_i with A_i = [t^i](I - tJ)^(m),
J = v v^T, v = (0, 1, -1, 0). A_i factors integrally through R_i (row: pairs {g1, g2}, {g3, g4}) or C_i (column:
pairs {g1, g4}, {g2, g3}) channels with unimodular S, T; column for 2i <= m, row otherwise, so every piece has at
most 2^(c m) channels, c = log2(12)/2. The transforms are tensor products of 4 x 4 integer matrices (cost
O(m^2) per entry of the old pair table, no division). The new index theta (grade wx + wy) has at most
2^(max(a,b) + (c-1) min(a,b)) values, i.e. u_theta = 2 (1 + (c - 1)^2) = 3.256 mode units instead of 4. Costs in
mode units (2 per group pair): preparing a pair table 4r - 6 - 2 eta, absorbing it 4r - 8, then A, the r - 4
singletons one at a time and B last with the pair-pair edge 4 + u_theta (eps_bc1 below); the mode fraction of a
level is min(eps_mm(r), eps_bc1(r)). Needs shellmm.
profile (appd-profile, GPT-6 Astra, head-verified): the mode sum of a node with r groups is split by complete
shell profiles (one exact weight w_e per group pair e, sum w_e < p); inside a profile every edge has a finite
alphabet of N_e tuples and the network is contracted group by group, every step one shell matrix product, in the
order of the frozen menu PF_MENUS[r] with the smallest cost bound max_k 125 (X + Y + Z) - 24 min(X, Y, Z) (X, Y,
Z the sums of b_e = ceil(log2 N_e) over the edges A-B, A-C, B-C of step k; alpha' = 101/125 > log2(7) - 2).
With l_e <= beta sqrt(2 w_e) and Cauchy-Schwarz, min over the menu <= the mean of its envelopes gives the mode
fraction eps <= (2 / r^2) max over one form per order of ||mean of the forms||^2 = PF_EPS[r] (exact, 101/125);
rounding to whole bits adds at most C(r, 2) bits, profiles and steps are polynomially many. The mode fraction of
a level is then min(eps_mm, eps_bc1 with bincross, PF_EPS[r]) for 4 <= r <= 16. Needs shellmm.
"""
import math
from fractions import Fraction

LN2 = math.log(2.0)
TARGETS = {"third": 1.0 / 3.0, "occ3sat": math.log2(1.1199), "ppsz": math.log2(1.307031578),
           "sharp3sat": math.log2(1.637), "brute": 1.0}
PRIMARY = "third"            # fair counting competitor for this class: 2^(n/3) (pathwidth DP)
FLAGS = ("memo", "modedp", "chunks", "perclass", "leafvars", "labelfree", "newton", "ldu", "shellmm", "guarded",
         "bincross", "profile")
CODE_FLAGS = ("memo", "modedp", "chunks", "leafvars", "newton", "ldu", "shellmm", "guarded", "bincross", "profile")


def eps_dp(r):
    return max(2 * (k - 1) * (r - k) + 2 * (r - 1) for k in range(1, r + 1)) / float(r * r)


def f_all(r):
    return (r - 1.0) / r


def f_grp(r):
    return 2.0 * (r - 1) / (r * r)


OMEGA_MM = math.log2(7.0)          # Strassen's algorithm


def f_mm(a, b, c, r, om=OMEGA_MM):
    """Mode fraction of one shell-matrix-product merge: blocks A (a groups) and B (b groups), c groups left,
    x = A-C, y = A-B, z = B-C. Naive cost 2(ab + ac + bc)/r^2; with square-block fast matrix multiplication of
    exponent om the smallest of the three index sets enters with weight (om - 2)^2 (appd-shellmm)."""
    s = a * b + a * c + b * c
    return 2.0 * (s - (1.0 - (om - 2.0) ** 2) * min(a * b, a * c, b * c)) / float(r * r)


def eps_mm(r):
    """appd-shellmm: eps of the group-by-group order (one group per step) with Strassen; eps_mm(r) <= eps_dp(r)."""
    return max([2.0 * (r - 1) / float(r * r)] + [f_mm(k - 1, 1, r - k, r) for k in range(2, r + 1)])


ETA_MM = 1.0 - (OMEGA_MM - 2.0) ** 2           # Strassen discount of the smallest index set
C_BIN = 0.5 * math.log2(12.0)                 # appd-bincross: channels per paired bit, 2^c = sqrt(12)
U_THETA = 2.0 * (1.0 + (C_BIN - 1.0) ** 2)     # mode units of the crossing index theta (two plain edges: 4)


def eps_bc1(r):
    """appd-bincross: mode fraction with one binary crossing (r >= 5; eps_mm otherwise). Steps in mode units: the
    pair tables 4r - 6 - 2 eta (prepare) and 4r - 8 (absorb); then pair A takes the t = r - 4 singletons one at a
    time (step j: shared 2j + 2, A side h + 4(t - j) + 4(j - 1) + 2(j - 1)(t - j), singleton side 4 + 2(t - j),
    h = 4 + u_theta) and pair B comes last (h + 4t). Equal to the exact best merge tree for r <= 16 (head DP)."""
    if r < 5:
        return eps_mm(r)
    t, h = r - 4, 4.0 + U_THETA
    qs = [4.0 * r - 6.0 - 2.0 * ETA_MM, 4.0 * r - 8.0, h + 4.0 * t]
    for j in range(1, t + 1):
        X = 2.0 * j + 2.0
        Y = h + 4.0 * (t - j) + 4.0 * (j - 1) + 2.0 * (j - 1) * (t - j)
        Z = 4.0 + 2.0 * (t - j)
        qs.append(X + Y + Z - ETA_MM * min(X, Y, Z))
    return min(eps_mm(r), max(qs) / float(r * r))


PF_ALPHA_NUM, PF_ALPHA_DEN = 101, 125          # appd-profile: alpha' = 101/125 > log2(7) - 2 (7^125 < 2^351)
# appd-profile menus: orders of the positions 0..r-1 of a node's groups (first the identity) and their certificates
PF_MENUS = {
    4: ((0, 1, 2, 3), (0, 2, 3, 1), (1, 2, 0, 3)),
    5: ((0, 1, 2, 3, 4), (0, 3, 4, 2, 1), (2, 3, 0, 4, 1), (3, 1, 0, 4, 2)),
    6: ((0, 1, 2, 3, 4, 5), (2, 3, 0, 5, 4, 1), (4, 2, 0, 5, 1, 3), (5, 2, 1, 4, 3, 0)),
    7: ((0, 1, 2, 3, 4, 5, 6), (1, 0, 4, 5, 3, 6, 2), (0, 3, 5, 4, 2, 1, 6), (5, 2, 0, 6, 1, 4, 3)),
    8: ((0, 1, 2, 3, 4, 5, 6, 7), (7, 0, 6, 2, 1, 4, 3, 5), (3, 4, 0, 6, 2, 1, 7, 5), (2, 7, 3, 4, 6, 5, 1, 0)),
    9: ((0, 1, 2, 3, 4, 5, 6, 7, 8), (3, 7, 4, 0, 6, 5, 1, 8, 2), (5, 3, 4, 8, 1, 2, 7, 0, 6), (7, 2, 4, 5, 0, 1, 6, 3, 8)),
    10: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), (3, 7, 2, 5, 6, 1, 9, 4, 0, 8), (5, 8, 1, 0, 7, 2, 6, 9, 4, 3), (6, 0, 9, 5, 3, 8, 7, 2, 1, 4)),
    11: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), (8, 9, 1, 5, 2, 6, 10, 7, 0, 3, 4), (9, 2, 4, 10, 5, 3, 8, 7, 1, 0, 6)),
    12: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), (7, 8, 0, 6, 3, 4, 9, 2, 1, 11, 10, 5), (4, 6, 11, 1, 7, 2, 8, 0, 3, 9, 5, 10)),
    13: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), (8, 6, 5, 3, 2, 10, 11, 0, 12, 1, 7, 9, 4), (11, 0, 5, 7, 2, 10, 1, 8, 6, 9, 12, 4, 3)),
    14: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13), (10, 4, 6, 13, 9, 2, 3, 8, 1, 5, 11, 7, 12, 0), (6, 8, 0, 10, 13, 1, 4, 12, 5, 7, 3, 9, 2, 11)),
    15: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), (13, 10, 8, 1, 4, 11, 0, 2, 14, 12, 9, 6, 3, 5, 7), (6, 5, 11, 2, 14, 8, 0, 10, 12, 7, 4, 1, 9, 13, 3)),
    16: ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), (1, 0, 8, 7, 2, 14, 12, 10, 3, 13, 15, 4, 6, 9, 11, 5), (3, 10, 2, 9, 0, 14, 5, 13, 4, 15, 6, 1, 7, 8, 11, 12)),
}
PF_EPS = {
    4: Fraction(185701, 375000),
    5: Fraction(185239, 390625),
    6: Fraction(332719, 750000),
    7: Fraction(651717, 1531250),
    8: Fraction(102163, 250000),
    9: Fraction(338143, 843750),
    10: Fraction(2440461, 6250000),
    11: Fraction(2280778, 5671875),
    12: Fraction(265501, 675000),
    13: Fraction(620336, 1584375),
    14: Fraction(1781599, 4593750),
    15: Fraction(1349194, 3515625),
    16: Fraction(6840623, 18000000),
}


def pf_forms(order, r):
    """The cost forms of a sequential order, scaled by PF_ALPHA_DEN: per step k and per cross class, coefficient
    PF_ALPHA_NUM on that class, PF_ALPHA_DEN on the other two, 0 on the edges inside A, B or C (edge order: (i, j),
    i < j, lexicographic)."""
    E = [(i, j) for i in range(r) for j in range(i + 1, r)]
    pos = {g: k for k, g in enumerate(order)}
    out = set()
    for k in range(1, r):
        cls = lambda g: 0 if pos[g] < k else (1 if pos[g] == k else 2)
        lab = [tuple(sorted((cls(i), cls(j)))) for i, j in E]
        for weak in ((0, 1), (0, 2), (1, 2)):
            out.add(tuple(0 if a == b else (PF_ALPHA_NUM if (a, b) == weak else PF_ALPHA_DEN) for a, b in lab))
    return sorted(out)


def pf_certificate(r):
    """Exact recomputation of PF_EPS[r] from PF_MENUS[r] (slow for r >= 14; used by the head's validation)."""
    menu = PF_MENUS[r]
    k = len(menu)
    fs = [pf_forms(o, r) for o in menu]
    best = 0
    acc = [tuple(f) for f in fs[0]]
    for j in range(1, k):
        acc = [tuple(x + y for x, y in zip(a, b)) for a in acc for b in fs[j]]
    best = max(sum(x * x for x in a) for a in acc)
    return Fraction(2 * best, k * k * PF_ALPHA_DEN * PF_ALPHA_DEN * r * r)


def eps_mode(r, fl):
    """Mode fraction of the root sum: all modes without modedp, eps_dp with it, eps_mm with shellmm as well,
    min(eps_mm, eps_bc1) with bincross and at most PF_EPS[r] with profile (4 <= r <= 16)."""
    if not fl.get("modedp"):
        return f_all(r)
    if not fl.get("shellmm"):
        return eps_dp(r)
    e = eps_bc1(r) if fl.get("bincross") else eps_mm(r)
    if fl.get("profile") and r in PF_EPS:
        e = min(e, float(PF_EPS[r]))
    return e


def gamma_guard(a, kappa):
    """appd-guarded: ratio of K p (modes times precision) of the guarded Lagrange deformation against the paper's
    Cauchy deformation at packet size Q = 2^a and guards H = Q/2^kappa: (Q/N)(Q - 1)/(2Q - 1), N = Q - 2H."""
    return (1.0 / (1.0 - 2.0 ** (1 - kappa))) * (1.0 - 2.0 ** -a) / (2.0 - 2.0 ** -a)


def need_bits(lp_in, M, a, u, kappa=None):
    """log2 of the headroom 2^(a - u + 1) a level needs. Paper (appd-tight item 2): (p_in + 4M(a + 1))(u - 1)/(4M).
    Guarded: h' ~ 2^(kappa-1) M/(u - 1) and D'_t = 2^(a - kappa - i), so the slack 4 h' D' at the top level is
    2M 2^(a-u+1)/(u - 1) and the evaluation loss is below 2^(kappa+1) M (a + 1): (p_in + 2^(kappa+1) M (a+1))(u-1)/(2M)."""
    if kappa is None:
        return lse2(lp_in, math.log2(4 * M * (a + 1))) + math.log2(u - 1) - math.log2(4 * M)
    return lse2(lp_in, math.log2(2.0 ** (kappa + 1) * M * (a + 1))) + math.log2(u - 1) - math.log2(2 * M)


def next_lp(lp_in, M, a, u, kappa=None):
    """log2(precision of the worst node (i = 0) / n), which the children of this level inherit: 6 h 2^a for the
    paper's deformation, 6 h' H = 3 h 2^a for the guarded one (h = M/(u - 1)), plus p_in and the evaluation loss."""
    h = M / (u - 1)
    if kappa is None:
        return lse2(math.log2(6 * h) + a, lp_in, math.log2(4 * M * (a + 1)))
    return lse2(math.log2(3 * h) + a, lp_in, math.log2(2.0 ** (kappa + 1) * M * (a + 1)))


def tup(u, f, dm, cK=2.0):
    """log2(#mode tuples) per unit M/sqrt(u): K = cK (1 + dm) f M / D per weight, p_t = 6 h D, h = M/(u - 1)."""
    return (2 * math.pi / LN2) * math.sqrt((u / (u - 1.0)) * cK * (1 + dm) * f)


def ck(params):
    """Mode-tuple constant: modes per unit weight per crossing class times 2 (power 2, Newton 1, Cauchy LDU 1/2)."""
    if params.get("ldu"):
        return 0.5
    return 1.0 if params.get("newton") else 2.0


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
    kap = [None] * L
    if fl["guarded"]:
        try:
            kap = [int(x) for x in params["kappa"]]
        except Exception:  # noqa: BLE001
            return _bad("guarded needs an integer kappa per level")
        if len(kap) != L or not all(2 <= k <= 30 for k in kap):
            return _bad("kappa must have one value in [2, 30] per level")
        if not fl["ldu"]:
            return _bad("guarded needs ldu (its modes are the Cauchy LDU modes)")
    if lam < 64:
        return _bad("lam below 64")
    if not 1e-4 <= delta <= 0.25:
        return _bad("delta outside [1e-4, 0.25]")
    if not 1 <= Lc <= 2 ** 20:
        return _bad("Lc outside [1, 2^20]")
    if fl["modedp"] and not fl["memo"]:
        return _bad("modedp needs memo (the elimination reads the child tables)")
    if fl["shellmm"] and not fl["modedp"]:
        return _bad("shellmm needs modedp (it reorganises the group-by-group sum)")
    if fl["bincross"] and not fl["shellmm"]:
        return _bad("bincross needs shellmm (its merges are shell matrix products)")
    if fl["profile"] and not fl["shellmm"]:
        return _bad("profile needs shellmm (its steps are shell matrix products)")
    if fl["leafvars"] and not (fl["chunks"] and L == 1):
        return _bad("leafvars applies to the depth-1 leaves of a chunk network only")
    for d in range(L):
        if not 2 <= r[d] <= 4096:
            return _bad("r outside [2, 4096]")
        if not 2 <= u[d] < a[d] <= 10 ** 7:
            return _bad("need 2 <= u < a <= 1e7")
        if kap[d] is not None and u[d] > a[d] - kap[d]:
            return _bad("guards need u <= a - kappa")
    M = [(1.0 + 1.0 / Lc) if fl["chunks"] else 3.0]
    for d in range(L):
        M.append(M[-1] * 2 * (1 + delta) / r[d])
    checks = {}
    # precision headroom (appd-tight item 2): (p_in + 4 M (a + 1)) / 2^(a - u + 1) <= 4 h, h = M/(u - 1)
    lp_in = 0.0                                   # root: P0 = n + 1, i.e. log2(p_in / n) = 0
    for d in range(L):
        need = need_bits(lp_in, M[d], a[d], u[d], kap[d])
        if need > a[d] - u[d] + 1:
            return _bad(f"precision headroom at depth {d}: need a - u + 1 >= {need:.2f}")
        # this call's worst-node precision (i = 0) is what its children inherit
        lp_in = next_lp(lp_in, M[d], a[d], u[d], kap[d])
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
    cK = ck(fl)
    worst, cum, terms = 0.0, 0.0, []
    for d in range(L):
        x = M[d] / math.sqrt(u[d])
        cKd = cK * (gamma_guard(a[d], kap[d]) if kap[d] else 1.0)
        arith = tup(u[d], eps_mode(r[d], fl), dm, cKd) * x
        if fl["memo"]:
            child = tup(u[d], f_grp(r[d]), dm, cKd) * x
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
