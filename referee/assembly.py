"""FROZEN cost assembly (the referee scores with this file). Version 6 (2026-10-31).

Bounded-occurrence 3-CNF (each variable in <= 3 clauses, clauses of width <= 3), counting
algorithm of Sections 2-4 of "The Exponential Time Hypothesis is False". Exponents per
variable: cost ~ 2^(e*n), polynomial factors ignored. lam = log2 n, M = 3n occurrences.

algo/assembly_public.py is an identical copy for the agent's own parameter search.
Editing that copy does not change any score.

Parameters: r, Q, R, A0 (integers), regime (string), delta (float, regime "tight" only),
designated (bool, default False), aggregated (bool, default False), units (bool, default False),
chunk (int, default 0), modular (bool, default False), perpoint (bool, default False),
eadapt (bool, default False).

Regimes (how many labels D may be balanced and how many copies h that costs):
  paper    : Lemma 4.5 as written, 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M)
  freedman : head extension, audited as a sketch (state/claims.md): D ln D <= M / (4(b+w)),
             h <= 2M/D (two-stage Bernstein + union bound + pessimistic estimators)
  tight    : same argument with a tolerance 0 < delta <= 1: D ln D <= delta^2 M / (4(b+w)),
             h <= (1 + delta) M/D. delta = 1 is the freedman regime.

designated = True (head extension v2, audited as a sketch): every clause reads one of its own
variables v(f) directly and is assigned to that variable's group at every point
(sigma_f = a_v(f)). That occurrence needs no replacement row (degree 0) and no auxiliary
position, so only M_eff = sum_f (|I_f| - 1) <= 2n occurrences carry replacements: the
T-degree of H drops to M_eff (R-1) and hD <= (1 + delta) M_eff. Balancing uses one stage per
colour of a proper colouring of the variable conflict graph (<= 7 colours), which costs the
constant in D ln(28 D) <= delta^2 M_eff / 128. Allowed with regimes freedman and tight only.

aggregated = True (head extension v3, audited as a sketch): Lemma 3.3 summed per U,
  sum_x p_x(x) = sum_{|U| <= q} sum_{x >= U} w_K(|x minus U|) p_x(1_U),  K = q - |U|,
  w_K(m) = sum_{j <= K} (-1)^j C(m, j).
The weight depends on x only through an additive count, so each U needs one group-separated
sum (with a polynomial in y for the count). Restricted sums per point: sum_{s <= q} C(n, s)
instead of sum_s C(n, s) 2^s, so the Taylor exponent is H2(x) (1 for x >= 1/2) instead of
x + H2(x) (log2 3 for x >= 2/3).

units = True (head extension v4, needs designated = True; audited as a sketch): a greedy maximal
matching of clauses; the variables of every matched clause share one schedule (a unit), so all
occurrences of a matched clause are read directly, and so are the occurrences of other clauses
whose variables lie in the clause's designated unit. Each matched width-3 clause blocks at most
6 others, so at least |F|/7 clauses are matched and M_eff <= 2n - 2n/7 = 12n/7. Balancing by
colour classes of units (<= 19 colours, <= 18 edges per unit): D ln(76 D) <= delta^2 M_eff / 1024.

chunk = L (integer 2 <= L <= 4096, needs units = True; head extension v5, audited as a sketch):
units are connected chunks of the variable graph (variables sharing a clause) instead of matched
clauses. A component with fewer than L variables is one unit; a larger component is cut along a
BFS spanning tree into connected units of L to 1 + 7(L - 1) variables. A unit of s variables that
is connected through its clauses saves at least s - 1 replacement rows, and M - |F| <= 2n, so
M_eff <= n + (number of units) <= (1 + 1/L) n. Balancing by colour classes of units (<= 42L + 1
colours, unit weight <= 63L): D ln(172 L D) <= delta^2 M_eff / (8192 L^2).

modular = True (head extension v5, audited as a sketch; decision targets only, needs designated
and aggregated): decide satisfiability through S(w) = sum over satisfying x of prod_v w_v^(x_v)
for random weights in GF(p^k), with a prime Q < p <= 2Q (Bertrand) and p^k >= 2n
(Schwartz-Zippel, one-sided error, like PPSZ). Every coefficient is p-integral (denominators
divide 2^a (R-1)! Q!), so the measurement values are needed only mod p: log2(2Q) bits per value
instead of Q. Only exponent_of(..., decision=True) accepts it.

perpoint = True (head extension v6, audited as a sketch): Definition 3.1 lets the measurement
matrices C_{t,g} depend on the point t. At point t every operator X_t (x) Y whose position-t factor
maps (a, d) only to pairs (d, *) is free (it stays in the clause's group at t); the swap of position
t and the projector on equal pairs at t are such operators. Writing I - M0 by the type of its
position-t factor (equal, symmetric, antisymmetric) and dropping the free parts leaves
P^N_t (x) (m(j) + m(j + 1)) with m(j) = 1 - (-1)^j f(j), j the antisymmetric count of the other
Q - 1 positions. So the measured fraction at every point is
   out = p * Pr[Bin(Q - 1, p) outside [A0, A0 + R - 2]],  p = (r - 1)/(2r),
instead of Pr[Bin(Q, p) outside [A0, A0 + R - 1]]. Rows still touch 2 groups, entries {0, +-1}.

eadapt = True (head extension v6, audited as a sketch): the interpolation window depends on the
number e of equal pairs (a(s) = d(s)) of a label: M0 = V0 sum_e f_e(W0) [N_E = e], where N_E
counts equal pairs (diagonal, commutes with every swap, so the support bound |E_uw| <= R - 1 is
unchanged). Given e the antisymmetric count is Bin(Q - e, 1/2), and I_e = [a_e, a_e + R - 1] with
a_e = (Q - e - R + 1) // 2 (no rows at all when Q - e < R). A0 is ignored. Measured fraction
   out = sum_e Pr[Bin(Q, 1/r) = e] Pr[Bin(Q - e, 1/2) outside I_e],
and with perpoint as well
   out = p * sum_e Pr[Bin(Q - 1, 1/r) = e] Pr[Bin(Q - 1 - e, 1/2) outside [a_e, a_e + R - 2]].
For large Q (more than 300 values of e within 14 sd) the sum over e is taken in about 65 blocks
between fixed quantiles, each scored at its smallest e (the largest spread), and the mass beyond 14 standard deviations counts as fully measured, so the
value is an upper bound.

A score that uses designated, aggregated, units, chunk, modular, perpoint or eadapt also needs
the verify gate to confirm exact counts with that construction or counting routine (see referee.py).
"""
import math
import numpy as np
from scipy.stats import binom

LOG2_3 = math.log2(3)
TARGETS = {"brute": 1.0, "sharp3sat": math.log2(1.637), "ppsz": math.log2(1.307031578)}
DECISION_TARGETS = ("ppsz",)  # yes/no targets: modular parameters are allowed only for these
REPORT_TARGETS = {"occ3sat": math.log2(1.1199)}  # reported, decision: best SAT bound, <= 3 occurrences
M_PER_N = 3        # occurrences per variable, worst case of the class
M_EFF_DESIG = 2    # occurrences per variable that still need replacement rows when designated
M_EFF_UNITS = 12 / 7  # ... with units (greedy matching of clauses), worst case of the class
BW = 6             # b + w (b <= 3 factors per variable, w <= 3 variables per factor)
DESIG_CONST = 128  # D ln(28 D) <= delta^2 M_eff / 128 for the designated construction
UNITS_CONST = 1024  # D ln(76 D) <= delta^2 M_eff / 1024 with units
CHUNK_CONST = 8192  # D ln(172 L D) <= delta^2 M_eff / (8192 L^2) with chunk units
CHUNK_OFFSET = 172
CHUNK_MAX = 4096
REGIMES = ("paper", "freedman", "tight")
DELTA_MIN = 1e-6


def H2(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def taylor(x):
    """(1/n) log2 sum_{j <= x n} 2^j C(n, j) for large n (pairs U <= E, |E| <= q)."""
    if x <= 0:
        return 0.0
    if x >= 2.0 / 3.0:
        return LOG2_3
    return x + H2(x)


def taylor_aggregated(x):
    """(1/n) log2 sum_{j <= x n} C(n, j) for large n (one restricted sum per U, |U| <= q)."""
    if x <= 0:
        return 0.0
    if x >= 0.5:
        return 1.0
    return H2(x)


def taylor_term(x, aggregated=False):
    return taylor_aggregated(x) if aggregated else taylor(x)


# ---- accounting facts (each one is also checked on real designs by the measure gate)
def groups_per_row(r):            # a measurement row is non-zero on at most 2 groups
    return min(2, r)


def row_abs_sum_bound(Q):         # restricted row: absolute coefficient sum <= 2^(Q-1)
    return 2 ** (Q - 1) if Q < 1000 else math.inf


def value_bits(Q, modular=False):  # values in [-2^(Q-1), 2^(Q-1)]: log2(2^Q + 1); mod p <= 2Q: log2(2Q)
    if modular:
        return math.log2(2.0 * Q)
    return float(Q) if Q > 60 else math.log2(2.0 ** Q + 1)


def row_degree(R):                # replacement rows have degree <= R - 1
    return R - 1


def occ_per_n(designated=False, units=False, chunk=0):
    """Occurrences per variable that carry replacement rows."""
    if units and chunk:
        return 1.0 + 1.0 / chunk
    if units:
        return M_EFF_UNITS
    return M_EFF_DESIG if designated else M_PER_N


def _solve(target, offset):
    """Solve x + log2(offset + x ln 2) = target for x = log2 D."""
    x = max(target, 1.0)
    for _ in range(200):
        x = target - math.log2(max(offset + x * math.log(2), 1e-12))
    return x


def max_log2D(lam, regime, delta=1.0, designated=False, units=False, chunk=0):
    """Largest allowed log2 D for n = 2^lam."""
    if regime not in REGIMES:
        raise ValueError(regime)
    d = delta if regime == "tight" else 1.0
    if units and chunk:
        if not designated or regime == "paper":
            raise ValueError("chunk units need designated = True and regime freedman or tight")
        target = lam + math.log2(1.0 + 1.0 / chunk) + 2 * math.log2(d) - math.log2(CHUNK_CONST * chunk * chunk)
        return _solve(target, math.log(CHUNK_OFFSET * chunk))  # D ln(172 L D) <= d^2 M_eff / (8192 L^2)
    if units:
        if not designated or regime == "paper":
            raise ValueError("units needs designated = True and regime freedman or tight")
        target = lam + math.log2(M_EFF_UNITS) + 2 * math.log2(d) - math.log2(UNITS_CONST)
        return _solve(target, math.log(76))     # D ln(76 D) <= d^2 M_eff / 1024
    if designated:
        if regime == "paper":
            raise ValueError("designated needs regime freedman or tight")
        target = lam + math.log2(M_EFF_DESIG) + 2 * math.log2(d) - math.log2(DESIG_CONST)
        return _solve(target, math.log(28))     # D ln(28 D) <= d^2 M_eff / 128
    logM = lam + math.log2(M_PER_N)
    if regime == "paper":                       # 2(b+w) D^2 <= M
        return (logM - math.log2(2 * BW)) / 2
    target = logM - math.log2(4 * BW) + 2 * math.log2(d)
    return _solve(target, 0.0)                  # D ln D <= d^2 M / (4(b+w))


def hD_over_n(lam, log2D, regime, delta=1.0, designated=False, units=False, chunk=0):
    """h * D / n with the regime's bound on the number of copies h."""
    m = occ_per_n(designated, units, chunk)
    if regime == "paper":                       # (M/D + sqrt(2(b+w)M)) D / n = 3 + 6 D / sqrt(n)
        ex = log2D - lam / 2.0
        return M_PER_N + math.sqrt(2 * BW * M_PER_N) * (2.0 ** ex if ex > -300 else 0.0)
    if regime == "freedman":                    # 2 M/D * D / n
        return 2.0 * m
    if regime == "tight":                       # (1 + delta) M/D * D / n
        return (1.0 + delta) * m
    raise ValueError(regime)


QUANTILES = tuple(sorted({1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02}
                         | {round(x, 4) for x in np.linspace(0.03, 0.97, 48)}
                         | {0.98, 0.99, 0.995, 0.999, 1 - 1e-4, 1 - 1e-6, 1 - 1e-9, 1 - 1e-12}))


def _eadapt_terms(Q, r, perpoint):
    """Blocks of the equal-pair count e over the positions other than the point's own (all Q of
    them without perpoint): block starts, block weights, and the mass beyond 14 sd."""
    n = Q - 1 if perpoint else Q
    if n <= 0:
        return np.zeros(1, np.int64), np.ones(1), 0.0
    q = 1.0 / r
    mu, sd = n * q, math.sqrt(n * q * (1.0 - q))
    lo = max(0, int(math.floor(mu - 14.0 * sd)))
    hi = min(n, int(math.ceil(mu + 14.0 * sd)))
    if hi - lo + 1 <= 300:
        es = np.arange(lo, hi + 1, dtype=np.int64)
        w = binom.pmf(es, n, q)
    else:                                       # blocks [edges[i], edges[i+1] - 1] at fixed quantiles
        qs = np.array(QUANTILES)
        mids = binom.ppf(qs, n, q).astype(np.int64)
        edges = np.unique(np.concatenate([[lo], mids, [hi + 1]]))
        edges = edges[(edges >= lo) & (edges <= hi + 1)]
        w = np.diff(binom.cdf(edges - 1, n, q))
        es = edges[:-1]
    rest = max(0.0, 1.0 - float(w.sum()))
    return es, w, rest


def _inner_eadapt(Q, es, Rs, perpoint):
    """Pr[antisymmetric count outside the window] for every block start e (rows) and R (columns)."""
    E = np.asarray(es, np.int64)[:, None]
    Rm = np.asarray(Rs, np.int64)[None, :]
    m = Q - E                                    # positions that are not equal pairs
    full = m < Rm                                # every count fits the window: nothing to measure
    a = np.where(full, 0, (m - Rm + 1) // 2)
    if perpoint:
        N = np.maximum(m - 1, 0)
        bad = binom.cdf(a - 1, N, 0.5) + binom.sf(a + Rm - 2, N, 0.5)
    else:
        N = np.maximum(m, 0)
        bad = binom.cdf(a - 1, N, 0.5) + binom.sf(a + Rm - 1, N, 0.5)
    return np.where(full, 0.0, bad)


def tail_out(Q, r, R, A0, perpoint=False, eadapt=False):
    """Fraction of basis vectors that must be measured (at every point, with perpoint):
    base Pr[Bin(Q, p) outside [A0, A0+R-1]]; see the module notes for perpoint and eadapt."""
    p = (r - 1) / (2 * r)
    if eadapt:
        es, w, rest = _eadapt_terms(Q, r, perpoint)
        val = float((w * _inner_eadapt(Q, es, np.array([R]), perpoint)[:, 0]).sum()) + rest
        return min(1.0, (p if perpoint else 1.0) * val)
    if perpoint:
        return float(p * (binom.cdf(A0 - 1, Q - 1, p) + binom.sf(A0 + R - 2, Q - 1, p)))
    return float(binom.cdf(A0 - 1, Q, p) + binom.sf(A0 + R - 1, Q, p))


def _flag(params, name):
    v = params.get(name, False)
    if v not in (True, False):
        raise ValueError(f"{name} must be true or false")
    return bool(v)


def _chunk(params):
    c = params.get("chunk", 0)
    if c is None or c is False:
        return 0
    if isinstance(c, bool) or float(c) != int(c):
        raise ValueError("chunk must be an integer")
    c = int(c)
    if c and not (2 <= c <= CHUNK_MAX):
        raise ValueError(f"chunk must be 0 or lie in [2, {CHUNK_MAX}]")
    return c


def _extra(params, regime):
    """delta, designated, aggregated, units, chunk, modular, perpoint and eadapt from params;
    ValueError when not allowed."""
    designated, aggregated, units = _flag(params, "designated"), _flag(params, "aggregated"), _flag(params, "units")
    chunk, modular = _chunk(params), _flag(params, "modular")
    perpoint, eadapt = _flag(params, "perpoint"), _flag(params, "eadapt")
    delta = 1.0
    if regime == "tight":
        delta = float(params["delta"])
        if not (DELTA_MIN <= delta <= 1.0):
            raise ValueError(f"delta must lie in [{DELTA_MIN}, 1]")
    if designated and regime == "paper":
        raise ValueError("designated needs regime freedman or tight")
    if units and not designated:
        raise ValueError("units needs designated = True")
    if chunk and not units:
        raise ValueError("chunk needs units = True")
    if modular and not (designated and aggregated):
        raise ValueError("modular needs designated = True and aggregated = True")
    return delta, designated, aggregated, units, chunk, modular, perpoint, eadapt


def exponent_of(params, lam, decision=False):
    """Score one parameter set at lam = log2 n. Returns dict(valid, e, terms, ...).
    decision=True scores a yes/no target (the only place where modular = True is allowed)."""
    try:
        r, Q, R, A0 = int(params["r"]), int(params["Q"]), int(params["R"]), int(params["A0"])
        regime = params.get("regime", "paper")
        if regime not in REGIMES:
            return dict(valid=False, reason=f"unknown regime {regime}")
        delta, designated, aggregated, units, chunk, modular, perpoint, eadapt = _extra(params, regime)
    except Exception as exc:                    # malformed parameters
        return dict(valid=False, reason=f"bad params: {exc}")
    if r < 2 or Q < 1 or R < 1 or A0 < 0 or A0 + R - 1 > Q:
        return dict(valid=False, reason="r>=2, Q>=1, R>=1 and 0 <= A0 <= A0+R-1 <= Q required")
    if modular and not decision:
        return dict(valid=False, reason="modular = True is only valid for decision targets")
    log2D = 2 * Q * math.log2(r)
    if log2D > max_log2D(lam, regime, delta, designated, units, chunk) + 1e-9:
        return dict(valid=False, reason="label space D too large for the balancing regime")
    out = tail_out(Q, r, R, A0, perpoint, eadapt)
    rows_per_n = groups_per_row(r) * hD_over_n(lam, log2D, regime, delta, designated, units, chunk) * out
    terms = dict(taylor=taylor_term(occ_per_n(designated, units, chunk) * row_degree(R) / Q, aggregated),
                 meas=rows_per_n * value_bits(Q, modular),
                 groups=1.0 / r)
    return dict(valid=True, e=sum(terms.values()), terms=terms, out=out, log2D=log2D, regime=regime,
                delta=delta, designated=designated, aggregated=aggregated, units=units, chunk=chunk,
                modular=modular, perpoint=perpoint, eadapt=eadapt, decision=bool(decision))


def exponents_over_R(lam, r, Q, Rs, regime, delta=1.0, designated=False, aggregated=False, units=False,
                     chunk=0, modular=False, perpoint=False, eadapt=False):
    """Vectorised helper for searches: e for each interval length in Rs (best A0 each; A0 = 0 with
    eadapt, where the windows are fixed by the equal-pair count)."""
    Rs = np.asarray(Rs, np.int64)
    p = (r - 1) / (2 * r)
    best = np.ones(Rs.shape)
    bestA = np.zeros(Rs.shape, np.int64)
    if eadapt:
        es, w, rest = _eadapt_terms(Q, r, perpoint)
        best = (w[:, None] * _inner_eadapt(Q, es, Rs, perpoint)).sum(axis=0) + rest
        best = np.minimum(1.0, best * (p if perpoint else 1.0))
    else:
        nq, wl = (Q - 1, Rs - 1) if perpoint else (Q, Rs)
        mu = nq * p
        for shift in (-1, 0, 1):
            A0 = np.floor(mu - (wl - 1) / 2).astype(np.int64) + shift
            A0 = np.clip(A0, 0, np.maximum(Q - Rs + 1, 0))
            out = binom.cdf(A0 - 1, nq, p) + binom.sf(A0 + wl - 1, nq, p)
            if perpoint:
                out = p * out
            better = out < best
            best = np.where(better, out, best)
            bestA = np.where(better, A0, bestA)
    log2D = 2 * Q * math.log2(r)
    rows = groups_per_row(r) * hD_over_n(lam, log2D, regime, delta, designated, units, chunk) * best
    m = occ_per_n(designated, units, chunk)
    tay = np.array([taylor_term(m * row_degree(int(R)) / Q, aggregated) for R in Rs])
    e = tay + rows * value_bits(Q, modular) + 1.0 / r
    return e, bestA
# VERSION6_MARK
