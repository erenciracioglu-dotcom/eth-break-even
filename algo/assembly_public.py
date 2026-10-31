"""FROZEN cost assembly (the referee scores with this file). Version 2 (2026-10-31).

Bounded-occurrence 3-CNF (each variable in <= 3 clauses, clauses of width <= 3), counting
algorithm of Sections 2-4 of "The Exponential Time Hypothesis is False". Exponents per
variable: cost ~ 2^(e*n), polynomial factors ignored. lam = log2 n, M = 3n occurrences.

algo/assembly_public.py is an identical copy for the agent's own parameter search.
Editing that copy does not change any score.

Parameters: r, Q, R, A0 (integers), regime (string), delta (float, regime "tight" only),
designated (bool, optional, default False).

Regimes (how many labels D may be balanced and how many copies h that costs):
  paper    : Lemma 4.5 as written, 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M)
  freedman : head extension, audited as a sketch (state/claims.md): D ln D <= M / (4(b+w)),
             h <= 2M/D (two-stage Bernstein + union bound + pessimistic estimators)
  tight    : same argument with a tolerance 0 < delta <= 1: D ln D <= delta^2 M / (4(b+w)),
             h <= (1 + delta) M/D. delta = 1 is the freedman regime.

designated = True (head extension v2, PROVISIONAL until audited): every clause reads one of its
own variables v(f) directly and is assigned to that variable's group at every point
(sigma_f = a_v(f)). That occurrence needs no replacement row (degree 0) and no auxiliary
position, so only M_eff = sum_f (|I_f| - 1) <= 2n occurrences carry replacements: the
T-degree of H drops to M_eff (R-1) and hD <= (1 + delta) M_eff. Balancing uses one stage per
colour of a proper colouring of the variable conflict graph (<= 7 colours), which costs the
constant in D ln(28 D) <= delta^2 M_eff / 128. Allowed with regimes freedman and tight only.
A score that uses designated = True also needs the verify gate to confirm exact counts with
the designated construction (see referee.py).
"""
import math
import numpy as np
from scipy.stats import binom

LOG2_3 = math.log2(3)
TARGETS = {"brute": 1.0, "sharp3sat": math.log2(1.637), "ppsz": math.log2(1.307031578)}
M_PER_N = 3       # occurrences per variable, worst case of the class
M_EFF_DESIG = 2   # occurrences per variable that still need replacement rows when designated
BW = 6            # b + w (b <= 3 factors per variable, w <= 3 variables per factor)
DESIG_CONST = 128  # D ln(28 D) <= delta^2 M_eff / 128 for the designated construction
REGIMES = ("paper", "freedman", "tight")
DELTA_MIN = 1e-6


def H2(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def taylor(x):
    """(1/n) log2 sum_{j <= x n} 2^j C(n, j) for large n."""
    if x <= 0:
        return 0.0
    if x >= 2.0 / 3.0:
        return LOG2_3
    return x + H2(x)


# ---- accounting facts (each one is also checked on real designs by the measure gate)
def groups_per_row(r):            # a measurement row is non-zero on at most 2 groups
    return min(2, r)


def row_abs_sum_bound(Q):         # restricted row: absolute coefficient sum <= 2^(Q-1)
    return 2 ** (Q - 1) if Q < 1000 else math.inf


def value_bits(Q):                # values in [-2^(Q-1), 2^(Q-1)]: log2(2^Q + 1)
    return float(Q) if Q > 60 else math.log2(2.0 ** Q + 1)


def row_degree(R):                # replacement rows have degree <= R - 1
    return R - 1


def occ_per_n(designated=False):
    """Occurrences per variable that carry replacement rows."""
    return M_EFF_DESIG if designated else M_PER_N


def _solve(target, offset):
    """Solve x + log2(offset + x ln 2) = target for x = log2 D."""
    x = max(target, 1.0)
    for _ in range(200):
        x = target - math.log2(max(offset + x * math.log(2), 1e-12))
    return x


def max_log2D(lam, regime, delta=1.0, designated=False):
    """Largest allowed log2 D for n = 2^lam."""
    if regime not in REGIMES:
        raise ValueError(regime)
    d = delta if regime == "tight" else 1.0
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


def hD_over_n(lam, log2D, regime, delta=1.0, designated=False):
    """h * D / n with the regime's bound on the number of copies h."""
    m = occ_per_n(designated)
    if regime == "paper":                       # (M/D + sqrt(2(b+w)M)) D / n = 3 + 6 D / sqrt(n)
        ex = log2D - lam / 2.0
        return M_PER_N + math.sqrt(2 * BW * M_PER_N) * (2.0 ** ex if ex > -300 else 0.0)
    if regime == "freedman":                    # 2 M/D * D / n
        return 2.0 * m
    if regime == "tight":                       # (1 + delta) M/D * D / n
        return (1.0 + delta) * m
    raise ValueError(regime)


def tail_out(Q, r, R, A0):
    """Fraction of basis vectors that must be measured: Pr[W outside [A0, A0+R-1]], W ~ Bin(Q, p)."""
    p = (r - 1) / (2 * r)
    return float(binom.cdf(A0 - 1, Q, p) + binom.sf(A0 + R - 1, Q, p))


def _extra(params, regime):
    """delta and designated from params; raises ValueError when they are not allowed."""
    designated = params.get("designated", False)
    if designated not in (True, False):
        raise ValueError("designated must be true or false")
    delta = 1.0
    if regime == "tight":
        delta = float(params["delta"])
        if not (DELTA_MIN <= delta <= 1.0):
            raise ValueError(f"delta must lie in [{DELTA_MIN}, 1]")
    if designated and regime == "paper":
        raise ValueError("designated needs regime freedman or tight")
    return delta, bool(designated)


def exponent_of(params, lam):
    """Score one parameter set at lam = log2 n. Returns dict(valid, e, terms, ...)."""
    try:
        r, Q, R, A0 = int(params["r"]), int(params["Q"]), int(params["R"]), int(params["A0"])
        regime = params.get("regime", "paper")
        if regime not in REGIMES:
            return dict(valid=False, reason=f"unknown regime {regime}")
        delta, designated = _extra(params, regime)
    except Exception as exc:                    # malformed parameters
        return dict(valid=False, reason=f"bad params: {exc}")
    if r < 2 or Q < 1 or R < 1 or A0 < 0 or A0 + R - 1 > Q:
        return dict(valid=False, reason="r>=2, Q>=1, R>=1 and 0 <= A0 <= A0+R-1 <= Q required")
    log2D = 2 * Q * math.log2(r)
    if log2D > max_log2D(lam, regime, delta, designated) + 1e-9:
        return dict(valid=False, reason="label space D too large for the balancing regime")
    out = tail_out(Q, r, R, A0)
    rows_per_n = groups_per_row(r) * hD_over_n(lam, log2D, regime, delta, designated) * out
    terms = dict(taylor=taylor(occ_per_n(designated) * row_degree(R) / Q),
                 meas=rows_per_n * value_bits(Q),
                 groups=1.0 / r)
    return dict(valid=True, e=sum(terms.values()), terms=terms, out=out, log2D=log2D, regime=regime,
                delta=delta, designated=designated)


def exponents_over_R(lam, r, Q, Rs, regime, delta=1.0, designated=False):
    """Vectorised helper for searches: e for each interval length in Rs (best A0 each)."""
    p = (r - 1) / (2 * r)
    mu = Q * p
    best = np.ones(Rs.shape)
    bestA = np.zeros(Rs.shape, np.int64)
    for shift in (-1, 0, 1):
        A0 = np.floor(mu - (Rs - 1) / 2).astype(np.int64) + shift
        A0 = np.clip(A0, 0, np.maximum(Q - Rs + 1, 0))
        out = binom.cdf(A0 - 1, Q, p) + binom.sf(A0 + Rs - 1, Q, p)
        better = out < best
        best = np.where(better, out, best)
        bestA = np.where(better, A0, bestA)
    log2D = 2 * Q * math.log2(r)
    rows = groups_per_row(r) * hD_over_n(lam, log2D, regime, delta, designated) * best
    m = occ_per_n(designated)
    tay = np.array([taylor(m * row_degree(int(R)) / Q) for R in Rs])
    e = tay + rows * value_bits(Q) + 1.0 / r
    return e, bestA
