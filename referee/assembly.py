"""FROZEN cost assembly (the referee scores with this file).

Bounded-occurrence 3-CNF (each variable in <= 3 clauses, clauses of width <= 3), counting
algorithm of Sections 2-4 of "The Exponential Time Hypothesis is False". Exponents per
variable: cost ~ 2^(e*n), polynomial factors ignored. lam = log2 n, M = 3n occurrences.

algo/assembly_public.py is an identical copy for the agent's own parameter search.
Editing that copy does not change any score.

Regimes (how many labels D may be balanced and how many copies h that costs):
  paper    : Lemma 4.5 as written, 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M)
  freedman : PROVISIONAL head extension (not in the paper): D ln D <= M / (4(b+w)),
             h <= 2M/D, by a Freedman/Bernstein bound per label + union bound,
             derandomised with pessimistic estimators. Must stay marked provisional.
"""
import math
import numpy as np
from scipy.stats import binom

LOG2_3 = math.log2(3)
TARGETS = {"brute": 1.0, "sharp3sat": math.log2(1.637), "ppsz": math.log2(1.307031578)}
M_PER_N = 3       # occurrences per variable, worst case of the class
BW = 6            # b + w (b <= 3 factors per variable, w <= 3 variables per factor)
REGIMES = ("paper", "freedman")


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


def max_log2D(lam, regime):
    """Largest allowed log2 D for n = 2^lam and M = 3n."""
    logM = lam + math.log2(M_PER_N)
    if regime == "paper":                       # 2(b+w) D^2 <= M
        return (logM - math.log2(2 * BW)) / 2
    if regime == "freedman":                    # D ln D <= M / (4(b+w))
        target = logM - math.log2(4 * BW)
        x = max(target, 1.0)
        for _ in range(200):                    # solve x + log2(x ln 2) = target
            x = target - math.log2(max(x * math.log(2), 1e-12))
        return x
    raise ValueError(regime)


def hD_over_n(lam, log2D, regime):
    """h * D / n with the regime's bound on the number of copies h."""
    if regime == "paper":                       # (M/D + sqrt(2(b+w)M)) D / n = 3 + 6 D / sqrt(n)
        ex = log2D - lam / 2.0
        return M_PER_N + math.sqrt(2 * BW * M_PER_N) * (2.0 ** ex if ex > -300 else 0.0)
    if regime == "freedman":                    # 2M/D * D / n
        return 2.0 * M_PER_N
    raise ValueError(regime)


def tail_out(Q, r, R, A0):
    """Fraction of basis vectors that must be measured: Pr[W outside [A0, A0+R-1]], W ~ Bin(Q, p)."""
    p = (r - 1) / (2 * r)
    return float(binom.cdf(A0 - 1, Q, p) + binom.sf(A0 + R - 1, Q, p))


def exponent_of(params, lam):
    """Score one parameter set at lam = log2 n. Returns dict(valid, e, terms, ...)."""
    try:
        r, Q, R, A0 = int(params["r"]), int(params["Q"]), int(params["R"]), int(params["A0"])
        regime = params.get("regime", "paper")
    except Exception as exc:                    # malformed parameters
        return dict(valid=False, reason=f"bad params: {exc}")
    if regime not in REGIMES:
        return dict(valid=False, reason=f"unknown regime {regime}")
    if r < 2 or Q < 1 or R < 1 or A0 < 0 or A0 + R - 1 > Q:
        return dict(valid=False, reason="r>=2, Q>=1, R>=1 and 0 <= A0 <= A0+R-1 <= Q required")
    log2D = 2 * Q * math.log2(r)
    if log2D > max_log2D(lam, regime) + 1e-9:
        return dict(valid=False, reason="label space D too large for the balancing regime")
    out = tail_out(Q, r, R, A0)
    rows_per_n = groups_per_row(r) * hD_over_n(lam, log2D, regime) * out
    terms = dict(taylor=taylor(M_PER_N * row_degree(R) / Q),
                 meas=rows_per_n * value_bits(Q),
                 groups=1.0 / r)
    return dict(valid=True, e=sum(terms.values()), terms=terms, out=out, log2D=log2D, regime=regime)


def exponents_over_R(lam, r, Q, Rs, regime):
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
    rows = groups_per_row(r) * hD_over_n(lam, log2D, regime) * best
    tay = np.array([taylor(M_PER_N * row_degree(int(R)) / Q) for R in Rs])
    e = tay + rows * value_bits(Q) + 1.0 / r
    return e, bestA
