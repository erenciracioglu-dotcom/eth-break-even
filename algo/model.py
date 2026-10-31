"""Parameter choice for the cost model (agent-editable).

The referee calls choose_params(lam) for many values of lam = log2 n and scores each result
with its own frozen copy of assembly_public.py. The score is the smallest lam at which the
cost exponent e drops below a target (best 3-SAT: log2(1.307031578) = 0.38630). Lower is better.

Rules: choose_params must return a dict with integer keys r, Q, R, A0 and a "regime" string
("paper" or "freedman"), must be deterministic, and must take under 3 seconds per call.

Two-stage search:
  stage 1: coarse scan over RSET x 36 geometric Q x 120 geometric R (always runs, so every
           lam gets a valid answer);
  stage 2: every integer r in [0.75 r_best, 1.35 r_best], Q = Qmax and a few values just
           below it, a dense R grid from exponents_over_R, then every integer R between the
           grid neighbours of the best R.
Every candidate is scored by exponents_over_R / exponent_of from assembly_public.
"""
import math
import time
import numpy as np
from assembly_public import exponents_over_R, exponent_of, max_log2D

RSET = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536,
        2048, 3072, 4096]
REGIME = "paper"

# stage-2 settings
R_LO, R_HI = 0.75, 1.35          # integer r window around the coarse r_best
Q_FRACS = (1.0, 0.995, 0.99, 0.98, 0.97, 0.95, 0.90)   # fractions of Qmax tried
N_R_FINE = 320                   # dense R grid per (r, Q)
MAX_INT_R = 220                  # cap on the integer-R refinement per (r, Q)
MAX_R_WINDOW = 40                # cap on the number of integer r values tried
DEADLINE = 1.8                   # seconds; stage 2 stops early past this (limit is 3 s)


def _Qladder(Qmax):
    """Qmax and a few values just below it (deduplicated, >= 2)."""
    Qs = {Qmax, Qmax - 1, Qmax - 2}
    for f in Q_FRACS:
        Qs.add(int(round(f * Qmax)))
    return sorted(q for q in Qs if 2 <= q <= Qmax)


def _scan(lam, regime, r, Q, Rs, best):
    """Score the interval lengths Rs for (r, Q); keep the best candidate in `best`."""
    Rs = np.unique(Rs.astype(np.int64))
    Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
    if Rs.size == 0:
        return best
    e, A0 = exponents_over_R(lam, r, Q, Rs, regime)
    i = int(np.argmin(e))
    if best is None or e[i] < best[0]:
        return (float(e[i]), dict(r=int(r), Q=int(Q), R=int(Rs[i]), A0=int(A0[i]),
                                  regime=regime))
    return best


def choose_params(lam, regime=None):
    t_start = time.time()
    regime = regime or REGIME
    mld = max_log2D(lam, regime)

    # ---------- stage 1: coarse scan (guarantees a valid answer) ----------
    best = None
    r_best = None
    for r in RSET:
        log2r = math.log2(r)
        Qmax = int(math.floor(mld / (2 * log2r)))
        if Qmax < 2:
            continue
        Qs = np.unique(np.concatenate([np.round(np.geomspace(2, Qmax, 36)), [Qmax]]).astype(np.int64))
        for Q in Qs:
            Q = int(Q)
            Rs = np.unique(np.round(np.geomspace(1, Q + 1, 120)).astype(np.int64))
            best = _scan(lam, regime, r, Q, Rs, best)
    if best is None:
        return dict(r=2, Q=1, R=1, A0=0, regime=regime)
    r_best = best[1]["r"]

    # ---------- stage 2: every integer r near the coarse optimum ----------
    r_lo = max(2, int(math.floor(R_LO * r_best)))
    r_hi = max(r_lo, int(math.ceil(R_HI * r_best)))
    rs = list(range(r_lo, r_hi + 1))
    if len(rs) > MAX_R_WINDOW:            # keep the window centred, strided
        stride = math.ceil(len(rs) / MAX_R_WINDOW)
        rs = rs[::stride]
    for r in rs:
        if time.time() - t_start > DEADLINE:
            break
        log2r = math.log2(r)
        Qmax = int(math.floor(mld / (2 * log2r)))
        if Qmax < 2:
            continue
        for Q in _Qladder(Qmax):
            if time.time() - t_start > DEADLINE:
                break
            # dense grid of interval lengths
            Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_FINE)).astype(np.int64))
            Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
            if Rs.size == 0:
                continue
            e, _A = exponents_over_R(lam, r, Q, Rs, regime)
            i = int(np.argmin(e))
            if best is None or e[i] < best[0]:
                best = (float(e[i]), dict(r=int(r), Q=int(Q), R=int(Rs[i]), A0=int(_A[i]),
                                          regime=regime))
            # every integer R between the grid neighbours of the best R
            lo = int(Rs[i - 1]) if i > 0 else 1
            hi = int(Rs[i + 1]) if i + 1 < Rs.size else int(Q + 1)
            lo = max(1, min(lo, int(Rs[i]) - MAX_INT_R // 2))
            hi = min(int(Q + 1), max(hi, int(Rs[i]) + MAX_INT_R // 2))
            hi = min(hi, lo + MAX_INT_R)
            if hi > lo:
                best = _scan(lam, regime, r, Q, np.arange(lo, hi + 1, dtype=np.int64), best)

    # ---------- final check: a couple of A0 shifts around the winner ----------
    p = dict(best[1])
    for d in (-2, -1, 0, 1, 2):
        cand = dict(p)
        cand["A0"] = int(p["A0"]) + d
        if cand["A0"] < 0 or cand["A0"] + cand["R"] - 1 > cand["Q"]:
            continue
        s = exponent_of(cand, lam)
        if s.get("valid") and s["e"] < best[0]:
            best = (float(s["e"]), cand)
    return best[1]
