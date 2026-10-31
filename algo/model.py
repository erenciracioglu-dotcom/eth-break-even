"""Parameter choice for the cost model (agent-editable).

The referee calls choose_params(lam) for many values of lam = log2 n and scores each result
with its own frozen copy of assembly_public.py. The score is the smallest lam at which the
cost exponent e drops below a target (best 3-SAT: log2(1.307031578) = 0.38630). Lower is better.

Rules: choose_params must return a dict with integer keys r, Q, R, A0 and a "regime" string
("paper" or "freedman"), must be deterministic, and must take under 3 seconds per call.

Regime used: "freedman" (provisional head extension of Lemma 4.5, D ln D <= M / (4(b+w)),
h <= 2M/D). It allows log2 D ~ log2 M instead of ~log2 M / 2, so Qmax roughly doubles and the
optimum moves to a somewhat larger r than in the paper regime (r ~ 15 at lam = 3e5, ~34 at
lam = 2.2e6, ~48 at lam = 5e6, ~456 at lam = 1e9), with Q at the very top of the ladder
(Q = floor(Qmax)). The search ranges below cover that: stage 1 has dense r values in 20..80,
stage 2 scans every integer r in a wide window around the coarse optimum and then refines the
best few r deeply (Q ladder up to Qmax, dense R grid, integer R sweep).

Three-stage search:
  stage 1: coarse scan over RSET x geometric Q x geometric R (always runs, so every lam gets
           a valid answer);
  stage 2: every integer r in [0.65 r_best, 1.5 r_best], Q = Qmax only, moderate R grid;
  stage 3: the best few r from stage 2, full Q ladder, dense R grid, then every integer R
           between the grid neighbours of the best R.
Every candidate is scored by exponents_over_R / exponent_of from assembly_public.
"""
import math
import time
import numpy as np
from assembly_public import exponents_over_R, exponent_of, max_log2D

# stage-1 r values: dense in 20..80 where the freedman optimum sits, sparse beyond
RSET = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 64, 72,
        80, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
REGIME = "freedman"

# stage-1 grid sizes (coarse; stages 2-3 do the fine work)
N_Q_COARSE = 24
N_R_COARSE = 90

# stage-2 settings (cheap r scan: Q = Qmax, moderate R grid)
R_LO, R_HI = 0.65, 1.50          # integer r window around the coarse r_best
MAX_R_WINDOW = 140               # cap on the number of integer r values scanned
N_R_SCAN = 300                   # R grid for the r scan
N_DEEP_R = 12                    # r values refined in stage 3

# stage-3 settings (deep refinement)
Q_FRACS = (1.0, 0.999, 0.998, 0.995, 0.99, 0.98, 0.97, 0.95, 0.90)   # fractions of Qmax tried
N_R_FINE = 480                   # dense R grid per (r, Q)
MAX_INT_R = 300                  # cap on the integer-R refinement per (r, Q)
DEADLINE = 2.0                   # seconds; stages 2-3 stop early past this (limit is 3 s)


def _Qladder(Qmax):
    """Qmax and a few values just below it (deduplicated, >= 1)."""
    Qs = {Qmax, Qmax - 1, Qmax - 2}
    for f in Q_FRACS:
        Qs.add(int(round(f * Qmax)))
    return sorted(q for q in Qs if 1 <= q <= Qmax)


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


def _Qmax(mld, r):
    return int(math.floor(mld / (2 * math.log2(r))))


def _deep(lam, regime, r, Qmax, best, t_start):
    """Stage 3: full Q ladder, dense R grid and an integer-R sweep for one r."""
    for Q in _Qladder(Qmax):
        if time.time() - t_start > DEADLINE:
            return best
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_FINE)).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A = exponents_over_R(lam, r, Q, Rs, regime)
        i = int(np.argmin(e))
        if best is None or e[i] < best[0]:
            best = (float(e[i]), dict(r=int(r), Q=int(Q), R=int(Rs[i]), A0=int(A[i]),
                                      regime=regime))
        # every integer R between the grid neighbours of the best R
        lo = int(Rs[i - 1]) if i > 0 else 1
        hi = int(Rs[i + 1]) if i + 1 < Rs.size else int(Q + 1)
        lo = max(1, min(lo, int(Rs[i]) - MAX_INT_R // 2))
        hi = min(int(Q + 1), max(hi, int(Rs[i]) + MAX_INT_R // 2))
        hi = min(hi, lo + MAX_INT_R)
        if hi > lo:
            best = _scan(lam, regime, r, Q, np.arange(lo, hi + 1, dtype=np.int64), best)
    return best


def choose_params(lam, regime=None):
    t_start = time.time()
    regime = regime or REGIME
    mld = max_log2D(lam, regime)

    # ---------- stage 1: coarse scan (guarantees a valid answer) ----------
    best = None
    for r in RSET:
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        if Qmax == 1:
            Qs = np.array([1], dtype=np.int64)
        else:
            Qs = np.unique(np.concatenate([np.round(np.geomspace(2, Qmax, N_Q_COARSE)),
                                           [Qmax]]).astype(np.int64))
        for Q in Qs:
            Q = int(Q)
            Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_COARSE)).astype(np.int64))
            best = _scan(lam, regime, r, Q, Rs, best)
    if best is None:
        # smallest possible label space: Q = 1 with r as small as the regime allows
        r0 = int(2.0 ** (mld / 2.0))
        if r0 < 2 or 2 * math.log2(r0) > mld + 1e-9:
            return dict(r=2, Q=1, R=1, A0=0, regime=regime)   # regime allows no labels at all
        return dict(r=int(r0), Q=1, R=1, A0=0, regime=regime)
    r_best = best[1]["r"]

    # ---------- stage 2: cheap scan of every integer r near the coarse optimum ----------
    r_lo = max(2, int(math.floor(R_LO * r_best)))
    r_hi = max(r_lo, int(math.ceil(R_HI * r_best)))
    rs = list(range(r_lo, r_hi + 1))
    if len(rs) > MAX_R_WINDOW:            # keep the window centred, strided
        stride = math.ceil(len(rs) / MAX_R_WINDOW)
        rs = rs[::stride]
    per_r = []                            # (e, r) for each scanned r
    for r in rs:
        if time.time() - t_start > DEADLINE:
            break
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        Q = Qmax
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_SCAN)).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A0 = exponents_over_R(lam, r, Q, Rs, regime)
        i = int(np.argmin(e))
        per_r.append((float(e[i]), int(r)))
        if best is None or e[i] < best[0]:
            best = (float(e[i]), dict(r=int(r), Q=int(Q), R=int(Rs[i]), A0=int(A0[i]),
                                      regime=regime))

    # ---------- stage 3: deep refinement of the best few r ----------
    per_r.sort()
    deep_rs = [r for _e, r in per_r[:N_DEEP_R]]
    if r_best not in deep_rs:
        deep_rs.append(r_best)
    for r in deep_rs:
        if time.time() - t_start > DEADLINE:
            break
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        best = _deep(lam, regime, r, Qmax, best, t_start)

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
