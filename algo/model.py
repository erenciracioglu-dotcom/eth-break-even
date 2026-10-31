"""Parameter choice for the cost model (agent-editable).

The referee calls choose_params(lam) for many values of lam = log2 n and scores each result
with its own frozen copy of assembly_public.py. The score is the smallest lam at which the
cost exponent e drops below a target (best 3-SAT: log2(1.307031578) = 0.38630). Lower is better.

Rules: choose_params must return a dict with integer keys r, Q, R, A0, a "regime" string and,
for regime "tight", a float "delta"; it must be deterministic and take under 3 seconds.

Regimes used: "tight" (head extension, audited) by default - the freedman argument with a
tolerance delta: D ln D <= delta^2 M / (4(b+w)) and h <= (1 + delta) M/D. delta = 1 is exactly
"freedman". Shrinking delta costs only 2*log2(1/delta) bits of log2 D (negligible at these
sizes, ~8 bits for delta = 1e-3 out of ~2.2e6) but multiplies the copy count h D / n from 3 to
1 + delta, which roughly halves the measurement term. The optimum sits at the very top of the
Q ladder (Q = floor(Qmax)) with r ~ 15-16 at lam = 3e5 and r ~ 35 at lam = 2.2e6, and the best
delta is small (e is flat below delta ~ 1e-3), so a small grid of deltas is searched.

Search (regime "tight"):
  stage 0: one coarse scan (RSET x geometric Q x geometric R) at a reference delta -> r window;
  stage 1: cheap scan of every integer r in that window at Q = Qmax, for every delta;
  stage 2: deep refinement (full Q ladder, dense R grid, integer-R sweep) for the best few r
           of the best N_DEEP_DELTAS deltas only, so the whole call stays under 3 s.
If choose_params is called with an explicit regime ("paper"/"freedman"), the old single-regime
three-stage search runs instead. designated is never used here.
"""
import math
import time
import numpy as np
from assembly_public import exponents_over_R, exponent_of, max_log2D

# stage-0 r values: dense in 12..80 where the optimum sits, sparse beyond
RSET = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 64, 72,
        80, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
REGIME = "tight"

# delta grid for the tight regime (delta = 1.0 is exactly "freedman")
DELTAS = (1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0001)
REF_DELTA = 0.01                   # delta of the stage-0 coarse scan that fixes the r window
N_DEEP_DELTAS = 2                  # deltas refined deeply (cheap pass ranks them)

# stage-0 grid sizes (coarse; the deep stages do the fine work)
N_Q_COARSE = 24
N_R_COARSE = 90

# stage-1 settings (cheap r scan: Q = Qmax, moderate R grid)
R_LO, R_HI = 0.65, 1.50            # integer r window around the coarse r_best
MAX_R_WINDOW = 140                 # cap on the number of integer r values scanned
N_R_SCAN = 300                     # R grid for the r scan
N_DEEP_R = 12                      # r values refined per deep delta

# stage-2 settings (deep refinement)
Q_FRACS = (1.0, 0.999, 0.998, 0.995, 0.99, 0.98, 0.97, 0.95, 0.90)   # fractions of Qmax tried
N_R_FINE = 480                     # dense R grid per (r, Q)
MAX_INT_R = 300                    # cap on the integer-R refinement per (r, Q)
DEADLINE = 1.6                     # seconds; later stages stop early past this (limit is 3 s)


def _Qladder(Qmax):
    """Qmax and a few values just below it (deduplicated, >= 1)."""
    Qs = {Qmax, Qmax - 1, Qmax - 2}
    for f in Q_FRACS:
        Qs.add(int(round(f * Qmax)))
    return sorted(q for q in Qs if 1 <= q <= Qmax)


def _Qmax(mld, r):
    return int(math.floor(mld / (2 * math.log2(r))))


def _mk(r, Q, R, A0, regime, delta):
    p = dict(r=int(r), Q=int(Q), R=int(R), A0=int(A0), regime=regime)
    if regime == "tight":
        p["delta"] = float(delta)
    return p


def _scan(lam, regime, delta, r, Q, Rs, best):
    """Score the interval lengths Rs for (r, Q); keep the best candidate in `best`."""
    Rs = np.unique(Rs.astype(np.int64))
    Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
    if Rs.size == 0:
        return best
    e, A0 = exponents_over_R(lam, r, Q, Rs, regime, delta)
    i = int(np.argmin(e))
    if best is None or e[i] < best[0]:
        return (float(e[i]), _mk(r, Q, int(Rs[i]), int(A0[i]), regime, delta))
    return best


def _coarse(lam, regime, delta, mld):
    """Stage 0: coarse scan over RSET x geometric Q x geometric R (guarantees a valid answer)."""
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
            best = _scan(lam, regime, delta, r, Q, Rs, best)
    return best


def _cheap_r_scan(lam, regime, delta, mld, rs):
    """Stage 1: Q = Qmax, moderate R grid, for every r in rs. Returns [(e, r, params)]."""
    out = []
    for r in rs:
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        Q = Qmax
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_SCAN)).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A0 = exponents_over_R(lam, r, Q, Rs, regime, delta)
        i = int(np.argmin(e))
        out.append((float(e[i]), int(r), _mk(r, Q, int(Rs[i]), int(A0[i]), regime, delta)))
    return out


def _deep(lam, regime, delta, r, Qmax, best, t_start):
    """Stage 2: full Q ladder, dense R grid and an integer-R sweep for one r."""
    for Q in _Qladder(Qmax):
        if time.time() - t_start > DEADLINE:
            return best
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, N_R_FINE)).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A = exponents_over_R(lam, r, Q, Rs, regime, delta)
        i = int(np.argmin(e))
        if best is None or e[i] < best[0]:
            best = (float(e[i]), _mk(r, Q, int(Rs[i]), int(A[i]), regime, delta))
        # every integer R between the grid neighbours of the best R
        lo = int(Rs[i - 1]) if i > 0 else 1
        hi = int(Rs[i + 1]) if i + 1 < Rs.size else int(Q + 1)
        lo = max(1, min(lo, int(Rs[i]) - MAX_INT_R // 2))
        hi = min(int(Q + 1), max(hi, int(Rs[i]) + MAX_INT_R // 2))
        hi = min(hi, lo + MAX_INT_R)
        if hi > lo:
            best = _scan(lam, regime, delta, r, Q, np.arange(lo, hi + 1, dtype=np.int64), best)
    return best


def _polish(lam, best, regime, delta):
    """Final check: a couple of A0 shifts around the winner."""
    p = dict(best[1])
    for d in (-2, -1, 0, 1, 2):
        cand = dict(p)
        cand["A0"] = int(p["A0"]) + d
        if cand["A0"] < 0 or cand["A0"] + cand["R"] - 1 > cand["Q"]:
            continue
        s = exponent_of(cand, lam)
        if s.get("valid") and s["e"] < best[0]:
            best = (float(s["e"]), cand)
    return best


def _r_window(r_best):
    rs = list(range(max(2, int(math.floor(R_LO * r_best))),
                    max(2, int(math.ceil(R_HI * r_best))) + 1))
    if len(rs) > MAX_R_WINDOW:                 # keep the window centred, strided
        stride = math.ceil(len(rs) / MAX_R_WINDOW)
        rs = rs[::stride]
    return rs


def _single_regime(lam, regime, t_start):
    """Legacy path: explicit regime ("paper" / "freedman"), three-stage search."""
    mld = max_log2D(lam, regime)
    best = _coarse(lam, regime, 1.0, mld)
    if best is None:
        r0 = int(2.0 ** (mld / 2.0))
        if r0 < 2 or 2 * math.log2(r0) > mld + 1e-9:
            return dict(r=2, Q=1, R=1, A0=0, regime=regime)   # regime allows no labels at all
        return dict(r=int(r0), Q=1, R=1, A0=0, regime=regime)
    rs = _r_window(best[1]["r"])
    per_r = [(e, r, p) for (e, r, p) in _cheap_r_scan(lam, regime, 1.0, mld, rs)
             if time.time() - t_start <= DEADLINE]
    for _e, _r, p in per_r:
        if best is None or _e < best[0]:
            best = (_e, p)
    deep_rs = [r for _e, r, _p in sorted(per_r, key=lambda t: t[0])[:N_DEEP_R]]
    if best[1]["r"] not in deep_rs:
        deep_rs.append(best[1]["r"])
    for r in deep_rs:
        if time.time() - t_start > DEADLINE:
            break
        best = _deep(lam, regime, 1.0, r, _Qmax(mld, r), best, t_start)
    return _polish(lam, best, regime, 1.0)[1]


def choose_params(lam, regime=None):
    t_start = time.time()
    if regime is not None and regime != "tight":
        return _single_regime(lam, regime, t_start)

    regime, delta_ref = "tight", REF_DELTA
    mld_ref = max_log2D(lam, regime, delta_ref)

    # ---------- stage 0: coarse scan at the reference delta (valid answer + r window) -------
    best = _coarse(lam, regime, delta_ref, mld_ref)
    if best is None:
        r0 = int(2.0 ** (mld_ref / 2.0))
        if r0 < 2 or 2 * math.log2(r0) > mld_ref + 1e-9:
            return dict(r=2, Q=1, R=1, A0=0, regime=regime, delta=1.0)  # no labels allowed
        return dict(r=int(r0), Q=1, R=1, A0=0, regime=regime, delta=delta_ref)
    rs = _r_window(best[1]["r"])

    # ---------- stage 1: cheap Q = Qmax r scan for every delta in the grid ----------
    per_delta = []                            # (e, delta, per_r list)
    for d in DELTAS:
        if time.time() - t_start > DEADLINE:
            break
        mld = max_log2D(lam, regime, d)
        per_r = _cheap_r_scan(lam, regime, d, mld, rs)
        if not per_r:
            continue
        per_delta.append((min(e for e, _r, _p in per_r), d, per_r))
        for e, _r, p in per_r:
            if e < best[0]:
                best = (e, p)
    if not per_delta:
        return _polish(lam, best, regime, delta_ref)[1]

    # ---------- stage 2: deep refinement of the best r for the best few deltas ----------
    per_delta.sort(key=lambda t: t[0])
    for _e0, d, per_r in per_delta[:N_DEEP_DELTAS]:
        mld = max_log2D(lam, regime, d)
        deep_rs = [r for _e, r, _p in sorted(per_r, key=lambda t: t[0])[:N_DEEP_R]]
        if best[1].get("delta") == d and best[1]["r"] not in deep_rs:
            deep_rs.append(best[1]["r"])
        for r in deep_rs:
            if time.time() - t_start > DEADLINE:
                break
            best = _deep(lam, regime, d, r, _Qmax(mld, r), best, t_start)

    return _polish(lam, best, regime, best[1].get("delta", delta_ref))[1]
