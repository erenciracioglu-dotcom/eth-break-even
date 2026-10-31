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
Q ladder (Q = floor(Qmax)) and the best delta is small (e is flat below delta ~ 1e-3), so a
small grid of deltas is searched.

Search (regime "tight"):
  stage 0: coarse scan (RSET x Q = Qmax x geometric R) at a reference delta -> r window;
  stage 1: cheap scan of every integer r in that window at Q = Qmax, for every delta;
  stage 2: deep refinement (full Q ladder, dense R grid, integer-R sweep) for the best few r
           of the best N_DEEP_DELTAS deltas.
If choose_params is called with an explicit regime ("paper"/"freedman"), the single-regime
three-stage search runs instead (same stages, no designated/aggregated/units extensions).

designated = True (head extension v2): every clause reads one of its own variables directly, so
only M_eff = 2n occurrences carry replacement rows; it enters through
max_log2D(lam, "tight", delta, True) and exponents_over_R(..., True). units = True (head
extension v4, needs designated): a greedy maximal matching of clauses turns whole clauses into
"units" sharing one group schedule, so M_eff <= 12n/7; it enters through
max_log2D(lam, "tight", delta, True, True) and exponents_over_R(..., True, True, True). The
optimum sits at r ~ 34 (lam = 3e5), 72 (lam = 2.2e6), 98 (lam = 5e6) and e = 0.33857 /
0.17813 / 0.13540 at those lam.

DETERMINISTIC WORK BUDGET (task T7). choose_params never reads a clock: every grid size is a
fixed integer chosen from lam alone (plus the exact r window that stage 0 produced). The sizes
come from a static cost model of assembly_public.exponents_over_R (per-call overhead a(lam) and
per-R-entry cost b(lam), measured offline and log-log interpolated) and a work budget
BUDGET_US: _cfg0 / _cfg12 shrink the grids, in order of increasing loss of accuracy, until the
*counted* work fits the budget. The same lam therefore always gives the same answer on any
machine, and the worst call over lam in [1e3, 1e13] stays under 2 s.
"""
import math
import numpy as np
from assembly_public import exponents_over_R, exponent_of, max_log2D

# stage-0 r values: dense in 12..80 where the optimum sits, sparse beyond
RSET = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64,
        68, 72, 80, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
REGIME = "tight"
AGGREGATED = True                    # head extension v3: Lemma 3.3 summed per U
UNITS = True                         # head extension v4: clause units (greedy matching), needs designated
CHUNK = 4096                         # head extension v5: connected chunk units of >= CHUNK vars

# delta grid for the tight regime (delta = 1.0 is exactly "freedman"), best deltas first so
# that truncating the grid (very large lam) keeps the ones that matter
DELTAS = (0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)
REF_DELTA = 0.01                   # delta of the stage-0 coarse scan that fixes the r window

# stage-1 settings (cheap r scan: Q = Qmax, geometric R grid)
R_LO, R_HI = 0.65, 1.50            # integer r window around the coarse r_best
# stage-2 settings (deep refinement): Q ladder = Qmax, Qmax-1, Qmax-2 + Q_FRACS of Qmax
Q_FRACS = (1.0, 0.999, 0.998, 0.995, 0.99, 0.98, 0.97, 0.95, 0.90)   # fractions of Qmax tried

# grid sizes actually used are chosen by _cfg0 / _cfg12; these are the largest (default) ones
DEFAULTS = dict(
    N_Q_COARSE=1,         # stage 0: geometric Q values per r (1 -> only Q = Qmax)
    N_R_COARSE=90,        # stage 0: geometric R grid per (r, Q)
    RSET_STRIDE=1,        # stage 0: thinning of the RSET r values
    MAX_R_WINDOW=140,     # stage 1: cap on the number of integer r values scanned
    N_R_SCAN=300,         # stage 1: geometric R grid per (r, Q)
    N_DELTAS=len(DELTAS),  # stage 1: how many deltas of the grid are scanned
    N_DELTA_R=8,          # stage 1: best r values that get scanned over the whole delta grid
    N_DEEP_DELTAS=1,      # stage 2: deltas refined deeply (cheap pass ranks them)
    N_DEEP_R=10,          # stage 2: r values refined per deep delta
    N_R_FINE=480,         # stage 2: dense geometric R grid per (r, Q)
    MAX_INT_R=300,        # stage 2: cap on the integer-R sweep per (r, Q)
)
# hard floors, reached only at very large lam (the search degrades gracefully, answers valid)
FLOORS = dict(N_Q_COARSE=1, N_R_COARSE=12, RSET_STRIDE=8, MAX_R_WINDOW=16, N_R_SCAN=24,
              N_DELTAS=1, N_DELTA_R=3, N_DEEP_DELTAS=1, N_DEEP_R=1, N_R_FINE=24, MAX_INT_R=8)
# order in which grids are shrunk: cheapest loss of accuracy first
SHRINK_ORDER = ("MAX_INT_R", "N_R_FINE", "N_DEEP_R", "N_R_SCAN", "N_DELTAS", "MAX_R_WINDOW",
                "N_R_COARSE", "RSET_STRIDE", "N_Q_COARSE")
LEGACY_N_DEEP_R = 12               # legacy ("paper"/"freedman") path: fixed counts, no clock

# ---- static cost model of assembly_public.exponents_over_R (microseconds) -------------
# measured offline on this machine: cost per call = A(lam) + B(lam) * (#R entries). B is the
# larger of the small-grid and in-search measurements (per-entry cost grows with the grid).
_COST_LAMS = (1e3, 1e4, 1e5, 3e5, 1e6, 2.2e6, 5e6, 1e7, 4.6e7, 1e8, 2.2e8, 1e9, 1e10, 1e11,
              1e12, 1e13)
_COST_A_US = (170.0, 170.0, 190.0, 190.0, 210.0, 210.0, 220.0, 230.0, 290.0, 300.0, 420.0,
              420.0, 550.0, 1250.0, 2150.0, 3850.0)
_COST_B_US = (1.0, 1.25, 1.7, 2.05, 2.7, 3.7, 5.9, 6.4, 11.8, 15.0, 17.8, 25.5, 58.0, 110.0,
              200.0, 435.0)
COST_MARGIN = 1.45                 # safety factor on the counted cost (slower machines etc.)
BUDGET_US = 1.95e6                 # counted work cap per call (~1.4 s of model work at most)
STAGE0_FRAC = 0.55                 # share of BUDGET_US that stage 0 may use


def _interp(x, xs, ys):
    """log-log interpolation, clamped (power-law extrapolation past the last point)."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        lx = (math.log(xs[-1]) - math.log(xs[-2])) / (math.log(ys[-1]) - math.log(ys[-2]))
        return ys[-1] * (x / xs[-1]) ** lx
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (math.log(x) - math.log(xs[lo])) / (math.log(xs[hi]) - math.log(xs[lo]))
    return math.exp((1 - t) * math.log(ys[lo]) + t * math.log(ys[hi]))


def _unit_cost(lam):
    """(a_us, b_us): per-call overhead and per-R-entry cost of exponents_over_R at this lam."""
    lam = max(1.0, float(lam))
    return (_interp(lam, _COST_LAMS, _COST_A_US) * COST_MARGIN,
            _interp(lam, _COST_LAMS, _COST_B_US) * COST_MARGIN)


def _Qmax(mld, r):
    return int(math.floor(mld / (2 * math.log2(r))))


def _nR(n, Q):
    """distinct R values scored by a geometric grid of n points over [1, Q + 1]"""
    return int(min(max(int(n), 1), Q + 1))


def _Qladder(Qmax):
    """Qmax and a few values just below it (deduplicated, >= 1)."""
    Qs = {Qmax, Qmax - 1, Qmax - 2}
    for f in Q_FRACS:
        Qs.add(int(round(f * Qmax)))
    return sorted(q for q in Qs if 1 <= q <= Qmax)


def _est0(lam, cfg):
    """counted work (us) of stage 0 with these grid sizes; no search is run"""
    a, b = _unit_cost(lam)
    mld = max_log2D(lam, REGIME, REF_DELTA, True, UNITS, CHUNK)
    calls = entries = 0
    for r in RSET[::cfg["RSET_STRIDE"]]:
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        nq = 1 if Qmax == 1 else min(cfg["N_Q_COARSE"] + 1, Qmax)
        calls += nq
        entries += nq * _nR(cfg["N_R_COARSE"], Qmax)
    return calls * a + entries * b


def _est12(lam, cfg, rs, designated=True, units=UNITS, chunk=CHUNK):
    """counted work (us) of stages 1-2 for the r window rs with these grid sizes"""
    a, b = _unit_cost(lam)
    calls = entries = 0
    for d in cfg["DELTAS"]:
        mld = max_log2D(lam, REGIME, d, designated, units, chunk)
        for r in rs:
            Qmax = _Qmax(mld, r)
            if Qmax < 1:
                continue
            calls += 1
            entries += _nR(cfg["N_R_SCAN"], Qmax)
    deep = rs[:cfg["N_DEEP_R"]] if cfg["N_DEEP_R"] < len(rs) else rs
    for d in cfg["DELTAS"][:cfg["N_DEEP_DELTAS"]]:
        mld = max_log2D(lam, REGIME, d, designated, units, chunk)
        for r in deep:
            Qmax = _Qmax(mld, r)
            if Qmax < 1:
                continue
            for Q in _Qladder(Qmax):
                calls += 2
                entries += _nR(cfg["N_R_FINE"], Q) + min(cfg["MAX_INT_R"], Q + 1) + 1
    return calls * a + entries * b


def _shrink(lam, cfg, fit, budget=None, order=SHRINK_ORDER):
    """shrink cfg (in `order`) until fit(lam, cfg) <= the budget share"""
    if budget is None:
        budget = STAGE0_FRAC * BUDGET_US
    for _ in range(400):
        if fit(lam, cfg) <= budget:
            return cfg
        moved = False
        for key in order:
            if cfg[key] > FLOORS[key]:
                if key == "RSET_STRIDE":
                    cfg[key] = min(FLOORS[key], cfg[key] * 2)
                elif key in ("N_DEEP_R", "N_DELTAS"):
                    cfg[key] -= 1
                else:
                    cfg[key] = max(FLOORS[key], int(math.ceil(0.8 * cfg[key])))
                moved = True
                break
        if not moved:
            return cfg
    return cfg


def _cfg0(lam):
    cfg = dict(DEFAULTS)
    return _shrink(lam, cfg, _est0)


def _cfg12(lam, rs, est0_us, designated=True, units=UNITS, deltas=DELTAS, chunk=CHUNK):
    """grid sizes for stages 1-2; rs is already fixed, so only the stage-1/2 grids shrink"""
    cfg = dict(DEFAULTS)
    cfg["DELTAS"] = deltas
    for key in ("N_Q_COARSE", "N_R_COARSE", "RSET_STRIDE", "MAX_R_WINDOW"):
        cfg[key] = FLOORS[key]          # these no longer cost anything here
    return _shrink(lam, cfg,
                   lambda L, c: est0_us + _est12(L, c, rs, designated, units, chunk),
                   budget=BUDGET_US,
                   order=("MAX_INT_R", "N_R_FINE", "N_DEEP_R", "N_R_SCAN", "N_DELTAS"))


def _mk(r, Q, R, A0, regime, delta, designated=False, aggregated=AGGREGATED, units=UNITS,
      chunk=CHUNK):
    p = dict(r=int(r), Q=int(Q), R=int(R), A0=int(A0), regime=regime)
    if regime == "tight":
        p["delta"] = float(delta)
    if designated:
        p["designated"] = True
    if aggregated:
        p["aggregated"] = True
    if units:
        p["units"] = True
    if units and chunk:
        p["chunk"] = int(chunk)
    return p


def _scan(lam, regime, delta, r, Q, Rs, best, designated=False, aggregated=AGGREGATED,
        units=UNITS, chunk=CHUNK):
    """Score the interval lengths Rs for (r, Q); keep the best candidate in `best`."""
    Rs = np.unique(Rs.astype(np.int64))
    Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
    if Rs.size == 0:
        return best
    e, A0 = exponents_over_R(lam, r, Q, Rs, regime, delta, designated, aggregated, units,
                              chunk)
    i = int(np.argmin(e))
    if best is None or e[i] < best[0]:
        return (float(e[i]), _mk(r, Q, int(Rs[i]), int(A0[i]), regime, delta, designated,
                                 aggregated, units, chunk))
    return best


def _coarse(lam, regime, delta, mld, cfg, designated=False, aggregated=AGGREGATED,
            units=UNITS, chunk=CHUNK):
    """Stage 0: coarse scan over RSET x geometric Q x geometric R (guarantees a valid answer)."""
    best = None
    for r in RSET[::cfg["RSET_STRIDE"]]:
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        if Qmax == 1:
            Qs = np.array([1], dtype=np.int64)
        else:
            Qs = np.unique(np.concatenate([np.round(np.geomspace(2, Qmax, cfg["N_Q_COARSE"])),
                                           [Qmax]]).astype(np.int64))
        for Q in Qs:
            Q = int(Q)
            Rs = np.unique(np.round(np.geomspace(1, Q + 1, cfg["N_R_COARSE"])).astype(np.int64))
            best = _scan(lam, regime, delta, r, Q, Rs, best, designated, aggregated, units,
                         chunk)
    return best


def _cheap_r_scan(lam, regime, delta, mld, rs, cfg, designated=False, aggregated=AGGREGATED,
                  units=UNITS, chunk=CHUNK):
    """Stage 1: Q = Qmax, moderate R grid, for every r in rs. Returns [(e, r, params)]."""
    out = []
    for r in rs:
        Qmax = _Qmax(mld, r)
        if Qmax < 1:
            continue
        Q = Qmax
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, cfg["N_R_SCAN"])).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A0 = exponents_over_R(lam, r, Q, Rs, regime, delta, designated, aggregated,
                                  units, chunk)
        i = int(np.argmin(e))
        out.append((float(e[i]), int(r),
                    _mk(r, Q, int(Rs[i]), int(A0[i]), regime, delta, designated, aggregated,
                        units, chunk)))
    return out


def _deep(lam, regime, delta, r, Qmax, best, cfg, designated=False,
          aggregated=AGGREGATED, units=UNITS, chunk=CHUNK):
    """Stage 2: full Q ladder, dense R grid and an integer-R sweep for one r."""
    for Q in _Qladder(Qmax):
        Rs = np.unique(np.round(np.geomspace(1, Q + 1, cfg["N_R_FINE"])).astype(np.int64))
        Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
        if Rs.size == 0:
            continue
        e, A = exponents_over_R(lam, r, Q, Rs, regime, delta, designated, aggregated,
                                units, chunk)
        i = int(np.argmin(e))
        if best is None or e[i] < best[0]:
            best = (float(e[i]),
                    _mk(r, Q, int(Rs[i]), int(A[i]), regime, delta, designated, aggregated,
                        units, chunk))
        # every integer R between the grid neighbours of the best R
        lo = int(Rs[i - 1]) if i > 0 else 1
        hi = int(Rs[i + 1]) if i + 1 < Rs.size else int(Q + 1)
        lo = max(1, min(lo, int(Rs[i]) - cfg["MAX_INT_R"] // 2))
        hi = min(int(Q + 1), max(hi, int(Rs[i]) + cfg["MAX_INT_R"] // 2))
        hi = min(hi, lo + cfg["MAX_INT_R"])
        if hi > lo:
            best = _scan(lam, regime, delta, r, Q, np.arange(lo, hi + 1, dtype=np.int64), best,
                         designated, aggregated, units, chunk)
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


def _r_window(r_best, cfg):
    rs = list(range(max(2, int(math.floor(R_LO * r_best))),
                    max(2, int(math.ceil(R_HI * r_best))) + 1))
    if len(rs) > cfg["MAX_R_WINDOW"]:      # keep the window centred, strided
        stride = math.ceil(len(rs) / cfg["MAX_R_WINDOW"])
        rs = rs[::stride]
    return rs


def _single_regime(lam, regime):
    """Legacy path: explicit regime ("paper" / "freedman"), three-stage search, fixed counts."""
    cfg = _cfg0(lam)
    mld = max_log2D(lam, regime)
    best = _coarse(lam, regime, 1.0, mld, cfg, False, False, False, 0)
    if best is None:
        r0 = int(2.0 ** (mld / 2.0))
        if r0 < 2 or 2 * math.log2(r0) > mld + 1e-9:
            return dict(r=2, Q=1, R=1, A0=0, regime=regime)   # regime allows no labels at all
        return dict(r=int(r0), Q=1, R=1, A0=0, regime=regime)
    rs = _r_window(best[1]["r"], cfg)
    cfg12 = _cfg12(lam, rs, _est0(lam, cfg), False, False, (1.0,), 0)
    per_r = _cheap_r_scan(lam, regime, 1.0, mld, rs, cfg12, False, False, False, 0)
    for _e, _r, p in per_r:
        if best is None or _e < best[0]:
            best = (_e, p)
    deep_rs = [r for _e, r, _p in sorted(per_r, key=lambda t: t[0])[:LEGACY_N_DEEP_R]]
    if best[1]["r"] not in deep_rs:
        deep_rs.append(best[1]["r"])
    for r in deep_rs:
        best = _deep(lam, regime, 1.0, r, _Qmax(mld, r), best, cfg12, False, False, False,
                   0)
    return _polish(lam, best, regime, 1.0)[1]


def choose_params(lam, regime=None):
    if regime is not None and regime != "tight":
        return _single_regime(lam, regime)

    regime, delta_ref = "tight", REF_DELTA
    designated = True                      # head extension v2: designated occurrences
    units = UNITS                          # head extension v4: clause units
    chunk = CHUNK                          # head extension v5: connected chunk units
    cfg = _cfg0(lam)
    mld_ref = max_log2D(lam, regime, delta_ref, designated, units, chunk)

    # ---------- stage 0: coarse scan at the reference delta (valid answer + r window) -------
    best = _coarse(lam, regime, delta_ref, mld_ref, cfg, designated, AGGREGATED, units,
                 chunk)
    if best is None:
        r0 = int(2.0 ** (mld_ref / 2.0))
        if r0 < 2 or 2 * math.log2(r0) > mld_ref + 1e-9:
            return dict(r=2, Q=1, R=1, A0=0, regime=regime, delta=1.0, designated=True,
                        aggregated=AGGREGATED, units=units, chunk=chunk)   # no labels
        return dict(r=int(r0), Q=1, R=1, A0=0, regime=regime, delta=delta_ref, designated=True,
                    aggregated=AGGREGATED, units=units, chunk=chunk)
    rs = _r_window(best[1]["r"], cfg)
    est0_us = _est0(lam, cfg)
    cfg = _cfg12(lam, rs, est0_us)     # stages 1-2 sized from the exact r window

    # ---------- stage 1: cheap Q = Qmax r scan for every delta in the grid ----------
    per_delta = []                            # (e, delta, per_r list)
    for d in DELTAS[:cfg["N_DELTAS"]]:
        mld = max_log2D(lam, regime, d, designated, units, chunk)
        per_r = _cheap_r_scan(lam, regime, d, mld, rs, cfg, designated, AGGREGATED, units,
                              chunk)
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
    for _e0, d, per_r in per_delta[:cfg["N_DEEP_DELTAS"]]:
        mld = max_log2D(lam, regime, d, designated, units, chunk)
        deep_rs = [r for _e, r, _p in sorted(per_r, key=lambda t: t[0])[:cfg["N_DEEP_R"]]]
        if best[1].get("delta") == d and best[1]["r"] not in deep_rs:
            deep_rs.append(best[1]["r"])
        for r in deep_rs:
            best = _deep(lam, regime, d, r, _Qmax(mld, r), best, cfg, designated, AGGREGATED,
                         units, chunk)

    return _polish(lam, best, regime, best[1].get("delta", delta_ref))[1]
