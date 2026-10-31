"""Parameter choice for the cost model (agent-editable).

The referee calls choose_params(lam) for many values of lam = log2 n and scores each result
with its own frozen copy of assembly_public.py. The score is the smallest lam at which the
cost exponent e drops below a target (best 3-SAT: log2(1.307031578) = 0.38630). Lower is better.

Rules: choose_params must return a dict with integer keys r, Q, R, A0 and a "regime" string
("paper" or "freedman"), must be deterministic, and must take under 3 seconds per call.
"""
import math
import numpy as np
from assembly_public import exponents_over_R, max_log2D

RSET = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536,
        2048, 3072, 4096]
REGIME = "paper"


def choose_params(lam, regime=None):
    regime = regime or REGIME
    best = None
    for r in RSET:
        log2r = math.log2(r)
        Qmax = int(math.floor(max_log2D(lam, regime) / (2 * log2r)))
        if Qmax < 2:
            continue
        Qs = np.unique(np.concatenate([np.round(np.geomspace(2, Qmax, 36)), [Qmax]]).astype(np.int64))
        for Q in Qs:
            Q = int(Q)
            Rs = np.unique(np.round(np.geomspace(1, Q + 1, 120)).astype(np.int64))
            Rs = Rs[(Rs >= 1) & (Rs <= Q + 1)]
            e, A0 = exponents_over_R(lam, r, Q, Rs, regime)
            i = int(np.argmin(e))
            if best is None or e[i] < best[0]:
                best = (float(e[i]), dict(r=r, Q=Q, R=int(Rs[i]), A0=int(A0[i]), regime=regime))
    if best is None:
        return dict(r=2, Q=1, R=1, A0=0, regime=regime)
    return best[1]
