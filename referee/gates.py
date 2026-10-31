"""FROZEN referee gates. Usage: python gates.py <algo_dir> <gate: verify|measure|score> <out_json>

Runs the agent's code from <algo_dir> in this separate process; only frozen code decides.
"""
import json
import math
import os
import sys
import time

import numpy as np

ALGO, GATE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
CASES = json.load(open(os.path.join(HERE, "cases.json"), encoding="utf-8"))
sys.path.insert(0, HERE)            # frozen assembly first
import assembly as A                # noqa: E402
sys.path.insert(1, ALGO)
os.chdir(ALGO)

LEMMA_KEYS = ("f_interp", "eigen", "orthogonal", "support", "P0_delta", "degree", "separation")


def brute(n, clauses):
    """Frozen brute-force count (independent of the agent's code)."""
    xs = np.arange(1 << n)
    sat = np.ones(xs.shape, bool)
    for cl in clauses:
        c = np.zeros(xs.shape, bool)
        for v, pos in cl:
            bit = (xs >> v) & 1
            c |= (bit == 1) if pos else (bit == 0)
        sat &= c
    return int(sat.sum())


def verify():
    from construction import Aux, random_formula, random_bounded_formula, build_design, check_design
    from counting import Counter
    out = {"ok": True, "items": []}
    for cfg in CASES["lemma_cfgs"]:
        rep = Aux(*cfg).report()
        ok = all(bool(rep[k]) for k in LEMMA_KEYS) and int(rep["min_kernel_dim"]) > 0
        out["items"].append({"kind": "lemma", "cfg": cfg, "ok": ok})
        out["ok"] &= ok
    for d in CASES["designs"]:
        rng = np.random.default_rng(d["seed"])
        aux = Aux(*d["cfg"])
        cl = random_bounded_formula(d["n"], rng)
        des = build_design(aux, d["n"], cl, rng)
        rep = check_design(aux, des)
        ok = bool(rep["kernel_condition"]) and bool(rep["P0_delta"])
        out["items"].append({"kind": "design", "cfg": d["cfg"], "n": d["n"], "ok": ok,
                             "min_kernel_dim": int(rep["min_kernel_dim"])})
        out["ok"] &= ok
    for c in CASES["counts"]:
        rng = np.random.default_rng(c["seed"])
        aux = Aux(*c["cfg"])
        cl = random_formula(c["n"], c["m"], rng) if c["type"] == "random" else random_bounded_formula(c["n"], rng)
        des = build_design(aux, c["n"], cl, rng)
        C = Counter(aux, des)
        data = {t: C.taylor_inclusion_exclusion(t, Zfun=lambda t_, S, U: C.Z_group_separated(t_, S, U))
                for t in range(1, aux.Q + 1)}
        h0, b = int(C.hermite_H0(data)), brute(c["n"], cl)
        ok = h0 == b
        out["items"].append({"kind": "count", "cfg": c["cfg"], "n": c["n"], "ok": ok, "H0": h0, "brute": b})
        out["ok"] &= ok
    return out


def measure():
    """The accounting facts used by the frozen assembly must hold on real designs."""
    from construction import Aux, random_bounded_formula, build_design, check_design
    out = {"ok": True, "items": []}
    for d in CASES["measure_designs"]:
        r, Q, R, A0 = d["cfg"]
        rng = np.random.default_rng(d["seed"])
        aux = Aux(r, Q, R, A0)
        cl = random_bounded_formula(d["n"], rng)
        des = build_design(aux, d["n"], cl, rng)
        rep = check_design(aux, des)
        M, D, h = int(des["M"]), int(aux.D), int(des["h"])
        Prow = des["Prow"]
        deg = max([k for k in range(Prow.shape[0]) if np.any(Prow[k])] or [0])
        if d["regime"] == "paper":
            in_regime = 2 * A.BW * D * D <= M
            hb = M / D + math.sqrt(2 * A.BW * M)
        else:
            in_regime = D * math.log(D) <= M / (4 * A.BW)
            hb = 2 * M / D
        checks = {
            "kernel_condition": bool(rep["kernel_condition"]),
            "design_in_regime": bool(in_regime),
            "groups_per_row": int(rep["max_groups_per_row"]) <= A.groups_per_row(r),
            "row_abs_sum": int(rep["max_restricted_abs"]) <= A.row_abs_sum_bound(Q),
            "row_degree": deg <= A.row_degree(R),
            "copies_h": h <= hb,
        }
        ok = all(checks.values())
        out["items"].append({"cfg": d["cfg"], "n": d["n"], "regime": d["regime"], "ok": ok, "checks": checks,
                             "M": M, "D": D, "h": h, "h_bound": round(hb, 1),
                             "max_groups_per_row": int(rep["max_groups_per_row"]),
                             "max_restricted_abs": int(rep["max_restricted_abs"]), "degree": deg})
        out["ok"] &= ok
    return out


def score():
    """Runs the agent's search; records every (lam, params, e) query. The parent process
    re-scores every query with a clean copy of the assembly, so this process is untrusted."""
    from model import choose_params
    calls = {"n": 0, "secs": 0.0}
    queries = []

    def f(lam):
        t0 = time.time()
        p = choose_params(int(lam))
        calls["n"] += 1
        calls["secs"] += time.time() - t0
        p = {k: (int(v) if k in ("r", "Q", "R", "A0") else v) for k, v in dict(p).items()}
        s = A.exponent_of(p, int(lam))
        e = s["e"] if s["valid"] else 9.0
        queries.append({"lam": int(lam), "params": p, "e": e})
        return e, p, s

    def first_lam(target, lo=10 ** 3, hi=10 ** 13):
        if f(hi)[0] >= target:
            return None
        if f(lo)[0] < target:
            return lo
        while hi / lo > 1.0005 and hi - lo > 1:
            mid = int(round(math.sqrt(lo * hi)))
            if mid in (lo, hi):
                break
            if f(mid)[0] < target:
                hi = mid
            else:
                lo = mid
        return hi

    res = {}
    for key in ("ppsz", "sharp3sat", "brute"):
        lam = first_lam(A.TARGETS[key])
        res[key] = lam
    e_ref, p_ref, s_ref = f(2224173)
    res["e_at_2224173"] = e_ref
    res["params_at_2224173"] = p_ref
    res["terms_at_2224173"] = s_ref.get("terms")
    res["choose_params_calls"] = calls["n"]
    res["choose_params_avg_s"] = calls["secs"] / max(calls["n"], 1)
    res["queries"] = queries
    res["ok"] = res["ppsz"] is not None and res["choose_params_avg_s"] < 3.0
    return res


if __name__ == "__main__":
    t0 = time.time()
    try:
        result = {"verify": verify, "measure": measure, "score": score}[GATE]()
    except Exception as exc:  # any crash in agent code fails the gate
        import traceback
        result = {"ok": False, "error": f"{type(exc).__name__}: {exc}",
                  "trace": traceback.format_exc()[-2000:]}
    result["seconds"] = round(time.time() - t0, 1)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, default=str)
