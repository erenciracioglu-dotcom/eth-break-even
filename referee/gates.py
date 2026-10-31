"""FROZEN referee gates, version 3. Usage: python gates.py <algo_dir> <verify|measure|score> <out_json>

Runs the agent's code from <algo_dir> in this separate process; only frozen code decides.
Version 2 adds the 'tight' regime and, when construction.build_design accepts designated=True,
frozen checks of the designated-variable construction (exact counts with the reduced Taylor
order, constant designated rows, group schedules, copies).
Version 3 adds, when counting.Counter offers taylor_aggregated and Zagg_group_separated,
exact counts with Lemma 3.3 aggregated per U and the number of restricted sums per point.
"""
import inspect
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


def supports_designated(fn):
    try:
        return "designated" in inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False


def designated_problems(aux, des):
    """Frozen checks of one designated design; returns a list of problems (empty = fine)."""
    if des.get("designated") is not True:
        return ["design is not marked designated=True"]
    occ = des["occ"]
    desig = np.asarray(des.get("desig", []), dtype=bool)
    if desig.shape != (len(occ),):
        return ["des['desig'] must hold one flag per occurrence"]
    Prow, av, df = des["Prow"], np.asarray(des["av"]), np.asarray(des["df"])
    probs, per_factor = [], {}
    for i, (f, v) in enumerate(occ):
        if not desig[i]:
            continue
        per_factor[f] = per_factor.get(f, 0) + 1
        if not np.array_equal(df[f], av[v]):
            probs.append(f"factor {f}: its group schedule differs from that of its designated variable {v}")
        row0 = np.zeros(Prow.shape[2], np.int64)
        row0[v] = aux.Lden
        if not np.array_equal(Prow[0, i], row0) or np.any(Prow[1:, i]):
            probs.append(f"occurrence {i}: the designated row is not the constant Lden * e_v")
    factors = {f for f, _ in occ}
    bad = [f for f in factors if per_factor.get(f, 0) != 1]
    if bad:
        probs.append(f"{len(bad)} factors do not have exactly one designated occurrence")
    if int(des.get("M_eff", -1)) != int((~desig).sum()):
        probs.append("des['M_eff'] is not the number of non-designated occurrences")
    return probs[:5]


def make_formula(c, rng):
    from construction import random_formula, random_bounded_formula
    if c["type"] == "random":
        return random_formula(c["n"], c["m"], rng)
    return random_bounded_formula(c["n"], rng)


def run_count(aux, des, cl, n):
    from counting import Counter
    C = Counter(aux, des)
    data = {t: C.taylor_inclusion_exclusion(t, Zfun=lambda t_, S, U: C.Z_group_separated(t_, S, U))
            for t in range(1, aux.Q + 1)}
    return int(C.hermite_H0(data)), brute(n, cl), int(C.J)


def verify():
    from construction import Aux, random_bounded_formula, build_design, check_design
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
        cl = make_formula(c, rng)
        des = build_design(aux, c["n"], cl, rng)
        h0, b, _ = run_count(aux, des, cl, c["n"])
        ok = h0 == b
        out["items"].append({"kind": "count", "cfg": c["cfg"], "n": c["n"], "ok": ok, "H0": h0, "brute": b})
        out["ok"] &= ok
    out["designated_supported"] = supports_designated(build_design)
    out["designated_verified"] = False
    if out["designated_supported"]:
        dok = True
        for c in CASES["designated_counts"]:
            rng = np.random.default_rng(c["seed"])
            aux = Aux(*c["cfg"])
            cl = make_formula(c, rng)
            des = build_design(aux, c["n"], cl, rng, designated=True)
            probs = designated_problems(aux, des)
            rep = check_design(aux, des)
            if not (rep["kernel_condition"] and rep["P0_delta"]):
                probs.append("kernel condition or P(0) = delta failed")
            h0, b, J = run_count(aux, des, cl, c["n"])
            j_exp = sum(max(len(x) - 1, 0) for x in cl) * (aux.R - 1) // aux.Q + 1
            if J != j_exp:
                probs.append(f"Taylor order J = {J}, expected {j_exp} from the designated degree")
            ok = (h0 == b) and not probs
            out["items"].append({"kind": "designated_count", "cfg": c["cfg"], "n": c["n"], "ok": ok,
                                 "H0": h0, "brute": b, "J": J, "J_expected": j_exp, "problems": probs})
            dok &= ok
        out["designated_verified"] = bool(dok)
        out["ok"] &= dok
    from counting import Counter
    out["aggregated_supported"] = all(hasattr(Counter, a) for a in ("taylor_aggregated", "Zagg_group_separated"))
    out["aggregated_verified"] = False
    if out["aggregated_supported"]:
        aok, ran = True, 0
        for c in CASES["aggregated_counts"]:
            if c.get("designated") and not out["designated_supported"]:
                continue
            ran += 1
            rng = np.random.default_rng(c["seed"])
            aux = Aux(*c["cfg"])
            cl = make_formula(c, rng)
            kw = {"designated": True} if c.get("designated") else {}
            des = build_design(aux, c["n"], cl, rng, **kw)
            C = Counter(aux, des)
            data, probs = {}, []
            for t in range(1, aux.Q + 1):
                data[t] = C.taylor_aggregated(t, Zagg=lambda t_, S, U: C.Zagg_group_separated(t_, S, U))
                exp_rs = sum(math.comb(c["n"], s) for s in range(min(int(C.q), c["n"]) + 1))
                got = getattr(C, "restricted_sums", None)
                if got != exp_rs:
                    probs.append(f"t={t}: restricted_sums = {got}, expected {exp_rs}")
            h0, b = int(C.hermite_H0(data)), brute(c["n"], cl)
            ok = (h0 == b) and not probs
            out["items"].append({"kind": "aggregated_count", "cfg": c["cfg"], "n": c["n"],
                                 "designated": bool(c.get("designated")), "ok": ok, "H0": h0, "brute": b,
                                 "problems": probs[:3]})
            aok &= ok
        out["aggregated_verified"] = bool(aok and ran > 0)
        out["ok"] &= aok
    return out


def measure():
    """The accounting facts used by the frozen assembly must hold on real designs."""
    from construction import Aux, random_bounded_formula, build_design, check_design
    out = {"ok": True, "items": []}

    def common_checks(rep, Prow, r, Q, R):
        deg = max([k for k in range(Prow.shape[0]) if np.any(Prow[k])] or [0])
        return deg, {
            "kernel_condition": bool(rep["kernel_condition"]),
            "groups_per_row": int(rep["max_groups_per_row"]) <= A.groups_per_row(r),
            "row_abs_sum": int(rep["max_restricted_abs"]) <= A.row_abs_sum_bound(Q),
            "row_degree": deg <= A.row_degree(R),
        }

    for d in CASES["measure_designs"]:
        r, Q, R, A0 = d["cfg"]
        rng = np.random.default_rng(d["seed"])
        aux = Aux(r, Q, R, A0)
        cl = random_bounded_formula(d["n"], rng)
        des = build_design(aux, d["n"], cl, rng)
        rep = check_design(aux, des)
        M, D, h = int(des["M"]), int(aux.D), int(des["h"])
        if d["regime"] == "paper":
            in_regime = 2 * A.BW * D * D <= M
            hb = M / D + math.sqrt(2 * A.BW * M)
        else:
            dl = float(d.get("delta", 1.0)) if d["regime"] == "tight" else 1.0
            in_regime = D * math.log(D) <= dl * dl * M / (4 * A.BW)
            hb = (1 + dl) * M / D
        deg, checks = common_checks(rep, des["Prow"], r, Q, R)
        checks.update(design_in_regime=bool(in_regime), copies_h=h <= hb)
        ok = all(checks.values())
        out["items"].append({"cfg": d["cfg"], "n": d["n"], "regime": d["regime"], "ok": ok, "checks": checks,
                             "M": M, "D": D, "h": h, "h_bound": round(hb, 1),
                             "max_groups_per_row": int(rep["max_groups_per_row"]),
                             "max_restricted_abs": int(rep["max_restricted_abs"]), "degree": deg})
        out["ok"] &= ok
    if supports_designated(build_design):
        for d in CASES["designated_measure"]:
            r, Q, R, A0 = d["cfg"]
            rng = np.random.default_rng(d["seed"])
            aux = Aux(r, Q, R, A0)
            cl = random_bounded_formula(d["n"], rng)
            des = build_design(aux, d["n"], cl, rng, designated=True)
            probs = designated_problems(aux, des)
            rep = check_design(aux, des)
            m_eff, D, h = int(des.get("M_eff", 0)), int(aux.D), int(des["h"])
            hb = 2 * m_eff / D
            deg, checks = common_checks(rep, des["Prow"], r, Q, R)
            checks.update(designated_rows=not probs, copies_h=h <= hb,
                          m_eff_bound=m_eff <= A.M_EFF_DESIG * d["n"])
            ok = all(checks.values())
            out["items"].append({"cfg": d["cfg"], "n": d["n"], "regime": "designated", "ok": ok,
                                 "checks": checks, "M_eff": m_eff, "D": D, "h": h, "h_bound": round(hb, 1),
                                 "problems": probs, "degree": deg})
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
