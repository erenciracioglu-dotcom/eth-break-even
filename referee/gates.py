"""FROZEN referee gates, version 8. Usage: python gates.py <algo_dir> <verify|measure|score> <out_json>

Runs the agent's code from <algo_dir> in this separate process; only frozen code decides.
Version 2 adds the 'tight' regime and, when construction.build_design accepts designated=True,
frozen checks of the designated-variable construction (exact counts with the reduced Taylor
order, constant designated rows, group schedules, copies).
Version 3 adds, when counting.Counter offers taylor_aggregated and Zagg_group_separated,
exact counts with Lemma 3.3 aggregated per U and the number of restricted sums per point.
Version 4 adds, when build_design accepts units=True, checks of clause units: units are the
variable sets of clauses (a maximal matching), share one schedule, and give exact counts.
Version 5 adds, when build_design accepts chunk=L, checks of chunk units (connected units of
>= L variables, small components whole, (1 + 1/L) n bound on real and adversarial formulas),
and, when Counter accepts residue_keys and weights, weighted counts mod a small prime with
measurement fibres keyed by residues; the score gate scores yes/no targets with
model.choose_params_decision when the model offers it.
Version 6 adds, when construction.Aux accepts perpoint and eadapt, frozen checks of per-point
measurement matrices (Aux.Cmeas_t, design key Chat_t) and of equal-pair windows: separation at every
point with the measurements used there (auxiliary rows and real designs), row counts against the
frozen formula, two groups per row, restricted sums, P(0) = delta, degree, and exact (or weighted,
mod p, residue-keyed) counts.
Version 7 adds, when counting.Counter offers Zagg_contracted, exact (and weighted) counts with the groups
summed one at a time, and a frozen check that clauses of group g see group h only through the rows of the
pair {g, h} (the bond condition behind assembly 'contract').
Version 8 adds, when construction.gf_rows exists, frozen checks of replacement rows over GF(p^2) with a
prime p < Q (assembly 'smallprime'): nodes, P(0) = delta, degree, Definition 3.1(4) at every node mod p,
and the weighted count from Taylor data at the nodes (direct evaluation, Hermite interpolation).
"""
import inspect
import math as _math_unused
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


def supports_param(fn, name):
    try:
        return name in inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False


def supports_designated(fn):
    return supports_param(fn, "designated")


def unit_problems(des, clauses):
    """Units must be variable sets of clauses (a maximal matching) sharing one schedule."""
    if des.get("units") is not True:
        return ["design is not marked units=True"]
    av = np.asarray(des["av"])
    n = av.shape[0]
    unit_of = np.asarray(des.get("unit_of", []))
    if unit_of.shape != (n,):
        return ["des['unit_of'] must hold one unit id per variable"]
    members = {}
    for v in range(n):
        members.setdefault(int(unit_of[v]), []).append(v)
    clause_sets = [frozenset(v for v, _ in cl) for cl in clauses]
    known = set(clause_sets)
    probs = []
    for u, vs in members.items():
        if len(vs) > 1:
            if len(vs) > 3 or frozenset(vs) not in known:
                probs.append(f"unit {u} = {vs} is not the variable set of one clause")
            if any(not np.array_equal(av[vs[0]], av[v]) for v in vs):
                probs.append(f"unit {u}: its variables do not share one schedule")
    single = {vs[0] for vs in members.values() if len(vs) == 1}
    free = [i for i, s in enumerate(clause_sets) if len(s) >= 2 and all(v in single for v in s)]
    if free:
        probs.append(f"the matching of clauses is not maximal: clause {free[0]} could still be matched")
    return probs


def var_graph(n, clauses):
    adj = [set() for _ in range(n)]
    for cl in clauses:
        vs = sorted({int(v) for v, _ in cl})
        for x in vs:
            adj[x].update(y for y in vs if y != x)
    return adj


def components(n, clauses):
    adj = var_graph(n, clauses)
    comp = -np.ones(n, np.int64)
    for s in range(n):
        if comp[s] >= 0:
            continue
        comp[s] = s
        stack = [s]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if comp[w] < 0:
                    comp[w] = s
                    stack.append(w)
    return comp, adj


def chunk_problems(des, clauses, L):
    """Chunk units: small components are one unit; in larger ones every unit has L to
    1 + (deg + 1)(L - 1) variables and is connected through clauses with >= 2 of its variables."""
    if des.get("units") is not True:
        return ["design is not marked units=True"]
    if int(des.get("chunk", 0)) != int(L):
        return [f"design is not marked chunk={L}"]
    av = np.asarray(des["av"])
    n = av.shape[0]
    unit_of = np.asarray(des.get("unit_of", []))
    if unit_of.shape != (n,):
        return ["des['unit_of'] must hold one unit id per variable"]
    comp, adj = components(n, clauses)
    deg = max((len(x) for x in adj), default=0)
    cap = 1 + (deg + 1) * (L - 1)
    csize = {}
    for v in range(n):
        csize[int(comp[v])] = csize.get(int(comp[v]), 0) + 1
    members = {}
    for v in range(n):
        members.setdefault(int(unit_of[v]), []).append(v)
    nbr = [set() for _ in range(n)]
    for cl in clauses:
        byu = {}
        for v in {int(v) for v, _ in cl}:
            byu.setdefault(int(unit_of[v]), []).append(v)
        for vs in byu.values():
            if len(vs) >= 2:
                for x in vs:
                    nbr[x].update(y for y in vs if y != x)
    probs = []
    for u, vs in members.items():
        cs = {int(comp[v]) for v in vs}
        if len(cs) != 1:
            probs.append(f"unit {u} spans {len(cs)} components")
            continue
        size_c = csize[cs.pop()]
        if any(not np.array_equal(av[vs[0]], av[v]) for v in vs):
            probs.append(f"unit {u}: its variables do not share one schedule")
        if size_c < L:
            if len(vs) != size_c:
                probs.append(f"unit {u}: a component of {size_c} < L variables must be a single unit")
        elif not (L <= len(vs) <= cap):
            probs.append(f"unit {u} has {len(vs)} variables, outside [{L}, {cap}]")
        if len(vs) > 1:
            seen, stack = {vs[0]}, [vs[0]]
            while stack:
                x = stack.pop()
                for y in nbr[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            if len(seen) != len(vs):
                probs.append(f"unit {u} is not connected through its clauses")
        if len(probs) >= 5:
            break
    return probs[:5]


def designation_problems(des):
    """Per factor: the designated occurrences are all occurrences in one unit with the most of them."""
    occ = des["occ"]
    desig = np.asarray(des.get("desig", []), dtype=bool)
    unit_of = np.asarray(des.get("unit_of", []))
    per_f = {}
    for i, (f, _v) in enumerate(occ):
        per_f.setdefault(f, []).append(i)
    probs = []
    for f, idx in per_f.items():
        cnt = {}
        for i in idx:
            u = int(unit_of[occ[i][1]])
            cnt[u] = cnt.get(u, 0) + 1
        du = {int(unit_of[occ[i][1]]) for i in idx if desig[i]}
        if len(du) != 1:
            probs.append(f"factor {f}: designated occurrences lie in {len(du)} units")
            continue
        u = du.pop()
        if cnt[u] != max(cnt.values()):
            probs.append(f"factor {f}: its designated unit does not hold the most of its occurrences")
        if sum(1 for i in idx if desig[i]) != cnt[u]:
            probs.append(f"factor {f}: not every occurrence in its designated unit is designated")
        if len(probs) >= 5:
            break
    return probs


def fano_formula(k, rng):
    """k disjoint Fano planes (7 variables, 7 clauses of width 3, every variable in 3 clauses)."""
    lines = [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]
    perm = rng.permutation(7 * k)
    clauses = []
    for b in range(k):
        for ln in lines:
            clauses.append([(int(perm[7 * b + x]), bool(rng.random() < 0.5)) for x in ln])
    return clauses


def weighted_brute(n, clauses, w, p):
    """Frozen sum over satisfying x of prod_v w_v^(x_v) mod p."""
    xs = np.arange(1 << n)
    sat = np.ones(xs.shape, bool)
    for cl in clauses:
        c = np.zeros(xs.shape, bool)
        for v, pos in cl:
            bit = (xs >> v) & 1
            c |= (bit == 1) if pos else (bit == 0)
        sat &= c
    tot = 0
    for x in xs[sat]:
        prod = 1
        for v in range(n):
            if (int(x) >> v) & 1:
                prod = prod * int(w[v]) % p
        tot = (tot + prod) % p
    return tot


def fibres_frozen(aux, des, t, p=None):
    """Number of distinct measurement keys per group at point t (exact, or residues mod p)."""
    av, Chat = np.asarray(des["av"]), np.asarray(des["Chat"])
    out = []
    for g in range(aux.r):
        vars_g = np.nonzero(av[:, t - 1] == g)[0]
        k = len(vars_g)
        if k == 0:
            out.append(1)
            continue
        lm = np.arange(1 << k)
        bits = ((lm[:, None] >> np.arange(k)[None, :]) & 1).astype(np.int64)
        vals = bits @ Chat[:, vars_g].T
        if p is not None:
            vals = vals % p
        out.append(len({tuple(int(x) for x in row) for row in vals}))
    return out


def designated_problems(aux, des, units=False, clauses=None, chunk=0):
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
    if units:
        bad = [f for f in factors if per_factor.get(f, 0) < 1]
        if bad:
            probs.append(f"{len(bad)} factors have no designated occurrence")
        probs += chunk_problems(des, clauses, chunk) if chunk else unit_problems(des, clauses)
    else:
        bad = [f for f in factors if per_factor.get(f, 0) != 1]
        if bad:
            probs.append(f"{len(bad)} factors do not have exactly one designated occurrence")
    if int(des.get("M_eff", -1)) != int((~desig).sum()):
        probs.append("des['M_eff'] is not the number of non-designated occurrences")
    return probs[:5]


# ------------------------------------------------ v6: per-point measurements, equal-pair windows
FROZEN_P = 33554393                     # prime below 2^25; residue products fit in int64


def _rref_mod(A, p):
    A = np.array(A % p, dtype=np.int64)
    nr, nc = A.shape
    piv, r = [], 0
    for c in range(nc):
        if r >= nr:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        piv.append(c)
        r += 1
    return A[:r], piv


def _in_rowspace_mod(C, B, p):
    """Every row of B in the row space of C (mod p)?"""
    Bm = np.array(B % p, dtype=np.int64)
    if C.shape[0] == 0:
        return not np.any(Bm)
    Rr, piv = _rref_mod(C, p)
    if not piv:
        return not np.any(Bm)
    acc = np.zeros(Bm.shape, np.int64)
    for k, c in enumerate(piv):
        acc = (acc + (Bm[:, c:c + 1] * Rr[k][None, :]) % p) % p
    return not np.any((Bm - acc) % p)


def meas_at(des, t):
    """The measurement matrix (rows over the variables) that the design uses at point t."""
    if des.get("perpoint"):
        return np.asarray(des["Chat_t"][t - 1], np.int64)
    return np.asarray(des["Chat"], np.int64)


def aux_pe_problems(aux, pp, ea):
    """Frozen checks of an auxiliary construction built with perpoint / eadapt."""
    r, Q, R, A0, D = aux.r, aux.Q, aux.R, aux.A0, aux.D
    a, d = np.asarray(aux.a), np.asarray(aux.d)
    Pc = np.asarray(aux.Pc, np.int64)
    probs = []
    if bool(getattr(aux, "perpoint", False)) != pp or bool(getattr(aux, "eadapt", False)) != ea:
        probs.append("Aux does not record the requested perpoint / eadapt flags")
    if Pc.shape[0] > R and np.any(Pc[R:]):
        probs.append("a replacement row has degree above R - 1")
    if not np.array_equal(Pc[0], aux.Lden * np.eye(D, dtype=np.int64)):
        probs.append("P(0) is not Lden * I")
    bound = D * A.tail_out(Q, r, R, A0, perpoint=pp, eadapt=ea) + 1e-6
    p = FROZEN_P
    PcM = Pc % p
    for t in range(1, Q + 1):
        Cm = np.asarray(aux.Cmeas_t[t - 1] if pp else aux.Cmeas, np.int64)
        if Cm.shape[0] > bound:
            probs.append(f"t={t}: {Cm.shape[0]} auxiliary measurement rows, frozen bound {bound:.2f}")
        grp = a[:, t - 1]
        touched = np.zeros(Cm.shape[0], np.int64)
        for g in range(r):
            sub = Cm[:, grp == g]
            touched += np.any(sub != 0, axis=1)
            if sub.size and int(np.abs(sub).sum(axis=1).max()) > A.row_abs_sum_bound(Q):
                probs.append(f"t={t}: a restricted measurement row has absolute sum above 2^(Q-1)")
        if touched.size and int(touched.max()) > A.groups_per_row(r):
            probs.append(f"t={t}: a measurement row touches {int(touched.max())} groups")
        Pt = np.zeros((D, D), np.int64)
        for k in range(Pc.shape[0]):
            Pt = (Pt + (pow(t, k, p) * PcM[k]) % p) % p
        for g in range(r):
            cols = np.nonzero(grp == g)[0]
            rows = np.nonzero(d[:, t - 1] != g)[0]
            if cols.size and rows.size and not _in_rowspace_mod(Cm[:, cols], Pt[np.ix_(rows, cols)], p):
                probs.append(f"t={t}, group {g}: replacement rows outside the span of the measurements")
                break
        if len(probs) >= 5:
            break
    return probs[:5]


def pe_design_problems(aux, des, pp, ea):
    """Frozen checks of a real design: Definition 3.1(4) at every point with the measurement matrix
    used there, row counts (h copies of the frozen bound), two groups per row, restricted sums."""
    probs = []
    if bool(des.get("perpoint", False)) != pp:
        probs.append("design does not record perpoint as requested")
        return probs
    Q, r = aux.Q, aux.r
    Prow = np.asarray(des["Prow"], np.int64)
    av, df, occ = np.asarray(des["av"]), np.asarray(des["df"]), des["occ"]
    h = int(des["h"])
    p = FROZEN_P
    PrM = Prow % p
    bound = aux.D * A.tail_out(Q, r, aux.R, aux.A0, perpoint=pp, eadapt=ea) + 1e-6
    for t in range(1, Q + 1):
        Cm = meas_at(des, t)
        if Cm.shape[0] > h * bound:
            probs.append(f"t={t}: {Cm.shape[0]} measurement rows, frozen bound h * {bound:.2f}")
        Pt = np.zeros(PrM.shape[1:], np.int64)
        for k in range(PrM.shape[0]):
            Pt = (Pt + (pow(t, k, p) * PrM[k]) % p) % p
        touched = np.zeros(Cm.shape[0], np.int64)
        for g in range(r):
            cols = np.nonzero(av[:, t - 1] == g)[0]
            sub = Cm[:, cols]
            touched += np.any(sub != 0, axis=1)
            if sub.size and int(np.abs(sub).sum(axis=1).max()) > A.row_abs_sum_bound(Q):
                probs.append(f"t={t}: a restricted measurement row has absolute sum above 2^(Q-1)")
            rows = [i for i, (f, _v) in enumerate(occ) if df[f][t - 1] != g]
            if cols.size and rows and not _in_rowspace_mod(sub, Pt[np.ix_(rows, cols)], p):
                probs.append(f"t={t}, group {g}: kernel condition (Definition 3.1(4)) fails")
        if touched.size and int(touched.max()) > A.groups_per_row(r):
            probs.append(f"t={t}: a measurement row touches {int(touched.max())} groups")
        if len(probs) >= 5:
            break
    return probs[:5]


def fibres_at(aux, des, t, p=None):
    """Distinct measurement keys per group at point t with the matrix used there."""
    av, Cm = np.asarray(des["av"]), meas_at(des, t)
    out = []
    for g in range(aux.r):
        vars_g = np.nonzero(av[:, t - 1] == g)[0]
        k = len(vars_g)
        if k == 0:
            out.append(1)
            continue
        lm = np.arange(1 << k)
        bits = ((lm[:, None] >> np.arange(k)[None, :]) & 1).astype(np.int64)
        vals = bits @ Cm[:, vars_g].T
        if p is not None:
            vals = vals % p
        out.append(len({tuple(int(x) for x in row) for row in vals}))
    return out


def flag_supported(Aux, name):
    """A construction flag counts as supported when Aux accepts it and records it (aux.<name> True)."""
    if not supports_param(Aux.__init__, name):
        return False
    try:
        return bool(getattr(Aux(2, 2, 2, 0, **{name: True}), name, False))
    except Exception:
        return False


def pe_kw(pp, ea):
    kw = {}
    if pp:
        kw["perpoint"] = True
    if ea:
        kw["eadapt"] = True
    return kw


def pe_verify(out, Aux, build_design, Counter):
    """v6 section of the verify gate."""
    out["perpoint_supported"] = flag_supported(Aux, "perpoint")
    out["eadapt_supported"] = flag_supported(Aux, "eadapt")
    out["perpoint_verified"] = False
    out["eadapt_verified"] = False
    if not (out["perpoint_supported"] or out["eadapt_supported"]):
        return
    okf = {"pp": True, "ea": True}
    ran = {"pp": 0, "ea": 0}

    def note(pp, ea, ok):
        for key, flag in (("pp", pp), ("ea", ea)):
            if flag:
                okf[key] &= bool(ok)
                ran[key] += 1

    for c in CASES["pe_lemma"]:
        pp, ea = bool(c.get("perpoint")), bool(c.get("eadapt"))
        if (pp and not out["perpoint_supported"]) or (ea and not out["eadapt_supported"]):
            continue
        aux = Aux(*c["cfg"], **pe_kw(pp, ea))
        rep = aux.report()
        probs = [k for k in LEMMA_KEYS if not bool(rep.get(k))]
        probs = ([f"report fails: {', '.join(probs)}"] if probs else []) + aux_pe_problems(aux, pp, ea)
        ok = not probs
        rows = [int(np.asarray(aux.Cmeas_t[t]).shape[0]) for t in range(aux.Q)] if pp else [int(np.asarray(aux.Cmeas).shape[0])]
        out["items"].append({"kind": "pe_lemma", "cfg": c["cfg"], "perpoint": pp, "eadapt": ea, "ok": ok,
                             "rows": rows, "bound": round(aux.D * A.tail_out(aux.Q, aux.r, aux.R, aux.A0, perpoint=pp, eadapt=ea), 2),
                             "problems": probs})
        note(pp, ea, ok)
    for c in CASES["pe_counts"]:
        pp, ea = bool(c.get("perpoint")), bool(c.get("eadapt"))
        L, desg, p = int(c.get("L", 0)), bool(c.get("designated")), c.get("p")
        if (pp and not out["perpoint_supported"]) or (ea and not out["eadapt_supported"]):
            continue
        if (desg and not out["designated_supported"]) or (L and not out["chunk_supported"]) \
                or (p and not out["modular_supported"]):
            continue
        rng = np.random.default_rng(c["seed"])
        aux = Aux(*c["cfg"], **pe_kw(pp, ea))
        cl = make_formula(c, rng)
        kw = {"designated": True} if desg else {}
        if L:
            kw.update(units=True, chunk=L)
        des = build_design(aux, c["n"], cl, rng, **kw)
        probs = pe_design_problems(aux, des, pp, ea)
        if p:
            p = int(p)
            extra = p * np.ones((1, c["n"]), np.int64)     # redundant row, 0 mod p: residue keys must merge
            if pp:
                des = dict(des, Chat_t=[np.vstack([np.asarray(M, np.int64), extra]) for M in des["Chat_t"]])
            else:
                des = dict(des, Chat=np.vstack([np.asarray(des["Chat"], np.int64), extra]))
            wrng = np.random.default_rng(c["seed"] + 1)
            for _ in range(16):
                w = wrng.integers(1, p, size=c["n"])
                if weighted_brute(c["n"], cl, w, p):
                    break
            C = Counter(aux, des, p=p, residue_keys=True, weights=w)
            data = {}
            for t in range(1, aux.Q + 1):
                data[t] = C.taylor_aggregated(t, Zagg=lambda t_, S, U: C.Zagg_group_separated(t_, S, U))
                got = [int(x) for x in C.fibre_counts(t)]
                exp = fibres_at(aux, des, t, p)
                if got != exp:
                    probs.append(f"t={t}: fibres per group {got}, expected {exp} (keys mod {p})")
            h0, ref = int(C.hermite_H0(data)) % p, weighted_brute(c["n"], cl, w, p)
        else:
            C = Counter(aux, des)
            data = {}
            for t in range(1, aux.Q + 1):
                if out.get("aggregated_supported"):
                    data[t] = C.taylor_aggregated(t, Zagg=lambda t_, S, U: C.Zagg_group_separated(t_, S, U))
                else:
                    data[t] = C.taylor_inclusion_exclusion(t, Zfun=lambda t_, S, U: C.Z_group_separated(t_, S, U))
                exp = fibres_at(aux, des, t, None)
                got = [int(x) for x in C.fibre_counts(t)] if hasattr(C, "fibre_counts") else exp
                if got != exp:
                    probs.append(f"t={t}: fibres per group {got}, expected {exp}")
            h0, ref = int(C.hermite_H0(data)), brute(c["n"], cl)
        ok = (h0 == ref) and not probs
        out["items"].append({"kind": "pe_count", "cfg": c["cfg"], "n": c["n"], "perpoint": pp, "eadapt": ea,
                             "L": L, "p": p, "ok": ok, "H0": h0, "reference": ref, "problems": probs[:5]})
        note(pp, ea, ok)
    out["perpoint_verified"] = bool(out["perpoint_supported"] and okf["pp"] and ran["pp"] > 0)
    out["eadapt_verified"] = bool(out["eadapt_supported"] and okf["ea"] and ran["ea"] > 0)
    out["ok"] &= okf["pp"] and okf["ea"]


# ------------------------------------------------ v7: contraction of the group sums
def bond_problems(aux, des):
    """Frozen check behind 'contract': at every point every measurement row touches at most 2 groups,
    and for every clause of group g and every other group h, its replacement rows restricted to S_h
    lie in the span of the rows of the pair {g, h} restricted to S_h (mod a prime)."""
    Q, r = aux.Q, aux.r
    Prow = np.asarray(des["Prow"], np.int64)
    av, df, occ = np.asarray(des["av"]), np.asarray(des["df"]), des["occ"]
    p = FROZEN_P
    PrM = Prow % p
    probs = []
    for t in range(1, Q + 1):
        Cm = meas_at(des, t)
        cols_of = [np.nonzero(av[:, t - 1] == g)[0] for g in range(r)]
        touch = np.zeros((r, Cm.shape[0]), bool)
        for g in range(r):
            if cols_of[g].size:
                touch[g] = np.any(Cm[:, cols_of[g]] != 0, axis=1)
        if Cm.shape[0] and int(touch.sum(axis=0).max()) > 2:
            probs.append(f"t={t}: a measurement row touches more than 2 groups")
        Pt = np.zeros(PrM.shape[1:], np.int64)
        for k in range(PrM.shape[0]):
            Pt = (Pt + (pow(t, k, p) * PrM[k]) % p) % p
        for g in range(r):
            rows = [i for i, (f, _v) in enumerate(occ) if df[f][t - 1] == g]
            if not rows:
                continue
            for h in range(r):
                if h == g or not cols_of[h].size:
                    continue
                pair = touch[h] & (touch.sum(axis=0) - touch[g] - touch[h] == 0)
                sub = Cm[pair][:, cols_of[h]]
                if not _in_rowspace_mod(sub, Pt[np.ix_(rows, cols_of[h])], p):
                    probs.append(f"t={t}: clauses of group {g} see group {h} outside the rows of the pair")
                    break
            if len(probs) >= 5:
                break
        if len(probs) >= 5:
            break
    return probs[:5]


def contract_verify(out, Aux, build_design, Counter):
    """v7 section of the verify gate: exact (or weighted) counts with Counter.Zagg_contracted."""
    out["contract_supported"] = bool(out.get("aggregated_supported")) and hasattr(Counter, "Zagg_contracted")
    out["contract_verified"] = False
    if not out["contract_supported"]:
        return
    ok_all, ran = True, 0
    for c in CASES["contract_counts"]:
        pp, ea = bool(c.get("perpoint")), bool(c.get("eadapt"))
        L, desg, p = int(c.get("L", 0)), bool(c.get("designated")), c.get("p")
        if (pp and not out.get("perpoint_supported")) or (ea and not out.get("eadapt_supported")):
            continue
        if (desg and not out["designated_supported"]) or (L and not out["chunk_supported"]) \
                or (p and not out["modular_supported"]):
            continue
        ran += 1
        rng = np.random.default_rng(c["seed"])
        kw_aux = {}
        if pp:
            kw_aux["perpoint"] = True
        if ea:
            kw_aux["eadapt"] = True
        aux = Aux(*c["cfg"], **kw_aux)
        cl = make_formula(c, rng)
        kw = {"designated": True} if desg else {}
        if L:
            kw.update(units=True, chunk=L)
        des = build_design(aux, c["n"], cl, rng, **kw)
        probs = bond_problems(aux, des)
        if p:
            p = int(p)
            wrng = np.random.default_rng(c["seed"] + 1)
            for _ in range(16):
                w = wrng.integers(1, p, size=c["n"])
                if weighted_brute(c["n"], cl, w, p):
                    break
            C = Counter(aux, des, p=p, residue_keys=True, weights=w)
            ref = weighted_brute(c["n"], cl, w, p)
        else:
            C = Counter(aux, des)
            ref = brute(c["n"], cl)
        data = {t: C.taylor_aggregated(t, Zagg=lambda t_, S, U: C.Zagg_contracted(t_, S, U))
                for t in range(1, aux.Q + 1)}
        h0 = int(C.hermite_H0(data)) % (p if p else (1 << 62))
        ok = (h0 == ref) and not probs
        out["items"].append({"kind": "contract_count", "cfg": c["cfg"], "n": c["n"], "perpoint": pp, "eadapt": ea,
                             "L": L, "p": p, "ok": ok, "H0": h0, "reference": ref, "problems": probs[:5]})
        ok_all &= ok
    out["contract_verified"] = bool(ok_all and ran > 0)
    out["ok"] &= ok_all


# ------------------------------------------------ v8: small prime with nodes in GF(p^2)
def gf_c(p):
    """Least quadratic non-residue mod an odd prime p (GF(p^2) = GF(p)[s] / (s^2 - c))."""
    return next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)


def gf_nodes(p, Q):
    """The first Q nonzero elements a + b s of GF(p^2), ordered by (b, a): 1..p-1 first."""
    out = []
    for b in range(p):
        for a in range(p):
            if (a, b) != (0, 0):
                out.append((a, b))
            if len(out) == Q:
                return np.array(out, np.int64)
    raise ValueError("GF(p^2) has fewer than Q nonzero elements")


def gmul(x, y, p, c):
    """Elementwise product in GF(p^2); x, y int64 arrays with a last axis (a, b)."""
    a = (x[..., 0] * y[..., 0] + c * ((x[..., 1] * y[..., 1]) % p)) % p
    b = (x[..., 0] * y[..., 1] + x[..., 1] * y[..., 0]) % p
    return np.stack([a, b], axis=-1)


def ginv(x, p, c):
    a, b = int(x[0]) % p, int(x[1]) % p
    nrm = (a * a - c * b * b) % p
    if nrm == 0:
        raise ZeroDivisionError("inverse of 0 in GF(p^2)")
    ni = pow(nrm, p - 2, p)
    return np.array([a * ni % p, (-b) * ni % p], np.int64)


def gpow(x, e, p, c):
    out = np.array([1, 0], np.int64)
    base = np.array(x, np.int64) % p
    while e:
        if e & 1:
            out = gmul(out, base, p, c)
        base = gmul(base, base, p, c)
        e >>= 1
    return out


def gsmul(A, B, p, c):
    """Truncated product of series over GF(p^2): A, B arrays (J, ..., 2), coefficient i at A[i]."""
    J = A.shape[0]
    out = np.zeros(np.broadcast_shapes(A.shape, B.shape), np.int64)
    for i in range(J):
        if not np.any(A[i]):
            continue
        out[i:] = (out[i:] + gmul(A[i][None, ...], B[:J - i], p, c)) % p
    return out


def gsolve(M, b, p, c):
    """Solve M x = b over GF(p^2) (M: (L, L, 2), b: (L, 2)); returns x (L, 2)."""
    L = M.shape[0]
    A = np.concatenate([M % p, (b % p)[:, None, :]], axis=1).astype(np.int64)
    for col in range(L):
        piv = next((r for r in range(col, L) if np.any(A[r, col])), None)
        if piv is None:
            raise ValueError("singular Hermite system")
        if piv != col:
            A[[col, piv]] = A[[piv, col]]
        A[col] = gmul(A[col], ginv(A[col, col], p, c)[None, :], p, c)
        for r in range(L):
            if r != col and np.any(A[r, col]):
                f = A[r, col].copy()
                A[r] = (A[r] - gmul(A[col], f[None, :], p, c)) % p
    return A[:, L]


def gf_direct_H0(aux, des, rows, w, p):
    """Weighted count S(w) from Taylor data of H_w(alpha_t + z) mod z^J at the GF(p^2) nodes, each
    computed directly over all assignments, then Hermite interpolation (frozen)."""
    c, nodes, Pg = rows["c"], np.asarray(rows["nodes"], np.int64), np.asarray(rows["Prow_gf"], np.int64) % p
    n, Q, R = des["n"], aux.Q, aux.R
    clauses, occ = des["clauses"], des["occ"]
    desig = np.asarray(des.get("desig", np.zeros(len(occ), bool)), bool)
    d = sum(1 for i, (f, _v) in enumerate(occ) if f < len(clauses) and not desig[i])
    J = d * (R - 1) // Q + 1
    occ_index = {fv: i for i, fv in enumerate(occ)}
    N = 1 << n
    X = ((np.arange(N)[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
    wX = np.ones(N, np.int64)
    for v in range(n):
        wX = np.where(X[:, v] == 1, wX * int(w[v]) % p, wX)
    K = Pg.shape[0]
    forms = np.stack([np.stack([(X @ Pg[k, :, :, comp].T) % p for comp in (0, 1)], axis=-1) for k in range(K)])
    data = {}
    for t in range(1, Q + 1):
        al = nodes[t - 1]
        S = np.zeros((J,) + forms.shape[1:], np.int64)      # (J, N, M, 2)
        for i in range(J):
            for k in range(i, K):
                coef = gpow(al, k - i, p, c) * (math.comb(k, i) % p) % p
                S[i] = (S[i] + gmul(forms[k], coef[None, None, :], p, c)) % p
        total = np.zeros((J, N, 2), np.int64)
        total[0, :, 0] = 1
        for f, cl in enumerate(clauses):
            prod = np.zeros((J, N, 2), np.int64)
            prod[0, :, 0] = 1
            for v, pos in cl:
                s = S[:, :, occ_index[(f, v)], :]
                if pos:
                    om = (-s) % p
                    om[0, :, 0] = (1 - s[0, :, 0]) % p
                else:
                    om = s
                prod = gsmul(prod, om, p, c)
            phi = (-prod) % p
            phi[0, :, 0] = (1 - prod[0, :, 0]) % p
            total = gsmul(total, phi, p, c)
        data[t] = np.stack([(total[:, :, comp] * wX[None, :]).sum(axis=1) % p for comp in (0, 1)], axis=-1)
    L = Q * J
    Mx = np.zeros((L, L, 2), np.int64)
    bx = np.zeros((L, 2), np.int64)
    row = 0
    for t in range(1, Q + 1):
        al = nodes[t - 1]
        for i in range(J):
            for u in range(i, L):
                Mx[row, u] = gpow(al, u - i, p, c) * (math.comb(u, i) % p) % p
            bx[row] = data[t][i]
            row += 1
    sol = gsolve(Mx, bx, p, c)
    return sol[0]


def gf_problems(aux, des, rows, p):
    """Frozen checks of GF(p^2) rows: field and nodes, P(0) = delta, degree, and Definition 3.1(4) at
    every node with the integer measurement rows reduced mod p (both coordinates of a GF(p^2) row must
    lie in the GF(p)-span, which is the GF(p^2)-span of rows over GF(p))."""
    probs = []
    c = int(rows.get("c", -1))
    if c != gf_c(p):
        probs.append(f"rows['c'] = {c} is not the least non-residue mod {p}")
        return probs
    nodes = np.asarray(rows.get("nodes"), np.int64)
    if nodes.shape != (aux.Q, 2) or not np.array_equal(nodes, gf_nodes(p, aux.Q)):
        probs.append("rows['nodes'] are not the first Q nonzero elements of GF(p^2) (frozen order)")
        return probs
    Pg = np.asarray(rows.get("Prow_gf"), np.int64) % p
    occ = des["occ"]
    if Pg.ndim != 4 or Pg.shape[1:] != (len(occ), des["n"], 2):
        probs.append("rows['Prow_gf'] must have shape (K, M, n, 2)")
        return probs
    E = np.zeros((len(occ), des["n"], 2), np.int64)
    E[np.arange(len(occ)), [v for _f, v in occ], 0] = 1
    if not np.array_equal(Pg[0], E):
        probs.append("P(0) is not the identity on the occurrences")
    if Pg.shape[0] > aux.R and np.any(Pg[aux.R:]):
        probs.append("a GF row has degree above R - 1")
    av, df = np.asarray(des["av"]), np.asarray(des["df"])
    for t in range(1, aux.Q + 1):
        al = nodes[t - 1]
        Pt = np.zeros(Pg.shape[1:], np.int64)
        for k in range(Pg.shape[0]):
            Pt = (Pt + gmul(Pg[k], gpow(al, k, p, c)[None, None, :], p, c)) % p
        Cm = meas_at(des, t) % p
        for g in range(aux.r):
            cols = np.nonzero(av[:, t - 1] == g)[0]
            rws = [i for i, (f, _v) in enumerate(occ) if df[f][t - 1] != g]
            if not cols.size or not rws:
                continue
            for comp in (0, 1):
                if not _in_rowspace_mod(Cm[:, cols], Pt[np.ix_(rws, cols)][..., comp], p):
                    probs.append(f"t={t}, group {g}: GF rows outside the span of the measurements mod {p}")
                    break
        if len(probs) >= 5:
            break
    return probs[:5]


def smallprime_verify(out, Aux, build_design):
    """v8 section of the verify gate: construction.gf_rows(aux, des, p) with p < Q."""
    try:
        import construction as _cons
        gf_rows = getattr(_cons, "gf_rows", None)
    except Exception:
        gf_rows = None
    out["smallprime_supported"] = bool(gf_rows is not None and out.get("modular_supported"))
    out["smallprime_verified"] = False
    if not out["smallprime_supported"]:
        return
    ok_all, ran = True, 0
    for c in CASES["smallprime_counts"]:
        pp, ea = bool(c.get("perpoint")), bool(c.get("eadapt"))
        L, p = int(c.get("L", 0)), int(c["p"])
        if (pp and not out.get("perpoint_supported")) or (ea and not out.get("eadapt_supported")):
            continue
        if L and not out["chunk_supported"]:
            continue
        ran += 1
        rng = np.random.default_rng(c["seed"])
        aux = Aux(*c["cfg"], **pe_kw(pp, ea))
        cl = make_formula(c, rng)
        kw = {"designated": True}
        if L:
            kw.update(units=True, chunk=L)
        des = build_design(aux, c["n"], cl, rng, **kw)
        rows = gf_rows(aux, des, p)
        probs = gf_problems(aux, des, rows, p)
        h0 = ref = None
        if not probs:
            wrng = np.random.default_rng(c["seed"] + 1)
            for _ in range(16):
                w = wrng.integers(1, p, size=c["n"])
                if weighted_brute(c["n"], cl, w, p):
                    break
            got = gf_direct_H0(aux, des, rows, w, p)
            ref = weighted_brute(c["n"], cl, w, p)
            h0 = [int(got[0]), int(got[1])]
            if h0 != [ref, 0]:
                probs.append(f"weighted count {h0} over GF({p}^2), expected [{ref}, 0]")
        ok = not probs
        out["items"].append({"kind": "smallprime_count", "cfg": c["cfg"], "n": c["n"], "p": p, "perpoint": pp,
                             "eadapt": ea, "L": L, "ok": ok, "H0": h0, "reference": ref, "problems": probs[:5]})
        ok_all &= ok
    out["smallprime_verified"] = bool(ok_all and ran > 0)
    out["ok"] &= ok_all


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
    out["units_supported"] = supports_param(build_design, "units") and out["designated_supported"]
    out["units_verified"] = False
    if out["units_supported"]:
        uok, ran = True, 0
        for c in CASES["units_counts"]:
            ran += 1
            rng = np.random.default_rng(c["seed"])
            aux = Aux(*c["cfg"])
            cl = make_formula(c, rng)
            des = build_design(aux, c["n"], cl, rng, designated=True, units=True)
            probs = designated_problems(aux, des, units=True, clauses=cl)
            rep = check_design(aux, des)
            if not (rep["kernel_condition"] and rep["P0_delta"]):
                probs.append("kernel condition or P(0) = delta failed")
            h0, b, J = run_count(aux, des, cl, c["n"])
            desig = np.asarray(des.get("desig", []), dtype=bool)
            nd = sum(1 for i, (f, _v) in enumerate(des["occ"]) if f < len(cl) and i < desig.size and not desig[i])
            j_exp = nd * (aux.R - 1) // aux.Q + 1
            if J != j_exp:
                probs.append(f"Taylor order J = {J}, expected {j_exp} from the non-designated clause occurrences")
            ok = (h0 == b) and not probs
            out["items"].append({"kind": "units_count", "cfg": c["cfg"], "n": c["n"], "ok": ok, "H0": h0,
                                 "brute": b, "J": J, "J_expected": j_exp, "problems": probs[:5]})
            uok &= ok
        out["units_verified"] = bool(uok and ran > 0)
        out["ok"] &= uok
    out["chunk_supported"] = bool(out["units_supported"]) and supports_param(build_design, "chunk")
    out["chunk_verified"] = False
    if out["chunk_supported"]:
        cok, ran = True, 0
        for c in CASES["chunk_counts"]:
            ran += 1
            rng = np.random.default_rng(c["seed"])
            aux = Aux(*c["cfg"])
            cl = make_formula(c, rng)
            L = int(c["L"])
            des = build_design(aux, c["n"], cl, rng, designated=True, units=True, chunk=L)
            probs = designated_problems(aux, des, units=True, clauses=cl, chunk=L) + designation_problems(des)
            rep = check_design(aux, des)
            if not (rep["kernel_condition"] and rep["P0_delta"]):
                probs.append("kernel condition or P(0) = delta failed")
            h0, b, J = run_count(aux, des, cl, c["n"])
            desig = np.asarray(des.get("desig", []), dtype=bool)
            nd = sum(1 for i, (f, _v) in enumerate(des["occ"]) if f < len(cl) and i < desig.size and not desig[i])
            j_exp = nd * (aux.R - 1) // aux.Q + 1
            if J != j_exp:
                probs.append(f"Taylor order J = {J}, expected {j_exp} from the non-designated clause occurrences")
            if c["type"] == "bounded" and int(des.get("M_eff", 0)) > (1 + 1 / L) * c["n"] + 1e-9:
                probs.append(f"M_eff = {des.get('M_eff')} above (1 + 1/L) n")
            ok = (h0 == b) and not probs
            out["items"].append({"kind": "chunk_count", "cfg": c["cfg"], "n": c["n"], "L": L, "ok": ok, "H0": h0,
                                 "brute": b, "J": J, "J_expected": j_exp, "M_eff": int(des.get("M_eff", -1)),
                                 "problems": probs[:5]})
            cok &= ok
        out["chunk_verified"] = bool(cok and ran > 0)
        out["ok"] &= cok
    out["modular_supported"] = bool(out["aggregated_supported"] and out["designated_supported"]
                                    and supports_param(Counter.__init__, "residue_keys")
                                    and supports_param(Counter.__init__, "weights")
                                    and hasattr(Counter, "fibre_counts"))
    out["modular_verified"] = False
    if out["modular_supported"]:
        mok, ran, merged = True, 0, False
        for c in CASES["modular_counts"]:
            L = int(c.get("L", 0))
            if L and not out["chunk_supported"]:
                continue
            ran += 1
            rng = np.random.default_rng(c["seed"])
            aux = Aux(*c["cfg"])
            cl = make_formula(c, rng)
            kw = {"designated": True}
            if L:
                kw.update(units=True, chunk=L)
            des = build_design(aux, c["n"], cl, rng, **kw)
            p = int(c["p"])
            # one extra measurement row p * (1, ..., 1): a valid, redundant measurement that is 0 mod p,
            # so exact integer keys split fibres that residue keys must keep together
            des = dict(des, Chat=np.vstack([np.asarray(des["Chat"], np.int64),
                                            p * np.ones((1, c["n"]), np.int64)]))
            wrng = np.random.default_rng(c["seed"] + 1)
            for _ in range(16):                     # weights with a nonzero weighted count when possible
                w = wrng.integers(1, p, size=c["n"])
                if weighted_brute(c["n"], cl, w, p):
                    break
            C = Counter(aux, des, p=p, residue_keys=True, weights=w)
            data, probs = {}, []
            for t in range(1, aux.Q + 1):
                data[t] = C.taylor_aggregated(t, Zagg=lambda t_, S, U: C.Zagg_group_separated(t_, S, U))
                got = [int(x) for x in C.fibre_counts(t)]
                exp = fibres_frozen(aux, des, t, p)
                if got != exp:
                    probs.append(f"t={t}: fibres per group {got}, expected {exp} (keys mod {p})")
                if sum(exp) < sum(fibres_frozen(aux, des, t, None)):
                    merged = True
            h0, ref = int(C.hermite_H0(data)) % p, weighted_brute(c["n"], cl, w, p)
            ok = (h0 == ref) and not probs
            out["items"].append({"kind": "modular_count", "cfg": c["cfg"], "n": c["n"], "p": p, "L": L, "ok": ok,
                                 "H0": h0, "weighted_brute": ref, "problems": probs[:3]})
            mok &= ok
        out["modular_merged"] = bool(merged)
        out["modular_verified"] = bool(mok and ran > 0 and merged)
        out["ok"] &= mok
    pe_verify(out, Aux, build_design, Counter)
    contract_verify(out, Aux, build_design, Counter)
    smallprime_verify(out, Aux, build_design)
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
    if supports_designated(build_design) and supports_param(build_design, "units"):
        for d in CASES["units_measure"]:
            r, Q, R, A0 = d["cfg"]
            rng = np.random.default_rng(d["seed"])
            aux = Aux(r, Q, R, A0)
            cl = random_bounded_formula(d["n"], rng)
            des = build_design(aux, d["n"], cl, rng, designated=True, units=True)
            probs = designated_problems(aux, des, units=True, clauses=cl)
            rep = check_design(aux, des)
            m_eff, D, h = int(des.get("M_eff", 0)), int(aux.D), int(des["h"])
            hb = 2 * m_eff / D
            deg, checks = common_checks(rep, des["Prow"], r, Q, R)
            checks.update(designated_rows=not probs, copies_h=h <= hb,
                          m_eff_bound=m_eff <= A.M_EFF_UNITS * d["n"] + 1e-9)
            ok = all(checks.values())
            out["items"].append({"cfg": d["cfg"], "n": d["n"], "regime": "units", "ok": ok,
                                 "checks": checks, "M_eff": m_eff, "D": D, "h": h, "h_bound": round(hb, 1),
                                 "problems": probs[:5], "degree": deg})
            out["ok"] &= ok
    if supports_designated(build_design) and supports_param(build_design, "units") and supports_param(build_design, "chunk"):
        chunk_measure(out, Aux, random_bounded_formula, build_design)
    return out


def chunk_measure(out, Aux, random_bounded_formula, build_design):
    for d in CASES["chunk_measure"]:
        r, Q, R, A0 = d["cfg"]
        rng = np.random.default_rng(d["seed"])
        aux = Aux(r, Q, R, A0)
        if d["family"] == "fano":
            n = 7 * (d["n"] // 7)
            cl = fano_formula(n // 7, rng)
        else:
            n = d["n"]
            cl = random_bounded_formula(n, rng)
        L = int(d["L"])
        des = build_design(aux, n, cl, rng, designated=True, units=True, chunk=L)
        probs = designated_problems(aux, des, units=True, clauses=cl, chunk=L) + designation_problems(des)
        m_eff, D, h = int(des.get("M_eff", 0)), int(aux.D), int(des["h"])
        checks = {"chunk_structure": not probs, "m_eff_bound": m_eff <= (1 + 1 / L) * n + 1e-9}
        if d.get("copies", True):
            checks["copies_h"] = h <= 2 * m_eff / D
        ok = all(checks.values())
        out["items"].append({"cfg": d["cfg"], "n": n, "family": d["family"], "L": L, "regime": "chunk", "ok": ok,
                             "checks": checks, "M_eff": m_eff, "bound": round((1 + 1 / L) * n, 1), "D": D, "h": h,
                             "problems": probs[:5]})
        out["ok"] &= ok


def score():
    """Runs the agent's search; records every (lam, params, e) query. The parent process
    re-scores every query with a clean copy of the assembly, so this process is untrusted."""
    from model import choose_params
    try:
        from model import choose_params_decision
    except ImportError:
        choose_params_decision = None
    calls = {"n": 0, "secs": 0.0}
    queries = []

    def f(lam, decision=False):
        fn = choose_params_decision if (decision and choose_params_decision is not None) else choose_params
        t0 = time.time()
        p = fn(int(lam))
        calls["n"] += 1
        calls["secs"] += time.time() - t0
        p = {k: (int(v) if k in ("r", "Q", "R", "A0", "chunk") else v) for k, v in dict(p).items()}
        s = A.exponent_of(p, int(lam), decision=decision)
        e = s["e"] if s["valid"] else 9.0
        queries.append({"lam": int(lam), "params": p, "e": e, "decision": bool(decision)})
        return e, p, s

    def first_lam(target, decision=False, lo=10 ** 3, hi=10 ** 13):
        if f(hi, decision)[0] >= target:
            return None
        if f(lo, decision)[0] < target:
            return lo
        while hi / lo > 1.0005 and hi - lo > 1:
            mid = int(round(math.sqrt(lo * hi)))
            if mid in (lo, hi):
                break
            if f(mid, decision)[0] < target:
                hi = mid
            else:
                lo = mid
        return hi

    res = {}
    for key in ("ppsz", "sharp3sat", "brute"):
        res[key] = first_lam(A.TARGETS[key], decision=key in A.DECISION_TARGETS)
    for key, target in A.REPORT_TARGETS.items():
        res[key] = first_lam(target, decision=True)
    res["decision_model"] = choose_params_decision is not None
    e_ref, p_ref, s_ref = f(2224173, decision=True)
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
# VERSION6_MARK
# VERSION7_MARK
# VERSION8_MARK
