"""FROZEN gates of referee_appd version 5 (round 8). Usage: python appd_gates.py <candidate_dir> <verify|score> <out_json> [verified_flags_json]

verify: the agent's appd_core.py against the frozen reference appd_ref.py on hidden cases (cases_appd.json):
  node   - node_value(method = paper, memo, dp) equals F(t + alpha z) mod (2^p, z^q) from the frozen direct
           evaluation of (D.8); memo reports the reference number of distinct children, dp the reference
           number of enumerated combinations; on chunk networks byvars_value equals the reference child tables
  e2e    - count_e2e equals the frozen brute-force count
  newton (version 2; these items only decide the newton flag, not the gate) - NodeExpansion(..., basis="newton")
           has exactly the reference modes (packet, tile y0, b,
           source, target, weight, source form) and the reference coefficient residues; node_value(..., basis=
           "newton") equals (D.8) for paper, memo and dp with the reference child and enumeration counts; byvars on
           chunk networks; count_e2e(..., a, u, basis="newton") equals the brute-force count at a = 2 and 3.
  ldu (version 3; decides only the ldu flag) - the same checks on the same hidden cases with basis="ldu" and
           weights="sharp": modes, weights, series source forms and coefficient residues equal the reference; node
           values by the three orders with the reference counts; byvars on chunk networks; end-to-end counts.
  shellmm (version 4; decides only the shellmm flag) - matmul equals the schoolbook product on 11 shapes and
           uses at most 7^(k-1) * 8 ring products (counted through ps_mul_mod) for 2^k x 2^k inputs;
           node_value(..., method="mm") equals (D.8) on the node cases (power basis) and the newton node
           cases (ldu basis), its per-step block shapes and the shapes of every matmul call equal the
           reference shell structure and its max_enum equals the reference dp count; count_e2e(...,
           basis="ldu", method="mm") equals the brute-force count on the newton end-to-end cases.
  guarded (version 5; decides only the guarded flag) - on hidden guarded cases (reference GuardedPacking with
           virtual size 2^a, guards 2^(a - kappa), labels on u levels) NodeExpansion(..., sharp, ldu) has exactly
           the reference modes, weights, series source forms and coefficient residues; node_value(method dp
           and mm) equals (D.8) with the guarded deformation and the reference dp count; byvars on chunk
           networks; count_e2e(..., guarded=True, kappa) equals the brute-force count.
  direct_taylor_mod and F_poly of the agent module are blocked while it runs.
  Verified flags: memo (every memo item), modedp (every dp item), chunks (every item on a chunk network and every
  chunk end-to-end count), leafvars (every byvars item), newton (every newton item), ldu (every ldu item),
  shellmm (every mm item), guarded (every gd item).
  With IMPLEMENTED = False
  nothing is tried and no code flag is verified (the gate still passes).
score: break-even sizes from the agent's appd_params.choose_params(lam, target), re-scored here with the frozen
  appd_assembly.py; parameters that use an unverified code flag count as invalid.
"""
import json
import math
import os
import random
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import appd_assembly as A  # noqa: E402
import appd_ref as R  # noqa: E402


def _blocked(*_a, **_k):
    raise RuntimeError("blocked by the referee during verification")


def _mode_table(ne):
    """(P, y0, b) -> (source, target, weight, sorted source form, {row wire: residue series})."""
    out = {}
    for k, md in enumerate(ne.modes):
        key = (int(md["P"]), int(md["y0"]), int(md["b"]))
        rows = {int(wr): [int(c) % ne.m for c in ser] for (wr, kk), ser in ne.coef.items() if kk == k}
        lin = sorted(((tuple(int(x) % ne.m for x in c) if isinstance(c, (list, tuple)) else int(c)), int(w2))
                     for c, w2 in md["lin"])
        out[key] = (md["src"], md["tgt"], int(md["w"]), lin, rows)
    return out


def newton_items(C, cases, basis="newton"):
    items = []
    for cs in cases.get("newton_node", []):
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = R.setup(rng, cs["n"], cs["m"], cs["a"], cs["r"], chunks=cs["chunks"], Lc=cs["Lc"])
        M = len(net.wires)
        if M < 2 or M > 8:
            continue
        C._WIREMAP.update({k: v for k, v in R._WIREMAP.items() if k[0] == id(net)})
        nodes = [s for s in pk.S if s and R.v2(s) < pk.a]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q, w = cs["p"], cs["q"], (cs["weights"] if basis == "newton" else "sharp")
        direct = R.direct_taylor_mod(net, pk, t, p, q)
        ref = R.NodeExpansion(net, pk, labels, t, p, q, weights=w, basis=basis)
        _, _, ref_calls = ref.sum_memo()
        _, ref_enum, _ = ref.sum_dp()
        base = dict(seed=cs["seed"], chunks=cs["chunks"], a=pk.a, M=M, t=t, p=p, q=q, w=w, modes=len(ref.modes))
        it = dict(base, kind=basis + "_modes")
        ne = None
        try:
            ne = C.NodeExpansion(net, pk, labels, t, p, q, weights=w, basis=basis)
            mine, want = _mode_table(ne), _mode_table(ref)
            it["ok"] = mine == want
            if not it["ok"]:
                it["diff"] = str(sorted(set(mine) ^ set(want)))[:200] or "same modes, different weights/forms/coefficients"
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
        for method in ("paper", "memo", "dp"):
            it = dict(base, kind=basis + "_" + method)
            try:
                ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method=method, basis=basis)
                ok = [int(x) % (1 << p) for x in ser] == direct
                if method == "memo":
                    ok = ok and int(st.get("child_calls", -1)) == ref_calls
                    it["child_calls"] = [st.get("child_calls"), ref_calls]
                if method == "dp":
                    ok = ok and int(st.get("max_enum", -1)) == ref_enum
                    it["max_enum"] = [st.get("max_enum"), ref_enum]
                it["ok"] = bool(ok)
            except Exception:  # noqa: BLE001
                it.update(ok=False, error=traceback.format_exc()[-400:])
            items.append(it)
        if cs["chunks"] and ne is not None:
            for g in ref.r_groups[:2]:
                tab = ref.group_table(g)
                for gm in sorted(tab.keys())[:3]:
                    it = dict(base, kind=basis + "_byvars", g=g)
                    try:
                        bv = C.byvars_value(ne, g, dict(gm))[0]
                        it["ok"] = [int(x) % (1 << p) for x in bv] == tab[gm]
                    except Exception:  # noqa: BLE001
                        it.update(ok=False, error=traceback.format_exc()[-400:])
                    items.append(it)
    modes_bad = any(i["kind"] == basis + "_modes" and not i["ok"] for i in items)
    for cs in cases.get("newton_e2e", []):
        rng = random.Random(cs["seed"])
        n = cs["n"]
        cls = R.random_formula(n, cs["m"], rng)
        it = dict(kind=basis + "_e2e", seed=cs["seed"], chunks=cs["chunks"], a=cs["a"], u=cs["u"])
        if modes_bad:      # a wrong basis can make these counts very slow; they cannot verify the flag anyway
            it.update(ok=False, error="not run: the candidate's " + basis + " modes differ from the reference")
            items.append(it)
            continue
        try:
            got = C.count_e2e(n, cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"], a=cs["a"],
                              u=cs["u"], basis=basis)
            want = R.count_sat(n, cls)
            it.update(ok=int(got) == want, got=int(got), want=want)
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
    return items


MM_SHAPES = ((1, 1, 1), (1, 3, 2), (2, 3, 1), (3, 1, 4), (5, 3, 7), (12, 4, 9), (7, 13, 6), (4, 4, 4), (8, 8, 8),
             (16, 16, 16), (32, 32, 32))
MM_MAX_PRODUCTS = {4: 7 * 8, 8: 7 ** 2 * 8, 16: 7 ** 3 * 8, 32: 7 ** 4 * 8}


def mm_items(C, cases):
    """Round 7 (shellmm): matmul exact against the schoolbook product, with at most 7^(k-1) * 8 ring products for
    2^k x 2^k inputs; node values through shell matrix products, with the reference block shapes in the stats and
    in every matmul call; end-to-end counts with method mm (ldu basis)."""
    if not callable(getattr(C, "matmul", None)):
        return [dict(kind=k, ok=False, error="appd_core has no matmul (round 7 not implemented)")
                for k in ("mm_matmul", "mm_node", "mm_e2e")]
    items = []
    rng = random.Random(7007)
    for (n, k, l) in MM_SHAPES:
        q, m = rng.choice([1, 2, 3]), 1 << rng.choice([3, 5, 8])
        Am = [[[rng.randrange(m) for _ in range(q)] for _ in range(k)] for _ in range(n)]
        Bm = [[[rng.randrange(m) for _ in range(q)] for _ in range(l)] for _ in range(k)]
        it = dict(kind="mm_matmul", shape=[n, k, l], q=q, m=m)
        cnt = [0]
        orig = C.ps_mul_mod

        def counting(a, b, qq, mm, _o=orig):
            cnt[0] += 1
            return _o(a, b, qq, mm)
        try:
            want = R.naive_matmul(Am, Bm, q, m)
            C.ps_mul_mod = counting
            try:
                got = C.matmul(Am, Bm, q, m)
            finally:
                C.ps_mul_mod = orig
            ok = [[[int(x) % m for x in s] for s in row] for row in got] == want
            if n == k == l and n in MM_MAX_PRODUCTS:
                ok = ok and cnt[0] <= MM_MAX_PRODUCTS[n]
            it.update(ok=bool(ok), ring_products=cnt[0])
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
    jobs = [(cs, "power") for cs in cases["node"]] + [(cs, "ldu") for cs in cases.get("newton_node", [])]
    for cs, basis in jobs:
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = R.setup(rng, cs["n"], cs["m"], cs["a"], cs["r"], chunks=cs["chunks"], Lc=cs["Lc"])
        M = len(net.wires)
        if M < 2 or M > 8:
            continue
        C._WIREMAP.update({k: v for k, v in R._WIREMAP.items() if k[0] == id(net)})
        nodes = [s for s in pk.S if s and R.v2(s) < pk.a]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q = cs["p"], cs["q"]
        w = "sharp" if basis == "ldu" else cs["weights"]
        direct = R.direct_taylor_mod(net, pk, t, p, q)
        ref = R.NodeExpansion(net, pk, labels, t, p, q, weights=w, basis=basis)
        _, ref_enum, _ = ref.sum_dp()
        want = ref.mm_shapes()
        it = dict(kind="mm_node", basis=basis, seed=cs["seed"], chunks=cs["chunks"], M=M, t=t, p=p, q=q,
                  modes=len(ref.modes), blocks=sum(len(s) for s in want))
        log = []
        orig_mm = C.matmul

        def logging_mm(Am, Bm, qq, mm, _o=orig_mm):
            log.append((len(Am), len(Bm), len(Bm[0]) if Bm else 0))
            return _o(Am, Bm, qq, mm)
        C.matmul = logging_mm
        try:
            ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method="mm", basis=basis)
            ok = [int(x) % (1 << p) for x in ser] == direct
            got_shapes = [sorted(tuple(int(v) for v in x) for x in s) for s in st.get("shapes", [])]
            ok = ok and got_shapes == want and sorted(log) == sorted(x for s in want for x in s)
            ok = ok and int(st.get("max_enum", -1)) == ref_enum
            it.update(ok=bool(ok), max_enum=[st.get("max_enum"), ref_enum], calls=len(log))
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        finally:
            C.matmul = orig_mm
        items.append(it)
    node_bad = any(i["kind"] == "mm_node" and not i["ok"] for i in items)
    for cs in cases.get("newton_e2e", []):
        rng = random.Random(cs["seed"])
        n = cs["n"]
        cls = R.random_formula(n, cs["m"], rng)
        it = dict(kind="mm_e2e", seed=cs["seed"], chunks=cs["chunks"], a=cs["a"], u=cs["u"])
        if node_bad:
            it.update(ok=False, error="not run: a shell matrix product node check failed")
            items.append(it)
            continue
        try:
            got = C.count_e2e(n, cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"], a=cs["a"],
                              u=cs["u"], basis="ldu", method="mm")
            want_c = R.count_sat(n, cls)
            it.update(ok=int(got) == want_c, got=int(got), want=want_c)
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
    return items


GUARD_KINDS = ("gd_modes", "gd_dp", "gd_mm", "gd_byvars", "gd_e2e")


def guarded_items(C, cases):
    """Round 8 (guarded): on the hidden guarded cases the agent's NodeExpansion with the reference GuardedPacking has
    exactly the reference modes (weights, series source forms with the odd units, coefficient residues with the row
    factors); node_value(method dp and mm) equals (D.8) with the guarded deformation and the reference dp count;
    byvars on chunk networks; count_e2e(..., guarded=True, kappa) equals the brute-force count."""
    if not callable(getattr(C, "GuardedPacking", None)):
        return [dict(kind=k, ok=False, error="appd_core has no GuardedPacking (round 8 not implemented)")
                for k in GUARD_KINDS]
    items = []
    for cs in cases.get("guarded_node", []):
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = R.setup_g(rng, cs["n"], cs["m"], cs["a"], cs["u"], cs["kappa"], cs["r"],
                                             chunks=cs["chunks"], Lc=cs["Lc"])
        M = len(net.wires)
        C._WIREMAP.update({k: v for k, v in R._WIREMAP.items() if k[0] == id(net)})
        nodes = [s for s in pk.S if s and R.v2(s) < pk.u]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q = cs["p"], cs["q"]
        direct = R.direct_taylor_mod_g(net, pk, t, p, q)
        ref = R.NodeExpansionG(net, pk, labels, t, p, q)
        _, ref_enum, _ = ref.sum_dp()
        base = dict(seed=cs["seed"], chunks=cs["chunks"], a=cs["a"], u=cs["u"], kappa=cs["kappa"], M=M, t=t, p=p, q=q,
                    modes=len(ref.modes))
        it = dict(base, kind="gd_modes")
        ne = None
        try:
            ne = C.NodeExpansion(net, pk, labels, t, p, q, weights="sharp", basis="ldu")
            mine, want = _mode_table(ne), _mode_table(ref)
            it["ok"] = mine == want
            if not it["ok"]:
                it["diff"] = str(sorted(set(mine) ^ set(want)))[:200] or "same modes, different weights/forms/coefficients"
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
        for method in ("dp", "mm"):
            it = dict(base, kind="gd_" + method)
            try:
                ser, st = C.node_value(net, pk, labels, t, p, q, weights="sharp", method=method, basis="ldu")
                ok = [int(x) % (1 << p) for x in ser] == direct and int(st.get("max_enum", -1)) == ref_enum
                it.update(ok=bool(ok), max_enum=[st.get("max_enum"), ref_enum])
            except Exception:  # noqa: BLE001
                it.update(ok=False, error=traceback.format_exc()[-400:])
            items.append(it)
        if cs["chunks"] and ne is not None:
            for g in ref.r_groups[:2]:
                tab = ref.group_table(g)
                for gm in sorted(tab.keys())[:3]:
                    it = dict(base, kind="gd_byvars", g=g)
                    try:
                        bv = C.byvars_value(ne, g, dict(gm))[0]
                        it["ok"] = [int(x) % (1 << p) for x in bv] == tab[gm]
                    except Exception:  # noqa: BLE001
                        it.update(ok=False, error=traceback.format_exc()[-400:])
                    items.append(it)
    modes_bad = any(i["kind"] == "gd_modes" and not i["ok"] for i in items)
    for cs in cases.get("guarded_e2e", []):
        rng = random.Random(cs["seed"])
        n = cs["n"]
        cls = R.random_formula(n, cs["m"], rng)
        it = dict(kind="gd_e2e", seed=cs["seed"], chunks=cs["chunks"], a=cs["a"], u=cs["u"], kappa=cs["kappa"],
                  method=cs["method"])
        if modes_bad:
            it.update(ok=False, error="not run: the candidate's guarded modes differ from the reference")
            items.append(it)
            continue
        try:
            got = C.count_e2e(n, cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"], a=cs["a"],
                              u=cs["u"], basis="ldu", method=cs["method"], guarded=True, kappa=cs["kappa"])
            want_c = R.count_sat(n, cls)
            it.update(ok=int(got) == want_c, got=int(got), want=want_c)
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
    return items


def verify(cand):
    sys.path.insert(0, cand)
    import appd_core as C  # noqa: E402
    cases = json.load(open(os.path.join(HERE, "cases_appd.json")))
    items = []
    if not getattr(C, "IMPLEMENTED", False):
        return {"ok": True, "implemented": False, "items": [], "flags": {k: False for k in A.CODE_FLAGS}}
    C.direct_taylor_mod = _blocked
    C.F_poly = _blocked
    for cs in cases["node"]:
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = R.setup(rng, cs["n"], cs["m"], cs["a"], cs["r"], chunks=cs["chunks"], Lc=cs["Lc"])
        M = len(net.wires)
        if M < 2 or M > 8:
            continue
        C._WIREMAP.update({k: v for k, v in R._WIREMAP.items() if k[0] == id(net)})
        nodes = [s for s in pk.S if s and R.v2(s) < pk.a]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q, w = cs["p"], cs["q"], cs["weights"]
        direct = R.direct_taylor_mod(net, pk, t, p, q)
        ref = R.NodeExpansion(net, pk, labels, t, p, q, weights=w)
        _, _, ref_calls = ref.sum_memo()
        _, ref_enum, _ = ref.sum_dp()
        base = dict(seed=cs["seed"], chunks=cs["chunks"], M=M, t=t, p=p, q=q, modes=len(ref.modes))
        for method in ("paper", "memo", "dp"):
            it = dict(base, kind=method)
            try:
                ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method=method)
                ok = [int(x) % (1 << p) for x in ser] == direct
                if method == "memo":
                    ok = ok and int(st.get("child_calls", -1)) == ref_calls
                    it["child_calls"] = [st.get("child_calls"), ref_calls]
                if method == "dp":
                    ok = ok and int(st.get("max_enum", -1)) == ref_enum
                    it["max_enum"] = [st.get("max_enum"), ref_enum]
                it["ok"] = bool(ok)
            except Exception:  # noqa: BLE001
                it.update(ok=False, error=traceback.format_exc()[-400:])
            items.append(it)
        if cs["chunks"]:
            ne = C.NodeExpansion(net, pk, labels, t, p, q, weights=w)
            for g in ref.r_groups[:2]:
                tab = ref.group_table(g)
                for gm in sorted(tab.keys())[:3]:
                    it = dict(base, kind="byvars", g=g)
                    try:
                        bv = C.byvars_value(ne, g, dict(gm))[0]
                        it["ok"] = [int(x) % (1 << p) for x in bv] == tab[gm]
                    except Exception:  # noqa: BLE001
                        it.update(ok=False, error=traceback.format_exc()[-400:])
                    items.append(it)
    for cs in cases["e2e"]:
        rng = random.Random(cs["seed"])
        n = cs["n"]
        cls = R.random_formula(n, cs["m"], rng)
        it = dict(kind="e2e", seed=cs["seed"], chunks=cs["chunks"])
        try:
            got = C.count_e2e(n, cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"])
            want = R.count_sat(n, cls)
            it.update(ok=int(got) == want, got=int(got), want=want)
        except Exception:  # noqa: BLE001
            it.update(ok=False, error=traceback.format_exc()[-400:])
        items.append(it)
    items += newton_items(C, cases)
    items += newton_items(C, cases, basis="ldu")
    items += mm_items(C, cases)
    items += guarded_items(C, cases)
    nk = ("newton_modes", "newton_paper", "newton_memo", "newton_dp", "newton_byvars", "newton_e2e")
    lk = tuple("ldu" + k[len("newton"):] for k in nk)
    mk = ("mm_matmul", "mm_node", "mm_e2e")
    gk = GUARD_KINDS
    core = ("paper", "memo", "dp", "byvars", "e2e")       # the newton items only decide the newton flag
    flags = {
        "memo": all(i["ok"] for i in items if i["kind"] == "memo") and any(i["kind"] == "memo" for i in items),
        "modedp": all(i["ok"] for i in items if i["kind"] == "dp") and any(i["kind"] == "dp" for i in items),
        "chunks": all(i["ok"] for i in items if i.get("chunks") and i["kind"] in core)
        and any(i.get("chunks") for i in items),
        "leafvars": all(i["ok"] for i in items if i["kind"] == "byvars") and any(i["kind"] == "byvars" for i in items),
        "newton": all(i["ok"] for i in items if i["kind"] in nk) and all(any(i["kind"] == k for i in items) for k in nk),
        "ldu": all(i["ok"] for i in items if i["kind"] in lk) and all(any(i["kind"] == k for i in items) for k in lk),
        "shellmm": all(i["ok"] for i in items if i["kind"] in mk) and all(any(i["kind"] == k for i in items) for k in mk),
        "guarded": all(i["ok"] for i in items if i["kind"] in gk) and all(any(i["kind"] == k for i in items) for k in gk),
    }
    return {"ok": all(i["ok"] for i in items if i["kind"] in core), "implemented": True, "items": items, "flags": flags,
            "counts": {k: sum(1 for i in items if i["kind"] == k) for k in ("paper", "memo", "dp", "byvars", "e2e") + nk + lk + mk + gk}}


def score(cand, verified):
    sys.path.insert(0, cand)
    import appd_params as P  # noqa: E402
    calls, spent = 0, 0.0

    def ev(lam, key):
        nonlocal calls, spent
        t0 = time.time()
        p = P.choose_params(int(lam), key)
        spent += time.time() - t0
        calls += 1
        s = A.exponent_of(p, int(lam))
        if not s["valid"]:
            return 9.0, p
        if any(s.get(f) and not verified.get(f) for f in A.CODE_FLAGS):
            return 9.0, p
        return s["e"], p

    out = {"thresholds": {}, "certificates": {}}
    for key, target in A.TARGETS.items():
        lo, hi = 64, 10 ** 9
        e_hi, p_hi = ev(hi, key)
        if e_hi >= target:
            out["thresholds"][key] = None
            continue
        e_lo, _ = ev(lo, key)
        if e_lo < target:
            out["thresholds"][key] = lo
            continue
        while hi - lo > 1:
            mid = int(math.sqrt(lo * hi)) if hi / lo > 1.01 else (lo + hi) // 2
            if mid in (lo, hi):
                mid = (lo + hi) // 2
            e, p = ev(mid, key)
            if e < target:
                hi, p_hi, e_hi = mid, p, e
            else:
                lo = mid
        out["thresholds"][key] = hi
        out["certificates"][key] = {"lam": hi, "params": p_hi, "e": e_hi}
    out["choose_params_avg_s"] = spent / max(calls, 1)
    out["calls"] = calls
    return out


if __name__ == "__main__":
    cand, gate, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    try:
        if gate == "verify":
            res = verify(cand)
        else:
            verified = json.loads(sys.argv[4]) if len(sys.argv) > 4 else {}
            res = score(cand, verified)
    except Exception:  # noqa: BLE001
        res = {"ok": False, "error": traceback.format_exc()[-1500:]}
    res["seconds"] = round(time.time() - t0, 1)
    with open(outp, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
