"""Quick self-check for the Appendix D workspace. Run from the workspace:
    D:/ETH-loop/.venv/Scripts/python.exe smoke_appd.py
Core checks (public seeds; they run once appd_core.IMPLEMENTED is True): node_value with methods paper, memo
and dp against direct_taylor_mod, the memo child count (sum of the table sizes), the dp combination count against
the paper tuple count; byvars_value against child_network_value on chunk networks; count_e2e against count_sat.
Newton checks (round 6): the same with NodeExpansion(..., basis="newton") at a = 2, 3 and 4, the weights
w_b = (b + 1 + v2(b!)) or (a + 1 - i + b + v2(b!)), and count_e2e(..., a, u, basis="newton") at a = 2 and 3.
Then e of choose_params at three lam values (lower is better; must be valid and under 3 s per call).
Prints SMOKE OK or SMOKE FAIL (core and parameters), then NEWTON OK or NEWTON FAIL."""
import math
import random
import time

import appd_assembly_public as A
import appd_core as C
import appd_params as P


def core_checks():
    if not C.IMPLEMENTED:
        print("appd_core.IMPLEMENTED is False: core checks skipped")
        return True
    ok_all = True
    pick = random.Random(7)
    for k in range(14):
        seed, chunks = 100 + k, k % 3 != 0
        n_, cls, net, labels, pk = C.setup(random.Random(seed), pick.randint(4, 6), pick.randint(3, 4),
                                           pick.choice([1, 2, 2, 3]), pick.choice([3, 4]), chunks=chunks, Lc=2)
        M = len(net.wires)
        nodes = [s for s in pk.S if s and C.v2(s) < pk.a]
        t, p, q, w = pick.choice(nodes), pick.choice([4, 5]), pick.choice([1, 2]), pick.choice(["paper", "sharp"])
        if M < 2 or M > 7:
            continue
        direct = C.direct_taylor_mod(net, pk, t, p, q)
        line, stats = [], {}
        for method in ("paper", "memo", "dp"):
            try:
                ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method=method)
                good = [int(x) % (1 << p) for x in ser] == direct
                stats[method] = st
            except Exception as exc:  # noqa: BLE001
                good = False
                stats[method] = repr(exc)[:160]
            line.append(f"{method}={'ok' if good else 'FAIL'}")
            ok_all &= good
        ne = C.NodeExpansion(net, pk, labels, t, p, q, weights=w)
        sizes = sum(len(ne.tuples(ne.inc[g], ne.p)) for g in ne.r_groups)
        if isinstance(stats.get("memo"), dict) and stats["memo"].get("child_calls") != sizes:
            line.append(f"memo child_calls {stats['memo'].get('child_calls')} != {sizes} FAIL")
            ok_all = False
        print(f"case {seed}: chunks={chunks} M={M} t={t} p={p} q={q} modes={len(ne.modes)} " + " ".join(line))
        if chunks:
            for g in ne.r_groups[:2]:
                for gm in ne.tuples(ne.inc[g], ne.p)[:3]:
                    ref = ne.child_network_value(g, dict(gm))
                    try:
                        good = [int(x) % (1 << p) for x in C.byvars_value(ne, g, dict(gm))[0]] == ref
                    except Exception as exc:  # noqa: BLE001
                        good = False
                        print("   byvars error:", repr(exc)[:160])
                    ok_all &= good
                    if not good:
                        print(f"   byvars group {g} gamma {gm}: FAIL")
    for k in range(4):
        seed, n = 200 + k, 3 + k % 2
        cls = C.random_formula(n, 2 + k % 2, random.Random(seed))
        try:
            got = C.count_e2e(n, cls, chunks=k % 2 == 1, Lc=2, r=3, seed=seed)
        except Exception as exc:  # noqa: BLE001
            got = repr(exc)[:160]
        want = C.count_sat(n, cls)
        print(f"count_e2e seed {seed}: got {got} want {want} {'ok' if got == want else 'FAIL'}")
        ok_all &= got == want
    return ok_all


def v2fact(b):
    return b - bin(b).count("1")


def newton_checks():
    ok_all = True
    pick = random.Random(11)
    done = 0
    for k in range(40):
        if done >= 10:
            break
        seed, chunks = 300 + k, k % 2 == 1
        a = pick.choice([2, 3, 3, 4])
        n_, cls, net, labels, pk = C.setup(random.Random(seed), pick.randint(4, 6), pick.randint(3, 4), a,
                                           pick.choice([3, 4]), chunks=chunks, Lc=2)
        M = len(net.wires)
        nodes = [s for s in pk.S if s and C.v2(s) < pk.a]
        t, p, q, w = pick.choice(nodes), pick.choice([5, 6]), pick.choice([1, 2]), pick.choice(["paper", "sharp"])
        if M < 2 or M > 7:
            continue
        done += 1
        direct = C.direct_taylor_mod(net, pk, t, p, q)
        line = []
        try:
            ne = C.NodeExpansion(net, pk, labels, t, p, q, weights=w, basis="newton")
        except Exception as exc:  # noqa: BLE001
            print(f"newton case {seed}: NodeExpansion(..., basis='newton') failed: {repr(exc)[:200]} FAIL")
            ok_all = False
            continue
        i = C.v2(t)
        bad_w = [md for md in ne.modes
                 if md["w"] != ((md["b"] + 1) if w == "paper" else (a + 1 - i + md["b"])) + v2fact(md["b"])]
        step = 1 << (i + 1)
        bad_lin = [md for md in ne.modes
                   if any(c != math.comb(pk.where[w2][1] // step, md["b"]) for c, w2 in md["lin"])]
        if bad_w:
            line.append(f"weights FAIL ({len(bad_w)} modes)")
            ok_all = False
        if bad_lin:
            line.append(f"source forms FAIL ({len(bad_lin)} modes)")
            ok_all = False
        for method in ("paper", "memo", "dp"):
            try:
                ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method=method, basis="newton")
                good = [int(x) % (1 << p) for x in ser] == direct
            except Exception as exc:  # noqa: BLE001
                good = False
                line.append(repr(exc)[:120])
            line.append(f"{method}={'ok' if good else 'FAIL'}")
            ok_all &= good
        pw = C.NodeExpansion(net, pk, labels, t, p, q, weights=w)
        print(f"newton case {seed}: a={a} chunks={chunks} M={M} t={t} p={p} q={q} modes newton {len(ne.modes)} "
              f"power {len(pw.modes)} " + " ".join(line))
        if chunks:
            for g in ne.r_groups[:2]:
                for gm in ne.tuples(ne.inc[g], ne.p)[:3]:
                    ref = ne.child_network_value(g, dict(gm))
                    try:
                        good = [int(x) % (1 << p) for x in C.byvars_value(ne, g, dict(gm))[0]] == ref
                    except Exception as exc:  # noqa: BLE001
                        good = False
                        print("   newton byvars error:", repr(exc)[:160])
                    ok_all &= good
                    if not good:
                        print(f"   newton byvars group {g} gamma {gm}: FAIL")
    for k, (a, u) in enumerate([(2, 1), (2, 2), (2, 1), (2, 2), (3, 1), (3, 1)]):
        seed, n = 400 + k, 3
        cls = C.random_formula(n, 1 + k % 2, random.Random(seed))
        try:
            got = C.count_e2e(n, cls, chunks=k % 2 == 1, Lc=2, r=3, seed=seed, a=a, u=u, basis="newton")
        except Exception as exc:  # noqa: BLE001
            got = repr(exc)[:160]
        want = C.count_sat(n, cls)
        print(f"newton count_e2e seed {seed} a={a} u={u}: got {got} want {want} {'ok' if got == want else 'FAIL'}")
        ok_all &= got == want
    return ok_all


def param_checks():
    ok_all = True
    for lam in (2700, 5300, 23056):
        t0 = time.time()
        p = P.choose_params(lam, "third")
        dt = time.time() - t0
        s = A.exponent_of(p, lam)
        good = s["valid"] and dt < 3.0
        ok_all &= good
        print(f"lam={lam}: e={s['e'] if s['valid'] else 'INVALID: ' + s.get('why', '')} ({dt:.2f} s) params={p}")
    return ok_all


if __name__ == "__main__":
    a = core_checks()
    b = param_checks()
    print("SMOKE OK" if a and b else "SMOKE FAIL")
    nw = newton_checks() if C.IMPLEMENTED else False
    print("NEWTON OK" if nw else "NEWTON FAIL")
