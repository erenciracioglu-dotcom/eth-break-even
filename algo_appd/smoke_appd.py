"""Quick self-check for the Appendix D workspace. Run from the workspace:
    D:/ETH-loop/.venv/Scripts/python.exe smoke_appd.py
Core checks (public seeds; they run once appd_core.IMPLEMENTED is True): node_value with methods paper, memo
and dp against direct_taylor_mod, the memo child count (sum of the table sizes), the dp combination count against
the paper tuple count; byvars_value against child_network_value on chunk networks; count_e2e against count_sat.
Then e of choose_params at three lam values (lower is better; must be valid and under 3 s per call).
Prints SMOKE OK or SMOKE FAIL at the end."""
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


def param_checks():
    ok_all = True
    for lam in (5404, 23056, 201458):
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
