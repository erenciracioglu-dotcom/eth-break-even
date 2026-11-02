"""T26 self-check (round 11, appd-profile): pf_choose on the public bit dicts of public_t26.json, node_value(...,
method="prof") on the public nodes against node_value(..., method="dp") and the head's per-profile records (NOTES.md,
'Round 11'), the shapes of the matmul calls, and count_e2e(..., basis="ldu", method="prof") against count_sat with
sum_prof used at exactly the nodes with 4 to 16 groups. Prints T26 OK when everything matches."""
import json
import os
import random
import time
import appd_core as C

HERE = os.path.dirname(os.path.abspath(__file__))


def canon(recs, mod):
    return [[list(map(int, d.get("weights", []))), list(map(int, d.get("sizes", []))), int(d.get("order", -1)),
             list(map(int, d.get("costs", []))), [list(map(int, s)) for s in d.get("shapes", [])],
             [int(x) % mod for x in d.get("value", [])]] for d in recs]


def main():
    data = json.load(open(os.path.join(HERE, "public_t26.json"), encoding="utf-8"))
    bad = 0
    for cs in data["choose"]:
        bits = {(i, j): b for i, j, b in cs["bits"]}
        got = C.pf_choose(bits, cs["r"])
        ok = int(got[0]) == cs["order"] and [int(x) for x in got[1]] == cs["costs"]
        if not ok:
            bad += 1
            print("pf_choose r", cs["r"], "got", got, "want", cs["order"], cs["costs"])
    print("pf_choose:", len(data["choose"]), "cases,", "ok" if bad == 0 else "%d wrong" % bad)
    for cs in data["node"]:
        t0 = time.time()
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = C.setup(rng, cs["n"], cs["m"], cs["a"], cs["r"], chunks=cs["chunks"], Lc=cs["Lc"])
        nodes = [s for s in pk.S if s and C.v2(s) < pk.a]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q, w, basis = cs["p"], cs["q"], cs["weights"], cs["basis"]
        mod = 1 << p
        log = []
        orig = C.matmul

        def logging(Am, Bm, qq, mm):
            log.append([len(Am), len(Bm), len(Bm[0]) if Bm else 0])
            return orig(Am, Bm, qq, mm)
        C.matmul = logging
        try:
            ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method="prof", basis=basis)
        finally:
            C.matmul = orig
        dp, _ = C.node_value(net, pk, labels, t, p, q, weights=w, method="dp", basis=basis)
        got = canon(st.get("records", []), mod)
        ok_v = [int(x) % mod for x in ser] == [int(x) % mod for x in dp] == cs["value"]
        ok_r = got == cs["records"]
        ok_m = sorted(log) == sorted(s for d in cs["records"] for s in d[4])
        if not (ok_v and ok_r and ok_m):
            bad += 1
            first = next((k for k in range(min(len(got), len(cs["records"]))) if got[k] != cs["records"][k]), None)
            print("node seed", cs["seed"], "value ok", ok_v, "records ok", ok_r, "matmul shapes ok", ok_m, "profiles",
                  len(got), "want", len(cs["records"]), "first difference",
                  None if first is None else (got[first], cs["records"][first]))
        else:
            print("node seed", cs["seed"], "ok:", len(got), "profiles,", round(time.time() - t0, 1), "s")
    for cs in data["e2e"]:
        rng = random.Random(cs["seed"])
        cls = C.random_formula(cs["n"], cs["m"], rng)
        calls = [0]
        orig = C.NodeExpansion.sum_prof

        def counting(self, *a, **k):
            calls[0] += 1
            return orig(self, *a, **k)
        C.NodeExpansion.sum_prof = counting
        try:
            got = C.count_e2e(cs["n"], cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"], a=cs["a"],
                              u=cs["u"], basis="ldu", method="prof")
        finally:
            C.NodeExpansion.sum_prof = orig
        want = C.count_sat(cs["n"], cls)
        ok = int(got) == want and calls[0] == cs["prof_nodes"]
        print("e2e seed", cs["seed"], "count", got, "want", want, "sum_prof calls", calls[0], "want", cs["prof_nodes"])
        bad += 0 if ok else 1
    print("T26 OK" if bad == 0 else "T26 FAIL (%d mismatches)" % bad)


if __name__ == "__main__":
    main()
