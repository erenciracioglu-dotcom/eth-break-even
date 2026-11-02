"""T24 self-check (round 10, appd-bincross part 2): node_value(..., method="cross") on the public nodes of
public_t24.json against node_value(..., method="dp") and the head's per-piece records (NOTES.md, 'Round 10'), and
count_e2e(..., basis="ldu", method="cross") against count_sat with sum_cross used at exactly the nodes with at
least four groups. Prints T24 OK when everything matches."""
import json
import os
import random
import time
import appd_core as C

HERE = os.path.dirname(os.path.abspath(__file__))
KEYS = ("wx", "wy", "nx", "ny", "a", "b", "m", "i", "pairing", "channels", "extra")


def canon(recs, mod):
    return [tuple(d.get(k) for k in KEYS) + (tuple(int(x) % mod for x in d.get("value", [])),) for d in recs]


def main():
    data = json.load(open(os.path.join(HERE, "public_t24.json"), encoding="utf-8"))
    bad = 0
    for cs in data["node"]:
        t0 = time.time()
        rng = random.Random(cs["seed"])
        n_, cls, net, labels, pk = C.setup(rng, cs["n"], cs["m"], cs["a"], cs["r"], chunks=cs["chunks"], Lc=cs["Lc"])
        nodes = [s for s in pk.S if s and C.v2(s) < pk.a]
        t = nodes[int(cs["tpick"] * len(nodes)) % len(nodes)]
        p, q, w, basis = cs["p"], cs["q"], cs["weights"], cs["basis"]
        mod = 1 << p
        ser, st = C.node_value(net, pk, labels, t, p, q, weights=w, method="cross", basis=basis)
        dp, _ = C.node_value(net, pk, labels, t, p, q, weights=w, method="dp", basis=basis)
        got = canon(st.get("records", []), mod)
        want = [tuple(x if not isinstance(x, list) else tuple(x) for x in row) for row in cs["records"]]
        ok_v = [int(x) % mod for x in ser] == [int(x) % mod for x in dp] == cs["value"]
        ok_r = got == want
        if not (ok_v and ok_r):
            bad += 1
            first = next((k for k in range(min(len(got), len(want))) if got[k] != want[k]), None)
            print("node seed", cs["seed"], "value ok", ok_v, "records ok", ok_r, "pieces", len(got), "want", len(want),
                  "first difference", None if first is None else (got[first], want[first]))
        else:
            print("node seed", cs["seed"], "ok:", len(got), "pieces,", round(time.time() - t0, 1), "s")
    for cs in data["e2e"]:
        rng = random.Random(cs["seed"])
        cls = C.random_formula(cs["n"], cs["m"], rng)
        calls = [0]
        orig = C.NodeExpansion.sum_cross

        def counting(self, *a, **k):
            calls[0] += 1
            return orig(self, *a, **k)
        C.NodeExpansion.sum_cross = counting
        try:
            got = C.count_e2e(cs["n"], cls, chunks=cs["chunks"], Lc=cs["Lc"], r=cs["r"], seed=cs["seed"], a=cs["a"],
                              u=cs["u"], basis="ldu", method="cross")
        finally:
            C.NodeExpansion.sum_cross = orig
        want = C.count_sat(cs["n"], cls)
        ok = int(got) == want and calls[0] == cs["cross_nodes"]
        print("e2e seed", cs["seed"], "count", got, "want", want, "sum_cross calls", calls[0], "want", cs["cross_nodes"])
        bad += 0 if ok else 1
    print("T24 OK" if bad == 0 else "T24 FAIL (%d mismatches)" % bad)


if __name__ == "__main__":
    main()
