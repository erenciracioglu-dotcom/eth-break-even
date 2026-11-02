"""Round 12 (head solution, second pass): a fine analytic re-rank of the best coarse candidates over delta and Lc
neighbours, and a wider coordinate refinement. Patches work/head/r12/hs12/appd_params.py in place."""
P = "D:/ETH-loop/work/head/r12/hs12/appd_params.py"
s = open(P, encoding="utf-8").read()


def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:80], s.count(old))
    s = s.replace(old, new)


rep("REFINE_TOP = 3             # round 12: distinct best builds refined by coordinate search\n"
    "REFINE_BUILDS = 160        # round 12: cap on build() calls in the refinement\n",
    "REFINE_TOP = 4             # round 12: distinct best builds refined by coordinate search\n"
    "REFINE_BUILDS = 220        # round 12: cap on build() calls in the refinement\n"
    "FINE_TOP = 60              # round 12: coarse candidates re-ranked on finer delta and Lc neighbours\n"
    "FINE_DELTA = (0.6, 0.8, 1.25, 1.6)\n"
    "FINE_LC = (0.5, 2.0)\n")
rep("            tries = [(r, delta * f, Lc, kap) for f in (0.7, 0.85, 1.15, 1.4)]\n",
    "            tries = [(r, delta * f, Lc, kap) for f in (0.7, 0.85, 1.15, 1.4, 0.5, 2.0)]\n")
rep("            if kap is not None:\n                tries += [(r, delta, Lc, kap - 1), (r, delta, Lc, kap + 1)]\n",
    "            if kap is not None:\n                tries += [(r, delta, Lc, kap + dk) for dk in (-1, 1, -2, 2)]\n")
rep("""    cand.sort(key=lambda x: x[0])
    built = []
""", """    cand.sort(key=lambda x: x[0])
    fine = []
    for _, shape, delta, Lc in cand[:FINE_TOP]:          # round 12: finer delta and Lc around the best shapes
        lv = flags["chunks"] and len(shape) == 1
        for d2 in [delta * f for f in FINE_DELTA] + [delta]:
            for L2 in [int(Lc * g) for g in FINE_LC] + [Lc]:
                if (d2, L2) == (delta, Lc) or not (1e-4 <= d2 <= 0.25 and 2 <= L2 <= (1 << 20)):
                    continue
                e, _ = rank_e(list(shape), d2, L2, lam, lv, cK, flags["shellmm"], flags["bincross"],
                              flags["profile"])
                fine.append((e, shape, d2, L2))
    cand = sorted(cand + fine, key=lambda x: x[0])
    built = []
""")
open(P, "w", encoding="utf-8", newline="\n").write(s)
print("patched", P)
