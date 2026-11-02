"""Round 12 (head solution): appd_params.py with the matchbal flag (appd-matchbal) in the analytic replica and a
coordinate refinement of the best builds over delta, Lc, every r_d and the guard kappa. Patches
work/head/r12/hs12/appd_params.py in place. Usage: python patch_params_v8.py"""
P = "D:/ETH-loop/work/head/r12/hs12/appd_params.py"
s = open(P, encoding="utf-8").read().replace("\r\n", "\n")


def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:80], s.count(old))
    s = s.replace(old, new)


rep('_KAPPA = None              # set by choose_params (None = the paper\'s Cauchy deformation)\n',
    '_KAPPA = None              # set by choose_params (None = the paper\'s Cauchy deformation)\n'
    '_MATCHBAL = False          # set by choose_params: root label balance by matchings (appd-matchbal, round 12)\n'
    'REFINE_TOP = 3             # round 12: distinct best builds refined by coordinate search\n'
    'REFINE_BUILDS = 160        # round 12: cap on build() calls in the refinement\n')
rep('''    if d == 0:
        base = math.log2(8192.0 * Lc * Lc) + inv
        c1 = math.log(172.0 * Lc)
''', '''    if d == 0 and _MATCHBAL:                 # appd-matchbal: q = 126 Lc - 1 matchings at the root
        q = 126.0 * Lc - 1.0
        base = math.log2(8.0 * q) + inv
        c1 = math.log(8.0 * q)
    elif d == 0:
        base = math.log2(8192.0 * Lc * Lc) + inv
        c1 = math.log(172.0 * Lc)
''')
rep('''def choose_params(lam, target="third"):
    global _KAPPA
    flags = dict(perclass=True, labelfree=True, memo=False, modedp=False, chunks=False, leafvars=False,
                 newton=False, ldu=False, shellmm=False, guarded=False, bincross=False, profile=False)
    if USE_CODE_FLAGS:
        flags.update(memo=True, modedp=True, chunks=True, newton=True, ldu=True, shellmm=True, guarded=True,
                     bincross=True, profile=True)
    _KAPPA = GUARD_KAPPA if flags["guarded"] else None
''', '''def _refine(cands, flags, lam):
    """Round 12: coordinate search from the best builds. Moves delta (x 0.7, 0.85, 1.15, 1.4), Lc (x 1/2, 2), every
    r_d (+-1) and, with guards, kappa (+-1); every trial is a full build() (u by bisection, a = min_a) scored by
    the frozen model at lam; first improvement wins; at most REFINE_BUILDS builds in total."""
    global _KAPPA
    budget = [REFINE_BUILDS]

    def trial(r, delta, Lc, kap):
        global _KAPPA
        if budget[0] <= 0:
            return None
        if not (1e-4 <= delta <= 0.25 and 1 <= Lc <= (1 << 20) and all(2 <= x <= 4096 for x in r)):
            return None
        if kap is not None and not 2 <= kap <= 30:
            return None
        budget[0] -= 1
        _KAPPA = kap
        fl = dict(flags, leafvars=flags["chunks"] and len(r) == 1)
        res = build(list(r), delta, Lc, fl, lam)
        if res is None:
            return None
        s = A.exponent_of(res[1], lam)
        return (s["e"], res[1]) if s["valid"] else None

    best = None
    for e0, p0 in cands[:REFINE_TOP]:
        cur = (e0, p0)
        moved = True
        while moved and budget[0] > 0:
            moved = False
            r, delta, Lc = list(cur[1]["r"]), cur[1]["delta"], cur[1]["Lc"]
            kap = cur[1]["kappa"][0] if cur[1].get("kappa") else None
            tries = [(r, delta * f, Lc, kap) for f in (0.7, 0.85, 1.15, 1.4)]
            tries += [(r, delta, Lc // 2, kap), (r, delta, Lc * 2, kap)]
            for d in range(len(r)):
                for dr in (-1, 1):
                    rr = list(r)
                    rr[d] += dr
                    tries.append((rr, delta, Lc, kap))
            if kap is not None:
                tries += [(r, delta, Lc, kap - 1), (r, delta, Lc, kap + 1)]
            for t in tries:
                res = trial(*t)
                if res is not None and res[0] < cur[0] - 1e-15:
                    cur, moved = res, True
                    break
        if best is None or cur[0] < best[0]:
            best = cur
    return best


def choose_params(lam, target="third"):
    global _KAPPA, _MATCHBAL
    flags = dict(perclass=True, labelfree=True, memo=False, modedp=False, chunks=False, leafvars=False,
                 newton=False, ldu=False, shellmm=False, guarded=False, bincross=False, profile=False,
                 matchbal=False)
    if USE_CODE_FLAGS:
        flags.update(memo=True, modedp=True, chunks=True, newton=True, ldu=True, shellmm=True, guarded=True,
                     bincross=True, profile=True, matchbal=True)
    _KAPPA = GUARD_KAPPA if flags["guarded"] else None
    _MATCHBAL = flags["matchbal"]
''')
rep('''    best = None
    seen = set()
    for _, shape, delta, Lc in cand[:TOPK]:
        key = (shape, delta, Lc)
        if key in seen:
            continue
        seen.add(key)
        fl = dict(flags, leafvars=flags["chunks"] and len(shape) == 1)
        p = build(list(shape), delta, Lc, fl, lam)
        if p is None:
            continue
        s = A.exponent_of(p[1], lam)
        if s["valid"] and (best is None or s["e"] < best[0]):
            best = (s["e"], p[1])
    if best:
        return best[1]
''', '''    built = []
    seen = set()
    for _, shape, delta, Lc in cand[:TOPK]:
        key = (shape, delta, Lc)
        if key in seen:
            continue
        seen.add(key)
        fl = dict(flags, leafvars=flags["chunks"] and len(shape) == 1)
        p = build(list(shape), delta, Lc, fl, lam)
        if p is None:
            continue
        s = A.exponent_of(p[1], lam)
        if s["valid"]:
            built.append((s["e"], p[1]))
    built.sort(key=lambda x: x[0])
    best = _refine(built, flags, lam) if built else None
    _KAPPA = GUARD_KAPPA if flags["guarded"] else None
    if best:
        return best[1]
''')
rep('Search (attempt 3, round 9 / T22):', 'Search (round 12, head): matchbal in bal_u and a coordinate refinement of the best builds (_refine).\n\n'
    'Search (attempt 3, round 9 / T22):')
open(P, "w", encoding="utf-8", newline="\n").write(s)
print("patched", P)
