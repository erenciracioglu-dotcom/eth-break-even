"""Round 12 (head): referee_appd version 8 = version 7 plus the analysis flag matchbal (appd-matchbal, GPT-6 Astra's
matching-based root label balance, head-checked). Copies referee_appd to work/head/r12/referee_v8 and patches the
assembly (flag, dependency check, root slack, docstring) and the version strings of referee.py and appd_gates.py.
Usage: python patch_v8.py"""
import os
import shutil

ROOT = "D:/ETH-loop"
SRC = os.path.join(ROOT, "referee_appd")
DST = os.path.join(ROOT, "work", "head", "r12", "referee_v8")
if os.path.exists(DST):
    shutil.rmtree(DST)
os.makedirs(DST)
for name in ("appd_assembly.py", "appd_gates.py", "appd_ref.py", "referee.py", "cases_appd.json"):
    data = open(os.path.join(SRC, name), "rb").read().replace(b"\r\n", b"\n")
    open(os.path.join(DST, name), "wb").write(data)


def patch(name, pairs):
    p = os.path.join(DST, name)
    s = open(p, encoding="utf-8").read()
    for old, new in pairs:
        assert s.count(old) == 1, (name, old[:80], s.count(old))
        s = s.replace(old, new)
    open(p, "w", encoding="utf-8", newline="\n").write(s)


DOC = '''matchbal (appd-matchbal, GPT-6 Astra, head-checked; version 8): the root label balance of a chunk network by
matchings. A chunk unit has at most 1 + 7(Lc - 1) variables, at most 21 Lc designated clauses with at most two
outside incidences each and at most 21 Lc outside incidences of its own variables, so its wire degree is
Delta <= 63 Lc. A greedy edge colouring splits the crossing wires (a loopless multigraph on the units) into
q <= 2 Delta - 1 <= 126 Lc - 1 matchings. Inside one matching different wires have disjoint ends, so their pair
labels are independent and uniform over the D pairs; the load of a (matching j, pair) is binomial with mean mu_j.
Budgets t_j = (delta/2)(mu_j + mu/q) add up to delta mu, and Bernstein (or a rational pessimistic estimator, which
derandomises the labels) bounds each tail by 2 exp(-delta^2 mu / (5 q)) for delta <= 1/4. The union bound over the
2 q D events then holds once M / D >= (8 q / delta^2) ln(8 q D): the root slack becomes
log2(8 q) + 2 log2(1/delta) + log2(ln(8 q) + ln D), linear in Lc where the two-stage chunk balancing was
quadratic. Port bounds follow by summing the pair loads. An analysis flag (no code); needs chunks.
'''
patch("appd_assembly.py", [
    ('"""FROZEN Appendix D assembly, referee_appd version 7 (2026-11-02, round 11: profile flag).',
     '"""FROZEN Appendix D assembly, referee_appd version 8 (2026-11-02, round 12: matchbal flag; version 7 added profile).'),
    ("a level is then min(eps_mm, eps_bc1 with bincross, PF_EPS[r]) for 4 <= r <= 16. Needs shellmm.\n",
     "a level is then min(eps_mm, eps_bc1 with bincross, PF_EPS[r]) for 4 <= r <= 16. Needs shellmm.\n" + DOC),
    ('         "bincross", "profile")\nCODE_FLAGS',
     '         "bincross", "profile", "matchbal")\nCODE_FLAGS'),
    ('    if fl["leafvars"] and not (fl["chunks"] and L == 1):',
     '    if fl["matchbal"] and not fl["chunks"]:\n'
     '        return _bad("matchbal needs chunks (it balances the wires of the chunk units)")\n'
     '    if fl["leafvars"] and not (fl["chunks"] and L == 1):'),
    ('        if d == 0 and fl["chunks"]:\n'
     '            slack = math.log2(8192.0 * Lc * Lc) + 2 * math.log2(1 / delta) + math.log2(math.log(172.0 * Lc) + lnD)\n',
     '        if d == 0 and fl["chunks"] and fl["matchbal"]:\n'
     '            q = 126.0 * Lc - 1.0                  # appd-matchbal: q <= 2 Delta - 1 matchings, Delta <= 63 Lc\n'
     '            slack = math.log2(8.0 * q) + 2 * math.log2(1 / delta) + math.log2(math.log(8.0 * q) + lnD)\n'
     '        elif d == 0 and fl["chunks"]:\n'
     '            slack = math.log2(8192.0 * Lc * Lc) + 2 * math.log2(1 / delta) + math.log2(math.log(172.0 * Lc) + lnD)\n'),
])
patch("referee.py", [("FROZEN referee for the Appendix D route, version 7 (round 11: profile-adaptive order gates).",
                      "FROZEN referee for the Appendix D route, version 8 (round 12: matchbal analysis flag in the "
                      "assembly; gates as in version 7).")])
patch("appd_gates.py", [("FROZEN gates of referee_appd version 7 (round 11).",
                         "FROZEN gates of referee_appd version 8 (round 12; gates unchanged from version 7).")])
print("referee_v8 staged at", DST)
