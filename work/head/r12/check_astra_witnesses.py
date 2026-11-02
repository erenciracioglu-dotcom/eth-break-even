"""Head check of GPT-6 Astra's hostile review (2 November): her frozen-model witnesses, rebuilt from the review's
tables on top of the round-11 certificates, evaluated with the unchanged frozen v7 assembly; and bincross off."""
import copy, json, math, sys
sys.path.insert(0, "D:/ETH-loop/referee_appd")
import appd_assembly as A
cert = json.load(open("D:/ETH-loop/work/head/r11/final_referee_r11.json"))["certificates"]
W = {  # target: (lambda, r, u, a, delta, Lc, e reported by Astra)
    "brute": (100, [7], [11], [29], 0.0362, 32, 0.999534479794686),
    "sharp3sat": (157, [8], [19], [38], 0.0200, 64, 0.710254458410584),
    "third": (532, [8, 9, 11], [79, 63, 38], [101, 166, 206], 0.0029, 512, 0.333319277422392),
    "ppsz": (411, [8, 9, 11], [59, 47, 29], [81, 130, 161], 0.0021, 512, 0.386265410589750),
    "occ3sat": (2033, [8, 9, 16], [323, 257, 173], [349, 608, 783], 0.0004, 1048576, 0.163366065628536),
}
print("targets", A.TARGETS)
out = {}
for k, (lam, r, u, a, delta, Lc, e_ast) in W.items():
    p = copy.deepcopy(cert[k]["params"])
    p.update(r=r, u=u, a=a, delta=delta, Lc=Lc, kappa=[10] * len(r), levels=len(r))
    v = A.exponent_of(p, lam)
    v1 = A.exponent_of(p, lam - 1)
    tgt = A.TARGETS[k]
    ok = v["valid"] and v["e"] < tgt
    out[k] = dict(lam=lam, valid=v["valid"], e=v.get("e"), target=tgt, ok=ok, e_matches=abs((v.get("e") or 9) - e_ast) < 1e-12,
                  valid_at_lam_minus_1=v1["valid"], reason=v.get("reason") or v1.get("reason"))
    print(k, out[k])
print("--- bincross off on the round-11 certificates")
for k, c in cert.items():
    p = copy.deepcopy(c["params"]); p["bincross"] = False
    v0 = A.exponent_of(c["params"], c["lam"]); v = A.exponent_of(p, c["lam"])
    print(k, c["lam"], "e with bincross", v0["e"], "without", v.get("e"), "same" if abs(v0["e"] - v.get("e", 9)) < 1e-15 else "DIFFERENT")
json.dump(out, open("D:/ETH-loop/work/head/r12/astra_witnesses_check.json", "w"), indent=1)
