"""Round 12 (head): install referee_appd version 8, the public v8 model and the head's round-12 search into the repo,
record the referee's validation run as the new baseline, and append the round-12 ledger entry.
Usage: python install_r12.py"""
import hashlib
import json
import os
import shutil

ROOT = "D:/ETH-loop"
R12 = os.path.join(ROOT, "work", "head", "r12")
REF = os.path.join(ROOT, "referee_appd")
ALGO = os.path.join(ROOT, "algo_appd")
FROZEN = ("appd_assembly.py", "appd_gates.py", "appd_ref.py", "referee.py", "cases_appd.json")


def lf_copy(src, dst):
    data = open(src, "rb").read().replace(b"\r\n", b"\n")
    open(dst, "wb").write(data)


for name in FROZEN:
    lf_copy(os.path.join(R12, "referee_v8", name), os.path.join(REF, name))
lf_copy(os.path.join(R12, "referee_v8", "appd_assembly.py"), os.path.join(ALGO, "appd_assembly_public.py"))
lf_copy(os.path.join(R12, "hs12", "appd_params.py"), os.path.join(ALGO, "appd_params.py"))
# frozen fingerprints
sha_p = os.path.join(ROOT, "state", "referee_appd_sha.json")
old = json.load(open(sha_p, encoding="utf-8"))
new = {name: hashlib.sha256(open(os.path.join(REF, name), "rb").read()).hexdigest() for name in FROZEN}
assert set(old) == set(new), (sorted(old), sorted(new))
json.dump(new, open(sha_p, "w", encoding="utf-8"), indent=1)
# baseline = the referee v8 run on the head solution
good = json.load(open(os.path.join(R12, "val", "good", "referee.json"), encoding="utf-8"))
assert good["verify_ok"] and all(good["verified_flags"].values())
best_p = os.path.join(ROOT, "state", "best_appd.json")
best = json.load(open(best_p, encoding="utf-8"))
for k, v in good["thresholds"].items():
    best[k] = v
for k, c in good["certificates"].items():
    best["params_" + k] = c["params"]
best["verified_flags"] = good["verified_flags"]
best["referee_version"] = 8
best["note_r12"] = ("round 12 (head): referee v8 on the head search; the page reports the smallest dictionaries of three "
                    "searches, see work/head/r12/round12_witnesses.json and claims.md")
json.dump(best, open(best_p, "w", encoding="utf-8"), indent=1)
# ledger
cur = json.load(open(os.path.join(R12, "val", "current", "referee.json"), encoding="utf-8"))
bad = json.load(open(os.path.join(R12, "val", "bad", "referee.json"), encoding="utf-8"))
d426 = json.load(open(os.path.join(R12, "dens_v8_4.26.json"), encoding="utf-8"))["thresholds"]
d5 = json.load(open(os.path.join(R12, "dens_v8_5.0.json"), encoding="utf-8"))["thresholds"]


def fmt(t):
    return " / ".join(f"{t[k]:,}" if t[k] is not None else "none" for k in ("third", "occ3sat", "ppsz", "sharp3sat", "brute"))


vline = (f"Validation (2 November, 20:48-20:56): the head solution (work/head/r12/hs12) verifies all 11 code flags and "
         f"scores {fmt(good['thresholds'])} (third / occ3sat / ppsz / sharp3sat / brute) at "
         f"{good['score_summary']['choose_params_avg_s']:.2f} s per call, matchbal on in every certificate; the untouched "
         f"round-11 code keeps {fmt(cur['thresholds'])}; a variant that sets matchbal without chunks keeps its code flags "
         f"and scores {fmt(bad['thresholds'])} (every dictionary rejected).")
dline = (f"m = 4.26n gives {d426['sharp3sat']:,} / {d426['ppsz']:,} / {d426['brute']:,} against 1.637^n / 1.307^n / 2^n "
         f"(the withdrawn round-11 entries were 6,063 / 20,415 / 3,089); m = 5n gives {d5['sharp3sat']:,} / {d5['ppsz']:,} / "
         f"{d5['brute']:,} (withdrawn: 8,660 / 29,219 / 4,401), n about 10^{round(d5['brute'] * 0.30103):,} against 2^n.")
add = open(os.path.join(R12, "claims_r12.md"), encoding="utf-8").read().replace("\r\n", "\n")
add = add.replace("VALIDATION_V8", vline).replace("DENSITY_LINE", dline)
claims = os.path.join(ROOT, "state", "claims.md")
text = open(claims, encoding="utf-8").read()
assert "## Round 12 (2 November, evening)" not in text
open(claims, "w", encoding="utf-8", newline="\n").write(text.rstrip("\n") + "\n" + add)
print(vline)
print(dline)
