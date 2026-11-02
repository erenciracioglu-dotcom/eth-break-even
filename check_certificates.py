"""Re-check every round-12 break-even cell in the frozen cost model (a few seconds, standard library only).

For each competitor the certificate in work/head/r12/round12_witnesses.json must be valid at its log2 n,
with a cost exponent e below the competitor's exponent, and the same parameters must fail one size below,
where a size condition (label balance or padding) no longer holds. Run from the repository root:

    python check_certificates.py
"""
import hashlib
import importlib.util
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(ROOT, "referee_appd", "appd_assembly.py")
CERTS = os.path.join(ROOT, "work", "head", "r12", "round12_witnesses.json")
SHA_V8 = "67262791697b7f1ae477a7e303d99635b02656a5e84cc67787653b39ffe60514"
ORDER = [("brute", "2^n (brute force)"), ("sharp3sat", "1.637^n"), ("third", "2^(n/3)"),
         ("ppsz", "1.307^n (PPSZ)"), ("occ3sat", "1.1199^n")]


def main():
    sha = hashlib.sha256(open(MODEL, "rb").read()).hexdigest()
    print("cost model referee_appd/appd_assembly.py, SHA-256", sha,
          "(version 8)" if sha == SHA_V8 else "(NOT the version-8 file)")
    spec = importlib.util.spec_from_file_location("appd_assembly", MODEL)
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    data = json.load(open(CERTS, encoding="utf-8"))
    certs = data.get("certificates", data)
    ok = sha == SHA_V8
    print(f"{'competitor':<20}{'log2 n':>8}{'n about':>10}{'e':>11}{'target':>11}  valid  fails one below")
    for key, label in ORDER:
        c = certs[key]
        lam = c["lam"]
        here = A.exponent_of(c["params"], lam)
        below = A.exponent_of(c["params"], lam - 1)
        target = A.TARGETS[key]
        good = bool(here["valid"]) and here["e"] < target and not below["valid"]
        ok = ok and good
        print(f"{label:<20}{lam:>8,}{'10^' + str(round(lam * math.log10(2))):>10}{here['e']:>11.6f}{target:>11.6f}"
              f"  {str(bool(here['valid'])):<6} {str(not below['valid'])}")
    print("ALL CERTIFICATES OK" if ok else "CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
