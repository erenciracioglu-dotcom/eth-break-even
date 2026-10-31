"""FROZEN referee. Usage: python referee.py <algo_dir> <out_json>

Runs three gates on the agent's copy of algo/ in separate processes:
  verify  - exact counts on hidden formulas vs a frozen brute force, construction checks
  measure - the accounting facts of the frozen assembly hold on real designs
  score   - break-even sizes; every query is re-scored here with the clean assembly,
            so a patched assembly inside the agent's process cannot fake a score.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import assembly as A  # noqa: E402

TIMEOUTS = {"verify": 1200, "measure": 600, "score": 1200}


def run_gate(algo, gate, outdir):
    out = os.path.join(outdir, f"gate_{gate}.json")
    if os.path.exists(out):
        os.remove(out)
    t0 = time.time()
    try:
        proc = subprocess.run([sys.executable, "-I", os.path.join(HERE, "gates.py"), algo, gate, out],
                              timeout=TIMEOUTS[gate], capture_output=True, text=True,
                              encoding="utf-8", errors="replace")
        tail = (proc.stdout + proc.stderr)[-1500:]
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"{gate} gate timed out after {TIMEOUTS[gate]} s"}
    if not os.path.exists(out):
        return {"ok": False, "error": f"{gate} gate produced no result", "tail": tail}
    with open(out, encoding="utf-8") as fh:
        res = json.load(fh)
    res.setdefault("seconds", round(time.time() - t0, 1))
    return res


def rescore(score_res):
    """Clean re-scoring of every query; certified thresholds come only from these."""
    queries = score_res.get("queries") or []
    tamper, certs = [], {}
    for q in queries:
        s = A.exponent_of(q["params"], q["lam"])
        e = s["e"] if s["valid"] else 9.0
        if abs(e - float(q["e"])) > 1e-9:
            tamper.append({"lam": q["lam"], "claimed": q["e"], "clean": e})
        q["e_clean"] = e
    for key, target in A.TARGETS.items():
        good = [q for q in queries if q["e_clean"] < target]
        if good:
            best = min(good, key=lambda q: q["lam"])
            certs[key] = {"lam": best["lam"], "params": best["params"], "e": best["e_clean"]}
        else:
            certs[key] = None
    ref = [q for q in queries if q["lam"] == 2224173]
    e_ref = ref[-1]["e_clean"] if ref else None
    return certs, tamper, e_ref


def main():
    algo, out = os.path.abspath(sys.argv[1]), sys.argv[2]
    outdir = os.path.dirname(os.path.abspath(out))
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    verify = run_gate(algo, "verify", outdir)
    measure = run_gate(algo, "measure", outdir) if verify.get("ok") else {"ok": False, "skipped": True}
    score = run_gate(algo, "score", outdir) if measure.get("ok") else {"ok": False, "skipped": True}
    certs, tamper, e_ref = rescore(score) if score.get("queries") else ({}, [], None)
    result = {
        "verify_ok": bool(verify.get("ok")),
        "measure_ok": bool(measure.get("ok")),
        "score_ok": bool(score.get("ok")) and not tamper and bool(certs.get("ppsz")),
        "tamper": tamper[:5],
        "thresholds": {k: (v["lam"] if v else None) for k, v in certs.items()},
        "certificates": certs,
        "e_at_2224173": e_ref,
        "verify": {k: v for k, v in verify.items() if k != "trace"} | ({"trace": verify["trace"]} if "trace" in verify else {}),
        "measure": measure,
        "score_summary": {k: v for k, v in score.items() if k != "queries"},
        "seconds": round(time.time() - t0, 1),
    }
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, default=str)
    print(json.dumps({k: result[k] for k in ("verify_ok", "measure_ok", "score_ok", "thresholds", "e_at_2224173", "seconds")}))


if __name__ == "__main__":
    main()
