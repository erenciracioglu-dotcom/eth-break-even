"""FROZEN referee, version 5. Usage: python referee.py <algo_dir> <out_json>

Runs three gates on the agent's copy of algo/ in separate processes:
  verify  - exact counts on hidden formulas vs a frozen brute force, construction checks
            (and the designated-variable construction when the code offers it)
  measure - the accounting facts of the frozen assembly hold on real designs
  score   - break-even sizes; every query is re-scored here with the clean assembly,
            so a patched assembly inside the agent's process cannot fake a score.
Parameters with designated, aggregated, units, chunk or modular only count when the verify gate
confirmed that construction or counting routine; otherwise those queries are scored as invalid.
Yes/no targets (A.DECISION_TARGETS, A.REPORT_TARGETS) accept every query; counting targets accept
only queries scored with decision=False (modular parameters are never valid for them).
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


def rescore(score_res, ok):
    """Clean re-scoring of every query; certified thresholds come only from these.
    ok: dict of verified flags (designated, aggregated, units, chunk, modular)."""
    queries = score_res.get("queries") or []
    tamper, certs, blocked = [], {}, 0
    for q in queries:
        dec = bool(q.get("decision", False))
        s = A.exponent_of(q["params"], q["lam"], decision=dec)
        e = s["e"] if s["valid"] else 9.0
        if abs(e - float(q["e"])) > 1e-9:
            tamper.append({"lam": q["lam"], "claimed": q["e"], "clean": e})
        if s.get("valid") and ((s.get("designated") and not ok["designated"])
                               or (s.get("aggregated") and not ok["aggregated"])
                               or (s.get("units") and not ok["units"])
                               or (s.get("chunk") and not ok["chunk"])
                               or (s.get("modular") and not ok["modular"])):
            e = 9.0
            blocked += 1
        q["e_clean"] = e
    targets = dict(A.TARGETS)
    targets.update(A.REPORT_TARGETS)
    yes_no = set(A.DECISION_TARGETS) | set(A.REPORT_TARGETS)
    for key, target in targets.items():
        good = [q for q in queries if q["e_clean"] < target and (key in yes_no or not q.get("decision", False))]
        if good:
            best = min(good, key=lambda q: q["lam"])
            certs[key] = {"lam": best["lam"], "params": best["params"], "e": best["e_clean"]}
        else:
            certs[key] = None
    ref = [q for q in queries if q["lam"] == 2224173]
    e_ref = ref[-1]["e_clean"] if ref else None
    return certs, tamper, e_ref, blocked


def main():
    algo, out = os.path.abspath(sys.argv[1]), sys.argv[2]
    outdir = os.path.dirname(os.path.abspath(out))
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    verify = run_gate(algo, "verify", outdir)
    measure = run_gate(algo, "measure", outdir) if verify.get("ok") else {"ok": False, "skipped": True}
    score = run_gate(algo, "score", outdir) if measure.get("ok") else {"ok": False, "skipped": True}
    designated_ok = bool(verify.get("designated_verified"))
    aggregated_ok = bool(verify.get("aggregated_verified"))
    units_ok = bool(verify.get("units_verified"))
    chunk_ok = bool(verify.get("chunk_verified"))
    modular_ok = bool(verify.get("modular_verified"))
    flags = dict(designated=designated_ok, aggregated=aggregated_ok, units=units_ok, chunk=chunk_ok, modular=modular_ok)
    certs, tamper, e_ref, blocked = (rescore(score, flags) if score.get("queries") else ({}, [], None, 0))
    summary = {k: v for k, v in score.items() if k != "queries"}
    if blocked and not certs.get("ppsz") and not summary.get("error"):
        summary["error"] = ("parameters with designated, aggregated, units, chunk or modular were used, but the "
                            "verify gate did not confirm that construction or counting routine "
                            f"({', '.join(k + '_verified=' + str(v) for k, v in flags.items())})")
    result = {
        "version": 5,
        "verify_ok": bool(verify.get("ok")),
        "measure_ok": bool(measure.get("ok")),
        "score_ok": bool(score.get("ok")) and not tamper and bool(certs.get("ppsz")),
        "designated_verified": designated_ok,
        "aggregated_verified": aggregated_ok,
        "units_verified": units_ok,
        "chunk_verified": chunk_ok,
        "modular_verified": modular_ok,
        "designated_blocked": blocked,
        "tamper": tamper[:5],
        "thresholds": {k: (v["lam"] if v else None) for k, v in certs.items()},
        "certificates": certs,
        "e_at_2224173": e_ref,
        "verify": verify,
        "measure": measure,
        "score_summary": summary,
        "seconds": round(time.time() - t0, 1),
    }
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=1, default=str)
    print(json.dumps({k: result[k] for k in ("verify_ok", "measure_ok", "score_ok", "designated_verified", "aggregated_verified", "units_verified",
                                             "chunk_verified", "modular_verified",
                                             "thresholds", "e_at_2224173", "seconds")}))


if __name__ == "__main__":
    main()
