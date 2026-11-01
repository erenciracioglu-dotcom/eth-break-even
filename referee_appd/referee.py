"""FROZEN referee for the Appendix D route, version 2 (round 6: Newton basis gates). Usage: python referee.py <candidate_dir> <out_json>
Runs the verify gate, then the score gate with the flags verify confirmed (separate processes)."""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TIMEOUTS = {"verify": 1800, "score": 1800}


def run_gate(cand, gate, outdir, extra=None):
    out = os.path.join(outdir, f"gate_{gate}.json")
    if os.path.exists(out):
        os.remove(out)
    cmd = [sys.executable, "-I", os.path.join(HERE, "appd_gates.py"), cand, gate, out] + ([extra] if extra else [])
    try:
        proc = subprocess.run(cmd, timeout=TIMEOUTS[gate], capture_output=True, text=True, encoding="utf-8",
                              errors="replace")
        tail = (proc.stdout + proc.stderr)[-1500:]
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"{gate} gate timed out after {TIMEOUTS[gate]} s"}
    if not os.path.exists(out):
        return {"ok": False, "error": f"{gate} gate produced no result", "tail": tail}
    return json.load(open(out, encoding="utf-8"))


def main():
    cand, outp = sys.argv[1], sys.argv[2]
    outdir = os.path.dirname(os.path.abspath(outp))
    t0 = time.time()
    ver = run_gate(cand, "verify", outdir)
    flags = ver.get("flags") or {}
    sc = run_gate(cand, "score", outdir, json.dumps(flags))
    res = {"verify_ok": bool(ver.get("ok")), "implemented": ver.get("implemented"), "verified_flags": flags,
           "verify": {k: ver.get(k) for k in ("counts", "error", "seconds")},
           "verify_failures": [i for i in ver.get("items", []) if not i.get("ok")][:5],
           "thresholds": sc.get("thresholds"), "certificates": sc.get("certificates"),
           "score_summary": {k: sc.get(k) for k in ("choose_params_avg_s", "calls", "error", "seconds")},
           "seconds": round(time.time() - t0, 1)}
    with open(outp, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("verify_ok", "implemented", "verified_flags", "thresholds")}))


if __name__ == "__main__":
    main()
