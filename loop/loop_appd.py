"""ETH autoresearch loop, Appendix D route (round 5). Reuses loop.py (agent runner, Strata handling, helpers);
its own workspace (algo_appd), frozen referee (referee_appd), backlog (state/backlog_appd.json) and best
(state/best_appd.json). The score is the break-even against 2^(n/3), the fair counting competitor.

Usage:
  python loop_appd.py          run until STOP_APPD is touched, the backlog is done, or errors pile up
  python loop_appd.py --test   one dry iteration without the agent (checks the plumbing, changes no state)
"""
import argparse
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import loop as L  # noqa: E402

ROOT = L.ROOT
ALGO = os.path.join(ROOT, "algo_appd")
REFEREE = os.path.join(ROOT, "referee_appd")
BEST = os.path.join(L.STATE, "best_appd.json")
BACKLOG = os.path.join(L.STATE, "backlog_appd.json")
LOGJ = os.path.join(L.STATE, "log_appd.jsonl")
SHAS = os.path.join(L.STATE, "referee_appd_sha.json")
STOP = os.path.join(ROOT, "STOP_APPD")
LOCK = os.path.join(ROOT, "loop_appd.lock")
FROZEN = ("appd_assembly.py", "appd_gates.py", "appd_ref.py", "referee.py", "cases_appd.json")
MIN_GAIN = 0.0005
L.AGENT_WALL_S = 55 * 60       # Appendix D tasks are larger than the main-route ones
KEYS = ("third", "occ3sat", "ppsz", "sharp3sat", "brute")

# the shared helpers read these module globals
L.REFEREE = REFEREE
L.BACKLOG = BACKLOG


def tracked_files():
    out = L.git("ls-files", "algo_appd")
    return [line[len("algo_appd/"):] for line in out.splitlines() if line.startswith("algo_appd/")]


def repo_guard():
    dirty = L.git("status", "--porcelain", "--untracked-files=no", "--", "algo_appd", "referee_appd")
    if dirty:
        L.log("WARNING: algo_appd/referee_appd were modified outside the loop; restoring: " + dirty.replace("\n", " | "))
        L.git("checkout", "--", "algo_appd", "referee_appd")
    ref = L.jload(SHAS, {})
    for name in FROZEN:
        if L.sha(os.path.join(REFEREE, name)) != ref.get(name):
            raise RuntimeError(f"referee_appd file {name} does not match state/referee_appd_sha.json")
    return dirty


def flags_text(best):
    f = best.get("verified_flags") or {}
    on = [k for k, v in f.items() if v]
    return ", ".join(on) if on else "none"


def build_prompt(task, best, attempt):
    feedback = ""
    if attempt > 1 and task.get("last_feedback"):
        feedback = ("\nYour previous attempt was rejected: " + task["last_feedback"] +
                    "\nDo not repeat the same change; try a different approach.\n")
    return L.fill(L.read_text(os.path.join(L.LOOPDIR, "task_template_appd.txt")),
                  id=task["id"], attempt=attempt, max_attempts=task.get("max_attempts", 3), title=task["title"],
                  files=", ".join(task["files"]), spec=task["spec"], third=f"{best['third']:,}",
                  flags=flags_text(best), feedback=feedback, stop=task["stop"], steps=task.get("steps", L.AGENT_STEPS))


def problem(ref):
    if not ref or ref.get("error"):
        return (ref or {}).get("error", "referee produced no result")
    if not ref.get("verify_ok"):
        bad = ref.get("verify_failures") or []
        return "verify gate failed on " + "; ".join(str({k: v for k, v in i.items() if k != "ok"})[:300]
                                                     for i in bad[:3])
    s = ref.get("score_summary") or {}
    if s.get("error"):
        return "score gate crashed: " + str(s["error"])[:600]
    if (s.get("choose_params_avg_s") or 0) >= 3.0:
        return f"choose_params too slow: {s['choose_params_avg_s']:.2f} s per call (limit 3 s)"
    if not (ref.get("thresholds") or {}).get("third"):
        bad = ref.get("verify_failures") or []
        why = "; ".join(str({k: v for k, v in i.items() if k != "ok"})[:300] for i in bad[:3])
        return ("score gate: the 2^(n/3) target was never reached (parameters that use an unverified code flag "
                "are invalid)" + ("; failing checks: " + why if why else ""))
    return None


def decide(task, ref, best):
    p = problem(ref)
    if p:
        return False, p
    new, old = ref["thresholds"]["third"], best["third"]
    detail = f"third threshold {old:,} -> {new:,}"
    need = task.get("needs_flags") or []
    missing = [f for f in need if not (ref.get("verified_flags") or {}).get(f)]
    if missing:
        bad = ref.get("verify_failures") or []
        why = "; ".join(str({k: v for k, v in i.items() if k != "ok"})[:300] for i in bad[:3])
        return False, ("code flags not verified: " + ", ".join(missing) + "; " + detail
                       + ("; failing checks: " + why if why else ""))
    if task.get("kind", "score") == "score":
        if new < old * (1 - MIN_GAIN):
            return True, detail + f" ({100 * (1 - new / old):.2f}% lower threshold)"
        return False, "no gain: " + detail
    if new <= old * (1 + MIN_GAIN):
        return True, "gates pass, no regression: " + detail
    return False, "regression: " + detail


def accept(it, task, cand, ref, best):
    L.copy_files(cand, ALGO, task["files"])
    nb = dict(best)
    nb.update({k: (ref["thresholds"] or {}).get(k) for k in KEYS})
    nb.update(verified_flags=ref.get("verified_flags"), params_third=((ref.get("certificates") or {}).get("third") or {}).get("params"),
              iter=it, task=task["id"], updated=L.now())
    L.jsave(BEST, nb)
    L.git("add", "algo_appd", "state", check=True)
    msg = f"iter {it:03d} {task['id']} (Appendix D): third threshold {best['third']:,} -> {nb['third']:,}"
    L.git("commit", "-m", msg + L.TRAILER, check=True)
    return nb, L.git("rev-parse", "--short", "HEAD")


def pick_task(backlog):
    return L.pick_task(backlog)


def iteration(test=False):
    backlog = L.jload(BACKLOG, [])
    task = pick_task(backlog)
    if task is None:
        return "idle"
    best = L.jload(BEST)
    guard = repo_guard()
    if not test and not L.ensure_strata():
        return "strata_down"
    prefix = "test_appd" if test else "iter_appd"
    it = L.next_iter_no(prefix)
    itdir = os.path.join(L.WORK, f"{prefix}_{it:03d}")
    ws, cand = os.path.join(itdir, "ws"), os.path.join(itdir, "candidate")
    names = tracked_files()
    L.copy_files(ALGO, ws, names)
    attempt = int(task.get("attempts", 0)) + 1
    with open(os.path.join(ws, "AGENTS.md"), "w", encoding="utf-8") as fh:
        fh.write(L.read_text(os.path.join(L.LOOPDIR, "agents_template_appd.md")))
    prompt = build_prompt(task, best, attempt)
    with open(os.path.join(itdir, "prompt.txt"), "w", encoding="utf-8") as fh:
        fh.write(prompt)
    record = {"iter": it, "test": test, "task": task["id"], "attempt": attempt, "start": L.now(),
              "best_before": best["third"], "guard": guard or None}
    L.log(f"{prefix} {it:03d}: task {task['id']} attempt {attempt} - {task['title']}")
    if test:
        agent = {"rc": 0, "seconds": 0, "killed": None, "requests": 0, "tail": "(test: no agent)"}
    else:
        L.update_task(task["id"], status="running")
        L.set_status(state="agent_appd", iter=it, task=task["id"], attempt=attempt, since=L.now())
        agent = L.run_agent(itdir, ws, prompt, task.get("steps", L.AGENT_STEPS))
        L.log(f"agent finished: rc={agent['rc']} {agent['seconds']} s, {agent['requests']} requests"
              + (f", killed: {agent['killed']}" if agent["killed"] else ""))
    record["agent"] = {k: v for k, v in agent.items() if k != "tail"}
    record["agent_summary"] = agent["tail"][-800:]
    infra = (not test) and agent["rc"] != 0 and not agent["killed"] and not L.strata_up()
    L.copy_files(ALGO, cand, names)
    changed = [n for n in names if os.path.exists(os.path.join(ws, n))
               and L.sha(os.path.join(ws, n)) != L.sha(os.path.join(ALGO, n))]
    allowed = [n for n in changed if n in task["files"]]
    big = [n for n in allowed if os.path.getsize(os.path.join(ws, n)) > L.MAX_FILE_BYTES]
    L.copy_files(ws, cand, [n for n in allowed if n not in big])
    record.update(changed=allowed, ignored_edits=[n for n in changed if n not in task["files"]],
                  new_files=[f for f in L.ws_files(ws) if f not in names and f != "AGENTS.md"][:30])
    ref = None
    if infra:
        accepted, reason = False, "infrastructure: model server went down during the run (attempt not counted)"
    elif big:
        accepted, reason = False, "file above 1 MB: " + ", ".join(big)
    elif not allowed and not test:
        accepted, reason = False, "no change to " + ", ".join(task["files"])
    else:
        repo_guard()
        L.set_status(state="referee_appd", iter=it, task=task["id"], attempt=attempt, since=L.now())
        ref = L.run_referee(itdir, cand)
        accepted, reason = decide(task, ref, best)
    if ref:
        record["referee"] = {k: ref.get(k) for k in ("verify_ok", "implemented", "verified_flags", "thresholds",
                                                      "verify", "score_summary", "seconds", "error")}
    record.update(end=L.now(), reason=reason)
    if test:
        record["decision"] = "test-accept" if accepted else "test-reject"
        L.append_jsonl(os.path.join(L.STATE, "test_log.jsonl"), record)
        L.log(f"test iteration (appd): {record['decision']} - {reason}")
        return "test"
    if accepted:
        nb, commit = accept(it, task, cand, ref, best)
        record.update(decision="accepted", commit=commit, best_after=nb["third"])
        L.update_task(task["id"], status="done", attempts=attempt, result=reason, done_iter=it)
        L.bulletin_add(f"iter_appd {it:03d} {task['id']} attempt {attempt}: ACCEPTED - {reason}")
    elif infra:
        record["decision"] = "infra"
        L.update_task(task["id"], status="open")
        L.bulletin_add(f"iter_appd {it:03d} {task['id']}: infrastructure failure, will retry - {reason}")
    else:
        maxa = int(task.get("max_attempts", 3))
        record["decision"] = "rejected"
        feedback = reason + ". Agent's own summary: " + agent["tail"][-500:].replace("\n", " ")
        L.update_task(task["id"], status="failed" if attempt >= maxa else "open", attempts=attempt,
                      last_feedback=feedback[:1500], result=reason)
        L.bulletin_add(f"iter_appd {it:03d} {task['id']} attempt {attempt}/{maxa}: REJECTED - {reason}")
    L.append_jsonl(LOGJ, record)
    L.log(f"iter_appd {it:03d}: {record['decision']} - {reason}")
    return record["decision"]


def take_lock():
    old = L.jload(LOCK, {})
    if old and old.get("pid") and old.get("pid") != os.getpid() and L.pid_alive(old["pid"]):
        return False
    L.jsave(LOCK, {"pid": os.getpid(), "since": L.now()})
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    args = ap.parse_args()
    os.makedirs(L.LOGS, exist_ok=True)
    os.makedirs(L.WORK, exist_ok=True)
    if not take_lock():
        L.log("another Appendix D loop holds the lock; exiting")
        return
    L.log(f"Appendix D loop started ({'test' if args.test else 'loop'}), pid {os.getpid()}")
    errors = 0
    try:
        while not os.path.exists(STOP):
            try:
                res = iteration(test=args.test)
                errors = 0
            except Exception:  # noqa: BLE001
                errors += 1
                res = "error"
                L.log("ERROR in Appendix D iteration:\n" + traceback.format_exc())
                if errors >= L.MAX_ERRORS:
                    L.log("too many errors in a row; stopping")
                    break
            if args.test or res == "idle":
                break
            if res in ("strata_down", "error"):
                time.sleep(120)
        if os.path.exists(STOP):
            L.log("STOP_APPD requested")
    finally:
        if not args.test:
            L.stop_strata_if_ours()
            L.set_status(state="stopped", since=L.now())
        L.jsave(LOCK, {"pid": None, "released": L.now()})
        L.log("Appendix D loop exited")


if __name__ == "__main__":
    main()
