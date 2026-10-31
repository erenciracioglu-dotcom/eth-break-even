"""Finishes an iteration whose agent completed but whose bookkeeping crashed (iter 009, T9: a file named
'nul' in the workspace broke ws_files before the referee ran). Same steps as loop.iteration() after
the agent run, with the agent result taken from the iteration folder.
Usage: python rescue_iter.py <iter_no> <task_id> <agent_seconds> <agent_requests>"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import loop as L  # noqa: E402

it, task_id, secs, reqs = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
itdir = os.path.join(L.WORK, f"iter_{it:03d}")
ws, cand = os.path.join(itdir, "ws"), os.path.join(itdir, "candidate")
backlog = L.jload(L.BACKLOG, [])
task = next(t for t in backlog if t["id"] == task_id)
best = L.jload(L.BEST)
attempt = 1
names = L.tracked_algo_files()
agent = {"rc": 0, "seconds": secs, "killed": None, "requests": reqs,
         "tail": L.tail_text(os.path.join(itdir, "eragent.log"))}
record = {"iter": it, "test": False, "task": task_id, "attempt": attempt, "start": "(see loop log)",
          "best_before": best["ppsz"], "guard": None,
          "rescued": "bookkeeping crashed on a file named 'nul'; finished by loop/rescue_iter.py"}
record["agent"] = {k: v for k, v in agent.items() if k != "tail"}
record["agent_summary"] = agent["tail"][-800:]

L.copy_files(L.ALGO, cand, names)
changed = [n for n in names if os.path.exists(os.path.join(ws, n))
           and L.sha(os.path.join(ws, n)) != L.sha(os.path.join(L.ALGO, n))]
allowed = [n for n in changed if n in task["files"]]
big = [n for n in allowed if os.path.getsize(os.path.join(ws, n)) > L.MAX_FILE_BYTES]
L.copy_files(ws, cand, [n for n in allowed if n not in big])
record.update(changed=allowed, ignored_edits=[n for n in changed if n not in task["files"]],
              new_files=[f for f in L.ws_files(ws) if f not in names and f != "AGENTS.md"][:30])
L.log(f"iter {it:03d}: rescue - refereeing the finished agent run of {task_id} (changed: {', '.join(allowed)})")
ref = None
if big:
    accepted, reason = False, "file above 1 MB: " + ", ".join(big)
elif not allowed:
    accepted, reason = False, "no change to " + ", ".join(task["files"])
else:
    L.repo_guard()
    L.set_status(state="referee", iter=it, task=task_id, attempt=attempt, since=L.now())
    ref = L.run_referee(itdir, cand)
    accepted, reason = L.decide(task, ref, best)
if ref:
    s = ref.get("score_summary", {})
    record["referee"] = {k: ref.get(k) for k in ("verify_ok", "measure_ok", "score_ok", "thresholds",
                                                  "e_at_2224173", "seconds", "error")}
    record["referee"]["avg_s"] = s.get("choose_params_avg_s")
    record["referee"]["params_ppsz"] = ((ref.get("certificates") or {}).get("ppsz") or {}).get("params")
record.update(end=L.now(), reason=reason)
if accepted:
    nb, commit = L.accept(it, task, cand, ref, best)
    record.update(decision="accepted", commit=commit, best_after=nb["ppsz"])
    L.update_task(task_id, status="done", attempts=attempt, result=reason, done_iter=it)
    L.bulletin_add(f"iter {it:03d} {task_id} attempt {attempt}: ACCEPTED - {reason}")
else:
    record["decision"] = "rejected"
    L.update_task(task_id, status="open", attempts=0, result=reason)
    L.bulletin_add(f"iter {it:03d} {task_id}: rescue run REJECTED - {reason}")
L.append_jsonl(os.path.join(L.STATE, "log.jsonl"), record)
L.log(f"iter {it:03d}: {record['decision']} - {reason}")
L.set_status(state="idle", since=L.now())
print(record["decision"], reason)
