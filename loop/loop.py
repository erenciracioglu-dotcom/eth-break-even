"""ETH autoresearch loop driver.

Head (proposer / auditor): Claude. It writes tasks into state/backlog.json and audits
state/claims.md. Hands (implementer): the Eragent CLI driving Qwen3.8 Flash IQ3_S, served by
Strata on port 8080. Referee: referee/referee.py, frozen and outside the agent workspace.

Usage:
  python loop.py           run until STOP is touched, the backlog stays empty for 2 h, or errors pile up
  python loop.py --once    one real iteration, then exit (leaves Strata running)
  python loop.py --test    one dry iteration without the agent (checks the plumbing, changes no state)
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
import traceback
import urllib.request

ROOT = r"D:\ETH-loop"
ALGO = os.path.join(ROOT, "algo")
REFEREE = os.path.join(ROOT, "referee")
STATE = os.path.join(ROOT, "state")
WORK = os.path.join(ROOT, "work")
LOGS = os.path.join(ROOT, "logs")
LOOPDIR = os.path.join(ROOT, "loop")
BEST = os.path.join(STATE, "best.json")
BACKLOG = os.path.join(STATE, "backlog.json")
VENV_PY = os.path.join(ROOT, ".venv", "Scripts", "python.exe")
SYS_PY = r"C:\Users\Eren\AppData\Local\Programs\Python\Python312\python.exe"
ERAGENT = r"D:\Eragent\eragent.py"
GIT = r"C:\Program Files\Git\cmd\git.exe"
URL = "http://127.0.0.1:8080/v1"
STRATA_TASK = "ETH_Strata"
STOP = os.path.join(ROOT, "STOP")
LOCK = os.path.join(ROOT, "loop.lock")

AGENT_WALL_S = 40 * 60          # hard wall clock per agent run
AGENT_STEPS = 40                # default Eragent step cap (a task may set "steps")
WS_MAX_BYTES = 200 * 2 ** 20    # workspace size cap
WS_MAX_FILES = 3000             # workspace file-count cap
MAX_FILE_BYTES = 2 ** 20        # no accepted file above 1 MB
REFEREE_WALL_S = 60 * 60
STRATA_BOOT_S = 300
IDLE_EXIT_S = 2 * 3600          # empty backlog this long -> stop loop and the server it started
POLL_S = 300
MIN_GAIN = 0.0005               # a score task must lower the ppsz threshold by more than this
MAX_ERRORS = 5
NO_WINDOW = 0x08000000
FROZEN = ("assembly.py", "gates.py", "referee.py", "cases.json")
TRAILER = ("\n\nCo-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n"
           "Claude-Session: https://claude.ai/code/session_01Fmj5Vq76oDGgd8tqPeA1Go")
START = time.time()
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


# ------------------------------------------------------------------ small helpers
def now():
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    line = f"[{now()}] {msg}"
    print(line, flush=True)
    with open(os.path.join(LOGS, "loop.log"), "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def jload(path, default=None):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def jsave(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(tmp, path)


def append_jsonl(path, obj):
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, default=str) + "\n")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                          creationflags=NO_WINDOW, **kw)


def git(*args, check=False):
    p = run([GIT, "-C", ROOT, *args])
    if check and p.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {p.stderr.strip()[:500]}")
    return p.stdout.strip()


def kill_tree(pid):
    run(["taskkill", "/PID", str(pid), "/T", "/F"])


def pid_alive(pid):
    if not isinstance(pid, int):
        return False
    out = run(["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV", "/NH"]).stdout.lower()
    return f'"{pid}"' in out and "python" in out


def set_status(**kw):
    kw.update(updated=now(), pid=os.getpid())
    jsave(os.path.join(STATE, "status.json"), kw)


def stop_requested():
    try:
        return os.path.getmtime(STOP) > START
    except OSError:
        return False


def sleep_poll(seconds):
    t0 = time.time()
    while time.time() - t0 < seconds and not stop_requested():
        time.sleep(10)


def fill(template, **kw):
    for k, v in kw.items():
        template = template.replace("{{" + k + "}}", str(v))
    return template


def read_text(path, default=""):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except OSError:
        return default


def tail_text(path, chars=1500):
    return ANSI.sub("", read_text(path))[-chars:].strip()


# ------------------------------------------------------------------ Strata server
def strata_up():
    try:
        with urllib.request.urlopen(URL + "/models", timeout=5) as resp:
            return resp.status == 200
    except Exception:
        return False


def task_running(name):
    return "Running" in run(["schtasks", "/Query", "/TN", name, "/FO", "CSV", "/NH"]).stdout


def ensure_strata():
    if strata_up():
        return True
    if not task_running(STRATA_TASK):
        log("Strata is down: starting scheduled task " + STRATA_TASK)
        run(["schtasks", "/Run", "/TN", STRATA_TASK])
    t0 = time.time()
    while time.time() - t0 < STRATA_BOOT_S:
        time.sleep(10)
        if strata_up():
            log(f"Strata is up after {time.time() - t0:.0f} s")
            return True
    log("Strata did not answer within the boot limit")
    return False


def port_pids(port):
    pids = set()
    for line in run(["netstat", "-ano", "-p", "TCP"]).stdout.splitlines():
        parts = line.split()
        if len(parts) >= 5 and parts[1].endswith(f":{port}") and parts[3] == "LISTENING":
            pids.add(int(parts[4]))
    return pids


def stop_strata_if_ours():
    """Only a server launched through our scheduled task is stopped."""
    if not task_running(STRATA_TASK):
        return
    for pid in port_pids(8080):
        log(f"stopping the Strata server (PID {pid}) started by {STRATA_TASK}")
        kill_tree(pid)


# ------------------------------------------------------------------ repository guard
def tracked_algo_files():
    out = git("ls-files", "algo")
    return [line[len("algo/"):] for line in out.splitlines() if line.startswith("algo/")]


def repo_guard():
    """algo/ and referee/ change only through accepted iterations (committed by this loop)."""
    dirty = git("status", "--porcelain", "--untracked-files=no", "--", "algo", "referee")
    if dirty:
        log("WARNING: algo/referee were modified outside the loop; restoring: " + dirty.replace("\n", " | "))
        git("checkout", "--", "algo", "referee")
    ref = jload(os.path.join(STATE, "referee_sha.json"), {})
    for name in FROZEN:
        if sha(os.path.join(REFEREE, name)) != ref.get(name):
            raise RuntimeError(f"referee file {name} does not match state/referee_sha.json")
    return dirty


# ------------------------------------------------------------------ backlog
def pick_task(backlog):
    status = {t["id"]: t.get("status") for t in backlog}
    for t in sorted(backlog, key=lambda t: t.get("priority", 50)):
        if t.get("status") not in ("open", "running"):
            continue
        if any(status.get(dep) not in ("done", "failed") for dep in t.get("after", [])):
            continue
        return t
    return None


def update_task(task_id, **fields):
    """Re-read the backlog so edits the head made meanwhile are kept."""
    backlog = jload(BACKLOG, [])
    for t in backlog:
        if t["id"] == task_id:
            t.update(fields)
    jsave(BACKLOG, backlog)


def bulletin_add(line):
    path = os.path.join(STATE, "bulletin.md")
    lines = read_text(path).splitlines() or ["# Bulletin (latest first)"]
    body = [f"- {now()} {line}"] + lines[1:]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines[:1] + body[:60]) + "\n")


# ------------------------------------------------------------------ one iteration
def copy_files(src, dst, names):
    for name in names:
        d = os.path.join(dst, name)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        shutil.copy2(os.path.join(src, name), d)


def ws_stats(ws):
    files, total = 0, 0
    for root, _dirs, names in os.walk(ws):
        for n in names:
            try:
                total += os.path.getsize(os.path.join(root, n))
                files += 1
            except OSError:
                pass
    return files, total


def ws_files(ws):
    out = []
    for root, _dirs, names in os.walk(ws):
        if "__pycache__" in root:
            continue
        out += [os.path.relpath(os.path.join(root, n), ws).replace("\\", "/") for n in names]
    return out


def build_prompt(task, best, attempt):
    feedback = ""
    if attempt > 1 and task.get("last_feedback"):
        feedback = ("\nYour previous attempt was rejected: " + task["last_feedback"] +
                    "\nDo not repeat the same change; try a different approach.\n")
    return fill(read_text(os.path.join(LOOPDIR, "task_template.txt")),
                id=task["id"], attempt=attempt, max_attempts=task.get("max_attempts", 3),
                title=task["title"], files=", ".join(task["files"]), spec=task["spec"],
                ppsz=f"{best['ppsz']:,}", e_ref=f"{best['e_at_2224173']:.6f}",
                regime=best.get("regime", "paper"), feedback=feedback, stop=task["stop"],
                steps=task.get("steps", AGENT_STEPS))


def run_agent(itdir, ws, prompt, steps):
    logpath = os.path.join(itdir, "eragent.log")
    stats = os.path.join(itdir, "eragent_stats.jsonl")
    env = dict(os.environ, PYTHONIOENCODING="utf-8", ERAGENT_SPAWN="0", ERAGENT_REPLY_LANGUAGE="English")
    cmd = [SYS_PY, "-X", "utf8", ERAGENT, "--url", URL, "--workspace", ws, "--max-steps", str(steps),
           "--max-tokens", "16384", "--timeout", "600", "--stats", "--stats-log", stats, "-p", prompt]
    t0, killed = time.time(), None
    with open(logpath, "w", encoding="utf-8") as out:
        proc = subprocess.Popen(cmd, cwd=ws, env=env, stdout=out, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, creationflags=NO_WINDOW)
        while proc.poll() is None:
            time.sleep(15)
            files, total = ws_stats(ws)
            if time.time() - t0 > AGENT_WALL_S:
                killed = f"wall clock above {AGENT_WALL_S // 60} min"
            elif files > WS_MAX_FILES:
                killed = f"{files} files in the workspace"
            elif total > WS_MAX_BYTES:
                killed = f"workspace above {total / 2 ** 20:.0f} MB"
            if killed:
                log("killing the agent: " + killed)
                kill_tree(proc.pid)
                try:
                    proc.wait(60)
                except subprocess.TimeoutExpired:
                    pass
                break
    try:
        requests = sum(1 for _ in open(stats, encoding="utf-8"))
    except OSError:
        requests = 0
    return {"rc": proc.returncode, "seconds": round(time.time() - t0), "killed": killed,
            "requests": requests, "tail": tail_text(logpath)}


def run_referee(itdir, cand):
    out = os.path.join(itdir, "referee", "referee.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(os.path.join(itdir, "referee.log"), "w", encoding="utf-8") as fh:
        proc = subprocess.Popen([VENV_PY, "-I", os.path.join(REFEREE, "referee.py"), cand, out], cwd=itdir,
                                stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                creationflags=NO_WINDOW)
        try:
            proc.wait(REFEREE_WALL_S)
        except subprocess.TimeoutExpired:
            kill_tree(proc.pid)
            return {"error": f"referee above {REFEREE_WALL_S // 60} min"}
    return jload(out, {"error": "referee wrote no result"})


def gate_problem(ref):
    for gate, key in (("verify", "verify_ok"), ("measure", "measure_ok")):
        if not ref.get(key):
            g = ref.get(gate, {})
            if g.get("error"):
                return f"{gate} gate crashed: {g['error']}"
            bad = [i for i in g.get("items", []) if not i.get("ok")]
            return f"{gate} gate failed on " + "; ".join(json.dumps({k: v for k, v in i.items() if k != 'ok'})[:300]
                                                          for i in bad[:3])
    s = ref.get("score_summary", {})
    if ref.get("tamper"):
        return "score gate: reported exponents differ from the clean re-scoring (assembly was patched at run time)"
    if s.get("error"):
        return f"score gate crashed: {s['error']}"
    if s.get("choose_params_avg_s", 0) >= 3.0:
        return f"choose_params too slow: {s['choose_params_avg_s']:.2f} s per call (limit 3 s)"
    if not (ref.get("thresholds") or {}).get("ppsz"):
        return "score gate: the ppsz target was never reached"
    return None


def decide(task, ref, best):
    if not ref or ref.get("error"):
        return False, (ref or {}).get("error", "referee produced no result")
    problem = gate_problem(ref)
    if problem:
        return False, problem
    new, old = ref["thresholds"]["ppsz"], best["ppsz"]
    e_new, e_old = ref.get("e_at_2224173"), best.get("e_at_2224173")
    detail = f"ppsz threshold {old:,} -> {new:,}"
    if e_new is not None and e_old is not None:
        detail += f"; e at 2224173 {e_old:.6f} -> {e_new:.6f}"
    if task.get("kind", "score") == "score":
        if new < old * (1 - MIN_GAIN):
            return True, detail + f" ({100 * (1 - new / old):.2f}% lower threshold)"
        return False, "no gain: " + detail
    if new <= old * (1 + MIN_GAIN):
        return True, "gates pass, no regression: " + detail
    return False, "regression: " + detail


def accept(it, task, cand, ref, best):
    copy_files(cand, ALGO, task["files"])
    cert = ref["certificates"]["ppsz"]
    regime = cert["params"].get("regime", "paper")
    nb = dict(best)
    nb.update(ppsz=ref["thresholds"]["ppsz"], sharp3sat=ref["thresholds"].get("sharp3sat"),
              brute=ref["thresholds"].get("brute"), e_at_2224173=ref.get("e_at_2224173"),
              params_ppsz=cert["params"], regime=regime, iter=it, task=task["id"], updated=now())
    by = nb.setdefault("by_regime", {})
    if regime not in by or nb["ppsz"] < by[regime]["ppsz"]:
        by[regime] = {"ppsz": nb["ppsz"], "iter": it, "params": cert["params"]}
    jsave(BEST, nb)
    if regime != "paper":
        with open(os.path.join(STATE, "claims.md"), "a", encoding="utf-8") as fh:
            fh.write(f"- {now()} iter {it:03d} ({task['id']}): ppsz threshold {nb['ppsz']:,} uses regime "
                     f"'{regime}' -> PROVISIONAL until the head audits that regime.\n")
    git("add", "algo", "state", check=True)
    msg = f"iter {it:03d} {task['id']}: ppsz threshold {best['ppsz']:,} -> {nb['ppsz']:,} [{regime}]"
    git("commit", "-m", msg + TRAILER, check=True)
    return nb, git("rev-parse", "--short", "HEAD")


def next_iter_no(prefix):
    n = 0
    for name in os.listdir(WORK):
        m = re.fullmatch(prefix + r"_(\d+)", name)
        if m:
            n = max(n, int(m.group(1)))
    return n + 1


def iteration(test=False):
    backlog = jload(BACKLOG, [])
    task = pick_task(backlog)
    if task is None:
        return "idle"
    best = jload(BEST)
    guard = repo_guard()
    if not test and not ensure_strata():
        return "strata_down"
    prefix = "test" if test else "iter"
    it = next_iter_no(prefix)
    itdir = os.path.join(WORK, f"{prefix}_{it:03d}")
    ws, cand = os.path.join(itdir, "ws"), os.path.join(itdir, "candidate")
    names = tracked_algo_files()
    copy_files(ALGO, ws, names)
    attempt = int(task.get("attempts", 0)) + 1
    with open(os.path.join(ws, "AGENTS.md"), "w", encoding="utf-8") as fh:
        fh.write(read_text(os.path.join(LOOPDIR, "agents_template.md")))
    prompt = build_prompt(task, best, attempt)
    with open(os.path.join(itdir, "prompt.txt"), "w", encoding="utf-8") as fh:
        fh.write(prompt)
    record = {"iter": it, "test": test, "task": task["id"], "attempt": attempt, "start": now(),
              "best_before": best["ppsz"], "guard": guard or None}
    log(f"{prefix} {it:03d}: task {task['id']} attempt {attempt} - {task['title']}")
    if test:
        agent = {"rc": 0, "seconds": 0, "killed": None, "requests": 0, "tail": "(test: no agent)"}
    else:
        update_task(task["id"], status="running")
        set_status(state="agent", iter=it, task=task["id"], attempt=attempt, since=now())
        agent = run_agent(itdir, ws, prompt, task.get("steps", AGENT_STEPS))
        log(f"agent finished: rc={agent['rc']} {agent['seconds']} s, {agent['requests']} requests"
            + (f", killed: {agent['killed']}" if agent["killed"] else ""))
    record["agent"] = {k: v for k, v in agent.items() if k != "tail"}
    record["agent_summary"] = agent["tail"][-800:]
    infra = (not test) and agent["rc"] != 0 and not agent["killed"] and not strata_up()

    copy_files(ALGO, cand, names)
    changed = [n for n in names if os.path.exists(os.path.join(ws, n))
               and sha(os.path.join(ws, n)) != sha(os.path.join(ALGO, n))]
    allowed = [n for n in changed if n in task["files"]]
    big = [n for n in allowed if os.path.getsize(os.path.join(ws, n)) > MAX_FILE_BYTES]
    copy_files(ws, cand, [n for n in allowed if n not in big])
    record.update(changed=allowed, ignored_edits=[n for n in changed if n not in task["files"]],
                  new_files=[f for f in ws_files(ws) if f not in names and f != "AGENTS.md"][:30])

    ref = None
    if infra:
        accepted, reason = False, "infrastructure: model server went down during the run (attempt not counted)"
    elif big:
        accepted, reason = False, "file above 1 MB: " + ", ".join(big)
    elif not allowed and not test:
        accepted, reason = False, "no change to " + ", ".join(task["files"])
    else:
        repo_guard()
        set_status(state="referee", iter=it, task=task["id"], attempt=attempt, since=now())
        ref = run_referee(itdir, cand)
        accepted, reason = decide(task, ref, best)
    if ref:
        s = ref.get("score_summary", {})
        record["referee"] = {k: ref.get(k) for k in ("verify_ok", "measure_ok", "score_ok", "thresholds",
                                                      "e_at_2224173", "seconds", "error")}
        record["referee"]["avg_s"] = s.get("choose_params_avg_s")
        record["referee"]["params_ppsz"] = ((ref.get("certificates") or {}).get("ppsz") or {}).get("params")
    record.update(end=now(), reason=reason)

    if test:
        record["decision"] = "test-accept" if accepted else "test-reject"
        append_jsonl(os.path.join(STATE, "test_log.jsonl"), record)
        log(f"test iteration: {record['decision']} - {reason}")
        return "test"
    if accepted:
        nb, commit = accept(it, task, cand, ref, best)
        record.update(decision="accepted", commit=commit, best_after=nb["ppsz"])
        update_task(task["id"], status="done", attempts=attempt, result=reason, done_iter=it)
        bulletin_add(f"iter {it:03d} {task['id']} attempt {attempt}: ACCEPTED - {reason}")
    elif infra:
        record["decision"] = "infra"
        update_task(task["id"], status="open")
        bulletin_add(f"iter {it:03d} {task['id']}: infrastructure failure, will retry - {reason}")
    else:
        maxa = int(task.get("max_attempts", 3))
        record["decision"] = "rejected"
        feedback = reason + ". Agent's own summary: " + agent["tail"][-500:].replace("\n", " ")
        update_task(task["id"], status="failed" if attempt >= maxa else "open", attempts=attempt,
                    last_feedback=feedback[:1500], result=reason)
        bulletin_add(f"iter {it:03d} {task['id']} attempt {attempt}/{maxa}: REJECTED - {reason}")
    append_jsonl(os.path.join(STATE, "log.jsonl"), record)
    log(f"iter {it:03d}: {record['decision']} - {reason}")
    return record["decision"]


def render_report():
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("report", os.path.join(LOOPDIR, "report.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.render(ROOT)
    except Exception:
        log("report failed: " + traceback.format_exc()[-600:])


# ------------------------------------------------------------------ main
def take_lock():
    old = jload(LOCK, {})
    if old and old.get("pid") and old.get("pid") != os.getpid() and pid_alive(old["pid"]):
        return False
    jsave(LOCK, {"pid": os.getpid(), "since": now()})
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--test", action="store_true")
    args = ap.parse_args()
    os.makedirs(LOGS, exist_ok=True)
    os.makedirs(WORK, exist_ok=True)
    if not take_lock():
        log("another loop process holds the lock; exiting")
        return
    mode = "test" if args.test else ("once" if args.once else "loop")
    log(f"loop started ({mode}), pid {os.getpid()}")
    idle_since, errors = None, 0
    try:
        while not stop_requested():
            try:
                res = iteration(test=args.test)
                errors = 0
            except Exception:
                errors += 1
                res = "error"
                log("ERROR in iteration:\n" + traceback.format_exc())
                if errors >= MAX_ERRORS:
                    log("too many errors in a row; stopping")
                    break
            render_report()
            if args.once or args.test:
                break
            if res == "idle":
                idle_since = idle_since or time.time()
                set_status(state="idle", since=dt.datetime.fromtimestamp(idle_since).strftime("%Y-%m-%d %H:%M:%S"))
                if time.time() - idle_since > IDLE_EXIT_S:
                    log("backlog empty for 2 h; stopping")
                    break
                sleep_poll(POLL_S)
            else:
                idle_since = None
                if res in ("strata_down", "error"):
                    sleep_poll(120)
        if stop_requested():
            log("STOP requested")
    finally:
        if mode == "loop":
            stop_strata_if_ours()
        if mode != "test":
            set_status(state="stopped", since=now())
        render_report()
        jsave(LOCK, {"pid": None, "released": now()})
        log("loop exited")


if __name__ == "__main__":
    main()
