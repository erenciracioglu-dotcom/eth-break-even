"""Renders D:/ETH-loop/report.html from the state files (called by loop.py after each iteration)."""
import html
import json
import os
import sys
import time

CSS = """
:root { --bg:#fbfaf7; --fg:#1d1d1b; --mut:#6b6a66; --line:#e3e0d8; --card:#ffffff;
        --ok:#1f7a4d; --bad:#a3352b; --warn:#8a6100; --acc:#2f5f8f; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#161615; --fg:#ecebe6; --mut:#9c9a93; --line:#2f2e2b; --card:#1e1e1c;
          --ok:#5cc28d; --bad:#e2766b; --warn:#d9a640; --acc:#8ab4e0; } }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--fg);
       font:14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width:1180px; margin:0 auto; padding:24px 16px 48px; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:16px; margin:28px 0 8px; }
.sub { color:var(--mut); margin:0 0 18px; }
.cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:12px; }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
.card .k { color:var(--mut); font-size:12px; } .card .v { font-size:20px; font-variant-numeric:tabular-nums; }
.card .s { color:var(--mut); font-size:12px; overflow-wrap:anywhere; }
.wrap { overflow-x:auto; border:1px solid var(--line); border-radius:10px; background:var(--card); }
table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums; }
th, td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--line); vertical-align:top; }
th { color:var(--mut); font-weight:600; font-size:12px; }
tr:last-child td { border-bottom:0; }
td.reason { min-width:280px; color:var(--mut); }
.accepted { color:var(--ok); font-weight:600; } .rejected, .failed { color:var(--bad); }
.infra, .running, .open { color:var(--warn); } .done { color:var(--ok); }
pre { white-space:pre-wrap; background:var(--card); border:1px solid var(--line); border-radius:10px;
      padding:12px 14px; margin:0; font:12.5px/1.5 ui-monospace, Consolas, monospace; }
"""


def _j(path, default):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return default


def _t(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except Exception:
        return ""


def _n(x):
    return f"{x:,}" if isinstance(x, (int, float)) and not isinstance(x, bool) else "-"


def _e(x):
    return f"{x:.6f}" if isinstance(x, (int, float)) else "-"


def render(root):
    st = os.path.join(root, "state")
    esc = html.escape
    best = _j(os.path.join(st, "best.json"), {})
    backlog = _j(os.path.join(st, "backlog.json"), [])
    status = _j(os.path.join(st, "status.json"), {})
    recs = []
    for line in _t(os.path.join(st, "log.jsonl")).splitlines():
        try:
            recs.append(json.loads(line))
        except ValueError:
            pass
    accepted = sum(1 for r in recs if r.get("decision") == "accepted")

    base = 2224173
    gain = ""
    if isinstance(best.get("ppsz"), int):
        gain = f"{base / best['ppsz']:.3f}x lower than the start ({base:,})"
    by = best.get("by_regime", {})
    by_txt = "; ".join(f"{k}: {_n(v.get('ppsz'))}" for k, v in by.items()) or "-"
    cards = [
        ("Loop", esc(str(status.get("state", "-"))),
         esc(f"iter {status.get('iter', '-')} task {status.get('task', '-')} since {status.get('since', '-')}")),
        ("Best ppsz threshold (lam = log2 n)", _n(best.get("ppsz")), esc(gain)),
        ("e at lam = 2,224,173", _e(best.get("e_at_2224173")), esc(f"regime {best.get('regime', '-')}")),
        ("By regime", esc(by_txt), esc("freedman results stay PROVISIONAL until audited")),
        ("Iterations", str(len(recs)), f"{accepted} accepted"),
    ]
    card_html = "".join(f"<div class='card'><div class='k'>{k}</div><div class='v'>{v}</div>"
                        f"<div class='s'>{s}</div></div>" for k, v, s in cards)

    blog_rows = "".join(
        f"<tr><td>{esc(t['id'])}</td><td>{esc(t.get('title', ''))}</td>"
        f"<td class='{esc(str(t.get('status', '')))}'>{esc(str(t.get('status', '')))}</td>"
        f"<td>{t.get('attempts', 0)}/{t.get('max_attempts', 3)}</td>"
        f"<td class='reason'>{esc(str(t.get('result', '')))}</td></tr>" for t in backlog)

    it_rows = []
    for r in reversed(recs[-300:]):
        ref = r.get("referee") or {}
        ag = r.get("agent") or {}
        dec = str(r.get("decision", ""))
        it_rows.append(
            f"<tr><td>{r.get('iter', 0):03d}</td><td>{esc(str(r.get('start', '')))}</td>"
            f"<td>{esc(str(r.get('task', '')))}</td><td>{r.get('attempt', '')}</td>"
            f"<td>{ag.get('seconds', '-')}</td><td>{ag.get('requests', '-')}</td>"
            f"<td>{_n((ref.get('thresholds') or {}).get('ppsz'))}</td><td>{_e(ref.get('e_at_2224173'))}</td>"
            f"<td class='{esc(dec)}'>{esc(dec)}</td><td class='reason'>{esc(str(r.get('reason', '')))}</td></tr>")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="refresh" content="60">
<title>ETH loop</title><style>{CSS}</style></head>
<body><main>
<h1>ETH autoresearch loop</h1>
<p class="sub">Head: Claude (tasks, audits) - hands: Eragent + Qwen3.8 Flash IQ3_S - referee: frozen gates.
Page written {esc(time.strftime('%Y-%m-%d %H:%M:%S'))}, refreshes every minute.</p>
<div class="cards">{card_html}</div>
<h2>Backlog</h2>
<div class="wrap"><table><tr><th>Task</th><th>Title</th><th>Status</th><th>Attempts</th><th>Result</th></tr>
{blog_rows}</table></div>
<h2>Iterations</h2>
<div class="wrap"><table><tr><th>Iter</th><th>Start</th><th>Task</th><th>Try</th><th>Agent s</th>
<th>Requests</th><th>ppsz threshold</th><th>e at 2,224,173</th><th>Decision</th><th>Reason</th></tr>
{''.join(it_rows) or '<tr><td colspan="10">No iterations yet.</td></tr>'}</table></div>
<h2>Claims ledger</h2>
<pre>{esc(_t(os.path.join(st, 'claims.md')))}</pre>
<h2>Bulletin</h2>
<pre>{esc(_t(os.path.join(st, 'bulletin.md')))}</pre>
</main></body></html>
"""
    tmp = os.path.join(root, "report.html.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(page)
    os.replace(tmp, os.path.join(root, "report.html"))


if __name__ == "__main__":
    render(sys.argv[1] if len(sys.argv) > 1 else r"D:\ETH-loop")
