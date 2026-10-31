# Rules for this workspace (ETH research loop, Appendix D route)

You are the implementer in an automated research loop. A hidden referee scores the files you leave in this
workspace after you finish. A supervisor reads your final reply.

## Files
- appd_core.py: the Appendix D algorithm on small formulas (exact arithmetic). Functions that still raise
  NotImplementedError are yours to write; keep every signature and docstring contract.
- appd_params.py: choose_params(lam, target) -> parameter dict for the frozen cost model. The referee looks for
  the smallest lam at which e drops below 1/3 (the best known counting method for this class).
- appd_assembly_public.py: the frozen cost model. READ ONLY. The referee uses its own copy.
- smoke_appd.py: self-check. NOTES.md: background.

## How to test (run from the workspace, in bash)
    D:/ETH-loop/.venv/Scripts/python.exe smoke_appd.py
Run it once before you change anything and again after every change. It ends with SMOKE OK or SMOKE FAIL.

## Rules
1. Edit only the files the task names. Never edit appd_assembly_public.py.
2. Stay inside this workspace: never read or write D:/ETH-loop/referee_appd, D:/ETH-loop/state, D:/ETH-loop/algo
   or D:/ETH-loop/algo_appd, and do not look for the referee or its test cases.
3. No package installs, no network, no servers, never kill processes.
4. choose_params must be deterministic and take under 3 s per call.
5. All code must live in the files the task names; other files are discarded. Put throwaway experiments in
   scratch_*.py files. No file larger than 1 MB.
6. Use the venv python shown above for every Python command.
7. Never invent numbers: report only values that smoke_appd.py or your own scripts printed.

## Finish
When the stop condition of the task holds, reply with at most 6 lines: what you changed, the smoke_appd.py
result, and anything that did not work.
