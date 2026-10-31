# Rules for this workspace (ETH research loop)

You are the implementer in an automated research loop. A hidden referee scores the files you
leave in this workspace after you finish. A supervisor reads your final reply.

## Files
- model.py: choose_params(lam, regime=None) -> dict(r, Q, R, A0, regime). It chooses the
  parameters of a #3-SAT counting algorithm for problem size n = 2**lam. The referee looks for
  the smallest lam at which the cost exponent e drops below 0.386300 (best known 3-SAT exponent).
- assembly_public.py: the frozen cost model (exponent_of, exponents_over_R, max_log2D, ...).
  READ ONLY. The referee uses its own copy, so editing this file cannot change any score.
- construction.py, counting.py: the algorithm itself (keep their public APIs).
- smoke.py: quick self-check. NOTES.md: background.

## How to test (run from the workspace, in bash)
    D:/ETH-loop/.venv/Scripts/python.exe smoke.py
It prints e (lower is better) at three lam values and FAILS if choose_params returns invalid
parameters or needs more than 3 s per call. Run it once before you change anything, so you know
the baseline, and again after every change.

## Rules
1. Edit only the files the task names. Never edit assembly_public.py.
2. Stay inside this workspace: never read or write D:/ETH-loop/referee, D:/ETH-loop/state or
   D:/ETH-loop/algo, and do not look for the referee or its test cases.
3. No package installs, no network, no servers, never kill processes.
4. choose_params must be deterministic and take under 3 s per call. Import only numpy, scipy,
   math and assembly_public.
5. All code must live in the files the task names; any other file you create is discarded.
   Put throwaway experiments in scratch_*.py files. No file larger than 1 MB.
6. Use the venv python shown above for every Python command.
7. Never invent numbers: report only values that smoke.py or your own scripts printed.

## Finish
When the stop condition of the task holds, reply with at most 6 lines: what you changed, the
smoke.py e values before and after, and anything that did not work.
