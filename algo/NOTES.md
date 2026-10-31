# What this folder is

Code for a cost model of the counting algorithm in the paper "The Exponential Time
Hypothesis is False" (bounded-occurrence 3-CNF case). The goal of the project is to lower
the break-even size: the smallest lam = log2 n at which the algorithm's cost exponent e
drops below 0.38630 (the best known 3-SAT algorithm). Each task states the current best.

The cost exponent is e = taylor + meas + groups (see exponent_of in assembly_public.py):
- taylor = x + H2(x) with x = 3(R-1)/Q (number of inclusion-exclusion terms; wants R << Q),
- meas = rows * Q (measured basis vectors; wants the interval [A0, A0+R-1] to cover almost all
  of Bin(Q, (r-1)/(2r)), i.e. R of a few standard deviations sqrt(Q)/2),
- groups = 1/r (wants large r),
and the label space D = r^(2Q) must satisfy the balancing regime, which caps Q at
Qmax = floor(max_log2D(lam, regime) / (2 log2 r)).

## Files

- `model.py` - `choose_params(lam)`: picks the parameters r, Q, R, A0 and a regime.
  This is the main file to improve.
- `assembly_public.py` - copy of the frozen scoring formulas. Read it, use its functions,
  but do not edit it: the referee scores with its own copy, so edits here change nothing.
- `construction.py` - exact construction from Section 4 of the paper (used by the checks).
- `counting.py` - the counting algorithm on tiny formulas (used by the checks).
- `smoke.py` - quick public test. Run: `D:/ETH-loop/.venv/Scripts/python.exe smoke.py`
  (use exactly this python; the system python has no numpy).

## Rules

- Keep the public functions listed in the docstrings of construction.py and counting.py.
- choose_params must be deterministic and take under 3 seconds per call.
- Do not create files larger than 1 MB.
- A hidden referee checks correctness on other formulas and computes the real score.
