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

## Version 2 of the cost model (2026-10-31)

- Regime "tight": extra parameter delta in [1e-6, 1]; D ln D <= delta^2 M / (4(b+w)) and
  h <= (1 + delta) M/D. delta = 1 equals "freedman". Call max_log2D(lam, "tight", delta) and
  exponents_over_R(lam, r, Q, Rs, "tight", delta), and return "delta" in the dict.
- Parameter designated (bool, default False): every clause reads one of its own variables
  directly, so only M_eff = sum_f (|I_f| - 1) <= 2n occurrences carry replacement rows and the
  taylor term uses 2 instead of 3 occurrences per variable. Call max_log2D(lam, regime, delta,
  True) and exponents_over_R(lam, r, Q, Rs, regime, delta, True) and return designated=True.
  Allowed with regimes "freedman" and "tight" only. The hidden referee counts designated
  parameters only when its own hidden checks of build_design(..., designated=True) pass.

## Version 3 of the cost model (2026-10-31)

- Parameter aggregated (bool, default False): Lemma 3.3 summed per U, so each point needs one
  restricted sum per U with |U| <= q instead of one per pair U <= E. The taylor term becomes
  H2(x) (1 for x >= 1/2) instead of x + H2(x). Call exponents_over_R(lam, r, Q, Rs, regime,
  delta, designated, True) and return aggregated=True. The hidden referee counts aggregated
  parameters only when its own hidden checks of Counter.taylor_aggregated with
  Counter.Zagg_group_separated pass (exact counts and the number of restricted sums).

## Version 4 of the cost model (2026-10-31)

- Parameter units (bool, default False, needs designated=True): a greedy maximal matching of
  clauses; the variables of every matched clause share one schedule (a unit), so all of that
  clause's occurrences are read directly. Worst case M_eff <= 12n/7 instead of 2n. Call
  max_log2D(lam, regime, delta, True, True) and exponents_over_R(lam, r, Q, Rs, regime, delta,
  True, aggregated, True), and return units=True. The hidden referee counts units parameters only
  when its own hidden checks of build_design(..., designated=True, units=True) pass.
