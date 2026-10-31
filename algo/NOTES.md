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

## Version 5 of the cost model (2026-10-31)

- Parameter chunk (int L, 2 <= L <= 4096, needs units=True): units are connected chunks of the
  variable graph with at least L variables (a component with fewer than L variables is one unit).
  Worst case M_eff <= (1 + 1/L) n instead of 12n/7. Call max_log2D(lam, regime, delta, True, True, L)
  and exponents_over_R(lam, r, Q, Rs, regime, delta, True, aggregated, True, L) and return chunk=L.
  The referee counts chunk parameters only when its hidden checks of
  build_design(..., designated=True, units=True, chunk=L) pass.
- Parameter modular (bool, needs designated and aggregated; yes/no targets only): the yes/no
  question is answered by a weighted count mod a prime p with Q < p <= 2Q, so a measurement value
  costs log2(2Q) bits instead of Q. Only exponent_of(params, lam, decision=True) accepts it; the
  referee scores the 1.307^n target with model.choose_params_decision(lam) when model.py defines
  it (and with choose_params otherwise). Counting targets always use choose_params. The referee
  counts modular parameters only when its hidden checks of Counter(..., residue_keys=True,
  weights=w) pass.

## Version 6 of the cost model (2026-10-31)

- Parameter perpoint (bool): measurement matrices per evaluation point (Definition 3.1 allows
  C_{t,g} to depend on t). At point t the swap of position t is free, so only basis vectors with an
  antisymmetric factor at position t and m(j - 1) + m(j) != 0 are measured there,
  m(j) = 1 - (-1)^j f_I(j): out = p * Pr[Bin(Q - 1, p) outside [A0, A0 + R - 2]], about half of
  the fixed matrix. Pass perpoint=True to exponents_over_R / exponent_of (keyword) and return
  perpoint=True. Aux(..., perpoint=True) provides Cmeas_t (one row matrix per point); designs carry
  des['perpoint'] and des['Chat_t']. The referee counts perpoint parameters only when its hidden
  checks of Aux(..., perpoint=True) (separation at every point with Cmeas_t, row counts) and of such
  designs (exact and weighted counts) pass.
- Parameter eadapt (bool): the interpolation window depends on the number e of equal pairs of a
  label, I_e = [a_e, a_e + R - 1], a_e = (Q - e - R + 1) // 2 (a_e = 0 when Q - e < R); A0 is
  ignored. exponents_over_R with eadapt=True is about ten times slower (20 to 40 ms per call), so
  use it only in a final re-scan around the winner of the perpoint search. The referee counts
  eadapt parameters only when its hidden checks of Aux(..., eadapt=True) pass.
- assembly_public.tail_out(Q, r, R, A0, perpoint=False, eadapt=False) gives the measured fraction
  (argument order Q, r, R, A0).

## Version 7 of the cost model (2026-10-31)

- Parameter contract (bool, needs aggregated): the groups of a point are summed one at a time
  (variable elimination) by Counter.Zagg_contracted instead of enumerating every measurement tuple.
  A clause of group g sees group h only through the rows of the pair {g, h} restricted to h, so the
  carried table is keyed by the half-bonds between processed and remaining groups; the measurement
  term is multiplied by phi(r) = assembly_public.contract_factor(r) (about 1/2). Pass contract=True
  to exponents_over_R / exponent_of (keyword) and return contract=True. The referee counts contract
  parameters only when its hidden checks (the bond condition on real designs, exact counts with
  Zagg_contracted) pass.

## Version 8 of the cost model (2026-10-31)

- Parameter smallprime (bool, yes/no targets only, needs modular): the prime only has to exceed R - 1.
  The interpolation nodes are the first Q nonzero elements of GF(p^2) = GF(p)[s] / (s^2 - c) (c the
  least quadratic non-residue mod p) instead of the integers 1..Q, and the rows use the factors
  (1 - T / alpha_s). The measurement matrices stay integer, so a measured value is a residue mod p and
  costs log2 p bits, p the least odd prime above max(R - 1, 2) (assembly_public.value_bits), instead of
  log2(2Q). Pass smallprime=True to exponents_over_R / exponent_of (keyword) together with modular=True
  and return smallprime=True. The referee counts smallprime parameters only when its hidden checks of
  construction.gf_rows(aux, des, p) pass (Definition 3.1(4) at every GF(p^2) node and the weighted count
  from Taylor data at the nodes).

## Version 9 of the cost model (2026-10-31)

- Parameter edges (bool, needs contract): Counter.taylor_edges(t) computes the Taylor data H(t + z) mod z^J
  inside the group-by-group contraction. At point t a replaced occurrence (clause in group d, variable in
  group a != d) is an edge: its row meets only the groups d and a, and its z-part on group a has no
  z^0 term. Every clause is expanded multilinearly in these z-parts, so the contraction carries one bit per
  open edge (at most J - 1 bits set) instead of enumerating one restricted sum per set U. The Taylor term
  becomes E H2(x / E) with E = (1 + delta) m_eff eps(r), eps(r) = assembly_public.edges_factor(r) (about
  1/2), instead of H2(x). Pass edges=True to exponents_over_R / exponent_of (keyword) together with
  contract=True and return edges=True. The referee counts edges parameters only when its hidden checks pass
  (two-group support of every replaced row at every point, taylor_edges equal to a direct evaluation of the
  Taylor data at every point, exact and weighted counts).
