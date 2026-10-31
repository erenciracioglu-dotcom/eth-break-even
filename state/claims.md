# Claims ledger

PROVISIONAL = accepted by the referee under a head extension that is not yet audited.
AUDITED = the head checked the argument. REJECTED = the argument failed and results built on
it are reverted.

## Regimes

- paper [BASE]: Lemma 4.5 of the paper as written: 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M).
  Accounting refinements used by the frozen assembly (rows touch <= 2 groups, restricted row
  abs sum <= 2^(Q-1), replacement-row degree <= R-1) are checked on real designs by the
  measure gate on every referee run. The group term 1/r (the paper's (2.11) charges 2n/r per
  group) uses group sizes (1 + o(1)) n/r: the derandomised bound n/r + sqrt(2Qn) from the proof
  of Lemma 4.5 is (1 + o(1)) n/r at these sizes. Not checked by a gate (head argument).
- freedman [AUDITED as a head proof sketch, 2026-10-31 12:05; was PROVISIONAL]:
  D ln D <= M / (4(b+w)) with h <= 2M/D and |S_t,g| <= 2n/r replaces the paper's condition
  2(b+w) D^2 <= M in Lemma 4.5 (same conclusions, same polynomial construction time).
  The name is historical: the audited argument uses Bernstein's inequality in two stages, so
  every sum is linear in independent choices and the pessimistic estimators are explicit
  products (a one-stage bound on the quadratic form sum_(f,v) 1{a_v=alpha} 1{d_f=beta} would
  need an exponential moment of a quadratic form, which has no easy estimator).
  Notation: K = r^Q, D = K^2, b = max incidences per variable (3), w = max per factor (3).
  Stage 1 (factor schedules d_f): m_beta = sum_f |I_f| 1{d_f = beta} has independent terms in
    [0, w] and mean M/K. Bernstein + union bound over the K values of beta: every
    m_beta <= (1+eps) M/K once eps^2 M/K > 2w(1+eps/3) ln K. Under the regime
    M/K = K (M/D) >= 48 K ln K, so eps = 0.05 works for K >= 51 (in the scored range
    K = r^Q is astronomically large). Estimator: sum_beta exp(-s(1+eps)M/K) prod_f E[exp(s|I_f| 1{d_f=beta})].
  Stage 2 (variable schedules a_v, d fixed): c_(alpha,beta) = sum_v 1{a_v=alpha} N_v(beta) with
    N_v(beta) = #{f containing v : d_f = beta} in [0, b]; independent terms, mean m_beta/K <=
    (1+eps) M/D, variance <= b(1+eps) M/D. Bernstein with t = 2M/D - mean >= (1-eps) M/D:
    Pr[c_u > 2M/D] <= exp(-3(1-eps)^2 (M/D) / (4b(2+eps))) = exp(-0.110 M/D) at eps = 0.05, b = 3.
    Group sizes: Pr[n_t,g >= (1+d') n/r] <= exp(-d'^2 n/(3r)) (Chernoff); with d' = (n/r)^(-1/3)
    this is negligible, so |S_t,g| = (1 + o(1)) n/r, as the cost model's 1/r term assumes.
    Union bound: D exp(-0.110 M/D) + Qr exp(-(n/r)^(1/3)/3) <= D^(-1.64) + tiny < 1 when
    M/D >= 24 ln D (the regime) and n/r is large. Fix a_v one variable at a time with the
    product estimator (D + Qr terms, each a product over n variables): poly(n, D) time.
  Consequence: valid; it would even allow M/D >= 8.4 ln D, which moves log2 D by about 1.5 bits
  (no measurable gain). D cannot usefully exceed about M anyway, so this regime is close to the
  ceiling for D; further gains need a smaller label space per unit of Q (head research).
  Caveat: a head-level sketch, not peer reviewed. The measure gate also checks h <= 2M/D on a
  real freedman-regime design on every referee run.

- tight [AUDITED as a head proof sketch, 2026-10-31; referee v2]: the freedman argument with the
  target (1+delta) M/D instead of 2M/D, 0 < delta <= 1. Stage 2 Bernstein with t = (delta-eps) M/D:
  Pr[c_u > (1+delta) M/D] <= exp(-(delta-eps)^2 (M/D) / (6(1+eps) + 2(delta-eps))) for b = 3.
  With eps = delta/4 and M/D >= 24 ln D / delta^2 (the regime) the exponent is >= 1.5 ln D, so
  the union bound over the D labels is <= D^(-1/2) (+ negligible group-size terms). Stage 1
  needs eps^2 M/K > 2w(1+eps/3) ln K, true because M/K >= 48 K ln K / delta^2.
  Effect: hD/n drops from 6 to 3(1+delta); log2 D loses 2 log2(1/delta) + O(1) bits.
- designated [AUDITED as a head proof sketch, 2026-10-31 12:30; exact counts checked by
  referee v2 whenever a score uses it]: every clause f reads one of its own variables v(f)
  directly and is assigned to that variable's group at every point, sigma_f(t) = a_v(f)(t).
  1. Definition 3.1(4) holds for the constant row P_f,v(f)(T) = e_v(f): P(0) = delta, degree 0,
     and for g != sigma_f(t) the variable v(f) is not in S_t,g, so <P(t), xi~> = 0 for every xi
     supported on S_t,g (no measurement needed).
  2. Non-designated occurrences keep the paper's construction with label (a_v, d_f), d_f = a_v(f);
     Lemma 4.3 only needs d_u(t) = sigma_f(t), which holds.
  3. Clause factors are multilinear in their inputs, so deg_T H <= M_eff (R-1) with
     M_eff = sum_f (|I_f| - 1) <= (2/3) M <= 2n for widths <= 3. The Taylor order
     J = floor(M_eff (R-1)/Q) + 1 is enough for Hermite interpolation (QJ > deg H), and the
     inclusion-exclusion bound |E| <= q = J - 1 shrinks by the factor 2/3: taylor(2(R-1)/Q).
  4. Only the M_eff non-designated occurrences need auxiliary positions: hD <= (1+delta) M_eff.
  5. Balancing: d_f = a_v(f) makes label counts quadratic in the a's. Properly colour the
     conflict graph (variables sharing a clause; max degree 6, so <= 7 colours, greedy) and fix
     the a's one colour class per stage. Each non-designated occurrence (f, v) is an edge
     (v, v(f)) between different colours, completed at the later colour; given earlier stages
     its label indicator is linear in the current stage (each variable touches <= 6 edges).
     Per stage: Bernstein for the new label counts around their conditional means (deviation
     budget delta (E_i + M_eff/7) / (2D) per stage, which sums to delta M_eff / D), half-label
     balance events for later stages and group sizes (linear, huge means), union bound over
     28 D events, one product-form pessimistic estimator per stage. The worst stage exponent is
     >= delta^2 (M_eff/D) / 112, so D ln(28 D) <= delta^2 M_eff / 128 suffices (cost: a few bits
     of log2 D). Conclusions: h <= (1+delta) M_eff / D and |S_t,g| = (1+o(1)) n/r.
  Evidence in code: the head's reference implementation passes referee v2 with exact counts at
  the reduced Taylor order (for example J = 10 instead of 14 for n = 6, m = 9, Q = 2, R = 2).
  Caveat: head-level sketch, not peer reviewed.

- aggregated [AUDITED as a head proof sketch, 2026-10-31 12:45; exact counts and the number of
  restricted sums are checked by referee v3 whenever a score uses it]: Lemma 3.3 summed per U.
  For x >= U: sum over U <= E <= x with |E| <= q of (-1)^(|E|-|U|) equals
  w_K(|x minus U|) = sum_{k <= K} (-1)^k C(|x minus U|, k) with K = q - |U|, so
  sum_x p_x(x) = sum_{|U| <= q} sum_{x >= U} w_K(|x minus U|) p_x(1_U).
  The weight depends on x only through |x minus U| = sum_g |x_g minus U_g|. With the measurements
  fixed, the restricted sum factorises over groups (Lemma 3.2), so sum_x (1+y)^(|x minus U|) G_x(1_U, z)
  is a product of per-group polynomials in y (truncated at y^K): polynomial overhead only.
  Hence one group-separated sum per U with |U| <= q, i.e. sum_{s <= q} C(n, s) <= 2^(n H2(q/n))
  restricted sums per point (q <= n/2) instead of sum_s C(n, s) 2^s: taylor term H2(x).
  Evidence: in the cloud, 13 cases (truncated Taylor orders, plain and designated designs) give
  identical Taylor data and exact counts; the head's reference implementation passes referee v3.
  Caveat: head-level sketch, not peer reviewed.

## Results that depend on a provisional regime
- 2026-10-31 12:14:25 iter 002 (T2): ppsz threshold 1,117,111 uses regime 'freedman' -> PROVISIONAL until the head audits that regime. [head: regime audited 12:05, this result stands as audited]
- 2026-10-31 12:28:56 iter 003 (T3): ppsz threshold 1,070,239 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 12:42:32 iter 004 (T4): ppsz threshold 438,437 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
