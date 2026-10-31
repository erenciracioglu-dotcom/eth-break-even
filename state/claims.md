# Claims ledger

PROVISIONAL = accepted by the referee under a head extension that is not yet audited.
AUDITED = the head checked the argument. REJECTED = the argument failed and results built on
it are reverted.

## Regimes

- paper [BASE]: Lemma 4.5 of the paper as written: 2(b+w) D^2 <= M, h <= M/D + sqrt(2(b+w) M).
  Accounting refinements used by the frozen assembly (rows touch <= 2 groups, restricted row
  abs sum <= 2^(Q-1), replacement-row degree <= R-1, n/r groups) are checked on real designs by
  the measure gate on every referee run.
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
    Group sizes: Pr[n_t,g >= 2n/r] <= exp(-0.386 n/r) (Chernoff).
    Union bound: D exp(-0.110 M/D) + Qr exp(-0.386 n/r) <= D^(-1.64) + tiny < 1 when
    M/D >= 24 ln D (the regime) and n/r >> ln(Qr). Fix a_v one variable at a time with the
    product estimator (D + Qr terms, each a product over n variables): poly(n, D) time.
  Consequence: valid; it would even allow M/D >= 8.4 ln D, which moves log2 D by about 1.5 bits
  (no measurable gain). D cannot usefully exceed about M anyway, so this regime is close to the
  ceiling for D; further gains need a smaller label space per unit of Q (head research).
  Caveat: a head-level sketch, not peer reviewed. The measure gate also checks h <= 2M/D on a
  real freedman-regime design on every referee run.

## Results that depend on a provisional regime
