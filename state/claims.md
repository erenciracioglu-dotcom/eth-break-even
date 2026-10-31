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

- units [AUDITED as a head proof sketch, 2026-10-31 13:00; unit structure, maximality and
  exact counts are checked by referee v4 whenever a score uses it]: a greedy maximal matching of
  clauses (each clause matched only if none of its variables is in a unit yet). The variables of
  a matched clause share one schedule, so the clause is assigned to that schedule's group at every
  point and all its inputs are read directly (constant rows, as in 'designated'); any other clause
  is assigned to the unit holding most of its variables and reads all of those directly.
  1. Worst case (any mix of widths <= 3): matched clauses contribute 0 to M_eff, an unmatched
     clause at most |I_f| - 1 <= 2 and <= (2/3)|I_f|. By maximality every unmatched clause of
     width >= 2 shares a variable with a matched clause, and each such variable lies in <= 2
     other clauses, so there are <= 2 S_m of them (S_m = occurrences in matched clauses). Hence
     M_eff <= min(4 S_m, (2/3)(M - S_m)) <= min(4 S_m, (2/3)(3n - S_m)), maximised at
     S_m = 3n/7: M_eff <= 12n/7.
  2. Group sizes: units of <= 3 variables move together; Chernoff with weights <= 3 still gives
     |S_t,g| = (1 + o(1)) n/r.
  3. Balancing: the labels of non-free occurrences are (a_unit(v), a_unit(f)), pairs of unit
     schedules from different units. Colour the unit conflict graph (units sharing a clause;
     degree <= 18, so <= 19 colours) and fix the unit schedules colour by colour; each unit
     touches <= 18 edges. The 'designated' multi-stage Bernstein argument with b' = 18 and
     s = 19 stages gives D ln(76 D) <= delta^2 M_eff / 1024 (a few more bits of log2 D).
  Evidence: in the cloud, 10 cases (bounded and random formulas, plain and aggregated Taylor data)
  give exact counts with units; the head's reference implementation passes referee v4.
  Caveat: head-level sketch, not peer reviewed.

- chunk [AUDITED as a head proof sketch, 2026-10-31 15:00; unit structure, the (1 + 1/L) n bound
  and exact counts are checked by referee v5 whenever a score uses it]: units are connected chunks
  of the variable graph G (variables adjacent when they share a clause) instead of matched clauses.
  A component of G with fewer than L variables is one unit; a larger one is cut bottom-up along a
  BFS spanning tree (a vertex whose pending subtree reaches L vertices closes a unit; the leftover
  at the root, fewer than L vertices, joins a unit holding one of its tree children). Every unit is
  connected through clauses with >= 2 of its variables and, in a large component, has L to
  1 + 7(L - 1) variables (G has degree <= 6). Each clause is designated to the unit holding most of
  its variables, and all of those are read directly (constant rows, as in 'designated').
  1. Bound: let s_f = max_U |f cap U| - 1. A clause with >= 2 variables in U is internal to U and,
     for widths <= 3, internal to at most one unit. A unit connected through its internal clauses
     has sum (|f cap U| - 1) >= |U| - 1 over them (a connected hypergraph on |U| vertices). So in
     a large component C: M_eff(C) = sum_f (|f| - 1 - s_f) <= (M_C - |F_C|) - n_C + u_C, with
     M_C <= 3 n_C and |F_C| >= M_C / 3, hence M_C - |F_C| <= 2 n_C and M_eff(C) <= n_C + u_C
     <= (1 + 1/L) n_C. Small components contribute 0. Hence M_eff <= (1 + 1/L) n.
  2. Group sizes: units of <= 7L variables move together; Chernoff with weights <= 7L still gives
     |S_t,g| = (1 + o(1)) n/r at these n.
  3. Balancing: labels of non-designated occurrences are (a_unit(v), a_unit(f)), pairs of schedules
     of two different units that share a clause. Colour the unit conflict graph (degree <= 42L, so
     <= 42L + 1 colours) and fix the unit schedules colour by colour; one unit's choice moves at
     most 9 * 7L = 63L occurrences. The multi-stage Bernstein argument of 'designated' with
     b' = 63L and s = 42L + 1 stages gives a worst stage exponent >= delta^2 (M_eff/D) / (7224 L^2),
     so D ln(172 L D) <= delta^2 M_eff / (8192 L^2) suffices (about 29 bits of log2 D at L = 256).
  Evidence: the head's reference implementation passes referee v5 (exact counts with L = 2 and 3;
  M_eff below (1 + 1/L) n on random bounded formulas with n = 1200 to 1500 and L = 2 to 64, and on
  unions of Fano planes).
  Caveat: head-level sketch, not peer reviewed.

- modular [AUDITED as a head proof sketch, 2026-10-31 15:00; weighted counts mod small primes with
  residue-keyed fibres are checked by referee v5 whenever a score uses it; yes/no targets only]:
  decide satisfiability by S(w) = sum over satisfying x of prod_v w_v^(x_v) for w uniform in
  GF(p^k)^n, with a prime Q < p <= 2Q (Bertrand) and p^k >= 2n. S is multilinear with 0/1
  coefficients and degree <= n, so for a satisfiable formula S(w) = 0 with probability <= 1/2
  (Schwartz-Zippel), and S(w) = 0 always otherwise: one-sided error, amplified by repetition (PPSZ,
  the 1.307^n competitor, is randomized as well).
  1. The factor w_v^(x_v) depends on x_v only; at point t it joins the group of v, so the design
     conditions and deg_T H do not change.
  2. Every coefficient is p-integral: the rows are integers over Lden = (R-1)! Q!; the coefficients
     that express P_u(t) on S_h through C_h are (1 - mu_psi) psi_u / <psi, psi> with mu_psi in
     Z / (R-1)! and <psi, psi> a power of 2; Hermite interpolation needs the nodes 1..Q distinct
     mod p; the w_K weights are integers.
  3. With all coefficients in GF(p), a clause of group g sees another group h only through
     C_h x_h mod p. Keying the fibres X_t,g by residues still partitions {0,1}^S_g, and the
     identity of Lemma 3.2 (distributivity) holds over GF(p^k), so the sum over residue tuples of
     prod_g B_t,g is H_w(t). A measurement value costs log2 p <= log2(2Q) bits instead of
     log2(2^Q + 1).
  4. Exact counting cannot use this cheaply: recovering #F needs about n bits, i.e. CRT over primes
     up to about n, which costs about log2 n > Q bits per value.
  Evidence: the head's reference implementation passes referee v5 (weighted counts equal the frozen
  brute force mod 3, 5 and 7; an implementation with integer keys fails the fibre check).
  Caveat: head-level sketch, not peer reviewed.

- perpoint [AUDITED as a head proof sketch, 2026-10-31 16:40; Aux-level separation at every point,
  row counts and exact counts are checked by referee v6 whenever a score uses it]: Definition 3.1(3)
  lets C_{t,g} depend on t; the paper's Section 2.7 uses one C for every point.
  1. Free operators at t: a matrix X_t (x) Y (any Y on the other positions) whose position-t factor
     maps a pair (a, d) only to pairs (d, *) is supported on {(u, w): a_w(t) = d_u(t)}, i.e. it uses
     only the clause's own group at t. The swap tau_t (on unequal pairs) and the projector A_t on
     equal pairs are such factors.
  2. M0 = V0 f(W0) is diagonal in the product eigenbasis, so K = I - M0 = sum over the type of the
     position-t factor: (P^E_t + P^S_t) (x) m(W_rest) + P^N_t (x) m(W_rest + 1), m(j) = 1 - (-1)^j f(j).
     With P^E_t = A_t, P^S_t = (C_t + tau_t)/2, P^N_t = (C_t - tau_t)/2 (C_t the projector on unequal
     pairs), K is congruent modulo free operators to C_t (x) (m(W_rest) + m(W_rest + 1))/2 and hence to
     L_t = P^N_t (x) (m(W_rest) + m(W_rest + 1)). P(t) - L_t is supported on the own group (the
     paper's T-dependent part already is), so rows of L_t in rowspace(C_t) give Definition 3.1(4).
  3. rank L_t = r(r-1)/2 * #{rest basis vectors with m(j) + m(j+1) != 0} <= D p Pr[Bin(Q-1, p) outside
     [A0, A0 + R - 2]]. Given this P, it is optimal: on the pairs (a, d) with a > d the free part is
     strictly block triangular, so any representative has rank >= r(r-1)/2 rank(m(W)+m(W+1)).
  4. Rows are psi with an antisymmetric factor at t: entries {0, +-1}, exactly 2 groups at t,
     restricted absolute sums <= 2^(Q-1); h copies and the group split are as before.
  Evidence: exact checks (cloud, 6 cases): P(t) - L_t vanishes outside the free pattern (residual
  1e-16), rank L_t matches the formula (e.g. (r,Q,R,A0) = (2,4,3,0): 10 vs 13, (3,3,2,1): 135 vs 243);
  the head's reference implementation passes referee v6; a variant that drops m(j-1) fails the
  Aux-level separation check.
  Caveat: head-level sketch, not peer reviewed.

- eadapt [AUDITED as a head proof sketch, 2026-10-31 16:40; checked by referee v6 like perpoint]:
  M0 = V0 sum_e f_e(W0) [N_E = e], N_E = number of positions with a(s) = d(s). [N_E = e] is diagonal
  and commutes with every swap (a swap of an equal pair is the identity), so M0 is still a sum of
  partial swaps times diagonals: M0[u,w] != 0 needs w = swap_{S^c}(u) with |S| <= R - 1 and
  E_uw within S (Lemma 4.2 unchanged), and the rows keep degree <= R - 1. On a basis vector with e
  equal-pair factors and j antisymmetric factors M0 acts as (-1)^j f_e(j), so it is fixed when
  j lies in I_e; given e, j ~ Bin(Q - e, 1/2), and a centred window of R values is optimal.
  Evidence: exact checks (cloud, 5 cases): support violation 0, rank(I - M0) e.g. (2,4,3,0): 10 vs 13,
  (2,4,2,1): 46 vs 67, (3,3,2,1): 135 vs 189; these equal the best global min-rank values found by
  numerical search, which found nothing lower. Combined with perpoint: (2,4,2,1): 28, (2,5,3,1): 46.
  Caveat: head-level sketch, not peer reviewed.

- contract [AUDITED as a head proof sketch, 2026-10-31 17:05; the bond condition and exact counts
  with Zagg_contracted are checked by referee v7 whenever a score uses it]:
  1. Bond structure: every auxiliary measurement row (a basis vector psi) has a position-t factor
     supported on one unordered pair {a, b} of groups (or one group), so it touches at most 2 groups
     at t; after the substitution of Step 4 a row may meet fewer groups, never more.
  2. A clause of group g = d_u(t) with a replaced occurrence u (variable in group h = a_u(t)) has
     P_u(t) = L_t[u, :] + (own-group part), and L_t[u, :] is a combination of the psi with psi_u != 0,
     whose position-t factor sits on the pair {h, g}. So its form restricted to any group h' != g is 0
     unless h' = h, and then lies in the span of the rows of the pair {g, h} restricted to h. The
     T-dependent part uses only partial swaps of u (groups a_u(t), d_u(t)). Designated occurrences are
     read in their own group.
  3. Hence H(t) (and every restricted sum Z_{t,U}) is sum_x prod_g A_g(x_g, (pi_{h->g}(x_h))_h): a tensor
     network on the complete graph of the r groups. Eliminating the groups in a fixed order, the step
     that adds group k enumerates its assignments, the full bonds between the k - 1 earlier and the
     r - k later groups, the earlier halves of its bonds to earlier groups and the later halves of its
     bonds to later groups: 2(k-1)(r-k) + r - 1 half-bonds, against r(r-1) when every measurement tuple
     is enumerated. Each half-bond costs the same bits as before (values of its rows restricted to one
     group), so the measurement term is multiplied by phi(r) = max_k [2(k-1)(r-k) + r - 1] / (r(r-1)).
     The Taylor and group terms are unchanged; memory becomes exponential (time is what is compared).
  Evidence: the head's reference implementation passes referee v7 (exact and weighted counts with
  r = 2, 3, 4, with and without perpoint and eadapt; the bond condition holds on every design); a
  variant that does not enumerate the later half-bonds fails. Note: in small designs a pair row can
  meet only one group, so a half-bond must use the rows that touch g and no group other than g and h.
  Caveat: head-level sketch, not peer reviewed.

- smallprime [AUDITED as a head proof sketch, 2026-10-31 17:20; rows at every GF(p^2) node and the
  weighted count from Taylor data at the nodes are checked by referee v8 whenever a score uses it]:
  1. Where the prime enters: (R - 1)! must be invertible (the window polynomial f_I has denominator
     (R - 1)!), 2 must be invertible (the eigenprojectors psi psi^T / |psi|^2 have |psi|^2 a power of 2),
     and the Q interpolation nodes must be distinct and nonzero (Phi_uy(alpha_t) = prod_{s in S}
     (1 - alpha_t / alpha_s) vanishes exactly when t is in S). Integer nodes 1..Q force p > Q; nodes in
     GF(p^k) only need p odd and p > R - 1 (k = 2 when p^2 > Q).
  2. Lemma 4.3 over GF(p^2): for g != d_u(t), P_u(alpha_t) restricted to S_{t,g} equals (I - M0)
     restricted, a combination of the measured psi with p-integral coefficients (denominators are powers
     of 2 times divisors of (R - 1)!), so it lies in the GF(p^2)-span of the integer measurement rows
     reduced mod p; with perpoint the same holds for L_t = P^N_t (x) (m(W) + m(W + 1)).
  3. Taylor data: H(alpha_t + z) mod z^J uses Hasse derivatives (Taylor coefficients), which exist in
     every characteristic; the derivative terms are expanded in multilinear monomials of degree < J as
     before. Hermite interpolation at Q distinct nodes with J Taylor coefficients each is the Chinese
     remainder map for the pairwise coprime (T - alpha_t)^J, valid over any field.
  4. Yes/no: S(w) = sum over satisfying x of prod w_v^(x_v) is a nonzero multilinear polynomial of degree
     <= n when the formula is satisfiable; weights drawn from GF(p^m) with p^m >= 2n give error <= 1/2
     (Schwartz-Zippel). Field elements have O(log n) bits; only the keys (measured values mod p) enter
     the exponent: log2 p bits per measurement row instead of log2(2Q).
  Evidence: the head's reference implementation passes referee v8 (8 frozen designs, p = 3 and 5, with
  and without perpoint, eadapt and chunk units); three wrong variants (no division by (R - 1)!, factors
  1 - T alpha_s instead of 1 - T / alpha_s, constant identity rows) fail on every frozen design.
  Caveat: head-level sketch, not peer reviewed.

## Results that depend on a provisional regime
- 2026-10-31 12:14:25 iter 002 (T2): ppsz threshold 1,117,111 uses regime 'freedman' -> PROVISIONAL until the head audits that regime. [head: regime audited 12:05, this result stands as audited]
- 2026-10-31 12:28:56 iter 003 (T3): ppsz threshold 1,070,239 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 12:42:32 iter 004 (T4): ppsz threshold 438,437 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 13:05:42 iter 005 (T5): ppsz threshold 273,036 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 13:28:39 iter 006 (T6): ppsz threshold 193,841 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 14:10:18 iter 007 (T7): ppsz threshold 193,841 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 14:52:15 iter 008 (T8): ppsz threshold 57,882 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 15:24:03 iter 009 (T9): ppsz threshold 29,002 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 17:10:32 iter 011 (T10): ppsz threshold 25,699 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 17:36:43 iter 014 (T11): ppsz threshold 24,714 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 17:56:02 iter 015 (T13): ppsz threshold 22,118 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
