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

- edges [AUDITED as a head proof sketch, 2026-10-31 18:20; the two-group support of the rows and
  taylor_edges against a direct evaluation of the Taylor data are checked by referee v9 whenever a score uses it]:
  1. Two-group support: M0[u, y] != 0 only if y agrees with u or with u swapped at every position, so at point
     t the group a_y(t) is a_u(t) or d_u(t); hence every coefficient of P_u(T) is zero outside the clause's
     group d and the variable's group a (also on real designs: copies map labels to variables of the same group).
  2. On group a != d every column y has t in S(u, y), so Phi_uy(t + z) = (-z / t) Phi'_uy(t + z): the z-part
     B_i there has no z^0 term, while the constant part is the measured (I - M0) row (Lemma 4.3).
  3. Multilinearity: Phi_f(A + B) = sum over subsets E of the clause's edges of (prod_E B_i) d_E Phi_f(A), so
     H(t + z) mod z^J = sum over edge sets E with |E| < J of products over groups (A depends on its own group
     and on half-bond keys, B_i on group a_i only). Sequential elimination carries the bits of the edges
     between processed and unprocessed groups; the step for group k touches the edges of
     2(k - 1)(r - k) + 2(r - 1) ordered group pairs, each holding about m_eff n / r^2 edges, so at most
     sum_{j < J} C(E n, j) bit patterns with E = m_eff eps(r) (times 1 + delta for fluctuations); this
     replaces the sum over sets U (C(n, <= q)). The key states and the 2^{n/r} own assignments multiply as
     in 'contract'; the cost model adds the two maxima (an upper bound).
  Evidence: the head's reference Counter.taylor_edges equals the direct Taylor data at every point and gives
  exact and weighted (mod p, residue keys) counts on 63 designs (r = 2, 3, 4; with and without designation,
  chunk units, perpoint and eadapt); a variant without the edge terms fails. Caveat: head-level sketch,
  not peer reviewed.

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
- 2026-10-31 18:07:59 iter 016 (T12): ppsz threshold 18,817 uses regime 'tight' -> PROVISIONAL until the head audits that regime.
- 2026-10-31 20:13:04 iter 021 (T14b): ppsz threshold 10,086 uses regime 'tight' -> PROVISIONAL until the head audits that regime.

## External review notes (Gemini, 2026-10-31 21:40) and the head's checks
- edges: the per-pair link bound is deterministic given the label balance. At point t the projection
  u -> (a_u(t), d_u(t)) splits B into r^2 fibres of exactly D / r^2 labels, so with c_u <= (1 + delta) M_eff / D
  every ordered group pair carries at most (1 + delta) M_eff / r^2 replaced occurrences; the step for group k
  therefore touches at most (1 + delta) M_eff eps(r) edges. No adversarial formula can concentrate them.
- edges with chunk units: links between two chunks share one label and move as a block in the balancing;
  the chunk regime bound already pays for blocks (the 8192 L^2 term), so nothing changes.
- smallprime / modular: weights must come from an extension field GF(p^m) with p^m >= 2n (a satisfiable
  formula built from n/3 disjoint clauses would vanish mod p with weights in GF(p) alone). The modular and
  smallprime audits above already prescribe GF(p^k) weights with p^k >= 2n, so the stated error bound holds.
- Appendix D: Gemini's tile halving in D.4.5 does not hold (every tile of a packet with unequal
  symbols receives crossings from the other group). The first-pass cost model of 21:55 is superseded
  by the section "Appendix D route (tensor contraction): head extensions" below.

## Appendix D route (tensor contraction): head extensions
[AUDITED as head proof sketches, 2026-10-31 22:45; exact identities checked by
work/head/appD/appd_toy.py; cost model work/head/appD/appd_model2.py. Not run through the loop or the referee.]

Setting: Theorem D.1 on the clause-variable network of a 3-CNF in which every variable occurs at
most three times. A variable with d <= 3 occurrences is one equality vertex, so the network is
already bipartite (no subdivision): M <= 3n wires. Exponents are compared as in the main route
(polynomial factors ignored); memory may be exponential.

- appd-tight [accounting of the paper's own algorithm]:
  1. Multiplicities q_t = h max(0, u - i) with h = floor(M(2L - 1) / (L(2u - 2 + 2^(1-u)))) + 1 ~ M/(u - 1)
     already give sum_t q_t > M(2L - 1) (the paper takes h = ceil(3M/u) and counts odd points only).
     [Corrected 2026-11-01 after the GPT review: the earlier text said ceil(...), which gives only >=; the
     toy and the referee code always used floor(...) + 1.]
  2. Precision: -v2(E_t(0)) = h D_t sum_j (j/2^j) min(j, u - i) <= 6 h D_t (only 2 h D_t at i = u - 1).
     The evaluation loss (q_t - 1)(a + 1 - i) and the inherited precision fit under that slack once
     2^(a-u+1) >= (p_in + 4M(a+1))(u - 1)/(4M): a - u ~ 2 log2 u at the root, a_d - a_(d-1) - u_d ~
     log2 r below it. Hence p_t <= 6 h D_t (1 + o(1)) at the worst node (i = 0).
  3. Modes: (1 + delta) M / L packets with 2^(i+1) tiles each; a fraction 1 - 1/r of the packets has
     unequal symbols at coordinate i, so K = 2 (1 + delta)(1 - 1/r) M / D_t modes of each weight, and
     the paper's generating-function bound gives log2 N <= (2 pi / ln 2) sqrt(K p_t / 6).
  4. Sharper weights (optional; negligible in the model, large in toys): every crossing entry of
     W(t + alpha z) is divisible by 2^(a+1-i), because C(t + alpha z)/C'(0) is divisible by
     alpha = 2^(a+1) and v2(t - y0 + x) <= i; so w_Z = a + 1 - i + b is valid. The toy code asserts
     both weight rules coefficientwise.
- appd-bernstein [label balancing for Lemma D.3]: the main-route 'tight' argument (two-stage
  Bernstein with explicit product pessimistic estimators) applied to the vertex labels in [r]^a gives
  every pair-label load <= (1 + delta) M / D_arr once M / D_arr >= C Delta ln(D_arr) / delta^2
  (Delta = maximum degree: 3 at the root, 2L below it), and per coordinate and symbol at most
  (1 + delta) M / r real ports on each side. Padding every class to a multiple of L needs
  L D_arr <= delta M; at the break-even sizes this padding condition is the binding label constraint.
- appd-memo [every distinct child contracted once]: by D.4.6 the child of group g for tuple beta is
  built from the group-g vertices, the coefficients of the modes targeted at g and the powers
  L_Z(X)^beta_Z of the modes with source g. So Z_(g,beta) depends only on beta restricted to
  I_g = {Z : source(Z) = g or target(Z) = g}, and (D.16) reads
  F(t + alpha zeta) = sum_beta prod_g V_g[beta|I_g] (mod 2^p_t). Contract each table entry
  V_g[gamma] (gamma a tuple on I_g of weight < p_t) once. A call then costs (number of tuples) x r
  ring operations plus sum_g |V_g| child contractions, so along the recursion the exponents combine
  by maximum instead of by sum. The padding balances the modes of every ordered group pair exactly
  (each unordered coordinate pair {s, t} receives 2/r^2 of the packets, with 2^i tiles in each
  direction), so |V_g| <= exp(2 pi sqrt(K_g p_t / 6)) with K_g = (2(r - 1)/r^2) 2 (1 + delta) M / D_t.
- appd-modedp [mode sum contracted group by group]: V_g[gamma] is divisible by
  2^(sum over Z targeted at g of w_Z gamma_Z) (every use of Z by a row of group g carries a
  coefficient divisible by 2^w_Z), and every mode has exactly one target, so prod_g V_g[beta|I_g] is
  divisible by 2^w(beta). Any finite set of tuples that contains {w(beta) < p_t} therefore gives the
  same sum mod 2^p_t. Eliminate the groups in the order 1..r: Phi_k(beta|C_k) on the cut C_k (modes
  between {1..k} and {k+1..r}); step k+1 enumerates the tuples of weight < p_t on C_k plus N_(k+1)
  (modes between k+1 and {k+2..r}), a set that contains I_(k+1), multiplies by V_(k+1) and projects
  to C_(k+1). Every tuple of weight < p_t is counted exactly once along the chain; the extra tuples
  vanish mod 2^p_t. Step k+1 enumerates the modes of 2k(r - k - 1) + 2(r - 1) ordered group pairs,
  so its count uses K_dp = eps(r) 2 (1 + delta) M / D_t with
  eps(r) = max_k [2(k - 1)(r - k) + 2(r - 1)] / r^2 -> 1/2 (the main-route 'edges' factor).
- appd-chunks [only wires between chunks deformed]: build the main-route 'chunk' units (connected
  chunks of the variable graph with Lc to 1 + 7(Lc - 1) variables; each clause designated to the unit
  holding most of its variables) and merge every unit with its designated clauses into one vertex
  whose tensor sums over the unit's variables (2^O(Lc) work per entry). The remaining incidences
  (f, v) with v outside the unit of f are the wires; each has a row end (the unit of f) and a column
  end (the unit of v).
  1. The contraction is unchanged (sum the internal variables first).
  2. Bipartiteness is not needed at the root: use T(Y, X) = prod_v T_v(Y|v, X|v) in (D.8) instead of
     A(Y) B(X); F(0) is still the contraction. The locality of W(t) (D.4.5) uses only the slot
     orientation convention of the packets; the factorisation of (D.16) over groups needs
     T(Y, X) = prod_g T_g(Y_g, X_g), which holds because every vertex lies in one group. The children
     are bipartite (old vertices against packet vertices) as in the paper.
  3. M_root = number of crossing wires = M_eff <= (1 + 1/Lc) n (the main-route chunk bound).
  4. Labels: the pair label of a crossing wire is the pair of its super-vertex labels; the main-route
     chunk balancing (colour the unit conflict graph, fix the labels colour by colour, multi-stage
     Bernstein) needs D_arr ln(172 Lc D_arr) <= delta^2 M / (8192 Lc^2), a few dozen bits, far below
     the padding constraint L D_arr <= delta M at the break-even sizes.
  Children then have at most 2 (1 + delta) M_root / r wires.
Evidence (all four): appd_toy.py, seeds 1-4 in the cloud and seed 7 on this PC: plain and chunk
networks contract to #F (120 random formulas); tapered Hermite recovery from Taylor data reduced
mod 2^P_t (P_t = p + exact loss) returns F(0) mod 2^p in 60 cases with a = 2, 3 and u = 1..3, and
fails in 53 of them with 3 bits less; (D.16) computed per tuple, with memoised children and group
by group equals F(t + alpha z) from (D.8) mod (2^p, z^q) in 153 node checks (r = 3..5, a = 1..3,
up to 23 modes and 873 tuples, plain and chunk networks, both weight rules); the child networks of
D.4.6 contract to the table values; end to end with a = 1 the mode sums at the Lemma D.5 precision
give the exact count on 17 formulas. Not covered: recursion below the first level (children are
summed directly, the optimum in the cost model) and the label balancing (labels in the toys are
random). Caveat: head-level sketches, not peer reviewed.
Cost model (appd_model2.py), log2 n break-even vs 2^n / 1.637^n / 1.307^n / 1.1199^n:
  paper's constants 1.24e10 / 2.45e10 / 8.29e10 / 4.63e11 (its invariant needs log2 n >= 3.9e25);
  appd-tight with the paper's labels 43,460 / 93,569 / 359,349 / 2,269,327;
  + appd-bernstein 23,673 / 50,334 / 191,516 / 1,200,418;
  + appd-memo 11,420 / 25,203 / 97,912 / 643,153;
  + appd-modedp 7,700 / 15,010 / 58,640 / 375,303;
  + appd-chunks 780 / 1,527 / 5,393 / 34,371 (without appd-modedp: 972 / 2,182 / 8,561 / 60,973).
- Gemini review (2026-10-31 23:30) of appd-memo, appd-modedp, appd-chunks, the appd_model2 arithmetic and the
  sharper weights: all hold; it now agrees that its earlier tile halving does not hold. It proposed three more
  steps; the head checked them:
- appd-perclass [AUDITED as a head sketch, 23:45; the toy packing already uses it]: an unordered-pair pattern
  with a' unequal coordinates takes 2^(a-a') slots from each of its 2^a' orientation classes per packet, so its
  classes only need a common load that is a multiple of 2^(a-a'): at most one dummy packet (L slots) per pattern.
  With two-sided Bernstein loads the orientation classes differ by <= 2 delta M / D_arr, so the dummies total
  <= 2 delta M + L (r(r+1)/2)^a; the padding condition becomes (r(r+1))^a <= delta M instead of (2 r^2)^a.
  Packets per coordinate pair are balanced only to (1 + 2 delta) (the per-pattern rounding adds
  <= L (r(r+1)/2)^(a-1) slots per pair, a fraction < delta of 2M/r^2), and the per-symbol port bound comes
  directly from Bernstein instead of from the padding.
- appd-leafvars [AUDITED as a head sketch, 23:45; checked by byvars_value in appd_toy.py]: a leaf child of
  group g in a chunk network is summed over the variables x_g of its chunks. Fixing x_g fixes every column wire
  of the group (a column bit is the value of its variable). A designated clause then depends only on its
  external row bits: constant if its internal literals are true, otherwise 1 - prod_j [literal j false], two
  product terms (a clause always keeps at least one variable in its own chunk, so it has at most two external
  wires). After expanding these terms every row bit carries a unary weight, which moves into its packet factor:
  sum_y w(y) l_y(form) = w(0) + (w(1) - w(0)) form, so each packet tensor is still one coefficient extraction.
  Cost <= 2^(n_g + c2_g) with c2_g the group's clauses with two external wires, so the leaf exponent is
  <= (1 + delta)(n + c2)/r <= 1.5 (1 + delta)(1 + 1/Lc) n / r in the worst case (M_eff = n, c2 = M_eff/2),
  instead of 2 (1 + delta) M_eff / r edge bits. Evidence: equals the table value V_g[gamma] in 298 checks
  (seeds 11-13 in the cloud; seed 11 reproduced on this PC), including clauses with two external wires.
- Second recursion level (Gemini item 3) [model only, twolevel.py]: a depth-1 call inherits its parent's
  precision (a_1 >= a_0 + u_1 + log2 r_0 + O(1)) and has packet vertices of degree 2^(a_0+1), which costs
  label space at depth 1. With free r_0, r_1 and u_1 the model gains only against 1.1199^n: 28,442 (r_0 = 16,
  r_1 = 4, u_1 = 0.7 u_0) instead of 30,115 for one level; Gemini's estimate of 20,000 to 23,000 assumed
  u_1 of about 90, which makes the depth-1 tuple count far too large.
  Cost model v3 (appd_model3.py), log2 n break-even vs 2^n / 1.637^n / 1.307^n / 1.1199^n:
    appd-chunks + appd-perclass 734 / 1,396 / 4,856 / 31,498; + appd-leafvars 689 / 1,343 / 4,601 / 30,115;
    two levels against 1.1199^n 28,442.

## External review notes (Fable, 2026-10-31 23:40) and the head's checks
- Counting competitor [ACCEPTED]: for 3-CNF with every variable in at most three clauses the clause-variable
  graph has maximum degree 3, so its pathwidth is at most (1/6 + eps) of its vertices of degree 3 (Fomin and
  Hoie, IPL 97 (2006), Theorem 5 for maximum degree 3), at most n/3 + o(n) here; dynamic programming over a
  path decomposition counts models in 2^(pw) poly. So the fair counting competitor is 2^(n/3) ~ 1.26^n, not
  1.637^n or 2^n. Head re-runs against e = 1/3: loop counting version 43,518 (th_third.py; Fable 43,519);
  first audit: as written 1.66e12, re-tuned 5,710,562, tighter accounting 3,389,953; Appendix D rows
  1.11e11 / 496,177 / 264,399 / 143,651 / 78,641 / 7,705 (chunks) / 6,345 (perclass + leafvars; two levels
  give 6,461). The page now leads with the two class-specific competitors (2^(n/3) counting, 1.1199^n yes/no).
- Status of the paper [ACCEPTED]: not yet independently confirmed (review led by R. Williams with R. Impagliazzo
  and R. Paturi had not reported by 31 October; OpenAI's Lean check covers an unstated version). Every number
  is conditional on Theorem 1.1; stated on the page.
- Memory [ACCEPTED]: T12, T14 and the Appendix D memo/modedp steps use exponential memory (about 2^(0.36n)
  table entries at the break-even); PPSZ and the 1.1199^n algorithm use polynomial space, the 2^(n/3) DP
  exponential space. Stated on the page.
- '<= 1,000' cells [ACCEPTED as a caveat]: the model's regime conditions are not re-checked near log2 n = 1,000.
- Soundness of T1-T14b: no error found. Fable's re-runs if the least certain steps failed (yes/no vs 1.307^n):
  group sizes 2n/r 11,701; no contraction saving 11,880; no small prime 11,500; no edge bits 18,826.
- Appendix D (Fable saw only the first-pass model): agrees every tile carries modes and that Z_(g,beta) depends
  only on modes touching g (the appd-memo argument, found independently); sees no use of bipartiteness in
  D.4.3-D.4.6 beyond bookkeeping; the linear taper of Lemma D.5 is within 1.5% of the LP optimum at u ~ 700.
  Its three objections to the first-pass model (Bernstein packing, root-only chunk saving, smaller headroom)
  are answered by appd-bernstein, appd-chunks (children of 2(1 + delta) M_root / r wires) and appd-tight item 2.
- Dead ends Fable checked: re-tapering Lemma D.5, conjugate nodes over GF(p^2), CRT with small primes for exact
  counting, restricted per-position label alphabets. Sections 2-4: no large lever left in its view.

## Appendix D route: label-free headroom and recursion levels (head, 2026-11-01)
- appd-labelfree [AUDITED as a head sketch, 00:45; checked end to end on toys]: a call with u retained valuation
  levels evaluates only nodes t with v2(t) < u (q_t = h max(0, u - i) vanishes above), and the grouping at such a
  node reads only label coordinate i = v2(t) (D.4.5). The slot bits u .. a-1 serve only the precision headroom of
  Lemma D.5, so every vertex can take the same symbol there: those coordinates are equal at both ends of every
  wire and become copy bits of the packing (Lemma D.3 holds with them as equal coordinates). F(T), Lemma D.4 and
  the Hermite recovery do not depend on the labels. Labels then live in [r]^u: the Bernstein balance needs
  2 u log2 r <= log2 M - slack - log2(max degree) instead of 2 a log2 r, and the per-pattern padding is
  L (r(r+1)/2)^u <= delta M (with appd-perclass) instead of L (r(r+1)/2)^a. Evidence (work/head/appD,
  labelfree_check.py and labelfree_packets.py): exact counts end to end at a = 2, u = 1 with constant headroom
  labels in 7 formulas (up to 63 modes and 1.24 million enumerated mode combinations at one node); at a = 3,
  u = 1 the packing needs 579 packets instead of 1,141 with random headroom labels (197 formulas).
- appd-levels [model only, appd_model4.py and three_level.py]: with label-free headroom a deeper recursion pays.
  Break-even log2 n vs 2^n / 1.637^n / 1.307^n / 1.1199^n / 2^(n/3): 598 / 1,226 / 4,033 / 23,062 / 5,405
  (one level, r = 5 and 7, for the first two; two levels r = (7, 12), u_1 = 0.6 u_0 for 1.307^n and 2^(n/3);
  r = (8, 16), u_1 = 0.5 u_0 for 1.1199^n). A third level only changes the 1.1199^n column (22,884).

## Round 5: the Appendix D route in the loop (2026-11-01)
- Frozen referee referee_appd v1: appd_assembly.py (the cost model with levels 1-3 and the steps above, as a
  function of explicit parameters), appd_gates.py, appd_ref.py (the head's toy implementation, = appd_toy.py)
  and cases_appd.json. Score: the smallest log2 n at which the exponent drops below 1/3 (2^(n/3) counting),
  plus the four other targets.
- Code flags memo, modedp, chunks and leafvars enter the score only when the candidate's code reproduces the
  reference on hidden formulas (node values by the paper/memo/dp sums, child-call and enumeration counts,
  byvars on chunk cases, exact end-to-end counts). perclass and labelfree are head-audited model steps and are
  on from the start.
- Baseline (skeleton, no code flags): 2^(n/3) 201,458; 1.1199^n 1,475,104; 1.307^n 140,745; 1.637^n 34,975;
  2^n 16,738. The head's own implementation scores 5,404 / 23,056 / 4,032 / 1,213 / 579 with all four flags.
- Tasks: T15 (Qwen implements the core: group_table, sum_paper, sum_memo, sum_dp, byvars_value, node_value,
  count_e2e), T16 (Qwen widens the parameter search: delta, Lc, three levels, free u_d, wider r).
- iter_appd 001, T15 attempt 1 [ACCEPTED 01:29, head audit 01:45]: all four code flags verified on hidden cases
  (21 node values each by the paper, memo and dp sums, 91 byvars table entries on chunk networks, 10 exact
  end-to-end counts); 2^(n/3) 201,458 -> 5,404, 1.1199^n 23,056, 1.307^n 4,032, 1.637^n 1,213, 2^n 579 (equal to the
  head's own implementation). Audit of commit 09d787a: the sums follow the contracts; byvars_value is Qwen's own
  construction from the description (assignment of the group's variables, two-term split of the designated
  clauses not satisfied internally, unary row weights folded into the packet factors, coefficient of the target
  modes' tuple, source-mode linear forms raised to their exponents); no blocked function, no file access. The
  agent met its stop condition at about 25 minutes, kept running its own extra tests and was stopped by the
  55-minute wall; the referee scored the files it left.

## Appendix D: levers examined after round 5 and set aside (head, 2026-11-01 02:00)
Where the exponent sits: at every used node level K p_t / 6 = 2 (1 + delta) f M^2 / (u - 1) (independent of i;
corrected from 'K p_t = ...' after the GPT review, the model code was right), so a
call costs about (2 pi / ln 2) sqrt(2 (1 + delta) f / (u - 1)) M bits, and u is capped by the label budget,
about log2 M / (2 log2 r) (Bernstein over the r^(2u) wire types) or log2 M / log2(r(r + 1)) (padding). The binding
term is the root's mode sum (Qwen's T16 analysis finds the same). Ideas checked and not pursued:
- Low-cut partitions instead of random labels (fewer crossing modes, smaller children): the packing needs every
  orientation class of a pattern filled, i.e. crossing wires whose directions vary independently across the label
  coordinates; one partition reused (even relabelled) fills only 2 of 2^a' classes. Independent low-cut partitions
  would be needed per coordinate, and with chunks the wires (about n) are already the boundary incidences, so a
  k-level bisection (at most k n / 3 cut incidences for max-degree-3 incidence graphs) gives no smaller crossing
  fraction than random labels.
- Symbols that depend on lower slot bits (a trie instead of a product per packet): locality at node t still holds
  (row x and column x + t share bits below i), and the padding term improves from log2(r(r+1)/2) to log2(r^2/2)
  per coordinate (about 3% in log2 n), but a tile's target rows then lie in several groups, so modes need one copy
  per target group and K grows. Net loss.
- Per-row weights: rows with v2(t - y0 + x) = v gain i + 1 - v per moment order, but splitting the modes by row
  class multiplies the mode count (effective K grows by the harmonic number H_(i+1)). Finer tiles (step 2^(i+1+k))
  give 2^k (k + 1)^-1 times the effective K. No gain.
- Base-B slots (B = 3, 4) would cut the precision factor 6 to B(B + 1)/(B - 1)^2, but locality at a level needs
  col_sym(d + t_i) = row_sym(d) for every nonzero digit t_i, which forces one group unless every digit value gets
  its own label coordinate; the doubled label cost cancels the gain.
- Precision-aware children: V_g[gamma] is needed only mod 2^(p_t - w_src,g(gamma)), but the headroom a_1 - u_1 that
  absorbs the inherited precision changes by O(1) bits only.
- Decision via isolation (minimum-weight solution unique, S mod 2^(2 n^2)): the target precision is negligible
  against the recovery loss 6 h D, so yes/no costs the same as counting; no separate yes/no route.
- Leaves summed over variables below depth 1: the deeper children also contain packet vertices whose port bits
  remain, about the edge-bit count already charged.
Conclusion: after T16's parameter search, no further head step with a sound and checkable gain is in sight;
round 5 ends there.

- iter_appd 003, T16 attempt 2 [ACCEPTED 02:47, head audit 03:30]: a deterministic search (imports math and the
  public assembly only): shapes pre-ranked with an analytic replica of the model, the ten best refined with the
  real model (bisection on u_0, the other u_d maximised per level, a_d = min_a), 0.50 s per call. Best for 2^(n/3):
  two levels, r = (7, 10), u = (907, 642), a = (926, 1569), delta = 1e-4, Lc = 2^20 (delta and Lc sit at the edges of
  the model's allowed ranges; going further changes only the factors (1 + 2 delta) and 1 + 1/Lc, below 1%).
  Break-evens: 2^(n/3) 5,404 -> 5,300, 1.1199^n 23,056 -> 21,942, 1.307^n 4,032 -> 3,958, 1.637^n 1,213 -> 1,212,
  2^n 579 -> 572. Attempt 1 explored for the full 55 minutes without editing the file (rejected: no change);
  attempt 2 finished in 21 minutes with the original task text (the sharpened text prepared for a third attempt was
  not needed). The loop then found the backlog empty, stopped the Strata server it had started and exited.

## External review notes (GPT-6 Astra, 2026-11-01 15:30) and the head's checks
- Two documentation errors [ACCEPTED, fixed above]: the multiplicity must be floor(d/S) + 1 (the code always was),
  and the 'where the exponent sits' line meant K p_t / 6.
- The best scores use two recursion levels, which no gate executes [ACCEPTED]: head audit T19 planned.
- Newton (binomial) moments [checked by the head 15:45, work in progress]: G/(A - c j) = sum_(b<=j)
  [G c^b b! / prod_(l<=b)(A - c l)] binom(j, b) with c = 2^(i+1), A = T - y0 + x, so the source forms become
  sum_j binom(j, b) X_(y0+cj) and the coefficient of mode b has valuation >= a + 1 - i + b + v2(b!) (worst rows
  v2(A) = i; every factor l c / (A - c l) has valuation >= v2(l) + 1). Exact checks (gpt_checks.py): identity and
  bound coefficientwise mod z^q in 2,038 cases at a = 3..6, no failure; the extra valuation over the power rule
  is >= v2(b!) as predicted.
- Cauchy LDU of the crossing block [checked by the head 15:45, work in progress]: K_kl = 1/(T + x_k - x_l),
  x_k = 2^(i+1) floor(k/2^i) + (k mod 2^i). Closed-form pivots give v2(d_k) >= 2b + 2 v2(b!) - i (b = floor(k/2^i))
  before the factor G of valuation >= a + 1; in 570 exact factorisations (a = 3..6, z in {0, 1, 3, 5}) no pivot
  vanished, L and U were 2-adically integral and the diagonal bound held.
- What-if in the cost model (appd_model4, weight step s acts like cK = 2/s): break-even vs 2^(n/3) / 1.307^n /
  1.1199^n / 1.637^n: power basis 5,405 / 4,033 / 23,062 / 1,226; Newton (s = 2) 2,719 / 2,033 / 11,550 / 656;
  LDU (s = 4) 1,383 / 1,052 / 5,793 / 371. Not yet scored by a referee.

## Round 6 (1 November): binomial (Newton) mode basis [appd-newton]
- Lemma [head proof, suggested by GPT-6 Astra]. Node t with i = v2(t) < u, c = 2^(i+1), alpha = 2^(a+1),
  T = t + alpha z; a packet whose pattern differs at coordinate i; tile y0 with columns y_j = y0 + c j
  (0 <= j < L/c); a crossing row x, so A0 = t - y0 + x has e = v2(A0) <= i. With A = A0 + alpha z:
    W_(x, y_j)(T) = sum_(b = 0..j) N_(x,b)(T) binom(j, b),  N_(x,b) = W_(x, y0)(T) prod_(l = 1..b) c l / (A - c l),
  and coefficientwise in z, v2(N_(x,b)) >= v2(W_(x,y0)) + b + v2(b!), i.e. weight b + 1 + v2(b!) (paper) or
  a + 1 - i + b + v2(b!) (sharp).
  Proof: W_(x,y)(T) = G(T)/(T - (y - x)) with G = C(T)/C'(0), so W_(x, y_j) = G f(j), f(j) = 1/(A - c j). The divided
  difference of s -> 1/(A - s) at s = 0, c, ..., b c is 1/prod_(l<=b)(A - c l), so Delta^b f(0) = b! c^b /
  prod_(l<=b)(A - c l), and Newton's forward formula f(j) = sum_(b<=j) binom(j, b) Delta^b f(0) is exact for integers
  j >= 0. For l >= 1, v2(A0 - c l) = e because v2(c l) >= i + 1 > e; as a series in z, 1/(A - c l) has coefficients of
  valuation >= -e + k (a + 1 - e) >= -e, so c l/(A - c l) has coefficients of valuation >= i + 1 + v2(l) - e >= 1 + v2(l).
  Valuation bounds multiply under series products. QED.
- Consequence: the crossing part of every row factor is sum_b N_(x,b) L_b(X) with integer source forms
  L_b(X) = sum_(j>=b) binom(j, b) X_(y0 + c j) that read only the source group's column bits: everything (D.13)-(D.16),
  D.4.6, memo, modedp and the leaves need, with new weights. v2(b!) = b - s_2(b) >= b - log2(b + 1), so a tile has
  at most (W + log2(W + 1))/2 + 1 modes of weight <= w0 + W instead of W + 1. In the partition bound
  (2 pi/ln 2) sqrt(K p/6) for the mode tuples the density K halves: the model constant cK = 2 becomes cK = 1. Exact
  tuple counts (newton_partition.py): log N(Newton)/log N(power) = 0.738 at one tile and p = 4,000, falling
  towards 1/sqrt 2 = 0.707; the s_2 correction is lower order.
- Head checks (appd_toy.py with basis = 'newton', md5 61944a04): 178 node instances at a = 2..4 with both weight
  rules (D.12 and the weights asserted entrywise; node values by the three summation orders equal (D.8) computed
  directly; 1,172 child networks of D.4.6 and 980 by-variables leaves equal the tables); 88 end-to-end exact counts
  at a = 2 (u = 1, 2) and a = 3 (u = 1, 2), all correct; on the same instances the power basis needed a median
  36 times more mode tuples (largest: 5.9 million against 1.05 trillion). With basis = 'power' the toy's report is
  byte-identical to round 5.
- referee_appd version 2 (round 6): the assembly gains the code flag newton (cK = 1); the gates add 18 hidden node
  cases (a = 2..4, the candidate's Newton modes, weights, forms and coefficient residues must equal the reference;
  node values by the three orders with the reference child and enumeration counts; by-variables leaves) and 8 hidden
  end-to-end counts at a = 2, 3; these items decide only the newton flag. Version 1 results are reproduced exactly
  (third 5,300, occ3sat 21,942, ppsz 3,958, sharp3sat 1,212, brute 572 for the current code). The head's own
  Newton implementation verifies every flag and scores third 2,670, occ3sat 10,993, ppsz 1,998, sharp3sat 652,
  brute 320; a deliberately wrong candidate (Newton forms with the power weights: right values, wrong weights)
  loses the newton flag on 18 checks.
- T17 audit (head, 17:30): Qwen's change matches the specification; N_b is recomputed as the full product for each
  b (quadratic in b, same values), the D.12 self-check uses bfun, count_e2e takes a, u and basis; appd_params passes
  A.ck(flags) to its ranking replica. The agent met the stop condition within minutes, then spent the rest of the
  hour on extra stress runs and was stopped at the 55-minute wall clock; the referee verified every flag including
  newton and scored third 2,670, occ3sat 10,993, ppsz 1,998, sharp3sat 652, brute 320, equal to the head's own
  Newton implementation.

## Round 6: Cauchy LDU of the crossing block [appd-ldu]
- Lemma [head proof, suggested by GPT-6 Astra]. Node t with i = v2(t) < u, c = 2^(i+1), alpha = 2^(a+1),
  T = t + alpha z. For a packet whose pattern differs at coordinate i and a half h, order all L/2 slots with bit i = h
  as x_k = c floor(k/2^i) + (k mod 2^i) + h 2^i and write b_k = floor(k/2^i), r_k = k mod 2^i. The crossing block is
  W_(x_k, x_l)(T) = G(T) K_kl with K_kl = 1/(T + x_k - x_l), a Cauchy matrix (a_k = T + x_k, b_l = x_l), and
  K = L diag(d) U with
    d_j  = (1/T) prod_(m<j) -D_jm^2 / ((T + D_jm)(T - D_jm)),   D_jm = x_j - x_m,
    L_kj = T/(T + x_k - x_j) prod_(m<j) [(x_k - x_m)/(x_j - x_m)] [(T + x_j - x_m)/(T + x_k - x_m)]   (k > j),
    U_jl = T/(T + x_j - x_l) prod_(m<j) [(x_l - x_m)/(x_j - x_m)] [(T + x_m - x_j)/(T + x_m - x_l)]   (l > j).
  Coefficientwise in z: L, U are in Z_(2)[[z]] and v2(d_j) >= -i + 2 b_j + 2 v2(b_j!).
  Proof: v2(x_k - x_m) is < i when r_k != r_m and equals i + 1 + v2(b_k - b_m) when r_k = r_m, never i; so
  v2(t + D) = min(i, v2(D)) for every difference D of the block, and v2(D) - min(i, v2(D)) = E(k, m) :=
  [r_k = r_m] (1 + v2(b_k - b_m)). Hence v2(d_j(0)) = -i + 2 sum_(m<j) E(j, m) = -i + 2 sum_(b' < b_j) (1 + v2(b_j - b'))
  = -i + 2 b_j + 2 v2(b_j!), and v2(L_kj(0)) = [v2(T) - min(i, v2(x_k - x_j))] + sum_(m<j) (E(k, m) - E(j, m)) >= 0:
  for r_k = r_j the sum is v2(binom(b_k, b_j)) >= 0; for r_k != r_j the predecessors of j in class r_k include all
  b' < b_j, giving b_j + v2(b_k!/(b_k - b_j)!) >= b_j + v2(b_j!). U is the same with k replaced by l. Every factor
  (t + D) + alpha z has its constant term of valuation <= i < a + 1, so its inverse series has coefficients of
  valuation >= -v2(t + D); the bounds at z = 0 therefore hold for every z-coefficient. QED.
- Consequence: the crossing part of row k is G sum_(m<=k) L_km d_m (U X)_m, so mode m has the source form
  (U X)_m = sum_(l>=m) U_ml X_(x_l) (2-adic integer z-series coefficients, reading only the source group's column
  bits) and weight a + 1 - i + 2b + 2 v2(b!) with b = floor(m/2^i): the same multiplicities as the Newton tiles,
  twice the growth, so cK = 1/2. The factorisation is of the full block (all L/2 positions); factoring only the
  occupied sub-block loses integrality (a missing predecessor can make L non-integral), so sparse packets carry
  modes for unoccupied positions. Asymptotically the packets are full (padding <= delta M), so this costs nothing
  there; on the toys, whose packets are nearly empty, LDU can have more tuples than Newton at low precision. The
  series-valued source forms enter the child networks of D.4.6 as series-valued packet tensors, which the recursion
  already handles (see appd-twolevel).
- Head checks: ldu_closed.py (closed forms equal exact Doolittle over Q[[z]] mod z^q and every pivot valuation equals
  -i + 2b + 2 v2(b!), 76 blocks at a = 2..5, no mismatch); appd_toy.py basis = 'ldu' (md5 51d503b7): 71 node
  instances at a = 2..4 (D.12, the weights and the integrality of L and U asserted coefficientwise; node values by
  the three orders equal (D.8); 630 child networks and 519 by-variables leaves over the three bases), 63 end-to-end
  exact counts at a = 2, 3 (u = 1, 2), all correct; power and Newton outputs byte-identical. Tuple counts for the
  three weight rules on full tiles (ldu_partition.py): log N(LDU)/log N(Newton) = 0.720 at p = 4,000, falling
  towards 1/sqrt 2.
- referee_appd version 3: the newton items run a second time with basis = 'ldu' and sharp weights (same hidden
  cases); they decide only the new flag ldu (cK = 1/2). The head's LDU implementation verifies every flag and scores
  third 1,370, occ3sat 5,518, ppsz 1,044, sharp3sat 368, brute 197.

## Round 6: two recursion levels executed on toys [appd-twolevel, audit T19]
- appd_twolevel.py builds every child network of D.4.6 that a root node needs as an explicit network (the group's
  old vertices plus one packet vertex per packet, whose entries are the parent's z-series residues lifted to
  integers) and contracts it with the algorithm itself over R1 = Z[z1]/(z1^q1): its own labels, packing, Cauchy
  deformation, nodes, Newton modes and exact Hermite recovery, at the precision the parent node needs (inherited
  precision), with level-2 leaves by brute force. Lifting is sound because the child value is a polynomial with
  integer coefficients in the table entries, so lifted residues give a value congruent mod 2^P.
- Result (seeds 31-36, 90 trials, no failure): 76 formulas (root a = 1 or 2, u = 1 or 2, r = 2 or 3; child a = 1 or 2, r = 2 or 3;
  30 of them on chunk networks) counted exactly through both levels, 236 root nodes whose values from recursively
  computed children equal (D.8), and 38,030 level-1 contractions each equal to the brute-force child value in every
  z1-coefficient mod 2^P; child networks up to 5 wires, packet vertices of degree up to 4 (at most 2 L_parent), child
  precision up to 77 bits (root nodes up to 58). 14 formulas were skipped by a memory guard (a level-1 node above 200,000 mode tuples).
  This executes what the cost model's second level assumes (child wire count, packet-vertex degree 2 L_parent,
  inherited precision); it does not test the asymptotic label balance.
- T18 audit (head, 19:05): Qwen's change matches the specification (inline LDU branch in __init__, the asserts on
  unknown bases and on integrality, the D.12 branch with series products, series-aware packet_tensor and
  byvars_value through the helpers; appd_params sets ldu with the other code flags). It also rewrote both files
  with CRLF line endings (no functional effect; 85 lines added, 13 removed when whitespace is ignored). One attempt
  of 11 minutes; the referee verified every flag including ldu and scored third 1,370, occ3sat 5,518, ppsz 1,044,
  sharp3sat 368, brute 197, equal to the head's own LDU implementation.
- Round 6 summary: against 2^(n/3) the Appendix D break-even went 5,300 -> 2,670 (Newton, T17) -> 1,370 (Cauchy
  LDU, T18), each step a halving as the cost model predicts (log N of the mode tuples scales with sqrt(cK)). Both
  lemmas are head proofs with exact checks; the second recursion level is now executed on toys (appd-twolevel).

## Round 7 (1 November, evening): three outside reviews (ChatGPT, Gemini, Fable), shell matrix products [T20]
- appd-shellmm [head proof sketch, suggested by ChatGPT]. In sum_dp the step for group g computes
  Phi'[x + z] = sum_y Phi[x + y] T_g[y + z] over w(x) + w(y) + w(z) < p_t, where x runs over the tuples on the cut modes
  that do not touch g (processed group to later group), y over the cut modes that touch g and z over the modes from g
  to later groups; every cut tuple of weight < p_t is a state (induction on the steps), so the enumerated
  combinations are exactly the triples of total weight < p_t. Grouping x, y and z by exact weight (shells), each
  shell triple (w1, w2, w3) with w1 + w2 + w3 < p_t is one full block product A[X_w1 x Y_w2] B[Y_w2 x Z_w3] and the
  triples partition the combinations, so the sum is unchanged term by term. Strassen's identities use only ring
  additions, subtractions and products, so they hold in (Z/2^P)[z]/(z^q); with S the largest power of two
  <= min(n, k, l) an n x k x l block product takes ceil(n/S) ceil(k/S) ceil(l/S) 7^(log2 S - 1) 8 ring products,
  O(n k l S^(log2 7 - 3)). Cost: there are fewer than p_t^3 shells (a polynomial factor). For a shell with sizes
  N_x, N_y, N_z the cost is N_x N_y N_z S^(omega - 3) = max over j of N_x N_y N_z N_j^(omega - 3). The generating
  function bound N_V(w) <= Z_V(s_V) e^(s_V w), with s_V = s for the two larger sets and s/(omega - 2) for set j,
  gives log cost <= s p_t + (c/s)(A_tot - (1 - (omega - 2)^2) A_j) for every s (A_V proportional to the modes of V,
  the same constants as the tuple count), hence 2 sqrt(c p_t (A_tot - (1 - (omega - 2)^2) A_j)), largest for the
  smallest A_j. With group-pair counts (k-1)(r-k), k-1 and r-k for x, y and z the mode fraction becomes
  eps_mm(r) = max_k 2[(k-1)(r-k) + (k-1) + (r-k) - (1 - (omega - 2)^2) min(k-1, (k-1)(r-k), r-k)] / r^2
  (omega = log2 7: 0.5696 at r = 7 against eps_dp 0.6122). Among all merge trees (blocks a, b, rest c, fraction
  2[ab + ac + bc - (1 - (omega - 2)^2) min(ab, ac, bc)]/r^2) the one-group-at-a-time order is already optimal for
  omega = 3 and log2 7 at every r <= 16 (merge_tree.py); trees that cut the groups into thirds only help near
  omega = 2.37.
  Head checks: shellmm.py + shellmm_check.py in the cloud (seeds 1-5): 108 node instances x 3 bases (power, Newton,
  LDU), 5,435 block products, up to 28 modes: values equal (D.8) and sum_dp, block shapes equal the shell structure,
  the largest step's block work sum n k l equals the dp enumeration; Strassen exact on 8 shapes with 7^(k-1) 8 ring
  products at 2^k (19,208 instead of 32,768 at 32 x 32 x 32). The same code is mm_blocks.MM_CODE (reference and head
  solution on this PC, work/head/r7).
  What-if with Qwen's choose_params on the frozen model (whatif.py): Strassen 1,280 / 5,137 / 977 / 347 / 187
  (third / occ3sat / ppsz / sharp3sat / brute; before 1,370 / 5,518 / 1,044 / 368 / 197); with the best known
  exponent omega = 2.371339 (galactic algorithms, not an implementable step) 1,177 / 4,766 / 892 / 303 / 176 one group
  at a time, 1,194 / 4,740 / 913 / 330 / 171 with the best trees (the parameter search is a heuristic).
- referee_appd version 4: mm items decide only the new flag shellmm (eps_mm instead of eps_dp): matmul exact on 11
  shapes with at most 7^(k-1) 8 ring products (counted through ps_mul_mod) at 4, 8, 16, 32; node_value(method="mm")
  equals (D.8) on the 21 node cases (power) and 18 newton node cases (LDU), with the reference block shapes both in
  its stats and in every matmul call (wrapped by the referee) and the reference dp count; count_e2e(method="mm",
  basis="ldu") exact on the 8 newton end-to-end cases. Dry runs: the head's solution verifies every flag and scores
  1,280 / 5,137 / 977 / 347 / 187 (= the what-if); the round-6 code keeps every old flag and its scores with shellmm
  off; a schoolbook matmul and a sum_mm that bypasses matmul both fail shellmm.
- appd-ldu finite bound [head proof, suggested by ChatGPT]. For k >= 0, B = 2^k, w_b = k + 2 + 4b - 2 s2(b) (b < B,
  s2 = binary digit sum; these are the LDU weights a + 1 - i + 2b + 2 v2(b!) with k = a - i - 1) and every s > 0:
  sum_b -log(1 - e^(-s w_b)) <= pi^2/(24 s). Proof: (1) W_k = {w_b} is majorized by V_k = {4b + 2}: by the top bit of
  b, W_k = (W_(k-1) + 1) u (W_(k-1) + 2^(k+1) - 1) and V_k = V_(k-1) u (V_(k-1) + 2^(k+1)), in both the lower half
  lies below the upper half (max W_(k-1) + 1 = 2^(k+1) - k < 2^(k+1) + k), so the sums of the j smallest elements
  satisfy P_W(j) - P_V(j) >= j for j <= H = 2^(k-1) and >= H - t >= 0 for j = H + t (induction), with equal totals
  2B^2. (2) g(x) = -log(1 - e^(-s x)) is convex, so Karamata gives sum g(w_b) <= sum g(4b + 2). (3) For convex g the
  midpoint value is below the mean over [4b, 4b + 4], so sum_b g(4b + 2) <= (1/4) int_0^inf g = pi^2/(24 s). QED.
  Consequence: a node with C full chains has at most exp(s p + C pi^2/(24 s)) tuples of weight < p for every s,
  i.e. log N <= pi sqrt(C p / 6), which is exactly the frozen model's count with cK = 1/2: the LDU cost is a
  finite-size upper bound, not only an asymptotic estimate (modes dropped at weight >= p or with empty forms only
  lower the count). Checks: majorization and the inequality for k <= 22 on 60 values of s (ldu_bound.py).
- Examined and set aside (with reasons):
  * Fable, half-packet row cap (the exponents of a packet half sum to at most its rows, since every row factor is
    affine in the modes; Fable estimated -4%): valid, but Chernoff with a second multiplier at the model's saddle
    (s 2^k = 2.19, D_t = 2^(a-i)) shows the uncapped mean exponent per mode exceeds 1 only for chains of about
    2^18 modes or more, and the cap costs a multiplier nu ~ e^(-4 s 2^k) ~ 1.6e-4: log N(capped)/log N = 0.99996 in the long-chain
    limit (rowcap2.py). The large exponents sit on a few light modes that carry almost no entropy. Effect < 0.01%.
  * Fable, balance constants 8192 Lc^2 and 172 halved (-1%) and partial packets (-2%, 1.1199^n column only): small;
    the first needs a sharper balancing proof. Not pursued.
  * Gemini 1 (fold the roots of C(T) into L and U for 6b-8b weight growth): impossible for this deformation. At z = 0
    any factorisation K = P diag(d) Q with P, Q invertible over Z_(2) has the Smith invariants of K(0) as its diagonal
    valuations; the LDU attains them (sum = v2 det K(0)), so no basis of the same block has heavier modes, and G is a
    common scalar factor that already adds a + 1 - i to every mode. A further halving of cK needs a different
    deformation W(T) (Fable's open question below), not a different factorisation.
  * Gemini 2 (pair valuation levels to cut the factor 6 to about 2): the 6 = sum j^2/2^j comes from the node counts per
    valuation level (half as many at each higher level), not from how levels are grouped; every multiplicity profile
    is inside the taper LP, where the linear taper of Lemma D.5 is within 1.5% of the optimum (ledger, round 5).
    No mechanism.
  * Gemini 3 (drop the modes of empty slots): modes with an empty source form are already dropped, and the padding
    keeps the packets full (dummies <= delta M). Already in the model.
  * Gemini 4 (general 3-SAT through sparsification and Tovey copying, 15,000-25,000): copying multiplies the wires and
    the sparsification density bound is far too large to give a meaningful number; the general-density entry below is
    the direct version.
- Open question (Fable, recorded as stated): can node values that are local at their own level make the 2-adic
  determinant excess of the crossing blocks larger than (L/2) D_t per half? Cauchy reaches it through the
  Vandermonde-square structure; XOR shifts, anti-shifts, inversion, free node values at unused levels, sparse packets
  and extra nodes all lose (Fable); no proof that Cauchy is optimal. By the Smith-invariant remark above, this is the
  only place a further halving of cK can come from. Fable's calibration: with labels, children and the cut free the
  exponent is still about 4.5/sqrt(lambda), a floor near log2 n = 200 against 2^(n/3).
- appd-density [model]: 3-CNF with m clauses and no occurrence bound. The chunk bound (appd-chunks item 1) gives
  M_eff <= sum_f (|f| - 1) - n + u <= 2m - n + n/Lc. The balancing constants assume at most three occurrences, so
  first replace each variable with more than B = 1024 occurrences by a chain of copies joined by equality constraints
  of arity 2: each copy adds one variable and one constraint with |f| - 1 = 1, so M_eff does not grow; the number of
  variables (and with it the leaf term) grows by at most 3m/(B - 2), and the balancing constant 8192 Lc^2 by (B/3)^4
  (unit size and conflict degree scale with B), about 34 bits. At m/n = 4.26 (the random 3-SAT threshold),
  conservatively also adding 3m/(B - 2) wires, the frozen v3 model with Qwen's choose_params gives 1.307^n 55,602,
  1.637^n 16,439, 2^n 8,334 (without the occurrence cap and its extra bits 55,416 / 16,381 / 8,301, matching Fable's
  55,420 / 16,384); with shellmm 51,728 / 15,295 / 7,755 (whatif3.py). The 2^(n/3) and 1.1199^n competitors are
  specific to three occurrences.
- T20 audit (head, 21:05): Qwen's change matches the specification: ps_sub_mod, naive_matmul (it also skips zero
  entries, which only lowers the product count), strassen_sq with Strassen's seven products (C22 written as
  M1 + M6 - M2 + M3), matmul on square blocks whose side is the largest power of two below the smallest dimension,
  with zero padding, NodeExpansion.shells and sum_mm as specified (through the module-level matmul), node_value
  method "mm" and count_e2e(method=...); appd_params passes shellmm to rank_e and uses A.eps_mode. One attempt of
  8 minutes (47 requests). The referee verified every flag including shellmm (11 matmul, 39 node and 8 end-to-end
  items, verify 148 s) and scored third 1,280, occ3sat 5,137, ppsz 977, sharp3sat 347, brute 187, equal to the
  head's own solution and to the what-if.
- Round 7 summary: against 2^(n/3) the Appendix D break-even went 1,370 -> 1,280 (-6.6%, shell matrix products with
  Strassen, T20); the other columns fell by 5 to 7%. With the best known matrix multiplication exponent the model
  would give about 1,180 (not implementable). The other levers in the three reviews are below 2% each or already in
  the model; a further halving needs a different deformation (the open question above).

## Round 8 (1 November, night): the guarded Lagrange deformation [appd-guarded, T21]
- Context: asked about Fable's open question (a better deformation than Cauchy), GPT-6 Astra proposed changing the
  global deformation while keeping the Cauchy crossing blocks: the Smith-invariant remark of round 7 fixes the mode
  weights of a given block, not the degree of the deformation, and the degree sets the Hermite precision.
- appd-guarded [head proof sketch of Astra's construction]. A packet is a virtual grid {0..Q-1}, Q = 2^a; real ports
  sit only on the central slots H <= x < Q - H, H = Q/2^kappa (kappa >= 2); labels live on the evaluated levels
  0..u-1 with u <= a - kappa, so H is a multiple of 2^u. Deformation W~_xy(T) = l_y(x + T), l_y the Lagrange basis
  polynomial of y on the grid.
  1. W~(0) = I and every entry has degree Q - 1, so F~(T) has degree <= M(Q - 1) instead of M(2Q - 1).
  2. Locality: for t in [-H, H) and a real row x, x + t is on the grid, so W~(t) is the shift x -> x + t; with
     i = v2(t) < u the two slots differ in bit i and the packet's orientation convention puts both ports in one
     group; a row whose shift lands on a guard or an empty slot has no real column at t.
  3. Integrality at T = t + alpha z, alpha = 2Q: every P_Q'(y) = +-y!(Q-1-y)! has valuation v = Q - 1 - a (row
     Q - 1 of Pascal's triangle is odd); with j = x + t the factor m = j of l_y(j + alpha z) is alpha z and the
     others have content 2^v2(j - m), so an off-shift entry has content >= 2^(a + 1 - v2(j - y)) (as the Cauchy
     entry) and the shift entry is 1 + O(2^(a+1-v2)).
  4. Modes: W~_xy = R_x K_xy c_y with the round-6 Cauchy kernel K_xy = 1/(t + alpha z + x - y), the row factor
     R_x = 2^-v P_Q(x + t + alpha z) (content >= 2^(a+1)) and the odd unit c_y = 2^v/P_Q'(y). The full-virtual-half
     LDU of K is unchanged (alpha is still 2Q), so the modes (U diag(c) X)_m with coefficients R_k d_m L_km have the
     round-6 weights a + 1 - i + 2b + 2 v2(b!), and the round-7 finite bound applies to every full chain.
  5. Recovery: the node grid [-H, H) has H/2^i nodes at level i, complete residue systems; with q_t = h(u - i) and
     h = floor(M(Q - 1)/sum_i (H >> i)(u - i)) + 1 the total multiplicity exceeds deg F~, and Lemma D.5's argument
     with Q replaced by H gives -v2(E_t(0)) <= 6 h H/2^i (hermite_loss.py: max ratio 0.992 for H <= 2^12, all u).
     As h ~ (Q/H) M/(2(u - 1)), the main loss is about 3MQ/((u - 1) 2^i): half of the paper's 6MQ/((u - 1) 2^i).
  6. Packing: [H, Q - H) is a union of whole periods of the low u bits, so every orientation class of a pattern with
     a' unequal coordinates gets N/2^a' slots, N = Q - 2H; padding and balance hold with capacity N <= Q, and child
     packet vertices have 2N <= 2Q ports.
  Cost per level: K p scales by gamma = (Q/N)(Q - 1)/(2Q - 1) -> 1/(2(1 - 2^(1-kappa))) (0.501 at kappa = 10); the
  multiplicity grows by about 2^(kappa-1), so the headroom needs 2^(a-u+1) >= (p_in + 2^(kappa+1) M (a+1))(u - 1)/(2M)
  (about kappa more bits of a at the root, 1 below it) and the children inherit 3h 2^a instead of 6h 2^a.
- Head checks: Astra's program reproduced (27 tensor contractions through 108 Taylor jets with exact recovery mod 2^9;
  8 full-block LDUs, 9,472 coefficients, 656 crossing entries). guarded.py on appd_toy (cloud, seeds 1-5): 75 node
  instances at a = 3, 4, 5, kappa = 2, 3, u = 1, 2 (up to 17 modes): the shift at every node, integrality, row
  factors >= 2^(a+1), the weights and (D.12) asserted, node values by paper / memo / dp / shell-mm equal (D.8)
  with the new deformation, 212 child networks of D.4.6 and 93 by-variables leaves equal the tables; 41 end-to-end
  exact counts (29 with wires), and the guarded recovery used 26 to 49% less precision than the paper's on the
  same formulas (e.g. 110 against 214 bits).
- What-if (whatif_guard.py: frozen v4 model, Qwen's T20 search with the guarded formulas): kappa 3: 883 / 3,438 /
  680 / 261 / 145; 4: 769; 6: 705; 8: 690 / 2,611 / 539 / 232 / 126; 10: 686 / 2,600 / 531 / 204 / 126; 12: 682;
  14: 686 (third / occ3sat / ppsz / sharp3sat / brute; before 1,280 / 5,137 / 977 / 347 / 187).
- referee_appd version 5: flag guarded with kappa per level (needs ldu; u <= a - kappa); gamma_guard, need_bits and
  next_lp in the assembly. Hidden cases v3: 21 guarded node cases (a = 3..5, kappa = 2, 3, u = 1, 2, at least two
  modes) and 8 guarded end-to-end cases. Items: the agent's NodeExpansion with the reference GuardedPacking has the
  reference modes, series source forms and coefficient residues; node_value (dp, mm) equals (D.8) with the guarded
  deformation and the reference dp count; byvars on chunk networks; count_e2e(guarded=True) exact. Dry runs: the
  head's solution verifies every flag (verify 143 s) and scores 686 / 2,600 / 531 / 204 / 126 (= the what-if at
  kappa = 10); the round-7 code keeps its flags and scores with guarded off; a single row factor for all rows and
  missing column units both fail guarded.
- T21 (iter_appd 007, attempt 1, 22:07 to 22:49, 38 requests): Qwen's smoke test printed all five OK lines at about
  22:21; it then kept running its own extra scratch tests (no fault found; some of its larger scratch cases ran out of
  memory in the existing tuple enumeration), so the head stopped the agent at 22:49 instead of waiting for the
  55-minute cap, and the loop ran the referee as usual. Referee: accepted at 22:53, every hidden check passed (verify
  179 s, among them 21 gd_modes, 21 gd_dp, 21 gd_mm, 42 gd_byvars and 8 gd_e2e items; all eight code flags verified).
- T21 audit (git diff -w 9d1006a HEAD -- algo_appd): exactly spec items 1-4. appd_core.py: NodeExpansion reads
  pk.guarded and asserts the LDU basis with the sharp weights for a guarded packing; it builds the W~ series from
  pk.Wfrac only for pairs of real slots; in the LDU branch it skips Cpoly and G, multiplies the U coefficients by
  pk.colunit(xk[k]) (integrality asserted on the product) and uses the row factor pk.rowfactor(xk[k], t, alpha, q) in
  place of G for the coefficient of mode m in row k; positions, weights, weight asserts and the (D.12) check are
  unchanged. direct_taylor_mod uses pk.Wfrac for a guarded packing. count_e2e(..., guarded=False, kappa=2) takes
  labels of u coordinates, GuardedPacking and guarded_multiplicities(pk, guarded_h(M, pk)); precision, node expansion,
  summation, recovery and the 2^unused factor are unchanged. appd_params.py: GUARD_KAPPA = 10 and _KAPPA (set in
  choose_params from the guarded flag); min_a starts at u + kappa and tests A.need_bits; every inherited precision
  (rank_e, make, build, full_ok) is A.next_lp; rank_e scales cK by A.gamma_guard for both the arithmetic and the child
  term; make adds kappa per level. One change beyond the spec, and a needed one: the fallback of choose_params now
  comes from make([2], [2], 1e-4, 4096, flags) instead of a fixed dict with a = 3, which would break u <= a - kappa
  under the guarded flag. GuardedPacking and the head helpers are untouched; no other file changed.
- Round 8 summary: GPT-6 Astra's guarded Lagrange deformation, proved and checked by the head, implemented by Qwen in
  one attempt and verified by referee_appd version 5: third 1,280 -> 686, occ3sat 5,137 -> 2,600, ppsz 977 -> 531,
  sharp3sat 347 -> 204, brute 187 -> 126 (-46% against 2^(n/3)), exactly the what-if at kappa = 10 and the head's dry
  run. At clause density 4.26 (the round-7 tree model with round 8's steps, whatif_dens_v5.py, kappa = 10): 25,946
  against 1.307^n, 7,689 against 1.637^n and 3,911 against 2^n (round 7: 51,728 / 15,295 / 7,755).

## Round 9 (2 November, night): Astra's precision reset and multiplicity balancing checked; a better search [T22]
- Context: asked once more, Fable found nothing new; GPT-6 Astra proposed two changes to the interfaces between
  recursive calls (note "Precision reset and multiplicity-sensitive label balancing", 1 November): (a) recover each
  child's high input precision from many low-precision contractions by CRT, (b) charge the parallel multiplicity
  instead of the maximum degree in the label balance of the bipartite children.
- appd-reset [head check of Astra's construction; model only, no code]: fix the canonical integer lifts of a child's
  entries at its input precision P (the creation-time representatives of D.4.1, never re-evaluated modulo q), reduce
  their coordinates modulo small odd primes q, contract each reduced network with the same 2-adic algorithm at a
  precision B_q with 2^(B_q - 1) > H_q = 2^M C^(nu - 1) (d (q - 1))^nu, lift centered (exact integer Z_q), and
  recombine by CRT until prod q > 2 H_P, H_P = 2^M C^(nu - 1) (d (2^P - 1))^nu. Exact, because reduction modulo q is
  a ring map of Z[U, V]/(f, V^l) and the contraction is a polynomial in the entries; polynomially many calls. In the
  model a child then inherits log2(B_q / n) ~ log2(M_child) + 16 (nu <= 2M, d = 4 lam, log2 C <= (4 lam - 1)
  log2(4 lam)) instead of a + log2(3h), so a_d ~ u_d + 23 at every depth. Astra's program rerun by the head (cloud): 18 CRT
  cases (1,494 recovery contractions) and 6 digit-polynomial cases (1,492), 16,384 packet-multiplicity labelings,
  no failure, every oracle precision as reported.
- appd-multbal [head check of Astra's lemma; model only]: in a bipartite child (old vertices against packet
  vertices) label the old side first: Bernstein for the per-label degree mass (summands <= Delta_L) and a k-collision
  cap on distinct neighbours of one label; then the packet side, with pair loads X_ab = sum_v N_v(a) 1{l(v) = b},
  summands N_v(a) < mu k, variance <= mu k (1 + eps) T: pair loads <= (1 + delta) M/S^2 under Astra's (L) and (R).
  Packet patterns fix the other endpoint's label, so the multiplicity of the next child is below mu k; mu_1 <= 63 Lc
  (the root unit degree, no cap at the root).
- What-if (r9/model9.py: the frozen v5 model re-implemented, equal to the assembly on the current certificate and
  on four new ones; r9/search9.py: greedy largest u per level with caps, shapes r_0 <= 12, r_1 <= 24, r_2 <= 64,
  delta 3e-4 to 1e-2, Lc = 2^6 to 2^20, k in 4..12; run on 22 cores): the current model gives 661 / 2,600 / 508 /
  190 / 117 (third / occ3sat / ppsz / sharp3sat / brute); with reset, multbal or both the same, except 189 against
  1.637^n with both. Reason: the binding constraints are all at the root. Its arithmetic needs u_0 ~ 108 at
  r_0 = 7, and then the root's chunk label balance (log2(8192 Lc^2) + 2 log2(1/delta) + log2(ln(172 Lc) + lnD),
  about 54 bits) and its padding (a_0 = u_0 + 23) fix lam; the children have slack, which is what the two steps
  enlarge. Even with 20 bits less root slack the break-even only moves to 654 (then the root padding binds), again
  equal with and without the two steps. Not adopted; the next levers are at the root (its balancing constant and,
  above all, the arithmetic constant cK gamma eps_mm: 1% less moves the 2^(n/3) break-even by about 6).
- The what-if also showed that Qwen's search missed the optimum of the current model: it tried only Lc = 2^12 and
  2^20 and pre-ranked shapes by the padding budget alone, which ignores the root balance (so Lc = 2^20 always won
  the pre-rank). T22: est_u also respects the assembly's balance test at each depth, a finer Lc and delta grid,
  20 refined candidates; no new math, referee v5 unchanged.
- Head dry run of the T22 search (r9/head_t22, the referee's own score gate): 663 / 2,600 / 508 / 190 / 114, 2.3 s per
  choose_params call.
- T22 (iter_appd 008, attempt 1, 00:27 to 00:33, 21 requests; check_t22.py printed T22 OK): accepted at 00:39, every
  hidden check passed (verify 187 s, all eight code flags verified), third threshold 686 -> 663; score gate 74 calls,
  2.5 s each.
- T22 audit (git diff -w 177fedd HEAD -- algo_appd): exactly spec items 1-3. appd_params.py only: the finer DELTAS,
  LCS and TOPK = 20; est_u renamed est_u_pad; bal_u copies the frozen balance test (root chunk slack, degree 2^(a_prev
  + 1) below) with the fixed-point start and a step-down; est_u is the minimum of the two; d, a_prev and Lc are
  threaded into rank_e, build and full_ok; a docstring paragraph. No other change (one unused local in bal_u).
- Round 9 summary: no new step adopted (Astra's two steps hold but do not bind); the better search gives 686 -> 663
  against 2^(n/3), 531 -> 508 against 1.307^n, 204 -> 190 against 1.637^n and 126 -> 114 against 2^n, 2,600 against
  1.1199^n unchanged, exactly the head's dry run. At clause density 4.26 the T22 search keeps 25,946 / 7,689 / 3,911
  (whatif_dens_v6.py: the optimum there already uses Lc = 2^20).

## Round 10 (2 November): binary crossing with fast integral absorption (GPT-6 Astra) [T23, T24, T25]
- appd-bincross [head-verified sketch; proposed by GPT-6 Astra, "Fast integral absorption of binary crossing
  pieces", 2 November, after the head rejected the runtime claim of the Schur crossing pieces].
  1. In the mode sum of a node the modes of the group pair {g1, g3} (edge x) and of {g2, g4} (edge y) are equality
     wires between the group tables. For fixed shells (exact weights wx, wy) the tuples get codes of a and b bits
     (absent codes are zero entries); on the top m = min(a, b) bit pairs the two wires are SWAP^(m) (digit
     2 x_k + y_k -> 2 u_k + v_k with u = y, v = x); the low |a - b| bits of the longer code pass through.
  2. SWAP = I - J, J = j j^T, j = (0, 1, -1, 0); SWAP^(m) = sum_{i=0..m} A_i with A_i = [t^i] (I - tJ)^(m).
  3. Row factorisation: S unimodular with first column j, A_i = (S^(m) E_i)(E_i^T (S^-1)^(m) A_i), E_i = words with
     >= i zero digits, R_i = sum_{k <= m-i} C(m, k) 3^k channels. Column flattening: colflat(I - tJ) = K (P0 - tI)
     with K^2 = I, K P0 = P0; T unimodular with first column (1, 0, 0, 1): colflat(A_i) =
     ((K T)^(m) F_i)(F_i^T (T^-1)^(m) B_i), F_i = words with >= m - i zero digits, C_i = sum_{k <= i} C(m, k) 3^k.
     Column pairing for 2i <= m, row pairing otherwise: every piece has at most sum_{k <= m/2} C(m, k) 3^k <= 2^(c m)
     channels, c = log2(12)/2 = 1.7925.
  4. Fast application: (S^T)^(m), ((K T)^T)^(m), (S^-1 - t S^-1 J)^(m) and (T^-1 P0 - t T^-1)^(m) act digit by
     digit: O(m^2 4^m) ring operations per assignment of the other indices of the old pair table, entries in
     {-1, 0, 1}, no division.
  5. The new index: log2 |theta_i| <= max(a, b) + (c - 1) min(a, b) <= beta sqrt(u_theta (wx + wy)) + O(1),
     u_theta = 2 (1 + (c - 1)^2) = 3.2561 mode units instead of 4 (Cauchy-Schwarz on a <= beta sqrt(2 wx),
     b <= beta sqrt(2 wy)); theta is graded by the original weight wx + wy and every cutoff uses original weights.
  6. Cost in mode units (2 per group pair; shell/Strassen merges q = X + Y + Z - eta min(X, Y, Z)): pair tables
     4r - 6 - 2 eta (prepare) and 4r - 8 (absorb); then pair A, the r - 4 singletons one at a time and pair B last,
     pair-pair edge 4 + u_theta: eps_bc1(r), equal to the exact best merge tree for r <= 16. Old -> new: r = 5
     0.5843 -> 0.5545, 6 0.5724 -> 0.5518, 7 0.5696 -> 0.5544, 8 0.5611 -> 0.5495; no gain at r = 4 (preparing a
     pair already costs the old peak).
  Head checks. Cloud: verify_binary.py (SWAP identity, both factorisations and their supports for m <= 4, ranks over
  Q (mod p) equal to R_i and C_i, fast = dense for m <= 5 with O(m^2 4^m) operation counts); verify_k6.py (216 exact
  K6 contractions with padded unequal alphabets: row-only, column-only and mixed pieces); crossing_dp.py (exact
  merge-tree DP: reproduces eps_mm and Astra's table for r = 4..10). PC: bc_head.py on 600+ random nodes of the toy
  algorithm (4 to 6 groups, power and LDU bases, chunks): total = (D.8) = sum_paper, per-piece values identical for
  row-only, column-only and mixed pairings and equal to a dense evaluation with the explicit A_i. Astra's checker
  and cost probe reproduced.
  What-if with the T22 search on the patched model: 646 / 2,531 / 497 / 186 / 112 (from 663 / 2,600 / 508 / 190 /
  114). Several disjoint crossings (k <= r/4, exact DP) lower eps only at r >= 9 and do not change the five values.
  Caveat: head-level sketch, not peer reviewed; the bound lives in the frozen unit model.
- referee_appd version 6: code flag bincross (needs shellmm; mode fraction min(eps_mm, eps_bc1)); non-scoring flag
  bc_fast. Gates: bc_fast (bc_mask, bc_left, bc_right equal the reference on 20 random series vectors, m = 0..4, both
  pairings), bc_node (15 hidden nodes with 4 to 6 groups, up to 129 pieces: the value and every per-piece record,
  with sum_paper, sum_memo, sum_dp and sum_mm blocked and the transforms called) and bc_e2e (an end-to-end count with
  sum_cross at exactly the nodes with four groups or more, 6 of them; a second case with 4 such nodes took over six
  minutes in the agent-side code path and was left out). Validation: the head solution verifies all
  ten flags (verify 179 s) and scores 646 / 2,531 / 497 / 186 / 112 (2.1 s per choose_params call); the current
  code keeps its eight flags and 663 / 2,600 / 508 / 190 / 114; wrong variants fail bincross: flatpiece (all piece
  values in piece 0: records), pairswap (row pairing for the low pieces: records), swapsign (one sign in -S^-1 J:
  bc_fast) and lazy_e2e (count_e2e ignores the crossing: sum_cross calls 0 of 6).
- Examined and set aside: Astra's Schur crossing pieces (head analysis on 1-2 November): integral and
  factorial-free, but theta's Schur degree is tied to its edge shells; in the coupled cost model the rewrite costs
  more at 2 x 2 (r = 4: 0.935 against 0.582) and 2 x 3 (r = 6: 0.638 against 0.572) and less only at 2 x 4
  (r = 8: 0.541 against 0.561), where it still needs explicit piece tensors and would not beat one binary crossing at
  r = 7 in the search.
- Incident: at 10:46 the PC rebooted after a bugcheck 0x10E (VIDEO_MEMORY_MANAGEMENT_INTERNAL) while the head ran
  about 16 search and referee processes at once; two of them reported MemoryError just before. The validation was
  rerun with six processes (peak use well below 1 GB each); round 10 keeps background work to a few processes.
- T23 (iter_appd 009, attempt 1, 11:37 to 11:41, 14 requests): accepted at 11:48, bc_fast verified and the eight
  earlier flags kept (referee 382 s); third threshold 663 -> 663 (code task).
- T23 audit (git diff -w 897fefe 51637b9 -- algo_appd): exactly the spec. appd_core.py only (137 lines): the six
  4 x 4 matrices of NOTES 'Round 10' (S, K T and N0, N1 for both pairings), the helpers _bc_pair, _bc_mats, _bc_rows,
  _bc_rows_T, _bc_apply, _bc_add_into, _bc_copy and bc_mask, bc_left, bc_right, digit by digit, with None as an
  internal zero and the inputs never modified. No other change.
- T24 attempt 1 (iter_appd 010, 11:48 to 12:20, 61 requests): rejected at 12:27, bincross not verified: 3 of the 15
  hidden nodes (all with 5 groups) wrong in value and records. Head diagnosis (diag_t24.py): sum_cross built the
  group-table keys by joining two tuples (the spectator part and the shell part) instead of self.key of the merged
  dict; self.key sorts by mode, so whenever a spectator mode at a group had a larger index than a crossing mode the
  lookup missed and the term was skipped silently. A copy with sorted keys passes the three nodes exactly. The public
  nodes of check_t24.py did not expose it; the agent's own summary blamed the pass-through bits. A spec note for a
  third attempt was added at 12:45 and was not needed.
- T24 attempt 2 (iter_appd 011, 12:27 to 13:22, 35 requests; killed at the 55-minute wall while it ran extra fuzz
  tests of its own, code final since 12:39): accepted at 13:29, all ten flags verified, bincross and bc_fast included
  (referee 412 s); 663 -> 663 (code task). Head pre-check at 13:08 on a snapshot: the three nodes of attempt 1 exact.
- T24 audit (git diff -w 51637b9 7bb1676 -- algo_appd): spec items 1-5. appd_core.py only: node_value(method =
  'cross'); count_e2e(method = 'cross') with sum_cross at nodes with four groups or more and sum_dp otherwise;
  _bc_spread and sum_cross (blocks, words, pass-through bits and records as specified, sorted keys restricted to
  inc[g], block vectors through bc_left / bc_right, pieces through bc_mask, no call to sum_paper, sum_memo, sum_dp or
  sum_mm). The commit also rewrote every line of appd_core.py with CRLF endings (cosmetic; the file is uniform
  again, like appd_params.py, so the head left it).
- T25 (iter_appd 012, attempt 1, 13:29 to 13:40, 29 requests): accepted at 13:46, all ten flags verified (referee
  393 s; score gate 75 choose_params calls, 2.3 s each); third threshold 663 -> 646 (2.56% lower). The loop then
  found the backlog done, stopped the Strata server and exited.
- T25 audit (git diff -w 7bb1676 7a79a65 -- algo_appd): spec items 1-3 in appd_params.py; besides them a cache for
  A.eps_mode keyed by (r, shellmm, bincross), a pure function of those arguments there (modedp is always on), so
  it changes the speed only. No other change.
- Round 10 summary: appd-bincross adopted. 663 -> 646 against 2^(n/3), 2,600 -> 2,531 against 1.1199^n, 508 -> 497
  against 1.307^n, 190 -> 186 against 1.637^n and 114 -> 112 against 2^n, exactly the head's what-if. Attempts:
  T23 1, T24 2, T25 1; loop 11:36 to 13:46. Best parameters against 2^(n/3): three levels, r = 7, 7, 10, u = 105,
  86, 59, a = 128, 216, 276, delta = 0.002, Lc = 256, kappa = 10.
- Head re-run below log2 n = 1,000 of the main loop's counting cells (Sections 2-4; 2 November, asked by the owner
  of the page). The main referee's score starts its search at log2 n = 1,000 and reported <= 1,000 for the
  counting version against 2^n from T12 and against 1.637^n from T14. Same frozen assembly and Qwen's choose_params
  of commits 67882a5 (T12, referee v8) and f4948c7 (T14b, referee v9), one grid point every 5% from 10 to 1,000 and
  bisection (work/head/below1000/scan.py); a cell is the size from which every grid point up to 1,000 is below the
  target. T12 beats 2^n from 940 (19 groups, Q = 103, R = 35; e(939) = 1.0002, e(940) = 0.9936). T14b beats 2^n
  from 57 and 1.637^n from 61, but there the search uses 2 and 12 groups with Q = 1 and R <= 2, i.e. without the
  interpolation, and from about 90 to 1,000 the counting exponent stays at 0.500 whatever Q, so these two cells
  measure the group-by-group contraction with its edge bits on its own (a dynamic program like the 2^(n/3)
  path-decomposition count). Below those sizes the search is erratic (e < 1 again at 13-34, invalid at 36-41,
  its fallback at 43-55). The cost model's regime conditions are not re-checked at these sizes.

## Verification of the cost model after round 10 (2 November, asked by the owner of the page; head only)
- Scan (work/head/verify_r10/scan_appd.py; frozen assembly v6, the T25 search, 286 evaluations from log2 n = 64 to
  3,200 on a 1.5% grid, then 4,000 to 10^9, integer bisection): e(lam) decreases monotonically, no grid point is
  invalid, and the five cells come out exactly as scored, 646 / 2,531 / 497 / 186 / 112, each a size from which the
  route stays ahead (no grid point below the target before it). Margins at the cells: 1/3 - e = 1.3e-4 (646),
  3.4e-5 (2,531), 3.9e-4 (497), 3.8e-3 (186), 2.1e-3 (112).
- Condition margins at the cells (bits): at 646 (r = 7, 7, 10; u = 105, 86, 59; a = 128, 216, 276; delta = 0.002;
  Lc = 256; kappa = 10) the root label balance has 0.82, the root padding 4.27; at 497 the root balance 0.84; at 186
  and 112 (one level) the balance 0.63 and 0.72; at 2,531 (Lc = 2^20) the root padding 0.94. The precision headroom
  is tight at every level by construction (the search takes the smallest a), and it does not depend on lam.
- The counting formulas are upper bounds at every size, not asymptotic estimates. (a) Tuples: if the sorted mode
  weights satisfy w_i >= i/K, then for every t > 0 #(tuples of weight <= N) <= e^(tN) prod_i (1 - e^(-t w_i))^(-1)
  <= exp(tN + K pi^2/(6t)) (sum of f(i/K) <= K int f for decreasing f), so ln # <= 2 pi sqrt(KN/6); tup() is exactly
  this with K = cK (1 + dm) f M/D and N = 6 h D (the mode-weight premise is the round 7 LDU lemma). (b) shellmm:
  X Y Z min(X, Y, Z)^(omega - 3) is nondecreasing in each argument, so the tuple bounds and Cauchy-Schwarz give
  (2 pi/sqrt 6) sqrt((K_X + K_Y + K_Z - eta K_min) N) for every shell triple, i.e. f_mm. (c) bincross: a and b are at
  most the shell bounds plus one bit, and Cauchy-Schwarz on sqrt(wx) + (c - 1) sqrt(wy) gives u_theta.
  Exact integer check (exact_counts.py, cloud; K = 1, 2, 3, 5, 8 modes per weight, budgets up to 2,000): ln count
  over the bound at most 0.948; the best shell triple with Strassen at 0.87 to 0.90 of the f_mm bound in bits (six
  unit mixes, every case below, by 34 bits or more); the crossing index below the u_theta bound everywhere (by
  1.79 bits or more). Dropped factors are polynomial (block padding, Strassen constants, the number of shells,
  steps and nodes); at n = 2^lam with lam >= 64 and these margins they cannot move a cell.
- What a computation cannot settle: the correctness of each step (head proof sketches; the referee checks their
  identities on toy nodes, and its shellmm and bincross gates check the block shapes and per-piece channels the model
  counts), and the label-balancing constant. The binding size condition at 646 is the root label balance: about
  2^56.5 wires per label type against the 2^55.6 that the multi-stage Bernstein sketch (8192 Lc^2 / delta^2 times
  ln(172 Lc D)) asks for. A single-stage Bernstein bound for one type (summands <= 63 Lc) would ask for about
  14 bits less (41.6); the constant comes from splitting the deviation budget over about 42 Lc stages because labels are pairs
  of unit schedules. A sharper balancing proof is the only size-related lever left; round 9's what-if (20 bits less
  slack: 661 -> 654, then the padding binds) puts it at about 1%. Not pursued.

## Round 11 (2 November): a contraction order per shell profile (GPT-6 Astra) [T26, T27]
- appd-profile [head-verified sketch; proposed by GPT-6 Astra, "Two contraction orders, chosen separately for each
  shell profile", 2 November, after the head's verification of the cost model; extended by the head to menus of
  three and four orders].
  1. Complete profiles: every group pair e = (i, j) of a node gets one exact weight w_e (one shell of its modes),
     sum_e w_e < p; every mode tuple of weight < p lies in exactly one profile, so the profile values add up to the
     mode sum (D.8) in any contraction order. There are at most (p + 1)^C(r, 2) profiles, polynomially many.
  2. Inside a profile every edge has a finite alphabet of N_e tuples, l_e = log2 N_e <= beta sqrt(2 w_e) (the
     tuple bound of the frozen model, an exact Chernoff bound with the round-7 LDU mode lemma). A sequential order
     costs, at step k (A the first k groups, B the next, C the rest), one shell matrix product with blocked
     Strassen: 2^(X + Y + Z - (3 - omega) min(X, Y, Z)) = 2^max(alpha X + Y + Z, X + alpha Y + Z, X + Y + alpha Z),
     alpha = omega - 2, X, Y, Z the sums of l_e over A-B, A-C, B-C; this also pays for reading and writing the
     tables. So an order has the cost envelope C_pi(l) = max over its 3 (r - 1) forms a . l (coefficient 0 on the
     edges inside A, B or C, alpha on one cross class, 1 on the other two).
  3. Per profile the algorithm runs the cheapest order of a fixed menu (estimates with whole bits b_e = (N_e -
     1).bit_length() and alpha' = 101/125 > alpha, since 7^125 < 2^351; ties to the smaller index). For any
     weights lam_j >= 0 with sum 1: min_j C_j <= sum_j lam_j C_j = max over one form per order of
     (sum_j lam_j a_j) . l <= beta sqrt(2P) ||sum_j lam_j a_j|| (Cauchy-Schwarz with sum_e w_e <= P). With uniform
     lam the root fraction is eps <= (2 / r^2) max over the form combinations of ||mean of the forms||^2, computed
     exactly with alpha' (integer arithmetic); whole bits add at most C(r, 2) bits per profile. Child tables are
     memoised across profiles and stay in the child term; the polynomial factors (profiles, steps, block padding,
     Strassen constants) cannot move a cell.
  4. One order reproduces eps_mm(r) exactly (normalisation check). Menus (frozen in appd_assembly v7, PF_MENUS,
     first order the identity): three orders for r = 4 and 11..16, four for r = 5..10, found by local search in the
     cloud (search_orders.py); certificates PF_EPS (exact fractions): r = 4 0.4952, 5 0.4742, 6 0.4436, 7 0.4256,
     8 0.4087, 9 0.4008, 10 0.3905, 11 0.4021, 12 0.3933, 13 0.3915, 14 0.3878, 15 0.3838, 16 0.3800 (one
     crossing: 0.5545 to 0.5327; at r = 7 0.5544 -> 0.4256). A uniform profile costs the same in every sequential
     order, so no menu of sequential orders goes below 2 (12 + 3 alpha)^2 / (21 * 49) = 0.404 at r = 7.
  Head checks. Cloud: cert_head.py (independent re-derivation; Astra's four two-order certificates reproduced
  exactly, 129201/250000, 196478/390625, 274603/562500, 51847/109375; one order = eps_mm(r)), finalize_menus.py
  (exact integer certificates of the frozen menus). PC: pf_certificate in the v7 assembly (pure-Python re-derivation)
  equals PF_EPS for r = 4..12; test_pf_head.py on 66 random nodes and find_pf_cases.py on 692 more (4 to 6
  groups, power and LDU bases, chunks): the total equals (D.8) and sum_dp, every profile value equals a direct sum
  over the profile's tuples and is the same in every order of the menu, and the chosen order follows the rule.
  Caveat: head-level sketch, not peer reviewed; the bound lives in the frozen unit model.
- referee_appd version 7: code flag profile (needs shellmm; mode fraction min(eps_mm, eps_bc1, PF_EPS[r]) for
  4 <= r <= 16). Gates: pf_choose (the module-level rule on 160 random edge-bit dicts, r = 4..16), pf_node (14 hidden
  nodes with 4 to 6 groups, 26 to 1,287 profiles, non-trivial choices in every one: the value and every per-profile
  record with the other sums blocked, every step one matmul call of the recorded shape, and again with a chooser
  that cycles through the menu, so the contraction must follow the chosen order) and pf_e2e (two end-to-end
  counts with sum_prof at exactly the nodes with 4 to 16 groups, 6 and 4 of them). Validation: Validation of referee_appd version 7 (2 November, 17:19-17:28): the head solution verifies all 11 flags and scores 533 / 2,034 / 413 / 159 / 103 (third / occ3sat / ppsz / sharp3sat / brute); the untouched round-10 code keeps its 10 flags with profile false and its thresholds 646 / 2,531 / 497 / 186 / 112 (baseline unchanged); four wrong variants each verify the old 10 flags and lose the profile flag at the intended gate: always the first order (pf_node: records, matrix shapes and cycled run fail, values right), ties to the higher index (pf_choose), contraction without the matrix product (pf_node: matmul and cycled fail, values and records right), end-to-end count without sum_prof (pf_e2e: right count, sum_prof called 0 of 6 times).
- Audit T26 (iter_appd 013, commit 4b7e03c, accepted 2 November 18:11 on the first attempt; agent 1,425 s, 35
  requests). appd_core.py gains 182 lines: the import of appd_assembly_public, pf_cost, pf_choose, sum_prof (profile
  enumeration depth first over the edges, contraction of each profile in the chosen order with one matmul call per
  step) and method 'prof' in node_value and count_e2e. The file's line endings changed from CRLF to LF (git diff
  --ignore-cr-at-eol shows only the 182 added lines). Read in full: no call of sum_paper, sum_memo, sum_dp, sum_mm or
  sum_cross inside sum_prof, no file or network access, no stored values; pf_choose and matmul are looked up by their
  module-level names at call time, so the referee's replacements took effect (pf_node with the cycling chooser
  passes on all 14 hidden nodes).
- Audit T27 (iter_appd 014, commit 26671db, accepted 18:19 on the first attempt; agent 90 s, 10 requests).
  appd_params.py threads the flag profile through eps_cached, rank_e and the all-flags branch of the search (9
  changed lines). Thresholds 646 / 2,531 / 497 / 186 / 112 -> 533 / 2,034 / 413 / 159 / 103 (third / occ3sat / ppsz /
  sharp3sat / brute), equal to the head solution's validation run with the same parameters in every cell; all 11
  flags verified (verify 273 s, score 153 s).
- What sets the round-11 cells (head, from the frozen model's size conditions with the certificates' parameters):
  a label balance condition in every cell. At the top level against 2^(n/3) (r = 8, u = 79: 8^158 = 2^474 label
  types, 58.17 bits of balance, 0.83 bits to spare), 1.307^n (2^360, 52.31, 0.69), 1.637^n (2^114, 44.75, 0.27) and
  2^n (7^22 = 2^61.76, 40.82, 0.44); at depth 2 against 1.1199^n (0.52 bits to spare). The exponent e does not depend
  on lambda in the model; lambda enters only through the size conditions, so a cell is the smallest size at which
  parameters with e below the competitor satisfy them.

## External review notes (GPT-6 Astra, 2026-11-02 18:09, relayed by the owner) and the head's answers
Astra commented on the head's plain-language message about practical meaning (in short: "this route cannot reach a
practical size; the threshold is the number of label types times the balance per type; the cost per level falls only
like 1/sqrt(u)").
1. Certainty. The count r^(2u) is the cost of the present labeling design, and the balance condition is a sufficient
   condition of the construction, not a proven lower bound for every construction; "at most 11 percent from a
   perfect balancing argument" holds for the present parameters and that margin only. Accepted: the head's claim
   holds for the present design, not as an impossibility result; the chance that another design becomes practical is
   low but not known to be zero.
2. "Mostly safety margins were removed" was wrong. The first drop (1.11e11 to 496,177 against 2^(n/3)) was
   accounting; from 496,177 to 533 (about 930 times) the steps were algorithmic: each distinct sub-contraction once,
   group-by-group contraction, chunks, the Newton basis, the Cauchy LDU, shell matrix products, the guarded
   deformation, the binary crossing, the profile orders. Accepted. Round 11 alone shrinks the variable count at the
   break-even by 2^113, about 10^34.
3. "A practical break-even means a few hundred variables" compared worst-case upper bounds; their crossing is not a
   running time of a program on a given formula, and dividing 2^67 operations by a FLOP rate ignores operation types,
   memory and data movement; practical solvers are judged on real instances. Accepted; the page now says that the
   cells compare worst-case bounds and gives the few hundred variables only as the scale at which such bounds can be
   run at all.
4. 533 in the head's message was the reference solution's score. The loop's accepted result is the same.
5. Proposals. (a) Write-up with separate statuses (checked on small instances by the referee; exact-arithmetic
   certificates; proof sketches; the paper's theorem assumed). Accepted. (b) A separate, time-boxed practicality
   experiment: do profile orders and fast transforms cut time or memory of concrete contractions at real sizes, and can
   a useful piece join an existing exact counter (for example GANAK in its deterministic mode; decision and counting
   judged separately)? Open, the owner's decision. Head's prior is low: the profile choice acts on the shells of the
   deformation's mode sums, which practical counters do not form, and even the route's smallest settings use
   packets of length L = 2^29. (c) If the theoretical threshold is pursued further, the target should be the labeling
   structure ("must the u coordinate carry every label pair?"; the obstacle is keeping the packets' orientation and
   locality while dropping labels). Accepted as the largest remaining lever. Head estimate with the round-11 exponent
   parameters fixed and only the label term of the size conditions scaled: labels halved / quartered / none give
   376 / 303 / 231 against 2^(n/3) (then the deepest packet's padding binds), 289 / 236 / 182 against 1.307^n,
   1,412 / 1,106 / 799 against 1.1199^n, 101 / 72 / 45 against 1.637^n and 71 / 55 / 35 against 2^n (with no labels
   the balance condition itself goes away; every no-label value is then set by padding). No such construction exists
   yet; galactic in every case.
- General clause density with round 11's steps (head what-if, work/head/r11/whatif_dens_v7.py: the round-9 density
  patch of the frozen v7 assembly with the accepted T27 search; occurrence cap 1,024 by chains of copies, as checked
  by Fable for round 6): m = 4.26n gives 6,063 / 20,415 / 3,089 against 1.637^n / 1.307^n / 2^n (round 8's steps:
  7,689 / 25,946 / 3,911); m = 5n gives 8,660 / 29,219 / 4,401. Anthropic's own paper (26 October 2026, tensor-network
  contraction, Lean files in anthropics/formal-math) is reported (lilting.ch) to compute its constants as n >
  10^10000 against 2^n at m = 5n; the Appendix D route here needs about 10^1325 there. Different algorithms and
  analyses: a comparison of two published routes, not of one algorithm.

## External review notes (GPT-6 Astra, hostile review, 2026-11-02 20:20, relayed by the owner) and the head's checks
Astra received the review package (brief, claims.md, frozen v7 assembly, round-11 certificates, results page PDF) and was asked
to refute. Her report (appendix_d_hostile_review.md, with a checker and certificates
in a zip) and the head's answers:
1. The cells are not minimal (changes cells). The unchanged frozen v7 model accepts 100 / 157 / 532 / 411 / 2,033
   (brute / sharp3sat / third / ppsz / occ3sat) with finer delta and Lc and, for third and ppsz, r = 11 and smaller
   u, a at the last level. Head check (work/head/r12/check_astra_witnesses.py): all five rebuilt from her tables on the
   round-11 certificates are valid in appd_assembly v7 with exactly her exponents (to 1e-12) and invalid one size
   lower. ACCEPTED. The cells are feasible thresholds found by a parameter search, valid from there on (e does not
   depend on lambda, every size condition relaxes as lambda grows), not proved minimal; the page's wording "above
   which" is to be replaced.
2. Steps: no counterexample to the algebra of appd-profile, appd-bincross, appd-shellmm, appd-guarded; all 13 PF_EPS
   constants reproduced exactly. Weakest interfaces (not false identities): (a) the chunk-root extension of
   appd-bernstein (proof incomplete; she supplies a replacement proof that supports the existing root condition),
   (b) full work accounting around appd-memo and recursive packet-entry evaluation, (c) the simultaneous weighted
   vertex balance needed by appd-leafvars (repairable). Head check: with bincross off the five round-11 exponents are
   unchanged (PF_EPS is below eps_bc1 at every level used), confirmed; with ldu on, newton does not set cK. ACCEPTED as
   statements of dependence; the write-up should say so.
3. Matching-based root balance (new, Astra): edge-colour the crossing wires into q <= 2 Delta - 1 <= 126 Lc - 1
   matchings, balance every matching and label pair separately (Bernstein / pessimistic estimator, derandomisable):
   sufficient condition M/D >= (8q/delta^2) log(8qD), linear in Lc where the ledger's is quadratic; the old condition
   dominates it, so the existing root test stands. With only that line replaced she reports 92 / 148 / 526 / 403 /
   2,033. PENDING: the head checks the proof and the witnesses before any adoption (new model version).
4. Work accounting: complete profiles are polynomial for fixed r <= 16; matrix padding is constant-factor; memo
   lookups need an indexed memory model; implicit tensor-entry evaluation (D.4.6: exp(O(Q log s)); chunk tensors
   2^O(Lc)) is not in the assembly. Her scale check: these are far below the margin n (target - e) at every cell.
   ACCEPTED: the write-up states the bound as polynomial x entry overhead x 2^(n max(worst, cum + leaf)) and gives the
   margins per cell explicitly instead of "cannot move a cell"; small child networks need an explicit direct rule.
   The root balance fixes the size; it is not the largest arithmetic term in every cell (third: the middle term).
5. Competitors and interfaces: 1.637^n (Kutzkov, IJCAI 2026), 1.307031578^n (Jiang and Cai) and 1.1199^n (Jain, Neoh
   and Stephan) confirmed; the counting competitor should be written 2^(n/3 + o(n)) (using it without the o(n) is
   conservative for us); "3-CNF" must mean width at most three (with exactly three distinct variables and at most
   three occurrences every formula is satisfiable, by Hall's theorem, so the decision column would be trivial).
   Definition 3.1 and Lemmas 4.5 / 4.7 are not the Appendix D interface; the right obligations are packet placement
   and real-port balance, F(0), integral jets, strict degree coverage, recovery precision, integral local modes and
   immutable tensor descriptions. ACCEPTED (wording).
6. General density (breaks the stated inequality): 2m - n + n/Lc omits small connected components (one clause on
   three variables, Lc = 4: M_eff = 0 > -1/4). ACCEPTED. Head repair to check: small components solved directly;
   with no isolated variables, every small component has at least (n_c - 1)/2 clauses, so the large components carry
   at most 2m - 2n/3 + n/Lc wires (single three-variable clauses are the worst case); the density patch's
   M_root = 2 rho - 1 + ... becomes 2 rho - 2/3 + ...; rerun and publish the density certificates with the patch.
   The B^4 scaling of the chunk-balance constant is a scaling argument, not an exact certificate (to be replaced, for
   example by the matching-based condition with the cap's Delta). Anthropic's 10^10000 (checked by Astra in its PDF)
   is a necessary threshold for its displayed bound to beat exhaustive search, not a certified crossover; a later
   preprint (Majdoddin, arXiv 2610.41622) lowers its crossover points. The comparison must state these statuses.
7. Recommended wording, ACCEPTED: "a feasible threshold for this conditional leading-exponent model, found by the
   stated parameter search; not proved minimal; polynomial and tensor-entry costs need explicit accounting for a
   full finite-runtime comparison."
- Page update after the review (version 23): the five frozen-model witnesses in the tables (100 / 157 / 532 / 411 /
  2,033), "a feasible threshold found by a parameter search, not proved minimal", width at most three, the counting
  competitor 2^(n/3 + o(n)), the general-density entries withdrawn, and the dropped costs with explicit margins: at
  the witnesses n (target - e) exceeds the review's scale 2^a (lambda + a) of the tensor-entry costs by 2^52.9 (2^n),
  2^101.1 (1.637^n), 2^300.4 (2^(n/3)), 2^225.7 (1.307^n) and 2^1220.6 (1.1199^n).

## Round 12 (2 November, evening): the hostile review's corrections and the matching-based root balance [head]
No loop task: the round changes one size condition of the model and the parameter search, no algorithm code.
- appd-matchbal [head-checked sketch; proposed by GPT-6 Astra in her hostile review, 2 November].
  1. A chunk unit has at most 1 + 7(Lc - 1) variables (appd-chunks), at most 21 Lc designated clauses with at most
     two outside incidences each, and its variables have at most 21 Lc outside incidences: wire degree
     Delta <= 63 Lc (the same count behind the old constant 8192 ~ 2 * 64^2).
  2. The crossing wires form a loopless multigraph on the units (parallel and opposite wires allowed); greedy edge
     colouring gives q <= 2 Delta - 1 <= 126 Lc - 1 matchings.
  3. Unit labels independent and uniform in [r]^u (K = r^u, D = K^2 ordered pairs): inside one matching the wires
     have disjoint ends, so the load X_j of a fixed pair is Binomial(|E_j|, 1/D), mean mu_j; sum_j mu_j = mu = M/D.
  4. Budgets t_j = (delta/2)(mu_j + mu/q) add up to delta mu. Bernstein: Pr[|X_j - mu_j| > t_j] <=
     2 exp(-t_j^2 / (2 (mu_j + t_j/3))); the exponent is smallest at mu_j + mu/q = 2 (mu/q) / (1 + delta/6), where
     it equals delta^2 mu / (2 q (1 + delta/6)^2), so every tail is <= 2 exp(-delta^2 mu / (5 q)) for delta <= 1/4
     (Astra's weaker form, with a rational pessimistic estimator that derandomises the labels by conditional
     expectations, polynomial in the network and in D <= M).
  5. Union bound over the 2 q D events: with M/D >= (8 q / delta^2) ln(8 q D) the total failure probability is at
     most (1/4) (8 q D)^(-0.6) < 1/4. Root slack: log2(8 q) + 2 log2(1/delta) + log2(ln(8 q) + ln D).
  6. Port bounds follow by summing the pair loads; the weighted vertex events of appd-leafvars (weights at most
     28 Lc) have ample slack under the same condition (Astra; not separately quantified by the head).
  Head checks: Delta <= 63 Lc against the appd-chunks definition; the Bernstein minimum by calculus; Astra's
  witnesses 92 / 148 / 526 / 403 / 2,033 reproduced exactly by an independently written patched copy of the
  frozen model (work/head/r12 and the cloud scratchpad) and by referee v8. The old two-stage condition needs 9 to
  13 bits more at every (r, u, Lc) used, so all earlier cells stay valid. Caveat: head-level sketch, not peer
  reviewed.
- referee_appd version 8: the assembly gains the analysis flag matchbal (needs chunks; not a code flag); gates and
  cases as in version 7. Without matchbal every round-11 certificate scores exactly as before.
  Validation (2 November, 20:48-20:56): the head solution (work/head/r12/hs12) verifies all 11 code flags and scores 525 / 2,033 / 404 / 147 / 93 (third / occ3sat / ppsz / sharp3sat / brute) at 2.86 s per call, matchbal on in every certificate; the untouched round-11 code keeps 533 / 2,034 / 413 / 159 / 103; a variant that sets matchbal without chunks keeps its code flags and scores none / none / none / none / none (every dictionary rejected).
- Round-12 cells (frozen v8 model, matchbal on): the smallest valid dictionary found by three searches, the
  head's online search (round-11 search plus matchbal in its replica and a coordinate refinement over delta, Lc,
  every r_d and kappa; 2.75 s per call), a head differential evolution over r_d, u_d, delta, Lc and kappa (two
  passes, seeded with the review's witnesses and the online search's best), and the review's own search. Every
  dictionary was checked on the PC with referee_v8/appd_assembly.py: valid at its size with the stated exponent,
  invalid one size lower. brute 92 (r 7, u 11, a 29, delta 0.0521, Lc 64, kappa 10; review and DE), sharp3sat 147
  (r 8, u 18, a 39, delta 0.00702, Lc 256, kappa 12; DE and online), ppsz 403 (r 8/9/10, u 60/46/36, a
  82/130/168, delta 0.0045, Lc 512, kappa 10; review and DE), third 524 (r 8/9/11, u 79/61/43, a 103/166/211,
  delta 0.00121, Lc 4096, kappa 12; DE pass 2), occ3sat 2,032 (r 8/9/16, u 323/256/172, a 350/608/782, delta
  0.000374, Lc 65536, kappa 11; DE, valid even without matchbal, so also a v7 witness). Binding condition: the
  root balance against 2^(n/3) (2^474 types and 49.78 bits, 0.22 to spare), 2^n (2^61.8, 30.25) and 1.637^n; the
  balance one level down against 1.307^n; the root padding against 1.1199^n. Dropped-cost margins: n (target - e)
  exceeds 2^a (lambda + a) by at least 2^46.3 (2^n), 2^89.5, 2^288.4, 2^211.7 and 2^1220.1. Label what-if with
  these exponent parameters fixed: half the labels 355 (2^(n/3)) and 61 (2^n); no labels 225 and 34 (padding).
  Status: feasible thresholds of the conditional model found by these searches, not proved minimal.
- appd-density repaired [head, after the review]: small connected components are solved directly; with no isolated
  variables a connected component with c clauses of width <= 3 has at most 2c + 1 variables, so the small ones
  remove at most n_s/3 from 2m - n and the large components carry at most 2m - 2n/3 + n/Lc wires (a single
  three-variable clause is the worst case); M_root = 2 rho - 2/3 + 3 rho/(B - 2) + 1/Lc with the chains of copies
  as before (equalities are arity-two constraints, so the copies add no wires). Root balance by matchings with the
  exact unit degree for occurrence cap B: units of at most (2B + 1) Lc variables (variable-graph degree <= 2B), at
  most B (2B + 1) Lc designated clauses with two outside incidences each and B (2B + 1) Lc outside incidences of
  their variables, Delta <= 3 B (2B + 1) Lc, q = 2 Delta - 1 (B = 3 gives 126 Lc - 1); this replaces the old
  (B/3)^4 scaling of the two-stage constant. work/head/r12/whatif_dens_v8.py (frozen v8 assembly and the head's
  round-12 search, both patched only in M_root, q and the leaf term): m = 4.26n gives 6,594 / 22,223 / 3,357 against 1.637^n / 1.307^n / 2^n (the withdrawn round-11 entries were 6,063 / 20,415 / 3,089); m = 5n gives 9,293 / 31,365 / 4,722 (withdrawn: 8,660 / 29,219 / 4,401), n about 10^1,421 against 2^n.
  Caveat: the weighted vertex events of the leaf step at high density are not separately quantified.
- Page version 24 (2 November, about 21:05): round-12 rows in both Appendix D tables (147 / 524 / 403 / 2,032 and 92 / 147 / 524 / 403 / 2,032), summary row, comparison and practical paragraphs, the round-12 paragraph, the density entries restored with the repaired count, dropped-cost margin at least 2^46.
