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
  1. Multiplicities q_t = h max(0, u - i) with h = ceil(M(2L - 1) / (L(2u - 2 + 2^(1-u)))) ~ M/(u - 1)
     already give sum_t q_t > M(2L - 1) (the paper takes h = ceil(3M/u) and counts odd points only).
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
Where the exponent sits: at every used node level K p_t = 2 (1 + delta) f M^2 / (u - 1) (independent of i), so a
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
