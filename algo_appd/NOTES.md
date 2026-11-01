# Appendix D route: background for the implementer

The paper's Appendix D counts the solutions of a CNF formula by contracting a Boolean tensor network. A call
deforms the wire identifications with a Cauchy matrix W(T) inside packets of L = 2^a slots, evaluates the
deformed contraction F(T) near interpolation nodes t, expands the interactions between groups of vertices in
"modes" Z (one per packet, tile and power b; each mode has one source group and one target group), and recovers
F(0) by tapered Hermite interpolation. The provided code (appd_core.py, everything up to NodeExpansion.__init__)
builds all of that; you implement the summation of the mode expansion:

    F(t + alpha z) = sum over tuples beta of weight < p of  prod over groups g of  V_g[beta restricted to I_g]

where I_g = self.inc[g] is the set of modes whose source or target is g, and V_g[gamma] is the contraction of the
child network of group g for gamma (self.child_network_value(g, gamma) computes one entry slowly).

- memo: every V_g[gamma] is computed once and kept in a table (group_table).
- dp: the sum over beta is eliminated group by group; the state is the tuple restricted to the cut (modes
  with one end processed and one not). Tuples of weight >= p may be skipped or included: they vanish mod 2^p.
- byvars (chunk networks): V_g[gamma] can also be computed by summing over the variables of the chunks of g,
  splitting each clause with external wires into 1 - prod [literal false] and pushing unary weights on the
  row bits into the packet factors.
- count_e2e: the whole algorithm with a = 1, u = 1 gives the exact count.

All series are lists of q integers mod 2^p (coefficients of z^0 .. z^(q-1)); use ps_mul_mod and ps_add_mod.
The referee checks your functions against frozen reference values on hidden formulas and blocks
direct_taylor_mod and F_poly while it does so. appd_params.py chooses the parameters of the cost model; the
referee re-scores them with its own frozen copy of appd_assembly_public.py.

## Round 6: binomial (Newton) moments
At a node t with i = v2(t), c = 2^(i+1) and alpha = 2^(a+1), a packet tile y0 has the columns y0 + c j
(0 <= j < L/c). For a crossing row x, W_(x, y0 + c j)(T) = G(T) / (A - c j) with A = T - y0 + x and
G = C(T)/C'(0). The power basis expands 1/(A - c j) = sum_b (1/A) (c/A)^b j^b: source forms sum_j j^b X_y.
Newton's forward differences give an exact, finite expansion instead:
    G / (A - c j) = sum_(b = 0..j) N_b * binom(j, b),    N_b = (G/A) * prod_(l = 1..b) c l / (A - c l),
so the source forms become sum_j binom(j, b) X_(y0 + c j), and N_b = N_(b-1) * c b / (A - c b).
Because v2(A - c l) = v2(t - y0 + x) <= i for every l, each factor c l / (A - c l) is divisible by 2^(1 + v2(l)),
so N_b carries v2(b!) more factors of 2 than the power coefficient: the weight of mode b grows by
v2(b!) = b - (number of ones in the binary form of b). Weights then grow about twice as fast, a class has about
half as many modes per unit of weight, and the frozen model counts mode tuples with the constant cK = 1 instead of
cK = 2 when the newton flag is on (appd_assembly_public.ck). The sums, the child tables and byvars_value need no
change: they only read self.modes (dict P, y0, b, src, tgt, w, lin) and self.coef.

## Round 6 (LDU): Cauchy LDU of the crossing block
At a node t (i = v2(t), c = 2^(i+1), alpha = 2^(a+1)) the crossing entries of a packet are the pairs of slots
(x, y) with the same bit i. Take one half h (slots with bit i = h) in the order x_k = c (k >> i) + (k mod 2^i)
+ h 2^i, k = 0 .. L/2 - 1, all L/2 positions whether occupied or not. On that block W_(x_k, x_l)(T) = G(T) K_kl
with K_kl = 1/(T + x_k - x_l), a Cauchy matrix. Factor K = L diag(d) U over power series in z (Doolittle without
pivoting: series_ldu). Then the crossing part of row k is G sum_(m <= k) L_km d_m (U X)_m, so mode m has the
source form (U X)_m = sum_(l >= m) U_ml X_(x_l) (its coefficients are z-series) and, in row k >= m, the
coefficient G L_km d_m. Closed forms (checked exactly): d_j = (1/T) prod_(m<j) -D^2/((T + D)(T - D)) with
D = x_j - x_m, and L, U have entries in Z_(2)[[z]]. Hence v2(d_j) >= -i + 2b + 2 v2(b!) with b = j >> i, and with
v2(G) >= a + 1 the weight of mode m is a + 1 - i + 2b + 2 v2(b!): twice the Newton growth. The frozen model uses
cK = 1/2 when the ldu flag is on. Only the sharp weights exist for this basis. Because the source forms have series
coefficients, every place that evaluates a form must multiply series (helpers lin_eval, src_mul, src_scale,
src_zero in appd_core.py keep int-only forms on the exact old code path).
