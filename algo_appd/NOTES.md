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


## Round 7: the mode sum as shell matrix products (appd-shellmm)
sum_dp eliminates the groups in the order r_groups. At the step for group g the state tuples live on the cut
modes. Split them into X (cut modes that do not touch g: processed group to later group) and Y (cut modes that
touch g: processed group to g), and let Z be the new modes (g to later groups). The step computes
Phi'[x + z] = sum over y of Phi[x + y] * T_g[y + z] over w(x) + w(y) + w(z) < p: a matrix product with rows x,
inner index y and columns z, restricted by the weight budget. Grouping each of X, Y and Z by exact weight
(shells) makes every allowed shell triple (w1, w2, w3), w1 + w2 + w3 < p, one full block product
A[X_w1 x Y_w2] * B[Y_w2 x Z_w3]; summed over the triples this is the same sum as sum_dp, term by term. A fast
matrix multiplication algorithm then multiplies the blocks. Strassen's algorithm uses only ring additions,
subtractions and multiplications, so it works over (Z/2^p)[z]/(z^q); on square blocks of side S (the smallest
dimension) it needs about n k l S^(log2(7) - 3) ring products instead of n k l. In the cost model the dominant
shell then costs 2^(log N_x + log N_z + (omega - 2) log N_y) (y the smallest index set) instead of
2^(log N_x + log N_y + log N_z); re-optimising the weight split gives the mode fraction
eps_mm(r) = max_k 2[(k-1)(r-k) + (k-1) + (r-k) - (1 - (omega - 2)^2) min(k-1, (k-1)(r-k), r-k)] / r^2 with
omega = log2(7) (appd_assembly_public.eps_mm, about 7 percent below eps_dp at r = 7). The frozen model uses it
when the shellmm flag is verified.


## Round 8: the guarded Lagrange deformation (appd-guarded)
Proposed by GPT-6 Astra, checked by the head. A packet is a virtual grid of Q = 2^a slots; real wires use only
the central slots H <= x < Q - H with H = Q/2^kappa (guards), and labels live on the evaluated levels 0..u-1
(u <= a - kappa). The deformation is W~_xy(T) = l_y(x + T), where l_y is the Lagrange basis polynomial of y on the
grid {0..Q-1}: W~(0) = I, every entry has degree Q - 1 (the paper's Cauchy deformation has 2Q - 1), and at a node
t in [-H, H) the matrix W~(t) is the shift x -> x + t (x + t stays on the grid thanks to the guards), so nodes with
v2(t) = i < u are local exactly as before. Near a node, T = t + alpha z with alpha = 2Q as before, and a crossing
entry factors as W~_xy = R_x(z) * 1/(t + alpha z + x - y) * c_y with the row factor
R_x = 2^-v P_Q(x + t + alpha z) (coefficientwise divisible by 2^(a+1), v = Q - 1 - a) and the odd unit
c_y = 2^v / P_Q'(y). So the Cauchy LDU of the full virtual half works unchanged: the mode coefficients become
R_k d_m L_km (instead of G d_m L_km) and the source forms get the units, (U diag(c) X)_m; the mode weights stay
a + 1 - i + 2b + 2 v2(b!). The node grid is [-H, H) (complete residue systems at every level), multiplicities
q_t = h (u - v2(t)) with h = floor(M (Q - 1) / sum_(i<u) (H >> i)(u - i)) + 1, and the Hermite loss is at most
6 h H / 2^i: about half of the paper's, because the degree halved. GuardedPacking (packing, Wfrac, rowfactor,
colunit), guarded_multiplicities and guarded_h are given in appd_core.py.
