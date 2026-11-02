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

## Round 9 (T22): a better parameter search (no new math, referee unchanged)
The frozen assembly's label-balance test (appd-bernstein, appd-chunks) often binds before the padding budget.
In exponent_of, section 'labels': lD = 2 u log2(r_d), logM = lam + log2(M_d), lnD = max(lD ln 2, 1); at the root
(chunks) slack = log2(8192 Lc^2) + 2 log2(1/delta) + log2(ln(172 Lc) + lnD); below the root
slack = log2(8) + a_(d-1) + 1 + 2 log2(1/delta) + log2(lnD); the test passes if logM - lD >= slack.
est_u (and so rank_e) looks only at the padding budget, so the pre-rank always prefers the chunk size Lc = 2^20
and never refines the smaller chunk sizes, where the root balance leaves room for a larger u_0. With est_u = the
smaller of the padding budget and the balance test, a finer Lc and delta grid and 20 refined candidates, the frozen
model gives about 663 against 2^(n/3) (head dry run; the current search needs 686).

## Round 10: one binary crossing with fast absorption (appd-bincross)
Proposed by GPT-6 Astra, checked by the head (exact checks on 600 random nodes). The mode sum of a node is a tensor
network on the groups: the modes whose two groups are {g, h} form one edge, and their exponent tuple is the index
that T_g and T_h share. Take g1, g2, g3, g4 = r_groups[0], r_groups[1], r_groups[2], r_groups[3] and the two edges
X = the modes of {g1, g3} and Y = the modes of {g2, g4}. Seen from the pairs {g1, g2} and {g3, g4} they cross.

Codes. Fix a shell of X (exact weight wx) and a shell of Y (wy), from self.shells(X) and self.shells(Y) (keep their
order). The code of a tuple is its position in its shell; a shell of N tuples has a = (N - 1).bit_length() bits
(a = 0 when N = 1), most significant bit first. Let b be the bits of the Y shell and m = min(a, b). The top m bits of
the two codes are paired; the low |a - b| bits of the longer code (call their value e) pass straight through.

Pieces. x and v are the X-tuples seen at g1 and g3, y and u the Y-tuples seen at g2 and g4; the two edges say
u = y and v = x. Per paired bit this is the 4 x 4 matrix SWAP = I - J with J = j j^T, j = (0, 1, -1, 0), on the digit
2 * (first bit) + (second bit). For m bits SWAP^(m) = sum over i = 0..m of A_i = [t^i] (I - tJ)^(m) (the pieces).
Each piece factors through few channels with integer matrices, in one of two pairings:
- row pairing: left block A = (g1 holding x, g2 holding y), right block B = (g4 holding u, g3 holding v);
- column pairing: left block A = (g1 holding x, g4 holding u), right block B = (g2 holding y, g3 holding v).
Piece i uses the column pairing when 2 i <= m, otherwise the row pairing. A block's word has m base-4 digits, digit k
(k = 0 the most significant) = 2 * (bit k of the first end's top bits) + (bit k of the second end's top bits).

The three module-level helpers (vectors are lists of 4^m series, None = zero, index = word):
- bc_mask(m, i, pairing): the words with at least i zero digits (row) or at least m - i zero digits (column), in
  increasing order.
- bc_left(vec, m, pairing, q, mod) = (M^T)^(m) vec with M = S (row) or K T (column): for k = 0..m-1 replace vec by
  new[w'] = sum over d of M^T[d'][d] * vec[w' with digit k set to d], where d' = digit k of w'.
- bc_right(vec, m, pairing, q, mod) = list R[0..m] with R[i] = [t^i] (N0 + t N1)^(m) vec: keep one vector per power of
  t; at digit k, N0 keeps the power and N1 raises it by one.
- piece i of a block pair = sum over w in bc_mask(m, i, pairing) of L[w] * R[i][w] (series products mod (2^p, z^q)).
Matrices (row index = output digit d', column index = input digit d):
  S    = [[0,1,0,0],[1,0,0,0],[-1,0,1,0],[0,0,0,1]]     K T = [[1,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]]
  row:    N0 = S^-1 = [[0,1,0,0],[1,0,0,0],[0,1,1,0],[0,0,0,1]]   N1 = -S^-1 J = [[0,-1,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
  column: N0 = T^-1 P0 = [[1,0,0,1],[0,0,0,0],[0,0,0,0],[0,0,0,0]] N1 = -T^-1 = [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[1,0,0,-1]]
(M^T)[d'][d] = M[d][d']. check_t23.py compares the helpers with dense Kronecker products of these matrices.

sum_cross (returns (series, records, child_calls)), pseudo-code:
    child_calls = 0; tables = group_table(g) for every group; g1..g4 = r_groups[:4]; singles = r_groups[4:]
    X, Y as above (sorted mode indices); others = every mode not in X or Y; total = zero series; records = []
    blocks: row -> ((g1, x, g2, y), (g4, u, g3, v)); col -> ((g1, x, g4, u), (g2, y, g3, v))
    for sx in shells(X), wx = its weight:  for sy in shells(Y), wy = its weight:
        skip when wx + wy >= p; budget = p - wx - wy; a, b, m as above; the shell of x and v is sx (a bits), of y, u sy
        for each pairing used by some i in 0..m, for each side (left, right) with block (hA, endA, hB, endB):
            bnd = sorted((inc[hA] | inc[hB]) - X - Y)            (the block's other modes, its inner edge included)
            for bt in tuples(bnd, budget), for code cA, tuple tA of endA's shell, for code cB, tuple tB of endB's shell:
                val = tables[hA][key((bt + tA) restricted to inc[hA])] * tables[hB][key((bt + tB) restricted to inc[hB])]
                word from the top m bits of cA (first) and cB (second); e = the low |a - b| bits of the longer code
                add val into vecs[(e, key(bt))][word]
            left side: L[(e, key)] = bc_left(vec, m, pairing, q, 2^p); right side: R[(e, key)] = bc_right(...)
        for i in 0..m: pairing = col if 2 i <= m else row; mask = bc_mask(m, i, pairing); piece = zero series
            for bt in tuples(others, budget):
                kA = key(bt restricted to the left block's bnd); kB = key(bt restricted to the right block's bnd)
                s = product over singles s of tables[s][key(bt restricted to inc[s])]
                for e in 0..2^|a-b| - 1 with (e, kA) in L and (e, kB) in R:
                    piece += s * (sum over w in mask of L[(e, kA)][w] * R[(e, kB)][i][w])
            total += piece; records.append(dict(wx, wy, nx = len(sx), ny = len(sy), a, b, m, i, pairing, channels =
            len(mask), extra = |a - b|, value = piece mod 2^p))
Every assignment of total weight < p appears once in every piece, so the pieces add up to sum_paper exactly; a
piece alone is the network with A_i in place of the two crossing edges. The referee checks the series and every
record, with the old sums (sum_paper, sum_memo, sum_dp, sum_mm) blocked.

Why it is cheaper (the cost model): a piece has at most 2^(c m) channels, c = log2(12)/2 = 1.79 per bit pair instead
of 2, so the new index theta of a piece has 2^(max(a,b) + (c - 1) min(a,b)) values: 3.256 mode units instead of the
4 of two edges. The transforms cost a factor m^2 per entry of the old pair table and use no division. With pair A,
the singletons one at a time and pair B last the mode fraction drops from 0.5696 to 0.5544 at r = 7
(appd_assembly_public.eps_bc1, flag bincross).

## Round 11: a contraction order per shell profile (appd-profile)
Proposed by GPT-6 Astra, checked and extended by the head (exact certificates; 66 random nodes where every profile
value equals a direct sum and is the same in every order). The mode sum of a node with r >= 4 groups is a tensor
network on the groups: the modes whose two groups are r_groups[i] and r_groups[j] (i < j) form the edge (i, j), and
their exponent tuple is the index that T_g and T_h share. Edge order: (0,1), (0,2), ..., (0,r-1), (1,2), ...,
(r-2,r-1). Every mode joins two different groups, so every mode lies on exactly one edge.

Profiles. Every edge has its shells: self.shells(sorted modes of the edge), in increasing weight (an edge without
modes has the single shell [{}] of weight 0). A profile picks one shell per edge with total weight < p. Enumerate
them depth first: edges in edge order, shells in increasing weight; at an edge stop as soon as a shell's weight
reaches the remaining budget (p minus the weights already chosen). Inside a profile every edge has a finite
alphabet, the tuples of its shell in shell order: N_e tuples and b_e = (N_e - 1).bit_length() bits.

Order choice, two module-level functions:
- pf_cost(order, bits, r): order = tuple of the positions 0..r-1; bits = dict {(i, j): b_e} over all edges (i < j).
  For k = 1..r-1: A = order[:k], B = order[k], C = order[k+1:]; X = sum of the bits of the edges between A and B,
  Y = between A and C, Z = between B and C; v_k = 125 * (X + Y + Z) - 24 * min(X, Y, Z). Return max over k of v_k
  (an integer: 125 times the Strassen cost exponent of the step with alpha = 101/125).
- pf_choose(bits, r, menu=None): menu = PF_MENUS[r] of appd_assembly_public when None; costs = [pf_cost(o, bits, r)
  for o in menu]; j = the index of the smallest cost, the smallest index among equal costs; return (j, costs).

sum_prof(self) (a NodeExpansion method; returns (series, records, child_calls)), pseudo-code:
    child_calls = 0; tables = group_table(g) for every group; r = len(r_groups); menu = PF_MENUS[r]
    edges and profiles as above; total = zero series; records = []
    for each profile, in enumeration order:
        alph[e] = the chosen shell of edge e (list of tuples); sizes[e] = len(alph[e])
        bits = {(i, j): (sizes[e] - 1).bit_length() for every edge e = (i, j)}
        j, costs = pf_choose(bits, r, menu)      # call it by this name, exactly once per profile
        value, shapes = the profile's network contracted in the order menu[j] (below)
        total += value
        records.append(dict(weights=[weight of the chosen shell, per edge], sizes=sizes, order=j, costs=costs,
                            shapes=shapes, value=value mod 2^p))
    return total, records, self.child_calls

Contraction of one profile in an order (positions):
    entry(g, one tuple per edge of g) = tables[g][self.key(union of those tuples)]
    state T: dict from (one alphabet index per open edge) to series; the open edges are the edges with exactly
    one end processed, kept in edge order.
    start: processed = {order[0]}; T[c] = entry(order[0], c) for every combination c of its edges' alphabets.
    for k = 1..r-1:
        b = order[k]; AB = the open edges touching b; AC = the other open edges (both in open-edge order);
        BC = the edges of b that are not open (to groups later in the order), in edge order
        rows = itertools.product over the AC alphabets (index tuples), inner = over AB, cols = over BC
        Am[row][inner] = T[the AC and AB indices placed in open-edge order]
        Bm[inner][col] = entry(b, the AB and BC indices)
        Cm = matmul(Am, Bm, q, 2^p)          # the module-level matmul, one call per step
        shapes.append((len(rows), len(inner), len(cols)))
        processed += {b}; the new open edges are AC and BC in edge order; T[their indices] = Cm[row][col]
    value = T[()] after the last step.

node_value(..., method="prof"): ser, records, calls = ne.sum_prof(); stats = {"records": records, "child_calls":
calls, "modes": len(ne.modes)}. count_e2e(..., method="prof"): at every node ne.sum_prof() when len(ne.r_groups) is
a key of PF_MENUS (4 to 16 groups), otherwise ne.sum_dp() (as for the default method).
Every mode tuple of weight < p lies in exactly one profile, so the profile values add up to sum_paper exactly, in
any order. The referee checks pf_choose on random bits, every per-profile record and the shape of every matmul
call, again after replacing pf_choose by a chooser that cycles through the menu, and end-to-end counts, with
sum_paper, sum_memo, sum_dp, sum_mm and sum_cross blocked.

Why it is cheaper (the cost model, flag profile): with one fixed order the worst profile sets the cost; choosing the
order per profile pays only min over the menu, and min <= the mean of the orders' cost envelopes. With
log2 N_e <= beta sqrt(2 w_e) and Cauchy-Schwarz this gives the mode fraction PF_EPS[r], an exact certificate over
all combinations of one cost form per order (appd_assembly_public.pf_certificate). At seven groups 0.554 (one
crossing) becomes 0.426 with the four orders of PF_MENUS[7]; whole bits add at most C(r, 2) bits per profile, and
profiles and steps are polynomially many.
