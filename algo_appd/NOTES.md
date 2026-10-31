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
