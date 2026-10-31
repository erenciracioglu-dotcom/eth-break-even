"""End-to-end run of the Section 2-3 counting algorithm on small formulas, mod a prime p.

H(T) = sum_x prod_f Phi_f(<P_{f,v}(T), x>),  H(0) = #F.
For t = 1..Q the Taylor data H(t+z) mod z^J are computed via Lemma 3.3 (inclusion-exclusion
over U <= E, |E| <= q); each restricted sum Z_{t,E,U} can be computed by Lemma 3.2 (fix
measurement tuples, sum every group separately using representatives), which is the
algorithm of the paper. Hermite interpolation of the Taylor data then gives H(0) mod p.

PUBLIC API (the referee calls these; keep names):
  Counter(aux, des, p=PRIMES[0], residue_keys=False, weights=None).taylor_inclusion_exclusion(t,
  Zfun=None), .Z_group_separated(t, S, U), .taylor_aggregated(t, Zagg=None), .Zagg_brute(t, S, U),
  .Zagg_group_separated(t, S, U), .restricted_sums, .hermite_H0(data), .fibre_counts(t), .Q;
  brute_count(n, clauses)

Modular form (cost model v5, yes/no targets): with residue_keys=True the fibres of every group
are keyed by the measurement values mod p (they only matter mod p), and with weights = w (int
array of length n, entries mod p) every assignment x carries the factor prod_v w_v^(x_v), so
Hermite interpolation of the aggregated Taylor data gives the weighted count
S(w) = sum over satisfying x of prod_v w_v^(x_v) mod p.

Aggregated form (cost model v3): Lemma 3.3 summed per U,
  sum_x p_x(x) = sum_{|U| <= q} sum_{x >= U} w_K(|x - U|) p_x(1_U),  K = q - |U|,
  w_K(m) = sum_{j <= K} (-1)^j C(m, j),
so one restricted sum per U instead of one per pair U <= E; the weight is additive in the
number of ones of x outside U, which a polynomial in y tracks inside the group separation.
"""
import itertools
from math import comb
import numpy as np

from construction import PRIMES


def smul(a, b, p):
    """Truncated product of power series stored as (J, N) arrays."""
    J = a.shape[0]
    c = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.int64)
    for i in range(J):
        ai = a[i]
        if not np.any(ai):
            continue
        c[i:] = (c[i:] + (ai[None, ...] * b[:J - i]) % p) % p
    return c


def mul_zy(A, B, p):
    """Product of two (J, Ky) tables, polynomials in z (rows, degree < J) and y (columns,
    degree < Ky), truncated in both; entries of A, B are kept below p so int64 is safe."""
    J, Ky = A.shape
    C = np.zeros((J, Ky), np.int64)
    for i in range(J):
        ai = A[i]
        if not np.any(ai):
            continue
        for j in range(J - i):
            bj = B[j]
            if not np.any(bj):
                continue
            C[i + j] = (C[i + j] + np.convolve(ai, bj, mode="full")[:Ky] % p) % p
    return C


def brute_count(n, clauses):
    xs = np.arange(1 << n)
    sat = np.ones(xs.shape, bool)
    for cl in clauses:
        c = np.zeros(xs.shape, bool)
        for v, pos in cl:
            bit = (xs >> v) & 1
            c |= (bit == 1) if pos else (bit == 0)
        sat &= c
    return int(sat.sum())


class Counter:
    def __init__(self, aux, des, p=PRIMES[0], residue_keys=False, weights=None):
        self.aux, self.des, self.p = aux, des, p
        self.residue_keys = bool(residue_keys)      # fibre keys = measurement values mod p
        self.n, self.Q = des["n"], aux.Q
        clauses = des["clauses"]
        self.clauses = clauses
        occ_index = {fv: i for i, fv in enumerate(des["occ"])}
        self.cl_occ = [[(occ_index[(f, v)], pos) for v, pos in cl] for f, cl in enumerate(clauses)]
        K = des["Prow"].shape[0]
        self.K = K
        invL = pow(aux.Lden % p, p - 2, p)
        self.P = (des["Prow"] % p) * invL % p                    # (K, M, n), unscaled mod p
        if des.get("designated"):
            # designated occurrences carry no replacement row: only the clause occurrences
            # that are NOT designated count (with units a clause may have none of them)
            d = sum(1 for i, (f, _v) in enumerate(des["occ"])
                    if f < len(clauses) and not des["desig"][i])
        else:
            d = sum(len(cl) for cl in clauses)
        self.dH = d * (aux.R - 1)                                # deg_T H <= d (R-1)
        self.J = self.dH // self.Q + 1
        self.q = self.J - 1
        N = 1 << self.n
        self.N = N
        xs = np.arange(N)
        self.X = ((xs[:, None] >> np.arange(self.n)[None, :]) & 1).astype(np.int64)
        self.pop = self.X.sum(axis=1)
        # self.wX[x] = prod over v with x_v = 1 of w_v (mod p); all ones when weights is None
        if weights is None:
            self.wX = np.ones(N, np.int64)
        else:
            w = np.asarray(weights, np.int64).reshape(-1) % p
            if w.size != self.n:
                raise ValueError("weights must have length n")
            self.wX = np.ones(N, np.int64)
            for v in range(self.n):
                if w[v] != 1:
                    sel = self.X[:, v] == 1
                    self.wX[sel] = (self.wX[sel] * int(w[v])) % p
        self.restricted_sums = 0
        self._grp = {}

    def series_vectors(self, t):
        """S[i] (M, n): coefficient of z^i in P_{f,v}(t+z)."""
        p, P, K = self.p, self.P, self.K
        S = np.zeros((self.J, P.shape[1], self.n), np.int64)
        for i in range(min(K, self.J)):
            acc = np.zeros(P.shape[1:], np.int64)
            for k in range(i, K):
                acc = (acc + (comb(k, i) * pow(t, k - i, p)) % p * P[k]) % p
            S[i] = acc
        return S

    def eval_clauses(self, inputs, which=None):
        """inputs: (J, N, M) series per occurrence; returns prod over clauses of Phi (J, N)."""
        p, J = self.p, self.J
        N = inputs.shape[1]
        total = np.zeros((J, N), np.int64)
        total[0] = 1
        for ci, lits in enumerate(self.cl_occ):
            if which is not None and ci not in which:
                continue
            prod = np.zeros((J, N), np.int64)
            prod[0] = 1
            for oi, pos in lits:
                s = inputs[:, :, oi]
                if pos:
                    oneminus = (-s) % p
                    oneminus[0] = (1 - s[0]) % p
                else:
                    oneminus = s
                prod = smul(prod, oneminus, p)
            phi = (-prod) % p
            phi[0] = (1 - prod[0]) % p
            total = smul(total, phi, p)
        return total

    def taylor_direct(self, t):
        p, S = self.p, self.series_vectors(t)
        inputs = np.stack([(self.X @ S[i].T) % p for i in range(self.J)])
        g = self.eval_clauses(inputs)
        return g.sum(axis=1) % p

    def g_U(self, S, U):
        p = self.p
        oneU = self.X[U]
        inputs = np.empty((self.J, self.N, S.shape[1]), np.int64)
        inputs[0] = (self.X @ S[0].T) % p
        for i in range(1, self.J):
            inputs[i] = ((oneU @ S[i].T) % p)[None, :]
        return self.eval_clauses(inputs)

    def taylor_inclusion_exclusion(self, t, Zfun=None):
        """Lemma 3.3: sum_{|E|<=q} sum_{U<=E} (-1)^{|E|-|U|} Z_{t,E,U}(z)."""
        p, N, n = self.p, self.N, self.n
        S = self.series_vectors(t)
        acc = np.zeros(self.J, np.int64)
        idx = np.arange(N)
        for U in range(N):
            if Zfun is None:
                g = self.g_U(S, U)
                Z = g.copy()
                for b in range(n):
                    m = idx[((idx >> b) & 1) == 0]
                    Z[:, m] = (Z[:, m] + Z[:, m | (1 << b)]) % p
            else:
                Z = Zfun(t, S, U)
            Es = idx[((idx & U) == U) & (self.pop <= self.q)]
            sign = np.where(((self.pop[Es] - self.pop[U]) % 2) == 0, 1, p - 1)
            acc = (acc + (Z[:, Es] * sign[None, :] % p).sum(axis=1)) % p
        return acc

    def Z_group_separated(self, t, S, U):
        """Z_{t,E,U} for all E at once, computed group by group with representatives."""
        p, aux, des = self.p, self.aux, self.des
        av, df, Chat = des["av"], des["df"], des["Chat"]
        if des.get("perpoint"):
            Chat = des["Chat_t"][t - 1]
        n, J = self.n, self.J
        groups = [np.nonzero(av[:, t - 1] == g)[0] for g in range(aux.r)]
        cl_group = [df[f][t - 1] for f in range(len(self.clauses))]
        oneU = self.X[U]
        corr = [((oneU @ S[i].T) % p) for i in range(1, J)]
        loc = []
        for g, vars_g in enumerate(groups):
            k = len(vars_g)
            lm = np.arange(1 << k)
            bits = ((lm[:, None] >> np.arange(k)[None, :]) & 1).astype(np.int64)
            gmask = (bits << vars_g[None, :]).sum(axis=1) if k else np.zeros(1, np.int64)
            Cg = Chat[:, vars_g]
            keys = [tuple(row) for row in (bits @ Cg.T)] if k else [()]
            fib = {}
            for i, key in enumerate(keys):
                fib.setdefault(key, []).append(int(gmask[i]))
            loc.append(fib)
        res = np.zeros((J, self.N), np.int64)
        stats = dict(combos=0, nontrivial=0)
        for combo in itertools.product(*[list(f.items()) for f in loc]):
            reps = [members[0] for _, members in combo]
            stats["combos"] += 1
            per_group = []
            for g, (key, members) in enumerate(combo):
                if len(members) > 1:
                    stats["nontrivial"] += 1
                others = sum(reps[h] for h in range(len(reps)) if h != g)
                xs = np.array(members, np.int64) + others
                Xs = ((xs[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
                inputs = np.empty((J, len(xs), S.shape[1]), np.int64)
                inputs[0] = (Xs @ S[0].T) % p
                for i in range(1, J):
                    inputs[i] = corr[i - 1][None, :]
                which = {ci for ci, gg in enumerate(cl_group) if gg == g}
                vals = self.eval_clauses(inputs, which)
                per_group.append((np.array(members, np.int64), vals))
            prodE = np.zeros((J, self.N), np.int64)
            prodE[0] = 1
            for g, (members, vals) in enumerate(per_group):
                gm = int(sum(1 << int(v) for v in groups[g]))
                Eg = np.arange(self.N) & gm
                sup = ((members[None, :] & Eg[:, None]) == Eg[:, None])
                sums = (vals @ sup.T.astype(np.int64)) % p
                prodE = smul(prodE, sums, p)
            res = (res + prodE) % p
        self.last_stats = stats
        return res

    # ------------------------------------------------ aggregated Lemma 3.3 (cost model v3)
    def group_structure(self, t):
        """Groups of point t, their masks/U masks, fibre members and fibre keys (cached)."""
        key = int(t)
        st = self._grp.get(key)
        if st is None:
            p, aux, des = self.p, self.aux, self.des
            av, df, Chat = des["av"], des["df"], des["Chat"]
            if des.get("perpoint"):
                Chat = des["Chat_t"][t - 1]
            n = self.n
            groups = [np.nonzero(av[:, t - 1] == g)[0] for g in range(aux.r)]
            cl_group = [df[f][t - 1] for f in range(len(self.clauses))]
            grp_mask, grp_bits, grp_Cg, loc = [], [], [], []
            for g, vars_g in enumerate(groups):
                k = len(vars_g)
                lm = np.arange(1 << k)
                bits = ((lm[:, None] >> np.arange(k)[None, :]) & 1).astype(np.int64)
                gmask = (bits << vars_g[None, :]).sum(axis=1) if k else np.zeros(1, np.int64)
                Cg = Chat[:, vars_g]
                if k:
                    vals = bits @ Cg.T
                    if self.residue_keys:
                        vals = vals % p
                    keys = [tuple(int(x) for x in row) for row in vals]
                else:
                    keys = [()]
                fib = {}
                for i, kk in enumerate(keys):
                    fib.setdefault(kk, []).append(int(gmask[i]))
                grp_mask.append(gmask)
                grp_bits.append(bits)
                grp_Cg.append(Cg)
                loc.append(fib)
            st = (groups, cl_group, grp_mask, grp_bits, grp_Cg, loc)
            self._grp[key] = st
        return st

    def fibre_counts(self, t):
        """Number of fibres (distinct measurement keys) of every group of point t."""
        return [len(fib) for fib in self.group_structure(t)[-1]]

    def taylor_aggregated(self, t, Zagg=None):
        """Aggregated Lemma 3.3: sum_{|U| <= q} sum_{x >= U} w_K(|x - U|) p_x(1_U), K = q - |U|."""
        p = self.p
        S = self.series_vectors(t)
        Zagg = self.Zagg_brute if Zagg is None else Zagg
        self.restricted_sums = 0
        acc = np.zeros(self.J, np.int64)
        for U in range(self.N):
            if self.pop[U] > self.q:
                continue
            self.restricted_sums += 1
            acc = (acc + Zagg(t, S, U)) % p
        return acc % p

    def Zagg_brute(self, t, S, U):
        """sum_{x >= U} w_K(popcount(x) - popcount(U)) g[:, x] mod p, with g = g_U(S, U)."""
        p, q = self.p, self.q
        K = q - int(self.pop[U])
        if K < 0:
            return np.zeros(self.J, np.int64)
        g = self.g_U(S, U)
        xs = np.arange(self.N, dtype=np.int64)
        sel = (xs & U) == U
        m = self.pop[sel] - int(self.pop[U])
        w = np.zeros(int(m.size), np.int64)
        for j in range(K + 1):
            w = (w + ((-1) ** j) * np.array([comb(int(mm), j) for mm in m], np.int64)) % p
        w = (w * self.wX[sel]) % p          # weight prod_v w_v^(x_v) of each assignment x
        return (g[:, sel] * w[None, :]).sum(axis=1) % p

    def Zagg_group_separated(self, t, S, U):
        """Same value as Zagg_brute, computed group by group like Z_group_separated, but one
        length-J vector: the count of ones of x_g outside U_g is tracked by a polynomial in y."""
        p, aux = self.p, self.aux
        n, J, q = self.n, self.J, self.q
        K = q - int(self.pop[U])
        if K < 0:
            return np.zeros(J, np.int64)
        groups, cl_group, grp_mask, grp_bits, grp_Cg, loc = self.group_structure(t)
        oneU = self.X[U]
        corr = [((oneU @ S[i].T) % p) for i in range(1, J)]
        Ky = K + 1
        Ugrp = [int(U & int(sum(1 << int(v) for v in vars_g))) if len(vars_g) else 0
                for vars_g in groups]
        popg = [int(bin(u).count("1")) for u in Ugrp]
        binom_tab = np.array([[comb(c, j) % p for j in range(Ky)] for c in range(n + 1)], np.int64)
        total = np.zeros(J, np.int64)
        for combo in itertools.product(*[list(f.items()) for f in loc]):
            reps = [members[0] for _, members in combo]
            per_group = []
            for g, (key, members) in enumerate(combo):
                xs = np.array(members, np.int64) + sum(reps[h] for h in range(len(reps)) if h != g)
                Xs = ((xs[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
                inputs = np.empty((J, len(xs), S.shape[1]), np.int64)
                inputs[0] = (Xs @ S[0].T) % p
                for i in range(1, J):
                    inputs[i] = corr[i - 1][None, :]
                which = {ci for ci, gg in enumerate(cl_group) if gg == g}
                vals = self.eval_clauses(inputs, which)
                ug = Ugrp[g]
                mem = np.array(members, np.int64)
                keep = (mem & ug) == ug
                if not keep.any():
                    per_group.append(None)
                    break
                mk = mem[keep]
                # weight of x_g: members only carry bits of group g
                vals = (vals[:, keep] * self.wX[mk][None, :]) % p
                # members only carry bits of group g, so this is the number of ones outside U_g
                c = np.array([int(bin(int(mm)).count("1")) for mm in mk], np.int64) - popg[g]
                Pg = np.einsum("im,mj->ij", vals, binom_tab[c]) % p  # (J, Ky)
                per_group.append(Pg)
            if any(Pg is None for Pg in per_group):
                continue
            prod = np.zeros((J, Ky), np.int64)
            prod[0, 0] = 1
            for Pg in per_group:
                prod = mul_zy(prod, Pg, p)
            signs = np.array([(-1) ** j for j in range(Ky)], np.int64)
            total = (total + (prod * signs[None, :]).sum(axis=1)) % p
        return total % p

    # ------------------------------------------------ contraction (cost model v7)
    def contract_structure(self, t):
        """Half-bond keys of every ordered pair (g, h) of groups of point t, and the
        representative mask of every half-bond value (cached). keys[g][h][i] is the value of
        the rows of the pair {g, h} restricted to group g, read on mask i of group g;
        reps[g][h][value] is the first mask of group g with that value."""
        cache = getattr(self, "_ctr", None)
        if cache is None:
            cache = self._ctr = {}
        ck = int(t)
        st = cache.get(ck)
        if st is None:
            p, aux, des = self.p, self.aux, self.des
            av = des["av"]
            Cm = np.asarray(des["Chat_t"][t - 1] if des.get("perpoint") else des["Chat"],
                            np.int64)
            groups, cl_group, grp_mask, grp_bits, grp_Cg, loc = self.group_structure(t)
            r = aux.r
            touch = np.zeros((r, Cm.shape[0]), bool)
            for g in range(r):
                if len(groups[g]):
                    touch[g] = np.any(Cm[:, groups[g]] != 0, axis=1)
            nt = touch.astype(np.int64).sum(axis=0)
            keys = []
            reps = []
            for g in range(r):
                kg, rg = [], []
                for h in range(r):
                    if h == g:
                        kg.append(None)
                        rg.append(None)
                        continue
                    if not len(groups[g]):                      # k = 0: the only mask is 0
                        kg.append([()])
                        rg.append({(): int(grp_mask[g][0])})
                        continue
                    both = (nt - touch[g].astype(np.int64) - touch[h].astype(np.int64)) == 0
                    rows_gh = np.nonzero(touch[g] & both)[0]
                    sub = Cm[np.ix_(rows_gh, groups[g])]
                    if self.residue_keys:
                        sub = sub % p
                    vals = grp_bits[g] @ sub.T                  # (2^k, #rows_gh)
                    kk = [tuple(int(x) for x in row) for row in vals]
                    rr = {}
                    for i, v in enumerate(kk):
                        if v not in rr:
                            rr[v] = int(grp_mask[g][i])
                    kg.append(kk)
                    rg.append(rr)
                keys.append(kg)
                reps.append(rg)
            st = (keys, reps)
            cache[ck] = st
        return st

    def Zagg_contracted(self, t, S, U):
        """Same length-J vector as Zagg_group_separated, computed by summing the groups one at
        a time (variable elimination): the table is keyed by the half-bonds between the
        processed and the unprocessed groups."""
        p, aux = self.p, self.aux
        n, J, q = self.n, self.J, self.q
        K = q - int(self.pop[U])
        if K < 0:
            return np.zeros(J, np.int64)
        groups, cl_group, grp_mask, grp_bits, grp_Cg, loc = self.group_structure(t)
        keys, reps = self.contract_structure(t)
        r = aux.r
        oneU = self.X[U]
        corr = [((oneU @ S[i].T) % p) for i in range(1, J)]
        Ky = K + 1
        Ugrp = [int(U & int(sum(1 << int(v) for v in vars_g))) if len(vars_g) else 0
                for vars_g in groups]
        popg = [int(bin(u).count("1")) for u in Ugrp]
        binom_tab = np.array([[comb(c, j) % p for j in range(Ky)] for c in range(n + 1)],
                             np.int64)
        cl_of_group = [{ci for ci, gg in enumerate(cl_group) if gg == g} for g in range(r)]
        unit = np.zeros((J, Ky), np.int64)
        unit[0, 0] = 1
        T = {(): unit}
        for g in range(r):
            future = list(range(g + 1, r))
            gm = grp_mask[g]
            keep = [i for i in range(len(gm)) if (int(gm[i]) & Ugrp[g]) == Ugrp[g]]
            which = cl_of_group[g]
            newT = {}
            for key, val in T.items():
                kd = {(a, b): (pi_ab, pi_ba) for ((a, b), pi_ab, pi_ba) in key}
                cons = [i for i in keep
                        if all(keys[g][a][i] == kd[(a, g)][1] for a in range(g))]
                if not cons:
                    continue
                past = 0
                for a in range(g):
                    past += reps[a][g][kd[(a, g)][0]]
                rest = tuple(e for e in key if e[0][1] != g)
                comb_list = [list(reps[h][g].items()) for h in future]
                for combo in itertools.product(*comb_list):
                    xs = gm[np.asarray(cons, np.int64)] + past
                    for _vh, mh in combo:
                        xs = xs + int(mh)
                    Xs = ((xs[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
                    inputs = np.empty((J, len(xs), S.shape[1]), np.int64)
                    inputs[0] = (Xs @ S[0].T) % p
                    for i in range(1, J):
                        inputs[i] = corr[i - 1][None, :]
                    vals = self.eval_clauses(inputs, which)
                    for jj, i in enumerate(cons):
                        m = int(gm[i])
                        col = (vals[:, jj] * self.wX[m]) % p
                        c = int(bin(m).count("1")) - popg[g]
                        Pg = (col[:, None] * binom_tab[c][None, :]) % p
                        newkey = tuple(sorted(rest + tuple(
                            ((g, h), keys[g][h][i], vh)
                            for h, (vh, _m) in zip(future, combo))))
                        add = mul_zy(val, Pg, p)
                        cur = newT.get(newkey)
                        newT[newkey] = add if cur is None else (cur + add) % p
            T = newT
            if not T:
                return np.zeros(J, np.int64)
        total = T[()]
        signs = np.array([(-1) ** j for j in range(Ky)], np.int64)
        return (total * signs[None, :]).sum(axis=1) % p

    # ------------------------------------------------ edge bits (cost model v9)
    def taylor_edges(self, t):
        """Same length-J vector as taylor_aggregated(t), computed inside the group-by-group
        contraction of 'contract': one bit per edge instead of one restricted sum per set U.

        At point t a replaced occurrence oi (des['desig'][oi] False, or no 'desig') whose variable
        v = des['occ'][oi][1] lies in group a = des['av'][v][t - 1] different from the clause group
        d = des['df'][f][t - 1] is an edge: the series rows S[k][oi] are zero outside the variables
        of the groups d and a. A-part: SA[0][oi] = S[0][oi], SA[k][oi] = S[k][oi] restricted to the
        variables of group d for k >= 1; B-part: SB[0][oi] = 0, SB[k][oi] = S[k][oi] restricted to
        the variables of group a for k >= 1. Every clause is multilinear in its occurrences, so
        H(t + z) mod z^J is the sum over the sets E of edges with |E| < J of sum_x w(x) prod_f D_f
        prod_{oi in E} B_oi, with D_f = 1 - prod_l u_l(A_l) if f has no edge in E and
        D_f = -(prod_{oi in E cap f} u'_oi) prod_{l in f, l not in E} u_l(A_l) otherwise, where
        u(y) = 1 - y for a positive literal and y for a negative one, u' = -1 (positive) / +1
        (negative). Every summand factorises over the groups, so the contraction carries one bit
        per open edge: a state is (key, O, c) with key the half-bond key of Zagg_contracted, O the
        edges with bit 1 that have exactly one processed end, c the number of bit-1 edges decided.
        """
        p, aux, des = self.p, self.aux, self.des
        n, J = self.n, self.J
        groups, cl_group, grp_mask, grp_bits, grp_Cg, loc = self.group_structure(t)
        keys, reps = self.contract_structure(t)
        r = aux.r
        S = self.series_vectors(t)
        # ---------------- edges, SA, SB (masked rows of length n) ----------------
        desig = des.get("desig")
        edges = []                       # e -> (oi, d, a), a != d
        e_of_occ = {}
        for oi, (f, v) in enumerate(des["occ"]):
            if f >= len(self.clauses):
                continue                 # singleton factors carry no clause factor
            if desig is not None and desig[oi]:
                continue                 # designated: constant row, no replacement
            d = int(des["df"][f][t - 1])
            a = int(des["av"][v][t - 1])
            if a == d:
                continue
            e_of_occ[oi] = len(edges)
            edges.append((oi, d, a))
        gvars = [np.zeros(n, bool) for _ in range(r)]
        for g in range(r):
            if len(groups[g]):
                gvars[g][groups[g]] = True
        SA, SB = [], []
        for (oi, d, a) in edges:
            A = np.zeros((J, n), np.int64)
            B = np.zeros((J, n), np.int64)
            A[0] = S[0][oi]
            for k in range(1, J):
                A[k] = S[k][oi] * gvars[d]
                B[k] = S[k][oi] * gvars[a]
            SA.append(A)
            SB.append(B)
        cl_edges = [[] for _ in self.clauses]
        for e, (oi, _d, _a) in enumerate(edges):
            cl_edges[int(des["occ"][oi][0])].append(e)
        new_e, close_e, b_side = [], [], []
        for g in range(r):
            fut, past = set(range(g + 1, r)), set(range(g))
            ng, cg, bg = [], [], []
            for e, (oi, d, a) in enumerate(edges):
                if (d == g and a in fut) or (a == g and d in fut):
                    ng.append(e)
                if (d == g and a in past) or (a == g and d in past):
                    cg.append(e)
                if a == g:
                    bg.append(e)
            new_e.append(ng)
            close_e.append(cg)
            b_side.append(bg)
        cl_of_group = [{ci for ci, gg in enumerate(cl_group) if gg == g} for g in range(r)]
        # rows of S / SA / SB against the own masks are recomputed per (key, cons, combo)
        unit = np.zeros(J, np.int64)
        unit[0] = 1
        T = {((), frozenset(), 0): unit.copy()}
        for g in range(r):
            future = list(range(g + 1, r))
            gm = grp_mask[g]
            cl_g = sorted(cl_of_group[g])
            ng, cg, bg = new_e[g], close_e[g], b_side[g]
            csets = [[]]
            for k in range(1, len(ng) + 1):
                csets += [list(cn) for cn in itertools.combinations(ng, k)]
            cgset = set(cg)
            newT = {}
            ctxcache = {}
            for (key, O, c), val in T.items():
                kd = {(hh, bb): (pi_ab, pi_ba) for ((hh, bb), pi_ab, pi_ba) in key}
                cons = [i for i in range(len(gm))
                        if all(keys[g][a][i] == kd[(a, g)][1] for a in range(g))]
                if not cons:
                    continue
                past = 0
                for a in range(g):
                    past += reps[a][g][kd[(a, g)][0]]
                rest = tuple(e2 for e2 in key if e2[0][1] != g)
                ck = (key, tuple(cons), past)
                ctx = ctxcache.get(ck)
                if ctx is None:
                    ctx = []
                    for cmb in itertools.product(
                            *[list(reps[h][g].items()) for h in future]):
                        xs = gm[np.asarray(cons, np.int64)] + past
                        for _vh, mh in cmb:
                            xs = xs + int(mh)
                        Xs = ((xs[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
                        ctx.append([cmb, Xs, self.wX[gm[cons]].copy(), {}, {}, {}, {}])
                    ctxcache[ck] = ctx
                lim = J - 1 - c
                if lim < 0:
                    continue
                for cmb, Xs, wsel, serA, serB, Dc, Tc in ctx:
                    N = Xs.shape[0]
                    nk = [tuple(sorted(rest + tuple(
                        ((g, h), keys[g][h][cons[jj]], vh)
                        for h, (vh, _m) in zip(future, cmb)))) for jj in range(len(cons))]
                    for N1 in [s for s in csets if len(s) <= lim]:
                        n1 = set(N1)
                        pats = tuple(tuple(1 if ((e in O) or (e in n1)) else 0
                                          for e in cl_edges[ci]) for ci in cl_g)
                        bset = frozenset(e for e in bg if ((e in O) or (e in n1)))
                        total = Tc.get((pats, bset))
                        if total is None:
                            total = np.zeros((J, N), np.int64)
                            total[0] = 1
                            for k, ci in enumerate(cl_g):
                                pat = pats[k]
                                D = Dc.get((ci, pat))
                                if D is None:
                                    prod = np.zeros((J, N), np.int64)
                                    prod[0] = 1
                                    for oi, pos in self.cl_occ[ci]:
                                        e = e_of_occ.get(oi)
                                        on = (e is not None) and ((e in O) or (e in n1))
                                        if on:
                                            s = np.zeros((J, N), np.int64)
                                            s[0] = p - 1 if pos else 1  # u' = -1 / +1
                                        else:
                                            y = serA.get(oi)
                                            if y is None:
                                                if e is None:
                                                    y = (S[:, oi] @ Xs.T) % p  # full rows
                                                else:
                                                    y = (SA[e] @ Xs.T) % p     # A-part rows
                                                serA[oi] = y
                                            if pos:
                                                s = (-y) % p
                                                s[0] = (1 - y[0]) % p
                                            else:
                                                s = y
                                        prod = smul(prod, s, p)
                                    D = (-prod) % p
                                    if not any(pat):
                                        D[0] = (1 - prod[0]) % p
                                    Dc[(ci, pat)] = D
                                total = smul(total, D, p)
                            for e in bset:
                                B = serB.get(e)
                                if B is None:
                                    B = (SB[e] @ Xs.T) % p
                                    serB[e] = B
                                total = smul(total, B, p)
                            total = (total * wsel[None, :]) % p
                            Tc[(pats, bset)] = total
                        O2 = (O - cgset) | frozenset(N1)
                        c2 = c + len(N1)
                        for jj in range(len(cons)):
                            add = np.convolve(val, total[:, jj])[:J] % p
                            st = (nk[jj], O2, c2)
                            cur = newT.get(st)
                            newT[st] = add if cur is None else (cur + add) % p
            T = newT
            if not T:
                return np.zeros(J, np.int64)
        acc = np.zeros(J, np.int64)
        for (key, O, c), val in T.items():
            if not O:
                acc = (acc + val.reshape(J)) % p
        return acc % p

    def hermite_H0(self, data):
        p, Q, J = self.p, self.Q, self.J
        L = Q * J
        A, b = [], []
        for t in range(1, Q + 1):
            for i in range(J):
                A.append([(comb(u, i) * pow(t, u - i, p)) % p if u >= i else 0 for u in range(L)])
                b.append(int(data[t][i]) % p)
        M = [row + [bb] for row, bb in zip(A, b)]
        for c in range(L):
            piv = next(r for r in range(c, L) if M[r][c] % p)
            M[c], M[piv] = M[piv], M[c]
            inv = pow(M[c][c], p - 2, p)
            M[c] = [(x * inv) % p for x in M[c]]
            for r in range(L):
                if r != c and M[r][c]:
                    f = M[r][c]
                    M[r] = [(x - f * y) % p for x, y in zip(M[r], M[c])]
        return M[0][L]
