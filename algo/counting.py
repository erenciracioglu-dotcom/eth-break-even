"""End-to-end run of the Section 2-3 counting algorithm on small formulas, mod a prime p.

H(T) = sum_x prod_f Phi_f(<P_{f,v}(T), x>),  H(0) = #F.
For t = 1..Q the Taylor data H(t+z) mod z^J are computed via Lemma 3.3 (inclusion-exclusion
over U <= E, |E| <= q); each restricted sum Z_{t,E,U} can be computed by Lemma 3.2 (fix
measurement tuples, sum every group separately using representatives), which is the
algorithm of the paper. Hermite interpolation of the Taylor data then gives H(0) mod p.

PUBLIC API (the referee calls these; keep names):
  Counter(aux, des).taylor_inclusion_exclusion(t, Zfun=None), .Z_group_separated(t, S, U),
  .hermite_H0(data), .Q;  brute_count(n, clauses)
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
    def __init__(self, aux, des, p=PRIMES[0]):
        self.aux, self.des, self.p = aux, des, p
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
            # designated occurrences carry no replacement row: only |I_f| - 1 per clause count
            d = sum(len(cl) - 1 for cl in clauses)
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
