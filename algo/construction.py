"""Exact implementation of the replacement-form construction (Section 4) of
"The Exponential Time Hypothesis is False", plus checks of Lemmas 4.1-4.3,
the transfer to actual formulas (Sections 4.4-4.5) and Definition 3.1(4).

Integer arithmetic throughout: M0 is stored as (R-1)! * M0, the rows P_u(T)
as Lden * P_u(T) with Lden = (R-1)! * Q!.

PUBLIC API (the referee calls these; keep names and return fields):
  Aux(r, Q, R, A0).report() -> dict with f_interp, eigen, orthogonal, support, P0_delta,
                               degree, separation, min_kernel_dim
  random_formula(n, m, rng), random_bounded_formula(n, rng)
  build_design(aux, n, clauses, rng) -> dict (keys used: Prow, Chat, av, df, occ, h, n, M, clauses)
  check_design(aux, des) -> dict (keys used: kernel_condition, P0_delta, min_kernel_dim,
                               max_groups_per_row, max_restricted_abs)
"""
import itertools, math
from math import factorial
import numpy as np


# ------------------------------------------------------------ small helpers
def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def primes_below(limit, k):
    out, x = [], limit
    while len(out) < k:
        x -= 1
        if is_prime(x):
            out.append(x)
    return out


PRIMES = primes_below(1 << 25, 2)      # products of two residues fit easily in int64


def gen_binom(x, i):
    """Generalised binomial C(x, i) for integer x (possibly negative)."""
    num = 1
    for a in range(i):
        num *= (x - a)
    return num // factorial(i)


def f_I(j, R, A0):
    """f_I(j) from (4.7)."""
    return (-1) ** A0 * sum((-2) ** i * gen_binom(j - A0, i) for i in range(R))


def rref_mod(A, p):
    A = np.array(A % p, dtype=np.int64)
    nr, nc = A.shape
    piv, r = [], 0
    for c in range(nc):
        if r >= nr:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + int(nz[0])
        if i != r:
            A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        piv.append(c)
        r += 1
    return A[:r], piv


def matmul_mod(A, B, p, block=2048):
    out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    for s in range(0, A.shape[1], block):
        out = (out + (A[:, s:s + block] @ B[s:s + block]) % p) % p
    return out


def in_rowspace_mod(C, B, p):
    """Every row of B lies in the row space of C (mod p)? Returns (ok, rank(C))."""
    Rr, piv = rref_mod(C, p)
    Bm = np.array(B % p, dtype=np.int64)
    if not piv:
        return (not np.any(Bm)), 0
    res = (Bm - matmul_mod(Bm[:, piv], Rr, p)) % p
    return (not np.any(res)), len(piv)


# ------------------------------------------------------- auxiliary construction
def local_basis(r):
    """Eigenbasis of the swap on [r]x[r]: (vector over r*r, eigenvalue sign)."""
    out = []
    for al in range(r):
        v = np.zeros(r * r, np.int64)
        v[al * r + al] = 1
        out.append((v, +1))
    for al in range(r):
        for de in range(al + 1, r):
            v = np.zeros(r * r, np.int64); v[al * r + de] = 1; v[de * r + al] = 1
            w = np.zeros(r * r, np.int64); w[al * r + de] = 1; w[de * r + al] = -1
            out.append((v, +1))
            out.append((w, -1))
    return out


class Aux:
    """B = ([r]x[r])^Q, M0 = V0 f_I(W0), C^B and rows P_u(T) (Sections 4.1-4.3)."""

    def __init__(self, r, Q, R, A0):
        self.r, self.Q, self.R, self.A0 = r, Q, R, A0
        D = r ** (2 * Q)
        self.D = D
        idx = np.arange(D, dtype=np.int64)
        a = np.empty((D, Q), np.int64)
        d = np.empty((D, Q), np.int64)
        rem = idx.copy()
        for s in range(Q):
            pr = rem % (r * r)
            rem //= (r * r)
            a[:, s] = pr // r
            d[:, s] = pr % r
        self.a, self.d = a, d
        self.base = (r * r) ** np.arange(Q, dtype=np.int64)
        perms = []
        for s in range(Q):
            A2, D2 = a.copy(), d.copy()
            A2[:, s], D2[:, s] = d[:, s], a[:, s]
            perms.append(self.index_of(A2, D2))
        self.perms = perms
        self.full = self.index_of(d, a)                 # V0: swap every position

        # (R-1)! * f_I(W0) = sum_j (-1)^(A0+j) ((R-1)!/j!) prod_{a<j} (2W0 - 2(A0+a) I)
        fac = factorial(R - 1)
        self.fac = fac
        Pj = np.eye(D, dtype=np.int64)
        F = np.zeros((D, D), np.int64)
        for j in range(R):
            F += ((-1) ** (A0 + j)) * (fac // factorial(j)) * Pj
            if j < R - 1:
                twoW = Q * Pj - sum(Pj[pm, :] for pm in perms)   # (2 W0) Pj
                Pj = twoW - 2 * (A0 + j) * Pj
        self.M0s = F[self.full, :]                      # (R-1)! * M0
        assert np.abs(self.M0s).max() < 2 ** 40

        # product eigenbasis psi and the measurement matrix C^B
        loc = local_basis(r)
        pairs = a * r + d
        Psi = np.empty((D, D), np.int64)
        js = np.empty(D, np.int64)
        for k, choice in enumerate(itertools.product(range(len(loc)), repeat=Q)):
            v = np.ones(D, np.int64)
            j = 0
            for s, c in enumerate(choice):
                vec, sg = loc[c]
                v *= vec[pairs[:, s]]
                j += sg < 0
            Psi[k] = v
            js[k] = j
        self.Psi, self.js = Psi, js
        self.inI = (js >= A0) & (js <= A0 + R - 1)
        self.Cmeas = Psi[~self.inI]

        # rows P_u(T), scaled by Lden = (R-1)! Q!   -- equation (4.10)
        QF = factorial(Q)
        self.Lden = fac * QF
        nz_u, nz_y = np.nonzero(self.M0s)
        mis = self.a[nz_y] != self.d[nz_u]
        self.max_mismatch = int(mis.sum(axis=1).max())
        Pc = np.zeros((max(R, self.max_mismatch + 1), D, D), np.int64)
        Pc[0] = self.Lden * np.eye(D, dtype=np.int64) - QF * self.M0s
        masks = (mis * (1 << np.arange(Q))).sum(axis=1)
        for m in np.unique(masks):
            S = [s + 1 for s in range(Q) if (m >> s) & 1]
            coeffs = [1]
            for s in S:                                   # multiply by (s - T)
                new = [0] * (len(coeffs) + 1)
                for i, c in enumerate(coeffs):
                    new[i] += c * s
                    new[i + 1] -= c
                coeffs = new
            scale = QF // math.prod(S) if S else QF       # QF * prod (1 - T/s)
            sel = masks == m
            uu, yy = nz_u[sel], nz_y[sel]
            vals = self.M0s[uu, yy]
            for k, c in enumerate(coeffs):
                if c:
                    np.add.at(Pc[k], (uu, yy), vals * (c * scale))
        self.Pc = Pc

    def index_of(self, A, Dd):
        return ((A * self.r + Dd) * self.base).sum(axis=-1)

    # ------------------------------------------------------------ checks
    def check_f(self):
        return all(f_I(j, self.R, self.A0) == (-1) ** j for j in range(self.A0, self.A0 + self.R))

    def check_eigen(self):
        """M0 psi = (-1)^j f_I(j) psi for every basis vector; basis orthogonal."""
        lhs = self.M0s @ self.Psi.T
        fv = np.array([(-1) ** int(j) * f_I(int(j), self.R, self.A0) for j in self.js], dtype=np.int64)
        rhs = self.Psi.T * (self.fac * fv)[None, :]
        G = self.Psi @ self.Psi.T
        return bool(np.array_equal(lhs, rhs)), bool(np.count_nonzero(G - np.diag(np.diag(G))) == 0)

    def check_support(self):
        """Lemma 4.2: nonzero M0[u,y] => #{s: a_y(s) != d_u(s)} <= R-1."""
        return self.max_mismatch <= self.R - 1

    def check_rows(self):
        """P_u(0) = delta_u and deg P_u <= R-1."""
        ok0 = bool(np.array_equal(self.Pc[0], self.Lden * np.eye(self.D, dtype=np.int64)))
        okdeg = all(not np.any(self.Pc[k]) for k in range(self.R, self.Pc.shape[0]))
        return ok0, okdeg

    def Pt(self, t):
        return sum((t ** k) * self.Pc[k] for k in range(self.Pc.shape[0]))

    def check_separation(self):
        """Lemma 4.3: for g != d_u(t), P_u(t)|S_{t,g} lies in rowspace(C^B restricted to S_{t,g})."""
        ok, kdims = True, []
        for p in PRIMES:
            for t in range(1, self.Q + 1):
                P = self.Pt(t)
                for g in range(self.r):
                    cols = np.nonzero(self.a[:, t - 1] == g)[0]
                    rows = np.nonzero(self.d[:, t - 1] != g)[0]
                    good, rank = in_rowspace_mod(self.Cmeas[:, cols], P[np.ix_(rows, cols)], p)
                    ok &= good
                    kdims.append(len(cols) - rank)
        return ok, min(kdims)

    def report(self):
        e, orth = self.check_eigen()
        r0, rdeg = self.check_rows()
        sep, kmin = self.check_separation()
        return dict(r=self.r, Q=self.Q, R=self.R, A0=self.A0, D=self.D, cB=int(self.Cmeas.shape[0]),
                    f_interp=self.check_f(), eigen=e, orthogonal=orth, support=self.check_support(),
                    P0_delta=r0, degree=rdeg, separation=sep, min_kernel_dim=int(kmin))


# ------------------------------------------------------ formulas and designs
def random_bounded_formula(n, rng, occ=3):
    """Clauses of width <= 3, every variable in exactly `occ` clauses (repeats in a clause dropped)."""
    slots = [v for v in range(n) for _ in range(occ)]
    rng.shuffle(slots)
    clauses = []
    for i in range(0, len(slots), 3):
        vs = list(dict.fromkeys(slots[i:i + 3]))
        clauses.append([(int(v), bool(rng.random() < 0.5)) for v in vs])
    return clauses


def random_formula(n, m, rng, k=3):
    clauses = []
    for _ in range(m):
        vs = rng.choice(n, size=k, replace=False)
        clauses.append([(int(v), bool(rng.random() < 0.5)) for v in vs])
    return clauses


def build_design(aux, n, clauses, rng):
    """Sections 4.4-4.5: schedules, labels, copies, rows P_{f,v}(T), matrices C_{t,g}."""
    r, Q, R = aux.r, aux.Q, aux.R
    used = {v for cl in clauses for v, _ in cl}
    factors = [[v for v, _ in cl] for cl in clauses] + [[v] for v in range(n) if v not in used]
    av = rng.integers(0, r, size=(n, Q))
    df = rng.integers(0, r, size=(len(factors), Q))
    occ = [(f, v) for f, sc in enumerate(factors) for v in sc]
    labels = np.array([int(aux.index_of(av[v], df[f])) for f, v in occ], dtype=np.int64)
    copy = np.zeros(len(occ), np.int64)
    cnt = {}
    for i, u in enumerate(labels):
        c = cnt.get(int(u), 0)
        copy[i] = c
        cnt[int(u)] = c + 1
    h = max(cnt.values())
    pos_y = [[] for _ in range(h)]
    pos_v = [[] for _ in range(h)]
    for i, (f, v) in enumerate(occ):
        pos_y[copy[i]].append(labels[i])
        pos_v[copy[i]].append(v)
    pos_y = [np.array(x, np.int64) for x in pos_y]
    pos_v = [np.array(x, np.int64) for x in pos_v]
    K = aux.Pc.shape[0]
    Prow = np.zeros((K, len(occ), n), np.int64)
    for i in range(len(occ)):
        c = copy[i]
        for k in range(K):
            np.add.at(Prow[k, i], pos_v[c], aux.Pc[k][labels[i], pos_y[c]])
    cB = aux.Cmeas.shape[0]
    Chat = np.zeros((h * cB, n), np.int64)
    for c in range(h):
        block = np.zeros((n, cB), np.int64)
        np.add.at(block, pos_v[c], aux.Cmeas[:, pos_y[c]].T)
        Chat[c * cB:(c + 1) * cB] = block.T
    return dict(factors=factors, av=av, df=df, occ=occ, labels=labels, copy=copy, h=h,
                Prow=Prow, Chat=Chat, n=n, M=len(occ), clauses=clauses)


def check_design(aux, des):
    """Definition 3.1(4) on the actual variables, plus the refined accounting claims."""
    n, Q, r = des["n"], aux.Q, aux.r
    Prow, Chat, av, df, occ = des["Prow"], des["Chat"], des["av"], des["df"], des["occ"]
    ok = True
    kdims = []
    for p in PRIMES:
        for t in range(1, Q + 1):
            Pt = sum((t ** k) * Prow[k] for k in range(Prow.shape[0]))
            for g in range(r):
                cols = np.nonzero(av[:, t - 1] == g)[0]
                rows = [i for i, (f, v) in enumerate(occ) if df[f][t - 1] != g]
                good, rank = in_rowspace_mod(Chat[:, cols], Pt[np.ix_(rows, cols)], p)
                ok &= good
                kdims.append(len(cols) - rank)
    # P_{f,v}(0) = delta_v
    E = np.zeros((len(occ), n), np.int64)
    E[np.arange(len(occ)), [v for _, v in occ]] = aux.Lden
    p0 = bool(np.array_equal(Prow[0], E))
    # refined accounting: groups touched per row, restricted abs sums
    max_groups, max_restr, nonzero_rows = 0, 0, 0
    for t in range(Q):
        touched = np.zeros(Chat.shape[0], np.int64)
        for g in range(r):
            cols = av[:, t] == g
            sub = Chat[:, cols]
            nzr = np.any(sub != 0, axis=1)
            touched += nzr
            nonzero_rows += int(nzr.sum())
            if sub.size:
                max_restr = max(max_restr, int(np.abs(sub).sum(axis=1).max()))
        max_groups = max(max_groups, int(touched.max()))
    sizes = [int((av[:, t] == g).sum()) for t in range(Q) for g in range(r)]
    return dict(n=n, M=des["M"], D=aux.D, h=des["h"], cB=int(aux.Cmeas.shape[0]),
                kernel_condition=bool(ok), P0_delta=p0, min_kernel_dim=int(min(kdims)),
                max_group_size=max(sizes),
                rows_paper_count=int(r * des["h"] * aux.Cmeas.shape[0]),       # r*h*cB per point
                rows_nonzero_avg=nonzero_rows / Q,                               # actual, per point
                max_groups_per_row=max_groups, max_restricted_abs=max_restr,
                bound_restricted=2 ** (Q - 1), paper_bound_abs=aux.D,
                max_full_abs=int(np.abs(Chat).sum(axis=1).max()))
