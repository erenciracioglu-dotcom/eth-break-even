#!/usr/bin/env python3
"""Toy implementation of the Appendix D contraction (Theorem D.1) and the head's extensions.

Everything is exact (Python integers and Fractions). Checks:
  1. networks: the clause-variable network and the chunk super-vertex network contract to #F
  2. packing (Lemma D.3 (1)-(2)); W(0) = I and W(t) is local for coordinate v2(t) (D.4.5)
  3. F(T) (D.8) has degree <= M(2L-1), F(0) = contraction, F(t + alpha z) in Z_(2)[z] (Lemma D.4)
  4. tapered Hermite recovery (Lemma D.5): Taylor data reduced mod 2^P_t give F(0) mod 2^p
  5. modes (D.12)-(D.16): F(t + alpha z) = sum_beta prod_g Z_{g,beta} mod (2^p, z^q), computed
     (a) per tuple as in the paper, (b) with each distinct child once ('memo'), (c) group by group
     ('modedp'); child networks built as in D.4.6 contract to the same values
  6. end to end (a = 1): #F from the mode sums at the precision of Lemma D.5
Weights: 'paper' w = b + 1, or 'sharp' w = a + 1 - i + b (every crossing entry of W(t + alpha z)
is divisible by 2^(a+1-i): C(t + alpha z)/C'(0) is divisible by alpha = 2^(a+1)).
Bases (round 6): 'power' (the paper: source forms sum_j j^b X_(y0 + c j), coefficients W_(x, y0) (c / A)^b) or
'newton' (binomial moments, suggested by GPT-6 Astra: source forms sum_j binom(j, b) X_(y0 + c j), coefficients
W_(x, y0) prod_(l=1..b) c l / (A - c l), the b-th forward difference in j; every weight grows by v2(b!)), or
'ldu' (Cauchy LDU, also suggested by GPT-6 Astra; sharp weights only): the crossing block of a packet half (slots
with bit i = h, ordered x_k = 2^(i+1) floor(k/2^i) + (k mod 2^i) + h 2^i) is K_kl = 1/(T + x_k - x_l) = L D U with unit
triangular L, U over Z_(2)[[z]]; mode m has the source form (U X)_m (series coefficients) and coefficient
G L_km d_m in row k, weight a + 1 - i + 2b + 2 v2(b!) with b = floor(m/2^i).
"""
import itertools
import math
import random
import sys
from fractions import Fraction as Fr


# ----------------------------------------------------------------------------- helpers
def v2(x):
    x = Fr(x)
    if x == 0:
        return math.inf
    n, d, c = x.numerator, x.denominator, 0
    while n % 2 == 0:
        n //= 2
        c += 1
    while d % 2 == 0:
        d //= 2
        c -= 1
    return c


def red(x, P):
    """2-adic integer Fraction -> residue mod 2^P."""
    x = Fr(x)
    n, d = x.numerator, x.denominator
    assert d % 2 == 1, "not 2-adically integral"
    m = 1 << P
    return (n * pow(d, -1, m)) % m


def poly_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                c[i + j] += ai * bj
    return c


def ser_mul(a, b, q):
    c = [Fr(0)] * q
    for i in range(q):
        if a[i]:
            for j in range(q - i):
                c[i + j] += a[i] * b[j]
    return c


def ser_inv(a, q):
    """1/a as a power series (a[0] != 0)."""
    inv = [Fr(0)] * q
    inv[0] = Fr(1) / a[0]
    for k in range(1, q):
        s = sum(a[j] * inv[k - j] for j in range(1, min(k, len(a) - 1) + 1))
        inv[k] = -s / a[0]
    return inv


def v2fact(b):
    """v2(b!) = b - (number of ones in b)."""
    return b - bin(b).count("1")


def lin_eval(lin, bit, q, m):
    """Value of a source form at column bits: an int when every coefficient is an int, else a series mod m."""
    if all(isinstance(c, int) for c, _ in lin):
        return sum(c * bit(w2) for c, w2 in lin)
    out = [0] * q
    for c, w2 in lin:
        if bit(w2):
            if isinstance(c, int):
                out[0] = (out[0] + c) % m
            else:
                out = [(x + y) % m for x, y in zip(out, c)]
    return out


def src_mul(f, v, e, q, m):
    """f * v^e for f, v ints or series mod m (the result stays an int while both are ints)."""
    if isinstance(f, int) and isinstance(v, int):
        return f * v ** e
    fs = [f % m] + [0] * (q - 1) if isinstance(f, int) else f
    vs = [v % m] + [0] * (q - 1) if isinstance(v, int) else v
    for _ in range(e):
        fs = ps_mul_mod(fs, vs, q, m)
    return fs


def src_scale(f, ser, q, m):
    """f * ser for an int or series f."""
    if isinstance(f, int):
        return [(f * c) % m for c in ser]
    return ps_mul_mod(f, ser, q, m)


def src_zero(f):
    return (f == 0) if isinstance(f, int) else not any(f)


def series_ldu(K, q):
    """K = L diag(d) U over Q[[z]] mod z^q (unit triangular L, U); entries are lists of q Fractions."""
    n = len(K)
    A = [[list(K[r][s]) for s in range(n)] for r in range(n)]
    unit = [Fr(1)] + [Fr(0)] * (q - 1)
    Lm = [[None] * n for _ in range(n)]
    Um = [[None] * n for _ in range(n)]
    d = []
    for k in range(n):
        piv = A[k][k]
        assert piv[0] != 0, "zero pivot"
        pinv = ser_inv(piv, q)
        d.append(piv)
        Lm[k][k] = Um[k][k] = unit
        for r in range(k + 1, n):
            Lm[r][k] = ser_mul(A[r][k], pinv, q)
        for s in range(k + 1, n):
            Um[k][s] = ser_mul(A[k][s], pinv, q)
        for r in range(k + 1, n):
            for s in range(k + 1, n):
                A[r][s] = [x - y for x, y in zip(A[r][s], ser_mul(Lm[r][k], A[k][s], q))]
    return Lm, d, Um


def ps_mul_mod(a, b, q, m):
    c = [0] * q
    for i in range(q):
        ai = a[i]
        if ai:
            for j in range(q - i):
                if b[j]:
                    c[i + j] = (c[i + j] + ai * b[j]) % m
    return c


def ps_add_mod(a, b, m):
    return [(x + y) % m for x, y in zip(a, b)]


# ----------------------------------------------------------------------------- formulas, networks
def count_sat(n, clauses):
    c = 0
    for x in range(1 << n):
        if all(any(((x >> (abs(l) - 1)) & 1) == (l > 0) for l in cl) for cl in clauses):
            c += 1
    return c


def clause_sat(cl, val):
    return any(val[abs(l)] == (1 if l > 0 else 0) for l in cl)


class Net:
    """wires: list of (row_vertex, col_vertex); vertex: dict(ports=[(wire, 'r'|'c')], table={bits: int})."""

    def __init__(self):
        self.wires, self.verts, self.unused = [], [], 0

    def contract(self):
        return self.contract_core() << self.unused

    def contract_core(self):
        M = len(self.wires)
        tot = 0
        for x in range(1 << M):
            bit = [(x >> w) & 1 for w in range(M)]
            p = 1
            for v in self.verts:
                p *= v["table"].get(tuple(bit[w] for w, _ in v["ports"]), 0)
                if p == 0:
                    break
            tot += p
        return tot


def formula_net(n, clauses):
    net = Net()
    occ = {v: [] for v in range(1, n + 1)}
    for ci, cl in enumerate(clauses):
        vs = sorted(set(abs(l) for l in cl))
        net.verts.append(dict(ports=[], table={}, kind=("clause", ci)))
        for v in vs:
            w = len(net.wires)
            net.wires.append([ci, None])
            net.verts[ci]["ports"].append((w, "r"))
            occ[v].append(w)
        tab = {}
        for bits in itertools.product((0, 1), repeat=len(vs)):
            if clause_sat(cl, dict(zip(vs, bits))):
                tab[bits] = 1
        net.verts[ci]["table"] = tab
    for v in range(1, n + 1):
        if not occ[v]:
            net.unused += 1
            continue
        k = len(net.verts)
        d = len(occ[v])
        net.verts.append(dict(ports=[(w, "c") for w in occ[v]],
                              table={(b,) * d: 1 for b in (0, 1)}, kind=("var", v)))
        for w in occ[v]:
            net.wires[w][1] = k
    net.wires = [tuple(w) for w in net.wires]
    return net


def make_units(n, clauses, Lc, rng):
    """Connected chunks of the variable graph with about Lc variables (greedy BFS)."""
    adj = {v: set() for v in range(1, n + 1)}
    for cl in clauses:
        vs = [abs(l) for l in cl]
        for a in vs:
            for b in vs:
                if a != b:
                    adj[a].add(b)
    used = sorted(set(abs(l) for cl in clauses for l in cl))
    left = set(used)
    units = []
    order = used[:]
    rng.shuffle(order)
    for s in order:
        if s not in left:
            continue
        U, frontier = [s], [s]
        left.discard(s)
        while frontier and len(U) < Lc:
            nxt = []
            for a in frontier:
                for b in sorted(adj[a]):
                    if b in left and len(U) < Lc:
                        U.append(b)
                        left.discard(b)
                        nxt.append(b)
            frontier = nxt
        units.append(set(U))
    return units


_WIREMAP = {}


def chunk_net_mapped(n, clauses, units):
    """chunk_net with the wire -> (clause, variable) map registered before the tables are built."""
    unit_of = {v: k for k, U in enumerate(units) for v in U}
    desig = []
    for cl in clauses:
        vs = sorted(set(abs(l) for l in cl))
        cnt = {}
        for v in vs:
            cnt[unit_of[v]] = cnt.get(unit_of[v], 0) + 1
        desig.append(min(cnt, key=lambda k: (-cnt[k], k)))
    net = Net()
    w = 0
    for ci, cl in enumerate(clauses):
        for v in sorted(set(abs(l) for l in cl)):
            if unit_of[v] != desig[ci]:
                _WIREMAP[(id(net), w)] = (ci, v)
                _WIREMAP[(id(net), "cv", ci, v)] = w
                w += 1
    # build with the same numbering
    tmp = chunk_net_build(n, clauses, units, unit_of, desig, net)
    return tmp


def chunk_net_build(n, clauses, units, unit_of, desig, net):
    net.units, net.unit_of, net.desig, net.clauses = units, unit_of, desig, clauses
    net.unused = n - len(unit_of)
    net.verts = [dict(ports=[], table={}, kind=("unit", k)) for k in range(len(units))]
    for ci, cl in enumerate(clauses):
        for v in sorted(set(abs(l) for l in cl)):
            if unit_of[v] != desig[ci]:
                w = len(net.wires)
                net.wires.append((desig[ci], unit_of[v]))
                net.verts[desig[ci]]["ports"].append((w, "r"))
                net.verts[unit_of[v]]["ports"].append((w, "c"))
    for k, U in enumerate(units):
        Us = sorted(U)
        ports = net.verts[k]["ports"]
        cls = [ci for ci in range(len(clauses)) if desig[ci] == k]
        tab = {}
        for pbits in itertools.product((0, 1), repeat=len(ports)):
            tot = 0
            for xb in itertools.product((0, 1), repeat=len(Us)):
                val = dict(zip(Us, xb))
                ok = True
                rowval = {}
                for (w, end), b in zip(ports, pbits):
                    if end == "c":
                        if val[_WIREMAP[(id(net), w)][1]] != b:
                            ok = False
                            break
                    else:
                        rowval[w] = b
                if not ok:
                    continue
                for ci in cls:
                    vv = {}
                    for v in sorted(set(abs(l) for l in clauses[ci])):
                        vv[v] = val[v] if v in U else rowval[_WIREMAP[(id(net), "cv", ci, v)]]
                    if not clause_sat(clauses[ci], vv):
                        ok = False
                        break
                if ok:
                    tot += 1
            if tot:
                tab[pbits] = tot
        net.verts[k]["table"] = tab
    return net


# ----------------------------------------------------------------------------- packing (Lemma D.3)
class Packing:
    def __init__(self, net, labels, a):
        self.a, self.L = a, 1 << a
        self.S = list(range(-self.L, self.L))
        by = {}
        for w, (rv, cv) in enumerate(net.wires):
            rl, cl = labels[rv], labels[cv]
            pat = tuple((min(rl[i], cl[i]), max(rl[i], cl[i])) for i in range(a))
            ori = tuple(0 if rl[i] <= cl[i] else 1 for i in range(a) if rl[i] != cl[i])
            by.setdefault(pat, {}).setdefault(ori, []).append(w)
        self.packets = []
        self.where = {}
        for pat, d in sorted(by.items()):
            uneq = [i for i in range(a) if pat[i][0] != pat[i][1]]
            eq = [i for i in range(a) if pat[i][0] == pat[i][1]]
            cap = 1 << len(eq)
            npk = max((len(ws) + cap - 1) // cap for ws in d.values())
            for k in range(npk):
                pk = dict(pat=pat, wire_at={})
                for ori, ws in sorted(d.items()):
                    for c, w in enumerate(ws[k * cap:(k + 1) * cap]):
                        s = 0
                        for idx, i in enumerate(uneq):
                            s |= ori[idx] << i
                        for idx, i in enumerate(eq):
                            s |= ((c >> idx) & 1) << i
                        assert s not in pk["wire_at"]
                        pk["wire_at"][s] = w
                        self.where[w] = (len(self.packets), s)
                self.packets.append(pk)
        # Lemma D.3 (1)-(2): every wire at the same slot at both ends, symbols match the convention
        for w, (rv, cv) in enumerate(net.wires):
            P, s = self.where[w]
            for i in range(a):
                assert self.row_sym(P, s, i) == labels[rv][i]
                assert self.col_sym(P, s, i) == labels[cv][i]
        self.Cp0 = 1
        for s in self.S:
            if s != 0:
                self.Cp0 *= -s

    def row_sym(self, P, x, i):
        lo, hi = self.packets[P]["pat"][i]
        return lo if lo == hi or not (x >> i) & 1 else hi

    def col_sym(self, P, y, i):
        lo, hi = self.packets[P]["pat"][i]
        return lo if lo == hi else (hi if not (y >> i) & 1 else lo)

    def Wpoly(self, d):
        """C'(0) W_xy(T) for y - x = d: prod_{s in S, s != d} (T - s) (integer polynomial in T)."""
        p = [1]
        for s in self.S:
            if s != d:
                p = poly_mul(p, [-s, 1])
        return p


# ----------------------------------------------------------------------------- F(T) (D.8)
def F_poly(net, pk):
    """Exact F(T) as a list of Fraction coefficients, by evaluation at deg+1 integer points."""
    M, L = len(net.wires), pk.L
    deg = M * (2 * L - 1)
    nonzero = []
    for x in range(1 << (2 * M)):
        X = [(x >> w) & 1 for w in range(M)]
        Y = [(x >> (M + w)) & 1 for w in range(M)]
        p = 1
        for v in net.verts:
            p *= v["table"].get(tuple(Y[w] if e == "r" else X[w] for w, e in v["ports"]), 0)
            if p == 0:
                break
        if p:
            nonzero.append((p, X, Y))
    Wp = {}
    pts = list(range(deg + 1))
    vals = []
    for T0 in pts:
        wv = {}
        for P, pkd in enumerate(pk.packets):
            for x in pkd["wire_at"]:
                for y in pkd["wire_at"]:
                    d = y - x
                    if d not in Wp:
                        Wp[d] = pk.Wpoly(d)
                    wv[(P, x, y)] = sum(c * T0 ** k for k, c in enumerate(Wp[d]))
        tot = 0
        for p, X, Y in nonzero:
            prod = p
            for w in range(M):
                P, x = pk.where[w]
                form = sum(wv[(P, x, y)] * X[w2] for y, w2 in pk.packets[P]["wire_at"].items())
                prod *= form if Y[w] else (pk.Cp0 - form)
                if prod == 0:
                    break
            tot += prod
        vals.append(Fr(tot, pk.Cp0 ** M))
    # Newton interpolation -> monomial coefficients
    n = len(pts)
    coef = vals[:]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (pts[i] - pts[i - j])
    poly = [coef[n - 1]]
    for i in range(n - 2, -1, -1):
        poly = poly_mul(poly, [Fr(-pts[i]), Fr(1)])
        poly[0] += coef[i]
    return poly


def taylor_at(poly, t, alpha, q):
    """Coefficients of F(t + alpha z) mod z^q (exact)."""
    out = [Fr(0)] * q
    for k, c in enumerate(poly):
        if c:
            for j in range(min(k, q - 1) + 1):
                out[j] += c * math.comb(k, j) * Fr(t) ** (k - j) * Fr(alpha) ** j
    return out


# ----------------------------------------------------------------------------- Lemma D.5
def multiplicities(a, u, h):
    L = 1 << a
    S = list(range(-L, L))
    q = {}
    for t in S:
        if t != 0:
            i = v2(t)
            q[t] = h * max(0, u - i)
    return q


def recover(data, qmult, a, exact=True, P=None):
    """Hermite recovery (D.9). data[t] = list of q_t coefficients (Fractions, or residues mod 2^P[t])."""
    L = 1 << a
    alpha = 2 ** (a + 1)
    S = list(range(-L, L))
    total = Fr(0)
    for t, qt in qmult.items():
        if qt == 0:
            continue
        # E_t(t + alpha z) = prod_{q != 0, t} (1 + alpha z / (t - q))^{q_q}, E_t(0) = prod ((0 - q)/(t - q))^{q_q}
        Et0 = Fr(1)
        ser = [Fr(1)] + [Fr(0)] * (qt - 1)
        for s in S:
            if s in (0, t) or qmult.get(s, 0) == 0:
                continue
            Et0 *= Fr(-s, t - s) ** qmult[s]
            fac = [Fr(1), Fr(alpha, t - s)] + [Fr(0)] * max(0, qt - 2)
            for _ in range(qmult[s]):
                ser = ser_mul(ser, fac[:qt], qt)
        inv = ser_inv(ser, qt)
        d = [Fr(c) for c in data[t]]
        quo = ser_mul(d, inv, qt)
        zval = Fr(-t, alpha)
        total += Et0 * sum(c * zval ** j for j, c in enumerate(quo))
    return total


def loss_bits(qmult, a, t):
    """-v2(E_t(0)) + (q_t - 1)(a + 1 - i): precision lost by the recovery at node t."""
    L = 1 << a
    Et0 = Fr(1)
    for s in range(-L, L):
        if s in (0, t) or qmult.get(s, 0) == 0:
            continue
        Et0 *= Fr(-s, t - s) ** qmult[s]
    i = v2(t)
    return -v2(Et0) + (qmult[t] - 1) * (a + 1 - i)


# ----------------------------------------------------------------------------- modes (D.4.5-D.4.6)
class NodeExpansion:
    """Mode expansion of F(t + alpha z) mod (2^p, z^q) at node t with i = v2(t)."""

    def __init__(self, net, pk, labels, t, p, q, weights="paper", basis="power"):
        assert basis in ("power", "newton", "ldu")
        assert basis != "ldu" or weights == "sharp", "the ldu basis has sharp weights only"
        self.basis = basis
        self.net, self.pk, self.t, self.p, self.q = net, pk, t, p, q
        self.m = 1 << p
        a, L = pk.a, pk.L
        self.alpha = 2 ** (a + 1)
        i = v2(t)
        self.i = i
        self.grp = [labels[v][i] for v in range(len(net.verts))]
        M = len(net.wires)
        # exact series of C'(0) W_xy(t + alpha z) / C'(0)
        self.Wser = {}
        for P, pkd in enumerate(pk.packets):
            for x in range(L):
                for y in range(L):
                    d = y - x
                    poly = pk.Wpoly(d)
                    s = taylor_at([Fr(c, pk.Cp0) for c in poly], t, self.alpha, q)
                    self.Wser[(P, x, y)] = s
        # groups of the ports
        self.rowg = {w: self.grp[net.wires[w][0]] for w in range(M)}
        self.colg = {w: self.grp[net.wires[w][1]] for w in range(M)}
        # modes
        self.modes = []          # (P, y0, b, source, target, weight)
        self.coef = {}           # (w_row, mode index) -> residue series
        self.local = {}          # w_row -> list of (w_col, residue series)
        for w in range(M):
            P, x = pk.where[w]
            self.local[w] = []
            for y, w2 in pk.packets[P]["wire_at"].items():
                if self.colg[w2] == self.rowg[w]:
                    self.local[w].append((w2, [red(c, p) for c in self.Wser[(P, x, y)]]))
        step = 1 << (i + 1)
        if basis == "ldu":
            self._ldu_modes(i, step)
        for P, pkd in enumerate(pk.packets):
            lo, hi = pkd["pat"][i]
            if lo == hi or basis == "ldu":
                continue
            for y0 in range(step):
                cols = [(j, pkd["wire_at"][y0 + step * j]) for j in range(L // step)
                        if (y0 + step * j) in pkd["wire_at"]]
                src = pk.col_sym(P, y0, i)
                tgt = hi if src == lo else lo
                rows = [(x, w) for x, w in pkd["wire_at"].items() if pk.row_sym(P, x, i) == tgt]
                if not cols or not rows:
                    continue
                cache = {}           # (x, b) -> exact coefficient series of mode b in row x
                for b in range(p):
                    extra = v2fact(b) if basis == "newton" else 0
                    wgt = (b + 1 + extra) if weights == "paper" else (a + 1 - i + b + extra)
                    if wgt >= p:
                        break
                    # L_Z(X) = sum_j j^b X_y (power) or sum_j binom(j, b) X_y (newton); trivial if all vanish
                    lin = [(self.bfun(j, b), w2) for j, w2 in cols]
                    lin = [(c, w2) for c, w2 in lin if c]
                    k = len(self.modes)
                    if lin:
                        self.modes.append(dict(P=P, y0=y0, b=b, src=src, tgt=tgt, w=wgt, lin=lin))
                    for x, wr in rows:
                        A = Fr(t - y0 + x)
                        assert v2(A) <= i
                        if b == 0:
                            c = taylor_at([Fr(c, pk.Cp0) for c in pk.Wpoly(y0 - x)], t, self.alpha, q)
                        else:
                            prev = cache[(x, b - 1)]
                            if basis == "power":       # ratio c / (A + alpha z)
                                base, num = [A, Fr(self.alpha)], Fr(step)
                            else:                      # ratio c b / (A - c b + alpha z)
                                base, num = [A - step * b, Fr(self.alpha)], Fr(step * b)
                            base = (base + [Fr(0)] * q)[:q]
                            ratio = [num * cc for cc in ser_inv(base, q)]
                            c = ser_mul(prev, ratio, q)
                        cache[(x, b)] = c
                        assert all(v2(cc) >= wgt for cc in c), ("weight", basis, b, wgt, [v2(cc) for cc in c])
                        if lin:
                            self.coef[(wr, k)] = [red(cc, p) for cc in c]
        # (D.12) check: every crossing entry equals its mode expansion mod 2^p
        for P, pkd in enumerate(pk.packets):
            for x, wr in pkd["wire_at"].items():
                for y, wc in pkd["wire_at"].items():
                    if self.colg[wc] == self.rowg[wr]:
                        continue
                    y0, j = y % step, y // step
                    acc = [0] * q
                    for k, md in enumerate(self.modes):
                        if basis == "ldu":
                            if md["P"] == P and (wr, k) in self.coef:
                                cy = dict((w2, c) for c, w2 in md["lin"]).get(wc)
                                if cy is not None:
                                    acc = ps_add_mod(acc, ps_mul_mod(list(cy), self.coef[(wr, k)], q, self.m), self.m)
                        elif md["P"] == P and md["y0"] == y0 and (wr, k) in self.coef:
                            jb = self.bfun(j, md["b"])
                            acc = ps_add_mod(acc, [(jb * c) % self.m for c in self.coef[(wr, k)]], self.m)
                    direct = [red(c, p) for c in self.Wser[(P, x, y)]]
                    assert acc == direct, ("D.12 failed", P, x, y)
        self.r_groups = sorted(set(self.grp))
        self.tg = {g: [k for k, md in enumerate(self.modes) if md["tgt"] == g] for g in self.r_groups}
        self.sg = {g: [k for k, md in enumerate(self.modes) if md["src"] == g] for g in self.r_groups}
        self.inc = {g: sorted(set(self.tg[g]) | set(self.sg[g])) for g in self.r_groups}
        self.child_calls = 0

    def _ldu_modes(self, i, step):
        """Cauchy LDU modes (see the module docstring): one block per packet half h (slots with bit i = h)."""
        pk, p, q, t, m = self.pk, self.p, self.q, self.t, self.m
        a, L, alpha = pk.a, pk.L, self.alpha
        half = 1 << i
        n = L // 2
        Cpoly = [1]
        for s in pk.S:
            Cpoly = poly_mul(Cpoly, [-s, 1])
        G = taylor_at([Fr(c, pk.Cp0) for c in Cpoly], t, alpha, q)
        self.ldu_blocks = []
        for P, pkd in enumerate(pk.packets):
            lo, hi = pkd["pat"][i]
            if lo == hi:
                continue
            for h in (0, 1):
                pos = [step * (k >> i) + (k & (half - 1)) + h * half for k in range(n)]
                src = pk.col_sym(P, pos[0], i)
                tgt = hi if src == lo else lo
                occ = [(k, pkd["wire_at"][pos[k]]) for k in range(n) if pos[k] in pkd["wire_at"]]
                rows = [(k, w) for k, w in occ if pk.row_sym(P, pos[k], i) == tgt]
                if not occ or not rows:
                    continue
                K = [[ser_inv(([Fr(t + pos[r] - pos[s]), Fr(alpha)] + [Fr(0)] * q)[:q], q) for s in range(n)]
                     for r in range(n)]
                Lm, d, Um = series_ldu(K, q)
                self.ldu_blocks.append(dict(P=P, h=h, n=n, d=d))
                for mm in range(n):
                    b = mm >> i
                    wgt = a + 1 - i + 2 * b + 2 * v2fact(b)
                    if wgt >= p:
                        break
                    lin = []
                    for k, w2 in occ:
                        if k >= mm:
                            cser = tuple(red(c, p) for c in Um[mm][k])
                            assert all(v2(c) >= 0 for c in Um[mm][k]), "U not 2-adically integral"
                            if any(cser):
                                lin.append((cser, w2))
                    if not lin:
                        continue
                    kk = len(self.modes)
                    self.modes.append(dict(P=P, y0=h, b=mm, src=src, tgt=tgt, w=wgt, lin=lin))
                    gd = ser_mul(G, d[mm], q)
                    for k, wr in rows:
                        if k < mm:
                            continue
                        assert all(v2(c) >= 0 for c in Lm[k][mm]), "L not 2-adically integral"
                        c = ser_mul(gd, Lm[k][mm], q)
                        assert all(v2(cc) >= wgt for cc in c), ("ldu weight", mm, wgt, [v2(cc) for cc in c])
                        self.coef[(wr, kk)] = [red(cc, p) for cc in c]

    def bfun(self, j, b):
        """Value at column index j of the b-th source form."""
        if self.basis == "newton":
            return math.comb(j, b)
        return j ** b if (j or b == 0) else 0

    # -- tuples
    def tuples(self, idx, budget):
        """All exponent tuples (dict mode -> exp > 0) on modes idx with total weight < budget."""
        out = []

        def rec(j, left, cur):
            if j == len(idx):
                out.append(dict(cur))
                return
            k = idx[j]
            wk = self.modes[k]["w"]
            e = 0
            while e * wk < left:
                if e:
                    cur[k] = e
                rec(j + 1, left - e * wk, cur)
                e += 1
            cur.pop(k, None)
        rec(0, budget, {})
        return out

    def weight(self, beta):
        return sum(self.modes[k]["w"] * e for k, e in beta.items())

    # -- children
    def group_table(self, g):
        """V_g[gamma] for all tuples gamma on I_g of weight < p (each child contracted once)."""
        net, m, q, p = self.net, self.m, self.q, self.p
        verts = [v for v in range(len(net.verts)) if self.grp[v] == g]
        rows = [w for w in range(len(net.wires)) if self.rowg[w] == g]
        cols = [w for w in range(len(net.wires)) if self.colg[w] == g]
        ports = sorted(set(rows) | set(cols))
        gammas = self.tuples(self.inc[g], p)
        tg = set(self.tg[g])
        table = {self.key(gm): [0] * q for gm in gammas}
        for bits in itertools.product((0, 1), repeat=len(rows) + len(cols)):
            Y = dict(zip(rows, bits[:len(rows)]))
            X = dict(zip(cols, bits[len(rows):]))
            tv = 1
            for v in verts:
                tv *= net.verts[v]["table"].get(
                    tuple(Y[w] if e == "r" else X[w] for w, e in net.verts[v]["ports"]), 0)
                if tv == 0:
                    break
            if tv == 0:
                continue
            # product of the row factors as a polynomial in the target modes
            polyd = {(): [tv % m] + [0] * (q - 1)}
            for w in rows:
                loc = [0] * q
                for w2, ser in self.local[w]:
                    if X[w2]:
                        loc = ps_add_mod(loc, ser, m)
                terms = [((), loc if Y[w] else [(1 - loc[0]) % m] + [(-c) % m for c in loc[1:]])]
                for k in self.tg[g]:
                    if (w, k) in self.coef:
                        c = self.coef[(w, k)]
                        terms.append((k, c if Y[w] else [(-cc) % m for cc in c]))
                new = {}
                for key, val in polyd.items():
                    for k, c in terms:
                        if k == ():
                            nk = key
                        else:
                            dd = dict(key)
                            dd[k] = dd.get(k, 0) + 1
                            if self.weight(dd) >= p:
                                continue
                            nk = tuple(sorted(dd.items()))
                        prod = ps_mul_mod(val, c, q, m)
                        if nk in new:
                            new[nk] = ps_add_mod(new[nk], prod, m)
                        else:
                            new[nk] = prod
                polyd = new
            Lval = {k: lin_eval(self.modes[k]["lin"], X.__getitem__, q, m) for k in self.sg[g]}
            for gm in gammas:
                tkey = tuple(sorted((k, e) for k, e in gm.items() if k in tg))
                if tkey not in polyd:
                    continue
                f = 1
                for k, e in gm.items():
                    if k not in tg:
                        f = src_mul(f, Lval[k], e, q, m)
                if src_zero(f):
                    continue
                kk = self.key(gm)
                table[kk] = ps_add_mod(table[kk], src_scale(f, polyd[tkey], q, m), m)
        self.child_calls += len(gammas)
        return table

    @staticmethod
    def key(beta):
        return tuple(sorted(beta.items()))

    def sum_paper(self):
        """(D.16) as in the paper: every tuple beta, every group (children recomputed per tuple)."""
        tables = {g: self.group_table(g) for g in self.r_groups}
        tot = [0] * self.q
        nb = 0
        for beta in self.tuples(list(range(len(self.modes))), self.p):
            nb += 1
            val = [1] + [0] * (self.q - 1)
            for g in self.r_groups:
                sub = {k: e for k, e in beta.items() if k in self.inc[g]}
                val = ps_mul_mod(val, tables[g][self.key(sub)], self.q, self.m)
            tot = ps_add_mod(tot, val, self.m)
        return tot, nb, nb * len(self.r_groups)

    def sum_memo(self):
        self.child_calls = 0
        tables = {g: self.group_table(g) for g in self.r_groups}
        tot = [0] * self.q
        nb = 0
        for beta in self.tuples(list(range(len(self.modes))), self.p):
            nb += 1
            val = [1] + [0] * (self.q - 1)
            for g in self.r_groups:
                sub = {k: e for k, e in beta.items() if k in self.inc[g]}
                val = ps_mul_mod(val, tables[g][self.key(sub)], self.q, self.m)
            tot = ps_add_mod(tot, val, self.m)
        return tot, nb, self.child_calls

    def sum_dp(self):
        """Group-by-group elimination: states are tuples on the cut, steps enumerate cut + new modes."""
        self.child_calls = 0
        tables = {g: self.group_table(g) for g in self.r_groups}
        order = self.r_groups
        pos = {g: n for n, g in enumerate(order)}
        ends = [(md["src"], md["tgt"]) for md in self.modes]
        Phi = {(): [1] + [0] * (self.q - 1)}
        max_enum = 0
        for n_, g in enumerate(order):
            N_new = [k for k in self.inc[g] if max(pos[ends[k][0]], pos[ends[k][1]]) > n_
                     and min(pos[ends[k][0]], pos[ends[k][1]]) == n_]
            new_tuples = self.tuples(N_new, self.p)
            newPhi = {}
            enum = 0
            for skey, sval in Phi.items():
                s = dict(skey)
                ws = self.weight(s)
                for nt in new_tuples:
                    if ws + self.weight(nt) >= self.p:
                        continue
                    enum += 1
                    full = dict(s)
                    full.update(nt)
                    gamma = {k: e for k, e in full.items() if k in self.inc[g]}
                    val = ps_mul_mod(sval, tables[g][self.key(gamma)], self.q, self.m)
                    # project: drop modes whose both ends are now processed
                    proj = {k: e for k, e in full.items()
                            if max(pos[ends[k][0]], pos[ends[k][1]]) > n_}
                    pk_ = self.key(proj)
                    newPhi[pk_] = ps_add_mod(newPhi[pk_], val, self.m) if pk_ in newPhi else val
            Phi = newPhi
            max_enum = max(max_enum, enum)
        return Phi.get((), [0] * self.q), max_enum, self.child_calls

    def child_network_value(self, g, gamma):
        """Build the child network of D.4.6 explicitly and contract it by brute force."""
        net, pk, m, q = self.net, self.pk, self.m, self.q
        verts = [v for v in range(len(net.verts)) if self.grp[v] == g]
        child_ports = sorted({w for w in range(len(net.wires)) if self.rowg[w] == g} |
                             {w for w in range(len(net.wires)) if self.colg[w] == g})
        # edges of the child: (wire, end) for every retained port; each joins its old vertex and its packet
        edges = [(w, "r") for w in child_ports if self.rowg[w] == g] + \
                [(w, "c") for w in child_ports if self.colg[w] == g]
        tg = set(self.tg[g])
        tot = [0] * q
        for bits in itertools.product((0, 1), repeat=len(edges)):
            val = dict(zip(edges, bits))
            tv = 1
            for v in verts:
                tv *= net.verts[v]["table"].get(
                    tuple(val[(w, e)] for w, e in net.verts[v]["ports"]), 0)
                if tv == 0:
                    break
            if tv == 0:
                continue
            ser = [tv % m] + [0] * (q - 1)
            for P in range(len(pk.packets)):
                ser = ps_mul_mod(ser, self.packet_tensor(g, P, gamma, val), q, m)
            tot = ps_add_mod(tot, ser, m)
        return tot

    def packet_tensor(self, g, P, gamma, val):
        m, q = self.m, self.q
        pkd = self.pk.packets[P]
        rows = [w for w in pkd["wire_at"].values() if self.rowg[w] == g]
        tmodes = [k for k in self.tg[g] if self.modes[k]["P"] == P]
        smodes = [k for k in self.sg[g] if self.modes[k]["P"] == P]
        want = tuple(sorted((k, gamma[k]) for k in tmodes if k in gamma))
        polyd = {(): [1] + [0] * (q - 1)}
        for w in rows:
            loc = [0] * q
            for w2, ser in self.local[w]:
                if val[(w2, "c")]:
                    loc = ps_add_mod(loc, ser, m)
            yb = val[(w, "r")]
            terms = [((), loc if yb else [(1 - loc[0]) % m] + [(-c) % m for c in loc[1:]])]
            for k in tmodes:
                if (w, k) in self.coef:
                    c = self.coef[(w, k)]
                    terms.append((k, c if yb else [(-cc) % m for cc in c]))
            new = {}
            for key, v in polyd.items():
                for k, c in terms:
                    if k == ():
                        nk = key
                    else:
                        dd = dict(key)
                        dd[k] = dd.get(k, 0) + 1
                        if dd[k] > gamma.get(k, 0):
                            continue
                        nk = tuple(sorted(dd.items()))
                    prod = ps_mul_mod(v, c, q, m)
                    new[nk] = ps_add_mod(new[nk], prod, m) if nk in new else prod
            polyd = new
        ser = polyd.get(want, [0] * q)
        f = 1
        for k in smodes:
            e = gamma.get(k, 0)
            if e:
                f = src_mul(f, lin_eval(self.modes[k]["lin"], lambda w2: val[(w2, "c")], q, m), e, q, m)
        return src_scale(f, ser, q, m)


# ----------------------------------------------------------------------------- tests
def random_formula(n, m, rng, width=3, maxocc=3):
    for _ in range(1000):
        occ = {v: 0 for v in range(1, n + 1)}
        cls = []
        ok = True
        for _ in range(m):
            cand = [v for v in range(1, n + 1) if occ[v] < maxocc]
            if len(cand) < 1:
                ok = False
                break
            k = min(width, len(cand), rng.randint(2, width))
            vs = rng.sample(cand, k)
            for v in vs:
                occ[v] += 1
            cls.append(tuple(v if rng.random() < 0.5 else -v for v in vs))
        if ok:
            return cls
    raise RuntimeError


def check_networks(rng, trials=30):
    for _ in range(trials):
        n = rng.randint(3, 7)
        m = rng.randint(2, 5)
        cls = random_formula(n, m, rng)
        c = count_sat(n, cls)
        assert formula_net(n, cls).contract() == c
        units = make_units(n, cls, rng.randint(1, 3), rng)
        assert chunk_net_mapped(n, cls, units).contract() == c
    return True


def setup(rng, n, m, a, r, chunks=False, Lc=2):
    cls = random_formula(n, m, rng)
    if chunks:
        net = chunk_net_mapped(n, cls, make_units(n, cls, Lc, rng))
    else:
        net = formula_net(n, cls)
    labels = [tuple(rng.randrange(r) for _ in range(a)) for _ in net.verts]
    return n, cls, net, labels, Packing(net, labels, a)


def direct_taylor_mod(net, pk, t, p, q):
    """F(t + alpha z) mod (2^p, z^q) straight from (D.8), without modes."""
    m, alpha, M = 1 << p, 2 ** (pk.a + 1), len(net.wires)
    Wres = {}
    for P, pkd in enumerate(pk.packets):
        for x in pkd["wire_at"]:
            for y in pkd["wire_at"]:
                Wres[(P, x, y)] = [red(c, p) for c in
                                   taylor_at([Fr(c, pk.Cp0) for c in pk.Wpoly(y - x)], t, alpha, q)]
    tot = [0] * q
    for xx in range(1 << (2 * M)):
        X = [(xx >> w) & 1 for w in range(M)]
        Y = [(xx >> (M + w)) & 1 for w in range(M)]
        tv = 1
        for v in net.verts:
            tv *= v["table"].get(tuple(Y[w] if e == "r" else X[w] for w, e in v["ports"]), 0)
            if tv == 0:
                break
        if tv == 0:
            continue
        ser = [tv % m] + [0] * (q - 1)
        for w in range(M):
            P, x = pk.where[w]
            form = [0] * q
            for y, w2 in pk.packets[P]["wire_at"].items():
                if X[w2]:
                    form = ps_add_mod(form, Wres[(P, x, y)], m)
            f = form if Y[w] else [(1 - form[0]) % m] + [(-c) % m for c in form[1:]]
            ser = ps_mul_mod(ser, f, q, m)
        tot = ps_add_mod(tot, ser, m)
    return tot


def byvars_value(ne, g, gamma):
    """Leaf child of group g for a chunk network, summed over the variables of its chunks (Gemini, item 2):
    fixing the variables fixes every column wire of the group; each designated clause then depends only on
    its external row bits and splits into at most two product terms; the row bits enter the packet factors
    as unary weights. Must equal V_g[gamma]."""
    net, pk, m, q = ne.net, ne.pk, ne.m, ne.q
    units_g = [k for k in range(len(net.verts)) if ne.grp[k] == g]
    vars_g = sorted(set().union(*[net.units[k] for k in units_g])) if units_g else []
    cls_g = [ci for ci in range(len(net.clauses)) if net.desig[ci] in units_g]
    rows_g = [w for w in range(len(net.wires)) if ne.rowg[w] == g]
    tg = set(ne.tg[g])
    total = [0] * q
    branch_max = 0
    for xb in itertools.product((0, 1), repeat=len(vars_g)):
        val = dict(zip(vars_g, xb))
        X = {w: val[_WIREMAP[(id(net), w)][1]] for w in range(len(net.wires)) if ne.colg[w] == g}
        # term lists per clause: (coefficient, {row wire: (weight at y=0, weight at y=1)})
        termlists = []
        dead = False
        for ci in cls_g:
            U = net.units[net.desig[ci]]
            internal_true = False
            ext = []
            for l in net.clauses[ci]:
                v = abs(l)
                if v in U:
                    internal_true |= (val[v] == (1 if l > 0 else 0))
                else:
                    ext.append((_WIREMAP[(id(net), "cv", ci, v)], l > 0))
            if internal_true:
                continue
            if not ext:
                dead = True
                break
            # sat = 1 - prod_j [literal j false]: two product terms
            falsew = {}
            for w, pos in ext:
                falsew[w] = (1, 0) if pos else (0, 1)      # weight of y when literal j is false
            termlists.append([(1, {}), (-1, falsew)])
        if dead:
            continue
        branch_max = max(branch_max, sum(1 for tl in termlists if len(tl) > 1))
        for choice in itertools.product(*termlists) if termlists else [()]:
            coef = 1
            wts = {}
            for c, wd in choice:
                coef *= c
                for w, wv in wd.items():
                    old = wts.get(w, (1, 1))
                    wts[w] = (old[0] * wv[0], old[1] * wv[1])
            ser = [coef % m] + [0] * (q - 1)
            for P in range(len(pk.packets)):
                pkd = pk.packets[P]
                prow = [w for w in pkd["wire_at"].values() if ne.rowg[w] == g]
                tmodes = [k for k in ne.tg[g] if ne.modes[k]["P"] == P]
                smodes = [k for k in ne.sg[g] if ne.modes[k]["P"] == P]
                want = tuple(sorted((k, gamma[k]) for k in tmodes if k in gamma))
                polyd = {(): [1] + [0] * (q - 1)}
                for w in prow:
                    loc = [0] * q
                    for w2, sr in ne.local[w]:
                        if X[w2]:
                            loc = ps_add_mod(loc, sr, m)
                    w0, w1 = wts.get(w, (1, 1))
                    # w0 (1 - form) + w1 form = w0 + (w1 - w0) form
                    base = [(w0 + (w1 - w0) * loc[0]) % m] + [((w1 - w0) * c) % m for c in loc[1:]]
                    terms = [((), base)]
                    for k in tmodes:
                        if (w, k) in ne.coef:
                            terms.append((k, [((w1 - w0) * c) % m for c in ne.coef[(w, k)]]))
                    new = {}
                    for key, v in polyd.items():
                        for k, c in terms:
                            if k == ():
                                nk = key
                            else:
                                dd = dict(key)
                                dd[k] = dd.get(k, 0) + 1
                                if dd[k] > gamma.get(k, 0):
                                    continue
                                nk = tuple(sorted(dd.items()))
                            pr = ps_mul_mod(v, c, q, m)
                            new[nk] = ps_add_mod(new[nk], pr, m) if nk in new else pr
                    polyd = new
                pv = polyd.get(want, [0] * q)
                f = 1
                for k in smodes:
                    e = gamma.get(k, 0)
                    if e:
                        f = src_mul(f, lin_eval(ne.modes[k]["lin"], X.__getitem__, q, m), e, q, m)
                ser = ps_mul_mod(ser, src_scale(f, pv, q, m), q, m)
            total = ps_add_mod(total, ser, m)
    return total, len(vars_g), branch_max


BYVARS = []


def stress(rng, trials):
    """Mode identities on larger toys: direct (D.8) = per tuple = memo = group-by-group DP."""
    out = []
    for trial in range(trials):
        a = rng.choice([1, 2, 2, 3])
        r = rng.choice([3, 4, 5])
        chunks = trial % 3 == 2
        n_, cls, net, labels, pk = setup(rng, rng.randint(4, 8 if chunks else 5), rng.randint(3, 5), a, r,
                                         chunks=chunks, Lc=rng.choice([2, 3]))
        M = len(net.wires)
        if M < 3 or M > 8:
            continue
        t = rng.choice([s for s in pk.S if s and v2(s) < a])
        p, q = rng.choice([4, 5, 6]), rng.choice([1, 2, 3])
        wts = rng.choice(["paper", "sharp"])
        ne = NodeExpansion(net, pk, labels, t, p, q, weights=wts)
        direct = direct_taylor_mod(net, pk, t, p, q)
        s1, nb, _ = ne.sum_paper()
        s2, _, calls2 = ne.sum_memo()
        s3, enum, _ = ne.sum_dp()
        assert s1 == direct and s2 == direct and s3 == direct, ("stress D.16", trial)
        if chunks:
            for g in ne.r_groups:
                tab = ne.group_table(g)
                for gm in rng.sample(list(tab.keys()), min(4, len(tab))):
                    bv, nv, br = byvars_value(ne, g, dict(gm))
                    assert bv == tab[gm], ("byvars leaf", trial, g)
                    BYVARS.append((len([w for w in range(M) if ne.rowg[w] == g or ne.colg[w] == g]), nv, br))
        out.append((a, r, chunks, M, t, p, q, wts, len(ne.modes), nb, nb * len(ne.r_groups), calls2, enum))
    return out


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    report = {}
    check_networks(rng)
    report["networks"] = "plain and chunk networks contract to #F on 30 random formulas"

    # 3 + 4: F(T), Lemma D.4, tapered Hermite recovery with truncated Taylor data
    rec = []
    for trial in range(8):
        a = rng.choice([2, 3])
        n_, cls, net, labels, pk = setup(rng, rng.randint(3, 4), rng.randint(2, 3), a, 3,
                                     chunks=(trial % 2 == 1))
        M = len(net.wires)
        if M == 0 or M > 6:
            continue
        F = F_poly(net, pk)
        deg = max(k for k, c in enumerate(F) if c) if any(F) else 0
        assert deg <= M * (2 * pk.L - 1)
        assert F[0] == net.contract_core()
        alpha = 2 ** (a + 1)
        for t in pk.S:
            if t:
                assert all(v2(c) >= 0 for c in taylor_at(F, t, alpha, 6))      # Lemma D.4
        for u in range(1, a + 1):
            Lh = pk.L
            hden = sum((Lh >> i) * (u - i) for i in range(u))
            h = (M * (2 * Lh - 1)) // hden + 1
            qm = multiplicities(a, u, h)
            assert sum(qm.values()) > M * (2 * Lh - 1)
            data = {t: taylor_at(F, t, alpha, qt) for t, qt in qm.items() if qt}
            assert recover(data, qm, a) == F[0]
            p = 6
            Pt = {t: p + loss_bits(qm, a, t) for t in data}
            trunc = {t: [red(c, Pt[t]) for c in data[t]] for t in data}
            got = recover(trunc, qm, a)
            assert v2(got - F[0]) >= p, ("precision", v2(got - F[0]))
            # with 3 bits less at every node the recovery may fail: record how often
            trunc2 = {t: [red(c, max(1, Pt[t] - 3)) for c in data[t]] for t in data}
            ok2 = v2(recover(trunc2, qm, a) - F[0]) >= p
            bound = {t: p + 6 * h * (Lh >> v2(t)) + 4 * M * (a + 1) for t in data}
            rec.append(dict(a=a, u=u, M=M, h=h, Pt_max=max(Pt.values()),
                            bound_max=max(bound.values()), ok_minus3=ok2))
    report["hermite"] = rec

    # 5: mode expansion identities at small precision
    modes_res = []
    for trial in range(10):
        a = rng.choice([1, 2, 2, 3])
        r = rng.choice([3, 4])
        chunks = trial % 3 == 2
        n_, cls, net, labels, pk = setup(rng, rng.randint(3, 5), rng.randint(2, 4), a, r, chunks=chunks)
        M = len(net.wires)
        if M == 0 or M > 6:
            continue
        F = F_poly(net, pk)
        alpha = 2 ** (a + 1)
        for t in rng.sample([s for s in pk.S if s and v2(s) < a], 2):
            p, q = rng.choice([4, 5]), rng.choice([1, 2, 3])
            for wts in ("paper", "sharp"):
                ne = NodeExpansion(net, pk, labels, t, p, q, weights=wts)
                direct = [red(c, p) for c in taylor_at(F, t, alpha, q)]
                s1, nb, calls1 = ne.sum_paper()
                s2, nb2, calls2 = ne.sum_memo()
                s3, enum, calls3 = ne.sum_dp()
                assert s1 == direct and s2 == direct and s3 == direct, ("D.16", t, wts, direct, s1, s2, s3)
                # child networks of D.4.6 contract to the table values
                for g in ne.r_groups:
                    tab = ne.group_table(g)
                    for gm in rng.sample(list(tab.keys()), min(3, len(tab))):
                        assert ne.child_network_value(g, dict(gm)) == tab[gm]
                modes_res.append(dict(a=a, r=r, chunks=chunks, M=M, t=t, p=p, q=q, w=wts,
                                      modes=len(ne.modes), tuples=nb, paper_children=nb * len(ne.r_groups),
                                      memo_children=calls2, dp_max_enum=enum))
    report["modes"] = modes_res

    # 6: end to end at a = 1 (u = 1): #F from mode sums at the precision of Lemma D.5
    e2e = []
    for trial in range(6):
        chunks = trial % 2 == 1
        r = 3
        n_, cls, net, labels, pk = setup(rng, rng.randint(3, 4), rng.randint(2, 3), 1, r, chunks=chunks)
        M = len(net.wires)
        if M == 0 or M > 5:
            continue
        nvars = n_
        target = count_sat(nvars, cls)
        p = nvars + 1
        a, u = 1, 1
        Lh = 2
        h = (M * (2 * Lh - 1)) // Lh + 1
        qm = multiplicities(a, u, h)
        data = {}
        stats = []
        for t, qt in qm.items():
            if not qt:
                continue
            Pt = p + loss_bits(qm, a, t)
            ne = NodeExpansion(net, pk, labels, t, Pt, qt, weights="sharp")
            s3, enum, calls = ne.sum_dp()
            data[t] = s3
            stats.append((t, Pt, qt, len(ne.modes), enum, calls))
        got = recover(data, qm, a)
        assert v2(got - Fr(net.contract_core())) >= p, "end to end"
        cnt = red(got, p) << net.unused
        assert cnt == target, (cnt, target)
        e2e.append(dict(chunks=chunks, M=M, count=target, nodes=stats))
    report["end_to_end"] = e2e
    st = stress(rng, int(sys.argv[2]) if len(sys.argv) > 2 else 40)
    report["stress"] = dict(cases=len(st), max_modes=max(x[8] for x in st),
                            max_tuples=max(x[9] for x in st),
                            dp_below_paper=sum(1 for x in st if x[12] < x[9]),
                            memo_below_paper=sum(1 for x in st if x[11] < x[10]),
                            byvars_checks=len(BYVARS), byvars_rows=BYVARS[:20],
                            rows=st)
    return report


if __name__ == "__main__":
    import json
    rep = main()
    print(json.dumps(rep, indent=1, default=str))


# ---- round 7 (appd-shellmm): the group-by-group mode sum as shell matrix products with Strassen's algorithm


def ps_sub_mod(a, b, m):
    return [(x - y) % m for x, y in zip(a, b)]


def naive_matmul(A, B, q, m):
    """Schoolbook product of an n x k and a k x l matrix whose entries are series mod (m, z^q)."""
    n, k, l = len(A), len(B), len(B[0])
    out = [[[0] * q for _ in range(l)] for _ in range(n)]
    for i in range(n):
        for kk in range(k):
            a = A[i][kk]
            for j in range(l):
                out[i][j] = ps_add_mod(out[i][j], ps_mul_mod(a, B[kk][j], q, m), m)
    return out


def _mm_add(X, Y, m):
    return [[ps_add_mod(a, b, m) for a, b in zip(r1, r2)] for r1, r2 in zip(X, Y)]


def _mm_sub(X, Y, m):
    return [[ps_sub_mod(a, b, m) for a, b in zip(r1, r2)] for r1, r2 in zip(X, Y)]


def strassen_sq(X, Y, q, m):
    """S x S times S x S, S a power of two: Strassen's seven half-size products down to 2 x 2 blocks (schoolbook
    there), so 7^(k-1) * 8 ring products for S = 2^k >= 2."""
    S = len(X)
    if S <= 2:
        return naive_matmul(X, Y, q, m)
    h = S // 2

    def blk(M_, i, j):
        return [row[j * h:(j + 1) * h] for row in M_[i * h:(i + 1) * h]]
    A11, A12, A21, A22 = blk(X, 0, 0), blk(X, 0, 1), blk(X, 1, 0), blk(X, 1, 1)
    B11, B12, B21, B22 = blk(Y, 0, 0), blk(Y, 0, 1), blk(Y, 1, 0), blk(Y, 1, 1)
    M1 = strassen_sq(_mm_add(A11, A22, m), _mm_add(B11, B22, m), q, m)
    M2 = strassen_sq(_mm_add(A21, A22, m), B11, q, m)
    M3 = strassen_sq(A11, _mm_sub(B12, B22, m), q, m)
    M4 = strassen_sq(A22, _mm_sub(B21, B11, m), q, m)
    M5 = strassen_sq(_mm_add(A11, A12, m), B22, q, m)
    M6 = strassen_sq(_mm_sub(A21, A11, m), _mm_add(B11, B12, m), q, m)
    M7 = strassen_sq(_mm_sub(A12, A22, m), _mm_add(B21, B22, m), q, m)
    C11 = _mm_add(_mm_sub(_mm_add(M1, M4, m), M5, m), M7, m)
    C12 = _mm_add(M3, M5, m)
    C21 = _mm_add(M2, M4, m)
    C22 = _mm_add(_mm_add(_mm_sub(M1, M2, m), M3, m), M6, m)
    return [r1 + r2 for r1, r2 in zip(C11, C12)] + [r1 + r2 for r1, r2 in zip(C21, C22)]


def matmul(A, B, q, m):
    """n x k times k x l over (Z/m)[z]/(z^q): square blocks of side S = the largest power of two <= min(n, k, l)
    (zero series outside the matrices), Strassen on every block product; about n k l S^(log2(7) - 3) ring
    products."""
    n, k, l = len(A), len(B), len(B[0])
    s = min(n, k, l)
    S = 1
    while 2 * S <= s:
        S *= 2
    zero = [0] * q

    def get(M_, i, j, R_, C_):
        return M_[i][j] if i < R_ and j < C_ else zero
    out = [[zero for _ in range(l)] for _ in range(n)]
    for bi in range(-(-n // S)):
        for bj in range(-(-l // S)):
            acc = None
            for bk in range(-(-k // S)):
                X = [[get(A, bi * S + i, bk * S + j, n, k) for j in range(S)] for i in range(S)]
                Y = [[get(B, bk * S + i, bj * S + j, k, l) for j in range(S)] for i in range(S)]
                Pm = strassen_sq(X, Y, q, m)
                acc = Pm if acc is None else _mm_add(acc, Pm, m)
            for i in range(S):
                for j in range(S):
                    if bi * S + i < n and bj * S + j < l:
                        out[bi * S + i][bj * S + j] = acc[i][j]
    return out


def _mm_steps(ne):
    """Per step of the order ne.r_groups: (g, X, Y, Z) with X the cut modes not touching g, Y the cut modes
    touching g, Z the modes between g and later groups."""
    pos = {g: j for j, g in enumerate(ne.r_groups)}
    ends = [(min(pos[md["src"]], pos[md["tgt"]]), max(pos[md["src"]], pos[md["tgt"]])) for md in ne.modes]
    out = []
    for j, g in enumerate(ne.r_groups):
        cut = [k for k, (lo, hi) in enumerate(ends) if lo < j <= hi]
        inc = set(ne.inc[g])
        out.append((g, [k for k in cut if k not in inc], [k for k in cut if k in inc],
                    [k for k, (lo, hi) in enumerate(ends) if lo == j < hi]))
    return out


def _mm_shells(ne, idx):
    sh = {}
    for tpl in ne.tuples(idx, ne.p):
        tpl = dict(tpl)
        sh.setdefault(ne.weight(tpl), []).append(tpl)
    return dict(sorted(sh.items()))


def _mm_shapes(self):
    """Block shapes of sum_mm: per step, the sorted list of (|X_w1|, |Y_w2|, |Z_w3|) over w1 + w2 + w3 < p."""
    res = []
    for g, X, Y, Z in _mm_steps(self):
        xs, ys, zs = _mm_shells(self, X), _mm_shells(self, Y), _mm_shells(self, Z)
        res.append(sorted((len(a), len(b), len(c)) for w1, a in xs.items() for w2, b in ys.items()
                          for w3, c in zs.items() if w1 + w2 + w3 < self.p))
    return res


def _sum_mm(self):
    """appd-shellmm: the sum of sum_dp with every step done as shell matrix products through the module-level
    matmul. Return (series, shapes, child_calls)."""
    self.child_calls = 0
    tables = {g: self.group_table(g) for g in self.r_groups}
    q, m, p = self.q, self.m, self.p
    zero = [0] * q
    Phi = {(): [1] + [0] * (q - 1)}
    shapes = []
    for g, X, Y, Z in _mm_steps(self):
        xs, ys, zs = _mm_shells(self, X), _mm_shells(self, Y), _mm_shells(self, Z)
        new, st = {}, []
        for w1, xl in xs.items():
            for w2, yl in ys.items():
                if w1 + w2 >= p:
                    continue
                Am = [[Phi.get(self.key({**x, **y}), zero) for y in yl] for x in xl]
                for w3, zl in zs.items():
                    if w1 + w2 + w3 >= p:
                        continue
                    Bm = [[tables[g][self.key({**y, **z})] for z in zl] for y in yl]
                    Cm = matmul(Am, Bm, q, m)
                    st.append((len(xl), len(yl), len(zl)))
                    for i, x in enumerate(xl):
                        for j, z in enumerate(zl):
                            kk = self.key({**x, **z})
                            new[kk] = ps_add_mod(new[kk], Cm[i][j], m) if kk in new else Cm[i][j]
        shapes.append(sorted(st))
        Phi = new
    return Phi.get((), [0] * q), shapes, self.child_calls


NodeExpansion.mm_shapes = _mm_shapes
NodeExpansion.sum_mm = _sum_mm


# ---- round 8 (appd-guarded): Astra's guarded Lagrange deformation
# Virtual packet grid {0..Q-1}, Q = 2^a; real slots only in [H, Q - H), H = 2^(a - kappa); labels on the evaluated
# levels 0..u-1 (u <= a - kappa); W~_xy(T) = l_y(x + T), the Lagrange basis polynomial of y on the grid evaluated at
# x + T (degree Q - 1 instead of 2Q - 1); nodes t in [-H, H) with v2(t) < u; local coordinate T = t + alpha z with
# alpha = 2Q as before. At a node the crossing block of a packet half is diag(row factors) * Cauchy * diag(odd units):
# W~_xy(t + alpha z) = [2^-v P_Q(x + t + alpha z)] * 1/(t + alpha z + x - y) * [2^v / P_Q'(y)], v = Q - 1 - a.


class GuardedPacking:
    """Packing with guard slots: same interface as Packing (a, L, S, packets, where, row_sym, col_sym) plus
    guarded = True, u, kappa, H and the deformation Wfrac(x, y), rowfactor(x, t, alpha, q), colunit(y)."""
    guarded = True

    def __init__(self, net, labels, a, u, kappa):
        assert kappa >= 2 and 1 <= u <= a - kappa, "need kappa >= 2 and 1 <= u <= a - kappa"
        self.a, self.L, self.u, self.kappa = a, 1 << a, u, kappa
        Q = self.L
        H = self.H = 1 << (a - kappa)
        self.S = list(range(-H, H))
        highs = list(range(H >> u, (Q - H) >> u))
        by = {}
        for w, (rv, cv) in enumerate(net.wires):
            rl, cl = labels[rv], labels[cv]
            pat = tuple((min(rl[i], cl[i]), max(rl[i], cl[i])) for i in range(u))
            ori = tuple(0 if rl[i] <= cl[i] else 1 for i in range(u) if rl[i] != cl[i])
            by.setdefault(pat, {}).setdefault(ori, []).append(w)
        self.packets, self.where = [], {}
        for pat, d in sorted(by.items()):
            uneq = [i for i in range(u) if pat[i][0] != pat[i][1]]
            eq = [i for i in range(u) if pat[i][0] == pat[i][1]]
            ne_ = len(eq)
            cap = (1 << ne_) * len(highs)
            npk = max((len(ws) + cap - 1) // cap for ws in d.values())
            for k in range(npk):
                pkd = dict(pat=pat, wire_at={})
                for ori, ws in sorted(d.items()):
                    for c, w in enumerate(ws[k * cap:(k + 1) * cap]):
                        s = 0
                        for idx, i in enumerate(uneq):
                            s |= ori[idx] << i
                        ce, ch = c % (1 << ne_), c >> ne_
                        for idx, i in enumerate(eq):
                            s |= ((ce >> idx) & 1) << i
                        s |= highs[ch] << u
                        assert H <= s < Q - H and s not in pkd["wire_at"]
                        pkd["wire_at"][s] = w
                        self.where[w] = (len(self.packets), s)
                self.packets.append(pkd)
        for w, (rv, cv) in enumerate(net.wires):
            P, s = self.where[w]
            for i in range(u):
                assert self.row_sym(P, s, i) == labels[rv][i]
                assert self.col_sym(P, s, i) == labels[cv][i]
        self.dens = {}
        for y in range(Q):
            dd = 1
            for j in range(Q):
                if j != y:
                    dd *= (y - j)
            self.dens[y] = dd
        self.v = Q - 1 - a
        assert all(v2(dd) == self.v for dd in self.dens.values())
        self._wc, self._rc = {}, {}

    def row_sym(self, P, x, i):
        lo, hi = self.packets[P]["pat"][i]
        return lo if lo == hi or not (x >> i) & 1 else hi

    def col_sym(self, P, y, i):
        lo, hi = self.packets[P]["pat"][i]
        return lo if lo == hi else (hi if not (y >> i) & 1 else lo)

    def Wfrac(self, x, y):
        """W~_xy(T) = l_y(x + T) = prod_(j != y) (T + x - j) / P_Q'(y): Fraction coefficients, degree Q - 1."""
        if (x, y) not in self._wc:
            p = [1]
            for j in range(self.L):
                if j != y:
                    p = poly_mul(p, [x - j, 1])
            self._wc[(x, y)] = [Fr(c, self.dens[y]) for c in p]
        return self._wc[(x, y)]

    def rowfactor(self, x, t, alpha, q):
        """2^-v prod_(j < Q) (t + alpha z + x - j) mod z^q (exact; coefficientwise divisible by 2^(a+1))."""
        key = (x, t, q)
        if key not in self._rc:
            p = [1]
            for j in range(self.L):
                p = poly_mul(p, [x - j, 1])
            self._rc[key] = [c / (1 << self.v) for c in taylor_at([Fr(c) for c in p], t, alpha, q)]
        return self._rc[key]

    def colunit(self, y):
        """The odd unit 2^v / P_Q'(y)."""
        return Fr(1 << self.v, self.dens[y])


def guarded_multiplicities(pk, h):
    """Tapered multiplicities on the guarded node grid: q_t = h (u - v2(t)) for t in [-H, H), t != 0, v2(t) < u."""
    return {t: h * (pk.u - v2(t)) for t in pk.S if t and v2(t) < pk.u}


def guarded_h(M, pk):
    """Smallest h with sum_t q_t > M (Q - 1), the degree of F~."""
    return M * (pk.L - 1) // sum((pk.H >> i) * (pk.u - i) for i in range(pk.u)) + 1


def setup_g(rng, n, m, a, u, kappa, r, chunks=False, Lc=2):
    cls = random_formula(n, m, rng)
    net = chunk_net_mapped(n, cls, make_units(n, cls, Lc, rng)) if chunks else formula_net(n, cls)
    labels = [tuple(rng.randrange(r) for _ in range(u)) for _ in net.verts]
    return n, cls, net, labels, GuardedPacking(net, labels, a, u, kappa)


def direct_taylor_mod_g(net, pk, t, p, q):
    """F~(t + alpha z) mod (2^p, z^q) straight from (D.8) with the guarded deformation."""
    m, alpha, M = 1 << p, 2 ** (pk.a + 1), len(net.wires)
    Wres = {}
    for P, pkd in enumerate(pk.packets):
        for x in pkd["wire_at"]:
            for y in pkd["wire_at"]:
                ser = taylor_at(pk.Wfrac(x, y), t, alpha, q)
                assert all(v2(c) >= 0 for c in ser), "entry not 2-adically integral"
                Wres[(P, x, y)] = [red(c, p) for c in ser]
    tot = [0] * q
    for xx in range(1 << (2 * M)):
        X = [(xx >> w) & 1 for w in range(M)]
        Y = [(xx >> (M + w)) & 1 for w in range(M)]
        tv = 1
        for vtx in net.verts:
            tv *= vtx["table"].get(tuple(Y[w] if e == "r" else X[w] for w, e in vtx["ports"]), 0)
            if tv == 0:
                break
        if tv == 0:
            continue
        ser = [tv % m] + [0] * (q - 1)
        for w in range(M):
            P, x = pk.where[w]
            form = [0] * q
            for y, w2 in pk.packets[P]["wire_at"].items():
                if X[w2]:
                    form = ps_add_mod(form, Wres[(P, x, y)], m)
            f = form if Y[w] else [(1 - form[0]) % m] + [(-c) % m for c in form[1:]]
            ser = ps_mul_mod(ser, f, q, m)
        tot = ps_add_mod(tot, ser, m)
    return tot


class NodeExpansionG(NodeExpansion):
    """NodeExpansion for the guarded deformation (LDU basis, sharp weights); the sums are inherited."""

    def __init__(self, net, pk, labels, t, p, q, weights="sharp", basis="ldu"):
        assert basis == "ldu" and weights == "sharp"
        self.basis = basis
        self.net, self.pk, self.t, self.p, self.q = net, pk, t, p, q
        self.m = 1 << p
        a = pk.a
        self.alpha = 2 ** (a + 1)
        i = self.i = v2(t)
        assert t and -pk.H <= t < pk.H and i < pk.u
        self.grp = [labels[v][i] for v in range(len(net.verts))]
        M = len(net.wires)
        self.Wser = {}
        for P, pkd in enumerate(pk.packets):
            for x in pkd["wire_at"]:
                for y in pkd["wire_at"]:
                    s = taylor_at(pk.Wfrac(x, y), t, self.alpha, q)
                    assert s[0] == (1 if y == x + t else 0), "not a shift at the node"
                    self.Wser[(P, x, y)] = s
        self.rowg = {w: self.grp[net.wires[w][0]] for w in range(M)}
        self.colg = {w: self.grp[net.wires[w][1]] for w in range(M)}
        self.modes, self.coef, self.local = [], {}, {}
        for w in range(M):
            P, x = pk.where[w]
            self.local[w] = []
            for y, w2 in pk.packets[P]["wire_at"].items():
                if self.colg[w2] == self.rowg[w]:
                    self.local[w].append((w2, [red(c, p) for c in self.Wser[(P, x, y)]]))
        step = 1 << (i + 1)
        half = 1 << i
        n = pk.L // 2
        for P, pkd in enumerate(pk.packets):
            lo, hi = pkd["pat"][i]
            if lo == hi:
                continue
            for h in (0, 1):
                pos = [step * (k >> i) + (k & (half - 1)) + h * half for k in range(n)]
                src = pk.col_sym(P, pos[0], i)
                tgt = hi if src == lo else lo
                occ = [(k, pkd["wire_at"][pos[k]]) for k in range(n) if pos[k] in pkd["wire_at"]]
                rows = [(k, w) for k, w in occ if pk.row_sym(P, pos[k], i) == tgt]
                if not occ or not rows:
                    continue
                K = [[ser_inv(([Fr(t + pos[r_] - pos[s_]), Fr(self.alpha)] + [Fr(0)] * q)[:q], q) for s_ in range(n)]
                     for r_ in range(n)]
                Lm, d, Um = series_ldu(K, q)
                for mm in range(n):
                    b = mm >> i
                    wgt = a + 1 - i + 2 * b + 2 * v2fact(b)
                    if wgt >= p:
                        break
                    lin = []
                    for k, w2 in occ:
                        if k >= mm:
                            cs = [c * pk.colunit(pos[k]) for c in Um[mm][k]]
                            assert all(v2(c) >= 0 for c in cs), "U * unit not 2-adically integral"
                            cser = tuple(red(c, p) for c in cs)
                            if any(cser):
                                lin.append((cser, w2))
                    if not lin:
                        continue
                    kk = len(self.modes)
                    self.modes.append(dict(P=P, y0=h, b=mm, src=src, tgt=tgt, w=wgt, lin=lin))
                    for k, wr in rows:
                        if k < mm:
                            continue
                        Rk = pk.rowfactor(pos[k], t, self.alpha, q)
                        assert all(v2(c) >= a + 1 for c in Rk if c), "row factor below 2^(a+1)"
                        assert all(v2(c) >= 0 for c in Lm[k][mm]), "L not 2-adically integral"
                        c = ser_mul(ser_mul(Rk, d[mm], q), Lm[k][mm], q)
                        assert all(v2(cc) >= wgt for cc in c), ("guarded weight", mm, wgt)
                        self.coef[(wr, kk)] = [red(cc, p) for cc in c]
        for P, pkd in enumerate(pk.packets):             # (D.12) for every crossing entry
            for x, wr in pkd["wire_at"].items():
                for y, wc in pkd["wire_at"].items():
                    if self.colg[wc] == self.rowg[wr]:
                        continue
                    acc = [0] * q
                    for k, md in enumerate(self.modes):
                        if md["P"] == P and (wr, k) in self.coef:
                            cy = dict((w2, c) for c, w2 in md["lin"]).get(wc)
                            if cy is not None:
                                acc = ps_add_mod(acc, ps_mul_mod(list(cy), self.coef[(wr, k)], q, self.m), self.m)
                    assert acc == [red(c, p) for c in self.Wser[(P, x, y)]], ("D.12 failed", P, x, y)
        self.r_groups = sorted(set(self.grp))
        self.tg = {g: [k for k, md in enumerate(self.modes) if md["tgt"] == g] for g in self.r_groups}
        self.sg = {g: [k for k, md in enumerate(self.modes) if md["src"] == g] for g in self.r_groups}
        self.inc = {g: sorted(set(self.tg[g]) | set(self.sg[g])) for g in self.r_groups}
        self.child_calls = 0


def count_e2e_g(n, clauses, chunks=True, Lc=2, r=3, seed=0, a=3, u=1, kappa=2, method="dp"):
    """Exact #F through the guarded deformation (reference)."""
    rng = random.Random(seed)
    net = chunk_net_mapped(n, clauses, make_units(n, clauses, Lc, rng)) if chunks else formula_net(n, clauses)
    labels = [tuple(rng.randrange(r) for _ in range(u)) for _ in net.verts]
    pk = GuardedPacking(net, labels, a, u, kappa)
    M = len(net.wires)
    qm = guarded_multiplicities(pk, guarded_h(M, pk))
    p = n + 1
    data = {}
    for t, qt in sorted(qm.items()):
        Pt = p + loss_bits(qm, a, t)
        ne = NodeExpansionG(net, pk, labels, t, Pt, qt)
        ser = ne.sum_mm()[0] if method == "mm" else ne.sum_dp()[0]
        data[t] = ser
    got = red(recover(data, qm, a), p)
    return (got * pow(2, net.unused, 1 << p)) % (1 << p)
