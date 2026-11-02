#!/usr/bin/env python3
"""Agent-editable implementation of the Appendix D route (T15). Everything up to NodeExpansion.__init__ is
provided infrastructure from the head's toy implementation (exact arithmetic, networks, chunk super-vertices,
packing, the Cauchy deformation, tapered Hermite recovery, the mode construction of D.12). Implement the
functions that raise NotImplementedError: NodeExpansion.group_table, sum_paper, sum_memo, sum_dp, the module
functions byvars_value, node_value and count_e2e. The referee compares them with frozen reference values on
hidden formulas; direct_taylor_mod and F_poly are blocked while it does so."""
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


def v2fact(b):
    """v2(b!) = b - (number of ones in the binary form of b)."""
    return b - bin(b).count("1")


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
    """K = L diag(d) U over Q[[z]] mod z^q (unit triangular L, U, no pivoting); entries are lists of q Fractions."""
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


def ps_sub_mod(a, b, m):
    return [(x - y) % m for x, y in zip(a, b)]


# ----------------------------------------------------------------------------- shell matrix products (appd-shellmm)
def naive_matmul(A, B, q, m):
    """Schoolbook product of an n x k matrix A of series and a k x l matrix B of series (list of q ints mod m):
    n rows of l series; every ring product goes through the module-level ps_mul_mod."""
    n, l = len(A), len(B[0]) if B else 0
    C = [[[0] * q for _ in range(l)] for _ in range(n)]
    for i in range(n):
        Ai, Ci = A[i], C[i]
        for t in range(len(B)):
            a = Ai[t]
            if not any(a):
                continue
            Bt = B[t]
            for j in range(l):
                b = Bt[j]
                if not any(b):
                    continue
                Ci[j] = ps_add_mod(Ci[j], ps_mul_mod(a, b, q, m), m)
    return C


def _mat_add(X, Y, m):
    return [[ps_add_mod(x, y, m) for x, y in zip(rowx, rowy)] for rowx, rowy in zip(X, Y)]


def _mat_sub(X, Y, m):
    return [[ps_sub_mod(x, y, m) for x, y in zip(rowx, rowy)] for rowx, rowy in zip(X, Y)]


def _quad(X, r0, r1, c0, c1, q, m):
    """The rows r0:r1 and columns c0:c1 of X as a matrix; zero series where X has no entry."""
    out = []
    for r in range(r0, r1):
        row = X[r] if 0 <= r < len(X) else None
        out.append([row[c] if (row is not None and 0 <= c < len(row)) else [0] * q for c in range(c0, c1)])
    return out


def strassen_sq(X, Y, q, m):
    """S x S matrix product with Strassen's algorithm (S a power of two, entries series mod m); falls back to
    naive_matmul for S <= 2. Only ring additions, subtractions and ps_mul_mod products are used."""
    S = len(X)
    if S <= 2:
        return naive_matmul(X, Y, q, m)
    h = S // 2
    A11 = _quad(X, 0, h, 0, h, q, m)
    A12 = _quad(X, 0, h, h, S, q, m)
    A21 = _quad(X, h, S, 0, h, q, m)
    A22 = _quad(X, h, S, h, S, q, m)
    B11 = _quad(Y, 0, h, 0, h, q, m)
    B12 = _quad(Y, 0, h, h, S, q, m)
    B21 = _quad(Y, h, S, 0, h, q, m)
    B22 = _quad(Y, h, S, h, S, q, m)
    M1 = strassen_sq(_mat_add(A11, A22, m), _mat_add(B11, B22, m), q, m)
    M2 = strassen_sq(_mat_add(A21, A22, m), B11, q, m)
    M3 = strassen_sq(A11, _mat_sub(B12, B22, m), q, m)
    M4 = strassen_sq(A22, _mat_sub(B21, B11, m), q, m)
    M5 = strassen_sq(_mat_add(A11, A12, m), B22, q, m)
    M6 = strassen_sq(_mat_sub(A21, A11, m), _mat_add(B11, B12, m), q, m)
    M7 = strassen_sq(_mat_sub(A12, A22, m), _mat_add(B21, B22, m), q, m)
    C11 = _mat_add(_mat_sub(_mat_add(M1, M4, m), M5, m), M7, m)
    C12 = _mat_add(M3, M5, m)
    C21 = _mat_add(M2, M4, m)
    C22 = _mat_add(_mat_sub(_mat_add(M1, M6, m), M2, m), M3, m)
    return [C11[r] + C12[r] for r in range(h)] + [C21[r] + C22[r] for r in range(h)]


def _slice(row, c0, c1, q):
    """The entries of row in the columns c0:c1 as series; zero series where row has no entry."""
    return [row[c] if (row is not None and 0 <= c < len(row) and row[c] is not None) else [0] * q
            for c in range(c0, c1)]


def matmul(A, B, q, m):
    """Product of an n x k matrix A and a k x l matrix B of series mod m, by Strassen on S x S blocks with
    S the largest power of two <= min(n, k, l); zero series pad the blocks outside A and B."""
    n, k = len(A), len(B)
    l = len(B[0]) if k else 0
    S = 1 << (min(n, k, l).bit_length() - 1)
    out = [[None] * ((l + S - 1) // S) for _ in range((n + S - 1) // S)]
    for bi in range(len(out)):
        r0, r1 = bi * S, min(bi * S + S, n)
        Arows = [row[0:k] for row in A[r0:r1]]
        Arows += [[None] * k for _ in range(S - (r1 - r0))]
        for bj in range(len(out[bi])):
            c0, c1 = bj * S, min(bj * S + S, l)
            acc = None
            for bk in range((k + S - 1) // S):
                g0, g1 = bk * S, min(bk * S + S, k)
                Ab = [_slice(row, g0, g1, q) + [[0] * q] * (S - (g1 - g0)) for row in Arows]
                Bb = [_slice(row, c0, c1, q) for row in B[g0:g1]] + \
                     [[[0] * q] * (c1 - c0) for _ in range(S - (g1 - g0))]
                P = strassen_sq(Ab, Bb, q, m)
                if acc is None:
                    acc = P
                else:
                    acc = _mat_add(acc, P, m)
            if acc is None:
                acc = [[[0] * q for _ in range(S)] for _ in range(S)]
            out[bi][bj] = [row[c0 - bj * S:c1 - bj * S] for row in acc[:r1 - r0]]
    res = []
    for bi in range(len(out)):
        for r in range(len(out[bi][0])):
            row = []
            for bj in range(len(out[0])):
                row.extend(out[bi][bj][r])
            res.append(row)
    return res


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
        assert basis in ("power", "newton", "ldu"), "unknown basis " + str(basis)
        assert not (basis == "ldu" and weights != "sharp"), "the ldu basis has only the sharp weights"
        self.net, self.pk, self.t, self.p, self.q = net, pk, t, p, q
        self.basis = basis
        self.guarded = getattr(pk, "guarded", False)
        if self.guarded:
            assert basis == "ldu" and weights == "sharp", "the guarded deformation uses the LDU basis only"
        self.m = 1 << p
        a, L = pk.a, pk.L
        self.alpha = 2 ** (a + 1)
        i = v2(t)
        self.i = i
        self.grp = [labels[v][i] for v in range(len(net.verts))]
        M = len(net.wires)
        # exact series of C'(0) W_xy(t + alpha z) / C'(0)
        self.Wser = {}
        if self.guarded:
            # guarded: W~_xy(t + alpha z) itself, only for the pairs of real slots
            for P, pkd in enumerate(pk.packets):
                for x in pkd["wire_at"]:
                    for y in pkd["wire_at"]:
                        self.Wser[(P, x, y)] = taylor_at(pk.Wfrac(x, y), t, self.alpha, q)
        else:
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
        for P, pkd in enumerate(pk.packets if basis != "ldu" else []):
            lo, hi = pkd["pat"][i]
            if lo == hi:
                continue
            for y0 in range(step):
                cols = [(j, pkd["wire_at"][y0 + step * j]) for j in range(L // step)
                        if (y0 + step * j) in pkd["wire_at"]]
                src = pk.col_sym(P, y0, i)
                tgt = hi if src == lo else lo
                rows = [(x, w) for x, w in pkd["wire_at"].items() if pk.row_sym(P, x, i) == tgt]
                if not cols or not rows:
                    continue
                for b in range(p):
                    wgt = ((b + 1) if weights == "paper" else (a + 1 - i + b))
                    if basis == "newton":
                        wgt += v2fact(b)
                    if wgt >= p:
                        break
                    # L_Z(X) = sum_j j^b X_y (power) or sum_j binom(j, b) X_y (Newton); trivial if all vanish
                    lin = [(self.bfun(j, b), w2) for j, w2 in cols]
                    lin = [(c, w2) for c, w2 in lin if c]
                    if not lin:
                        continue
                    k = len(self.modes)
                    self.modes.append(dict(P=P, y0=y0, b=b, src=src, tgt=tgt, w=wgt, lin=lin))
                    for x, wr in rows:
                        A = Fr(t - y0 + x)
                        assert v2(A) <= i
                        first = taylor_at([Fr(c, pk.Cp0) for c in pk.Wpoly(y0 - x)], t, self.alpha, q)
                        if basis == "newton":
                            # N_b = N_(b-1) * c*b / (A - c*b + alpha z), N_0 = first
                            cc = first
                            for bb in range(1, b + 1):
                                den = [A - Fr(step * bb), Fr(self.alpha)] + [Fr(0)] * max(0, q - 2)
                                rb = [Fr(step * bb) * d for d in ser_inv(den[:q], q)]
                                cc = ser_mul(cc, rb, q)
                        else:
                            base = [A, Fr(self.alpha)] + [Fr(0)] * max(0, q - 2)
                            invb = ser_inv(base[:q], q)
                            ratio = [Fr(2 ** (i + 1)) * c for c in invb]
                            c = first
                            for _ in range(b):
                                c = ser_mul(c, ratio, q)
                            cc = c
                        assert all(v2(d) >= wgt for d in cc), ("weight", wgt, [v2(d) for d in cc])
                        self.coef[(wr, k)] = [red(d, p) for d in cc]
        if basis == "ldu":
            # Cauchy LDU of the crossing block: W_(x_k, x_l)(t + alpha z) = G * K_kl, K_kl = 1/(T + x_k - x_l),
            # K = L diag(d) U over Z_(2)[[z]]; mode m has source form (U X)_m = sum_(l >= m) U_ml X_(x_l) and
            # coefficient G L_km d_m in row k >= m.
            half = 1 << i
            n = L // 2
            if not self.guarded:
                Cpoly = [1]
                for sv in pk.S:
                    Cpoly = poly_mul(Cpoly, [-sv, 1])
                G = taylor_at([Fr(c, pk.Cp0) for c in Cpoly], t, self.alpha, q)
            for P, pkd in enumerate(pk.packets):
                lo, hi = pkd["pat"][i]
                if lo == hi:
                    continue
                for h in (0, 1):
                    xk = [step * (k >> i) + (k & (half - 1)) + h * half for k in range(n)]
                    src = pk.col_sym(P, xk[0], i)
                    tgt = hi if src == lo else lo
                    occ = [(k, pkd["wire_at"][x]) for k, x in enumerate(xk) if x in pkd["wire_at"]]
                    rows = [(k, w) for (k, w) in occ if pk.row_sym(P, xk[k], i) == tgt]
                    if not occ or not rows:
                        continue
                    K = [[ser_inv(([Fr(t + xk[r] - xk[s]), Fr(self.alpha)] + [Fr(0)] * q)[:q], q)
                          for s in range(n)] for r in range(n)]
                    Lm, d, Um = series_ldu(K, q)
                    for m in range(n):
                        b = m >> i
                        wgt = a + 1 - i + 2 * b + 2 * v2fact(b)
                        if wgt >= p:
                            break
                        lin = []
                        for (k, w2) in occ:
                            if k < m:
                                continue
                            if self.guarded:
                                cf = [c * pk.colunit(xk[k]) for c in Um[m][k]]
                            else:
                                cf = Um[m][k]
                            assert all(v2(c) >= 0 for c in cf), ("U not integral", P, h, m, k)
                            tup = tuple(red(c, p) for c in cf)
                            if any(tup):
                                lin.append((tup, w2))
                        if not lin:
                            continue
                        mi = len(self.modes)
                        self.modes.append(dict(P=P, y0=h, b=m, src=src, tgt=tgt, w=wgt, lin=lin))
                        for (k, wr) in rows:
                            if k < m:
                                continue
                            assert all(v2(c) >= 0 for c in Lm[k][m]), ("L not integral", P, h, m, k)
                            if self.guarded:
                                R_k = pk.rowfactor(xk[k], t, self.alpha, q)
                                cc = ser_mul(ser_mul(R_k, d[m], q), Lm[k][m], q)
                            else:
                                cc = ser_mul(ser_mul(G, d[m], q), Lm[k][m], q)
                            assert all(v2(c) >= wgt for c in cc), ("weight", wgt, [v2(c) for c in cc])
                            self.coef[(wr, mi)] = [red(c, p) for c in cc]
        # (D.12) check: every crossing entry equals its mode expansion mod 2^p
        for P, pkd in enumerate(pk.packets):
            for x, wr in pkd["wire_at"].items():
                for y, wc in pkd["wire_at"].items():
                    if self.colg[wc] == self.rowg[wr]:
                        continue
                    acc = [0] * q
                    for k, md in enumerate(self.modes):
                        if md["P"] != P or (wr, k) not in self.coef:
                            continue
                        if basis == "ldu":
                            cb = None
                            for c, w2 in md["lin"]:
                                if w2 == wc:
                                    cb = c
                                    break
                            if cb is None:
                                continue
                            acc = ps_add_mod(acc, ps_mul_mod(list(cb), self.coef[(wr, k)], q, self.m), self.m)
                        else:
                            y0, j = y % step, y // step
                            if md["y0"] != y0:
                                continue
                            jb = self.bfun(j, md["b"])
                            acc = ps_add_mod(acc, [(jb * c) % self.m for c in self.coef[(wr, k)]], self.m)
                    direct = [red(c, p) for c in self.Wser[(P, x, y)]]
                    assert acc == direct, ("D.12 failed", P, x, y)
        self.r_groups = sorted(set(self.grp))
        self.tg = {g: [k for k, md in enumerate(self.modes) if md["tgt"] == g] for g in self.r_groups}
        self.sg = {g: [k for k, md in enumerate(self.modes) if md["src"] == g] for g in self.r_groups}
        self.inc = {g: sorted(set(self.tg[g]) | set(self.sg[g])) for g in self.r_groups}
        self.child_calls = 0

    # -- basis functions
    def bfun(self, j, b):
        """Value of the basis function of index b at the column index j: j**b for the power basis,
        math.comb(j, b) for the Newton (binomial) basis; not used by the Cauchy LDU basis, whose source forms
        carry z-series coefficients in self.modes[k]["lin"]."""
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
        """V_g[gamma] for every tuple gamma on I_g = self.inc[g] with weight < self.p (appd-memo): the contraction of
        the child network of group g (D.4.6) for that gamma, computed once per distinct gamma. Return a dict
        {self.key(gamma): series} (series = list of self.q residues mod 2^self.p). self.child_network_value(g, gamma)
        is a slow reference for one gamma. Add the number of entries to self.child_calls."""
        gammas = self.tuples(self.inc[g], self.p)
        tbl = {self.key(gm): self.child_network_value(g, gm) for gm in gammas}
        self.child_calls += len(gammas)
        return tbl

    @staticmethod
    def key(beta):
        return tuple(sorted(beta.items()))

    def sum_paper(self):
        """(D.16) as in the paper: F(t + alpha z) = sum over all tuples beta (weight < p) of prod_g Z_{g,beta}
        mod (2^p, z^q). Return (series, number of tuples, number of tuples * number of groups)."""
        tables = {g: self.group_table(g) for g in self.r_groups}
        total = [0] * self.q
        nb = 0
        for beta in self.tuples(list(range(len(self.modes))), self.p):
            nb += 1
            val = [1] + [0] * (self.q - 1)
            for g in self.r_groups:
                sub = {k: e for k, e in beta.items() if k in self.inc[g]}
                val = ps_mul_mod(val, tables[g][self.key(sub)], self.q, self.m)
            total = ps_add_mod(total, val, self.m)
        return total, nb, nb * len(self.r_groups)

    def sum_memo(self):
        """Same sum, with the child tables of group_table computed once per group (appd-memo).
        Return (series, number of tuples, self.child_calls) with self.child_calls reset to 0 first."""
        self.child_calls = 0
        tables = {g: self.group_table(g) for g in self.r_groups}
        total = [0] * self.q
        nb = 0
        for beta in self.tuples(list(range(len(self.modes))), self.p):
            nb += 1
            val = [1] + [0] * (self.q - 1)
            for g in self.r_groups:
                sub = {k: e for k, e in beta.items() if k in self.inc[g]}
                val = ps_mul_mod(val, tables[g][self.key(sub)], self.q, self.m)
            total = ps_add_mod(total, val, self.m)
        return total, nb, self.child_calls

    def sum_dp(self):
        """Same sum, eliminated group by group in the order self.r_groups (appd-modedp): states are tuples on the
        cut (modes with one end processed and one not); the step for group g enumerates the state tuples together
        with every tuple on the new modes between g and later groups, skipping combinations of weight >= p,
        multiplies by the table entry of g and projects out the modes whose both ends are processed.
        Return (series, max_enum, self.child_calls) with self.child_calls reset to 0 first; max_enum is the largest
        number, over the steps, of (state, new tuple) combinations of total weight < p that the step multiplies."""
        self.child_calls = 0
        tables = {g: self.group_table(g) for g in self.r_groups}
        pos = {g: j for j, g in enumerate(self.r_groups)}
        endpos = {}
        for k, md in enumerate(self.modes):
            endpos[k] = (min(pos[md["src"]], pos[md["tgt"]]), max(pos[md["src"]], pos[md["tgt"]]))
        Phi = {(): [1] + [0] * (self.q - 1)}
        max_enum = 0
        for j, g in enumerate(self.r_groups):
            new = sorted(k for k in self.inc[g] if endpos[k][0] == j and endpos[k][1] > j)
            newset = set(new)
            keep = sorted(k for k, (lo, hi) in endpos.items() if hi > j)
            keepset = set(keep)
            newPhi = {}
            enum = 0
            for s, ser in Phi.items():
                sd = dict(s)
                for nt in self.tuples(new, self.p - self.weight(sd)):
                    ntd = dict(nt)
                    full = dict(sd)
                    full.update(ntd)
                    if self.weight(full) >= self.p:
                        continue
                    enum += 1
                    gamma = {k: e for k, e in full.items() if k in self.inc[g]}
                    val = ps_mul_mod(ser, tables[g][self.key(gamma)], self.q, self.m)
                    proj = tuple(sorted((k, e) for k, e in full.items() if k in keepset))
                    newPhi[proj] = ps_add_mod(newPhi[proj], val, self.m) if proj in newPhi else val
            max_enum = max(max_enum, enum)
            Phi = newPhi
        return Phi.get((), [0] * self.q), max_enum, self.child_calls

    def shells(self, idx):
        """self.tuples(idx, self.p) put into shells by their exact weight self.weight(...); the shells are returned
        in increasing weight (each a list of dicts)."""
        byw = {}
        for t in self.tuples(idx, self.p):
            byw.setdefault(self.weight(t), []).append(t)
        return [byw[w] for w in sorted(byw)]

    def sum_mm(self):
        """Same sum as sum_dp (appd-shellmm), group by group, as shell matrix products: at the step of group g the
        cut modes split into X (not touching g), Y (touching g) and the new modes Z (g to later groups); each of X,
        Y and Z is grouped by exact weight (shells) and every shell triple (w1, w2, w3) with w1 + w2 + w3 < p is one
        full block product matmul(A[X_w1 x Y_w2], B[Y_w2 x Z_w3]) whose entries are added into the new state.
        Return (series, shapes, self.child_calls) with self.child_calls reset to 0 first; shapes is the list, in step
        order, of the sorted lists of the triples (len(xl), len(yl), len(zl)) of the step."""
        self.child_calls = 0
        tables = {g: self.group_table(g) for g in self.r_groups}
        pos = {g: j for j, g in enumerate(self.r_groups)}
        endpos = {}
        for k, md in enumerate(self.modes):
            endpos[k] = (min(pos[md["src"]], pos[md["tgt"]]), max(pos[md["src"]], pos[md["tgt"]]))
        zero = [0] * self.q
        Phi = {(): [1] + [0] * (self.q - 1)}
        shapes = []
        for j, g in enumerate(self.r_groups):
            incg = set(self.inc[g])
            cut = sorted(k for k, (lo, hi) in endpos.items() if lo < j <= hi)
            X = [k for k in cut if k not in incg]
            Y = [k for k in cut if k in incg]
            Z = sorted(k for k, (lo, hi) in endpos.items() if lo == j < hi)
            shX, shY, shZ = self.shells(X), self.shells(Y), self.shells(Z)
            new = {}
            step = []
            for xl in shX:
                w1 = self.weight(xl[0])
                for yl in shY:
                    w2 = self.weight(yl[0])
                    if w1 + w2 >= self.p:
                        continue
                    Am = []
                    for x_i in xl:
                        xd = dict(x_i)
                        Am.append([Phi.get(self.key({**xd, **dict(y_s)}), zero) for y_s in yl])
                    for zl in shZ:
                        w3 = self.weight(zl[0])
                        if w1 + w2 + w3 >= self.p:
                            continue
                        Bm = [[tables[g][self.key({**dict(y_s), **dict(z_j)})] for z_j in zl] for y_s in yl]
                        Cm = matmul(Am, Bm, self.q, self.m)
                        step.append((len(xl), len(yl), len(zl)))
                        for i, x_i in enumerate(xl):
                            xd = dict(x_i)
                            for jj, z_j in enumerate(zl):
                                key = self.key({**xd, **dict(z_j)})
                                val = Cm[i][jj]
                                new[key] = ps_add_mod(new[key], val, self.m) if key in new else val
            shapes.append(sorted(step))
            Phi = new
        return Phi.get((), [0] * self.q), shapes, self.child_calls

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
        if src_zero(f):
            return [0] * q
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
    gd = getattr(pk, "guarded", False)
    for P, pkd in enumerate(pk.packets):
        for x in pkd["wire_at"]:
            for y in pkd["wire_at"]:
                if gd:
                    Wres[(P, x, y)] = [red(c, p) for c in taylor_at(pk.Wfrac(x, y), t, alpha, q)]
                else:
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
    """appd-leafvars: V_g[gamma] of a chunk network computed by summing over the variables of the chunks of group g
    (net.units, net.unit_of, net.desig, net.clauses; wire -> (clause, variable) in _WIREMAP). Fixing the variables
    fixes every column wire of the group; each clause designated to a chunk of g then depends only on its external
    row bits: constant 1 if one of its internal literals is true, otherwise 1 - prod_j [external literal j false],
    two product terms. After expanding these terms every row bit y carries a unary weight (w0, w1), and
    sum_y w(y) l_y(form) = w0 + (w1 - w0) form enters that row's packet factor. Must equal ne.group_table(g)[key].
    Return (series, number of variables summed, max number of clauses split into two terms)."""
    net, pk, m, q = ne.net, ne.pk, ne.m, ne.q
    kg = [k for k in range(len(net.verts)) if ne.grp[k] == g]
    kset = set(kg)
    varset = set()
    for k in kg:
        varset |= set(net.units[k])
    vars_ = sorted(varset)
    cls = [ci for ci in range(len(net.clauses)) if net.desig[ci] in kset]
    lit = {}
    for ci in cls:
        for l in net.clauses[ci]:
            lit[(ci, abs(l))] = 1 if l > 0 else -1
    # column wires of the group: their bit is the value of the variable they carry
    colvar = {w: _WIREMAP[(id(net), w)][1] for w in range(len(net.wires)) if ne.colg[w] == g}
    # row wires of the group: one per external literal of a clause designated to a chunk of g
    rowinfo, extof = {}, {ci: [] for ci in cls}
    for w in range(len(net.wires)):
        if ne.rowg[w] == g:
            ci, v = _WIREMAP[(id(net), w)]
            rowinfo[w] = (ci, v, lit[(ci, v)])
            extof[ci].append(w)
    rows = sorted(rowinfo)
    tmodes = [k for k in ne.tg[g]]
    smodes = [k for k in ne.sg[g]]
    want = tuple(sorted((k, gamma[k]) for k in tmodes if k in gamma))
    rowmodes = {w: [k for k in tmodes if (w, k) in ne.coef] for w in rows}
    total = [0] * q
    nsplit_max = 0
    for xb in itertools.product((0, 1), repeat=len(vars_)):
        val = dict(zip(vars_, xb))
        split = [ci for ci in cls
                 if not any(val[abs(l)] == (1 if l > 0 else 0)
                            for l in net.clauses[ci] if abs(l) in net.units[net.desig[ci]])]
        nsplit_max = max(nsplit_max, len(split))
        for choice in itertools.product((0, 1), repeat=len(split)):
            w01 = {w: (1, 1) for w in rows}
            sign = 1
            for ci, ch in zip(split, choice):
                if not ch:
                    continue
                sign = -sign
                for w in extof[ci]:
                    w01[w] = (1, 0) if rowinfo[w][2] > 0 else (0, 1)
            polyd = {(): [1] + [0] * (q - 1)}
            for w in rows:
                loc = [0] * q
                for w2, ser in ne.local[w]:
                    if val[colvar[w2]]:
                        loc = ps_add_mod(loc, ser, m)
                a0 = (w01[w][1] - w01[w][0]) % m
                base = [(a0 * c) % m for c in loc]
                base[0] = (base[0] + w01[w][0]) % m
                terms = [((), base)] + [(k, [(a0 * c) % m for c in ne.coef[(w, k)]]) for k in rowmodes[w]]
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
                    f = src_mul(f, lin_eval(ne.modes[k]["lin"], lambda w2: val[colvar[w2]], q, m), e, q, m)
            if src_zero(f):
                continue
            term = src_scale(f, ser, q, m)
            if sign < 0:
                term = [(-c) % m for c in term]
            total = ps_add_mod(total, term, m)
    return total, len(vars_), nsplit_max


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


IMPLEMENTED = True      # set to True once every function below works


def node_value(net, pk, labels, t, p, q, weights="paper", method="dp", basis="power"):
    """F(t + alpha z) mod (2^p, z^q) from the mode expansion at node t (v2(t) < pk.a): method "paper", "memo"
    or "dp" selects NodeExpansion.sum_paper / sum_memo / sum_dp; basis "power", "newton" or "ldu" selects the
    mode basis (ldu needs weights="sharp"). Return (series, stats) where stats is {"tuples" or "max_enum": int, "child_calls": int,
    "modes": len(ne.modes)}."""
    ne = NodeExpansion(net, pk, labels, t, p, q, weights=weights, basis=basis)
    if method == "paper":
        ser, tup_n, calls = ne.sum_paper()
        stats = {"tuples": tup_n, "child_calls": calls, "modes": len(ne.modes)}
    elif method == "memo":
        ser, tup_n, calls = ne.sum_memo()
        stats = {"tuples": tup_n, "child_calls": calls, "modes": len(ne.modes)}
    elif method == "dp":
        ser, max_enum, calls = ne.sum_dp()
        stats = {"max_enum": max_enum, "child_calls": calls, "modes": len(ne.modes)}
    elif method == "cross":
        ser, records, calls = ne.sum_cross()
        stats = {"records": records, "child_calls": calls, "modes": len(ne.modes)}
    elif method == "mm":
        ser, shapes, calls = ne.sum_mm()
        stats = {"shapes": shapes,
                 "max_enum": max([sum(n * k * l for (n, k, l) in s) for s in shapes] or [0]),
                 "child_calls": calls, "modes": len(ne.modes)}
    else:
        raise ValueError("unknown method " + str(method))
    return ser, stats


def count_e2e(n, clauses, chunks=True, Lc=2, r=3, seed=0, a=1, u=1, basis="power", method="dp",
              guarded=False, kappa=2):
    """Exact #F by the whole algorithm with packet coordinates a and retained valuation levels u: build the
    network (chunk super-vertices when chunks, units from make_units(n, clauses, Lc, random.Random(seed))),
    labels of a coordinates tuple(rng.randrange(r)) * a per vertex from the same rng, packing,
    h = M (2L - 1) // sum_{i < u} (L >> i)(u - i) + 1 with L = 2^a, multiplicities(a, u, h); at every node t
    with q_t > 0 use precision P_t = p + loss_bits(...) with p = n + 1 and the 'sharp' weights with the given
    basis ("power", "newton" or "ldu"), F(t + alpha z) mod (2^P_t, z^q_t) from sum_dp, then recover(...) and reduce
    mod 2^p; multiply by 2^(unused variables). With a = 1, u = 1 this is the old a = 1, u = 1 call and must
    equal count_sat. With guarded=True the packing is GuardedPacking(net, labels, a, u, kappa) with labels of u
    coordinates, the nodes are guarded_multiplicities(pk, guarded_h(M, pk)) and the basis must be "ldu"."""
    rng = random.Random(seed)
    if chunks:
        net = chunk_net_mapped(n, clauses, make_units(n, clauses, Lc, rng))
    else:
        net = formula_net(n, clauses)
    lab_n = u if guarded else a
    labels = [tuple(rng.randrange(r) for _ in range(lab_n)) for _ in net.verts]
    pk = GuardedPacking(net, labels, a, u, kappa) if guarded else Packing(net, labels, a)
    M = len(net.wires)
    L = pk.L
    if guarded:
        qm = guarded_multiplicities(pk, guarded_h(M, pk))
    else:
        h = M * (2 * L - 1) // sum((L >> i) * (u - i) for i in range(u)) + 1
        qm = multiplicities(a, u, h)
    p = n + 1
    data = {}
    for t, qt in sorted(qm.items()):
        if qt == 0:
            continue
        Pt = p + loss_bits(qm, a, t)
        ne = NodeExpansion(net, pk, labels, t, Pt, qt, weights="sharp", basis=basis)
        if method == "mm":
            ser, _sh, _cc = ne.sum_mm()
        elif method == "cross" and len(ne.r_groups) >= 4:
            ser, _rec, _cc = ne.sum_cross()
        else:
            ser, _me, _cc = ne.sum_dp()
        data[t] = ser
    got = red(recover(data, qm, a), p)
    return (got * pow(2, net.unused, 1 << p)) % (1 << p)


# ------------------------------------------------------------------ binary crossing (round 10, appd-bincross)
# The 4 x 4 matrices of NOTES.md, section 'Round 10' (row index = output digit, column = input digit).
_BC_S = [[0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 1]]                    # S (row pairing)
_BC_KT = [[1, 0, 0, 1], [0, 0, -1, 0], [0, -1, 0, 0], [1, 0, 0, 0]]                  # K T (column pairing)
_BC_N0_ROW = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]]                # S^-1
_BC_N1_ROW = [[0, -1, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]               # -S^-1 J
_BC_N0_COL = [[1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]                # T^-1 P0
_BC_N1_COL = [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [1, 0, 0, -1]]            # -T^-1


def _bc_pair(pairing):
    """'row' or 'col' ('column' accepted) -> the pairing name used below."""
    if pairing == "row":
        return "row"
    if pairing in ("col", "column"):
        return "col"
    raise ValueError("unknown pairing " + str(pairing))


def _bc_mats(pairing):
    """(M, N0, N1) of NOTES.md for the given pairing; M = S (row) or K T (column)."""
    if _bc_pair(pairing) == "row":
        return _BC_S, _BC_N0_ROW, _BC_N1_ROW
    return _BC_KT, _BC_N0_COL, _BC_N1_COL


def _bc_rows(A):
    """A -> for every output digit d2 the (input digit d, coefficient) pairs with A[d2][d] != 0."""
    return [[(d, A[d2][d]) for d in range(4) if A[d2][d]] for d2 in range(4)]


def _bc_rows_T(A):
    """Rows of A^T: for every output digit d2 the (d, coefficient) pairs with A^T[d2][d] = A[d][d2] != 0."""
    return [[(d, A[d][d2]) for d in range(4) if A[d][d2]] for d2 in range(4)]


def _bc_apply(old, m, k, rows, q, mod):
    """One digit step: new[w2] = sum over d of rows[d2][d] * old[w2 with digit k replaced by d], where d2 is
    digit k of w2 (digit 0 is the most significant). None entries mean zero; words that come out zero are
    returned as None. The input list and its series are not modified."""
    n = 4 ** m
    shift = 2 * (m - 1 - k)
    clear = ~ (3 << shift)
    out = [None] * n
    for w in range(n):
        base = w & clear
        acc = None
        for d, c in rows[(w >> shift) & 3]:
            src = old[base | (d << shift)]
            if src is None:
                continue
            if acc is None:
                acc = [c * x for x in src]
            else:
                acc = [a + c * x for a, x in zip(acc, src)]
        if acc is not None:
            out[w] = [a % mod for a in acc]
    return out


def _bc_add_into(dst, src, mod):
    """dst[w] += src[w] (None = zero) mod mod, in place on the list dst; returns dst."""
    for w, s in enumerate(src):
        if s is None:
            continue
        cur = dst[w]
        if cur is None:
            dst[w] = s
        else:
            dst[w] = [(a + b) % mod for a, b in zip(cur, s)]
    return dst


def _bc_copy(vec):
    """A copy of a vector of series (None kept as None); the input list and its series are not touched."""
    return [None if v is None else [x for x in v] for v in vec]


def bc_mask(m, i, pairing):
    """The words w in range(4 ** m), increasing, with at least i zero base-4 digits when pairing == 'row', or
    at least m - i zero digits when pairing == 'col'. Digit k of w (k = 0 the most significant) is
    (w >> (2 * (m - 1 - k))) & 3."""
    need = i if _bc_pair(pairing) == "row" else m - i
    n = 4 ** m
    if need <= 0:
        return list(range(n))
    if need > m:
        return []
    out = []
    for w in range(n):
        z = 0
        t = w
        for _ in range(m):
            t, d = t >> 2, t & 3
            if d == 0:
                z += 1
        if z >= need:
            out.append(w)
    return out


def bc_left(vec, m, pairing, q, mod):
    """(M^T)^(x)m vec mod mod as a list of 4 ** m series with no None, M = S for 'row' and M = K T for 'col';
    one digit at a time with A = M^T, i.e. A[d2][d] = M[d][d2]. The input is never modified."""
    M = _bc_mats(pairing)[0]
    rows = _bc_rows_T(M)
    cur = _bc_copy(vec)
    if m == 0:
        s = cur[0] if cur[0] is not None else [0] * q
        return [[int(x) % mod for x in s]]
    for k in range(m):
        cur = _bc_apply(cur, m, k, rows, q, mod)
    return [[x for x in s] if s is not None else [0] * q for s in cur]


def bc_right(vec, m, pairing, q, mod):
    """R = a list of m + 1 lists of 4 ** m series with R[i] = [t^i] (N0 + t N1)^(x)m vec mod mod, (N0, N1) of
    NOTES.md for the pairing; here A = N0 or N1 directly (no transpose). The input is never modified."""
    _M, N0, N1 = _bc_mats(pairing)
    r0, r1 = _bc_rows(N0), _bc_rows(N1)
    cur = _bc_copy(vec)
    if m == 0:
        s = cur[0] if cur[0] is not None else [0] * q
        return [[[int(x) % mod for x in s]]]
    powers = [cur]
    for k in range(m):
        nxt = [None] * (m + 1)
        for i, v in enumerate(powers):
            if v is None:
                continue
            for j, rows in ((i, r0), (i + 1, r1)):
                res = _bc_apply(v, m, k, rows, q, mod)
                nxt[j] = res if nxt[j] is None else _bc_add_into(nxt[j], res, mod)
        powers = nxt
    return [[s if s is not None else [0] * q for s in p] for p in powers]


# ------------------------------------------------------------------ sum_cross (round 10, appd-bincross part 2)
def _bc_spread(bits, m):
    """Bits of `bits` (bit m-1 the most significant) spread into the base-4 digits of a word: digit k
    (k = 0 the most significant) = bit m-1-k; the digits are 0/1, so 2 * spread(first) + spread(second) is the
    word of a block whose digit k is 2 * (bit k of the first end's top bits) + (bit k of the second end's)."""
    w = 0
    for k in range(m):
        if (bits >> (m - 1 - k)) & 1:
            w |= 1 << (2 * (m - 1 - k))
    return w


def sum_cross(self):
    """Same sum as sum_dp (appd-bincross), one binary crossing with fast absorption: g1..g4 = r_groups[:4],
    X = the modes of the pair {g1, g3}, Y = the modes of {g2, g4} (they cross seen from {g1, g2} and {g3, g4});
    per shell pair the codes of the tuples give a = (len(sx) - 1).bit_length() and b = (len(sy) - 1).bit_length()
    bits, m = min(a, b) paired bits and extra = |a - b| free bits. Piece i of SWAP^(m) = sum_i A_i is one block
    pair product through bc_mask / bc_left / bc_right (column pairing when 2 i <= m, row pairing otherwise),
    summed over the block boundaries and the singletons. Return (series, records, self.child_calls) with
    self.child_calls reset to 0 first; records has one dict per (sx, sy, i) in loop order with keys wx, wy, nx,
    ny, a, b, m, i, pairing, channels, extra, value (the piece's series, entries mod 2^p)."""
    assert len(self.r_groups) >= 4, "sum_cross needs at least four groups"
    self.child_calls = 0
    p, q, mod = self.p, self.q, self.m
    g1, g2, g3, g4 = self.r_groups[:4]
    singles = self.r_groups[4:]
    tables = {g: self.group_table(g) for g in self.r_groups}
    X = sorted(k for k, md in enumerate(self.modes) if {md["src"], md["tgt"]} == {g1, g3})
    Y = sorted(k for k, md in enumerate(self.modes) if {md["src"], md["tgt"]} == {g2, g4})
    XY = set(X) | set(Y)
    others = [k for k in range(len(self.modes)) if k not in XY]
    egrp = {"x": g1, "y": g2, "u": g4, "v": g3}          # the four ends and the group that holds them
    blocks = {"row": (("x", "y"), ("u", "v")), "col": (("x", "u"), ("y", "v"))}
    one = [1] + [0] * (q - 1)
    shX = self.shells(X)
    shY = self.shells(Y)
    total = [0] * q
    records = []
    for sx in shX:
        wx = self.weight(sx[0])
        nx = len(sx)
        a = (nx - 1).bit_length()
        for sy in shY:
            wy = self.weight(sy[0])
            if wx + wy >= p:
                continue
            budget = p - wx - wy
            ny = len(sy)
            b = (ny - 1).bit_length()
            m = min(a, b)
            extra = abs(a - b)
            nE = 1 << extra
            nW = 1 << (2 * m)
            shell_t = {"x": sx, "v": sx, "y": sy, "u": sy}
            ebits = {"x": a, "v": a, "y": b, "u": b}
            used = sorted(set("col" if 2 * i <= m else "row" for i in range(m + 1)))
            side = {}
            for pairing in used:
                for which, (eA, eB) in enumerate(blocks[pairing]):
                    hA, hB = egrp[eA], egrp[eB]
                    bnd = sorted((set(self.inc[hA]) | set(self.inc[hB])) - XY)
                    bndset = set(bnd)
                    incA, incB = set(self.inc[hA]), set(self.inc[hB])
                    preA = []
                    for cA, tA in enumerate(shell_t[eA]):
                        sub = tuple(sorted((k, e) for k, e in tA.items() if k in incA))
                        top = _bc_spread(cA >> (ebits[eA] - m), m) if m else 0
                        ex = cA & ((1 << (ebits[eA] - m)) - 1) if ebits[eA] > m else 0
                        preA.append((top, ex, sub))
                    preB = []
                    for cB, tB in enumerate(shell_t[eB]):
                        sub = tuple(sorted((k, e) for k, e in tB.items() if k in incB))
                        top = _bc_spread(cB >> (ebits[eB] - m), m) if m else 0
                        ex = cB & ((1 << (ebits[eB] - m)) - 1) if ebits[eB] > m else 0
                        preB.append((top, ex, sub))
                    vecs = {}
                    for bt in self.tuples(bnd, budget):
                        items = tuple(sorted(bt.items()))          # sorted pairs, usable as a dict key
                        sA = tuple(x for x in items if x[0] in incA)
                        sB = tuple(x for x in items if x[0] in incB)
                        tA_of = []
                        for top, ex, sub in preA:
                            v = tables[hA].get(tuple(sorted(sA + sub)))
                            if v is not None:
                                tA_of.append((top, ex, v))
                        tB_of = []
                        for top, ex, sub in preB:
                            v = tables[hB].get(tuple(sorted(sB + sub)))
                            if v is not None:
                                tB_of.append((top, ex, v))
                        if not tA_of or not tB_of:
                            continue
                        for topA, exA, vA in tA_of:
                            for topB, exB, vB in tB_of:
                                word = 2 * topA + topB
                                ee = exA if ebits[eA] >= ebits[eB] else exB
                                vec = vecs.get((ee, items))
                                if vec is None:
                                    vec = vecs[(ee, items)] = [None] * nW
                                pr = ps_mul_mod(vA, vB, q, mod)
                                cur = vec[word]
                                vec[word] = pr if cur is None else [(x + y) % mod for x, y in zip(cur, pr)]
                    tr = []
                    for vec in vecs.values():
                        tr.append(bc_left(vec, m, pairing, q, mod) if which == 0
                                  else bc_right(vec, m, pairing, q, mod))
                    if pairing not in side:
                        side[pairing] = [None, None]
                    side[pairing][which] = (bndset, dict(zip(list(vecs), tr)))
            otup = self.tuples(others, budget)
            prep = {}
            for pairing in used:
                rows = []
                bndL, bndR = side[pairing][0][0], side[pairing][1][0]
                incs = [set(self.inc[sg]) for sg in singles]
                for bt in otup:
                    items = sorted(bt.items())
                    kA = tuple(x for x in items if x[0] in bndL)
                    kB = tuple(x for x in items if x[0] in bndR)
                    s = one
                    ok = True
                    for sg, ic in zip(singles, incs):
                        v = tables[sg].get(tuple(x for x in items if x[0] in ic))
                        if v is None:
                            ok = False
                            break
                        s = ps_mul_mod(s, v, q, mod)
                    if ok:
                        rows.append((kA, kB, s))
                prep[pairing] = rows
            for i in range(m + 1):
                pairing = "col" if 2 * i <= m else "row"
                mask = bc_mask(m, i, pairing)
                Ls = side[pairing][0][1]
                Rs = side[pairing][1][1]
                piece = [0] * q
                for kA, kB, s in prep[pairing]:
                    for e in range(nE):
                        lv = Ls.get((e, kA))
                        if lv is None:
                            continue
                        rv = Rs.get((e, kB))
                        if rv is None:
                            continue
                        acc = None
                        for w in mask:
                            a1, b1 = lv[w], rv[i][w]
                            if not any(a1) or not any(b1):
                                continue
                            pr = ps_mul_mod(a1, b1, q, mod)
                            acc = pr if acc is None else [(x + y) % mod for x, y in zip(acc, pr)]
                        if acc is None:
                            continue
                        pr = ps_mul_mod(s, acc, q, mod)
                        piece = [(x + y) % mod for x, y in zip(piece, pr)]
                total = ps_add_mod(total, piece, mod)
                records.append(dict(wx=wx, wy=wy, nx=nx, ny=ny, a=a, b=b, m=m, i=i, pairing=pairing,
                                    channels=len(mask), extra=extra, value=piece))
    return total, records, self.child_calls


NodeExpansion.sum_cross = sum_cross
