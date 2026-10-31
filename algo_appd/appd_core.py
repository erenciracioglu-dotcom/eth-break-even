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

    def __init__(self, net, pk, labels, t, p, q, weights="paper"):
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
        for P, pkd in enumerate(pk.packets):
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
                    wgt = (b + 1) if weights == "paper" else (a + 1 - i + b)
                    if wgt >= p:
                        break
                    # L_Z(X) = sum_j j^b X_y ; trivial if all j^b vanish
                    lin = [(j ** b if (j or b == 0) else 0, w2) for j, w2 in cols]
                    lin = [(c, w2) for c, w2 in lin if c]
                    if not lin:
                        continue
                    k = len(self.modes)
                    self.modes.append(dict(P=P, y0=y0, b=b, src=src, tgt=tgt, w=wgt, lin=lin))
                    for x, wr in rows:
                        A = Fr(t - y0 + x)
                        assert v2(A) <= i
                        first = taylor_at([Fr(c, pk.Cp0) for c in pk.Wpoly(y0 - x)], t, self.alpha, q)
                        base = [A, Fr(self.alpha)] + [Fr(0)] * max(0, q - 2)
                        invb = ser_inv(base[:q], q)
                        ratio = [Fr(2 ** (i + 1)) * c for c in invb]
                        c = first
                        for _ in range(b):
                            c = ser_mul(c, ratio, q)
                        assert all(v2(cc) >= wgt for cc in c), ("weight", wgt, [v2(cc) for cc in c])
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
                        if md["P"] == P and md["y0"] == y0 and (wr, k) in self.coef:
                            jb = j ** md["b"] if (j or md["b"] == 0) else 0
                            acc = ps_add_mod(acc, [(jb * c) % self.m for c in self.coef[(wr, k)]], self.m)
                    direct = [red(c, p) for c in self.Wser[(P, x, y)]]
                    assert acc == direct, ("D.12 failed", P, x, y)
        self.r_groups = sorted(set(self.grp))
        self.tg = {g: [k for k, md in enumerate(self.modes) if md["tgt"] == g] for g in self.r_groups}
        self.sg = {g: [k for k, md in enumerate(self.modes) if md["src"] == g] for g in self.r_groups}
        self.inc = {g: sorted(set(self.tg[g]) | set(self.sg[g])) for g in self.r_groups}
        self.child_calls = 0

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
        raise NotImplementedError("T15: implement group_table")

    @staticmethod
    def key(beta):
        return tuple(sorted(beta.items()))

    def sum_paper(self):
        """(D.16) as in the paper: F(t + alpha z) = sum over all tuples beta (weight < p) of prod_g Z_{g,beta}
        mod (2^p, z^q). Return (series, number of tuples, number of tuples * number of groups)."""
        raise NotImplementedError("T15: implement sum_paper")

    def sum_memo(self):
        """Same sum, with the child tables of group_table computed once per group (appd-memo).
        Return (series, number of tuples, self.child_calls) with self.child_calls reset to 0 first."""
        raise NotImplementedError("T15: implement sum_memo")

    def sum_dp(self):
        """Same sum, eliminated group by group in the order self.r_groups (appd-modedp): states are tuples on the
        cut (modes with one end processed and one not); the step for group g enumerates the state tuples together
        with every tuple on the new modes between g and later groups, skipping combinations of weight >= p,
        multiplies by the table entry of g and projects out the modes whose both ends are processed.
        Return (series, max_enum, self.child_calls) with self.child_calls reset to 0 first; max_enum is the largest
        number, over the steps, of (state, new tuple) combinations of total weight < p that the step multiplies."""
        raise NotImplementedError("T15: implement sum_dp")

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
                f *= sum(c * val[(w2, "c")] for c, w2 in self.modes[k]["lin"]) ** e
        return [(f * c) % m for c in ser]


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
    """appd-leafvars: V_g[gamma] of a chunk network computed by summing over the variables of the chunks of group g
    (net.units, net.unit_of, net.desig, net.clauses; wire -> (clause, variable) in _WIREMAP). Fixing the variables
    fixes every column wire of the group; each clause designated to a chunk of g then depends only on its external
    row bits: constant 1 if one of its internal literals is true, otherwise 1 - prod_j [external literal j false],
    two product terms. After expanding these terms every row bit y carries a unary weight (w0, w1), and
    sum_y w(y) l_y(form) = w0 + (w1 - w0) form enters that row's packet factor. Must equal ne.group_table(g)[key].
    Return (series, number of variables summed, max number of clauses split into two terms)."""
    raise NotImplementedError("T15: implement byvars_value")


IMPLEMENTED = False      # set to True once every function below works


def node_value(net, pk, labels, t, p, q, weights="paper", method="dp"):
    """F(t + alpha z) mod (2^p, z^q) from the mode expansion at node t (v2(t) < pk.a): method "paper", "memo"
    or "dp" selects NodeExpansion.sum_paper / sum_memo / sum_dp. Return (series, stats) where stats is
    {"tuples" or "max_enum": int, "child_calls": int, "modes": len(ne.modes)}."""
    raise NotImplementedError("T15: implement node_value")


def count_e2e(n, clauses, chunks=True, Lc=2, r=3, seed=0):
    """Exact #F by the whole algorithm with a = 1 (L = 2), u = 1: build the network (chunk super-vertices when
    chunks, units from make_units(n, clauses, Lc, random.Random(seed))), labels tuple(rng.randrange(r)) per vertex
    from the same rng, packing, h = M(2L - 1) // L + 1, multiplicities(1, 1, h); at every node t with q_t > 0 use
    precision P_t = p + loss_bits(...) with p = n + 1 and the 'sharp' weights, F(t + alpha z) mod (2^P_t, z^q_t)
    from sum_dp, then recover(...) and reduce mod 2^p; multiply by 2^(unused variables). Must equal count_sat."""
    raise NotImplementedError("T15: implement count_e2e")
