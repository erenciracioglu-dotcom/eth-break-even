"""T23 self-check (round 10, appd-bincross part 1): bc_mask, bc_left and bc_right of appd_core against dense
Kronecker products of the 4 x 4 matrices of NOTES.md, section 'Round 10', for m = 0..3, both pairings, random
series vectors with some None entries. Prints T23 OK when everything matches."""
import itertools
import random
import appd_core as C

S_ROW = [[0, 1, 0, 0], [1, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 1]]           # S (row pairing)
KT_COL = [[1, 0, 0, 1], [0, 0, -1, 0], [0, -1, 0, 0], [1, 0, 0, 0]]         # K T (column pairing)
N0_ROW = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]]           # S^-1
N1_ROW = [[0, -1, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]          # -S^-1 J
N0_COL = [[1, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]           # T^-1 P0
N1_COL = [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [1, 0, 0, -1]]       # -T^-1


def kron_entry(mats, row, col, m):
    """Entry (row, col) of mats[0] (x) mats[1] (x) ... (x) mats[m-1]; digit 0 is the most significant."""
    v = 1
    for k in range(m):
        sh = 2 * (m - 1 - k)
        v *= mats[k][(row >> sh) & 3][(col >> sh) & 3]
        if v == 0:
            return 0
    return v


def dense_apply(mats, vec, m, q, mod):
    out = []
    for row in range(4 ** m):
        acc = [0] * q
        for col in range(4 ** m):
            c = kron_entry(mats, row, col, m)
            if c and vec[col] is not None:
                acc = [(x + c * y) % mod for x, y in zip(acc, vec[col])]
        out.append(acc)
    return out


def main():
    rng = random.Random(2323)
    bad = 0
    for m in range(4):
        for pairing in ("row", "col"):
            M = S_ROW if pairing == "row" else KT_COL
            MT = [[M[d][d2] for d in range(4)] for d2 in range(4)]
            N0, N1 = (N0_ROW, N1_ROW) if pairing == "row" else (N0_COL, N1_COL)
            for trial in range(3):
                q, mod = rng.choice([1, 2, 3]), 1 << rng.choice([3, 6, 9])
                vec = [None if rng.random() < 0.2 else [rng.randrange(mod) for _ in range(q)] for _ in range(4 ** m)]
                want_mask = [[w for w in range(4 ** m)
                              if sum(1 for k in range(m) if (w >> (2 * k)) & 3 == 0) >= (i if pairing == "row" else m - i)]
                             for i in range(m + 1)]
                got_mask = [list(C.bc_mask(m, i, pairing)) for i in range(m + 1)]
                want_left = dense_apply([MT] * m, vec, m, q, mod)
                got_left = [[int(x) % mod for x in s] for s in C.bc_left([None if v is None else list(v) for v in vec],
                                                                          m, pairing, q, mod)]
                want_right = [[[0] * q for _ in range(4 ** m)] for _ in range(m + 1)]
                for sel in itertools.product((0, 1), repeat=m):
                    part = dense_apply([N1 if s else N0 for s in sel], vec, m, q, mod)
                    i = sum(sel)
                    want_right[i] = [[(x + y) % mod for x, y in zip(a, b)] for a, b in zip(want_right[i], part)]
                got_right = [[[int(x) % mod for x in s] for s in piece]
                             for piece in C.bc_right([None if v is None else list(v) for v in vec], m, pairing, q, mod)]
                ok = got_mask == want_mask and got_left == want_left and got_right == want_right
                if not ok:
                    bad += 1
                    print("MISMATCH m", m, pairing, "mask", got_mask == want_mask, "left", got_left == want_left,
                          "right", got_right == want_right)
    # the pieces of the crossing: sum_i sum_{w in mask_i} L[w] R[i][w] = f^T SWAP^(m) g
    for m in range(1, 4):
        for pairing in ("row", "col"):
            mod, q = 1 << 9, 1
            f = [[rng.randrange(mod)] for _ in range(4 ** m)]
            g = [[rng.randrange(mod)] for _ in range(4 ** m)]
            L = C.bc_left(f, m, pairing, q, mod)
            R = C.bc_right(g, m, pairing, q, mod)
            tot = sum(int(L[w][0]) * int(R[i][w][0]) for i in range(m + 1) for w in C.bc_mask(m, i, pairing)) % mod
            # row: words (x, y) -> (u, v) = (y, x) digit by digit; column: f on (x, u), g on (y, v), [u = y][v = x]
            want = 0
            for wf in range(4 ** m):
                wg = 0
                for k in range(m):
                    d = (wf >> (2 * (m - 1 - k))) & 3
                    hi, lo = d >> 1, d & 1
                    wg = 4 * wg + (2 * lo + hi if pairing == "row" else 2 * lo + hi)
                want += f[wf][0] * g[wg][0]
            if tot != want % mod:
                bad += 1
                print("PIECES DO NOT SUM TO THE CROSSING m", m, pairing)
    print("T23 OK" if bad == 0 else "T23 FAIL (%d mismatches)" % bad)


if __name__ == "__main__":
    main()
