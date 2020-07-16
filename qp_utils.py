from tableaux import Tableau
from permutations import Permutation
from signed import SignedPermutation


def rsk(w):
    p, q = {}, {}
    for a in w:
        i = 1
        while True:
            j = 1
            while (i, j) in p and p[(i, j)] <= a:
                j += 1
            if (i, j) not in p:
                p[(i, j)] = a
                q[(i, j)] = len(p)
                break
            else:
                a, p[(i, j)] = p[(i, j)], a
            i += 1
    ans = Tableau(p), Tableau(q)
    assert ans[0].is_standard() and ans[1].is_standard()
    return ans


def truncate_a(w, n):
    return tuple(i + 1 if a > n else a for i, a in enumerate(w[:n]))


def truncate_bc(w, n):
    return tuple(i + 1 if a > n else -i - 1 if a < -n else a for i, a in enumerate(w[:n]))


def gelfand_rsk(v, n, sgn):
    assert all(v[v[i - 1] - 1] == i for i in range(1, len(v) + 1))

    w = truncate_a(v, n)
    cycles = sorted([(w[i - 1], i) for i in range(1, len(w) + 1) if i < w[i - 1]])
    if sgn:
        cycles += reversed(sorted([(i, i) for i in range(1, n + 1) if i == w[i - 1]]))
    else:
        cycles += sorted([(i, i) for i in range(1, n + 1) if i == w[i - 1]])
    p = {}
    for cyc in cycles:
        _, a = cyc
        i = 1
        while True:
            j = 1
            while (i, j) in p and p[(i, j)] <= a:
                j += 1
            if (i, j) not in p:
                p[(i, j)] = a
                if cyc[0] != cyc[1] and sgn:
                    while i > 1 and (i - 1, j + 1) not in p:
                        i = i - 1
                    p[(i, j + 1)] = cyc[0]
                elif cyc[0] != cyc[1] and not sgn:
                    while j > 1 and (i + 1, j - 1) not in p:
                        j = j - 1
                    p[(i + 1, j)] = cyc[0]
                break
            else:
                a, p[(i, j)] = p[(i, j)], a
            i += 1
    ans = Tableau(p)
    assert ans.is_standard()
    return ans


def gelfand_a_conjugate(ht, w, i):
    i += 1
    if w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if w(i) == i + 1 and w(i + 1) == i:
        return ht, (w, False)
    s = Permutation.s_i(i)
    return ht + 1, s * w * s


def gelfand_a_start(n, k):
    return Permutation(*[1 + i + (-1)**i for i in range(2 * k)])


def gelfand_a_printer(n, k):
    w = gelfand_a_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_a_conjugate(ht, v, i)
        return ht, v

    return printer


def gelfand_bc_conjugate(ht, w, i):
    s = SignedPermutation.s_i(i, w.rank)
    if i == 0 and w(1) in [-1, 1]:
        return ht + 1, w * s
    if i > 0 and abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
        return ht, (w, False)
    if i > 0 and w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if i > 0 and w(i) == -i and w(i + 1) == -i - 1:
        return ht, (w, True)
    return ht + 1, s * w * s


def gelfand_bc_start(n, k):
    return SignedPermutation(*(
        [1 + i + (-1)**i for i in range(2 * k)] +
        [i for i in range(2 * k + 1, n + 1)]
    ))


def gelfand_bc_printer(n, k):
    w = gelfand_bc_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_bc_conjugate(ht, v, i)
        return ht, v

    return printer


def gelfand_d_conjugate(ht, w, i):
    if i == 0:
        s = SignedPermutation.ds_i(-1, w.rank)
        t = SignedPermutation.ds_i(1, w.rank)
        if abs(w(1)) != 1 and abs(w(2)) != 2:
            if w * s == s * w:
                return ht, (w, False)
            return ht + 1, s * w * s
        if (w(1) == 1 and w(2) == 2) or (w(1) == -1 and w(2) == -2):
            return ht + 1, s * w * t
        if (w(1) == 1 and w(2) == -2) or (w(1) == -1 and w(2) == 2):
            return ht, (w, True)
        if (abs(w(1)) == 1 and abs(w(2)) != 2) or (abs(w(1)) != 1 and abs(w(2)) == 2):
            return ht + 1, (s * w * s).dstar()
        raise Exception

    if abs(w(i)) == i + 1 and abs(w(i + 1)) == i:
        return ht, (w, False)
    if w(i) == i and w(i + 1) == i + 1:
        return ht, (w, True)
    if w(i) == -i and w(i + 1) == -i - 1:
        return ht, (w, True)

    s = SignedPermutation.ds_i(i, w.rank)
    return ht + 1, s * w * s


def gelfand_d_start(n, k):
    return SignedPermutation(*(
        [1 + i + (-1)**i for i in range(2 * k)] +
        [i for i in range(2 * k + 1, n + 1)]
    ))


def gelfand_d_printer(n, k):
    w = gelfand_d_start(n, k)

    def printer(word):
        ht, v = 0, w
        for i in word:
            ht, v = gelfand_d_conjugate(ht, v, i)
        return ht, v

    return printer
