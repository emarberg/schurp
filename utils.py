from tableaux import Tableau


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
