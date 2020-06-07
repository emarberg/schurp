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
    return ans


def truncate_a(w, n):
    return tuple(i + 1 if a > n else a for i, a in enumerate(w[:n]))


def gelfand_rsk(w, n, sgn):
    assert all(w[w[i - 1] - 1] == i for i in range(1, len(w) + 1))

    w = truncate_a(w, n)
    cycles = sorted([(w[i - 1], i) for i in range(1, len(w) + 1) if i < w[i - 1]])
    if sgn:
        cycles += reversed(sorted([(i, i) for i in range(1, n + 1) if i == w[i - 1]]))
    else:
        cycles += sorted([(i, i) for i in range(1, n + 1) if i == w[i - 1]])
    p = {}
    for b, a in cycles:
        i = 1
        while True:
            j = 1
            while (i, j) in p and p[(i, j)] <= a:
                j += 1
            if (i, j) not in p:
                p[(i, j)] = a
                if a != b and sgn:
                    while i > 1 and (i - 1, j + 1) not in p:
                        i = i - 1
                    p[(i, j + 1)] = b
                elif a != b and not sgn:
                    while j > 1 and (i + 1, j - 1) not in p:
                        j = j - 1
                    p[(i + 1, j)] = b
                break
            else:
                a, p[(i, j)] = p[(i, j)], a
            i += 1
    return Tableau(p)
