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


def gelfand_rsk(w, sgn):
    assert all(w[w[i - 1] - 1] == i for i in range(1, len(w) + 1))
    cycles = sorted([(w[i - 1], i) for i in range(1, len(w) + 1) if i <= w[i - 1]])
    p = {}
    for b, a in cycles:
        if a == b and sgn:
            i = 1
            while (i, 1) in p:
                i += 1
            p[(i, 1)] = a
        elif a == b and not sgn:
            j = 1
            while (1, j) in p:
                j += 1
            p[(1, j)] = a
        else:
            i = 1
            while True:
                j = 1
                while (i, j) in p and p[(i, j)] <= a:
                    j += 1
                if (i, j) not in p:
                    p[(i, j)] = a
                    if sgn:
                        while i > 1 and (i - 1, j + 1) not in p:
                            i = i - 1
                        p[(i, j + 1)] = b
                    else:
                        while j > 1 and (i + 1, j - 1) not in p:
                            j = j - 1
                        p[(i + 1, j)] = b
                    break
                else:
                    a, p[(i, j)] = p[(i, j)], a
                i += 1
    return Tableau(p)
