from words import Word, Permutation
from vectors import Vector


def test_sums():
    assert Word(3, 2, 1) - Word(3, 2, 1) == Vector()
    assert Word(3, 2, 1) + Word(3, 2, 1) == Vector({Word(3, 2, 1): 2})


def test_shuffle():
    u = Word(2, 1)
    v = Word(4, 3)

    assert u * v == \
        Word(2, 1, 4, 3) + Word(2, 4, 1, 3) + Word(2, 4, 3, 1) + \
        Word(4, 2, 1, 3) + Word(4, 2, 3, 1) + Word(4, 3, 2, 1)

    assert u * u == 2 * Word(2, 1, 2, 1) + 4 * Word(2, 2, 1, 1)

    assert u * Word() == Vector.base(u)
    assert Word() * v == Vector.base(v)


def test_permutations():
    u = Permutation(1)
    v = Permutation(2, 1)
    w = Permutation(3, 2, 1)
    assert u * u == Vector.base(u)
    assert u * v == v * u == Vector.base(v)
    assert v * v == Permutation(2, 3, 1) + Permutation(3, 1, 2)
    assert w * v == Permutation(3, 2, 4, 1) + Permutation(3, 4, 1, 2) + Permutation(4, 2, 1, 3)

    assert Permutation(1, 2).oneline == (1, 2)
    assert Permutation(1, 2, 3).oneline == (1, 2, 3)


def test_exclude():
    u = Permutation(1, 2, 3)
    assert u.exclude(2) == Permutation(1, 2)

    u = Permutation(2, 1, 3)
    assert u.exclude(2) == Permutation(1, 2)

    u = Permutation(2, 1, 3, 5, 4, 6)
    assert u.exclude(5) == Permutation(2, 1, 3, 4, 5)
    assert u.exclude(6, 5, 3) == Permutation(2, 1, 3)

    u = Permutation(4, 1, 2, 3)
    assert u.exclude(2, 2) == Permutation(3, 1, 2)
    assert u.exclude(2, 3) == Permutation(2, 1)


def test_long_element_conjecture():
    def subset(m, n):
        if m == 0 or n == 0:
            if m == 1 or n == 1:
                return {Permutation()}
            else:
                return set()
        return set((Permutation.reverse(m) * Permutation.reverse(n)).keys())

    for i in range(1, 4):
        for j in range(1, 4):
            s = subset(i, j)
            a = set()
            b = set()
            for w in s:
                assert w(1) == i or w(i + j - 1) == i
                if w(1) == i:
                    a.add(w.exclude(i))
                if w(i + j - 1) == i:
                    b.add(w.exclude(i))
            print(s)
            print(a, subset(i - 1, j), i - 1, j)
            print(b, subset(i, j - 1), i, j - 1)
            assert a == subset(i - 1, j)
            assert b == subset(i, j - 1)


def test(k=4):
    for i in range(1, k):
        for j in range(1, k):
            for u in Permutation.all(i):
                for v in Permutation.all(j):
                    w = u * v

                    if u(i) == i:
                        oneline = u.oneline[:-1] + tuple(t + i - 1 for t in v.oneline)
                        q = Permutation(*oneline)
                        assert w == Vector.base(q)
                        continue

                    if v(1) == 1:
                        oneline = u.oneline + tuple(t + i - 1 for t in v.oneline[1:])
                        q = Permutation(*oneline)
                        assert w == Vector.base(q)
                        continue

                    a = u.find(i)
                    initial = tuple(u(t) for t in range(1, a + 1))

                    b = v.find(1)
                    final = tuple(v(t) + i - 1 for t in range(b, j + 1))
                    ufinal = tuple(v(t) for t in range(b, j + 1))

                    assert all(x.startswith(initial) or x.endswith(final) for x in w)

                    s, t = Vector(), Vector()
                    for x in w:
                        if x.startswith(initial):
                            s += x.exclude(*initial)
                        if x.endswith(final):
                            t += x.exclude(*final)
                    g, h = u.exclude(*initial) * v, u * v.exclude(*ufinal)

                    print(u, '*', v, '=', w)
                    print(initial, final)
                    print(s, ',', t)
                    print(g, ',', h)
                    print('')

                    assert g == s
                    assert h == t
