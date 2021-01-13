from words import (
    Word,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
    involution_insert
)
from marked import MarkedNumber
from hopf import HopfPermutation
from crystals import (
    OrthogonalCrystalGenerator,
    # SymplecticCrystalGenerator
)
from permutations import Permutation


def test_primed_crystals(n=4, k=2):
    for pi in Permutation.involutions(n):
        for w in pi.get_primed_involution_words():
            for f in Word.increasing_factorizations(w, k):
                f = tuple(tuple(_) for _ in f)
                p1, q1 = involution_insert(*f)

                for i in range(-1, k):
                    g = Word.incr_crystal_f(f, i)
                    assert g is None or Word.incr_crystal_e(g, i) == f

                    try:
                        if g is not None:
                            p2, q2 = involution_insert(*g)
                            assert p1 == p2
                            assert q1.shifted_crystal_f(i) == q2
                        else:
                            assert q1.shifted_crystal_f(i) is None
                    except AssertionError:
                        print(f, '-- %s -->' % i, g)
                        print('insertion tableaux:')
                        print(p1)
                        print(p2)
                        print('recording tableaux:')
                        print(q1)
                        print()
                        print(q2)
                        print()
                        print(q1.shifted_crystal_f(i, True))
                        print()
                        input('')
                        assert p1 == p2
                    # g = Word.incr_crystal_e(f, i)

                    # if g is not None:
                    #     p2, q2 = involution_insert(*g)
                    #     assert p1 == p2
                    #     assert q1.shifted_crystal_e(i) == q2
                    # else:
                    #     assert q1.shifted_crystal_e(i) is None


def invert_fac(f):
    w = []
    for i, a in enumerate(f):
        for j in a:
            while abs(j) >= len(w):
                w += [0]
            w[abs(j)] = i + 1 if j > 0 else -i - 1
    return Word(*w[1:])


def test_mixed_insert():
    w = Word(3, 3, 2, 3, 3, 2, 1, 2, 3)
    p = Tableau.from_string("1,2',2,3',3;2,3',3;3").shift()
    q = Tableau.from_string("1,2,4,5,9;3,6,8;7").shift()
    assert w.mixed_insert() == (p, q)

    testset = {
        Permutation.from_involution_word(*w) for w in Permutation.all(5)
    }
    for pi in list(testset):
        for k in [2, 3, 4]:
            fac = [
                tuple([Word(*i) for i in tup])
                for x in pi.get_primed_involution_words()
                for tup in OrthogonalCrystalGenerator.get_increasing_primed_factorizations(x, k)
            ]

            for f in fac:
                w = invert_fac(f)
                p, q = involution_insert(*f)
                pp, qq = w.mixed_insert()
                print('pi =', pi)
                print('f =', f)
                print(p, q)
                print('w =', w)
                print(pp, qq)
                print()
                assert p == qq
                assert q == pp

    # print(len(testset))


def test_specific_primed_insertion():
    f = (Word(1, 4, 5), Word(3, -4), Word(2,))
    g = (Word(1, 4, 5), Word(3,), Word(2, -4))

    p = Tableau.from_string("1,2,4,5;3,5'").shift()
    q = Tableau.from_string("1,1,1,3';2,2").shift()
    r = Tableau.from_string("1,1,1,3';2,3").shift()

    i, j = involution_insert(*f)
    print(f)
    print(p)
    print(q)
    print(i)
    print(j)
    print()
    assert (p, q) == (i, j)
    assert (p, r) == involution_insert(*g)


def doperator(tab, index):
    def find(t, a):
        for i, j in t:
            if abs(t.get(i, j)) == a:
                return i, j

    def locations(t, x):
        a, b, c = None, None, None
        w = t.shifted_reading_word()
        for i in range(len(w)):
            a = i if abs(w[i]) == x else a
            b = i if abs(w[i]) == (x + 1) else b
            c = i if abs(w[i]) == (x + 2) else c
        return a, b, c

    if index == 0:
        i, j = find(tab, 2)
        return tab.set(i, j, -tab.get(i, j))

    a, b, c = locations(tab, index)
    if a < b < c or c < b < a:
        return tab

    i1, j1 = find(tab, index)
    i2, j2 = find(tab, index + 1)
    i3, j3 = find(tab, index + 2)

    x1, x3 = tab.get(i1, j1), tab.get(i3, j3)

    if i1 == j1 and i3 == j3 and x1.number * x3.number < 0:
        tab = tab.set(i1, j1, -x1).set(i3, j3, -x3)

    if b < a < c or c < a < b:
        i1, j1, i2, j2 = i2, j2, i3, j3

    x1, x2 = tab.get(i1, j1), tab.get(i2, j2)

    if i1 == i2 or j1 == j2:
        tab = tab.set(i1, j1, -x1 if i1 != j1 else x1)
        tab = tab.set(i2, j2, -x2 if i2 != j2 else x2)
    else:
        tab = tab.set(i1, j1, x1 - 1 if x1.number < 0 else x1 + 1)
        tab = tab.set(i2, j2, x2 + 1 if x2.number < 0 else x2 - 1)

    return tab


def test_primed_ck(bound=5):
    records = {}

    def check(v, w, i):
        v = tuple(Word(a) for a in v)
        w = tuple(Word(a) for a in w)
        p, q = involution_insert(*v)
        _, r = involution_insert(*w)
        records[i] = records.get(i, []) + [(v, w, q, r)]
        return involution_insert(*w)[0] == involution_insert(*v)[0]

    for n in range(bound):
        for pi in Permutation.involutions(n):
            for w in pi.get_primed_involution_words():
                if len(w) < 2:
                    continue
                a, b = w[:2]
                if (a < 0 and b > 0) or (a > 0 and b < 0):
                    v = (-b, -a) + w[2:]
                else:
                    v = (b, a) + w[2:]
                check(v, w, 0)
                for i in range(len(w) - 2):
                    a, b, c = w[i:i + 3]
                    v = w
                    if a == c:
                        v = w[:i] + (b, a, b) + w[i + 3:]
                    elif -a == c > 0:
                        v = w[:i] + (b, c, -b) + w[i + 3:]
                    elif -a == c < 0:
                        v = w[:i] + (-b, a, b) + w[i + 3:]
                    elif abs(a) < abs(c) < abs(b) or abs(b) < abs(c) < abs(a):
                        v = w[:i] + (b, a, c) + w[i + 3:]
                    elif abs(b) < abs(a) < abs(c) or abs(c) < abs(a) < abs(b):
                        v = w[:i] + (a, c, b) + w[i + 3:]
                    if v == w:
                        continue
                    if not check(v, w, i + 1):
                        print(v)
                        print(involution_insert(Word(*v))[0])
                        print(w)
                        print(involution_insert(Word(*w))[0])
                        raise Exception
            for i in records:
                for v, w, q, r in records[i]:
                    if i > 0:
                        if doperator(q, i) != r:
                            print('\n\n\n')
                            print(q)
                            print(r)
                            print(i, ':', q.shifted_reading_word())
                            print(doperator(q, i))
                            input('\n?')
            records = {}


def test_primed_insertion(bound=5):
    for n in range(bound):
        for k in range(bound):
            seen = {}
            for pi in Permutation.involutions(n):
                fac = [
                    tuple([Word(*i) for i in tup])
                    for w in pi.get_primed_involution_words()
                    for tup in OrthogonalCrystalGenerator.get_increasing_primed_factorizations(w, k)
                ]
                for f in fac:
                    print('f =', f)
                    p, q = involution_insert(*f)
                    # print('\n', p, q)
                    # print(seen.get((p, q), None))
                    assert (p, q) not in seen
                    seen[(p, q)] = f
                    # print()


n = 4


def test_fpf_insertion():
    for w in get_fpf_involution_words((6, 5, 4, 3, 2, 1)):
        p, q, = Word(*w).fpf_insert()
        assert w == Tableau.inverse_fpf(p, q)


def test_involution_p_tableaux():
    k = 4
    for w in Permutation.involutions(n):
        for cg in OrthogonalCrystalGenerator.from_permutation(n, w, k):
            shapes = [
                {cg.insertion_tableau(i) for i in comp}
                for comp in cg.components
            ]
            assert all(len(sh) == 1 for sh in shapes)
            assert all(len(s & t) == 0 for s in shapes for t in shapes if s != t)


# def test_fpf_p_tableaux():
#     k = 4
#     for w in HopfPermutation.fpf_involutions(n):
#         cg = SymplecticCrystalGenerator(w.oneline, k)
#         shapes = [
#             {cg.insertion_tableau(i) for i in comp}
#             for comp in cg.components
#         ]
#         assert all(len(sh) == 1 for sh in shapes)
#         assert all(len(s & t) == 0 for s in shapes for t in shapes if s != t)


def test_involution_insertion():
    a = list(HopfPermutation.involutions(n))
    for w in a:
        for e in get_involution_words(w.oneline):
            print(e)
            assert Word(*e).shifted_hecke_insert() == Word(*e).involution_insert()


def test_shifted_hecke_insert():
    w = Word(3)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 3})
    assert q == Tableau({(1, 1): 1})

    w = Word(3, 5)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 3, (1, 2): 5})
    assert q == Tableau({(1, 1): 1, (1, 2): 2})

    w = Word(3, 5, 4)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 3, (1, 2): 4, (2, 2): 5})
    assert q == Tableau({(1, 1): 1, (1, 2): 2, (2, 2): 3})

    w = Word(3, 5, 4, 1)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 1, (1, 2): 3, (1, 3): 4, (2, 2): 5})
    assert q == Tableau({(1, 1): 1, (1, 2): 2, (2, 2): 3, (1, 3): -4})

    w = Word(3, 5, 4, 1, 2)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 1, (1, 2): 2, (1, 3): 4, (2, 2): 3, (2, 3): 5})
    assert q == Tableau({(1, 1): 1, (1, 2): 2, (2, 2): 3, (1, 3): -4, (2, 3): -5})

    w = Word(3, 5, 4, 1, 2, 3)
    p, q = w.shifted_hecke_insert()
    assert p == Tableau({(1, 1): 1, (1, 2): 2, (1, 3): 3, (2, 2): 3, (2, 3): 4, (3, 3): 5})
    assert q == Tableau({(1, 1): 1, (1, 2): 2, (2, 2): 3, (1, 3): -4, (2, 3): -5, (3, 3): 6})
