from partitions import Partition
from permutations import Permutation
from words import Word, eg_insert, weak_eg_insert
from vectors import Vector
from tableaux import Tableau


def pathtest(n):
    w = Permutation.random_reduced_word(n)
    last = None
    last_path = None
    for i in range(len(w)):
        p, q = eg_insert(*w[:i + 1])
        pp, q = weak_eg_insert(*w[:i + 1])
        print(p)
        print(pp)
        path = Word.lifting_path(p.row_reading_word())
        print('this:', path)
        if last is not None:
            box = [b for b in p.mapping if b not in last.mapping][0]
            print('last:', last_path)
            print() # print(box)
            print(list(last.partition()) + [0,])
            print(list(p.partition()))
            print()
        last = p
        last_path = path
        # print(p)
        input('')
    print(w)


def test_weak_eg_descent(n=5):
    for pi in Permutation.all(n):
        for w in pi.get_reduced_words():
            l = len(w)
            if l < 2:
                continue

            p, q = weak_eg_insert(*w)

            i1, j1, = [a for a, b in q.mapping.items() if abs(b) == 1][0]
            i2, j2, = [a for a, b in q.mapping.items() if abs(b) == 2][0]

            _p, _q = weak_eg_insert(*w)

            _i1, _j1, = [a for a, b in _q.mapping.items() if abs(b) == 1][0]
            _i2, _j2, = [a for a, b in _q.mapping.items() if abs(b) == 2][0]

            print()
            print(w)

            # q = Tableau({a: b if abs(b) in [1, 2] else 0 for a, b in q.mapping.items()})

            # p, q = eg_insert(*w)
            # q = Tableau({a: b for a, b in q.mapping.items() if abs(b) in [l - 1, l]})
            # print(p)
            # print(q)

            if w[-2] < w[-1]:
                assert j1 > j2
            #    assert i1 <= i2

            if w[-2] > w[-1]:
                assert j1 <= j2
            #    assert i1 > i2

            assert j1 == _j1
            assert j2 == _j2


def test_weak_eg_insertion(n=5):
    for flag in Partition.flags(n - 1):
        print('flag =', flag)
        for w in Permutation.all(n):
            for f in w.get_bounded_increasing_factorizations(flag=flag):
                p, q = weak_eg_insert(*f)
                if not q.is_key_flagged(flag=flag):
                    print(p)
                    print(q)
                    print(f)
                    print()
                    assert q.is_key_flagged(flag=flag)


def test_lift():
    assert Word.lift_alignment((6,), (3, 4, 7)) == [[None, 6, None], [3, 4, 7]]
    assert Word.lift_alignment((3, 5), (2, 4, 5)) == [[3, None, 5], [2, 4, 5]]

    assert Word.lift((3,), (1, 2)) == ((1, 3), (2,))
    assert Word.lift((2, 4), (1, 2, 3)) == ((1, 2, 4), (1, 3))
    assert Word.lift((6,), (3, 4, 7)) == ((3, 6, 7), (4,))
    assert Word.lift((3, 5), (2, 4, 5)) == ((3, 4, 5), (2, 4))

    assert Word.lift((6, 9, 3, 7, 8, 2, 3, 5, 9, 1, 2, 4, 5, 6)) == Tableau.from_string("1,2,4;2,3,4,5,6;3,8;6,7,8,9")


def test_drop():
    assert Word.drop_alignment((3, 6), (4, 7)) == [[3, 6, None], [None, 4, 7]]
    assert Word.drop_alignment((3, 4, 5), (2, 4)) == [[3, 4, 5], [2, None, 4]]

    assert Word.drop((1, 3), (2,)) == ((3,), (1, 2))
    assert Word.drop((1, 2, 4), (1, 3)) == ((2, 4), (1, 2, 3))
    assert Word.drop((3, 6), (4, 7)) == ((6,), (3, 4, 7))
    assert Word.drop((3, 4, 5), (2, 4)) == ((3, 5), (2, 4, 5))

    assert Word.drop((3, 6, 4, 7, 5, 2, 4)) == Tableau.from_string("2,4,5;3,5;6,7")


def test_drop_all(n=5):
    for g in Permutation.all(n):
        for w in g.get_reduced_words():
            tab = Word.drop(w)
            eg_tab = eg_insert(w)[0]
            if tab != eg_tab:
                print(w)
                print(tab)
                print(eg_tab)
            assert tab == eg_tab


def test_sums():
    s = {1, 2, 3}
    assert Word(3, 2, 1, subset=s) - Word(3, 2, 1, subset=s) == Vector()
    assert Word(3, 2, 1, subset=s) + Word(3, 2, 1, subset=s) == Vector({Word(3, 2, 1, subset=s): 2})


def test_shuffle():
    s = {1, 2, 3, 4}
    u = Word(2, 1, subset={1, 2})
    v = Word(4, 3, subset={3, 4})

    assert u * v == \
        Word(2, 1, 4, 3, subset=s) + Word(2, 4, 1, 3, subset=s) + Word(2, 4, 3, 1, subset=s) + \
        Word(4, 2, 1, 3, subset=s) + Word(4, 2, 3, 1, subset=s) + Word(4, 3, 2, 1, subset=s)

    assert u * Word() == Vector.base(u)
    assert Word() * v == Vector.base(v)


def test_coproduct():
    s = {1, 2, 3, 4}
    w = Word(2, 1, 2, 4, 4, subset=s)

    u = Word(2, 1, 2, subset={1, 2})
    v = Word(4, 4, subset={3, 4})
    assert w.coproduct({1, 2}, {3, 4}) == Vector.base((u, v))
    assert w.coproduct({3, 4}, {1, 2}) == Vector()

    h = Word(2, 1, 2, subset={1, 2, 3})
    x = Word(2, 1, 2, subset={1, 2})
    y = Word(subset={3})
    z = Word(4, 4, subset={4})
    assert w.coproduct({1, 2}, {3}, {4}) == Vector.base((x, y, z))
    assert w.coproduct({1, 2}, {4}, {3}) == Vector.base((x, z, y))
    assert w.coproduct({1, 2, 3}, {4}) == Vector.base((h, z))

    assert w.coproduct({1, 4}, {2, 3}) == Vector()
