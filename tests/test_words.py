from permutations import Permutation
from words import Word, eg_insert
from vectors import Vector
from tableaux import Tableau


def test_drop():
    assert Word.drop_alignment((3, 6), (4, 7)) == [[3, 6, None], [None, 4, 7]]
    assert Word.drop_alignment((3, 4, 5), (2, 4)) == [[3, 4, 5], [2, None, 4]]

    assert Word.drop((1, 3), (2,)) == ((3,), (1, 2))
    assert Word.drop((1,2,4),(1,3)) == ((2, 4), (1, 2, 3))
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
