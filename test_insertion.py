from words import (
    Word,
    HopfPermutation,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
)
from crystals import (
    ShiftedCrystalGenerator,
    FPFCrystalGenerator
)

n = 4


def test_fpf_insertion():
    for w in get_fpf_involution_words((6, 5, 4, 3, 2, 1)):
        p, q, = Word(*w).fpf_insert()
        assert w == Tableau.inverse_fpf(p, q)


def test_involution_p_tableaux():
    k = 4
    for w in HopfPermutation.involutions(n):
        cg = ShiftedCrystalGenerator(w.oneline, k)
        shapes = [
            {cg.insertion_tableau(i) for i in comp}
            for comp in cg.components
        ]
        assert all(len(sh) == 1 for sh in shapes)
        assert all(len(s & t) == 0 for s in shapes for t in shapes if s != t)


def test_fpf_p_tableaux():
    k = 4
    for w in HopfPermutation.fpf_involutions(n):
        cg = FPFCrystalGenerator(w.oneline, k)
        shapes = [
            {cg.insertion_tableau(i) for i in comp}
            for comp in cg.components
        ]
        assert all(len(sh) == 1 for sh in shapes)
        assert all(len(s & t) == 0 for s in shapes for t in shapes if s != t)


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
