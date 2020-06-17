from words import (
    Word,
    HopfPermutation,
    get_involution_words,
    get_fpf_involution_words,
    Tableau,
    involution_insert
)
from crystals import (
    OrthogonalCrystalGenerator,
    SymplecticCrystalGenerator
)
from permutations import Permutation


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


def test_primed_insertion():
    for n in range(5):
        for k in range(5):
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
                    print('\n', p, q)
                    print(seen.get((p, q), None))
                    assert (p, q) not in seen
                    seen[(p, q)] = f
                    print()


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
