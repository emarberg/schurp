from tableaux import Tableau
from marked import MarkedNumber
from partitions import Partition, StrictPartition
from words import Word, decomposition_insert, sk_insert
from crystals import AbstractQCrystal


def decomposition_monoid_maps(n=3, k=10, mu=None):
    partitions = [mu] if mu is not None else sorted({mu.transpose().tuple() for mu in Partition.all(k, max_part=n) if mu.transpose().is_strict()})
    print(partitions)
    print()

    for mu in partitions:
        b = AbstractQCrystal.decomposition_tableaux_from_strict_partition(mu, n)
        print('n =', n, 'mu =', mu, ':', len(b))
        for t in b:
            u = sk_insert(tuple(reversed(t.row_reading_word())))[0]
            v = decomposition_insert(tuple(reversed(u.row_reading_word())))[0]
            print(t.tex())
            print('& \\leftrightarrow ')
            print(u.tex())
            print('\\\\')
            assert t == v


def test_decomposition_insertion():
    t = Tableau()
    v = 2
    expected = Tableau.shifted_from_rows([[2]])
    assert t.decomposition_insert(v) == expected

    t = Tableau.shifted_from_rows([[6, 6, 1, 3, 5]])
    v = 2
    expected = Tableau.shifted_from_rows([[6, 6, 3, 2, 5], [1]])
    assert t.decomposition_insert(v) == expected

    t = Tableau.shifted_from_rows([[3, 2, 4]])
    v = 1
    expected = Tableau.shifted_from_rows([[4, 2, 1], [3]])
    assert t.decomposition_insert(v) == expected

    t = Tableau.shifted_from_rows([[6, 6, 1, 3, 5], [3, 2, 4]])
    v = 2
    expected = Tableau.shifted_from_rows([[6, 6, 3, 2, 5], [4, 2, 1], [3]])
    assert t.decomposition_insert(v) == expected

    p = Tableau.shifted_from_rows([[3, 2, 1], [2]])
    q = Tableau.shifted_from_rows([[1, 2, 4], [3]])
    assert (p, q) == decomposition_insert(2, 3, 2, 1)

    p = Tableau.shifted_from_rows([[4, 4, 4, 4, 4, 4], [3, 3, 3, 3], [2, 2], [1]])
    q = Tableau.shifted_from_rows([[1, 2, 4, 7, 8, 13], [3, 5, 9, 12], [6, 10], [11]])
    assert (p, q) == decomposition_insert(1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4)


def test_standardize(n=6):
    success = 0
    for mu in StrictPartition.all(n):
        mu = tuple(mu)
        for nu in StrictPartition.subpartitions(mu, strict=True):
            for t in Tableau.standard_shifted_marked(mu, nu, True):
                s = t.destandardize()
                assert s.standardize() == t
                success += 1
    assert success > 0


def test_shifted_crystal_word():
    tab = Tableau.from_string("1,1,4',4;,2,4',5';,,4,5")
    word, positions = tab.shifted_crystal_word()
    assert word == [MarkedNumber(i) for i in [-5, -4, -4, 4, 5, 2, 1, 1, 4]]
    assert all(tab[box] == word[i] for i, box in enumerate(positions))


def test_standard():
    p = Partition(2, 1)
    assert Tableau.get_standard(p) == {
        Tableau.from_string("1,2;3"),
        Tableau.from_string("1,3;2")
    }


def test_semistandard():
    p = Partition(2, 1)
    assert Tableau.get_semistandard(p) == {
        Tableau.from_string("1,1;2"),
        Tableau.from_string("1,2;2"),
        Tableau.from_string("1,3;2"),
        Tableau.from_string("1,2;3"),
    }


def test_standard_shifted():
    p = StrictPartition(2, 1)
    assert Tableau.get_standard_shifted(p) == {
        Tableau.from_string("1,2;3").shift(),
        Tableau.from_string("1,2';3").shift()
    }


def test_semistandard_shifted():
    p = StrictPartition(2, 1)
    assert Tableau.get_semistandard_shifted(p, n=3) == {
        Tableau.from_string("1,1;2").shift(),
        Tableau.from_string("1,1;3").shift(),
        Tableau.from_string("1,2;3").shift(),
        Tableau.from_string("1,3';3").shift(),
        Tableau.from_string("1,2';2").shift(),
        Tableau.from_string("1,2';3").shift(),
        Tableau.from_string("2,2;3").shift(),
        Tableau.from_string("2,3';3").shift(),
    }


def test_toggle():
    p = StrictPartition(4, 1)
    for i, t in enumerate(Tableau.get_semistandard_shifted(p)):
        assert t.toggle().toggle() == t


def test_bump():
    p = 4
    column_dir = False
    seq = ()
    assert Tableau.bump(p, column_dir, seq) == (None, False, (4,))

    seq = (5, 6, 7)
    assert Tableau.bump(p, column_dir, seq) == (5, True, (4, 6, 7))

    seq = (4, 6, 7)
    assert Tableau.bump(p, column_dir, seq) == (6, True, (4, 6, 7))

    seq = (3, 6, 7)
    assert Tableau.bump(p, column_dir, seq) == (6, False, (3, 4, 7))

    seq = (3, 4, 7)
    assert Tableau.bump(p, column_dir, seq) == (7, False, (3, 4, 7))


def test_inverse_inv():
    for w in [
        (5, 3, 4, 2, 3, 2, 1, 2, 3),
        (5, 4, 3, 4, 2, 3, 1, 2, 3),
        (1, 5, 3, 2, 4, 3, 1, 2, 1),
        (1, 3, 2, 5, 4, 1, 3, 2, 1),
        (3, 1, 2, 5, 1, 4, 3, 2, 1),
        (3, 5, 4, 2, 1, 3, 2, 1, 3),
        (3, 5, 1, 4, 2, 3, 1, 2, 1),
        (1, 5, 3, 2, 1, 4, 3, 2, 1),
        (5, 4, 3, 1, 2, 4, 1, 3, 2),
        (5, 1, 3, 2, 4, 3, 1, 2, 1),
        (1, 3, 5, 2, 1, 4, 3, 2, 1),
        (5, 3, 1, 2, 1, 4, 3, 2, 1),
        (3, 1, 5, 4, 2, 3, 1, 2, 1),
        (1, 5, 3, 2, 4, 3, 2, 1, 2),
        (5, 1, 3, 2, 1, 4, 3, 2, 1),
        (5, 4, 3, 1, 2, 1, 4, 3, 2),
        (3, 5, 1, 2, 4, 1, 3, 2, 1),
        (1, 3, 5, 4, 3, 2, 1, 3, 2),
        (3, 5, 2, 4, 3, 1, 2, 3, 1),
        (3, 5, 1, 4, 2, 3, 2, 1, 2),
        (3, 5, 4, 1, 2, 3, 1, 2, 1),
        (3, 5, 1, 4, 3, 2, 1, 3, 2),
        (5, 3, 1, 4, 3, 2, 1, 3, 2),
        (5, 1, 3, 2, 4, 3, 2, 1, 2),
        (3, 1, 5, 2, 4, 1, 3, 2, 1),
        (3, 1, 5, 4, 2, 3, 2, 1, 2),
        (3, 1, 2, 1, 5, 4, 3, 2, 1),
        (5, 4, 3, 4, 1, 2, 1, 3, 2),
        (1, 5, 3, 4, 3, 2, 1, 3, 2),
        (5, 3, 2, 4, 3, 2, 1, 2, 3),
        (5, 3, 4, 1, 3, 2, 1, 3, 2),
        (3, 5, 4, 1, 2, 3, 2, 1, 2),
        (5, 3, 4, 2, 1, 3, 2, 3, 1),
        (3, 5, 2, 4, 1, 3, 2, 1, 3),
        (5, 1, 3, 4, 3, 2, 1, 3, 2),
        (3, 5, 4, 3, 2, 1, 3, 2, 3),
        (3, 2, 1, 5, 4, 3, 2, 1, 3),
        (5, 3, 4, 2, 1, 2, 3, 2, 1),
        (5, 1, 4, 3, 2, 1, 4, 3, 2),
        (1, 5, 4, 3, 4, 2, 1, 3, 2),
        (5, 4, 3, 2, 1, 4, 2, 3, 2),
        (5, 4, 1, 3, 2, 4, 1, 3, 2),
        (5, 4, 1, 3, 4, 2, 3, 1, 2),
        (3, 2, 5, 4, 1, 3, 2, 3, 1),
        (3, 5, 4, 3, 1, 2, 3, 1, 2),
        (5, 1, 4, 3, 4, 2, 1, 3, 2),
        (1, 3, 2, 5, 4, 3, 1, 2, 1),
        (3, 5, 2, 1, 4, 3, 2, 1, 3),
        (5, 4, 3, 2, 4, 1, 2, 3, 2),
        (5, 4, 3, 2, 1, 2, 4, 3, 2),
        (5, 4, 3, 1, 2, 4, 3, 1, 2),
        (5, 4, 3, 4, 2, 1, 2, 3, 2),
        (1, 5, 4, 3, 2, 4, 3, 1, 2),
        (5, 3, 2, 4, 1, 3, 2, 3, 1),
        (1, 3, 2, 5, 4, 3, 2, 1, 2),
        (3, 2, 5, 1, 4, 3, 2, 3, 1),
        (3, 5, 1, 2, 4, 3, 1, 2, 1),
        (5, 3, 2, 4, 1, 2, 3, 2, 1),
        (5, 3, 1, 4, 3, 2, 3, 1, 2),
        (5, 1, 4, 3, 2, 4, 3, 1, 2),
        (1, 3, 5, 4, 3, 2, 3, 1, 2),
        (5, 3, 2, 1, 2, 4, 3, 2, 1),
        (3, 1, 5, 2, 4, 3, 1, 2, 1),
        (5, 4, 3, 4, 1, 2, 3, 1, 2),
        (5, 3, 1, 4, 2, 1, 3, 2, 1),
        (1, 3, 5, 4, 2, 1, 3, 2, 1),
        (3, 5, 1, 2, 4, 3, 2, 1, 2),
        (5, 3, 4, 3, 2, 1, 2, 3, 2),
        (1, 5, 3, 4, 3, 2, 3, 1, 2),
        (3, 2, 5, 4, 1, 2, 3, 2, 1),
        (5, 3, 2, 1, 4, 3, 2, 3, 1),
        (5, 4, 3, 1, 4, 2, 1, 3, 2),
        (5, 3, 4, 1, 3, 2, 3, 1, 2),
        (3, 1, 5, 2, 4, 3, 2, 1, 2),
        (3, 5, 4, 2, 3, 1, 2, 1, 3),
        (5, 3, 2, 1, 4, 2, 3, 2, 1),
        (1, 5, 3, 4, 2, 1, 3, 2, 1),
        (5, 1, 3, 4, 3, 2, 3, 1, 2),
        (5, 3, 4, 1, 2, 1, 3, 2, 1),
        (5, 1, 3, 4, 2, 1, 3, 2, 1),
        (3, 2, 5, 1, 4, 2, 3, 2, 1),
        (5, 3, 4, 3, 2, 3, 1, 2, 3),
        (5, 3, 4, 2, 3, 1, 2, 3, 1),
        (3, 5, 2, 4, 3, 1, 2, 1, 3),
        (5, 3, 4, 3, 1, 2, 1, 3, 2),
        (1, 5, 4, 3, 2, 4, 1, 3, 2),
        (1, 3, 2, 5, 1, 4, 3, 2, 1),
        (3, 5, 4, 2, 3, 2, 1, 2, 3),
        (3, 2, 5, 4, 3, 2, 1, 2, 3),
        (3, 2, 5, 4, 3, 1, 2, 1, 3),
        (5, 1, 4, 3, 2, 4, 1, 3, 2),
        (5, 3, 4, 2, 1, 3, 2, 1, 3),
        (1, 3, 5, 4, 2, 3, 1, 2, 1),
        (5, 3, 1, 4, 2, 3, 1, 2, 1),
        (3, 2, 1, 5, 4, 3, 2, 3, 1),
        (1, 5, 4, 3, 2, 1, 4, 3, 2),
        (3, 1, 2, 5, 4, 1, 3, 2, 1),
        (3, 5, 1, 2, 1, 4, 3, 2, 1),
        (5, 4, 3, 1, 4, 2, 3, 1, 2),
        (1, 3, 5, 2, 4, 1, 3, 2, 1),
        (5, 3, 1, 2, 4, 1, 3, 2, 1),
        (3, 1, 5, 2, 1, 4, 3, 2, 1),
        (1, 3, 2, 1, 5, 4, 3, 2, 1),
        (5, 3, 4, 1, 2, 3, 1, 2, 1),
        (1, 5, 3, 4, 2, 3, 1, 2, 1),
        (5, 3, 2, 4, 3, 1, 2, 3, 1),
        (1, 3, 5, 4, 2, 3, 2, 1, 2),
        (5, 1, 3, 4, 2, 3, 1, 2, 1),
        (5, 3, 1, 4, 2, 3, 2, 1, 2),
        (3, 1, 5, 4, 3, 2, 1, 3, 2),
        (1, 5, 3, 4, 2, 3, 2, 1, 2),
        (3, 5, 2, 4, 3, 2, 1, 2, 3),
    ]:
        p, q = Word(*w).involution_insert(verbose=True)
        assert Tableau.inverse_inv(p, q) == w


def test_from_composition():
    assert Tableau.from_composition(()) == Tableau()
    assert Tableau().weight() == ()

    alpha = (1, 0, 3, 2, 5, 0, 3)
    assert Tableau.from_composition(alpha) == Tableau({
        (1, 1): 1, (1, 2): 3, (1, 3): 3, (1, 4): 5, (1, 5): 5,
        (2, 1): 3, (2, 2): 4, (2, 3): 5,
        (3, 1): 4, (3, 2): 5, (3, 3): 7,
        (4, 1): 5, (4, 2): 7,
        (5, 1): 7
    })
    assert Tableau.from_composition(alpha).weight() == alpha
