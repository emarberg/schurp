from tableaux import Tableau
from marked import MarkedNumber
from partitions import Partition, StrictPartition
from words import Word


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
