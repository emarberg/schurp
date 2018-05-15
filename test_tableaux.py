from tableaux import Tableau
from partitions import Partition, StrictPartition


# def test_repr():
#     t = Tableau.from_string(" 1', 2',4;5 , 10,11' ;3,4'")
#     assert str(t) == "1'  2'  4  \n5   10  11'\n3   4'     "


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
    assert Tableau.get_semistandard_shifted(p) == {
        Tableau.from_string("1,1;2").shift(),
        Tableau.from_string("1,2;3").shift(),
        Tableau.from_string("1,2';2").shift(),
        Tableau.from_string("1,2';3").shift(),
    }


def test_toggle():
    p = StrictPartition(4, 2, 1)
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
