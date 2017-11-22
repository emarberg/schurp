from tableaux import Tableau
from partitions import Partition


def test_repr():
    t = Tableau.from_string(" 1', 2',4;5 , 10,11' ;3,4'")
    assert str(t) == "1'  2'  4  \n5   10  11'\n3   4'     "


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
