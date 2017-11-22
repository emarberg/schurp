from tableaux import Tableau


def test_repr():
    t = Tableau.from_string(" 1', 2',4;5 , 10,11' ;3,4'")
    assert str(t) == "1'  2'  4  \n5   10  11'\n3   4'     "
