from tableaux import ShiftedTableau


def test_repr():
    t = ShiftedTableau(string=" 1', 2',4,5 , 10,11' ;3,4'")
    assert str(t) == "1'  2'  4   5   10  11' \n" + \
                     "    3   4'  "
