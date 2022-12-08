from stable.boxes import BoxOperator


def test_multiply():
    assert BoxOperator() * BoxOperator() == BoxOperator()
    assert BoxOperator(1) * BoxOperator(2) == BoxOperator(3)
    assert BoxOperator(1, 1) * BoxOperator(2, -2) == BoxOperator(3, 1, -2)
    assert BoxOperator(1, 1) * 'r2' == BoxOperator(1, 1, 2)
    assert BoxOperator(1, 1) * 'C2' == BoxOperator(1, 1, -2)
    assert BoxOperator(1, 1) * 0 == 0


def test_add():
    assert BoxOperator() == 0
    assert 0 + BoxOperator(1, 1) == BoxOperator(1, 1) + 0 == BoxOperator(1, 1)
    assert 0 - BoxOperator(1, 1) == -BoxOperator(1, 1) - 0 == -BoxOperator(1, 1)
