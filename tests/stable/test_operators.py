from stable.operators import (
    operator_shifted_add,
    operator_shifted_row_remove,
    operator_shifted_column_remove,
    operator_AL,
    operator_R,
    operator_AR,
    operator_C,
)
from stable.vectors import Vector
from stable.polynomials import X


beta = X(0)


def test_shifted_add():
    mu = ()
    op = operator_shifted_add(0)
    assert op(mu) == Vector({(1,): 1})

    mu = (1, )
    assert op(mu) == Vector({(1,): beta})
    assert op(op(mu)) == Vector({(1,): beta**2})

    op = operator_shifted_add(1)
    assert op(mu) == Vector({(2,): 1})

    u = Vector({(1,): 2, (3, 1): -1})
    v = Vector({(2,): 2, (3, 2): -1})
    assert op(u) == v


def test_shifted_remove_row():
    mu = ()
    op = operator_shifted_row_remove(0)
    assert op(mu) == Vector()

    mu = (1,)
    assert op(mu) == Vector({(): 1})

    mu = (2,)
    op = operator_shifted_row_remove(1)
    assert op(mu) == Vector({(1,): 1, (): -beta})

    mu = (4, 1)
    op = operator_shifted_row_remove(3)
    assert op(mu) == Vector({(3, 1,): 1, (2, 1): -beta})


def test_shifted_remove_column():
    mu = ()
    op = operator_shifted_column_remove(0)
    assert op(mu) == Vector()

    mu = (1,)
    assert op(mu) == Vector({(): 1})

    mu = (2,)
    op = operator_shifted_column_remove(1)
    assert op(mu) == Vector({(1,): 1})

    mu = (4, 3, 2)
    op = operator_shifted_column_remove(1)
    assert op(mu) == Vector({(4, 3, 1,): 1, (4, 2, 1): -beta, (3, 2, 1): beta**2})
