from vectors import Vector
from partitions import StrictPartition


def test_repr():
    def printer(x):
        return 'P' + str(x.parts)

    mu = StrictPartition(3, 2, 1)
    nu = StrictPartition(2, 1)
    u = Vector.base(mu, printer)
    v = Vector.base(nu, printer)
    w = Vector()

    assert str(u) == 'P(3, 2, 1)'
    assert str(v) == 'P(2, 1)'
    assert str(w) == '0'
    assert str(3 * u + 2 * v) in ['3*P(3, 2, 1) + 2*P(2, 1)', '2*P(2, 1) + 3*P(3, 2, 1)']
    assert str(-3 * u + 2 * v) in ['-3*P(3, 2, 1) + 2*P(2, 1)', '2*P(2, 1) - 3*P(3, 2, 1)']
    assert str(-3 * u + -2 * v) in ['-3*P(3, 2, 1) - 2*P(2, 1)', '-2*P(2, 1) - 3*P(3, 2, 1)']


def test_eq():
    mu = StrictPartition(3, 2, 1)
    nu = StrictPartition(2, 1)
    u = Vector.base(mu)
    v = Vector.base(nu)
    w = Vector()

    assert w == u - u == v - v
    assert u == u + w
    assert v == v - w
    assert 2 * v + 3 * u == u + v + u + v + u
    assert -1 * v + 2 * u == u - v + u
