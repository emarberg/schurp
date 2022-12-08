from stable.vectors import Vector


def test_eq():
    mu = (3, 2, 1)
    nu = (2, 1)
    u = Vector.base(mu)
    v = Vector.base(nu)
    w = Vector()

    assert w == u - u == v - v
    assert u == u + w
    assert v == v - w
    assert 2 * v + 3 * u == u + v + u + v + u
    assert -1 * v + 2 * u == u - v + u
