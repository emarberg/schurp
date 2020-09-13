from words import Word
from vectors import Vector


def test_sums():
    s = {1, 2, 3}
    assert Word(3, 2, 1, subset=s) - Word(3, 2, 1, subset=s) == Vector()
    assert Word(3, 2, 1, subset=s) + Word(3, 2, 1, subset=s) == Vector({Word(3, 2, 1, subset=s): 2})


def test_shuffle():
    s = {1, 2, 3, 4}
    u = Word(2, 1, subset={1, 2})
    v = Word(4, 3, subset={3, 4})

    assert u * v == \
        Word(2, 1, 4, 3, subset=s) + Word(2, 4, 1, 3, subset=s) + Word(2, 4, 3, 1, subset=s) + \
        Word(4, 2, 1, 3, subset=s) + Word(4, 2, 3, 1, subset=s) + Word(4, 3, 2, 1, subset=s)

    assert u * Word() == Vector.base(u)
    assert Word() * v == Vector.base(v)


def test_coproduct():
    s = {1, 2, 3, 4}
    w = Word(2, 1, 2, 4, 4, subset=s)

    u = Word(2, 1, 2, subset={1, 2})
    v = Word(4, 4, subset={3, 4})
    assert w.coproduct({1, 2}, {3, 4}) == Vector.base((u, v))
    assert w.coproduct({3, 4}, {1, 2}) == Vector()

    h = Word(2, 1, 2, subset={1, 2, 3})
    x = Word(2, 1, 2, subset={1, 2})
    y = Word(subset={3})
    z = Word(4, 4, subset={4})
    assert w.coproduct({1, 2}, {3}, {4}) == Vector.base((x, y, z))
    assert w.coproduct({1, 2}, {4}, {3}) == Vector.base((x, z, y))
    assert w.coproduct({1, 2, 3}, {4}) == Vector.base((h, z))

    assert w.coproduct({1, 4}, {2, 3}) == Vector()
