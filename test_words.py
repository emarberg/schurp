from words import Word, Permutation
from vectors import Vector


def test_sums():
    assert Word(3, 2, 1) - Word(3, 2, 1) == Vector()
    assert Word(3, 2, 1) + Word(3, 2, 1) == Vector({Word(3, 2, 1): 2})


def test_shuffle():
    u = Word(2, 1)
    v = Word(4, 3)

    assert u * v == \
        Word(2, 1, 4, 3) + Word(2, 4, 1, 3) + Word(2, 4, 3, 1) + \
        Word(4, 2, 1, 3) + Word(4, 2, 3, 1) + Word(4, 3, 2, 1)

    assert u * u == 2 * Word(2, 1, 2, 1) + 4 * Word(2, 2, 1, 1)

    assert u * Word() == Vector.base(u)
    assert Word() * v == Vector.base(v)
