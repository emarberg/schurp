from stable.utils import G, GP, GS
from stable.permutations import Permutation


def test_grothendieck():
    n = 3
    w = Permutation(4, 3, 2, 1)
    f = G(n, (3, 2, 1))
    print(f)
    g = G(n, w)
    print(g)
    assert f == g


def test_grothendieck_p():
    n = 3
    w = Permutation(4, 3, 2, 1)
    f = GP(n, (2,))
    print(f)
    g = GP(n, w)
    print(g)
    assert f == g


def test_grothendieck_s():
    n = 3
    w = Permutation(-1)
    f = GS(n, (1,))
    print(f)
    g = GS(n, w)
    print(g)
    assert f == g
