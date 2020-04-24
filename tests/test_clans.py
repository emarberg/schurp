from clans import Clan
from permutations import Permutation


def test_init():
    w = Clan([1, True, 1, False])
    assert w.cycles() == [(1, 3)]
    assert w(1) == 3
    assert w(2) is True
    assert w(3) == 1
    assert w(4) is False


def test_multiply():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])
    assert 1 * w == w
    assert 3 * w == v


def test_multiply_b():
    w = Clan([True, True, False, False, False, True, True], Clan.TYPE_B)
    v = Clan([True, 1, 1, False, 2, 2, True], Clan.TYPE_B)
    u2 = Clan([1, True, 1, False, 2, True, 2], Clan.TYPE_B)
    u = Clan([True, 1, 2, False, 1, 2, True], Clan.TYPE_B)
    t = Clan([1, True, 2, False, 1, True, 2], Clan.TYPE_B)
    s = Clan([True, 1, 2, False, 2, 1, True], Clan.TYPE_B)
    r = Clan([1, True, 2, False, 2, True, 1], Clan.TYPE_B)
    q = Clan([1, 2, True, False, True, 2, 1], Clan.TYPE_B)
    p = Clan([1, 2, 3, True, 3, 2, 1], Clan.TYPE_B)

    assert 0 * w == w
    assert 1 * w == v
    assert 2 * w == w

    assert 0 * v == u
    assert 1 * v == v
    assert 2 * v == u2

    assert 0 * u == u
    assert 1 * u == s
    assert 2 * u == t

    assert 0 * s == s
    assert 1 * s == s
    assert 2 * s == r

    assert 0 * r == r
    assert 1 * r == q
    assert 2 * r == r

    assert 0 * q == p
    assert 1 * q == q
    assert 2 * q == q


def test_atoms():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])
    assert set(w.get_clan_words()) == {(2, 3)}
    assert set(v.get_clan_words()) == {(2,)}
    assert w.get_atoms() == {Permutation(1, 3, 4, 2)}
    assert v.get_atoms() == {Permutation(1, 3, 2)}


def test_all():
    w = Clan([1, True, 1, False])
    v = Clan([1, True, False, 1])
    a = set(Clan.all(2, 2))
    assert v in a and w in a
    assert len(a) == 6 + 12 + 3


def test_all_b():
    assert len(list(Clan.all_b(2, 1))) == 25


def test_atoms_b(n=4):
    for p in range(1, n):
        for q in [p - 1, p]:
            for clan in Clan.all_b(p, q):
                print(clan)
                z = -clan.richardson_springer_map()
                assert set(clan.get_atoms()).issubset(set(z.get_atoms()))
