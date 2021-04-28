from crystals import(
    AbstractGLCrystal,
    AbstractQCrystal,
    AbstractPrimedQCrystal,
)
from partitions import Partition
from symmetric import (
    SchurQ,
    SchurP,
    Schur
)
from words import Word


def _test_components(crystal, rank, factors, x):
    b = crystal.standard_object(rank)
    u = None
    uch = 1
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)

        u = b if u is None else u.tensor(b)

        groups = u.group_highest_weights()
        for g in groups.values():
            for i in range(len(g) - 1):
                assert u.isomorphic_subcrystals(g[i], g[i + 1])
        weights = list(groups)
        for i in range(len(weights) - 1):
            mu, nu = weights[i:i + 2]
            assert not u.isomorphic_subcrystals(groups[mu][0], groups[nu][0])

        uch = uch * x
        expected = {tuple(k.mu): v for k, v in uch.items() if len(k.mu) <= rank}
        actual = {Partition.trim(k): len(v) for k, v in groups.items()}
        try:
            assert expected == actual
        except:
            print('expected =', expected)
            print('  actual =', actual)
            print()


def test_gl_components(rank=3, factors=3):
    x = Schur(1)
    for rank in range(1, rank + 1):
        _test_components(AbstractGLCrystal, rank, factors, x)


def test_q_components(rank=3, factors=3):
    x = SchurP(1)
    for rank in range(2, rank + 1):
        _test_components(AbstractQCrystal, rank, factors, x)


def test_primed_q_components(rank=3, factors=3):
    x = SchurQ(1)
    for rank in range(2, rank + 1):
        _test_components(AbstractPrimedQCrystal, rank, factors, x)


def flatten(t, cl=tuple):
    if type(t) != cl:
        return (t,)
    if len(t) == 0:
        return cl()
    return (flatten(t[0]) if type(t[0]) == cl else (t[0],)) + flatten(t[1:])


def _test_operators_on_words(crystal, rank, factors):
    b = crystal.standard_object(rank)
    u = None
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)
        u = b if u is None else u.tensor(b)
        for vertex in u:
            flat = flatten(vertex)
            for i in u.indices:
                target = u.f_operator(i, vertex)
                if target:
                    assert flatten(target) == crystal.f_operator_on_words(i, flat)
                else:
                    assert crystal.f_operator_on_words(i, flat) is None

                target = u.e_operator(i, vertex)
                if target:
                    assert flatten(target) == crystal.e_operator_on_words(i, flat)
                else:
                    assert crystal.e_operator_on_words(i, flat) is None


def test_gl_operators_on_words(rank=3, factors=3):
    u = Word(1, 3, 3, 1, 2, 1, 2)
    v = Word(2, 3, 3, 1, 2, 1, 2)
    w = Word(3, 3, 3, 1, 2, 1, 2)

    assert AbstractGLCrystal.f_operator_on_words(2, u) is None
    assert AbstractGLCrystal.f_operator_on_words(2, v) == w
    assert AbstractGLCrystal.e_operator_on_words(2, w) == v

    for rank in range(1, rank + 1):
        _test_operators_on_words(AbstractGLCrystal, rank, factors)


def test_q_operators_on_words(rank=4, factors=4):
    for rank in range(2, rank + 1):
        _test_operators_on_words(AbstractQCrystal, rank, factors)


def test_primed_q_operators_on_words(rank=4, factors=4):
    for rank in range(2, rank + 1):
        _test_operators_on_words(AbstractPrimedQCrystal, rank, factors)
