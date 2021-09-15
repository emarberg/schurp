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
import random


def test_star_operator(rank=3, factors=3):
    # fails
    for n in range(1, rank + 1):
        crystal = AbstractGLCrystal
        b = crystal.standard_object(n)
        u = None
        for f in range(1, 1 + factors):
            print(crystal, 'rank =', n, 'factors =', f)
            u = b if u is None else u.tensor(b)
            for i in u.indices:
                for x in u:
                    try:
                        y = u.star_operator(u.e_operator(i, u.star_operator(x)))
                        z = u.f_operator(n - i, x)
                        assert y == z
                    except:
                        print('x =', x, 'i =', i, 'w_0 e_i w_0(x) =', y, 'f_{n-i}(x) =', z, '\n')
                        # u.draw()
                        # assert False
                    try:
                        y = u.star_operator(u.f_operator(i, u.star_operator(x)))
                        z = u.e_operator(n - i, x)
                        assert y == z
                    except:
                        print('x =', x, 'i =', i, 'w_0 f_i w_0(x) =', y, 'e_{n-i}(x) =', z, '\n')
                        # u.draw()
                        # assert False


def _test_components(crystal, rank, factors, x):
    b = crystal.standard_object(rank)
    u = None
    uch = 1
    for f in range(1, 1 + factors):
        print(crystal, 'rank =', rank, 'factors =', f)

        u = b if u is None else u.tensor(b)

        groups = u.group_highest_weights()
        lowgroups = u.group_lowest_weights()
        try:
            for mu in groups:
                rmu = tuple(reversed(mu))
                assert rmu in lowgroups
                assert len(groups[mu]) == len(lowgroups[rmu])
            for mu in lowgroups:
                rmu = tuple(reversed(mu))
                assert rmu in groups
                assert len(groups[rmu]) == len(lowgroups[mu])
        except:
            print('* highest:', groups)
            print('*  lowest:', lowgroups)
            print()
            assert False

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
            assert False


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


def flatten(t):
    if type(t) != tuple:
        return (t,)
    if len(t) == 0:
        return ()
    return (flatten(t[0]) if type(t[0]) == tuple else (t[0],)) + flatten(t[1:])


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
                target_word = crystal.f_operator_on_words(i, flat)
                if target:
                    assert flatten(target) == target_word
                else:
                    assert target_word is None

                target = u.e_operator(i, vertex)
                target_word = crystal.e_operator_on_words(i, flat)
                if target:
                    assert flatten(target) == target_word
                else:
                    assert target_word is None


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


def all_primed_words(max_letter, length):
    for sgn in range(2**length):
        for v in range(max_letter**length):
            a = []
            for _ in range(length):
                a.append((-1)**(sgn % 2) * ((v % max_letter) + 1))
                sgn = sgn >> 1
                v = v >> 1
            yield Word(a)


def random_primed_words(count):
    def generator(max_letter, length):
        for _ in range(count):
            a = []
            for _ in range(length):
                a.append((-1)**(random.randint(0, 1)) * random.randint(1, max_letter))
            yield Word(a)
    return generator


def _test_operators_on_shifted_tableaux(generator, rank, length, verbose=False):
    for word in generator(rank, length):
        tab = word.mixed_insert()[0]
        for i in range(-1, rank):
            target_word = AbstractPrimedQCrystal.f_operator_on_words(i, word)
            target_tab = tab.shifted_crystal_f(i)
            backtrack_tab = target_tab.shifted_crystal_e(i) if target_tab else None
            try:
                if target_word is None:
                    assert target_tab is None
                else:
                    assert target_tab == target_word.mixed_insert()[0]
                    assert backtrack_tab == tab
                if verbose:
                    print(word, '--', i, '-->', target_word)
                    print(tab)
                    print(target_tab)
                    print(backtrack_tab)
                    print()
            except:
                print(word, '--', i, '-->', target_word)
                print(tab)
                print(target_tab)
                print(backtrack_tab)
                print()
                assert False

            target_word = AbstractPrimedQCrystal.e_operator_on_words(i, word)
            target_tab = tab.shifted_crystal_e(i)
            backtrack_tab = target_tab.shifted_crystal_f(i) if target_tab else None
            try:
                if target_word is None:
                    assert target_tab is None
                else:
                    assert target_tab == target_word.mixed_insert()[0]
                    assert backtrack_tab == tab
                if verbose:
                    print(word, '<--', i, '--', target_word)
                    print(tab)
                    print(target_tab)
                    print(backtrack_tab)
                    print()
            except:
                print(word, '<--', i, '--', target_word)
                print(tab)
                print(target_tab)
                print(backtrack_tab)
                print()
                assert False


def test_primed_q_operators_on_tableaux(rank=4, length=4):
    _test_operators_on_shifted_tableaux(all_primed_words, rank, length)


def test_random_primed_q_operators_on_tableaux(rank=10, length=10):
    _test_operators_on_shifted_tableaux(random_primed_words(100), rank, length, True)
