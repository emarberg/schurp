from operators import (
    raising, lowering, diagonal, 
    i_raising, i_diagonal, i_lowering,
    f_raising, f_diagonal, f_lowering,
)
from vectors import Vector
from permutations import Permutation


def test_raising():
    n = 3
    u = Permutation()
    v = Permutation(2, 1, 3)
    w = Permutation(1, 3, 2)
    assert raising(u, n) == Vector({v: 1, w: 2})
    assert raising(Permutation(3, 2, 1), n) == 0


def test_lowering():
    n = 3
    t = Permutation(2, 3, 1)
    u = Permutation(3, 1, 2)
    v = Permutation(2, 1, 3)
    w = Permutation(1, 3, 2)
    assert lowering(t, n) == Vector({v: 1, w: 1})
    assert lowering(u, n) == Vector({v: 1, w: 3})
    assert lowering(Permutation(), n) == 0


def test_diagonal():
    n = 3
    t = Permutation(2, 3, 1)
    u = Permutation(3, 1, 2)
    v = Permutation(2, 1, 3)
    w = Permutation(1, 3, 2)
    assert diagonal(t, n) == Vector({t: 1})
    assert diagonal(u, n) == Vector({u: 1})
    assert diagonal(v, n) == Vector({v: -1})
    assert diagonal(w, n) == Vector({w: -1})


def test_sl_repn():
    for n in range(6):
        for w in Permutation.all(n):
            assert diagonal(raising(w, n), n) - raising(diagonal(w, n), n) == 2 * raising(w, n)
            assert diagonal(lowering(w, n), n) - lowering(diagonal(w, n), n) == -2 * lowering(w, n)
            assert raising(lowering(w, n), n) - lowering(raising(w, n), n) == diagonal(w, n)


def test_sl_i_repn():
    for n in range(4, 5):
        for w in Permutation.involutions(n):
            a = i_diagonal(i_raising(w, n), n)
            b = i_raising(i_diagonal(w, n), n)
            c = 2 * i_raising(w, n)

            print()
            print('n =', n, ', w =', w, ':', w.length(), w.involution_length())

            # print()
            # print('raising:')
            # print()
            # print(a)
            # print(b)
            # print(c)
            # print()
            # assert a - b == c

            a = i_diagonal(i_lowering(w, n), n)
            b = i_lowering(i_diagonal(w, n), n)
            c = -2 * i_lowering(w, n)

            # print()
            # print('lowering:')
            # print()
            # print(a)
            # print(b)
            # print(c)
            # print()
            # assert a - b == c

            a = i_raising(i_lowering(w, n), n)
            b = i_lowering(i_raising(w, n), n)
            c = i_diagonal(w, n)

            # print()
            # print('diagonal:')
            # print()
            # print(i_raising(w, n))
            # print(i_lowering(i_raising(w, n), n))
            # print()
            # print(a)
            # print(b)
            print()
            print(a - b, '=?=', c)
            print()


def test_sl_f_repn():
    for n in range(6, 8, 2):
        for w in Permutation.fpf_involutions(n):
            a = f_diagonal(f_raising(w, n), n)
            b = f_raising(f_diagonal(w, n), n)
            c = 2 * f_raising(w, n)

            print()
            print('n =', n, ', w =', w, ':', w.length(), w.fpf_involution_length())

            # print()
            # print('raising:')
            # print()
            # print(a)
            # print(b)
            # print(c)
            # print()
            assert a - b == c

            a = f_diagonal(f_lowering(w, n), n)
            b = f_lowering(f_diagonal(w, n), n)
            c = -2 * f_lowering(w, n)

            # print()
            # print('lowering:')
            # print()
            # print(a)
            # print(b)
            # print(c)
            # print()
            assert a - b == c

            a = f_raising(f_lowering(w, n), n)
            b = f_lowering(f_raising(w, n), n)
            c = f_diagonal(w, n)

            print()
            print('diagonal:')
            print()
            print(a)
            print(b)
            print()
            print(a - b, '=?=', c)
            print()
            print(a - b == c)
