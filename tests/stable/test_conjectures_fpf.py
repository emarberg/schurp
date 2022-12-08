from stable.permutations import Permutation
from stable.symmetric import SymmetricPolynomial
import pytest


def test_generate():
    assert set(Permutation.fpf_involutions(2)) == {Permutation(2, 1)}
    assert set(Permutation.fpf_involutions(3)) == set()
    assert set(Permutation.fpf_involutions(4)) == {
        Permutation(2, 1, 4, 3),
        Permutation(3, 4, 1, 2),
        Permutation(4, 3, 2, 1),
    }


def test_fpf_involution_shape():
    w = Permutation(2, 1)
    assert w.fpf_involution_shape() == ()
    assert w._fpf_grassmannian_shape() == ()

    w = Permutation(4, 3, 2, 1)
    print(w.fpf_involution_code())
    assert w.fpf_involution_shape() == (2,)
    assert w._fpf_grassmannian_shape() == (2,)

    w = Permutation(3, 4, 1, 2)
    assert w.fpf_involution_shape() == (1,)
    assert w._fpf_grassmannian_shape() == (1,)


def _test(it, upper):
    for w in it:
        nu = w.fpf_involution_shape()
        if nu:
            for n in range(upper):
                f = SymmetricPolynomial.stable_grothendieck_p(nu, n, n)
                g = w.symplectic_stable_grothendieck(degree_bound=n)
                print(w, nu)
                print()
                print(f)
                print()
                print(g)
                print()
                print(f - g)
                print()
                print()
                print()
                assert f == g


def test_one():
    _test(Permutation.fpf_involutions(1), 12)


def test_two():
    _test(Permutation.fpf_involutions(2), 8)


def test_three():
    _test(Permutation.fpf_involutions(3), 10)


@pytest.mark.slow
def test_four():
    _test(Permutation.fpf_involutions(4), 8)


@pytest.mark.slow
def test_six():
    _test(Permutation.fpf_involutions(6), 6)
