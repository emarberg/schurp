from stable.permutations import Permutation
from stable.symmetric import SymmetricPolynomial
import pytest


def _test(r, upper):
    w = Permutation.longest_element(r, signed=True)
    mu = tuple(i for i in range(r, 0, -1))
    for n in range(upper):
        f = SymmetricPolynomial.stable_grothendieck_s(n, mu, degree_bound=n)
        g = w.signed_involution_stable_grothendieck(n, degree_bound=n)
        print('f =', f)
        print()
        print('g =', g)
        print()
        print(f - g)
        print()
        print()
        print()
        assert f == g


def test_one():
    _test(1, 10)


def test_two():
    _test(2, 8)


@pytest.mark.slow
def test_three():
    _test(3, 10)


@pytest.mark.slow
def test_four():
    _test(4, 10)


@pytest.mark.slow
def test_five():
    _test(5, 8)


@pytest.mark.slow
def test_six():
    _test(6, 6)
