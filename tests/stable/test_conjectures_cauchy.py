from stable.polynomials import Polynomial, X
from stable.symmetric import SymmetricPolynomial
from stable.tableaux import Partition
from stable.polynomials import beta as BETA # noqa
import pytest


def _test(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            for lam1 in Partition.all(n):
                for lam2 in Partition.all(n):
                    print()
                    print()
                    print('* v =', v, ', n =', n, ', mu =', lam1, ', nu =', lam2)
                    print()

                    print('Computing LHS . . .')
                    print()

                    s = Polynomial()
                    for mu in Partition.all(n + max(sum(lam1), sum(lam2))):
                        a = SymmetricPolynomial.stable_grothendieck_doublebar(v, mu, lam1).truncate(n).polynomial('x')
                        b = SymmetricPolynomial.dual_stable_grothendieck(v, mu, lam2).truncate(n).polynomial('y')
                        s += (a * b).truncate(n)
                        print('  ', mu, ':', s, '|', a, '|', b)
                        print()
                    print('LHS =', s)
                    print()
                    print()

                    print('Computing RHS . . .')
                    print()

                    f = Polynomial.one()
                    x = Polynomial.x
                    y = Polynomial.y
                    for i in range(1, v + 1):
                        for j in range(1, v + 1):
                            a = x(i) * y(j)
                            term = Polynomial.one()
                            for e in range(1, n + 1):
                                term += a**e
                            f = (f * term).truncate(n)
                    print('  ', '   :', f)
                    print()

                    t = Polynomial()
                    for kappa in Partition.subpartitions(lam2):
                        a = SymmetricPolynomial.stable_grothendieck_doublebar(v, lam2, kappa).truncate(n)
                        b = SymmetricPolynomial.dual_stable_grothendieck(v, lam1, kappa).truncate(n)
                        t += (f * a.polynomial('x') * b.polynomial('y')).truncate(n)
                        print('  ', kappa, ':', t)
                        print()

                    print('RHS =', t)
                    print()
                    print()
                    print('diff =', s - t)
                    print()
                    assert s == t


def test_grothendieck_cauchy_fast():
    _test(3, 1)


@pytest.mark.slow
def test_grothendieck_cauchy_fast_two():
    _test(4, 2)


@pytest.mark.slow
def test_grothendieck_cauchy():
    _test(6, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_slow():
    _test(9, 6)


def kernel(n, v):
    t = Polynomial.one()
    x = Polynomial.x
    y = Polynomial.y
    for i in range(1, v + 1):
        for j in range(1, v + 1):
            a = x(i) * y(j)
            base = 1 + BETA * x(i) + x(i) * y(j)
            term = Polynomial.one()
            for e in range(1, n + 1):
                term += a**e
            base *= term
            term = Polynomial.one()
            for e in range(1, n + 1):
                term += (-BETA * x(i))**e
            base *= term
            t = (t * base).truncate(n)
    return t


def _test_shifted_q(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            for lam1 in Partition.all(n, strict=True):
                for lam2 in Partition.all(n, strict=True):
                    print()
                    print()
                    print('* v =', v, ', n =', n, ', mu =', lam1, ', nu =', lam2)
                    print()

                    print('Computing LHS . . .')
                    print()

                    s = Polynomial()
                    for mu in Partition.all(n + max(sum(lam1), sum(lam2)), strict=True):
                        a = SymmetricPolynomial.stable_grothendieck_q_doublebar(v, mu, lam1).truncate(n).polynomial('x')
                        b = SymmetricPolynomial.dual_stable_grothendieck_p(v, mu, lam2).truncate(n).polynomial('y')
                        s += (a * b).truncate(n)
                        print('  ', mu, ':', s, '|', a, '|', b)
                        print()
                    print('LHS =', s)
                    print()
                    print()

                    print('Computing RHS . . .')
                    print()

                    f = kernel(n, v)
                    print('  ', '   :', f)
                    print()

                    t = Polynomial()
                    for kappa in Partition.subpartitions(lam2, strict=True):
                        a = SymmetricPolynomial.stable_grothendieck_q_doublebar(v, lam2, kappa).truncate(n)
                        b = SymmetricPolynomial.dual_stable_grothendieck_p(v, lam1, kappa).truncate(n)
                        t += (f * a.polynomial('x') * b.polynomial('y')).truncate(n)
                        print('  ', kappa, ':', t)
                        print()

                    print('RHS =', t)
                    print()
                    print()
                    print('diff =', s - t)
                    print()
                    assert s == t


def test_grothendieck_cauchy_q_fast():
    _test_shifted_q(3, 1)


@pytest.mark.slow
def test_grothendieck_cauchy_q_fast_two():
    _test_shifted_q(4, 2)


@pytest.mark.slow
def test_grothendieck_cauchy_q():
    _test_shifted_q(5, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_q_slow():
    _test_shifted_q(8, 5)


def _test_shifted_p(n_max, v_max):
    for n in range(n_max + 1):
        for v in range(1, v_max + 1):
            for lam1 in Partition.all(n, strict=True):
                for lam2 in Partition.all(n, strict=True):
                    print()
                    print()
                    print('* v =', v, ', n =', n, ', mu =', lam1, ', nu =', lam2)
                    print()

                    print('Computing LHS . . .')
                    print()

                    s = Polynomial()
                    for mu in Partition.all(n + max(sum(lam1), sum(lam2)), strict=True):
                        a = SymmetricPolynomial.stable_grothendieck_p_doublebar(v, mu, lam1).truncate(n).polynomial('x')
                        b = SymmetricPolynomial.dual_stable_grothendieck_q(v, mu, lam2).truncate(n).polynomial('y')
                        s += (a * b).truncate(n)
                        print('  ', mu, ':', s, '|', a, '|', b)
                        print()
                    print('LHS =', s)
                    print()
                    print()

                    print('Computing RHS . . .')
                    print()

                    f = kernel(n, v)
                    print('  ', '   :', f)
                    print()

                    t = Polynomial()
                    for kappa in Partition.subpartitions(lam2, strict=True):
                        a = SymmetricPolynomial.stable_grothendieck_p_doublebar(v, lam2, kappa).truncate(n)
                        b = SymmetricPolynomial.dual_stable_grothendieck_q(v, lam1, kappa).truncate(n)
                        t += (f * a.polynomial('x') * b.polynomial('y')).truncate(n)
                        print('  ', kappa, ':', t)
                        print()

                    print('RHS =', t)
                    print()
                    print()
                    print('diff =', s - t)
                    print()
                    assert s == t


def test_grothendieck_cauchy_p_fast():
    _test_shifted_p(3, 1)


@pytest.mark.slow
def test_grothendieck_cauchy_p_fast_two():
    _test_shifted_p(4, 2)


@pytest.mark.slow
def test_grothendieck_cauchy_p():
    _test_shifted_q(5, 3)


@pytest.mark.slow
def test_grothendieck_cauchy_p_slow():
    _test_shifted_q(8, 5)
