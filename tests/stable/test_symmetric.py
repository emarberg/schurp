from stable.symmetric import SymmetricMonomial, SymmetricPolynomial
from stable.tableaux import Tableau, Partition
from stable.polynomials import Polynomial, X
from stable.utils import (
    G, GQ, GP, G_doublebar, GP_doublebar, GQ_doublebar, 
    jp, decomposition_jp, 
)
from stable.polynomials import beta as BETA # noqa
from stable.vectors import Vector

import itertools
import pytest


def test_doublebar():
    assert G(1, (), (1,)) == 0
    assert G_doublebar(1, (), (1,)) == 0

    assert GP(1, (), (1,)) == 0
    assert GP_doublebar(1, (), (1,)) == 0

    assert GQ(1, (), (1,)) == 0
    assert GQ_doublebar(1, (), (1,)) == 0


def test_multipeak_formula():
    """
    Tests that explicit summation formula for

      K^{(beta)}_{(2,1)}(x_1, x_2, ..., x_N)

    is the same as

      GQ^{(beta)}_{(2,1)}(x_1, x_2, ... ,x_N)

    """
    mu = (2, 1)
    for N in range(6):
        ans = 0
        for n in range(2, N + 1):
            for i in itertools.combinations(range(1, N + 1), n):
                term = BETA**(n - 3)
                for j in range(n):
                    x = X(i[j])
                    term *= (2 * x + BETA * x**2)
                sigma = (n - 1) * (n - 2) // 2
                for j in range(n):
                    for k in range(j + 1, n):
                        x = X(i[j])
                        y = X(i[k])
                        sigma += BETA * (x + y + BETA * x * y)
                print(i, term * sigma)
                ans += term * sigma
        bns = GQ(N, mu).polynomial()
        print(ans)
        print()
        print(bns)
        print()
        print(ans - bns)
        assert ans == bns


def test_overline_multipeak_formula():
    """
    Tests that explicit summation formula for

      oK^{(-1)}_{(2,1)}(x_1, x_2, ..., x_N)

    is the same as

      GP^{(-1)}_{(2,1)}(x_1, x_2, ... ,x_N)

    """
    mu = (2, 1)
    for N in range(7):
        ans = 0
        for n in range(2, N + 1):
            for i in itertools.combinations(range(1, N + 1), n):
                term = BETA**(n - 3)
                for t in range(n):
                    term *= X(i[t])
                for j in range(n):
                    x = X(i[j])
                    for k in range(j + 1, n):
                        y = X(i[k])
                        pi = 1
                        if j + 1 == k:
                            pi *= BETA * (x + y + BETA * x * y)
                        else:
                            pi *= (1 + BETA * x) * (1 + BETA * y)
                            for t in range(j + 1, k):
                                pi *= (2 + BETA * X(i[t]))
                        print(i, j, k, term * pi)
                        print()
                        ans += term * pi
        bns = GP(N, mu).polynomial()
        print()
        print(ans)
        print()
        print(bns)
        print()
        print(ans - bns)
        assert ans == bns


def test_destandardize():
    mu = (2, 1)
    n = 3
    assert set(SymmetricMonomial._destandardize(n, mu)) == {
        (1, 2, 0),
        (1, 0, 2),
        (0, 1, 2),
        (2, 1, 0),
        (2, 0, 1),
        (0, 2, 1),
    }
    assert len(list(SymmetricMonomial._destandardize(n, mu))) == 6

    mu = (2, 2)
    n = 3
    assert set(SymmetricMonomial._destandardize(n, mu)) == {
        (2, 2, 0),
        (2, 0, 2),
        (0, 2, 2),
    }
    assert len(list(SymmetricMonomial._destandardize(n, mu))) == 3


def test_basic():
    x = SymmetricMonomial(1, (1,))
    assert type(x) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * x
    print(x)
    print(f)
    print(type(f), f, SymmetricMonomial(1, (2,)), f - SymmetricMonomial(1, (2,)))
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == SymmetricMonomial(1, (2,))

    x = SymmetricMonomial(2, (1,))
    print(x)
    print(x * x)
    assert x * x == SymmetricMonomial(2, (2,)) + 2 * SymmetricMonomial(2, (1, 1))

    f = x + 1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * 1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == x

    f = x * 0
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * -1
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * -1 + x
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f.is_zero()

    x = SymmetricMonomial(2, (1,))
    y = SymmetricMonomial(2, (1, 1))
    assert type(y) in [SymmetricMonomial, SymmetricPolynomial]

    f = x * y
    assert type(f) in [SymmetricMonomial, SymmetricPolynomial]
    assert f == SymmetricMonomial(2, (2, 1))


def test_polynomials():
    m = SymmetricMonomial(5, ())
    assert m.polynomial() == Polynomial.one()

    m = SymmetricMonomial(5, (1,))
    x = Polynomial.x
    assert m.polynomial() == x(1) + x(2) + x(3) + x(4) + x(5)

    m = SymmetricMonomial(3, (2, 1,))
    x = Polynomial.y
    assert m.polynomial('y') == \
        x(1) * x(2)**2 + x(1)**2 * x(2) + \
        x(1) * x(3)**2 + x(1)**2 * x(3) + \
        x(3) * x(2)**2 + x(3)**2 * x(2)


def test_slow_symmetric_functions():
    for mu in [(), (1,), (1, 1), (2,), (2, 1)]:
        for n in range(4):
            f = SymmetricPolynomial._slow_schur(n, mu)
            g = SymmetricPolynomial._slow_stable_grothendieck(n, mu)
            h = SymmetricPolynomial._slow_schur_s(n, mu)
            k = SymmetricPolynomial._slow_stable_grothendieck_s(n, mu)
            assert g.lowest_degree_terms() == f
            assert k.lowest_degree_terms() == h


def test_slow_dual_stable_grothendieck():
    n = 3
    mu = (1,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(n, mu)
    print(g)
    assert g == SymmetricMonomial(n, (1,))

    n = 3
    mu = (2,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(n, mu)
    print(g)
    assert g == SymmetricMonomial(n, (1, 1)) + SymmetricMonomial(n, (2,))

    n = 3
    mu = (2, 1)
    for t in Tableau.semistandard_rpp(n, mu):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck(n, mu)
    print(g)
    assert g == -BETA * SymmetricMonomial(n, (1, 1)) - BETA * SymmetricMonomial(n, (2,)) + 2 * SymmetricMonomial(n, (1, 1, 1)) + SymmetricMonomial(n, (2, 1))


def test_slow_dual_stable_grothendieck_pq():
    n = 3
    mu = (1,)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(n, mu)
    print(g)
    print(h)
    print()
    assert g == SymmetricMonomial(n, (1,))
    assert h == 2 * SymmetricMonomial(n, (1,))

    n = 3
    mu = (2,)
    for t in Tableau.semistandard_marked_rpp(3, mu, diagonal_nonprimes=False):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu)
    print(g)
    assert g == \
        SymmetricPolynomial._slow_dual_stable_grothendieck(n, (1, 1)) + \
        SymmetricPolynomial._slow_dual_stable_grothendieck(n, (2,))

    n = 3
    mu = (2, 1)
    for t in Tableau.semistandard_marked_rpp(n, mu):
        print(t)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(n, mu)
    print(g)
    print(h)
    print()
    assert g == -BETA * SymmetricMonomial(n, (1, 1)) + -BETA * SymmetricMonomial(n, (2,)) + 2 * SymmetricMonomial(n, (1, 1, 1)) + SymmetricMonomial(n, (2, 1))

    n = 5
    mu = (3, 2)
    g = SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu)
    h = SymmetricPolynomial._slow_dual_stable_grothendieck_q(n, mu)

    n = 6
    mu = (3, 2, 1)
    assert SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu) == \
        SymmetricPolynomial._slow_dual_stable_grothendieck(n, mu)


def test_decompose():
    schur = SymmetricPolynomial.schur
    schur_P = SymmetricPolynomial.schur_p  # noqa
    dual_grothendieck = SymmetricPolynomial.dual_stable_grothendieck
    grothendieck = SymmetricPolynomial.stable_grothendieck
    GP = SymmetricPolynomial.stable_grothendieck_p  # noqa
    GQ = SymmetricPolynomial.stable_grothendieck_q  # noqa

    n = 6
    mu = (3, 2, 1)
    nu = (4, 2, 1)

    f = schur(n, mu) * schur(n, mu)
    exp = SymmetricPolynomial.schur_expansion(f)
    assert f == sum([schur(n, a) * coeff for a, coeff in exp.items()])

    f = schur_P(n, nu)
    exp = SymmetricPolynomial.schur_expansion(f)
    assert f == sum([schur(n, a) * coeff for a, coeff in exp.items()])

    f = GP(n, nu)
    exp = SymmetricPolynomial.grothendieck_expansion(f)
    assert f == sum([grothendieck(n, a) * coeff for a, coeff in exp.items()])

    f = GQ(n, nu)
    exp = SymmetricPolynomial.grothendieck_expansion(f)
    assert f == sum([grothendieck(n, a) * coeff for a, coeff in exp.items()])

    f = SymmetricPolynomial.dual_stable_grothendieck_p(n, nu)
    exp = SymmetricPolynomial.dual_grothendieck_expansion(f)
    assert f == sum([dual_grothendieck(n, a) * coeff for a, coeff in exp.items()])

    f = SymmetricPolynomial.dual_stable_grothendieck_q(n, nu)
    exp = SymmetricPolynomial.dual_grothendieck_expansion(f)
    assert f == sum([dual_grothendieck(n, a) * coeff for a, coeff in exp.items()])


@pytest.mark.slow
def test_symmetric_functions():
    nn = 6
    for mu in Partition.all(nn):
        for nu in Partition.all(nn, strict=True):
            for n in range(nn):
                print(n, mu, nu)
                print()

                f = SymmetricPolynomial.schur(n, mu, nu)
                g = SymmetricPolynomial.stable_grothendieck(n, mu, nu)
                h = SymmetricPolynomial.dual_stable_grothendieck(n, mu, nu)

                fs = SymmetricPolynomial._slow_schur(n, mu, nu)
                gs = SymmetricPolynomial._slow_stable_grothendieck(n, mu, nu)
                hs = SymmetricPolynomial._slow_dual_stable_grothendieck(n, mu, nu)

                if f != fs:
                    print(f)
                    print(fs)
                    print()
                if g != gs:
                    print(g)
                    print(gs)
                    print()
                if h != hs:
                    print(h)
                    print(hs)
                    print()
                    print()

                assert f == fs
                assert g == gs
                assert h == hs

                hh = SymmetricPolynomial.schur_s(n, mu, nu)
                kk = SymmetricPolynomial.stable_grothendieck_s(n, mu, nu)

                if mu == nu:
                    assert f == 1
                    assert g == 1
                    assert h == 1
                    assert fs == 1
                    assert gs == 1
                    assert hs == 1
                    assert hh == 1
                    assert kk == 1

                if not Partition.contains(mu, nu):
                    assert f == 0
                    assert g == 0
                    assert h == 0
                    assert fs == 0
                    assert gs == 0
                    assert hs == 0
                    assert hh == 0
                    assert kk == 0

                print(f)
                print(g)
                print()
                print(hh)
                print(kk)
                print()
                assert g.lowest_degree_terms() == f
                assert kk.lowest_degree_terms() == hh


@pytest.mark.slow
def test_strict_symmetric_functions():
    nn = 5
    for mu in Partition.all(nn, strict=True):
        for nu in Partition.all(nn, strict=True):
            for n in range(nn):
                print(n, mu, nu)
                print()

                # Schur-P and GP

                f = SymmetricPolynomial.schur_p(n, mu, nu)
                g = SymmetricPolynomial.stable_grothendieck_p(n, mu, nu)
                h = SymmetricPolynomial.dual_stable_grothendieck_p(n, mu, nu)

                fs = SymmetricPolynomial._slow_schur_p(n, mu, nu)
                gs = SymmetricPolynomial._slow_stable_grothendieck_p(n, mu, nu)
                hs = SymmetricPolynomial._slow_dual_stable_grothendieck_p(n, mu, nu)

                if f != fs:
                    print(f)
                    print(fs)
                    print()

                if g != gs:
                    print(g)
                    print(gs)
                    print()

                if h != hs:
                    print(h)
                    print(hs)
                    print()
                    print()
                    print()

                assert f == fs
                assert g == gs
                assert h == hs

                if mu == nu:
                    assert f == 1
                    assert g == 1
                    assert h == 1
                    assert fs == 1
                    assert gs == 1
                    assert hs == 1

                if not Partition.contains(mu, nu):
                    assert f == 0
                    assert g == 0
                    assert h == 0
                    assert fs == 0
                    assert gs == 0
                    assert hs == 0

                # Schur-Q and GQ

                f = SymmetricPolynomial.schur_q(n, mu, nu)
                g = SymmetricPolynomial.stable_grothendieck_q(n, mu, nu)
                h = SymmetricPolynomial.dual_stable_grothendieck_q(n, mu, nu)

                fs = SymmetricPolynomial._slow_schur_q(n, mu, nu)
                gs = SymmetricPolynomial._slow_stable_grothendieck_q(n, mu, nu)
                hs = SymmetricPolynomial._slow_dual_stable_grothendieck_q(n, mu, nu)

                if f != fs:
                    print(f)
                    print(fs)
                    print()

                if g != gs:
                    print(g)
                    print(gs)
                    print()

                if h != hs:
                    print(h)
                    print(hs)
                    print()
                    print()
                    print()

                assert f == fs
                assert g == gs
                assert h == hs

                if mu == nu:
                    assert f == 1
                    assert g == 1
                    assert h == 1
                    assert fs == 1
                    assert gs == 1
                    assert hs == 1

                if not Partition.contains(mu, nu):
                    assert f == 0
                    assert g == 0
                    assert h == 0
                    assert fs == 0
                    assert gs == 0
                    assert hs == 0
