from schubert import (
    Schubert,
    Grothendieck,
    InvSchubert,
    FPFSchubert,
    FPFGrothendieck,
    MPolynomial,
    Operator
)
from permutations import Permutation
import pytest


def test_formal_operators():
    D = Operator.create(1) # noqa
    x1 = MPolynomial.monomial(1)
    x2 = MPolynomial.monomial(2)
    x3 = MPolynomial.monomial(3)

    (x1 * x2).toggle(2) == x1 * x3
    (x1 * x2).toggle(1) == x1 * x2
    x1.toggle(1) == x2
    x2.toggle(1) == x1

    E = D * x1 # noqa
    F = x1 * D # noqa

    print('1:', D, type(D))
    print('2:', E, type(E))
    print('3:', F, type(F))
    print()

    D + D
    D - D
    D * D
    assert D * x1 != x1 * D
    assert D * x2 != x2 * D
    assert x3 * D == D * x3

    assert D * x1 == Operator.create() + x2 * D

    A = Operator.create(1) * x1 # noqa
    B = Operator.create(2) * x2 # noqa

    print(A)
    print(Operator.create(1) * x2 * Operator.create(1))
    print(A * Operator.create(1))
    print(A * Operator.create(1) - Operator.create(1))
    print()

    assert A * Operator.create() == A
    assert Operator.create() * A == A
    assert A * Operator.create(1) == Operator.create(1)
    assert Operator.create(1) * A == 0
    assert A * A == A
    assert B * B == B
    assert A * B * A == B * A * B

    a = Operator.create(1) * (1 - x2)
    b = Operator.create(2) * (1 - x3)
    assert a * a == a
    assert b * b == b
    assert a * b * a == b * a * b

    A = Operator.create(1) * x1 * (1 - x2) # noqa
    B = Operator.create(2) * x2 * (1 - x3) # noqa
    print(A * B * A)
    print(B * A * B)
    print()

    assert A * A == A
    assert B * B == B
    assert A * B * A == B * A * B


def test_formal_operators_experiments():
    x = lambda i: MPolynomial.monomial(i) # noqa
    A = lambda i: Operator.create(i) * x(i) * (1 - x(i + 1)) # noqa
    a = lambda i: Operator.create(i) * (1 - x(i + 1)) # noqa
    print(A(2) * A(1) - a(2) * a(1) * x(1)**2)
    print()
    print(a(1) * x(1))
    assert False


def test_divided_differences():
    x = MPolynomial.monomial(3)
    y = MPolynomial.monomial(4)
    f = x**3
    assert f.divided_difference(3) == x**2 + x * y + y**2
    assert (x * y).isobaric_divided_difference(3) == x * y


def test_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation.s_i(3)
    assert Schubert.get(w) == x + y + z

    w = Permutation(3, 4, 1, 2)
    assert Schubert.get(w) == x * x * y * y

    w = Permutation(1, 4, 3, 2)
    assert Schubert.get(w) == x**2 * y + x**2 * z + x * y**2 + x * y * z + y**2 * z


def test_inv_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation()
    assert InvSchubert.get(w) == 1

    w = Permutation(1, 4, 3, 2)
    assert InvSchubert.get(w) == x**2 + 2 * x * y + x * z + y**2 + y * z

    w = Permutation(3, 4, 1, 2)
    assert InvSchubert.get(w) == x**2 * y + x * y**2


def test_fpf_schubert():
    x = MPolynomial.monomial(1)
    y = MPolynomial.monomial(2)
    z = MPolynomial.monomial(3)

    w = Permutation(2, 1, 4, 3)
    assert FPFSchubert.get(w) == 1

    w = Permutation(3, 4, 1, 2)
    assert FPFSchubert.get(w) == x + y

    w = Permutation(4, 3, 2, 1)
    assert FPFSchubert.get(w) == x * x + x * y + x * z + y * z


@pytest.mark.slow
def test_product_grothendieck():
    n = 6
    g = list(Permutation.all(n))

    for w in g:
        a = Grothendieck.get(w) == Grothendieck.product(w)
        b = w.is_dominant()
        c = Grothendieck.get(w) == Schubert.get(w)
        assert a == b
        assert b == c


@pytest.mark.slow
def test_fpf_grothendieck():
    n = 6
    g = list(Permutation.fpf_involutions(n))

    for w in g:
        a = FPFGrothendieck.get(w) == FPFGrothendieck.product(w)
        b = w.is_fpf_dominant()
        assert a == b
        f = FPFGrothendieck.get(w)
        g = sum([
            (-1) ** (u.length() - w.fpf_involution_length()) * Grothendieck.get(u)
            for u in w.get_symplectic_hecke_atoms()
        ])
        print(w)
        print(f)
        print(g)
        print()
        print(list(w.get_symplectic_hecke_atoms()))
        print()
        assert f == g


def test_min():
    w = Permutation(3, 5, 1, 6, 2, 4)
    f = FPFSchubert.get(w)
    assert FPFSchubert.least_term(f) == {2: 1, 4: 1}
