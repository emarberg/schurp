from schubert import (
    Schubert,
    Grothendieck,
    InvSchubert,
    FPFSchubert,
    FPFGrothendieck,
    MPolynomial
)
import schubert
from permutations import Permutation


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


def test_product_grothendieck():
    n = 6
    g = list(Permutation.all(n))

    for w in g:
        a = Grothendieck.get(w) == Grothendieck.product(w)
        b = w.is_dominant()
        c = Grothendieck.get(w) == Schubert.get(w)
        assert a == b
        assert b == c


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
