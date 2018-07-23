from schubert import Schubert, InvSchubert, FPFSchubert, MPolynomial
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
