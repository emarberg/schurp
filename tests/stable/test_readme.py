from stable.utils import *


def test_functions():
    s = schur
    g = dual_grothendieck

    # number of variables
    n = 6

    # define a partition
    mu = (3, 2, 1)
    nu = (4, 2, 1)

    # compute Schur function in variables x_1, x_2, ..., x_n
    schur(n, mu)
    s(n, mu)

    # compute Schur P-function in variables x_1, x_2, ..., x_n
    schur_P(n, mu)
    P(n, mu)

    # compute Schur Q-function in variables x_1, x_2, ..., x_n
    schur_Q(n, mu)
    Q(n, mu)

    # compute Schur S-function in variables x_1, x_2, ..., x_n
    schur_S(n, mu)
    S(n, mu)

    # compute Schur decomposition
    f = schur(n, mu) * schur(n, mu)
    SymmetricPolynomial.schur_expansion(f)

    f = schur_P(n, nu)
    SymmetricPolynomial.schur_expansion(f)

    # compute stable Grothendieck polynomial in variables x_1, x_2, ..., x_n
    grothendieck(n, mu)
    G(n, mu)

    # compute K-theoretic Schur P-function in variables x_1, x_2, ..., x_n
    grothendieck_P(n, mu)
    GP(n, mu)

    # compute K-theoretic Schur Q-function in variables x_1, x_2, ..., x_n
    grothendieck_Q(n, mu)
    GQ(n, mu)

    # compute K-theoretic Schur S-function in variables x_1, x_2, ..., x_n
    grothendieck_S(n, mu)
    GS(n, mu)

    # compute stable Grothendieck decomposition
    f = GP(n, nu)
    SymmetricPolynomial.grothendieck_expansion(f)

    f = GQ(n, nu)
    SymmetricPolynomial.grothendieck_expansion(f)

    # compute dual stable Grothendieck polynomial in variables x_1, x_2, ..., x_n
    g(n, mu)
    dual_grothendieck(n, mu)

    # compute shifted dual stable Grothendieck polynomials in variables x_1, x_2, ..., x_n
    gp(n, mu)
    dual_grothendieck_P(n, mu)

    gq(n, mu)
    dual_grothendieck_Q(n, mu)

    # compute dual stable Grothendieck decomposition
    f = gp(n, nu)
    SymmetricPolynomial.dual_grothendieck_expansion(f)

    f = gq(n, nu)
    SymmetricPolynomial.dual_grothendieck_expansion(f)

    # add, subtract, and multiply polynomials
    f = P(n, mu)
    g = Q(n, mu)
    8 * f == g

    alpha, beta = (3,), (4,)
    f = GP(n, alpha) - GP(n, beta)
    g = GQ(n, alpha)
    f == g

    # serialize to change a SymmetricPolynomial into a dict representing its expansion into monomials
    f = GP(n, nu)
    f.serialize()


def test_tableaux():
    # generate tableaux and reverse plane partitions
    from tableaux import Tableau
    mu = (4, 2, 1)

    a = Tableau.semistandard(3, mu)
    b = Tableau.semistandard_setvalued(3, mu)

    c = Tableau.semistandard_shifted(3, mu, diagonal_primes=True)
    d = Tableau.semistandard_shifted_setvalued(4, mu, diagonal_primes=False)

    e = Tableau.semistandard_rpp(3, mu)
    f = Tableau.semistandard_shifted_rpp(3, mu, diagonal_nonprimes=False)

    # serialize to change a Tableau into a dictionary
    s = [t.serialize() for t in c]
