from schubert import GrothendieckB, GrothendieckC, GrothendieckD
from signed import SignedPermutation
from stable.utils import (
    GP,
    GP_expansion, 
    GQ, 
    GQ_expansion,
    SymmetricPolynomial,
)


def test_b(rank=4, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckB.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_c(rank=4, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckC.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GQ_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_d(rank=4, numvars=2):
    for w in SignedPermutation.all(rank, dtype=True):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckD.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_b_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_bc(rank):
        for n in range(numvars + 1):
            f = GrothendieckB.symmetric(n, w)
            g = GP(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g


def test_c_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_bc(rank):
        for n in range(numvars + 1):
            f = GrothendieckC.symmetric(n, w)
            g = GQ(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g


def test_d_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_d(rank):
        for n in range(numvars + 1):
            f = GrothendieckD.symmetric(n, w)
            g = GP(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g