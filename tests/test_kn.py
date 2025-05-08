from schubert import GrothendieckB, GrothendieckC, GrothendieckD
from signed import SignedPermutation
from stable.utils import GP, GQ


def test_b_grassmannian(rank=5, numvars=4):
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


def test_c_grassmannian(rank=5, numvars=4):
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


def test_d_grassmannian(rank=5, numvars=4):
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