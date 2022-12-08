from stable.utils import G, GP, GQ, GS, gp, gq, g
from stable.symmetric import SymmetricPolynomial
from stable.tableaux import Partition
from stable.vectors import Vector
from stable.polynomials import beta as BETA # noqa
import pytest


@pytest.mark.slow
def test_staircase_grothendieck_GP_positivity(): # noqa
    r = 6
    for k in range(r):
        delta = tuple(k - i for i in range(k))
        for nu in Partition.all(sum(delta)):
            if not Partition.contains(delta, nu):
                continue
            n = len(delta)
            f = G(n, delta, nu)
            expansion = SymmetricPolynomial.GP_expansion(f)
            normalized = Vector({
                lam: c * BETA**(sum(lam) - sum(delta) + sum(nu))
                for lam, c in expansion.items()
            })
            print('G_{%s/%s}(x_%s) =' % (delta, nu, n), normalized)
            print()
            assert all(c > 0 for c in normalized.values())


@pytest.mark.slow
def test_skew_GQ_positivity(): # noqa
    k = 10
    for mu in Partition.all(k, strict=True):
        for nu in Partition.all(k, strict=True):
            if not Partition.contains(mu, nu):
                continue
            n = len(mu)
            f = GQ(n, mu, nu)
            expansion = SymmetricPolynomial.GQ_expansion(f)
            normalized = Vector({
                lam: c * BETA**(sum(lam) - sum(mu) + sum(nu))
                for lam, c in expansion.items()
            })
            print('GQ_{%s/%s}(x_%s) =' % (mu, nu, n), normalized)
            print()
            assert all(c > 0 for c in normalized.values())


@pytest.mark.slow
def test_skew_GP_positivity(): # noqa
    k = 10
    for mu in Partition.all(k, strict=True):
        for nu in Partition.all(k, strict=True):
            if not Partition.contains(mu, nu):
                continue
            n = len(mu)
            f = GP(n, mu, nu)
            expansion = SymmetricPolynomial.GP_expansion(f)
            normalized = Vector({
                lam: c * BETA**(sum(lam) - sum(mu) + sum(nu))
                for lam, c in expansion.items()
            })
            print('GP_{%s/%s}(x_%s) =' % (mu, nu, n), normalized)
            print()
            assert all(c > 0 for c in normalized.values())


def is_binary_power(i):
    return len(list(filter(lambda x: x == '1', bin(i)))) == 1


def sgn(mu, nu):
    boxes = sorted({
        (i + 1, i + j) for i in range(len(mu)) for j in range(mu[i] + 1, nu[i] + 1)
    })
    sgn = 1
    for i in range(len(boxes)):
        _, j = boxes[i]
        sgn *= -1 if i > 0 and boxes[i - 1][1] == j else 1
    return sgn


def zero_one_tuples(n):
    if n == 0:
        yield ()
        return
    for i in range(2**n):
        ans = tuple(map(int, list(bin(i))[2:]))
        yield (n - len(ans)) * (0,) + ans


@pytest.mark.slow
def test_GQ_to_GP_expansion(): # noqa
    for mu in Partition.all(25, strict=True):
        print('mu =', mu)
        print()
        print(Partition.printable(mu, shifted=True))
        print()
        n = len(mu)
        q = GQ(n, mu)
        expansion = SymmetricPolynomial.GP_expansion(q)
        normalized = Vector({
            tuple(nu[i] - mu[i] for i in range(len(mu))):
            c * sgn(mu, nu) * BETA**(sum(nu) - sum(mu)) / 2**(len(mu) - sum(nu) + sum(mu))
            for nu, c in expansion.items()
        })
        unsigned = all(c > 0 for c in normalized.values())
        print('  mu =', mu, 'n =', n)
        print('  expansion =', expansion)
        print('  normalized expansion =', normalized)
        assert all(len(nu) == 0 or max(nu) <= 1 for nu in normalized)
        assert all(len(nu) == len(mu) for nu in expansion)
        assert all(Partition.contains(nu, mu) for nu in expansion)
        assert all(c % 2**(len(mu) - sum(nu) + sum(mu)) == 0 for nu, c in expansion.items())
        assert unsigned
        expected = {
            tuple(mu[i] + a[i] for i in range(len(a)))
            for a in zero_one_tuples(len(mu))
            if all(mu[i - 1] + a[i - 1] > mu[i] + a[i] for i in range(1, len(a)))
        }
        print('  expected =', expected)
        assert set(expansion) == expected
        print()
        print()


@pytest.mark.slow
def test_gq_to_gp_expansion(): # noqa
    for mu in Partition.all(15, strict=True):
        print('mu =', mu)
        print()
        print(Partition.printable(mu, shifted=True))
        print()
        n = len(mu)
        q = gq(n, mu)
        expansion = SymmetricPolynomial.gp_expansion(q)
        print('  mu =', mu, 'n =', n)
        print('  expansion =', expansion)
        assert all(len(nu) == len(mu) for nu in expansion)
        assert all(Partition.contains(mu, nu) for nu in expansion)
        # assert all(c % 2**(len(mu) - sum(nu) + sum(mu)) == 0 for nu, c in expansion.items())
        expected = {}
        for a in zero_one_tuples(len(mu)):
            if not all(mu[i - 1] - a[i - 1] > mu[i] - a[i] for i in range(1, len(a))):
                continue
            if not all(mu[i] - a[i] > 0 for i in range(len(a))):
                continue
            nu = Partition.trim(tuple(mu[i] - a[i] for i in range(len(a))))
            coeff = 2**(len(nu) - sum(a)) * sgn(nu, mu) * BETA**sum(a)
            assert coeff != 0
            expected[nu] = coeff
        print('  expected =', expected)
        assert expansion == Vector(expected)
        print()
        print()


def _expansion(n, function, expand, shifted=True, unsigned=True): # noqa
    for mu in Partition.all(n, strict=shifted):
        n = len(mu)
        p = function(n, mu)
        ansion = expand(p)
        if unsigned:
            expansion = {
                nu: coeff * (-1)**abs(sum(mu) - sum(nu))
                for nu, coeff in ansion.items()
            }
        else:
            expansion = ansion
        print('mu =', mu)
        print()
        print(Partition.printable(mu, shifted=shifted))
        print()
        print('  mu =', mu, 'n =', n)
        print('  expansion =', ansion)
        if unsigned:
            print('  unsigned expansion =', expansion)
        print()
        assert all(v > 0 for v in expansion.values())


@pytest.mark.slow
def test_refined_gq_to_gp_expansion(k=12): # noqa
    for mu in Partition.all(k, strict=True):
        print('mu =', mu)
        print()
        print(Partition.printable(mu, shifted=True))
        print()
        n = len(mu)
        q = SymmetricPolynomial._slow_refined_dual_stable_grothendieck_q(n, mu)
        expected = {}
        for a in zero_one_tuples(len(mu)):
            if not all(mu[i - 1] - a[i - 1] > mu[i] - a[i] for i in range(1, len(a))):
                continue
            if not all(mu[i] - a[i] > 0 for i in range(len(a))):
                continue
            nu = Partition.trim(tuple(mu[i] - a[i] for i in range(len(a))))
            coeff = 2**(len(nu) - sum(a)) * sgn(nu, mu) * BETA**sum(a)
            assert coeff != 0
            expected[nu] = coeff
        print('  expected =', expected)
        expected = sum([
            coeff * SymmetricPolynomial._slow_refined_dual_stable_grothendieck_p(n, nu)
            for (nu, coeff) in expected.items()
        ])
        print('           =', expected)
        assert q == expected
        print()
        print()


def _schur_expansion(n, function, shifted=True): # noqa
    _expansion(n, function, SymmetricPolynomial.schur_expansion, shifted)


def _dual_grothendieck_expansion(n, function, shifted=True, unsigned=True): # noqa
    _expansion(n, function, SymmetricPolynomial.dual_grothendieck_expansion, shifted, unsigned)


def test_gp_to_schur_expansion(): # noqa
    _schur_expansion(8, gp)


def test_gp_to_schur_expansion(): # noqa
    _schur_expansion(8, gq)


def test_G_to_g_expansion(): # noqa
    _dual_grothendieck_expansion(8, G, shifted=False, unsigned=False)


def test_GQ_to_g_expansion(): # noqa
    _dual_grothendieck_expansion(8, GQ, unsigned=False)


def test_GP_to_g_expansion(): # noqa
    _dual_grothendieck_expansion(8, GP, unsigned=False)


def test_gq_to_g_expansion(): # noqa
    assert gq(2, (3, 1)) == \
        -4 * (BETA**2) * (-1)**2 * g(2, (2,)) + \
        2 * BETA * (-1) * g(2, (2, 1)) + \
        4 * g(2, (2, 2)) + \
        4 * g(2, (3, 1))

def test_gp_to_g_expansion(): # noqa
    assert gp(2, (3, 1)) == \
        -(BETA**2) * (-1)**2 * g(2, (2,)) + \
        BETA * (-1) * g(2, (2, 1)) + \
        g(2, (2, 2)) + \
        g(2, (3, 1))


@pytest.mark.slow
def test_dewitt_conjecture(n=5):

    def mu(m, k):
        return k * (m,)

    def nu(m, k):
        ans = [m + k - 1]
        while ans[-1] - 2 >= abs(m - k) + 1:
            ans += [ans[-1] - 2]
        return tuple(ans)

    for m in range(1, 10):
        for k in range(1, 10):
            a = mu(m, k)
            b = nu(m, k)
            print(' n =', n)
            print('mu =', a)
            print('nu =', b)
            print()

            if n < len(a) and n < len(b):
                continue

            gs = GS(n, a)
            gq = GQ(n, b)
            print('GS =', gs)
            print()
            print('GQ =', gq)
            print()
            print()
            assert gs == gq
