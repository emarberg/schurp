from stable.symmetric import SymmetricPolynomial
from stable.partitions import Partition
from stable.polynomials import beta as BETA
from stable.tableaux import Tableau
from stable.utils import gq, gp, beta, jp, jq, gp_expansion, gq_expansion, jp_expansion, jq_expansion, shifted_ribbon
from stable.vectors import Vector
import itertools


def subsets(s):
    for k in range(len(s) + 1):
        for t in itertools.combinations(s, k):
            yield t


def nchoosek(n, k):
    return len(list(itertools.combinations(range(n), k)))


def test_shifted_setvalued_copieri(n=5): # noqa
    for nu in Partition.all(n, strict=True):
        for lam in Partition.shifted_ribbon_complements(nu):
            n = 1 + max(len(nu), Partition.get(nu, 1))

            print()
            # hat b^{\nu}_{\lambda,(r)}
            ans = Tableau.shifted_setvalued_copieri(nu, lam, diagonal_primes=False)
            for r in ans:
                f = gq(n, lam)
                onerow = Partition.trim((r,))
                g = gq(n, onerow)
                expansion = gq_expansion(f * g)
                expected = expansion[nu]
                expected = expected.substitute(0, 1) if type(expected) != int else expected
                print('b: nu =', nu, 'lambda =', lam, 'r =', r, ':', expected, '==', len(ans[r]))
                if expected != len(ans[r]):
                    print()
                    print(expansion)
                    print(ans)
                assert expected == len(ans[r])

            print()
            # hat c^{\nu}_{\lambda,(r)}
            ans = Tableau.shifted_setvalued_copieri(nu, lam, diagonal_primes=True)
            for r in ans:
                f = gp(n, lam)
                onerow = Partition.trim((r,))
                g = gq(n, onerow)
                expansion = gp_expansion(f * g)
                expected = expansion[nu]
                expected = expected.substitute(0, 1) if type(expected) != int else expected
                print('c: nu =', nu, 'lambda =', lam, 'r =', r, ':', expected, '==', len(ans[r]))
                if expected != len(ans[r]):
                    print()
                    print(expansion)
                    print(ans)
                assert expected == len(ans[r])


def test_shifted_ribbon_jq(n=5):
    for alpha in Partition.compositions(n):
        mu, nu = shifted_ribbon(alpha)
        shape = Partition.shifted_shape(mu, nu)
        p = len([
            (i, j) for (i, j) in shape
            if ((i, j + 1) in shape and (i + 1, j) in shape) or
            ((i, j - 1) in shape and (i - 1, j) in shape)
        ])
        expansion = jq_expansion(jq(1, mu, nu))
        expected = Vector({(n - j,): nchoosek(p, j) * (-beta)**j for j in range(p + 1)})
        print(alpha)
        print(Tableau({a: 0 for a in shape}))
        print("corners =", p)
        print()
        print(expansion)
        print(expected)
        print()
        print()
        assert expansion == expected


def test_shifted_ribbon_jp(n=5):
    for alpha in Partition.compositions(n):
        mu, nu = shifted_ribbon(alpha)
        shape = Partition.shifted_shape(mu, nu)
        p = len([
            (i, j) for (i, j) in shape
            if ((i, j + 1) in shape and (i + 1, j) in shape) or
            ((i, j - 1) in shape and (i - 1, j) in shape) or
            (i == j and (i - 1, j) in shape)
        ])
        expansion = jp_expansion(jp(1, mu, nu))
        expected = Vector({(n - j,): nchoosek(p, j) * (-beta)**j for j in range(p + 1)})
        print(alpha)
        print(Tableau({a: 0 for a in shape}))
        print("corners =", p)
        print()
        print(expansion)
        print(expected)
        print()
        print()
        assert expansion == expected


def gp_pieri(mu, p):
    ans = {}
    shape = Partition.shifted_shape(mu)
    corners = {(i, j) for (i, j) in shape if (i + 1, j) not in shape and (i, j + 1) not in shape}
    nu = ((mu[0] if mu else 0) + p,) + mu
    outer = Partition.shifted_shape(nu, mu)
    for a in subsets(corners):
        for b in subsets(outer):
            c = a + b
            if any(i1 < i2 and j1 < j2 for (i1, j1) in c for (i2, j2) in c):
                continue
            lam = Partition.from_shape(shape | set(c))
            if not Partition.is_strict_partition(lam):
                continue
            if (Partition.shifted_shape(lam) - shape) | set(a) != set(c):
                continue
            free = len([(i, j) for (i, j) in c if (i, j - 1) not in c and (i + 1, j) not in c])
            if len(c) <= p <= len(c) + free:
                diff = p - len(c)
                bet = beta**(sum(mu) + p - sum(lam))
                coeff = nchoosek(free, diff) * bet * 2**(free - diff)
                ans[lam] = ans.get(lam, 0) + coeff
    return ans


def test_gp_pieri(n=4, m=5, l=5):
    for mu in Partition.all(m, strict=True):
        f = gp(n, mu)
        for p in range(1, l + 1):
            g = gq(n, (p,))

            actual = f * g

            expected = 0
            pieri = gp_pieri(mu, p)
            for (lam, c) in pieri.items():
                expected += gp(n, lam) * c

            # print('gp_%s(x_%s)' % (mu, n), '*', 'gp_%s(x_%s)' % (p, n))
            # print('  ', actual - expected)
            # print()
            assert actual == expected

            q = mu[0] if mu else 0
            nu = (q + p,) + mu
            print()
            print(Partition.printable(nu, mu, shifted=True))
            print()
            print('mu =', mu, 'p =', p, 'q =', q)
            print()

            rho1 = (q + p,) + mu[1:]
            rho2 = (q + p - 1,) + mu[1:]

            print(pieri)
            print()

            for lam in pieri:
                co = pieri[lam].substitute(0, 1)
                co = co - 1 if lam in [rho1, rho2] else co
                klg = Tableau.KLG_counts(nu, lam, q)
                # print(lam, co, klg)
                assert co == klg
            print()


def gq_pieri(mu, p):
    ans = {}
    shape = Partition.shifted_shape(mu)
    corners = {(i, j) for (i, j) in shape if (i + 1, j) not in shape and (i, j + 1) not in shape}
    nu = ((mu[0] if mu else 0) + p,) + mu
    outer = Partition.shifted_shape(nu, mu)
    for a in subsets(corners):
        for b in subsets(outer):
            c = a + b
            if any(i1 < i2 and j1 < j2 for (i1, j1) in c for (i2, j2) in c):
                continue
            lam = Partition.from_shape(shape | set(c))
            if not Partition.is_strict_partition(lam):
                continue
            if (Partition.shifted_shape(lam) - shape) | set(a) != set(c):
                continue
            free = len([(i, j) for (i, j) in c if i != j and (i, j - 1) not in c and (i + 1, j) not in c])
            if len(c) <= p <= len(c) + free:
                diff = p - len(c)
                bet = beta**(sum(mu) + p - sum(lam))
                coeff = nchoosek(free, diff) * bet * 2**(free - diff)
                ans[lam] = ans.get(lam, 0) + coeff
    return ans


def test_gq_pieri(n=4, m=5, l=5):
    for mu in Partition.all(m, strict=True):
        f = gq(n, mu)
        for p in range(1, l + 1):
            g = gq(n, (p,))

            actual = f * g

            expected = 0
            pieri = gq_pieri(mu, p)
            for (lam, c) in pieri.items():
                expected += gq(n, lam) * c

            # print('gq_%s(x_%s)' % (mu, n), '*', 'gq_%s(x_%s)' % (p, n))
            # print('  ', actual - expected)
            # print()
            assert actual == expected

            q = mu[0] if mu else 0
            nu = (q + p,) + mu
            print()
            print(Partition.printable(nu, mu, shifted=True))
            print()
            print('mu =', mu, 'p =', p, 'q =', q)
            print()

            rho1 = (q + p,) + mu[1:]
            rho2 = (q + p - 1,) + mu[1:]

            for lam in pieri:
                co = pieri[lam].substitute(0, 1)
                co = co - 1 if mu and lam in [rho1, rho2] else co
                kog = Tableau.KOG_counts(nu, lam, q)
                # print(lam, co, kog)
                assert co == kog
            print()


def test_simple():
    f = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_q(4, (4,), beta=BETA)
    g = SymmetricPolynomial.dual_stable_grothendieck_q(4, (4,))
    assert f.omega_schur_expansion(f) == g.schur_expansion(g)

    f = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_p(4, (4,), beta=BETA)
    g = SymmetricPolynomial.dual_stable_grothendieck_p(4, (4,))
    assert f.omega_schur_expansion(f) == g.schur_expansion(g)


def test_p(n=4):
    counted = 0
    lam = tuple([n - i for i in range(n)])
    for mu in Partition.subpartitions(lam, strict=True):
        for nu in [()]: # Partition.subpartitions(mu, strict=True):
            n = Partition.get(mu, 1)
            f = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_p(n, mu, nu, beta=BETA)
            g = SymmetricPolynomial.dual_stable_grothendieck_p(n, mu, nu)
            ef, eg = f.omega_schur_expansion(f), g.schur_expansion(g)
            print('n =', n, 'mu =', mu)
            print(ef)
            print(eg)
            print()
            assert ef == eg
            counted += 1
            print('tested =', counted)
            print()


def test_q(n=4):
    counted = 0
    lam = tuple([n - i for i in range(n)])
    for mu in Partition.subpartitions(lam, strict=True):
        for nu in [()]: # Partition.subpartitions(mu, strict=True):
            n = Partition.get(mu, 1)
            f = SymmetricPolynomial.slow_transposed_dual_stable_grothendieck_q(n, mu, nu, beta=BETA)
            g = SymmetricPolynomial.dual_stable_grothendieck_q(n, mu, nu)
            ef, eg = f.omega_schur_expansion(f), g.schur_expansion(g)
            print('n =', n, 'mu =', mu)
            print(ef)
            print(eg)
            print()
            assert ef == eg
            counted += 1
            print('tested =', counted)
            print()
