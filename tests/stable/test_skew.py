from stable.utils import GP, GQ, SymmetricPolynomial, GP_doublebar, GQ_doublebar, beta, gp, gq
from stable.partitions import Partition
from stable.vectors import Vector
from tests.stable.test_conjectures_skew import zero_one_tuples


gp_skew = {}
gq_skew = {}


def substitute(vec, b=None):
    if b is None:
        return vec
    return Vector({key: value.substitute(0, b).constant_term() for (key, value) in vec.items()})


def gp_expand(f, n, mu, nu):
    return substitute(SymmetricPolynomial.GP_expansion(f(n, mu, nu)), -2)


def search(n, k):
    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for nu in partitions:
            if Partition.contains(mu, nu) and len(nu) > 0 and nu[0] != mu[0] and not any(nu[i] > mu[i + 1] for i in range(min(len(nu), len(mu) - 1))):
                target = gp_expand(GQ, n, mu, nu)
                print('\ntarget:', target, '     ', mu, nu)
                print()

                seen = set()
                s = []
                for aa in Partition.successors(mu, strict=True):
                    for a in Partition.successors(aa, strict=True):
                        for bb in Partition.successors(nu, strict=True):
                            for b in Partition.successors(bb, strict=True):
                                if Partition.contains(a, b) and Partition.skew_key(a, b, True) not in seen:
                                    s += [(gp_expand(GP, n, a, b), a, b)]
                                    seen.add(Partition.skew_key(a, b, True))

                display = [(str(x), a, b) for (x, a, b) in s]
                width = 48
                display = [(x if len(x) < width else x[:width] + '...', a, b) for (x, a, b) in display]
                m = max([len(x) for x, _, _ in display])
                for i, (x, a, b) in enumerate(display):
                    print('       ', x, ' ' * (m - len(x)), '       ', a, b, ':', i + 1)

                assert target.is_expressable(*[x[0] for x in s])


def expand_search(n, k):
    if n not in gp_skew:
        gp_skew[n] = {}
    if n not in gq_skew:
        gq_skew[n] = {}

    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for kappa in partitions:
            if Partition.contains(mu, nu) and not (nu and nu[0] == mu[0]):
                key = (nu, mu)
                gp_skew[n][key] = substitute(SymmetricPolynomial.GP_expansion(GP(n, mu, nu)), None)
                gq_skew[n][key] = substitute(SymmetricPolynomial.GQ_expansion(GQ(n, mu, nu)), None)


def overlap(bigger, smaller):
    ans = 0
    skew = Partition.shifted_shape(bigger, smaller)
    for (i, j) in skew:
        if (i + 1, j) in skew:
            ans += 1
    return ans


def cols(bigger, smaller):
    ans = Partition.shifted_shape(bigger) - Partition.shifted_shape(smaller)
    return (-1) ** len({j for (i, j) in ans})


def test_undo_doublebar(n=5, k=5):
    partitions = list(Partition.all(k, strict=True))
    for lam in partitions:
        for mu in partitions:
            if not Partition.contains(lam, mu):
                continue
            lhs = GP(n, lam, mu)
            rhs = sum([(-beta)**(sum(mu) - sum(nu)) * GP_doublebar(n, lam, nu) for nu in partitions if Partition.contains(mu, nu)])
            print('GP:', n, lam, mu)
            assert lhs == rhs
            lhs = GQ(n, lam, mu)
            rhs = sum([(-beta)**(sum(mu) - sum(nu)) * GQ_doublebar(n, lam, nu) for nu in partitions if Partition.contains(mu, nu)])
            print('GQ:', n, lam, mu)
            assert lhs == rhs


def test(n=5, k=5):
    partitions = list(Partition.all(k, strict=True))
    for mu in partitions:
        for kappa in partitions:
            rhs = 0
            expected = {
                tuple(mu[i] + a[i] for i in range(len(a)))
                for a in zero_one_tuples(len(mu))
                if all(mu[i - 1] + a[i - 1] > mu[i] + a[i] for i in range(1, len(a)))
            }
            for lam in expected:
                rhs += 2**(len(mu) - sum(lam) + sum(mu)) * cols(lam, mu) * (-beta) ** (sum(lam) - sum(mu)) * GP_doublebar(n, lam, kappa)
            lhs = 0
            expected = {
                tuple(kappa[i] - a[i] for i in range(len(a)))
                for a in zero_one_tuples(len(kappa))
                if all(kappa[i - 1] - a[i - 1] > kappa[i] - a[i] for i in range(1, len(a))) and all(0 < kappa[i] - a[i] <= (mu[i] if i < len(mu) else 0) for i in range(len(a)))
            }
            for nu in expected:
                lhs += 2**(len(kappa) - sum(kappa) + sum(nu)) * cols(kappa, nu) * (-beta) ** (sum(kappa) - sum(nu)) * GQ_doublebar(n, mu, nu)
            if lhs != 0 or rhs != 0:
                print('n =', n, 'mu =', mu, 'kappa =', kappa)
                # print()
                # print(lhs)
                # print()
            assert lhs == rhs

            if Partition.contains(mu, kappa):
                lhs = 2**(sum(kappa) + len(kappa)) * GQ_doublebar(n, mu, kappa)
                rhs = 0
                expected = {
                    tuple(mu[i] + a[i] for i in range(len(a)))
                    for a in zero_one_tuples(len(mu))
                    if all(mu[i - 1] + a[i - 1] > mu[i] + a[i] for i in range(1, len(a)))
                }
                for nu in Partition.subpartitions(kappa, strict=True):
                    if len(nu) == len(kappa):
                        for lam in expected:
                            rhs += 2**(len(mu) + overlap(kappa, nu) + sum(nu) - sum(lam) + sum(mu)) * cols(lam, mu) * (-beta) ** (sum(lam) - sum(mu) + sum(kappa) - sum(nu)) * GP_doublebar(n, lam, nu)
                print('n =', n, 'mu =', mu, 'kappa =', kappa)
                print()
                print(lhs)
                print(rhs)
                print()
                assert lhs == rhs


def test_dual(n=5, k=5):
    partitions = list(Partition.all(k, strict=True))
    for lam in partitions:
        for kappa in partitions:
            if Partition.contains(lam, kappa):
                lhs = 2**sum(lam) * gq(n, lam, kappa)
                rhs = 0
                expected = {
                    tuple(lam[i] - a[i] for i in range(len(a)))
                    for a in zero_one_tuples(len(lam))
                    if all(lam[i] - a[i] > 0 for i in range(len(a))) and all(lam[i - 1] - a[i - 1] > lam[i] - a[i] for i in range(1, len(a)))
                }
                for mu in expected:
                    for nu in Partition.subpartitions(mu, strict=True):
                        if len(nu) == len(kappa) and Partition.contains(nu, kappa):
                            rhs += 2**(len(mu) - len(nu) + overlap(nu, kappa) + sum(kappa) + sum(mu) - sum(nu)) * cols(lam, mu) * (-beta) ** (sum(lam) - sum(mu) + sum(nu) - sum(kappa)) * gp(n, mu, nu)
                print('n =', n, 'lambda =', lam, 'kappa =', kappa)
                if lhs != rhs:
                    print()
                    print(lhs)
                    print(rhs)
                    print()
                assert lhs == rhs


def test_other_dual(n=5, k=5):
    partitions = list(Partition.all(k, strict=True))
    for lam in partitions:
        for kappa in partitions:
            if Partition.contains(lam, kappa):
                lhs = 2**sum(lam) * gp(n, lam, kappa)
                rhs = 0
                expected = {
                    tuple(lam[i] - a[i] for i in range(len(a)))
                    for a in zero_one_tuples(len(lam))
                    if all(lam[i] - a[i] > 0 for i in range(len(a))) and all(lam[i - 1] - a[i - 1] > lam[i] - a[i] for i in range(1, len(a)))
                }
                for mu in expected:
                    for nu in Partition.subpartitions(mu, strict=True):
                        if len(nu) == len(kappa) and Partition.contains(nu, kappa):
                            rhs += 2**(len(nu) - len(mu) + overlap(lam, mu) + sum(kappa) + sum(mu) - sum(nu)) * cols(nu, kappa) * (-beta) ** (sum(lam) - sum(mu) + sum(nu) - sum(kappa)) * gq(n, mu, nu)
                print('n =', n, 'lambda =', lam, 'kappa =', kappa)
                if lhs != rhs:
                    print()
                    print(lhs)
                    print(rhs)
                    print()
                assert lhs == rhs


def test_inclusion_excluson(k=5):
    partitions = list(Partition.all(k, strict=True))
    for kappa in partitions:
        ans, queue = Vector(), {kappa: 1}
        while queue:
            pop = next(iter(queue))
            coeff = queue[pop]
            ans += Vector({pop: coeff})
            del queue[pop]
            additions = {
                tuple(pop[i] - a[i] for i in range(len(a)))
                for a in zero_one_tuples(len(kappa))
                if not all(a[i] == 0 for i in range(len(a))) and all(pop[i - 1] - a[i - 1] > pop[i] - a[i] for i in range(1, len(a))) and all(0 < pop[i] - a[i] for i in range(len(a)))
            }
            for nu in additions:
                queue[nu] = queue.get(nu, 0) - cols(pop, nu) * coeff
        print('kappa =', kappa)
        print()
        for nu in ans:
            print(Partition.printable(kappa, nu, True))
            print()
            print('coefficient =', ans[nu])
            print()
        print()
        print(ans)
        print()
        print()
        assert all(ans[nu] == 2**overlap(kappa, nu) for nu in ans)
        assert set(ans) == {nu for nu in Partition.subpartitions(kappa, strict=True) if len(nu) == len(kappa)}
#         ans, queue = Vector(), ans.dictionary


def test_dual_inclusion_excluson(k=5):
    partitions = list(Partition.all(k, strict=True))
    for lam in partitions:
        for nu in partitions:
            if not Partition.contains(lam, nu):
                continue
            ans, queue = Vector(), {nu: 1}
            while queue:
                kappa = next(iter(queue))
                coeff = queue[kappa]
                ans += Vector({kappa: coeff})
                del queue[kappa]
                additions = {
                    tuple(kappa[i] + a[i] for i in range(len(a)))
                    for a in zero_one_tuples(len(kappa))
                    if not all(a[i] == 0 for i in range(len(a))) and all(kappa[i - 1] + a[i - 1] > kappa[i] + a[i] for i in range(1, len(a))) and all(lam[i] >= kappa[i] + a[i] for i in range(len(a)))
                }
                for gam in additions:
                    queue[gam] = queue.get(gam, 0) - cols(gam, kappa) * coeff
            print('lambda =', lam, 'nu =', nu)
            print()
            for gam in ans:
                print(Partition.printable(gam, nu, True))
                print()
                print('coefficient =', ans[gam])
                print()
            print()
            print(ans)
            print()
            print()
            assert all(ans[g] == 2**overlap(g, nu) for g in ans)
            assert set(ans) == {g for g in Partition.subpartitions(lam, strict=True) if len(g) == len(nu) and Partition.contains(g, nu)}
