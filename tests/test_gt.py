from gt import GTPattern
from permutations import Permutation
from keys import key, p_key, get_composition, restrict_variables
from tests.test_keys import try_to_decompose_p
from stable.partitions import Partition
from stable.utils import schur, schur_p


def test_tableaux(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        patterns = GTPattern.from_partition(mu, n)
        print('n =', n, 'mu =', mu, 'patterns =', len(patterns))
        for gt in patterns:
            t = gt.tableau()
            nt = GTPattern.from_tableau(t, n)
            assert gt == nt
            assert nt.tableau() == t
        mu_plus_delta = tuple(n - i - 1 + (mu[i] if i < len(mu) else 0) for i in range(n))
        patterns = GTPattern.strict_from_partition(mu, n)
        print('n =', n, 'mu + delta =', mu_plus_delta, 'patterns =', len(patterns))
        for gt in patterns:
            t = gt.shifted_tableau()
            assert gt == GTPattern.from_shifted_tableau(t, n)
            assert GTPattern.from_shifted_tableau(t).shifted_tableau() == t
        print()


def test_schur(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        patterns = GTPattern.from_partition(mu, n)
        print('n =', n, 'mu =', mu, 'patterns =', len(patterns))
        expected = schur(n, mu).polynomial()
        actual = 0
        for gt in patterns:
            actual += gt.weight
        assert expected == actual


def test_key(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        patterns = GTPattern.from_partition(mu, n)
        for w in Permutation.all(n):
            flipped_w = Permutation.longest_element(n) * w
            alpha = get_composition(w, mu)
            print('n =', n, 'mu =', mu, 'w =', w, 'alpha =', alpha, 'patterns =', len(patterns))
            expected = key(alpha)
            actual = 0
            for gt in patterns:
                dkp = gt.dual_kogan_permutation()
                if flipped_w.strong_bruhat_less_equal(dkp):
                    actual += gt.weight
            if expected != actual:
                print()
                print('expected =', expected)
                print('  actual =', actual)
            assert expected == actual


def test_schur_p(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        mu_plus_delta = tuple(n - i - 1 + (mu[i] if i < len(mu) else 0) for i in range(n))
        patterns = GTPattern.strict_from_partition(mu, n)
        print('n =', n, 'mu + delta =', mu_plus_delta, 'patterns =', len(patterns))
        for gt in patterns:
            print(gt)
            print(gt.weight)
        expected = schur_p(n, mu_plus_delta).polynomial()
        actual = 0
        for gt in patterns:
            actual += gt.weight
        if expected != actual:
            print('  ', expected - actual)
            print()
            print('  ', expected)
            print()
            print('  ', actual)
        assert expected == actual


def test_pkey(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        mu_plus_delta = tuple(n - i - 1 + (mu[i] if i < len(mu) else 0) for i in range(n))
        patterns = GTPattern.strict_from_partition(mu, n)
        for w in Permutation.all(n):
            flipped_w = Permutation.longest_element(n) * w
            alpha = get_composition(w, Partition.skew_symmetric_double(mu_plus_delta))
            print('n =', n, 'mu =', mu, 'w =', w, 'alpha =', alpha, 'patterns =', len(patterns))
            expected = restrict_variables(p_key(alpha), n)
            actual = 0
            for gt in patterns:
                skp = gt.strict_kogan_permutation()
                if flipped_w.strong_bruhat_less_equal(skp):
                    actual += gt.weight
            # if expected != actual:
            #     print()
            #     print('  ', expected)
            #     print()
            #     print('  ', actual)
            #     print()
            #     print('  ', actual - expected)
            #     input('')
            print('  ', actual == expected)
            print('  ', try_to_decompose_p(actual))


