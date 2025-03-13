from gt import GTPattern
from permutations import Permutation
from keys import key, get_composition
from stable.partitions import Partition
from stable.utils import schur


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
                    actual += gt.monomial
            if expected != actual:
                print()
                print('expected =', expected)
                print('  actual =', actual)
            assert expected == actual


def test_schur(nrows=3, size=6):
    n = nrows
    mus = list(Partition.generate(size, max_row=n))
    for mu in mus:
        patterns = GTPattern.from_partition(mu, n)
        print('n =', n, 'mu =', mu, 'patterns =', len(patterns))
        expected = schur(n, mu).polynomial()
        actual = 0
        for gt in patterns:
            actual += gt.monomial
        assert expected == actual


