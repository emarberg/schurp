from gt import GTPattern
from stable.partitions import Partition
from stable.utils import schur


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


