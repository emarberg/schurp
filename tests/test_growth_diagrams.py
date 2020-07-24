from partitions import Partition
import itertools


def test_basic():
    m = {(1, 1): 1, (1, 2): 1, (3, 1): 1, (4, 1): 1, (2, 2): 1}
    assert Partition.growth_diagram(m) == [
        [Partition(), Partition(), Partition()],
        [Partition(), Partition(1), Partition(2)],
        [Partition(), Partition(1), Partition(3)],
        [Partition(), Partition(2), Partition(3, 1)],
        [Partition(), Partition(3), Partition(3, 2)],
    ]


def test_symmetric(n=3):
    a = {(i, j) for i in range(1, n + 1) for j in range(1, i)}
    for k in range(len(a)):
        for _subset in itertools.combinations(a, k):
            subset = set(_subset) | {(j, i) for (i, j) in _subset}
            Partition.print_growth_diagram(subset)
