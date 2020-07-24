from partitions import Partition
import itertools
import random

def test_basic():
    m = {(1, 1): 1, (1, 2): 1, (3, 1): 1, (4, 1): 1, (2, 2): 1}
    assert Partition.growth_diagram(m) == [
        [Partition(), Partition(), Partition()],
        [Partition(), Partition(1), Partition(2)],
        [Partition(), Partition(1), Partition(3)],
        [Partition(), Partition(2), Partition(3, 1)],
        [Partition(), Partition(3), Partition(3, 2)],
    ]


def test_symmetric(n=4):
    a = {(i, j) for i in range(1, n + 1) for j in range(1, i + 1)}
    for k in range(len(a)):
        for _subset in itertools.combinations(a, k):
            subset = set(_subset) | {(j, i) for (i, j) in _subset}
            dictionary = {}
            for (i, j) in subset:
                x = random.randint(1, 10)
                dictionary[(i, j)] = x
                dictionary[(j, i)] = x
            g = Partition.growth_diagram(dictionary, n, n)
            # Partition.print_growth_diagram(g)
            for i in range(n + 1):
                mu = Partition.transpose(g[i][i])
                sigma = sum([dictionary.get((j, j), 0) for j in range(i + 1)])
                emp = sum([abs(v % 2) for v in mu])
                assert emp == sigma
