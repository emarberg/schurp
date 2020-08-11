from partitions import Partition
from permutations import Permutation
from tableaux import Tableau
from words import Word
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


def test_shifted_growth_diagram():
    w = (4, 2, 1, 1, 2, 3, 2)
    g, e, c = Partition.shifted_growth_diagram(w)

    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)

    gtest = [[[], [], [], [], [], [], [], []], [[], [], [], [1], [1], [1], [1], [1]], [[], [], [1], [2], [2], [2], [2], [2]], [[], [], [1], [2], [2], [2], [3], [3, 1]], [[], [1], [2], [3], [3], [3, 1], [3, 1], [3, 2]]]
    gtest = [[Partition(*x) for x in row] for row in gtest]
    assert g == gtest

    etest = [[False, False, False, False, False, False, False, False], [False, False, False, False, False, False, False, False], [False, False, False, True, True, False, False, False], [False, False, False, True, True, False, False, False], [False, False, True, True, True, False, False, True]]
    assert e == etest

    ctest = [[None, None, None, None, None, None, None, None], [None, None, None, None, 1, None, None, None], [None, None, None, None, 2, 1, None, 1], [None, None, None, None, 2, 1, None, None], [None, None, None, None, 3, None, 2, None]]
    assert c == ctest

    w = (4, 5, 1, 2, 3, 4, 6, 5, 6, 4)
    g, e, c = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    print(p)
    print(q)
    pp, qq = Word(*w).involution_insert()
    assert p == pp and q == qq

    w = (1, 3, 2, 5, 6, 4, 3, 5, 2, 4, 5, 6)
    g, e, c = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    print(p)
    print(q)
    pp, qq = Word(*w).involution_insert()
    assert p == pp and q == qq


def test_shifted_growth_words(n=5):
    for a in Permutation.involutions(n):
        for w in a.get_involution_words():
            p, q = Word(*w).involution_insert()
            g, e, c = Partition.shifted_growth_diagram(w)
            pp, qq = Tableau.from_shifted_growth_diagram(g, e, c)
            # print(w)
            # print(p)
            # print(pp)
            # print(q)
            # print(qq)
            # Partition.print_growth_diagram(g)
            # Partition.print_growth_diagram(e)
            # Partition.print_growth_diagram(c)
            assert p == pp and q == qq
