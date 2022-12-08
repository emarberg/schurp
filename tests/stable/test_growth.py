from stable.partitions import Partition
from stable.permutations import Permutation
from stable.tableaux import Tableau
from stable.insertion import InsertionAlgorithm
import itertools
import random


def test_basic():
    m = {(1, 1): 1, (1, 2): 1, (3, 1): 1, (4, 1): 1, (2, 2): 1}
    assert Partition.growth_diagram(m) == [
        [(), (), ()],
        [(), (1,), (2,)],
        [(), (1,), (3,)],
        [(), (2,), (3, 1)],
        [(), (3,), (3, 2)],
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
    g, e, c, r = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)

    gtest = [[[], [], [], [], [], [], [], []], [[], [], [], [1], [1], [1], [1], [1]], [[], [], [1], [2], [2], [2], [2], [2]], [[], [], [1], [2], [2], [2], [3], [3, 1]], [[], [1], [2], [3], [3], [3, 1], [3, 1], [3, 2]]]
    gtest = [[tuple(x) for x in row] for row in gtest]
    assert g == gtest

    etest = [[False, False, False, False, False, False, False, False], [False, False, False, False, False, False, False, False], [False, False, False, True, True, False, False, False], [False, False, False, True, True, False, False, False], [False, False, True, True, True, False, False, True]]
    assert e == etest

    ctest = [[None, None, None, None, None, None, None, None], [None, None, None, None, 1, None, None, None], [None, None, None, None, 2, 1, None, 1], [None, None, None, None, 2, 1, None, None], [None, None, None, None, 3, None, 2, None]]
    assert c == ctest

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    pp, qq = InsertionAlgorithm.orthogonal_hecke(w)
    print(p)
    print(q)
    assert p == pp and q == qq

    w = (4, 5, 1, 2, 3, 4, 6, 5, 6, 4)
    g, e, c, r = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)
    Partition.print_growth_diagram(r)

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    pp, qq = InsertionAlgorithm.orthogonal_hecke(w)
    print(p)
    print(q)
    assert p == pp and q == qq

    w = (1, 3, 2, 5, 6, 4, 3, 5, 2, 4, 5, 6)
    g, e, c, r = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)
    Partition.print_growth_diagram(r)

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    pp, qq = InsertionAlgorithm.orthogonal_hecke(w)
    print(p)
    print(q)
    assert p == pp and q == qq

    w = (3, 2, 1)
    g, e, c, r = Partition.shifted_growth_diagram(w)
    Partition.print_growth_diagram(g)
    Partition.print_growth_diagram(e)
    Partition.print_growth_diagram(c)
    Partition.print_growth_diagram(r)

    p, q = Tableau.from_shifted_growth_diagram(g, e, c)
    pp, qq = InsertionAlgorithm.orthogonal_hecke(w)
    print(p)
    print(q)
    assert p == pp and q == qq


def test_shifted_growth_words(n=4):
    c = True
    reasons = set()
    for a in Permutation.involutions(n):
        for w in a.get_involution_words():
            print(w)
            p, q = InsertionAlgorithm.orthogonal_hecke(w)
            g, e, c, r = Partition.shifted_growth_diagram(w)
            reasons |= set([x for row in r for x in row])
            # Partition.print_growth_diagram(g)
            # Partition.print_growth_diagram(e)
            # Partition.print_growth_diagram(c)
            # print(p)
            # print(q)
            pp, qq = Tableau.from_shifted_growth_diagram(g, e, c)
            # print(pp)
            # print(qq)
            assert p == pp and q == qq
    print(sorted(reasons))
