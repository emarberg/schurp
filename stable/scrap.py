from .tableaux import Partition
from .permutations import *


def skew(w):
    n = w.rank
    f = tuple(i for i in range(1, n + 1) if w(i) > i)
    r = len(f)
    lam = tuple(w(f[-1]) - r - f[i] + i + 1 for i in range(r))
    mu = tuple(w(f[-1]) - r - w(f[i]) + i + 1 for i in range(r))
    print(Partition.printable(lam))
    print()
    print(Partition.printable(mu))
    return lam, mu


def askew(r, mu):
    n = 2 * r
    mu = mu + (r - len(mu)) * (0,)
    b = tuple(r - mu[i] + i + 1 for i in range(r))
    a = sorted(tuple(i for i in range(1, n + 1) if i not in b))
    y = Permutation()
    for i, j in zip(a, b):
        if i != j:
            y *= Permutation.transposition(i, j)
    return y


def test():
    a = (1, 2, 5, 6, 8, 11)
    b = (3, 4, 7, 9, 10, 12)
    y = Permutation()
    for i, j in zip(a, b):
        y *= Permutation.transposition(i, j)
    w = y.get_min_atom()
    return y, w


from polynomials import *
from utils import *
import itertools


    def weak_k_knuth_class(cls, word, length_bound=-1):

        def cab_acb(w):
            for i in range(len(w) - 2):
                c, a, b = w[i:i + 2]
                if a < b < c:
                    yield w[:i] + (a, c, b) + w[i + 2:]

        def acb_cab(w):
            for i in range(len(w) - 2):
                a, c, b = w[i:i + 2]
                if a < b < c:
                    yield w[:i] + (c, a, b) + w[i + 2:]

        def bca_bac(w):
            for i in range(len(w) - 2):
                b, c, a = w[i:i + 2]
                if a < b < c:
                    yield w[:i] + (b, a, c) + w[i + 3:]

        def bac_bca(w):
            for i in range(len(w) - 2):
                b, a, c = w[i:i + 2]
                if a < b < c:
                    yield w[:i] + (b, c, a) + w[i + 3:]

        def aba_bab(w):
            for i in range(len(w) - 2):
                a, b, c = w[i:i + 2]
                if c == a < b:
                    yield w[:i] + (b, a, b) + w[i + 3:]

        def bab_aba(w):
            for i in range(len(w) - 2):
                b, a, c = w[i:i + 2]
                if a < b == c:
                    yield w[:i] + (a, b, a) + w[i + 3:]

        def ab_ba(w):
            if len(w) >= 2 and w[0] != w[1]:
                yield (w[1], w[0]) + w[2:]

        def a_aa(w):
            for i in range(len(w)):
                a = w[i]
                if i == len(w) - 1 or w[i + 1] != a:
                    yield w[:i] + (a,) + w[i + 1:]

        def aa_a(w):
            for i in range(len(w) - 1):
                a, b = w[i], w[i + 1]
                if i == len(w) - 2 or w[i + 2] != a:
                    yield w[:i] + w[i + 1:]

        gens = [cab_acb, acb_cab, bac_bca, bac_bca, aba_bab, bab_aba, ab_ba, a_aa, aa_a]

        assert length_bound == -1 or len(word) <= length_bound
        seen = set()
        add = {word}
        while add:
            for w in add:
                    yield w
                    seen.add(w)
            add = {v for f in gens for w in seen for v in f(w) if (length_bound == -1 or len(v) <= length_bound) and w not in seen}


def test(N):
    """
    Tests that explicit summaation formula for

      K^{(-1)}_{(2,1)}(x_1, x_2, ..., x_N)

    is the same as

      GQ^{(-1)}_{(2,1)}(x_1, x_2, ... ,x_N)

    """
    ans = 0
    for n in range(2, N + 1):
        for i in itertools.combinations(range(1, N + 1), n):
            term = 1 if (n - 3) % 2 == 0 else -1
            for j in range(n):
                x = X(i[j])
                term *= (2 * x - x**2)
            sigma = (n - 1) * (n - 2) // 2
            for j in range(n):
                for k in range(j + 1, n):
                    x = X(i[j])
                    y = X(i[k])
                    sigma -= (x + y - x * y)
            print(i, term * sigma)
            ans += term * sigma
    bns = GQ((2, 1), N).polynomial()
    print(ans)
    print()
    print(bns)
    print()
    print(ans - bns)
    assert ans == bns
