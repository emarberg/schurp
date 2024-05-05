from partitions import Partition
from ogroth import (
    grothendieck_transitions,
    grothendieck_double_transitions,
    ogroth_expansion,
)
from schubert import Grothendieck, InvGrothendieck, Permutation
import itertools
import time


def read_cplusplus_ogroth(mu):
    DIRECTORY = "/Users/emarberg/examples/ogroth/"
    file = DIRECTORY + "(" + "".join([str(a) + "," for a in mu])  + ").txt"
    ans = []
    with open(file, 'r') as f:
        for line in f.readlines():
            c, w = line.split(',', 1)
            ans.append((eval(w), int(c)))
    return sorted(ans)


def test_cplusplus_ogroth(n=5, verbose=False):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    t0 = t1 = time.time()
    for i, mu in enumerate(mus):
        ans = sorted(ogroth_expansion(mu))
        bns = Grothendieck.decompose(InvGrothendieck.get(Permutation.from_involution_shape(*mu)))
        bns = sorted([(tuple(x.oneline), bns.dictionary[x]) for x in bns.dictionary])
        cns = read_cplusplus_ogroth(mu)

        if verbose:
            print()
            print('mu =', mu)
            print()
            for z, c in ans:
                print('   ', c, '*', z)
            print()
            for z, c in bns:
                print('   ', c, '*', z)
            print()
            for z, c in cns:
                print('   ', c, '*', z)
            print()
        assert ans == bns
        assert ans == cns
        ###
        print('  #', i + 1, 'of', len(mus), 'mu =', mu, ':', time.time() - t1)
        t1 = time.time()
        ###
    print()
    print('n =', n, 'time =', time.time() - t0)


def test_grothendieck_transitions():
    w = (1, 3, 4, 5, 2)
    j = 3
    assert set(grothendieck_transitions(w, j)) == {
        ((1, 3, 4, 5, 2), 1),
        ((1, 3, 5, 4, 2), 1),
        ((1, 4, 3, 5, 2), -1),
        ((1, 4, 5, 3, 2), -1),
        ((3, 4, 1, 5, 2), 1),
        ((3, 4, 5, 1, 2), 1),
        ((3, 4, 2, 5, 1), 1),
        ((3, 4, 5, 2, 1), 1)
    }


def test_all(n=7, verbose=False):
    total = 1
    for i in range(1, n + 1):
        total *= i
    t = total // 100
    ###
    t0 = t1 = time.time()
    w0 = tuple(i + 1 for i in range(n))
    for i, w in enumerate(itertools.permutations(w0)):
        for j in range(1, n + 2):
            ans = grothendieck_transitions(w, j)

            if verbose:
                print()
                print('w =', w, 'j =', j)
                print()
                for y, c in ans:
                    print('  ', c, '*', y)
                print()
        ###
        if n > 7 and (i + 1) % t == 0:
            print('  ', (i + 1) // t, '%', time.time() - t1)
            t1 = time.time()
        ###
    print('n =', n, 'time =', time.time() - t0)


def test_double_all(n=7, verbose=False):
    total = 1
    for i in range(1, n + 1):
        total *= i
    t = total // 100
    ###
    t0 = t1 = time.time()
    w0 = tuple(i + 1 for i in range(n))
    for i, w in enumerate(itertools.permutations(w0)):
        for j in range(1, n + 2):
            for k in range(j, n + 2):
                ans = grothendieck_double_transitions(w, j, k)

                if verbose:
                    print()
                    print('w =', w, 'j =', j, 'k =', k)
                    print()
                    for y, c in ans:
                        print('  ', c, '*', y)
                    print()

        ###
        if n > 6 and (i + 1) % t == 0:
            print('  ', (i + 1) // t, '%', time.time() - t1)
            t1 = time.time()
        ###
    print('n =', n, 'time =', time.time() - t0)


GT_CACHE = {}

def test_ogroth_expansion(n=6, gtcheck=True, verbose=False):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    ###
    t0 = t1 = time.time()
    for i, mu in enumerate(mus):
        ans = ogroth_expansion(mu)
        if gtcheck:
            if mu not in GT_CACHE:
                bns = Grothendieck.decompose(InvGrothendieck.get(Permutation.from_involution_shape(*mu)))
                bns = sorted([(tuple(x.oneline), bns.dictionary[x]) for x in bns.dictionary])
                GT_CACHE[mu] = bns
            assert sorted(ans) == GT_CACHE[mu]
        if verbose:
            print()
            print('mu =', mu)
            print()
            for z, c in sorted(ans):
                print('   ', c, '*', z)
            print()
        ###
        print('  #', i + 1, 'of', len(mus), 'mu =', mu, ':', time.time() - t1)
        t1 = time.time()
        ###
    print()
    print('n =', n, 'time =', time.time() - t0)
