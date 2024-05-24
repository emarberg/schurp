from partitions import Partition
from ogroth import (
    grothendieck_transitions,
    grothendieck_double_transitions,
    ogroth_expansion,
)
from schubert import Grothendieck, AltInvGrothendieck, InvGrothendieck, Permutation, X
from vectors import Vector
from stable.tableaux import nchoosek
import itertools
import time
from collections import defaultdict


def test_k_pieri_chains(n=3):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    
    for mu in [delta]:#mus:
        expected = {Permutation(*w).inverse(): c for w, c in read_cplusplus_ogroth(mu)}
        k = len(mu)

        mapping = {}
        seen = set()
        z = Permutation.from_involution_shape(*mu)
        for v in z.get_involution_hecke_atoms():
            v = v.inverse()
            print(mu, v)
            for w, forced, prohibited, path in v.inverse_k_pieri_chains(k, k):
                length = len(path)
                print('  ', w, length - forced - prohibited, forced, path)
                d = length - forced - prohibited
                seen.add(d)
                mapping[w] = mapping.get(w, []) + [(forced, prohibited, length)]
                assert w in expected
            print()
            print('  ', 'seen', seen)
            print()

        for w in sorted(mapping, key=lambda x: (x.rank, expected[x])):
            coeff = []
            for forced, prohibited, length in sorted(mapping[w], key=lambda x:-x[-1]):
                d = length - forced - prohibited
                f = forced
                for p in range(f, min(k, d + f) + 1):
                    coeff += [2**(k - p) * nchoosek(d, p - f) * (-1 if p >= 2 and p % 2 == 0 else 1)]
            print(mu, (' ' if w.rank ==n else '') + str(w), coeff, '==', expected[w])
            assert expected[w] == sum(coeff)
        print()
    print()


def test_alt_inv_grothendieck(n=5):
    f = {
        w: Vector({x:2**w.number_two_cycles() for x in w.get_involution_hecke_atoms()})
        for w in Permutation.involutions(n)
    }
    g = {
        w: Grothendieck.decompose(AltInvGrothendieck.get(w))
        for w in Permutation.involutions(n)
    }
    assert True not in {x.is_vexillary() for x in f if f[x] != g[x]}
    assert False not in {x.is_vexillary() for x in f if f[x] == g[x]}


def chinese_class(w, n):
    def span(w, seen):
        if len(w) == n:
            for i in range(n // 2, len(w)):
                if max((0,) + w[i:]) <= n // 2:
                    v = w[:i] + (n + 1,) + w[i:]
                    if v not in seen:
                        seen.add(v)
                        yield v

        if len(w) == n + 1:
            v = tuple(i for i in w if i <= n)
            if v not in seen:
                seen.add(v)
                yield v

        for i in range(len(w) - 2):
            c, a, b = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (b, c, a) + w[i + 3:],
                    w[:i] + (c, b, a) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

            b, c, a = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (c, a, b) + w[i + 3:],
                    w[:i] + (c, b, a) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

            c, b, a = w[i: i + 3]
            if a < b < c and c != n + 1:
                for v in [
                    w[:i] + (b, c, a) + w[i + 3:],
                    w[:i] + (c, a, b) + w[i + 3:],
                ]:
                    if v not in seen:
                        seen.add(v)
                        yield v

    seen = {w}
    add = {w}
    while add:
        nextadd = set()
        for w in add:
            yield w
            nextadd |= set(span(w, seen))
        add = nextadd


def test_longest_grothendieck(n):
    w0 = Permutation.longest_element(n)
    s = InvGrothendieck.top(w0)
    m = w0.involution_length()
    d = Grothendieck.decompose(s)
    f = {
        tuple(w.inverse().oneline): c * X(0)**(w.length() - m) for w, c in d.dictionary.items()
    }
    a = set(f)
    classes = []
    for w in a:
        if any(w in cl for cl in classes):
            continue
        c = set(chinese_class(w, n))
        classes.append(c)
        assert c.issubset(a)
    assert len(classes) == 1
    print()
    d = defaultdict(int)
    a = sorted(a, key=lambda x: (f[x].substitute(0, 1), f[x].degree(), f[x], x))
    for w in a:
        if len(w) > n:
            continue
        d[f[w]] += 1
        print('  ', w, ':', f[w].set(0, 1))
    print()
    print(d)
    print()
    return a


def test_longest_grothendieck_indices(n):
    w = Permutation.longest_element(n)
    mu = tuple(range(n - 1, 0, -2))
    assert mu == w.involution_shape().tuple()
    tup = tuple(w.get_min_atom().inverse().oneline)
    ans = {tuple(Permutation(*w).inverse().oneline) for w in chinese_class(tup, n)}
    bns = {w for w,_ in read_cplusplus_ogroth(mu)}
    assert ans == bns


def read_cplusplus(mu, directory):
    assert directory in ['ogroth/', 'spgroth/']
    DIRECTORY = "/Users/emarberg/examples/test/" + directory
    file = DIRECTORY + "(" + "".join([str(a) + "," for a in mu])  + ").txt"
    ans = []
    with open(file, 'r') as f:
        for line in f.readlines():
            c, w = line.split(',', 1)
            ans.append((eval(w), int(c)))
    return sorted(ans)


def read_cplusplus_ogroth(mu):
    return read_cplusplus(mu, 'ogroth/')


def read_cplusplus_spgroth(mu):
    return read_cplusplus(mu, 'spgroth/')


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


def test_cplusplus_spgroth(n=6, verbose=False):
    delta = tuple(range(n - 2, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    t0 = t1 = time.time()
    for i, mu in enumerate(mus):
        z = Permutation.from_fpf_involution_shape(*mu)
        ans = sorted([(tuple(w.oneline), 1) for w in z.get_symplectic_hecke_atoms()])
        bns = read_cplusplus_spgroth(mu)

        if verbose:
            print()
            print('mu =', mu)
            print()
            for z, c in ans:
                print('   ', c, '*', z)
            print()
        assert ans == bns
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
