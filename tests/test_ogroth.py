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
from collections import defaultdict, deque
import math
import subprocess


OGROTH_HECKE_ATOMS_CACHE = {}
OGROTH_EXPAND_CACHE = {}


def attempt_factor(p, roots=(1,2)):
    ans = []
    for i in p.variables():
        for c in roots:
            if p.set(i, -c) == 0:
                ans.append(X(i) + c)
                p = p.divide_linear(i, c)
    assert p == p.one() or -p == p.one()
    if p != p.one():
        ans.append(p)
    return ans


def print_factors(p):
    return ''.join(map(lambda x: '(' + str(x) + ')', attempt_factor(p)))


def ogroth_expand(w):
    D = lambda i: lambda f: (f * (1 + X(i + 1))).divided_difference(i)
    swap = lambda i: lambda f: f.substitute(i, X(-1)).substitute(i + 1, X(i)).substitute(-1, X(i + 1))

    assert w.is_involution()
    assert w.is_vexillary()
    
    if w in OGROTH_EXPAND_CACHE:
        return OGROTH_EXPAND_CACHE[w]
    
    if w.is_dominant():
        k = w.number_two_cycles()
        f = X(0)**0
        for i in range(1, k + 1):
            f *= (2 + X(i))
        OGROTH_EXPAND_CACHE[w] = Vector({w: f})
    else:
        n = w.rank
        for i in range(1, n):
            s = Permutation.s_i(i)
            z = s * w * s
            if len(z) == len(w) + 2 and z.is_vexillary():
                bns = ogroth_expand(z)
                ans = {}
                for v, f in bns.items():
                    df = D(i)(f)
                    ans[v] = ans.get(v, 0) + df
                    if v(i) > v(i + 1):
                        sf = swap(i)(f)
                        ans[v] += sf
                        sv = s * v if s * v == v * s else s * v * s
                        assert s * v != v * s
                        ans[sv] = ans.get(sv, 0) + sf
                OGROTH_EXPAND_CACHE[w] = Vector(ans)
                break
    return ogroth_expand(w)

def get_shiftable(w):
    left = [a for a in range(1, w.rank + 1) if a < w(a)]
    segments = []
    for a in left:
        if not segments or a > segments[-1][-1] + 1:
            segments.append([a])
        else:
            segments[-1].append(a)

    crossing_pairs = []
    for seg in segments if w(1) == 1 else segments[1:]:
        for i, b in enumerate(seg):
            a = [a for a in seg[:i] if w(a) < w(b)]
            if a:
                a = a[-1]
                crossing_pairs.append((a, b))

    def get_cycle(s):
        ans = Permutation()
        for seg in segments[1:] if 1 < w(1) else segments:
            cycle = [seg[0] - 1] + [c for c in seg if c in s]
            ans *= Permutation.from_cycles(cycle)
        return ans

    def get_coefficient(s):
        ans = X(0)**0
        for seg in segments:
            for i, a in enumerate(seg):
                if any(b in s for (x, b) in crossing_pairs if x == a):
                    ans *= -1
                elif a not in s:
                    ans *= 2 + X(a)
                else:
                    ans *= 1 + X(a)
        return ans

    mobile = (set(left) - set(segments[0])) if 1 < w(1) else set(left)

    for k in range(len(mobile) + 1):
        for s in itertools.combinations(mobile, k):
            if not any(a not in s and b in s for (a, b) in crossing_pairs):
                sigma = get_cycle(s)
                coeff = get_coefficient(s)
                yield s, sigma, coeff


def _test_ogroth_expand(w, ans):
    expected = {}
    for s, sigma, f in get_shiftable(w):
        v = sigma.inverse() * w * sigma
        expected[v] = f
    # for v, f in expected.items():
    #     print('  ', v.cycle_repr(), ':', print_factors(f))
    # print()
    assert ans == Vector(expected)


def print_ogroth_expand(w):
    ans = ogroth_expand(w)
    print('z =', w.cycle_repr(), 'is dominant:', w.is_dominant(), w.involution_shape())
    print()
    w.print_rothe_diagram(sep='.')
    print()
    w.print_involution_rothe_diagram(sep='.')
    print()
    for v, f in ans.items():
        s = {a for a in range(1, w.rank + 1) if a < w(a) != v(a)}
        print('  ', v.cycle_repr(), ':', print_factors(f), s, v.involution_shape())
    print()
    _test_ogroth_expand(w, ans)


def print_all_ogroth_expand(n, k=None):
    count = 1 
    for w in Permutation.involutions(n):
        w = w.shift(1)
        if k and w.number_two_cycles() != k:
            continue
        if w.is_vexillary():
            print('c =', count)
            print_ogroth_expand(w)
            count += 1



def kbruhat_extended_hecke_graph(w):
    u = set(w.get_extended_hecke_atoms())
    n = w.rank + 1

    j = min([i for (i, j) in w.rothe_diagram() if i == j], default=1) - 1
    k = max([i for (i, j) in w.rothe_diagram() if i == j], default=1)

    e = set()
    for x in u:
        for y in u:
            ab =  jk_bruhat_cover(x, y, j, k)
            if ab is not None:
                a, b = ab
                if a < w(a) or b > w(b):
                    e.add((x, y, ab))
    return (u, e)


def ranked_extended_hecke_graph(w):
    u = set(w.get_extended_hecke_atoms())
    n = w.rank + 1
    e = set()
    for i in range(1, n + 1):
        for x in u:
            y = Permutation.s_i(i) * x
            if len(y) == len(x) + 1 and y in u:
                e.add((x, y, i))
    return (u, e)


def subword_extended_hecke_graph(w):
    n = w.rank + 1

    def undo(x):
        return Permutation(*([1] + [i for i in x.inverse().oneline if i not in [1, n]] + [n])).inverse()
    
    u = set(w.get_extended_hecke_atoms())
    h = set(w.get_involution_hecke_atoms())
    e = set()
    for i in range(1, n + 1):
        # for x in h:
        #     y = Permutation.s_i(i) * x
        #     if len(y) == len(x) + 1 and y in h:
        #         e.add((x, y, i))
        for y in u - h:
            x = undo(y)
            if x in h:
                e.add((x, y, 0))
    return (u, e)


def test_reverse_extended_hecke(n=3):
    assert n > 2
    w = Permutation.longest_element(n - 1).shift(1)
    
    u = set(w.get_extended_hecke_atoms())
    h = set(w.get_involution_hecke_atoms())

    p = 1 + n // 2
    q = (1 + n // 2) if n % 2 == 0 else (2 + n // 2)

    m = w.rank + 1
    expected = set()
    for x in h:
        s = list(x.inverse().oneline)[1:]
        for i in range(m):
            for j in range(i + 1, m):
                a = s[:i]
                b = s[j - 1:]
                if all(t > p for t in a) and all(t < q for t in b):
                    ss = a + [1] + s[i:j - 1] + [m] + b
                    # print(s, i, j, ss)
                    expected.add(Permutation(*ss).inverse())
    assert expected == u


def is_connected(v, e):
    if len(v) == 0:
        return True
    component = {next(iter(v))}
    size = len(component)
    while True:
        component |= {b for (a, b, i) in e if a in component and b in v} | {a for (a, b, i) in e if b in component and a in v}
        if len(component) == size:
            break
        size = len(component)
    return len(component) == len(v)


def test_ranked_graph_edges(n=3):
    dominant_equal = 0
    dominant_total = 0
    dominant_set = set()
    
    vexillary_equal = 0
    vexillary_total = 0
    
    overcount = 0
    totalcount = 0

    invol = [w for w in Permutation.involutions(n)]
    for i, w in enumerate(invol):
        dominant_total += int(w.is_dominant())
        vexillary_total += int(w.is_vexillary())

        left = [a for a in range(1, n + 1) if a < w(a)]
        right = [a for a in range(1, n + 1) if a > w(a)]

        v, e = ranked_extended_hecke_graph(w)
        ogroth = ogroth_hecke_atoms(w)

        j = min([i for (i, j) in w.rothe_diagram() if i == j], default=1) - 1

        a1 = is_connected(v, e)
        a2 = is_connected(ogroth, e)
        b1 = all(x.inverse()(i + 1) in right or x.inverse()(i) in left for (x, y, i) in e)
        b2 = all(x.inverse()(i + 1) in right or x.inverse()(i) in left for (x, y, i) in e if x in ogroth and y in ogroth)
        # for (x, y, i) in e:
        #    print('..', x,y,i,':',x.inverse()(i + 1) in right, x.inverse()(i) in left)
        c = all(
            (x, x.s_i(i) * x, i) in e
            for x in v
            for i in range(max(1, j), w.rank + 1)
            if x.inverse()(i) < x.inverse()(i + 1) and (x.inverse()(i) in left or x.inverse()(i+1) in right)
        )
        d = all(
            (y.s_i(i) * y, y, i) in e
            for y in v
            for i in range(max(1, j), w.rank + 1)
            if y.inverse()(i) > y.inverse()(i + 1) and (y.inverse()(i+1) in left or y.inverse()(i) in right) and w(y.inverse()(i)) != y.inverse()(i+1)
        )
        print('[', i + 1, 'of', len(invol), ']', w.cycle_repr(), ':', c, d, len(ogroth), len(v), 'vex:', w.is_vexillary()) #, ': up has all', c, ': down has all', d)
        assert a1
        assert a2
        assert b1
        assert b2
        if w.is_vexillary():
            assert set(ogroth).issubset(set(v))
            assert set(w.get_involution_hecke_atoms()).issubset(set(ogroth))
        
        if w.is_dominant():
            if len(ogroth) == len(v):
                dominant_equal += 1
            else:
                dominant_set |= {w.shape()}
        vexillary_equal += int(w.is_vexillary() and len(ogroth) == len(v))
        if w.is_vexillary():
            overcount += len(v) - len(ogroth)
            totalcount += len(ogroth)
        
    print()
    print('vexillary:', vexillary_equal, 'of', vexillary_total)
    print(' dominant:', dominant_equal, 'of', dominant_total, dominant_set)
    print('overcount:', overcount, 'of', totalcount)


def jk_bruhat_cover(x, y, j, k):
    if len(y) != len(x) + 1:
        return None
    for a in range(max(1, j), k + 1):
        for b in range(k + 1, x.rank + 2):
            if y == x * x.transposition(a, b):
                return (a, b)


def draw_subword_extended_hecke_graph(n, neato=False, extended=True, labels=True):
    w = Permutation.longest_element(n - 1).shift(1)
    vertices, edges = subword_extended_hecke_graph(w)
    ogroth = ogroth_hecke_atoms(w)
    hecke = w.get_involution_hecke_atoms()
    draw_graph(vertices, edges, ogroth, hecke, neato, extended, labels)


def draw_ranked_extended_hecke_graph(w, neato=False, extended=True, labels=True):
    vertices, edges = ranked_extended_hecke_graph(w)
    ogroth = ogroth_hecke_atoms(w)
    hecke = w.get_involution_hecke_atoms()
    draw_graph(vertices, edges, ogroth, hecke, neato, extended, labels)


def draw_kbruhat_extended_hecke_graph(w, neato=False, extended=True, labels=True):
    vertices, edges = kbruhat_extended_hecke_graph(w)
    ogroth = ogroth_hecke_atoms(w)
    hecke = w.get_involution_hecke_atoms()
    draw_graph(vertices, edges, ogroth, hecke, neato, extended, labels)


def draw_graph(vertices, edges, ogroth, hecke, neato, extended, labels):
    if not extended:
        vertices = ogroth

    s = ['digraph G {']
    s += ['    overlap=false;']
    s += ['    splines=true;']
    s += ['    node [shape=box; style=filled];']

    def prt(x):
        # ans = x.cycle_repr()
        # ans = ''.join(map(str, x.oneline))
        ans = ''.join(map(str, x.inverse().oneline)) + '⁻¹'
        ans += ' : ' + (str(ogroth[x]) if x in ogroth else '0')
        return ans

    for x in vertices:
        s += ['    "%s" [fillcolor=%s];' % (prt(x), 'lightskyblue1' if x in hecke else 'white' if x in ogroth else 'grey85')]

    for x, y, label in edges:
        if x in vertices and y in vertices:
            if labels:
                s += ['    "%s" -> "%s" [label="%s"];' % (prt(x), prt(y), label)]
            else:
                s += ['    "%s" -> "%s";' % (prt(x), prt(y))]

    s += ['}']
    s = '\n'.join(s)

    filename = 'tmp' + str(len(vertices))
    BASE_DIRECTORY = '/Users/emarberg/examples/test/'
    dot_filename = BASE_DIRECTORY + 'graphs/dot/' + '%s.dot' % filename
    png_filename = BASE_DIRECTORY + 'graphs/png/' + '%s.png' % filename
    with open(dot_filename, 'w') as f:
        f.write(s)
    subprocess.run(["neato" if neato else "dot", "-Tpng", dot_filename, "-o", png_filename])
    subprocess.run(["open", png_filename])


def extended_hecke_atoms(w):
    return w.get_extended_hecke_atoms()


def ogroth_hecke_atoms(w):
    if w not in OGROTH_HECKE_ATOMS_CACHE:
        OGROTH_HECKE_ATOMS_CACHE[w] = {w: c for (w, c) in Grothendieck.decompose(InvGrothendieck.get(w)).items()}
    return OGROTH_HECKE_ATOMS_CACHE[w]


def get_paths(target, v, e):
    edges = {}
    for x, y, t in e:
        edges[y] = edges.get(y, []) + [(x, t)]
    q = deque([(target, ())])
    while q:
        (x, path) = q.popleft()
        yield x, tuple(reversed(path))
        for y, t in edges.get(x, []):
            q.append((y, path + (t,)))


def test_upward_k_pieri_tree(n=3, verbose=True):
    for w in Permutation.all(n):
        k = (w.inverse() % w).number_two_cycles()
        f = InvGrothendieck.diagonal_product(k)
        for i in range(1, k + 1):
            f *= X(i)**-1
        expected = Grothendieck.get(w) * f

        v, e = w.upward_k_pieri_tree()
        result = sum([c * Grothendieck.get(x) for x, c in v])
        if verbose:
            print('w =', w)
            print('  ', Grothendieck.decompose(expected))
            print('  ', Grothendieck.decompose(result))
            print('  ', Grothendieck.decompose(expected - result))
        assert result == expected


def test_downward_k_pieri_tree(n=3):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    for mu in mus:
        expected = {Permutation(*w): c for w, c in read_cplusplus_ogroth(mu)}

        k = len(mu)
        z = Permutation.from_involution_shape(*mu)
        n = z.rank
        
        for w in Permutation.all(n + 1):
            v, e = w.downward_k_pieri_tree(k)
            coeff = [c for x, c in v if x.inverse() % x == z]
            assert sum(coeff) == expected.get(w, 0)


def test_partial_ogroth(n=3):
    delta = tuple(range(n - 1, 0, -2))
    mus = [delta] # sorted(Partition.subpartitions(delta, strict=True), key=sum)
    for mu in mus:
        k = len(mu)
        coeff = {}

        start = AltInvGrothendieck.get(Permutation.from_involution_shape(*mu))
        dec = Grothendieck.decompose(start)
        for w in sorted(dec):
            coeff[w] = [dec[w]]
            print(mu, 'i =', 0, 'w =', (' ' if w.rank == n else '') + str(w.inverse()), ':', coeff[w])
        print()

        for i in range(1, k + 1):
            start *= 2 + X(i)
            dec = Grothendieck.decompose(start)
            for w in sorted(dec, key=lambda x:(x.rank, len(x))):
                coeff[w] = coeff.get(w, i * [0]) + [dec[w]]
                print(mu, 'i =', i, 'w =', (' ' if w.rank == n else '') + str(w.inverse()), ':', coeff[w])
                assert dec[w] > 0
                # fails # assert [a for a in coeff[w] if a][0] == 1
                assert all(coeff[w][i] - coeff[w][i - 1] > 0 for i in range(1, len(coeff[w])) if coeff[w][i] != 0 or coeff[w][i - 1] != 0)
            print()


def test_bplus(n=3):
    total = 0
    equal = 0
    good = set()
    t = time.time()
    for z in Permutation.involutions(n):
        if z.is_vexillary():
            total += 1
            h = set(z.get_involution_hecke_atoms())
            d = set(Grothendieck.decompose(InvGrothendieck.get(z)))
            k = max([i for (i, j) in z.rothe_diagram() if i == j], default=1)
            s = set(z.get_extended_hecke_atoms())
            b = d.issubset(s)
            print(total, ':', 'z =', z, '=', z.cycle_repr(), time.time() - t)
            print('  ', 'k =', k, ':', d == s)
            print()
            equal += 1 if d == s else 0
            if d == s:
                good.add(z)
            #print()
            #z.print_involution_rothe_diagram(sep='.')
            #print()
            assert b
            assert h.issubset(d)
            t = time.time()
    print('equal:', equal, 'of', total)
    print()
    return good


def test_unexpected_grothendieck_terms(n=3, verbose=False):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    unexpected = {}
    t = time.time()
    for i, mu in enumerate(mus):
        expected = {Permutation(*w).inverse() for w, c in read_cplusplus_ogroth(mu)}
        k = len(mu)
        z = Permutation.from_involution_shape(*mu)
        for v in z.get_involution_hecke_atoms():
            v = v.inverse()
            for w, forced, prohibited, path in v.inverse_k_pieri_chains(k, k):
                if w not in expected:
                    unexpected[mu] = unexpected.get(mu, []) + [(v, w)]
        print(i + 1, 'of', len(mus), ':', mu, time.time() - t)
        if mu in unexpected:
            print('  ', unexpected[mu])
            print()
        t = time.time()
    print()
    print(unexpected)


def test_k_pieri_chains(n=3, verbose=False):
    delta = tuple(range(n - 1, 0, -2))
    mus = sorted(Partition.subpartitions(delta, strict=True), key=sum)
    unexpected = {}
    for mu in mus:
        expected = {Permutation(*w).inverse(): c for w, c in read_cplusplus_ogroth(mu)}
        k = len(mu)

        mapping = {}
        seen = set()
        z = Permutation.from_involution_shape(*mu)
        for v in z.get_involution_hecke_atoms():
            v = v.inverse()
            # print(mu, v)
            for w, forced, prohibited, path in v.inverse_k_pieri_chains(k, k):
                length = len(path)
                # print('  ', w, length - forced - prohibited, forced, path)
                d = length - forced - prohibited
                seen.add(d)
                mapping[w] = mapping.get(w, []) + [(forced, prohibited, path)]
                if w not in expected:
                    unexpected[mu] = unexpected.get(mu, []) + [(v, w)]
            # print()
            # print('  ', 'seen', seen)
            # print()

        for w in sorted(mapping, key=lambda x: (x.rank, expected.get(x, math.inf))):
            coeff = []
            paths = []
            for forced, prohibited, path in sorted(mapping[w], key=lambda x:x[-1]):
                length = len(path)
                d = length - forced - prohibited
                f = forced
                coeff += [2**(k - length + prohibited) * (-1 if f >= 2 and f % 2 == 0 else 1)]
                paths.append((coeff[-1], tuple(reversed(path))))
                # for p in range(f, min(k, d + f) + 1):
                #     coeff += [2**(k - p) * nchoosek(d, p - f) * (-1 if p >= 2 and p % 2 == 0 else 1)]
            if w.rank == n:
                print(mu, (' ' if w.rank == n else '') + str(w), '= w :', length, forced, prohibited, ':', expected.get(w, 0), '==', max(coeff), '+', sum([c for c in coeff if c > 0]) - max(coeff), '-', sum([-c for c in coeff if c < 0]))
                if verbose:
                    print()
                    for c, path in sorted(paths, key=lambda p: p[1]):
                        print('  ' if c > 0 else ' ', c, ':', path)
                    print()
            assert expected.get(w, 0) == sum(coeff)
    print()
    print(unexpected)


def test_alt_inv_grothendieck(n=5):
    f = {
        w: Vector({x:1 for x in w.get_involution_hecke_atoms()})
        for w in Permutation.involutions(n)
    }
    g = {
        w: Grothendieck.decompose(AltInvGrothendieck.get(w))
        for w in Permutation.involutions(n)
    }
    assert f == g


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


def test_longest_grothendieck(n=3):
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


def test_longest_grothendieck_indices(n=3):
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
