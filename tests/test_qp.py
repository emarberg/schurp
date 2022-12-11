from qp import QPWGraph, QPModule, gelfand_d_printer, b_print
import polynomials
from polynomials import q_coeff
from permutations import Permutation
from signed import SignedPermutation
from even import EvenSignedPermutation
from tableaux import Tableau
from qp_utils import rsk, gelfand_rsk
import random
import pytest


def test_gelfand_recurrences():
    for sgn in [True, False]:
        modules = []
        for n in range(5):
            for k in range(0, n + 2, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_a, QPModule.create_gelfand_a)
                modules.append(w)
        for n in range(5):
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_bc, QPModule.create_gelfand_bc)
                modules.append(w)
        for n in range(2, 5):
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_d, QPModule.create_gelfand_d)
                modules.append(w)
        for w in modules:
            m = w.qpmodule
            print(m.family, m.rank, m.layer)

            def p(i, j):
                return w.get_cbasis_polynomial(i, j, False)

            def des(i):
                if sgn:
                    return set(m.weak_ascents(i)) | set(m.strict_descents(i))
                else:
                    return set(m.weak_descents(i)) | set(m.strict_descents(i))

            def asc(i):
                if sgn:
                    return set(m.weak_descents(i)) | set(m.strict_ascents(i))
                else:
                    return set(m.weak_ascents(i)) | set(m.strict_ascents(i))

            def sigma(i, j, s):
                ans = 0
                sj = m.operate(j, s)
                for k in range(i, sj):
                    if s in des(k):
                        ans += polynomials.q(m.height(j) - m.height(k)) * w.get_cbasis_leading(k, sj) * p(i, k)
                return ans

            for t in m:
                for u in m:
                    for s in des(t) | asc(t):
                        su = m.operate(u, s)
                        st = m.operate(t, s)

                        if s in des(t):
                            assert p(u, t) == p(su, t)
                        if s in m.strict_descents(t) and s in des(u):
                            assert p(u, t) == p(u, st) * polynomials.q(2) + p(su, st) - sigma(u, t, s)
                        if s in m.strict_descents(t) and s in m.strict_ascents(u):
                            assert p(u, t) == p(u, st) + p(su, st) * polynomials.q(2) - sigma(u, t, s)
                        if s in des(t) and s in (m.weak_descents(u) if sgn else m.weak_ascents(u)):
                            assert p(u, t) == 0


def test_a_descents(n=6):
    def elems(n):
        ans = set()
        for k in range(0, n, 2):
            try:
                m = QPModule.read_gelfand_a(n - 1, k // 2)
            except FileNotFoundError:
                m = QPModule.create_gelfand_a(n - 1, k // 2)
                m.write()
            for i in m:
                oneline = m.permutation(i)
                oneline += tuple(i - (-1) ** i for i in range(len(oneline) + 1, 2 * n + 1))
                ans.add(Permutation(*oneline))
        return ans

    def varphi(v):
        a = [i for i in range(1, n + 1) if abs(v(i)) <= n]

        forward = {i + 1: a[i] for i in range(len(a))}
        backward = {}
        for i in range(len(a)):
            backward[a[i]] = i + 1
            backward[-a[i]] = -i - 1

        b = [i for i in range(-n, n + 1) if v(i) > n]
        w = Permutation(*(a + b))
        z = Permutation(*[backward[v(forward[i])] for i in range(1, len(a) + 1)])
        return w, z

    # test descents
    s = Permutation.s_i

    def is_min_coset_rep(w, z):
        k = z.rank
        return all(w(i) < w(i + 1) for i in range(1, k)) and all(w(i) < w(i + 1) for i in range(k + 1, n))

    def sigma(w, z, i):
        a = s(i)
        if is_min_coset_rep(a * w, z):
            return 0
        b = w**-1 * a * w
        if b * z * b != z:
            return 0
        if b in {s(i) for i in range(1, z.rank)}:
            return -1
        return 1

    def height(w, z, i):
        a = s(i)
        if is_min_coset_rep(a * w, z):
            return len(a * w) - len(w)
        b = w**-1 * a * w
        return (len(b * z * b) - len(z)) // 2

    for v in elems(n):
        w, z = varphi(v)

        # des^=
        vdes_weak = {i for i in range(1, n) if s(i) * v * s(i) == v}
        wdes_weak = {i for i in range(1, n) if sigma(w, z, i) == -1}
        # print(v, w, z, z.rank, vdes_weak, wdes_weak)
        assert vdes_weak == wdes_weak

        # asc^=
        vasc_weak = {i for i in range(1, n) if v * s(i) * v in {s(i) for i in range(n + 1, 2 * n)}}
        wasc_weak = {i for i in range(1, n) if sigma(w, z, i) == 1}
        assert vasc_weak == wasc_weak

        # des^<
        vdes_strict = {i for i in range(1, n) if len(s(i) * v) < len(v)} - vdes_weak - vasc_weak
        wdes_strict = {i for i in range(1, n) if height(w, z, i) == -1}
        assert vdes_strict == wdes_strict

        # asc^<
        vasc_strict = {i for i in range(1, n) if len(s(i) * v) > len(v)} - vdes_weak - vasc_weak
        wasc_strict = {i for i in range(1, n) if height(w, z, i) == 1}
        assert vasc_strict == wasc_strict


def test_bc_descents(n=5):
    def elems(n):
        ans = set()
        for k in range(0, n + 1, 2):
            try:
                m = QPModule.read_gelfand_bc(n, k // 2)
            except FileNotFoundError:
                m = QPModule.create_gelfand_bc(n, k // 2)
                m.write()
            for i in m:
                oneline = m.permutation(i)
                oneline += tuple(i - (-1) ** i for i in range(len(oneline) + 1, 2 * n + 1))
                ans.add(SignedPermutation(*oneline))
        return ans

    def varphi(v):
        a = [i for i in range(1, n + 1) if abs(v(i)) <= n]

        forward = {i + 1: a[i] for i in range(len(a))}
        backward = {}
        for i in range(len(a)):
            backward[a[i]] = i + 1
            backward[-a[i]] = -i - 1

        b = [i for i in range(-n, n + 1) if v(i) > n]
        w = SignedPermutation(*(a + b))
        z = SignedPermutation(*[backward[v(forward[i])] for i in range(1, len(a) + 1)])
        return w, z

    # test descents
    s = SignedPermutation.s_i

    def is_min_coset_rep(w, z):
        k = z.rank
        if k > 0 and w(1) < 0:
            return False
        return all(w(i) < w(i + 1) for i in range(1, k)) and all(w(i) < w(i + 1) for i in range(k + 1, n))

    def sigma(w, z, i):
        a = s(i, n)
        if is_min_coset_rep(a * w, z):
            return 0
        b = w**-1 * a * w
        y = SignedPermutation(*[z(i) if i <= z.rank else i for i in range(1, n + 1)])
        if b * y * b != y:
            return 0
        if b in {s(i, n) for i in range(z.rank)}:
            return -1
        return 1

    def height(w, z, i):
        a = s(i, n)
        if is_min_coset_rep(a * w, z):
            return len(a * w) - len(w)
        b = w**-1 * a * w
        y = SignedPermutation(*[z(i) if i <= z.rank else i for i in range(1, n + 1)])
        return (len(b * y * b) - len(y)) // 2

    for v in elems(n):
        w, z = varphi(v)

        # des^=
        vdes_weak = {i for i in range(n) if s(i, 2 * n) * v * s(i, 2 * n) == v}
        wdes_weak = {i for i in range(n) if sigma(w, z, i) == -1}
        # print(v, w, z, z.rank, vdes_weak, wdes_weak)
        assert vdes_weak == wdes_weak

        # asc^=
        vasc_weak = {i for i in range(n) if v * s(i, 2 * n) * v in {s(i, 2 * n) for i in range(n + 1, 2 * n)}}
        wasc_weak = {i for i in range(n) if sigma(w, z, i) == 1}
        assert vasc_weak == wasc_weak

        # des^<
        vdes_strict = {i for i in range(n) if len(s(i, 2 * n) * v) < len(v)} - vdes_weak - vasc_weak
        wdes_strict = {i for i in range(n) if height(w, z, i) == -1}
        assert vdes_strict == wdes_strict

        # asc^<
        vasc_strict = {i for i in range(n) if len(s(i, 2 * n) * v) > len(v)} - vdes_weak - vasc_weak
        wasc_strict = {i for i in range(n) if height(w, z, i) == 1}
        assert vasc_strict == wasc_strict


def test_d_descents(n=5):
    def elems(n):
        ans = set()
        for k in range(0, n + 1, 2):
            try:
                m = QPModule.read_gelfand_d(n, k // 2)
            except FileNotFoundError:
                m = QPModule.create_gelfand_d(n, k // 2)
                m.write()
            for i in m:
                oneline = m.permutation(i)
                oneline += tuple(i - (-1) ** i for i in range(len(oneline) + 1, 2 * n + 1))
                ans.add(EvenSignedPermutation(*oneline))
        return ans

    def varphi(v):
        a = [i for i in range(1, n + 1) if abs(v(i)) <= n]

        forward = {i + 1: a[i] for i in range(len(a))}
        backward = {}
        for i in range(len(a)):
            backward[a[i]] = i + 1
            backward[-a[i]] = -i - 1

        b = [i for i in range(-n, n + 1) if v(i) > n]
        e = len([i for i in b if i < 0])
        if e % 2 == 0:
            w = EvenSignedPermutation(*(a + b))
            z = EvenSignedPermutation(*[backward[v(forward[i])] for i in range(1, len(a) + 1)])
        else:
            w = EvenSignedPermutation(*([-a[0]] + a[1:] + b))
            z = EvenSignedPermutation(*[backward[v(forward[i])] for i in range(1, len(a) + 1)]).star()
        return w, z, len(a) // 2

    def iota(v):
        e = n + len([i for i in range(1, n + 1) if abs(v(i)) > n])
        x = EvenSignedPermutation(*([-i for i in v.oneline[:e]] + [i for i in v.oneline[e:]]))
        s = EvenSignedPermutation(*(list(range(1, n + 1)) + list(range(e, n, -1)) + list(range(e + 1, 2 * n + 1))))
        x = s * x * s
        if e % 4 != 0:
            x = x.star()
        return x

    # test descents
    s = EvenSignedPermutation.s_i

    def is_min_coset_rep(w, z):
        k = z.rank
        if k > 0 and abs(w(1)) > w(2):
            return False
        return all(w(i) < w(i + 1) for i in range(2, k)) and all(w(i) < w(i + 1) for i in range(k + 1, n))

    def sigma(w, z, i):
        a = s(i, n)
        if is_min_coset_rep(a * w, z):
            return 0
        b = w**-1 * a * w
        y = EvenSignedPermutation(*[z(i) if i <= z.rank else i for i in range(1, n + 1)])
        if b * y * b != y:
            return 0
        if b in {s(i, n) for i in range(z.rank)}:
            return -1
        return 1

    def height(w, z, i):
        a = s(i, n)
        if is_min_coset_rep(a * w, z):
            return len(a * w) - len(w)
        b = w**-1 * a * w
        y = EvenSignedPermutation(*[z(i) if i <= z.rank else i for i in range(1, n + 1)])
        return (len(b * y * b) - len(y)) // 2

    for v in elems(n):
        w, z, _ = varphi(v)

        # des^=
        vdes_weak = {i for i in range(n) if s(i, 2 * n) * v * s(i, 2 * n) == v}
        wdes_weak = {i for i in range(n) if sigma(w, z, i) == -1}
        # print(v, w, z, z.rank, vdes_weak, wdes_weak)
        assert vdes_weak == wdes_weak

        # asc^=
        vasc_weak = {i for i in range(n) if v * s(i, 2 * n) * v in {s(i, 2 * n) for i in range(n + 1, 2 * n)}}
        wasc_weak = {i for i in range(n) if sigma(w, z, i) == 1}
        assert vasc_weak == wasc_weak

        # des^<
        vdes_strict = {i for i in range(n) if len(s(i, 2 * n) * v) < len(v)} - vdes_weak - vasc_weak
        wdes_strict = {i for i in range(n) if height(w, z, i) == -1}
        assert vdes_strict == wdes_strict

        # asc^<
        vasc_strict = {i for i in range(n) if len(s(i, 2 * n) * v) > len(v)} - vdes_weak - vasc_weak
        wasc_strict = {i for i in range(n) if height(w, z, i) == 1}
        assert vasc_strict == wasc_strict

    # test duality
    dictionary = {}
    for v in elems(n):
        w1, z1, k = varphi(v)
        w2, z2, _ = varphi(iota(v))
        if (n + k) % 2 != 0:
            w2 = w2.star()
            z2 = z2.star()
        if n % 2 != 0:
            z1 = z1.star()
        dictionary[k] = dictionary.get(k, []) + [(w1, z1, w2, z2)]

    for k in sorted(dictionary):
        for w1, z1, w2, z2 in dictionary[k]:
            print('n =', n, '| k =', k, '|', v, '|', w1.inverse() * w2, z1.inverse() * z2)
            u = EvenSignedPermutation(*([-i for i in range(1, 1 + 2 * k)] + list(range(n, 2 * k, -1))))
            assert z2 == z1 * EvenSignedPermutation.longest_element(2 * k)
            assert w2 == w1 * u * EvenSignedPermutation.longest_element(n)
        print()


def read_or_create(rank, layer, sgn, read, create):
    try:
        m = read(rank, layer)
        w = QPWGraph.read(m.get_directory(), sgn=sgn)
        if not w.is_cbasis_computed:
            w.compute_cbasis()
            w.write()
    except FileNotFoundError:
        m = create(rank, layer)
        w = QPWGraph(m, sgn=sgn)
        w.compute_cbasis()
        m.write()
        w.write()
    return w


def read_or_compute_wgraph(w, check=False):
    if not w.is_wgraph_computed:
        w.compute_wgraph()
        w.write()
    if check:
        edges = w.slow_compute_wgraph()
        for x in w.qpmodule:
            if sorted(w.get_wgraph_edges(x, True)) != sorted(edges.get(x, [])):
                print(x)
                print('* ', list(w.get_wgraph_edges(x, True)))
                print('* ', edges.get(x, []))
                print()
                raise Exception
    w.compute_cells()


@pytest.mark.slow
def test_gelfand_a_positivity(n=8):
    # fails for n >= 8
    success = True
    for sgn in [False, True]:
        print('type A, n =', n, 'sgn =', sgn)
        values = set()
        for k in range(0, n + 2, 2):
            w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_a, QPModule.create_gelfand_a)
            read_or_compute_wgraph(w)
            for i in w.qpmodule:
                for j, mu in w.get_wgraph_edges(i, True):
                    values.add(mu)
                    if mu < 0:
                        print('k =', k, (i, j), values)
                        print(w.get_cbasis_polynomial(j, i))
                        print(list(w.get_wgraph_edges(i, True)))
                        success = False
        print('* success:', values)
        print()
    assert success


@pytest.mark.slow
def test_gelfand_bc_positivity(n=4):
    # fails n >= 4
    success = True
    for sgn in [False, True]:
        print('type BC, n =', n, 'sgn =', sgn)
        values = set()
        for k in range(0, n + 1, 2):
            w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_bc, QPModule.create_gelfand_bc)
            read_or_compute_wgraph(w)
            for i in w.qpmodule:
                for j, mu in w.get_wgraph_edges(i, True):
                    values.add(mu)
                    if mu < 0:
                        print('k =', k, (i, j), values)
                        print(w.get_cbasis_polynomial(j, i))
                        print(list(w.get_wgraph_edges(i, True)))
                        success = False
        print('* success:', values)
        print()
    assert success


@pytest.mark.slow
def test_gelfand_d_positivity(n=7):
    # fails n >= 7
    success = True
    for sgn in [False, True]:
        print('type D, n =', n, 'sgn =', sgn)
        values = set()
        for k in range(0, n + 1, 2):
            w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_d, QPModule.create_gelfand_d)
            read_or_compute_wgraph(w)
            for i in w.qpmodule:
                for j, mu in w.get_wgraph_edges(i, True):
                    values.add(mu)
                    if mu < 0:
                        print('k =', k, (i, j), values)
                        print(w.get_cbasis_polynomial(j, i))
                        print(list(w.get_wgraph_edges(i, True)))
                        success = False
        print('* success:', values)
        print()
    assert success


def b_toggle(n):
    def toggle(c):
        x = SignedPermutation(*[-a for a in c])
        s = SignedPermutation(*(list(range(1, n + 1)) + list(range(x.rank, n, -1))))
        return tuple((s * x * s).oneline)
    return toggle


def test_gelfand_bc_duality(nn=4):
    read, create = QPModule.read_gelfand_bc, QPModule.create_gelfand_bc
    for n in range(2, nn + 1):
        toggle = b_toggle(n)
        for k in range(0, n + 1, 2):
            w = read_or_create(n, k // 2, True, read, create)
            read_or_compute_wgraph(w)

            v = read_or_create(n, k // 2, False, read, create)
            read_or_compute_wgraph(v)

            duality, index, jndex = {}, {}, {}
            for i in w.qpmodule:
                index[w.permutation(i)] = i
            for j in v.qpmodule:
                jndex[v.permutation(j)] = j
            for c in index:
                duality[c] = jndex[toggle(c)]

            print('n =', n, 'k =', k)
            for x, i in index.items():
                ii = duality[x]

                asc_i = set(w.qpmodule.strict_ascents(i)) | set(w.qpmodule.weak_ascents(i))
                asc_ii = set(v.qpmodule.strict_ascents(ii)) | set(v.qpmodule.weak_descents(ii))
                gen = set(range(n))
                print(asc_i, gen - asc_ii)
                assert asc_i == gen - asc_ii

                for j, mu in w.get_wgraph_edges(i, True):
                    y = w.permutation(j)
                    jj = duality[y]
                    assert (ii, mu) in v.get_wgraph_edges(jj, True)
            print('* success\n')

            # for y, j in index.items():
            #     for x, i in index.items():
            #         summation = []
            #         for z, k in index.items():
            #             f = w.get_cbasis_polynomial(i, k)
            #             g = v.get_cbasis_polynomial(duality[y], duality[z])
            #             c = (-1) ** (w.qpmodule.height(i) + w.qpmodule.height(k))
            #             term = c * f * g
            #             if term != 0:
            #                 summation += [term]
            #         print('  ', i, j, ':', summation)
            #         assert sum(summation) in [0, 1]
            #         assert (sum(summation) == 1) == (i == j)

            # for y, j in index.items():
            #     for x, i in index.items():
            #         if i == j:
            #             continue
            #         f = w.get_cbasis_polynomial(i, j)
            #         ii = duality[x]
            #         jj = duality[y]
            #         g = v.get_cbasis_polynomial(jj, ii)
            #         if f != g and (q_coeff(f, -1) != 0 or q_coeff(g, -1) != 0):
            #             print('  ', i, j, ':', f, '<->', b_print(w)(i), b_print(w)(j))
            #             print('  ', ii, jj, ':', g, '<->', b_print(w)(ii), b_print(w)(jj))
            #             print()
            #             assert q_coeff(f, -1) == q_coeff(g, -1)


def d_toggle(n):
    def toggle(c):
        x = EvenSignedPermutation(*[-a for a in c])
        s = EvenSignedPermutation(*(list(range(1, n + 1)) + list(range(x.rank, n, -1))))
        x = s * x * s
        if x.rank % 4 != 0:
            x = x.star()
        return tuple(x.oneline)
    return toggle


def test_gelfand_d_duality(nn=4):
    read, create = QPModule.read_gelfand_d, QPModule.create_gelfand_d
    for n in range(2, nn + 1):
        toggle = d_toggle(n)
        for k in range(0, n + 1, 2):
            w = read_or_create(n, k // 2, True, read, create)
            read_or_compute_wgraph(w)

            v = read_or_create(n, k // 2, False, read, create)
            read_or_compute_wgraph(v)

            print('\nn =', n, 'k =', k)
            duality, index, jndex = {}, {}, {}
            for i in w.qpmodule:
                index[w.permutation(i)] = i
            for j in v.qpmodule:
                jndex[v.permutation(j)] = j
            for c in index:
                duality[c] = jndex[toggle(c)]

            for x, i in index.items():
                ii = duality[x]

                asc_i = set(w.qpmodule.strict_ascents(i)) | set(w.qpmodule.weak_ascents(i))
                asc_ii = set(v.qpmodule.strict_ascents(ii)) | set(v.qpmodule.weak_descents(ii))
                if len(x) % 4 != 0:
                    asc_ii = {(j - 1 if j == 1 else j + 1 if j == 0 else j) for j in asc_ii}
                gen = set(range(n))
                assert asc_i == gen - asc_ii

                for j, mu in w.get_wgraph_edges(i, True):
                    y = w.permutation(j)
                    jj = duality[y]
                    assert (ii, mu) in v.get_wgraph_edges(jj, True)
            print('* success\n')

            # for y, j in index.items():
            #     for x, i in index.items():
            #         summation = []
            #         for z, k in index.items():
            #             f = w.get_cbasis_polynomial(i, k)
            #             g = v.get_cbasis_polynomial(duality[y], duality[z])
            #             c = (-1) ** (w.qpmodule.height(i) + w.qpmodule.height(k))
            #             term = c * f * g
            #             if term != 0:
            #                 summation += [term]
            #         print('  ', i, j, ':', summation)
            #         assert sum(summation) in [0, 1]
            #         assert (sum(summation) == 1) == (i == j)

            # for x, i in index.items():
            #     for y, j in index.items():
            #         if i == j:
            #             continue
            #         f = w.get_cbasis_polynomial(i, j)
            #         ii = duality[x]
            #         jj = duality[y]
            #         xx = v.permutation(ii)
            #         yy = v.permutation(jj)
            #         g = v.get_cbasis_polynomial(jj, ii)
            #         if f != g and (q_coeff(f, -1) != 0 or q_coeff(g, -1) != 0):
            #             print('  ', i, j, ':', f, '<->', x, y)
            #             print('  ', ii, jj, ':', g, '<->', xx, yy)
            #             print()
            #             assert q_coeff(f, -1) == q_coeff(g, -1)


def test_gelfand_rsk():
    assert gelfand_rsk((4, 3, 2, 1), 3, False) == Tableau({(1, 1): 1, (2, 1): 2, (3, 1): 3})


def test_bytes():
    assert QPWGraph._bytes(0) == 1
    assert QPWGraph._bytes(1) == 1
    assert QPWGraph._bytes(7) == 1
    assert QPWGraph._bytes(255) == 1
    assert QPWGraph._bytes(256) == 2
    assert QPWGraph._bytes(256 * 256 - 1) == 2
    assert QPWGraph._bytes(256 * 256) == 3


def test_hecke_cells_a(nn=4):
    for n in range(nn + 1):
        w = read_or_create(n, None, None, QPModule.read_hecke_a, QPModule.create_hecke_a)
        read_or_compute_wgraph(w)
        cells = w.get_cells_as_permutations()
        molecules = w.get_molecules_as_permutations()
        for c in cells:
            r = {rsk(w)[0] for w in c}
            # print(c)
            # print(r)
            # print()
            assert len(r) == 1
        print('n =', n, '::', '#edges =', w.get_wgraph_size(), '#cells =', len(cells), '#molecules =', len(molecules))
        print()
        assert sorted([sorted(c) for c in cells]) == sorted([sorted(c) for c in molecules])


def test_hecke_cells_bc(nn=4):
    for n in range(2, nn + 1):
        w = read_or_create(n, None, None, QPModule.read_hecke_bc, QPModule.create_hecke_bc)
        read_or_compute_wgraph(w)
        cells = w.cells
        molecules = w.molecules
        print('n =', n, '::', '#edges =', w.get_wgraph_size(), '#cells =', len(cells), '#molecules =', len(molecules))
        print('* cells == molecules:', sorted([sorted(c) for c in cells]) == sorted([sorted(c) for c in molecules]))
        print()


def test_hecke_cells_d(nn=4):
    for n in range(2, nn + 1):
        w = read_or_create(n, None, None, QPModule.read_hecke_d, QPModule.create_hecke_d)
        read_or_compute_wgraph(w)
        cells = w.cells
        molecules = w.molecules
        print('n =', n, '::', '#edges =', w.get_wgraph_size(), '#cells =', len(cells), '#molecules =', len(molecules))
        print('* cells == molecules:', sorted([sorted(c) for c in cells]) == sorted([sorted(c) for c in molecules]))
        print()


def test_two_sided_cells_a(nn=4):
    for n in range(nn + 1):
        w = read_or_create(n, None, None, QPModule.read_two_sided_hecke_a, QPModule.create_two_sided_hecke_a)
        read_or_compute_wgraph(w)
        cells = w.cells
        molecules = w.molecules
        seen = set()
        for c in w.get_cells_as_permutations():
            r = {rsk(w)[0].partition() for w in c}
            # print(c)
            # print(r)
            # print()
            assert len(r) == 1
            assert not r.issubset(seen)
            seen |= r
        print('n =', n, '::', '#edges =', w.get_wgraph_size(), '#cells =', len(cells), '#molecules =', len(molecules))
        assert sorted([sorted(c) for c in cells]) == sorted([sorted(c) for c in molecules])

        # cells = [set(c) for c in w.get_cells_as_permutations()]
        #
        # print()
        # for sgn in [True, False]:
        #     for k in range(0, n + 2, 2):
        #         w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_a, QPModule.create_gelfand_a)
        #         read_or_compute_wgraph(w)
        #         gelfand_cells = [{truncate_a(t, n + 1) for t in c} for c in w.get_cells_as_permutations()]
        #         for c in gelfand_cells:
        #             b = all(c & sup == c or c & sup == set() for sup in cells)
        #             print('sgn =', sgn, '::', 'k =', k, '::', b, '::', c)
        #             r = {gelfand_rsk(x, n + 1, sgn=sgn).partition() for x in c}
        #             print()
        # print()
        # print()


def test_two_sided_cells_bc(nn=4):
    for n in range(2, nn + 1):
        w = read_or_create(n, None, None, QPModule.read_two_sided_hecke_bc, QPModule.create_two_sided_hecke_bc)
        read_or_compute_wgraph(w)
        print('rank =', n, '#edges =', w.get_wgraph_size(), '#cells =', len(w.cells), '#molecules =', len(w.molecules))
        # cells = [set(c) for c in w.get_cells_as_permutations()]

        # print()
        # for sgn in [True, False]:
        #     for k in range(0, n + 1, 2):
        #         w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_bc, QPModule.create_gelfand_bc)
        #         read_or_compute_wgraph(w)
        #         gelfand_cells = [{truncate_bc(t, n) for t in c} for c in w.get_cells_as_permutations()]
        #         for c in gelfand_cells:
        #             b = all(c & sup == c or c & sup == set() for sup in cells)
        #             print('sgn =', sgn, '::', 'k =', k, '::', b, '::', c)
        #             print()
        # print()
        # print()


def test_two_sided_cells_d(nn=4):
    for n in range(2, nn + 1):
        w = read_or_create(n, None, None, QPModule.read_two_sided_hecke_d, QPModule.create_two_sided_hecke_d)
        read_or_compute_wgraph(w)
        print('rank =', n, '#edges =', w.get_wgraph_size(), '#cells =', len(w.cells), '#molecules =', len(w.molecules))

        # cells = [set(c) for c in w.get_cells_as_permutations()]

        # print()
        # for sgn in [True, False]:
        #     for k in range(0, n + 1, 2):
        #         w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_d, QPModule.create_gelfand_d)
        #         read_or_compute_wgraph(w)
        #         gelfand_cells = [{truncate_bc(t, n) for t in c} for c in w.get_cells_as_permutations()]
        #         for c in gelfand_cells:
        #             b = all(c & sup == c or c & sup == set() for sup in cells)
        #             print('sgn =', sgn, '::', 'k =', k, '::', b, '::', c)
        #             print()
        # print()
        # print()


def test_gelfand_cells_a(nn=5, s=None):
    for n in range(nn + 1):
        for sgn in ([False, True] if s is None else [s]):
            vertices = []
            cells = []
            molecules = []
            edges = []
            comp = 0
            for k in range(0, n + 2, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_a, QPModule.create_gelfand_a)
                read_or_compute_wgraph(w)
                cells += w.get_cells_as_permutations()
                molecules += w.get_molecules_as_permutations()
                edges += [w.get_wgraph_size()]
                vertices += [w.qpmodule.size]
                comp += w.count_weakly_connected_components()
                # w.print_wgraph()
                # w.print_cells()
            print(
                ('N' if sgn else 'M') + str(n),
                'vertices =', vertices, '(', sum(vertices), ')',
                'edges =', edges, '(', sum(edges), ')',
                'components =', comp,
                'cells =', len(cells),
                'molecules =', len(molecules),
                'weights =', QPWGraph.get_wgraph_weights_gelfand_a(n, sgn),
                'cell weights =', QPWGraph.get_cell_weights_gelfand_a(n, sgn))
            seen = {}
            for c in cells:
                r = {gelfand_rsk(x, n + 1, sgn=sgn).partition() for x in c}
                if not sgn:
                    assert all(gelfand_rsk(x, n + 1, sgn=sgn) == rsk(x)[0].restrict(n + 1) for x in c)
                mu = next(iter(r))
                if len(r) != 1 or mu in seen:
                    print(c)
                    print()
                    print(seen[mu])
                    print()
                    print({gelfand_rsk(x, n + 1, sgn=sgn) for x in c})
                    print()
                    print()
                    print(r)
                    print()
                assert len(r) == 1
                assert not r.issubset(seen)
                seen[mu] = c
            assert sorted([sorted(c) for c in cells]) == sorted([sorted(c) for c in molecules])
        print()


def test_gelfand_cells_bc(nn=5, s=None):
    cellmap = {}
    for n in range(2, nn + 1):
        cellmap[n] = {}
        for sgn in ([False, True] if s is None else [s]):
            cellmap[n][sgn] = set()
            comp = 0
            cells = []
            molecules = []
            edges = []
            vertices = []
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_bc, QPModule.create_gelfand_bc)
                read_or_compute_wgraph(w)
                # print(n, k, len(w.cells))
                cells += w.cells
                molecules += w.molecules
                edges += [w.get_wgraph_size()]
                vertices += [w.qpmodule.size]
                comp += w.count_weakly_connected_components()
                cellmap[n][sgn] |= w.get_cells_as_permutations()
                # w.print_wgraph()
                # w.print_cells()
            print(
                ('N' if sgn else 'M') + str(n),
                'vertices =', vertices, '(', sum(vertices), ')',
                'edges =', edges, '(', sum(edges), ')',
                'components =', comp,
                'cells =', len(cells),
                'molecules =', len(molecules),
                'weights =', QPWGraph.get_wgraph_weights_gelfand_bc(n, sgn),
                'cell weights =', QPWGraph.get_cell_weights_gelfand_bc(n, sgn))
            if (not sgn) in cellmap[n]:
                toggle = b_toggle(n)
                negated = {tuple(sorted(toggle(c) for c in cell)) for cell in cellmap[n][not sgn]}
                if negated != cellmap[n][sgn]:
                    for c in sorted(cellmap[n][sgn], key=len):
                        print('1 ', c)
                    print()
                    for c in sorted(negated, key=len):
                        print('2 ', c)
                assert negated == cellmap[n][sgn]
                # print('* isomorphism checked\n')
        print()


def test_gelfand_cells_d(nn=5, s=None):
    cellmap = {}
    for n in range(2, nn + 1):
        cellmap[n] = {}
        for sgn in ([False, True] if s is None else [s]):
            cellmap[n][sgn] = set()
            comp = 0
            cells = []
            molecules = []
            edges = []
            vertices = []
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_d, QPModule.create_gelfand_d)
                read_or_compute_wgraph(w)
                # print(n, k, len(w.cells))
                cells += w.cells
                molecules += w.molecules
                comp += w.count_weakly_connected_components()
                edges += [w.get_wgraph_size()]
                vertices += [w.qpmodule.size]
                cellmap[n][sgn] |= w.get_cells_as_permutations()
                # w.print_wgraph()
                # w.print_cells()
            print(
                ('N' if sgn else 'M') + str(n),
                'vertices =', vertices, '(', sum(vertices), ')',
                'edges =', edges, '(', sum(edges), ')',
                'components =', comp,
                'cells =', len(cells),
                'molecules =', len(molecules),
                'weights =', QPWGraph.get_wgraph_weights_gelfand_d(n, sgn),
                'cell weights =', QPWGraph.get_cell_weights_gelfand_d(n, sgn))
            if (not sgn) in cellmap[n]:
                toggle = d_toggle(n)
                negated = {tuple(sorted(toggle(c) for c in cell)) for cell in cellmap[n][not sgn]}
                if negated != cellmap[n][sgn]:
                    for c in sorted(cellmap[n][sgn], key=len):
                        print('1 ', c)
                    print()
                    for c in sorted(negated, key=len):
                        print('2 ', c)
                    assert negated == cellmap[n][sgn]
                # print('* isomorphism checked\n')
        print()


def test_qpwgraph(n=6, k=3):
    m = QPModule.create_gelfand_a(n, k)
    w = QPWGraph(m, False)
    w.frame = bytearray(w.size)

    random.seed(12345)
    expected = {}
    for j in m:
        for i in range(j - 1, -1, -1):
            if m.height(i) == m.height(j):
                continue
            e = random.randint(0, 2**(8 * w.nbytes - 1) - 1)
            w.set_cbasis(i, j, e, set_bytes=False)
            expected[(i, j)] = e

    for j in m:
        for i in m:
            if m.height(i) < m.height(j):
                assert m._int(w.get_cbasis(i, j), signed=True) == expected[(i, j)]


def test_qpwgraph_cbasis(n=3):
    m = QPModule.create_hecke_a(n)
    w = QPWGraph(m)
    w.compute_cbasis()
    for j in m:
        for i in m:
            f = w.get_cbasis_polynomial(i, j)
            if f:
                print(m.permutation(i), m.permutation(j), f, w.get_cbasis_leading(i, j), m.height(j) - m.height(i))
            assert f.is_zero() or f.is_positive()
        print()


def test_qpwgraph_speed(n=4):
    m = QPModule.create_hecke_a(n)
    w1 = QPWGraph(m)
    w1.compute_cbasis()
    w2 = QPWGraph(m)
    w2._slowcompute()
    assert w1.frame == w2.frame


def test_eq():
    m = QPModule.create_gelfand_a(1, 0)
    n = QPModule.create_gelfand_a(1, 0)
    assert m == n

    n = QPModule.create_gelfand_a(1, 1)
    assert m != n
    assert len(m) == len(n)


def test_create():
    m = QPModule.create_gelfand_a(0, 0)
    m = QPModule.create_gelfand_a(1, 0)

    m = QPModule.create_gelfand_a(1, 0)
    assert set(m.weak_descents(0)) == set()
    assert set(m.weak_ascents(0)) == {0}

    m = QPModule.create_gelfand_a(1, 1)
    assert set(m.weak_descents(0)) == {0}
    assert set(m.weak_ascents(0)) == set()

    m = QPModule.create_gelfand_a(1, 1)
    m = QPModule.create_gelfand_a(2, 1)


def test_slow_gelfand_a(nin=6):
    for n in range(nin + 1):
        for k in range(0, n + 2, 2):
            m = QPModule.create_gelfand_a(n, k // 2)
            slow = QPModule.slow_create_gelfand_a(n, k // 2)
            print(m)
            print(slow)
            assert m == slow


def test_slow_gelfand_bc(nin=6):
    for n in range(2, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_bc(n, k // 2)
            slow = QPModule.slow_create_gelfand_bc(n, k // 2)
            assert m == slow


def test_slow_gelfand_d(nin=6):
    for n in range(3, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            slow = QPModule.slow_create_gelfand_d(n, k // 2)
            assert m == slow


def _test_wgraph(m, sgn=None):
    w = QPWGraph(m, sgn=sgn)
    w.write()
    w.compute_cbasis(verbose=True)
    w.write()
    read = QPWGraph.read(m.get_directory(), sgn=sgn)
    assert w == read



def test_io_two_sided_hecke(n=3):
    m = QPModule.create_two_sided_hecke_a(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)

    m = QPModule.create_two_sided_hecke_bc(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)

    m = QPModule.create_two_sided_hecke_d(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)


def test_io_hecke(n=3):
    m = QPModule.create_hecke_a(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)

    m = QPModule.create_hecke_bc(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)

    m = QPModule.create_hecke_d(n)
    m.write()
    read = QPModule.read(m.get_directory())
    assert m == read
    _test_wgraph(m)


def test_io_a(n=3):
    for k in range(0, n + 2, 2):
        m = QPModule.create_gelfand_a(n, k // 2)
        m.write()
        read = QPModule.read(m.get_directory())
        assert m == read
        _test_wgraph(m, True)
        _test_wgraph(m, False)


def test_io_bcd(n=3):
    for k in range(0, n + 1, 2):
        m = QPModule.create_gelfand_bc(n, k // 2)
        m.write()
        read = QPModule.read(m.get_directory())
        assert m == read
        _test_wgraph(m, True)
        _test_wgraph(m, False)
    for k in range(0, n + 1, 2):
        m = QPModule.create_gelfand_d(n, k // 2)
        m.write()
        read = QPModule.read(m.get_directory())
        assert m == read
        _test_wgraph(m, True)
        _test_wgraph(m, False)


def test_gelfand_a(nin=5):
    for n in range(nin + 1):
        for k in range(0, n + 2, 2):
            m = QPModule.create_gelfand_a(n, k // 2)
            assert len({m.printer(m.reduced_word(i)) for i in m}) == m.size
            for e in m:
                for i in range(n - 1):
                    for j in range(n - 1):
                        if abs(i - j) > 1:
                            assert m.operate(e, i, j) == m.operate(e, j, i)
                        elif abs(i - j) == 1:
                            assert m.operate(e, i, j, i) == m.operate(e, j, i, j)
                        else:
                            assert m.operate(e, i, i) == e


def test_hecke_gelfand_a(nin=4):
    q = polynomials.q
    for n in range(nin + 1):
        for k in range(0, n + 2, 2):
            m = QPModule.create_gelfand_a(n, k // 2)
            for e in m:
                for i in range(n - 1):
                    for j in range(n - 1):
                        print(n, k // 2, ':', e, i, j)
                        if abs(i - j) > 1:
                            assert m.element(e).operate(i, j) == m.element(e).operate(j, i)
                        elif abs(i - j) == 1:
                            assert m.element(e).operate(i, j, i) == m.element(e).operate(j, i, j)
                        else:
                            f = m.operate(e, i)
                            if e == f:
                                assert m.element(e).operate(i) in [m.element(e) * q(1), m.element(e) * -q(-1)]
                            else:
                                assert m.element(e).operate(i, i) == m.element(e) + m.element(e).operate(i) * (q(1) - q(-1))


def test_gelfand_bc(nin=4):
    for n in range(2, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_bc(n, k // 2)
            assert len({m.printer(m.reduced_word(i)) for i in m}) == m.size
            for e in m:
                for i in range(n):
                    for j in range(n):
                        print(n, k // 2, ':', e, i, j)
                        if {i, j} == {0, 1}:
                            assert m.operate(e, i, j, i, j) == m.operate(e, j, i, j, i)
                        elif abs(i - j) > 1:
                            assert m.operate(e, i, j) == m.operate(e, j, i)
                        elif abs(i - j) == 1:
                            assert m.operate(e, i, j, i) == m.operate(e, j, i, j)
                        else:
                            assert m.operate(e, i, i) == e


def test_hecke_gelfand_bc(nin=3):
    q = polynomials.q
    for n in range(2, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_bc(n, k // 2)
            for e in m:
                for i in range(n):
                    for j in range(n):
                        print(n, k // 2, ':', e, i, j)
                        if {i, j} == {0, 1}:
                            assert m.element(e).operate(i, j, i, j) == m.element(e).operate(j, i, j, i)
                        elif abs(i - j) > 1:
                            assert m.element(e).operate(i, j) == m.element(e).operate(j, i)
                        elif abs(i - j) == 1:
                            assert m.element(e).operate(i, j, i) == m.element(e).operate(j, i, j)
                        else:
                            f = m.operate(e, i)
                            if e == f:
                                assert m.element(e).operate(i) in [m.element(e) * q(1), m.element(e) * -q(-1)]
                            else:
                                assert m.element(e).operate(i, i) == m.element(e) + m.element(e).operate(i) * (q(1) - q(-1))


def test_gelfand_d(nin=4):
    for n in range(2, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            assert len({m.printer(m.reduced_word(i)) for i in m}) == m.size
            for e in m:
                print('n =', n, 'k =', k // 2, ':', 'element =', m.element(e))
                for i in range(n):
                    for j in range(n):
                        if {i, j} == {0, 2} or (abs(i - j) == 1 and {i, j} != {0, 1}):
                            x = m.operate(e, i, j, i)
                            y = m.operate(e, j, i, j)
                            print('  R1:', 'i =', i, 'j =', j, '->', m.element(x), '==', m.element(y))
                            assert x == y
                        elif {i, j} == {0, 1} or abs(i - j) > 1:
                            x = m.operate(e, i, j)
                            y = m.operate(e, j, i)
                            print('  R2:', 'i =', i, 'j =', j, '->', m.element(x), '==', m.element(y))
                            assert x == y
                        elif i == j:
                            x = e
                            y = m.operate(e, i, i)
                            print('  R3', 'i =', i, 'j =', j, '->', m.element(x), '==', m.element(y))
                            assert x == y
                print()


def test_hecke_gelfand_d(nin=4):
    q = polynomials.q
    for n in range(2, nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            for e in m:
                # print('n =', n, 'k =', k // 2, ':', 'element =', m.element(e), 'of', m.size)
                for i in range(n):
                    for j in range(n):
                        if {i, j} == {0, 2} or (abs(i - j) == 1 and {i, j} != {0, 1}):
                            x = m.element(e).operate(i, j, i)
                            y = m.element(e).operate(j, i, j)
                            # print('  R1:', 'i =', i, 'j =', j, '->', x, '==', y)
                            assert x == y
                        elif {i, j} == {0, 1} or abs(i - j) > 1:
                            x = m.element(e).operate(i, j)
                            y = m.element(e).operate(j, i)
                            # print('  R2:', 'i =', i, 'j =', j, '->', x, '==', y)
                            assert x == y
                        elif i == j:
                            f = m.operate(e, i)
                            if e == f:
                                assert m.element(e).operate(i) in [m.element(e) * q(1), m.element(e) * -q(-1)]
                            else:
                                assert m.element(e).operate(i, i) == m.element(e) + m.element(e).operate(i) * (q(1) - q(-1))


def test_formula_gelfand_d(nin=5):
    def truncated(n, k, e):
        return gelfand_d_printer(n, k // 2)(m.reduced_word(e))[1]

    for n in range(3, nin + 1, 2):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            for e in m:
                v = truncated(n, k, e)
                a = (-1)**v.half_signs()
                for i in range(n):
                    f = m.operate(e, i)
                    w = truncated(n, k, f)
                    s = SignedPermutation.ds_i(i if i > 0 else -1, n)
                    if i > 0 and abs(v(i)) != i and abs(v(i + 1)) != i + 1:
                        if abs(v(i)) == i + 1 and abs(v(i + 1)) == i:
                            assert e == f
                            assert i in m.weak_descents(e)
                        elif v(i) < v(i + 1):
                            assert m.length(e) < m.length(f)
                        elif v(i) > v(i + 1):
                            assert m.length(e) > m.length(f)
                        assert w == s * v * s
                    if i > 0 and (abs(v(i)) != i or abs(v(i + 1)) != i + 1):
                        # print(m.height(e), ':', v, '*', i, '->', w, ':', m.height(f))
                        if i == -a * v(i) or i + 1 == a * v(i + 1):
                            assert m.length(e) < m.length(f)
                        if i == a * v(i) or i + 1 == -a * v(i + 1):
                            assert m.length(e) > m.length(f)
                        assert w == s * v * s
                    if i > 0 and abs(v(i)) == i and abs(v(i + 1)) == i + 1:
                        if i == -a * v(i) and i + 1 == a * v(i + 1):
                            assert m.length(e) < m.length(f)
                        elif i == a * v(i) and i + 1 == -a * v(i + 1):
                            assert m.length(e) > m.length(f)
                        else:
                            assert e == f
                            assert i in m.weak_ascents(e)
                        assert w == s * v * s

                    if i == 0 and abs(v(1)) != 1 and abs(v(2)) != 2:
                        # print(m.height(e), ':', v, '*', i, '->', w, ':', m.height(f))
                        if abs(v(1)) == 2 and abs(v(2)) == 1:
                            assert e == f
                            assert i in m.weak_descents(e)
                        elif (v(1) < 0 and v(2) < 0) or 0 < v(2) < -v(1) or 0 < v(1) < -v(2):
                            assert m.length(e) > m.length(f)
                        else:
                            assert m.length(e) < m.length(f)
                        assert w == s * v * s
                    elif i == 0 and (abs(v(1)) != 1 or abs(v(2)) != 2):
                        # print(m.height(e), ':', v, '*', i, '->', w, ':', m.height(f))
                        if 1 == a * v(1) or 2 == a * v(2):
                            assert m.length(e) < m.length(f)
                        if 1 == -a * v(1) or 2 == -a * v(2):
                            assert m.length(e) > m.length(f)
                        assert w == (s * v * s).dstar()
                    elif i == 0 and abs(v(1)) == 1 and abs(v(2)) == 2:
                        t = SignedPermutation.ds_i(1, n)
                        if 1 == a * v(1) and 2 == a * v(2):
                            assert m.length(e) < m.length(f)
                        elif 1 == -a * v(1) and 2 == -a * v(2):
                            assert m.length(e) > m.length(f)
                        else:
                            assert e == f
                            assert i in m.weak_ascents(e)
                        assert w == s * v * t
