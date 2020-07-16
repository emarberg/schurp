from qp import QPWGraph, QPModule, gelfand_d_printer
import polynomials
from signed import SignedPermutation
from even import EvenSignedPermutation
from tableaux import Tableau
from qp_utils import rsk, gelfand_rsk, truncate_a, truncate_bc
import random


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
            cells = []
            molecules = []
            edges = []
            for k in range(0, n + 2, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_a, QPModule.create_gelfand_a)
                read_or_compute_wgraph(w)
                cells += w.get_cells_as_permutations()
                molecules += w.get_molecules_as_permutations()
                edges += [w.get_wgraph_size()]
            print('sgn =', sgn, 'n =', n, '::', '#edges =', edges, '#cells =', len(cells), '#molecules =', len(molecules))
            print()
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


def test_gelfand_cells_bc(nn=5, s=None):
    cellmap = {}
    for n in range(2, nn + 1):
        cellmap[n] = {}
        for sgn in ([False, True] if s is None else [s]):
            cellmap[n][sgn] = set()
            cells = []
            molecules = []
            edges = []
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_bc, QPModule.create_gelfand_bc)
                read_or_compute_wgraph(w)
                print(n, k, len(w.cells))
                cells += w.cells
                molecules += w.molecules
                edges += [w.get_wgraph_size()]
                cellmap[n][sgn] |= w.get_cells_as_permutations()
            print()
            print('sgn =', sgn, 'n =', n, '::', '#edges =', edges, '#cells =', len(cells), '#molecules =', len(molecules))
            print()
            if (not sgn) in cellmap[n]:
                def toggle(c):
                    x = SignedPermutation(*[-a for a in c])
                    s = SignedPermutation(*(list(range(1, n + 1)) + list(range(x.rank, n, -1))))
                    return tuple((s * x * s).oneline)

                negated = {tuple(sorted(toggle(c) for c in cell)) for cell in cellmap[n][not sgn]}
                if negated != cellmap[n][sgn]:
                    for c in sorted(cellmap[n][sgn], key=len):
                        print('1 ', c)
                    print()
                    for c in sorted(negated, key=len):
                        print('2 ', c)
                assert negated == cellmap[n][sgn]
                print('* isomorphism checked\n')


def test_gelfand_cells_d(nn=5, s=None):
    cellmap = {}
    for n in range(2, nn + 1):
        cellmap[n] = {}
        for sgn in ([False, True] if s is None else [s]):
            cellmap[n][sgn] = set()
            cells = []
            molecules = []
            edges = []
            for k in range(0, n + 1, 2):
                w = read_or_create(n, k // 2, sgn, QPModule.read_gelfand_d, QPModule.create_gelfand_d)
                read_or_compute_wgraph(w)
                print(n, k, len(w.cells))
                cells += w.cells
                molecules += w.molecules
                edges += [w.get_wgraph_size()]
                cellmap[n][sgn] |= w.get_cells_as_permutations()
            print()
            print('sgn =', sgn, 'n =', n, '::', '#edges =', edges, '#cells =', len(cells), '#molecules =', len(molecules))
            print()
            if (not sgn) in cellmap[n]:
                def toggle(c):
                    x = EvenSignedPermutation(*[-a for a in c])
                    s = EvenSignedPermutation(*(list(range(1, n + 1)) + list(range(x.rank, n, -1))))
                    x = s * x * s
                    if x.rank % 4 != 0:
                        x = x.star()
                    return tuple(x.oneline)

                negated = {tuple(sorted(toggle(c) for c in cell)) for cell in cellmap[n][not sgn]}
                if negated != cellmap[n][sgn]:
                    for c in sorted(cellmap[n][sgn], key=len):
                        print('1 ', c)
                    print()
                    for c in sorted(negated, key=len):
                        print('2 ', c)
                    assert negated == cellmap[n][sgn]
                print('* isomorphism checked\n')



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
