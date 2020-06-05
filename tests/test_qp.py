from qp import QPWGraph, QPModule, gelfand_d_printer
import polynomials
from signed import SignedPermutation
from utils import rsk
import random


def test_hecke_cells_a(n=4):
    m = QPModule.create_hecke_a(n)
    w = QPWGraph(m)
    w.compute_wgraph()
    cells = w.cells
    w.print_wgraph()
    for c in cells:
        r = {rsk(w)[0] for w in c}
        print(c)
        print(r)
        print()
        assert len(r) == 1


def test_gelfand_cells_a(n=4):
    cells = []
    for k in range(0, n + 2, 2):
        m = QPModule.create_gelfand_a(n, k // 2)
        w = QPWGraph(m, sgn=False)
        w.compute_wgraph()
        cells += w.cells

    seen = set()
    for c in cells:
        r = {rsk(w)[1].restrict(n + 1).partition() for w in c}
        print(c)
        print({rsk(w)[1].restrict(n + 1) for w in c})
        print()
        assert len(r) == 1
        assert not r.issubset(seen)
        seen |= r



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
