from modules import QPModule, q
from signed import SignedPermutation


def test_eq():
    m = QPModule.create_gelfand_a(2, 0)
    n = QPModule.create_gelfand_a(2, 0)
    assert m == n

    n = QPModule.create_gelfand_a(2, 1)
    assert m != n
    assert len(m) == len(n)


def test_create():
    m = QPModule.create_gelfand_a(0, 0)
    m = QPModule.create_gelfand_a(1, 0)

    m = QPModule.create_gelfand_a(2, 0)
    assert set(m.weak_descents(0)) == set()
    assert set(m.weak_ascents(0)) == {0}

    m = QPModule.create_gelfand_a(2, 1, False)
    assert set(m.weak_descents(0)) == set()
    assert set(m.weak_ascents(0)) == {0}

    m = QPModule.create_gelfand_a(2, 1)
    assert set(m.weak_descents(0)) == {0}
    assert set(m.weak_ascents(0)) == set()

    m = QPModule.create_gelfand_a(2, 0, False)
    assert set(m.weak_descents(0)) == {0}
    assert set(m.weak_ascents(0)) == set()

    m = QPModule.create_gelfand_a(2, 1)
    m = QPModule.create_gelfand_a(3, 1)


def test_gelfand_a(nin=5):
    for n in range(nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_a(n, k // 2)
            assert len({m.printer(m.reduced_word(i)) for i in m}) == m.size
            for e in m:
                for i in range(n - 1):
                    for j in range(n - 1):
                        print(n, k, e, i, j)
                        if abs(i - j) > 1:
                            assert m.operate(e, i, j) == m.operate(e, j, i)
                        elif abs(i - j) == 1:
                            assert m.operate(e, i, j, i) == m.operate(e, j, i, j)
                        else:
                            assert m.operate(e, i, i) == e


def test_hecke_gelfand_a(nin=4):
    for n in range(nin + 1):
        for k in range(0, n + 1, 2):
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
    for n in range(nin + 1):
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
    for n in range(nin + 1):
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


def test_gelfand_d(nin=5):
    for n in range(3, nin + 1, 2):
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


def test_hecke_gelfand_d(nin=5):
    for n in range(3, nin + 1, 2):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            for e in m:
                print('n =', n, 'k =', k // 2, ':', 'element =', e, 'of', m.size)
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
    for n in range(3, nin + 1, 2):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_d(n, k // 2)
            for e in m:
                v = m.permutation(e)
                a = (-1)**v.dkappa()
                for i in range(n):
                    f = m.operate(e, i)
                    w = m.permutation(f)
                    s = SignedPermutation.ds_i(i if i > 0 else -1, n)
                    if i > 0 and abs(v(i)) != i and abs(v(i + 1)) != i + 1:
                        if abs(v(i)) == i + 1 and abs(v(i + 1)) == i:
                            assert e == f
                            assert i in m.weak_ascents(e)
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
                            assert i in m.weak_descents(e)
                        assert w == s * v * s

                    if i == 0 and abs(v(1)) != 1 and abs(v(2)) != 2:
                        # print(m.height(e), ':', v, '*', i, '->', w, ':', m.height(f))
                        if abs(v(1)) == 2 and abs(v(2)) == 1:
                            assert e == f
                            assert i in m.weak_ascents(e)
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
                            assert i in m.weak_descents(e)
                        assert w == s * v * t
