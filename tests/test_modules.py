from modules import QPModule, q


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
