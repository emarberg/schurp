from modules import QPModule


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


def test_gelfand_a(nin=7):
    for n in range(nin + 1):
        for k in range(0, n + 1, 2):
            m = QPModule.create_gelfand_a(n, k // 2)
            for e in m:
                for i in range(n - 1):
                    for j in range(n - 1):
                        if abs(i - j) > 1:
                            assert m.act(e, i, j) == m.act(e, j, i)
                        elif abs(i - j) == 1:
                            assert m.act(e, i, j, i) == m.act(e, j, i, j)
                        else:
                            assert m.act(e, i, i) == e
