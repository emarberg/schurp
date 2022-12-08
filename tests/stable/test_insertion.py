from stable.tableaux import Tableau
from stable.insertion import FState, InsertionAlgorithm
from stable.permutations import Permutation


def test_next_f_r1():
    # transition (R1)
    state = FState().add(1)

    t = Tableau().add(1, 2, 1)
    assert state == FState(t, (1, 2))

    state, box = state.next()
    i, j = box

    t = Tableau().add(1, 1, 1)
    assert state == FState(t)
    assert i == 1
    assert j == 1


def test_next_f_r2():
    # transition (R2)
    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 5,
        (2, 2): 4, (2, 3): 6, (2, 4): 7, (2, 6): 5
    })
    state = FState(t, (2, 6))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 5,
        (2, 2): 4, (2, 3): 5, (2, 4): 7,
        (3, 6): 6
    })
    assert state == FState(t, (3, 6))
    assert i == 2
    assert j == 3


def test_next_f_r3():
    # transition (R3)
    #
    #   1 2 3  1
    # 5 5 3 4
    #
    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 1, (2, 3): 2, (2, 4): 3, (2, 6): 1
    })
    state = FState(t, (2, 6))

    state, box = state.next()
    i, j = box

    #     2
    #
    #   1 2 3
    # 5 5 5 4
    #
    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 1, (2, 3): 2, (2, 4): 3, (4, 3): 2
    })
    assert state == FState(t, (4, 3))
    assert i == 2
    assert j == 3


def test_next_f_d1():
    # transition (D1)
    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 5,
        (2, 2): 6, (2, 3): 7, (2, 4): 8, (2, 6): 4
    })
    state = FState(t, (2, 6))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 5,
        (2, 2): 4, (2, 3): 7, (2, 4): 8, (4, 3): 6
    })
    assert state == FState(t, (4, 3))
    assert i == 2
    assert j == 2


def test_next_f_d2():
    # transition (D2)
    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 4, (2, 3): 5, (2, 4): 6, (2, 6): 3
    })
    state = FState(t, (2, 6))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 4, (2, 3): 5, (2, 4): 6, (4, 3): 5
    })
    assert state == FState(t, (4, 3))
    assert i == 2
    assert j == 2


def test_next_f_c1():
    # transition (C1)
    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 4, (4, 3): 6
    })
    state = FState(t, (4, 3))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 5, (1, 2): 5, (1, 3): 5, (1, 4): 4,
        (2, 2): 4, (2, 3): 6
    })
    assert state == FState(t)
    assert i == 2
    assert j == 3


def test_next_f_c2():
    # transition (C2)
    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 5, (1, 4): 6,
        (2, 2): 4, (2, 3): 6, (2, 4): 7,
        (3, 3): 8,
        (5, 3): 4
    })
    state = FState(t, (5, 3))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 6,
        (2, 2): 4, (2, 3): 6, (2, 4): 7,
        (3, 3): 8,
        (5, 4): 5
    })
    assert state == FState(t, (5, 4))
    assert i == 1
    assert j == 3


def test_next_f_c3():
    # transition (C3)
    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 5, (1, 4): 6,
        (2, 2): 4, (2, 3): 6, (2, 4): 7,
        (3, 3): 8,
        (5, 3): 5
    })
    state = FState(t, (5, 3))

    state, box = state.next()
    i, j = box

    t = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 5, (1, 4): 6,
        (2, 2): 4, (2, 3): 6, (2, 4): 7,
        (3, 3): 8,
        (5, 4): 6
    })
    assert state == FState(t, (5, 4))
    assert i == 2
    assert j == 3


def test_fstate_insertion():
    a = Tableau()
    b = Tableau({(1, 1): 4})
    c = Tableau({(1, 1): 2, (1, 2): 4})
    d = Tableau({(1, 1): 2, (1, 2): 3, (2, 2): 4})
    e = Tableau({(1, 1): 2, (1, 2): 3, (2, 2): 4, (1, 3): 4})
    f = Tableau({(1, 1): 2, (1, 2): 3, (2, 2): 4, (1, 3): 4, (2, 3): 5})

    state = FState()
    assert state.is_terminal() and state.tableau == a

    state, path = state.insert(4)
    assert state.is_terminal() and state.tableau == b
    assert path == [(1, 1)]

    state, path = state.insert(2)
    assert state.is_terminal() and state.tableau == c
    assert path == [(1, 1), (1, 2)]

    state, path = state.insert(3)
    assert state.is_terminal() and state.tableau == d
    assert path == [(1, 2), (2, 2)]

    state, path = state.insert(1)
    assert state.is_terminal() and state.tableau == e
    assert path == [(1, 1), (2, 2), (1, 3)]

    state, path = state.insert(2)
    assert state.is_terminal() and state.tableau == f
    assert path == [(1, 2), (2, 2), (2, 3)]

    p, q = FState.insertion_tableaux(4, 2, 3, 1, 2)
    assert p == f


def test_reverse_insertion():
    w = (2, 1, 4, 1, 4, 3, 4, 1, 2)
    p, q = FState.insertion_tableaux(*w)
    r = (2, 3, {-9})

    assert q.values() == {1, -2, 3, -4, 5, 6, 7, -8, -9}

    assert p == Tableau({(1, 1): 2, (1, 2): 3, (1, 3): 4, (2, 2): 4, (2, 3): 5})

    s, box = FState(p).previous(r)
    assert box == (2, 3)
    assert s.tableau == Tableau({(1, 1): 2, (1, 2): 3, (1, 3): 4, (2, 2): 4, (4, 3): 5})
    assert s.outer == (4, 3)

    s, box = s.previous()
    assert box == (2, 2)
    assert s.tableau == Tableau({(1, 1): 2, (1, 2): 3, (1, 3): 4, (2, 2): 4, (2, 5): 3})
    assert s.outer == (2, 5)

    s, box = s.previous()
    assert box == (1, 2)
    assert s.tableau == Tableau({(1, 1): 2, (1, 2): 3, (1, 3): 4, (2, 2): 4, (1, 5): 2})
    assert s.outer == (1, 5)
    assert s.is_initial()

    assert FState.inverse_insertion(p, q) == w

    #
    #     22
    #   6  8 22
    # 2 4  6 20
    #
    w = (2, 6, 8, 22, 4, 6, 20, 6)
    p, q = FState.insertion_tableaux(*w)
    assert p == Tableau({
        (3, 3): 22,
        (2, 2): 6, (2, 3): 8, (2, 4): 22,
        (1, 1): 2, (1, 2): 4, (1, 3): 6, (1, 4): 20
    })
    assert FState.inverse_insertion(p, q) == w


def test_reverse_insertion_complete():
    n = 8
    count = 200
    for w in Permutation.symplectic_hecke_words(n):
        count -= 1
        if count == 0:
            break
        print(w)
        p, q = FState.insertion_tableaux(*w)
        print()
        print(p)
        print()
        print(q)
        v = FState.inverse_insertion(p, q)
        print()
        print(v)
        print()
        assert v == w


def test_symplectic_hecke_insertion():
    w = (4, 2, 6, 1, 7, 5, 3, 4, 2, 1, 3, 2)
    i = (1, 2, 2, 3, 3, 4, 5, 5, 6, 8, 8, 9)

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w)

    p = Tableau({
        (1, 1): 2, (1, 2): 3, (1, 3): 4, (1, 4): 5, (1, 5): 6, (1, 6): 7,
        (2, 2): 4, (2, 3): 5, (2, 4): 6, (2, 5): 7,
        (3, 3): 6, (3, 4): 7,
    })
    q = Tableau({
        (1, 1): 1, (1, 2): -2, (1, 3): 3, (1, 4): -4, (1, 5): 5, (1, 6): -10,
        (2, 2): 6, (2, 3): -7, (2, 4): -9, (2, 5): -12,
        (3, 3): 8, (3, 4): -11,
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, tuple(range(1, 13)))

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w, i)

    q = Tableau({
        (1, 1): 1, (1, 2): -2, (1, 3): 2, (1, 4): -3, (1, 5): 3, (1, 6): -8,
        (2, 2): 4, (2, 3): -5, (2, 4): -6, (2, 5): -9,
        (3, 3): 5, (3, 4): -8,
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, i)


def test_symplectic_hecke_insertion_setvalued():
    w = (2, 2, 4, 3)
    i = (1, 1, 1, 4)

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w)

    p = Tableau({
        (1, 1): 2, (1, 2): 3,
        (2, 2): 4
    })
    q = Tableau({
        (1, 1): (1, 2), (1, 2): 3,
        (2, 2): 4
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, (1, 2, 3, 4))

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w, i)

    q = Tableau({
        (1, 1): (1, 1), (1, 2): 1,
        (2, 2): 4
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, i)

    i = (1, 2, 3, 3)
    insertion, recording = InsertionAlgorithm.symplectic_hecke(w, i, False)

    q = Tableau({
        (1, 1): (1, 2), (1, 2): -3,
        (2, 2): 3
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q, False) == (w, i)

    w = (4, 2, 2, 3)
    i = (2, 4, 4, 4)

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w)

    p = Tableau({
        (1, 1): 2, (1, 2): 3,
        (2, 2): 4
    })
    q = Tableau({
        (1, 1): 1, (1, 2): (-3, -2),
        (2, 2): 4
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, (1, 2, 3, 4))

    insertion, recording = InsertionAlgorithm.symplectic_hecke(w, i)

    q = Tableau({
        (1, 1): 2, (1, 2): (-4, -4),
        (2, 2): 4
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q) == (w, i)

    i = (1, 1, 2, 3)
    insertion, recording = InsertionAlgorithm.symplectic_hecke(w, i, False)

    q = Tableau({
        (1, 1): 1, (1, 2): (1, 2),
        (2, 2): 3
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_symplectic_hecke(p, q, False) == (w, i)


def test_orthogonal_hecke_insertion():
    w = (4, 5, 1, 1, 3, 2)
    i = (1, 1, 3, 4, 4, 6)

    insertion, recording = InsertionAlgorithm.orthogonal_hecke(w)

    p = Tableau({
        (1, 1): 1, (1, 2): 2, (1, 3): 4, (1, 4): 5,
        (2, 2): 3
    })
    q = Tableau({
        (1, 1): 1, (1, 2): 2, (1, 3): (-3, -4), (1, 4): -6,
        (2, 2): 5
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_orthogonal_hecke(p, q) == (w, (1, 2, 3, 4, 5, 6))

    insertion, recording = InsertionAlgorithm.orthogonal_hecke(w, i)

    q = Tableau({
        (1, 1): 1, (1, 2): 1, (1, 3): (-3, -4), (1, 4): -6,
        (2, 2): 4
    })
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_orthogonal_hecke(p, q) == (w, i)


def test_hecke_insertion():
    w = (3, 2, 1, 1)
    i = (1, 1, 1, 2)
    insertion, recording = InsertionAlgorithm.hecke(w, i, multisetvalued=False)
    p = Tableau({
        (1, 1): 1, (1, 2): 2, (1, 3): 3
    })
    q = Tableau({
        (1, 1): 1, (1, 2): 1, (1, 3): (1, 2)
    })
    print(insertion)
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_hecke(p, q, multisetvalued=False) == (w, i)

    p = Tableau({
        (1, 1): 1, (1, 2): 2,
        (2, 1): 2
    })
    q = Tableau({
        (1, 1): 2, (1, 2): 3,
        (2, 1): 3
    })
    assert InsertionAlgorithm.inverse_hecke(p, q, multisetvalued=False) == ((1, 2, 1), (2, 3, 3))

    i = (1, 2, 3, 3)
    insertion, recording = InsertionAlgorithm.hecke(w, i)
    p = Tableau({
        (1, 1): 1, (2, 1): 2, (3, 1): 3
    })
    q = Tableau({
        (1, 1): 1, (2, 1): 2, (3, 1): (3, 3)
    })
    print(insertion)
    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_hecke(p, q) == (w, i)

    w = (1, 5, 1, 3, 3)
    i = (1, 1, 3, 3, 3)
    insertion, recording = InsertionAlgorithm.hecke(w)

    p = Tableau({
        (1, 1): 1, (1, 2): 3,
        (2, 1): 5
    })
    q = Tableau({
        (1, 1): 1, (1, 2): (2, 5),
        (2, 1): (3, 4)
    })

    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_hecke(p, q) == (w, (1, 2, 3, 4, 5))

    insertion, recording = InsertionAlgorithm.hecke(w, i)

    q = Tableau({
        (1, 1): 1, (1, 2): (1, 3),
        (2, 1): (3, 3)
    })

    assert insertion == p
    assert recording == q
    assert InsertionAlgorithm.inverse_hecke(p, q) == (w, i)

    insertion, q = InsertionAlgorithm.hecke((3, 4, 1, 2, 4))

    p = Tableau({
        (1, 1): 1, (1, 2): 2, (1, 3): 4,
        (2, 1): 3, (2, 2): 4,
    })
    assert insertion == p
    assert InsertionAlgorithm.inverse_hecke(p, q) == ((3, 4, 1, 2, 4), (1, 2, 3, 4, 5))

    insertion, q = InsertionAlgorithm.hecke((3, 1, 2, 4))

    p = Tableau({
        (1, 1): 1, (1, 2): 2, (1, 3): 4,
        (2, 1): 3
    })
    assert insertion == p
    assert InsertionAlgorithm.inverse_hecke(p, q) == ((3, 1, 2, 4), (1, 2, 3, 4))
