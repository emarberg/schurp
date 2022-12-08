from stable.permutations import Permutation
from stable.tableaux import Partition


def test_init():
    w = Permutation(1)
    v = Permutation([1])
    u = Permutation(*[1])
    assert u == v == w


def test_constructor():
    try:
        Permutation(3)
        assert False
    except:
        pass

    # valid constructions
    w = Permutation(1, 2)
    assert all(w(i) == i for i in range(10))

    w = Permutation(-2, 1)
    assert w(1) == -2
    assert w(-1) == 2
    assert w(0) == 0
    assert w(2) == 1
    assert w(-2) == -1

    # input args must represent each congruence class exactly once
    try:
        Permutation(3, 2, 5, 15)
        assert False
    except:
        pass


def test_identity():
    w = Permutation.identity()
    assert w == Permutation(1, 2)
    assert w.is_identity()
    assert w.is_involution()


def test_s_i():
    assert Permutation.s_i(1) == Permutation(2, 1)
    assert Permutation.s_i(1) == Permutation(2, 1, 3)
    assert Permutation.s_i(2) == Permutation(1, 3, 2)

    w = Permutation.s_i(5)
    assert w(5) == 6
    assert w(6) == 5
    assert w(55) == 55
    assert w(56) == 56


def test_inverse():
    r = Permutation.s_i(1)
    s = Permutation.s_i(2)
    t = Permutation.s_i(3)

    assert r.inverse() == r
    assert (r * s).inverse() == s * r
    assert (r * s * t).inverse() == t * s * r
    assert (r * s * t * r).inverse() == r * t * s * r
    assert (r * s * t * r * s).inverse() == s * r * t * s * r
    assert (r * s * t * r * s * t).inverse() == t * s * r * t * s * r

    w = r * s * t * r * s * t
    assert w.inverse() == w**-1


def test_multiplication():
    s = Permutation.s_i(1)
    t = Permutation.s_i(2)
    u = Permutation.s_i(3)

    assert s * s == Permutation.identity()
    assert s % s == s
    assert s * t * u == s % t % u
    assert s * t * s * t == t * s
    assert s % t % s % t == s * t * s == t * s * t


def test_lt():
    assert Permutation(1, 2, 3) < Permutation(1, 3, 2)
    assert not (Permutation(1, 2, 3) < Permutation(1, 2, 3))
    assert not (Permutation(1, 3, 2) < Permutation(1, 2, 3))


def test_length():
    w = Permutation.identity()
    assert len(w) == w.length == 0
    assert w.involution_length == 0

    r = Permutation.s_i(1)
    s = Permutation.s_i(2)
    t = Permutation.s_i(3)

    assert len(r) == 1
    assert len(r * s) == 2
    assert len(r * s * t) == 3
    assert len(r * s * t * r) == 4
    assert len(r * s * t * r * s) == 5
    assert len(r * s * t * r * s * t) == 4

    assert r.involution_length == 1
    assert (s * r * s).involution_length == 2
    assert (t * s * r * s * t).involution_length == 3
    assert (r * t * s * r * s * t * r).involution_length == 2
    assert (s * r * t * s * r * s * t * r * s).involution_length == 1
    assert (t * s * r * t * s * r * s * t * r * s * t).involution_length == 1


def test_descents():
    w = Permutation(3, 5, -6, 1, -2, 4)
    assert w.right_descent_set == {2, 4}
    assert w.left_descent_set == {1, 4, 5}

    w = Permutation.identity()
    assert w.right_descent_set == set()
    assert w.left_descent_set == set()


def test_all():
    ans = 1
    for n in range(2, 6):
        ans *= n
        assert len({w for w in Permutation.all(n, signed=False)}) == ans


def test_reduced_words():
    w = Permutation.identity()
    assert set(w.get_reduced_words()) == {tuple()}
    assert set(w.reduced_words) == {tuple()}

    r = Permutation.s_i(1)
    s = Permutation.s_i(2)
    t = Permutation.s_i(3)

    assert set(r.reduced_words) == {(1,)}
    assert set((r * s).reduced_words) == {(1, 2)}
    assert set((r * s * t).reduced_words) == {(1, 2, 3)}
    assert set((r * s * t * r).reduced_words) == {(1, 2, 3, 1), (1, 2, 1, 3), (2, 1, 2, 3)}
    assert set((r * s * r).reduced_words) == {(1, 2, 1), (2, 1, 2)}

    n = 5
    w = Permutation(*reversed(range(1, n + 1)))
    words = set(w.get_reduced_words())
    assert len(words) == 768
    for e in words:
        v = Permutation.identity()
        for i in e:
            v *= Permutation.s_i(i)
        assert v == w


def test_involution_words():
    w = Permutation.identity()
    assert set(w.get_involution_words()) == {tuple()}
    assert set(w.involution_words) == {tuple()}

    r = Permutation.s_i(1)
    s = Permutation.s_i(2)
    t = Permutation.s_i(3)

    # can only access this property for affine permutations which are involutions
    try:
        set((r * s).involution_words)
        assert False
    except:
        pass

    assert set(r.involution_words) == {(1,)}
    assert set((s * r * s).involution_words) == {(1, 2), (2, 1)}
    assert set((t * s * r * s * t).involution_words) == {(1, 2, 3), (2, 1, 3), (2, 3, 1), (3, 2, 1)}

    n = 5
    w = Permutation(*reversed(range(1, n + 1)))
    words = set(w.get_involution_words())
    assert len(words) == 80
    for e in words:
        v = Permutation.identity()
        for i in e:
            s = Permutation.s_i(i)
            v = s % v % s
        assert v == w


def test_involutions():
    assert set(Permutation.involutions(1)) == {Permutation()}
    assert set(Permutation.involutions(1, True)) == {Permutation(), Permutation(-1)}
    assert set(Permutation.involutions(2)) == {Permutation(), Permutation(2, 1)}
    assert set(Permutation.involutions(2, True)) == {
        Permutation(),
        Permutation(-1),
        Permutation(2, 1),
        Permutation(-2, -1),
        Permutation(-1, -2),
        Permutation(1, -2)
    }

    n = 6
    assert set(Permutation.involutions(n)) == {
        w for w in Permutation.all(n) if w.inverse() == w
    }
    assert set(Permutation.involutions(n, True)) == {
        w for w in Permutation.all(n, True) if w.inverse() == w
    }
    assert len(set(Permutation.involutions(n))) == len(list(Permutation.involutions(n)))

    assert set(Permutation.fpf_involutions(4)) == {
        Permutation(2, 1, 4, 3),
        Permutation(3, 4, 1, 2),
        Permutation(4, 3, 2, 1),
    }


def test_atoms():
    y = Permutation(3, 4, 1, 2, 7, 9, 5, 10, 6, 8, 12, 11)
    assert y.get_min_atom() in y.get_atoms()

    y = Permutation(3, 4, 1, 2)
    assert Permutation(3, 1, 4, 2).inverse() in y.get_atoms()
    assert Permutation(3, 1, 4, 2) not in y.get_atoms()


def test_fpf_atoms():
    y = Permutation(4, 3, 2, 1)
    s = Permutation.s_i(1)
    t = Permutation.s_i(2)
    u = Permutation.s_i(3)
    assert set(y.get_fpf_atoms()) == {t * s, t * u}

    n = 6
    for w in Permutation.fpf_involutions(n):
        atoms = set()
        for word in w.get_fpf_involution_words():
            atoms.add(Permutation.from_word(word))
        assert atoms == set(w.get_fpf_atoms())


def test_mod():
    s = Permutation.s_i(1)
    t = Permutation.s_i(2)

    assert (s * t) % (t * s) == s * t * s
    assert (s * t) % (t * s * t) == s * t * s


def test_transitions_finite():
    n = 4
    for w in Permutation.all(n):
        for r in range(1, n + 1):
            phi_plus = {w * Permutation.transposition(r, j) for j in w.upper_transitions(r)}
            phi_minus = {w * Permutation.transposition(i, r) for i in w.lower_transitions(r)}

            assert all(w in v.bruhat_covers() for v in phi_plus)
            assert all(w in v.bruhat_covers() for v in phi_minus)


def test_transitions_all():
    n = 4
    for w in Permutation.all(n):
        for r in range(1, n + 1):
            phi_plus = {w * Permutation.transposition(r, j, n) for j in w.upper_transitions(r)}
            phi_minus = {w * Permutation.transposition(i, r, n) for i in w.lower_transitions(r)}

            assert all(w in v.bruhat_covers() for v in phi_plus)
            assert all(w in v.bruhat_covers() for v in phi_minus)

            print (w.get_reduced_word())
            print(r)
            print([x.get_reduced_word() for x in phi_plus])
            print([x.get_reduced_word() for x in phi_minus])


def test_fpf_transitions_simple():
    y = Permutation(2, 1, 4, 3)
    p = 1
    q = 2

    def t(i, j):
        return Permutation.transposition(i, j)

    phi_plus = {t(q, j) * y * t(q, j) for j in y.upper_fpf_involution_transitions(q)}
    assert phi_plus == {Permutation(3, 4, 1, 2)}
