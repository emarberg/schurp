from signed import SignedPermutation
from polynomials import X


def test_get_minimal_fpf_involution():
    assert SignedPermutation.get_minimal_fpf_involution(1) == SignedPermutation(-1)
    assert SignedPermutation.get_minimal_fpf_involution(2) == SignedPermutation(2, 1)
    assert SignedPermutation.get_minimal_fpf_involution(3) == SignedPermutation(-1, 3, 2)
    assert SignedPermutation.get_minimal_fpf_involution(4) == SignedPermutation(2, 1, 4, 3)


def test_ncsp_matchings():
    m = ((-5, -1), (-4, -3), (-2, 2), (1, 5), (3, 4))
    base = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    ncsp = list(SignedPermutation.ncsp_matchings(base))
    assert m not in ncsp

    base = {-8, -6, -4, -2, 2, 4, 6, 8}
    assert set(SignedPermutation.ncsp_matchings(base)) == {
        ((-8, 8), (-6, 6), (-4, 4), (-2, 2)),
        ((-8, -6), (-4, 4), (-2, 2), (6, 8)),
        ((-8, -2), (-6, -4), (2, 8), (4, 6)),
        ((-8, 8), (-6, 6), (-4, -2), (2, 4)),
        ((-8, 8), (-6, -4,), (-2, 2), (4, 6)),
        ((-8, -6), (-4, -2), (2, 4), (6, 8))
    }


def test_get_involution_word():
    s = SignedPermutation.s_i(0, 3)
    t = SignedPermutation.s_i(1, 3)
    assert (t * s * t).get_involution_word() == (0, 1)
    assert (s * t * s).get_involution_word() == (1, 0)
    assert set((t * s * t).get_involution_words()) == {(0, 1)}
    assert set((s * t * s).get_involution_words()) == {(1, 0)}


def test_fpf_involution_words(n=4):
    for w in SignedPermutation.fpf_involutions(n):
        words = set(w.get_fpf_involution_words())
        if words:
            print(w, '->', len(words))


def test_get_atoms(n=4):
    i = list(SignedPermutation.involutions(n))
    for w in i:
        words = set()
        for a in w.get_atoms():
            assert len(a) == w.involution_length()
            assert a.inverse() % a == w
            words |= set(a.get_reduced_words())
        assert words == set(w.get_involution_words())


def test_get_abs_fpf_atoms(n=6):
    for y in SignedPermutation.abs_fpf_involutions(n):
        s = set()
        for w in y.get_atoms():
            winv = w.inverse()
            if all(winv(i - 1) > winv(i) for i in range(2, y.rank, 2)):
                s.add(SignedPermutation.get_minimal_fpf_involution(n) * w)
        assert s == set(y.get_fpf_atoms())


def test_brion_length(n=4):
    for y in SignedPermutation.involutions(n):
        ell = len(y.neg()) + len(y.pair())
        assert ell + len(y) == 2 * y.involution_length()
        for w in y.get_atoms():
            ndes = w.ndes()
            cdes = [(a, b) for a, b in ndes if 0 < a < -b]
            nfix = w.nfix()
            nneg = w.nneg()
            print('b:', w.brion_length_b(), 'c:', w.brion_length_c(), ':', ell, len(y.neg()), len(y.pair()), ':', 'cdes =', len(cdes), 'ndes =', len(ndes), 'nfix =', len(nfix), 'nneg =', len(nneg), w == y.get_min_atom(), w.shape())
            assert w.brion_length_b() == ell - len(cdes)
            assert w.brion_length_c() == len(y.pair()) + len(cdes)
            assert len(cdes) * 2 == len(y.neg()) - len(nneg)
        print()


def test_brion_weight_counts(n=4):
    x = X(0)
    for m in range(1, n + 1):
        y = SignedPermutation.longest_element(m)
        a = 0
        b = 0
        c = 0
        for w in y.get_atoms():
            a += len(w.get_reduced_words())
            b += x**w.brion_length_b() * len(w.get_reduced_words())
            c += x**w.brion_length_c() * len(w.get_reduced_words())
        print(m, ':', a, '\t\t\t', b, '\t\t\t', c)


def test_shape(n=4):
    cls = SignedPermutation
    w0 = cls.longest_element(n)
    for w in cls.involutions(n):
        shapes = {}
        for a in w.get_atoms():
            sh = tuple(sorted(a.shape()))
            shapes[sh] = shapes.get(sh, []) + [a]
        print(w, ' :: ', w * w0, ' :: ', (w * w0).fixed_points())
        print()
        for sh, atoms in shapes.items():
            print(' ', set(sh), '->', atoms)
            print(' ', len(str(set(sh))) * ' ', '  ', [a.inverse() for a in atoms])
            print()
        print()
        print()
