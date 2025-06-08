from signed import SignedPermutation, Permutation
from even import EvenSignedPermutation
from polynomials import X


def test_dshape(n=4):
    for z in SignedPermutation.involutions(n, dtype=True):
        for w in z.get_atoms_d():
            print(z, w.inverse())
            print(w.dshape())
            print()


def test_b_bruhat_order(n=4):
    length = lambda x: x.length()
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n + 1) for j in range(i, n + 1)]
    
    bruhat = {v: {v} for v in SignedPermutation.all(n)}
    for v in sorted(bruhat, key=lambda x: length(x) * -1):
        for t in reflections:
            vt = v * t
            if length(vt) > length(v):
                bruhat[v] |= bruhat[vt]

    for v in SignedPermutation.all(n):
        for u in SignedPermutation.all(n):
            expected = v.strong_bruhat_less_equal(u)
            actual = u in bruhat[v]
            assert expected == actual


def test_d_bruhat_order(n=4):
    length = lambda x: x.dlength()
    reflections = [SignedPermutation.reflection_t(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    reflections += [SignedPermutation.reflection_s(i, j, n) for i in range(1, n) for j in range(i + 1, n + 1)]
    
    bruhat = {v: {v} for v in SignedPermutation.all(n, dtype=True)}
    for v in sorted(bruhat, key=lambda x: x.dlength() * -1):
        for t in reflections:
            vt = v * t
            if length(vt) > length(v):
                bruhat[v] |= bruhat[vt]

    for v in SignedPermutation.all(n, dtype=True):
        for u in SignedPermutation.all(n, dtype=True):
            expected = v.dbruhat_less_equal(u)
            actual = u in bruhat[v]
            print(v, '<=', u, ':', expected, actual)
            assert expected == actual


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
        pos =  y.count_nontrivial_positive_cycles()
        assert ell + len(y) == 2 * y.involution_length()
        for w in y.get_atoms():
            ndes = w.ndes()
            cdes = [(a, b) for a, b in ndes if 0 < a < -b]
            nfix = w.nfix()
            nneg = w.nneg()
            print('b:', w.brion_length_b(), 'c:', w.brion_length_c(), ':', ell, len(y.neg()), len(y.pair()), ':', 'cdes =', len(cdes), 'ndes =', len(ndes), 'nfix =', len(nfix), 'nneg =', len(nneg), w == y.get_min_atom(), w.shape())
            
            # (nontrivial blocks of w.shape()) / 2 == len(cdes) 
            # y.get_min_atom().ell_zero() == len(y.neg()) + |{a : a < 0 < -a < y(a)}|
            #                             == ell - |{a : 0 < a < y(a)}|

            # (nontrivial blocks of w.shape()) / 2 == ell - pos - w.ell_zero()

            sh = w.shape()
            nb = len([(a, b) for (a, b) in sh if a + b != 0])
            assert nb % 2 == 0
            assert nb // 2 == y.get_min_atom().ell_zero() - w.ell_zero()
            assert y.get_min_atom().ell_zero() == ell - pos

            assert w.brion_length_b() == ell - len(cdes)
            assert w.brion_length_b() == w.ell_zero() + pos

            assert w.brion_length_c() == len(y.pair()) + len(cdes)
            assert w.brion_length_c() == y.count_positive_ascents() - w.ell_zero()
            assert len(cdes) * 2 == len(y.neg()) - len(nneg)
        print()


def test_brion_length_d(n=4):
    for twisted in [False, True]:
        offset = 1 if twisted else 0
        u = SignedPermutation.s_i(0, n) if twisted else SignedPermutation.identity(n)
        for y in SignedPermutation.involutions(n, dtype=True, twisted=twisted):
            ell = len(y.neg()) + len(y.pair())
            # assert ell + len(y) == 2 * y.involution_length()
            for w in y.get_atoms_d(twisted=twisted):
                sh = w.dshape(offset)
                
                ndes = w.ndes()
                cdes = [(a, b) for a, b in ndes if 0 < a < -b]
                nfix = w.nfix()
                nneg = w.nneg()
                
                nb = len([(a, b) for (a, b) in sh if a + b != 0])
                pair = len((u*y).pair())
                
                assert w.brion_length_d(twisted=twisted) == nb // 2 + pair
                
                # (nontrivial blocks of w.shape()) / 2 == len(cdes) 
                # y.get_min_atom().ell_zero() == len(y.neg()) + |{a : a < 0 < -a < y(a)}|
                #                             == ell - |{a : 0 < a < y(a)}|

                # (nontrivial blocks of w.shape()) / 2 == ell - pos - w.ell_zero()

                #assert nb % 2 == 0
                #assert nb // 2 == y.get_min_atom().ell_zero() - w.ell_zero()
                #assert y.get_min_atom().ell_zero() == ell - pos

                #assert w.brion_length_b() == ell - len(cdes)
                #assert w.brion_length_b() == w.ell_zero() + pos

                #assert w.brion_length_c() == len(y.pair()) + len(cdes)
                #assert w.brion_length_c() == y.count_positive_ascents() - w.ell_zero()
                #assert len(cdes) * 2 == len(y.neg()) - len(nneg)


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


def test_nneg(n=4):
    for w in SignedPermutation.involutions(n):
        for atom in w.get_atoms():
            o = atom.inverse().oneline

            ndes, fix, neg = atom.ndes(), atom.nfix(), atom.nneg()

            sh = {(a, -a) for a in neg}
            sh |= {(a, -b) for a, b in ndes if 0 < a < -b}
            sh |= {(b, -a) for a, b in ndes if 0 < a < -b}
            sh = tuple(sorted(sh))

            for a, b in ndes:
                if a < -b:
                    neg += (a, b)
            neg = tuple(sorted([abs(i) for i in neg]))

            cfix = tuple(i for i in w.fixed_points() if i > 0)
            cneg = tuple(i for i in w.negated_points() if i > 0)
            csh = tuple(sorted(atom.shape()))

            print(w, ':', o, '->', csh, '==', sh)
            print('  fixed points:', cfix, '==', fix)
            print('negated points:', cneg, '==', neg)
            print()
            assert fix == cfix
            assert neg == cneg
            assert csh == sh
