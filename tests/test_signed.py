from signed import SignedPermutation
from polynomials import X


def test_fpf_involution_words(n=4):
    for w in SignedPermutation.involutions(n):
        words = set(w.get_fpf_involution_words())
        if words:
            print(w, '->', len(words))


def test_get_atoms(n=4):
    i = list(SignedPermutation.involutions(n))
    for w in i:
        for a in w.get_atoms():
            assert len(a) == w.involution_length()
            assert a.inverse() % a == w


def test_get_abs_fpf_atoms(n=6):
    for y in SignedPermutation.abs_fpf_involutions(n):
        s = set()
        for w in y.get_atoms():
            winv = w.inverse()
            if all(winv(i - 1) > winv(i) for i in range(2, y.rank, 2)):
                s.add(SignedPermutation.get_minimal_fpf_involution(n) * w)
        assert s == set(y.get_fpf_atoms())


def test_brion_length(n=6):
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
