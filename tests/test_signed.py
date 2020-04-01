from signed import SignedPermutation


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


def test_get_fpf_atoms(n=6):
    for y in SignedPermutation.fpf_involutions(n):
        s = set()
        for w in y.get_atoms():
            winv = w.inverse()
            if all(winv(i - 1) > winv(i) for i in range(2, y.rank, 2)):
                s.add(SignedPermutation.get_minimal_fpf_involution(n) * w)
        assert s == set(y.get_fpf_atoms())
