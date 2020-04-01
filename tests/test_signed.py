from signed import SignedPermutation


def test_fpf_involution_words(n=4):
    for w in SignedPermutation.involutions(n):
        words = set(w.get_fpf_involution_words())
        if words:
            print(w, '->', len(words))


def test_get_atoms(n=5):
    i = list(SignedPermutation.involutions(n))
    for w in i:
        for a in w.get_atoms():
            assert len(a) == w.involution_length()
            assert a.inverse() % a == w
