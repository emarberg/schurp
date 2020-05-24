from even import EvenSignedPermutation


def test_get_minimal_fpf_involution():
    cls = EvenSignedPermutation
    assert cls.get_minimal_fpf_involution(1) == cls(1)
    assert cls.get_minimal_fpf_involution(2) == cls(2, 1)
    assert cls.get_minimal_fpf_involution(3) == cls(1, 3, 2)
    assert cls.get_minimal_fpf_involution(4) == cls(2, 1, 4, 3)


def test_get_involution_word():
    s = EvenSignedPermutation.s_i(0, 3)
    t = EvenSignedPermutation.s_i(2, 3)
    assert (t * s * t).get_involution_word() in {(0, 2), (2, 0)}
    assert (s * t * s).get_involution_word() in {(0, 2), (2, 0)}
    assert set((t * s * t).get_involution_words()) == {(0, 2), (2, 0)}
    assert set((s * t * s).get_involution_words()) == {(0, 2), (2, 0)}


def test_fpf_involution_words(n=4):
    for w in EvenSignedPermutation.fpf_involutions(n):
        words = set(w.get_fpf_involution_words())
        if words:
            print(w, '->', len(words))


def test_get_atoms(n=4):
    i = list(EvenSignedPermutation.involutions(n))
    for w in i:
        words = set()
        for a in w.get_atoms():
            assert len(a) == w.involution_length()
            assert a.inverse() % a == w
            words |= set(a.get_reduced_words())
        assert words == set(w.get_involution_words())


def test_get_twisted_atoms(n=4):
    i = list(EvenSignedPermutation.twisted_involutions(n))
    for w in i:
        words = set()
        for a in w.get_twisted_atoms():
            assert len(a) == w.twisted_involution_length()
            assert a.inverse().star() % a == w
            words |= set(a.get_reduced_words())
        assert words == set(w.get_twisted_involution_words())
