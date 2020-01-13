from signed import SignedPermutation


def test_fpf_involution_words(n=4):
    for w in SignedPermutation.involutions(n):
        words = set(w.get_fpf_involution_words())
        if words:
            print(w, '->', len(words))
