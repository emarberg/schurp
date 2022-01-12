from pipedreams import BumplessPipedream, SymmetricBumplessPipedream
from permutations import Permutation
from schubert import DoubleSchubert, FPFSchubert


def schubert_via_bumpless(w):
    ans = 0
    for bpd in BumplessPipedream.from_permutation(w):
        ans += bpd.weight()
    return ans


def fpf_schubert_via_bumpless(w):
    ans = 0
    for bpd in BumplessPipedream.from_fpf_involution(w):
        ans += bpd.fpf_weight()
    return ans


def test_schubert(n=4):
    for w in Permutation.all(n):
        expected = DoubleSchubert.get(w)
        frombpd = schubert_via_bumpless(w)
        if expected != frombpd:
            print(w)
            print(BumplessPipedream.from_permutation(w))
            print()
            print(' S_w =', expected)
            print()
            print('      ', frombpd)
            print()
            print('diff =', expected - frombpd)
            raise Exception


def test_fpf_schubert(n=6):
    for w in Permutation.fpf_involutions(n):
        expected = FPFSchubert.get(w)
        frombpd = fpf_schubert_via_bumpless(w)
        if expected != frombpd:
            print(w)
            print(BumplessPipedream.from_fpf_involution(w))
            print()
            print(' S_w =', expected)
            print()
            print('      ', frombpd)
            print()
            print('diff =', expected - frombpd)
            raise Exception


def test_symmetric(n=8):
    for w in Permutation.fpf_involutions(n):
        one = BumplessPipedream.from_fpf_involution_slow(w, n)
        two = SymmetricBumplessPipedream.from_fpf_involution(w, n)
        if one != two:
            print(w, n, ':', len(one), '=?=', len(two))
            print()
            print(one)
            print()
            print(two)



