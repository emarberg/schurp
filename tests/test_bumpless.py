from pipedreams import BumplessPipedream, SymmetricBumplessPipedream
from permutations import Permutation
from schubert import (
    DoubleSchubert, 
    DoubleGrothendieck,
    FPFSchubert, 
    InvSchubert,
    AltInvGrothendieck,
)


def inv_schubert_via_bumpless(w, strict):
    ans = 0
    for bpd in BumplessPipedream.from_involution(w, strict=strict):
        ans += bpd.inv_weight()
    return ans


def inv_grothendieck_via_bumpless(w, strict):
    ans = 0
    for bpd in BumplessPipedream.from_involution(w, reduced=False, strict=strict):
        if AltInvGrothendieck.beta in [-1, 1]:
            sgn = AltInvGrothendieck.beta**w.involution_length()
        else:
            sgn = AltInvGrothendieck.beta**(-w.involution_length())
        ans += sgn * bpd.inv_kweight()    
    return ans


def test_inv_bumpless(n=6):
    for strict in [True, False]:
        for w in Permutation.involutions(n):
            print(w, strict)
            expected = InvSchubert.get(w)
            frombpd = inv_schubert_via_bumpless(w, strict) 
            assert expected == frombpd

            expected = AltInvGrothendieck.get(w)
            frombpd = inv_grothendieck_via_bumpless(w, strict)
            assert expected == frombpd
        
        dreams = list(BumplessPipedream.from_involution(w))
        test = sum(d.inv_weight() for d in dreams)
        actual = InvSchubert.get(w)


def schubert_via_bumpless(w, strict):
    ans = 0
    for bpd in BumplessPipedream.from_permutation(w, strict=strict):
        ans += bpd.weight()
    return ans


def grothendieck_via_bumpless(w, strict):
    ans = 0
    for bpd in BumplessPipedream.from_permutation(w, reduced=False, strict=strict):
        if DoubleGrothendieck.beta in [-1, 1]:
            sgn = DoubleGrothendieck.beta**w.length()
        else:
            sgn = DoubleGrothendieck.beta**(-w.length())
        ans += sgn * bpd.kweight()    
    return ans


def fpf_schubert_via_bumpless(w):
    ans = 0
    for bpd in BumplessPipedream.from_fpf_involution_slow(w):
        ans += bpd.fpf_weight()
    return ans


def test_schubert(n=5):
    for strict in [True, False]:
        for w in Permutation.all(n):
            expected = DoubleSchubert.get(w)
            frombpd = schubert_via_bumpless(w, strict)    
            assert expected == frombpd

            expected = DoubleGrothendieck.get(w)
            frombpd = grothendieck_via_bumpless(w, strict)
            assert expected == frombpd


def test_fpf_schubert(n=6):
    for w in Permutation.fpf_involutions(n):
        expected = FPFSchubert.get(w)
        frombpd = fpf_schubert_via_bumpless(w)
        if expected != frombpd:
            print(w)
            print(BumplessPipedream.from_fpf_involution_slow(w))
            print()
            print(' S_w =', expected)
            print()
            print('      ', frombpd)
            print()
            print('diff =', expected - frombpd)
            raise Exception


def test_symmetric(n=6):
    for w in Permutation.fpf_involutions(n):
        one = BumplessPipedream.from_fpf_involution_slow(w, n)
        two = SymmetricBumplessPipedream.from_fpf_involution(w, n)
        if one != two:
            print(w, n, ':', len(one), '=?=', len(two))
            print()
            print(one)
            print()
            print(two)



