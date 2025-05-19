from permutations import Permutation
from signed import SignedPermutation


def get_inv_poset(z):
    n = z.rank
    atoms = z.get_atoms()
    upper_poset = {}
    for w in z.all(n):
        for a in atoms:
            if a.strong_bruhat_less_equal(w):
                upper_poset[w] = set()
    for x in upper_poset:
        for y in upper_poset:
            if x.strong_bruhat_less_equal(y):
                upper_poset[y].add(x)
    return upper_poset


def get_fpf_poset(z):
    n = z.rank
    atoms = z.get_fpf_atoms()
    upper_poset = {}
    for w in Permutation.all(n):
        for a in atoms:
            if a.strong_bruhat_less_equal(w):
                upper_poset[w] = set()
    for x in upper_poset:
        for y in upper_poset:
            if x.strong_bruhat_less_equal(y):
                upper_poset[y].add(x)
    return upper_poset


def get_twisted_poset(z, n):
    atoms = z.get_twisted_atoms(n)
    upper_poset = {}
    for w in Permutation.all(n):
        for a in atoms:
            if a.strong_bruhat_less_equal(w):
                upper_poset[w] = set()
    for x in upper_poset:
        for y in upper_poset:
            if x.strong_bruhat_less_equal(y):
                upper_poset[y].add(x)
    return upper_poset
                

def get_mobius(upper_poset):
    ans = {}

    def mobius(x):
        if x not in ans:
            mu = 0
            for z in upper_poset[x] - {x}:
                mu += mobius(z)
            ans[x] = 1 - mu
        return ans[x]

    for x in upper_poset:
        mobius(x)
    return ans


def test_twisted(n=4):
    for z in Permutation.twisted_involutions(n):
        print('z =', z)
        upper_poset = get_twisted_poset(z, n)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_twisted_hecke_atoms(n))
        actual = {x: (-1)**(x.length() - z.twisted_involution_length(n)) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_fpf(n=4):
    # tests formula in Theorem 3 of https://arxiv.org/pdf/0902.1930
    for z in Permutation.fpf_involutions(n):
        print('z =', z)
        upper_poset = get_fpf_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_symplectic_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.fpf_involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_inv(n=4):
    for z in Permutation.involutions(n):
        print('z =', z)
        upper_poset = get_inv_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
                # print('  ', x, expected[x], 'vs', actual.get(x, '*'))
        assert actual == expected


def test_fpf(n=4):
    for z in SignedPermutation.fpf_involutions(n):
        print('z =', z)
        upper_poset = get_fpf_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_symplectic_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.fpf_involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected


def test_binv(n=4):
    for z in SignedPermutation.involutions(n):
        print('z =', z)
        upper_poset = get_inv_poset(z)
        mobius = get_mobius(upper_poset)
        hecke = set(z.get_involution_hecke_atoms())
        actual = {x: (-1)**(x.length() - z.involution_length()) for x in hecke}
        expected = {}
        for x in upper_poset:
            if mobius[x] != 0:
                expected[x] = mobius[x]
        assert actual == expected