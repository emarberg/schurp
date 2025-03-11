from permutations import Permutation


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