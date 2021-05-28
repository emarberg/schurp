from permutations import Permutation


def test_essential_chains(n=4):
    print('n =', 4)
    print()
    for i, w in enumerate(Permutation.all(n)):
        print(i + 1, ':', w)
        chains = w.all_northeast_chains()
        ess = w.essential_northeast_chains()
        for b in chains:
            assert any(set(a).issubset(set(b)) for a in ess)
    print()


def test_reflected_chains(n=4):
    print('n =', 4)
    print()
    for i, w in enumerate(Permutation.fpf_involutions(n)):
        print(i + 1, ':', w)
        chains = w.all_reflected_northeast_chains()
        ess = w.essential_reflected_northeast_chains()
        for b in chains:
            assert any(set(a).issubset(set(b)) for a in ess)
    print()
