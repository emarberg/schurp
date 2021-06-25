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


def recover_permutation_from_rank_table(tbl):
    # determine permutation from rank table
    oneline = []
    for i in range(len(tbl)):
        for j in range(len(tbl[i])):
            a = tbl[i][j - 1] if j > 0 else 0
            b = tbl[i - 1][j] if i > 0 else 0
            c = tbl[i - 1][j - 1] if i > 0 and j > 0 else 0
            if a == b == c == tbl[i][j] - 1:
                # add j+1 since j is between 0 and n - 1
                oneline += [j + 1]
    return Permutation(*oneline)


def upper_bound_chains(p, q):
    # if (p, q) == (4, 3) then this method returns the list of chains
    #
    # . . . .  . . . .  . . . .
    # . . . .  . . . .  . . x .
    # . . . .  . . x .  . x . .
    # . . x .  . x . .  x . . .
    #
    # for example
    ans = []
    for k in range(min(p, q)):
        ans += [tuple((p - k + j, q - j) for j in range(k + 1))]
    return ans


def recover_permutation_from_northeast_chains(ch):
    # input `ch` is collection corresponding to \cA_w

    # determine n
    n = max([0] + [a for chain in ch for pair in chain for a in pair]) + 1
    # determine rank table from chains
    tbl = [[None for _ in range(n)] for _ in range(n)]
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            e = [chain for chain in upper_bound_chains(i, j) if chain in ch]
            tbl[i - 1][j - 1] = len(e[0]) - 1 if e else min(i, j)
    return recover_permutation_from_rank_table(tbl)


def test_recover_from_chains(n=5):
    for w in Permutation.all(n):
        ch = w.all_northeast_chains()
        assert w == recover_permutation_from_northeast_chains(ch)


def upper_bound_reflected_chains(p, q):
    # if (p, q) == (6, 4) then this method returns the *reflections* of the chains
    #
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . . x . .
    #
    # . . . x . .
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . x . . .
    #
    # . . . x . .
    # . . . . . .
    # . . . . . .
    # . . . . . .
    # . . x . . .
    # . x . . . .
    #
    # . . . x . .
    # . . x . . .
    # . . . . . .
    # . . . . . .
    # . x . . . .
    # x . . . . .
    #
    # for example
    assert p > q
    ans = []
    for k in range(min(p, q)):
        m = 1 + k // 2
        chain = []
        for j in range(k + 1):
            x = (p - j) if j < m else (k + 1 - j)
            y = q - k + j
            assert x != y
            chain += [(y, x) if x < y else (x, y)]
        ans += [tuple(sorted(chain))]
    return ans


def recover_involution_from_reflected_chains(ch):
    # input `ch` is collection corresponding to \ssA_z

    # determine n
    n = max([0] + [a for chain in ch for pair in chain for a in pair])
    n = n + 1 if n % 2 != 0 else n

    # determine below diagonal part of rank table from reflected chains
    tbl = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(1, n + 1):
        for j in range(1, i):
            e = [chain for chain in upper_bound_reflected_chains(i, j) if chain in ch]
            tbl[i - 1][j - 1] = len(e[0]) - 1 if e else min(i, j)
            # fill in transposed entry by symmetry
            tbl[j - 1][i - 1] = tbl[i - 1][j - 1]

    # fill is diagonal part of rank table
    for i in range(1, n):
        a = tbl[i - 1][i - 1]
        b = tbl[i][i - 1]
        tbl[i][i] = b + 1 if a < b else a

    return recover_permutation_from_rank_table(tbl)


def test_recover_from_reflected_chains(n=6):
    for z in Permutation.fpf_involutions(n):
        ch = z.all_reflected_northeast_chains()
        assert z.fpf_trim() == recover_involution_from_reflected_chains(ch)
