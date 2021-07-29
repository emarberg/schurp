from permutations import Permutation
from pipedreams import Pipedream
from schubert import InvSchubert, FPFSchubert
from partitions import Partition, StrictPartition
from words import fpf_insert, eg_insert
from tableaux import Tableau
import pytest


def sigma_k(mu, k):
    d = {(i + 1, j + 1) for i in range(len(mu)) for j in range(mu[i])}
    d = {(i, j) for (i, j) in d if i > k and j > k}
    code = []
    for (i, j) in d:
        while i - 1 >= len(code):
            code += [0]
        code[i - 1] += 1
    return Permutation.from_code(code)


def test_sigma(n=8, k=0):
    lam = tuple([n - i for i in range(0, n)])
    for mu in [[2 * k + n - i for i in range(0, 2 *k + n)]]: #Partition.all(n):
        mu = tuple(mu)
#        for k in range(0, n + 1):
        sigma = sigma_k(mu, k)
        if sigma.is_identity():
            continue
        n = sigma.rank
        c = 0
        for p in sigma.get_pipe_dreams():
            c += 1
            print(p)
            w = p.words()
            w = [[n - i for i in a] for a in w]
            print(w)
            Q = eg_insert(*w)[1]
            print(Q)
            assert Q.is_k_flagged(k)
        sigma.print_rothe_diagram()
        print()
        print('k =', k, 'mu =', mu, 'sigma =', sigma)
        print('pipe dreams:', c)
        print('  k flagged:', len(Tableau.k_flagged(k, lam)))
        input('\n')


def sigma_fpf(mu, k):
    d = {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(mu[i])}
    d = {(j + 1, i) for (i, j) in d if i > k}
    code = []
    for (i, j) in d:
        while i - 1 >= len(code):
            code += [0]
        code[i - 1] += 1
    return Permutation.from_fpf_involution_code(code)


def test_sigma_fpf(n=8, k=0):
    lam = tuple([n - i for i in range(0, n, 2)])
#    for mu in [[2 * k + n - i for i in range(0, n + 2 * k, 2)]]:  # StrictPartition.all(n):
#        mu = tuple(mu)
#        for k in range(2, n + 1, 2):
#        sigma = sigma_fpf(mu, k)
    oneline = [i - 1 if i % 2 == 0 else i + 1 for i in range(1, 2 * k + 1)] + [2 * n + 2 * k - i for i in range(2 * n)]
    sigma = Permutation(*oneline)
    # if sigma.is_identity():
    #    continue
    n = sigma.rank
    c = 0
    for p in sigma.get_fpf_involution_pipe_dreams():
        c += 1
        print(p)
        w = p.words()
        # print(w)
        # print(fpf_insert(*w)[1])
        v = [[n - i for i in a] for a in w[1 + k // 2:]]
        print(v)
        Q = fpf_insert(*v)[1]
        print(Q)
        assert Q.is_shifted_k_flagged(2 * k)
    sigma.print_fpf_rothe_diagram()
    print()
    print('k =', k, 'sigma =', sigma.cycle_repr())
    print('  fpf pipe dreams:', c)
    print('shifted k flagged:', len(Tableau.shifted_k_flagged(2 * k, lam)))
    print('mu =', lam)
    print(Tableau.shifted_k_flagged(2 * k, lam))
    input('\n')


def count_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    w = Permutation(*oneline)
    return w.count_pipe_dreams(), w


def count_involution_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    w = Permutation(*oneline)
    return w.count_involution_pipe_dreams(), w


def count_fpf_pipe_dreams(shift, w):
    if type(w) != Permutation:
        w = Permutation.longest_element(w)
    oneline = list(range(1, shift + 1)) + [i + shift for i in w.oneline]
    w = Permutation(*oneline)
    for i in range(1, shift, 2):
        w *= Permutation.s_i(i)
    return w.count_fpf_involution_pipe_dreams(), w


def test_pipe_dream_counts():
    def expected(k, n):
        ans = 1
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                ans *= 2 * k + i + j - 1
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                ans = ans // (i + j - 1)
        return ans

    for n in range(5):
        print()
        for k in range(5):
            e = expected(k, n)
            f, w = count_pipe_dreams(k, n)
            s = ''.join([str(i) for i in range(n, 0, -1)])
            print('1^%i x %s =' % (k, s), w, ':', e, '=', f)
            assert e == f


def test_inv_pipe_dream_counts():
    def expected(k, n):
        p = n // 2
        q = p if n % 2 == 0 else p + 1
        ans = 1
        for i in range(1, p + 1):
            for j in range(1, q + 1):
                ans *= i + j + k - 1
        for i in range(1, p + 1):
            for j in range(1, q + 1):
                assert ans % (i + j - 1) == 0
                ans = ans // (i + j - 1)
        return ans

    for n in range(6):
        print()
        print('n =', n)
        for k in range(6):
            e = expected(k, n)
            f, w = count_involution_pipe_dreams(k, n)
            s = ''.join([str(i) for i in range(n, 0, -1)])
            print('1^%i x %s =' % (k, s), w, ':', e, '=', f)
            assert e == f


def test_gr_inv_pipe_dream_counts():
    def expected(k, n):
        k = k // 2 if k % 2 == 0 else (k - 1) // 2
        n = n + 1
        return count_pipe_dreams(k, n)[0]

    for n in range(6):
        print()
        print('n =', n)
        for k in range(6):
            w = Permutation()
            for i in range(1, n + 1):
                w *= Permutation.transposition(i + k, n + i + k)
            e = expected(k, n)
            f = len(list(w.get_involution_pipe_dreams()))
            s = ''.join([str(i) for i in w.oneline[k:]])
            print('1^%i x %s =' % (k, s), w, ':', e, '=', f)
            assert e == f


def test_fpf_pipe_dream_counts():
    def expected(k, n):
        k = k // 2
        n = n // 2
        c, _ = count_pipe_dreams(k, n)
        ans = 1
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                if i != j:
                    ans *= 2 * k + i + j - 1
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                if i != j:
                    ans = ans // (i + j - 1)
        assert ans == c ** 2
        return ans

    up = 7
    for n in range(0, up, 2):
        print()
        print('n =', n)
        for k in range(0, up, 2):
            e = expected(k, n)
            f, w = count_fpf_pipe_dreams(k, n)
            s = ''.join([str(i) for i in range(n, 0, -1)])
            print('1^%i x %s =' % (k, s), w, ':', e, '=', f)
            assert e == f


def test_bottom_pipedream():
    w = Permutation(3, 1, 4, 6, 5, 2)
    p = w.get_bottom_pipe_dream()
    assert p.crossings == {(1, 1), (1, 2), (3, 1), (4, 1), (4, 2), (5, 1)}


def test_ladder_moves():
    w = Permutation(1, 4, 3, 2)
    p = w.get_bottom_pipe_dream()

    q = Pipedream({(1, 2), (2, 1), (2, 2)})
    r = Pipedream({(1, 3), (2, 1), (3, 1)})
    print(p)
    print(q)
    print(r)
    print()
    for x in p.ladder_moves():
        print(x)
    assert set(p.ladder_moves()) == {q, r}

    test = set(p.upper_ladder_interval())
    assert test == set(w.get_pipe_dreams())
    expected = {
        p,
        q,
        r,
        Pipedream({(1, 2), (1, 3), (3, 1)}),
        Pipedream({(1, 2), (1, 3), (2, 2)}),
    }
    assert test == expected


def test_ladder_moves_span():
    n = 5
    for w in Permutation.all(n):
        a = set(w.get_pipe_dreams())
        b = {
            Pipedream.from_word(*dream)
            for word in w.get_reduced_words()
            for dream in w._get_pipe_dreams_helper(word)
        }
        assert a == b


def test_involutions():
    n = 6
    for i, perm in enumerate(Permutation.involutions(n)):
        print('. . .', i, perm)
        dreams = perm.get_involution_pipe_dreams()
        test_s = sum([dream.inv_monomial() for dream in dreams])
        s = InvSchubert.get(perm)
        assert s == test_s


def test_fpf_involutions():
    n = 6
    for i, perm in enumerate(Permutation.fpf_involutions(n)):
        print('. . .', i, perm)
        dreams = perm.get_fpf_involution_pipe_dreams()
        test_s = sum([dream.fpf_monomial() for dream in dreams])
        s = FPFSchubert.get(perm)
        assert s == test_s


def test_involution_ladder_moves():
    w = Permutation(1, 2, 3, 6, 5, 4)
    p = w.get_bottom_involution_pipe_dream()
    assert p == Pipedream({(5, 1), (4, 1)})

    q = Pipedream({(4, 1), (4, 2)})
    print(p)
    print()
    for d in p.involution_ladder_moves():
        print(d)
        print()
    assert q in p.involution_ladder_moves()

    p = Pipedream({(3, 1), (3, 3), (4, 1)})
    q = Pipedream({(3, 1), (3, 2), (3, 3)})
    assert q not in p.involution_ladder_moves()


def print_discrepancy(a, b, w):
    if a != b:
        print(w)
        print()
        for p in a & b:
            print(p)
            print()
        print('Extra:\n')
        for p in b - a:
            print(p)
            print()
        print('Missing:\n')
        for p in a - b:
            print(p)
            print()
        print()


def test_involution_ladder_moves_span():
    n = 6
    for w in Permutation.involutions(n):
        a = set(w._get_involution_pipe_dreams_slow())
        b = set(w.get_involution_pipe_dreams())
        print_discrepancy(a, b, w)
        assert b.issubset(a)
        assert a == b

        a = set(w._get_involution_pipe_dreams_slow(True))
        b = set(w.get_involution_pipe_dreams(True))
        print_discrepancy(a, b, w)
        assert a == b


# def test_involution_chute_moves():
#     p = Pipedream({(1, 4), (1, 5)})
#     q = Pipedream({(1, 4), (2, 4)})
#     assert q in p.involution_chute_moves()


# def test_involution_chute_moves_span():
#     n = 4
#     for w in Permutation.involutions(n):
#         a = set(w.get_involution_pipe_dreams(extended=True))
#         b = set(w.get_bottom_involution_pipe_dream().involution_chute_span())
#         print_discrepancy(a, b, w)
#         assert a == b


def test_fpf_ladder_moves():
    p = Pipedream({(3, 2), (4, 2)})
    q = Pipedream({(3, 1), (3, 2)})
    print(p)
    print()
    for d in p.fpf_involution_ladder_moves():
        print(d)
        print()
    assert q in p.fpf_involution_ladder_moves()

    w = Permutation(2, 1, 6, 5, 4, 3)
    p = w.get_bottom_fpf_pipe_dream()
    assert p == Pipedream({(4, 1), (5, 1)})

    q = Pipedream({(3, 2), (5, 1)})
    print(p)
    print()
    for d in p.fpf_involution_ladder_moves():
        print(d)
        print()
    assert {q} == set(p.fpf_involution_ladder_moves())


def test_fpf_ladder_moves_span():
    # have checked up to n = 10
    n = 6
    for w in [
        Permutation.from_fpf_involution_word(4, 3, 6, 5, 4),
        Permutation.from_fpf_involution_word(2, 1, 4, 3, 2)
    ] + list(Permutation.fpf_involutions(n)):
        a = set(w._get_fpf_involution_pipe_dreams_slow())
        b = set(w.get_fpf_involution_pipe_dreams())
        print_discrepancy(a, b, w)
        assert b.issubset(a)
        assert a == b

        a = set(w._get_fpf_involution_pipe_dreams_slow(extended=True))
        b = set(w.get_fpf_involution_pipe_dreams(extended=True))
        print_discrepancy(a, b, w)
        assert a == b
