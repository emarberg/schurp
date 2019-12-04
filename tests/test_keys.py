from keys import (
    weak_compositions,
    monomial,
    sorting_permutation,
    strict_weak_compositions,
    key, atom,
    p_key, p_atom,
    q_key, q_atom,
    has_distinct_parts,
)
from symmetric import FPFStanleyExpander, InvStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert, X
from permutations import Permutation
from collections import defaultdict


def test_weak_compositions():
    assert set(weak_compositions(0, 0)) == {()}
    assert set(weak_compositions(4, 2)) == {(4, 0), (1, 3), (2, 2), (3, 1), (0, 4)}
    assert len(list(weak_compositions(4, 2))) == 5


def test_sorting_permutation():
    assert sorting_permutation((1, 0, 2, 1)) == (2, 1, 3)
    assert sorting_permutation((0, 0, 0, 0)) == tuple()
    assert sorting_permutation((1, 2, 3)) == (1, 2, 1)


def test_ordinary_key():
    weak_comp = (1, 0, 2, 1)
    expected_key = \
        monomial((2, 1, 1, 0)) + \
        monomial((1, 2, 1, 0)) + \
        monomial((1, 1, 2, 0)) + \
        monomial((2, 1, 0, 1)) + \
        monomial((2, 0, 1, 1)) + \
        monomial((1, 2, 0, 1)) + \
        monomial((1, 1, 1, 1)) + \
        monomial((1, 0, 2, 1))
    actual_key = key(weak_comp)
    print(expected_key)
    print()
    print(actual_key)
    assert expected_key == actual_key


def test_atom():
    for n in range(5):
        for k in range(5):
            for alpha in weak_compositions(n, k):
                kappa = atom(alpha)
                print(alpha, kappa)
                assert kappa.is_positive()


def test_p_atom():
    for n in range(8):
        for k in range(8):
            for alpha in strict_weak_compositions(n, k):
                kappa = p_atom(alpha)
                # print(alpha, max(kappa))
                assert kappa.is_positive()
                assert kappa.is_not_laurent_polynomial()
    print('success')


def test_q_atom():
    for n in range(8):
        for k in range(8):
            for alpha in strict_weak_compositions(n, k):
                kappa = q_atom(alpha)
                # print(alpha, max(kappa))
                assert kappa.is_positive()
                assert kappa.is_not_laurent_polynomial()
    print('success')


def schur(partition):
    assert all(partition[i] >= partition[i + 1] for i in range(len(partition) - 1))
    w = Permutation.get_grassmannian(*partition)
    n = len(partition)
    return Schubert.get(w).truncate(n)


def test_schur():
    assert schur((3, 1)) == key((1, 3, 0, 0))


def schurp(partition):
    assert all(partition[i] > partition[i + 1] for i in range(len(partition) - 1))
    w = Permutation.get_inv_grassmannian(*partition)
    w = w.shift(w.rank)
    return InvSchubert.get(w)


def icode(w):
    diagram = w.involution_rothe_diagram()
    ans = w.rank * [0]
    for i, j in diagram:
        ans[j - 1] += 1
    return tuple(ans)


def fcode(w):
    diagram = w.fpf_rothe_diagram()
    ans = w.rank * [0]
    for i, j in diagram:
        ans[j - 1] += 1
    return tuple(ans)


def tuplize(g):
    beta = []
    for i in g:
        while i > len(beta):
            beta += [0]
        beta[i - 1] = g[i]
    beta = tuple(beta)
    return beta


def get_exponents(kappa):
    betas = []
    for g in kappa:
        beta = tuplize(g)
        betas.append(beta)
    return betas


def test_leading_key(m=5):
    for n in range(m):
        for k in range(m):
            for alpha in weak_compositions(n, k):
                while alpha and alpha[-1] == 0:
                    alpha = alpha[:-1]
                kappa = key(alpha)
                betas = get_exponents(kappa)
                beta = min(betas)
                print(alpha, beta, betas, kappa)
                print()
                assert beta == alpha


def decompose_into_atoms(kappa):
    def dict_from_tuple(beta):
        mon = X(0)**0
        for i, a in enumerate(beta):
            mon *= X(i + 1)**a
        return max(mon)
    #
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * atom(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def test_p_atom_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in weak_compositions(n, k, reduced=True):
                kappa = p_atom(alpha)
                dec = decompose_into_atoms(kappa)
                if has_distinct_parts(alpha):
                    if kappa != 0:
                        print(alpha)
                        print(set(dec.values()))
                        print()
                    assert min({0} | set(dec.values())) >= 0
                else:
                    print(alpha)
                    print(set(dec.values()))
                    print()


def test_q_atom_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in weak_compositions(n, k, reduced=True):
                e = 2 ** len(set(alpha) - {0})
                kappa = q_atom(alpha)
                dec = decompose_into_atoms(kappa)
                if has_distinct_parts(alpha):
                    if kappa != 0:
                        print(alpha)
                        print({d // e for d in dec.values()})
                        print()
                        assert all(d % e == 0 for d in dec.values())
                    assert min({0} | set(dec.values())) >= 0
                else:
                    print(alpha)
                    print({d // e for d in dec.values()})
                    print()
                    assert all(d % e == 0 for d in dec.values())


def decompose_into_keys(kappa):
    def dict_from_tuple(beta):
        mon = X(0)**0
        for i, a in enumerate(beta):
            mon *= X(i + 1)**a
        return max(mon)
    #
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * key(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def test_p_key_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in weak_compositions(n, k, reduced=True):
                kappa = p_key(alpha)
                dec = decompose_into_keys(kappa)
                if has_distinct_parts(alpha):
                    if set(dec.values()) != {1}:
                        print(alpha)
                        print(set(dec.values()))
                        print()
                    assert min({0} | set(dec.values())) >= 0
                else:
                    v = set(dec.values())
                    assert not v or min(v) < 0
                if kappa == 0:
                    print(alpha, 'f =', kappa)
                    print(dec)
                    print()


def test_q_key_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in weak_compositions(n, k, reduced=True):
                kappa = q_key(alpha)
                dec = decompose_into_keys(kappa)
                if has_distinct_parts(alpha):
                    expected = {2 ** len(set(alpha) - {0})}
                    if set(dec.values()) != expected:
                        print(alpha)
                        print(set(dec.values()))
                        print()
                    assert set(dec.values()) == expected
                else:
                    print(alpha)
                    print(set(dec.values()))
                    print()
                if kappa == 0:
                    print(alpha, 'f =', kappa)
                    print(dec)
                    print()


def test_leading_p_key(m=5):
    for n in range(2 * m):
        for k in range(m):
            valuesdict = defaultdict(list)
            for alpha in strict_weak_compositions(n, k):
                while alpha and alpha[-1] == 0:
                    alpha = alpha[:-1]
                kappa = p_key(alpha)
                valuesdict[kappa].append(alpha)
                betas = get_exponents(kappa)
                beta = min([b for b in betas if len(b) == len(alpha)])
                assert alpha in betas
                # if beta != alpha:
                print(beta == alpha, alpha, beta, [b for b in betas if len(b) == k])
            print()
            print()
            print('repetitions')
            for k, v in valuesdict.items():
                if len(v) == 1:
                    continue
                print(v)
                print()
            print()
        #        assert beta == alpha


def test_leading_q_key(m=5):
    for n in range(2 * m):
        for k in range(m):
            valuesdict = defaultdict(list)
            for alpha in strict_weak_compositions(n, k):
                while alpha and alpha[-1] == 0:
                    alpha = alpha[:-1]
                kappa = q_key(alpha)
                betas = get_exponents(kappa)
                beta = min([b for b in betas if len(b) == len(alpha)])
                valuesdict[beta].append(alpha)
                assert alpha in betas
                print(q_key(beta) == kappa, beta == alpha, alpha, beta)#, [b for b in betas if len(b) == k])
            print()
            print()
            if any(len(v) > 1 for v in valuesdict.values()):
                print('repetitions')
                print()
            for k, v in valuesdict.items():
                if len(v) == 1:
                    continue
                print(v, len({q_key(x) for x in v}))
                print()
            print()


def test_inv_schubert(n=4):
    def isvex(w):
        f = InvStanleyExpander(w).expand()
        if len(f) > 1:
            return False
        mu = list(f.keys())[0].mu
        coeff = list(f.values())[0] * 2 ** w.number_two_cycles()
        return coeff == 2 ** len(mu)

    i = set(Permutation.involutions(n))
    s = {w: InvSchubert.get(w) * 2**w.number_two_cycles() for w in i}
    # x = {w for w in i if q_key(icode(w)) == s[w]}
    # return i, s, x
    m = max({w.involution_length() for w in i})
    e = {alpha: q_key(alpha) for l in range(m + 1) for alpha in strict_weak_compositions(l, n)}
    v = set(e.values())
    x = {w for w in s if s[w] in v}
    vex = {w for w in i if w.is_vexillary()}
    ivex = {w for w in i if isvex(w)}
    assert vex == ivex
    assert set(x).issubset(vex)
    # for w in vex:
    #     a = icode(w)
    #     if q_key(a) != s[w]:
    #         assert w not in x
    #         w.print_involution_rothe_diagram()
    #         # print([b for b in e if e[b] == s[w]])
    #         print()
    return i, s, x, e, vex


def test_fpf_schubert(n=4):
    def isvex(w):
        f = FPFStanleyExpander(w).expand()
        if len(f) > 1:
            return False
        if set(f.values()) != {1}:
            return False
        return True
    i = set(Permutation.fpf_involutions(n))
    s = {w: FPFSchubert.get(w) for w in i}
    # x = {w for w in i if p_key(fcode(w)) == s[w]}
    # return i, s, x
    m = max({w.fpf_involution_length() for w in i})
    e = {alpha: p_key(alpha) for l in range(m + 1) for alpha in strict_weak_compositions(l, n)}
    v = set(e.values())
    x = {w for w in s if s[w] in v}
    vex = {w for w in i if isvex(w)}
    for w in set(x) - vex:
        f = FPFStanleyExpander(w)
        v = f.expand()
        print(w)
        print(len(v))
        print()
        input('???')
    for w in x:
        a = fcode(w)
        if p_key(a) != s[w]:
            print(s[w])
            print(a)
            w.print_fpf_rothe_diagram()
            print([b for b in e if e[b] == s[w]])
            print()
            input('!!!')
    return i, s, x, e, vex
