from keys import (
    weak_compositions,
    monomial,
    sorting_permutation,
    symmetric_weak_compositions,
    key, atom,
    q_power,
    p_key, p_atom,
    q_key, q_atom,
    has_distinct_parts,
    symmetric_half, symmetric_double,
    get_exponents,
    decompose_into_keys, decompose_into_atoms, dict_from_tuple,
    symmetric_composition_from_row_column_counts, symmetric_halves
)
from symmetric import FPFStanleyExpander, InvStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert, X
from permutations import Permutation
from collections import defaultdict
from words import Word
from partitions import Partition, StrictPartition
from tableaux import Tableau


def test_symmetric_composition_from_row_column_counts():
    assert symmetric_composition_from_row_column_counts((), ()) == ()

    rc = (5, 0, 0, 2, 0)
    cc = (1, 1, 1, 2, 2)
    alpha = symmetric_composition_from_row_column_counts(rc, cc)
    assert alpha == (5, 1, 1, 3, 2)

    assert symmetric_composition_from_row_column_counts((2,), (1, 1)) == (2, 1)
    assert symmetric_composition_from_row_column_counts((2, 2), (1, 2, 0, 1)) == (2, 3, 0, 1)


def test_symmetric_halves(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                a, b = symmetric_halves(alpha)
                beta = symmetric_composition_from_row_column_counts(a, b)
                print(alpha)
                print(beta)
                print(a, b)
                print()
                assert alpha == beta


def partitions(n):
    for mu in Partition.all(n):
        yield tuple(mu.parts)


def strict_partitions(n):
    for mu in StrictPartition.all(n):
        yield tuple(mu.parts)


def symmetric_compositions(n):
    for mu in strict_partitions(n):
        m = mu[0] if mu else 0
        shape = {(i, j) for i in range(1, len(mu) + 1) for j in range(i, i + mu[i - 1])}
        shape |= {(j, i) for (i, j) in shape}
        shape = tuple(sorted(shape))
        #
        seen = set()
        toadd = {shape}
        while toadd:
            next_toadd = set()
            for mu in toadd:
                yield Tableau({b: 1 for b in mu})
                seen.add(mu)
                for i in range(1, m):
                    s = Permutation.s_i(i)
                    nu = tuple(sorted((s(a), s(b)) for (a, b) in mu))
                    if nu not in seen:
                        next_toadd.add(nu)
            toadd = next_toadd


def eg_insert(sequence):
    return Word(*sequence).eg_insert()


def o_eg_insert(sequence):
    return Word(*sequence).involution_insert()


def sp_eg_insert(sequence):
    return Word(*sequence).fpf_insert()


def test_insertion_definition(n=2):
    for w in Permutation.all(n):
        dictionary = {}
        pipedreams = list(w.get_pipe_dreams())
        for dream in pipedreams:
            word = dream.word()
            p, q = eg_insert(word)
            if p not in dictionary:
                dictionary[p] = dream.monomial()
            else:
                dictionary[p] += dream.monomial()
        for p in dictionary:
            f = dictionary[p]
            d = decompose_into_keys(f)
            alpha = tuple(reversed(p.weight()))
            beta = list(d)[0]
            if alpha != beta:
                print('w =', w)
                print()
                print(p)
                print(f, '=', d)
                print()
                print(alpha, '?=', beta)
                print()
                print()
            assert len(d) == 1
            assert d[beta] == 1


def test_p_insertion_definition(n=2):
    discovered = []
    dictionary = {}
    for w in Permutation.fpf_involutions(n):
        pipedreams = list(w.get_fpf_involution_pipe_dreams())
        for dream in pipedreams:
            word = dream.word()
            p, q = sp_eg_insert(word)
            if p not in dictionary:
                dictionary[(w, p)] = dream.fpf_monomial()
            else:
                dictionary[(w, p)] += dream.fpf_monomial()
    for w, p in sorted(dictionary, key=lambda t: t[1].partition()):
        f = dictionary[(w, p)]
        d = decompose_into_keys(f)
        betas = sorted(d, key=lambda x: (len(x), tuple(sorted(x))))
        betas = [b for b in betas if is_p_special(b)]
        success = False
        for b in betas:
            if p_key(b) == f:
                success = True
                discovered.append(b)
        if not success:
            print('FAILED')
        else:
            continue
            print('SUCCESS')
        print()
        print('w =', w)
        print()
        print(p)
        print('  key expansion =', d)
        print()
        print('trials:')
        for beta in betas:
            print()
            print('  key expansion(', beta, ') =', decompose_into_keys(p_key(beta)))
        print()
        print()
        print()
        print()
        print()
    print('discovered:')
    for beta in sorted(discovered):
        print(' *', beta, is_p_special(beta))
    print()
    print(len(set(discovered)))


def test_q_insertion_definition(n=2):
    discovered = []
    dictionary = {}
    for w in Permutation.involutions(n):
        pipedreams = list(w.get_involution_pipe_dreams())
        for dream in pipedreams:
            e = w.number_two_cycles()
            word = dream.word()
            p, q = o_eg_insert(word)
            if p not in dictionary:
                dictionary[(w, p)] = dream.inv_monomial() * 2**e
            else:
                dictionary[(w, p)] += dream.inv_monomial() * 2**e
    for w, p in sorted(dictionary, key=lambda t: t[1].partition()):
        f = dictionary[(w, p)]
        d = decompose_into_keys(f)
        betas = sorted(d, key=lambda x: (len(x), tuple(sorted(x))))
        betas = [b for b in betas if is_special(b)]
        success = False
        for b in betas:
            if q_key(b) == f:
                # print('success:', b)
                success = True
                discovered.append(b)
        if not success:
            print('FAILED')
        else:
            continue
            print('SUCCESS')
        print()
        print('w =', w)
        print()
        print(p)
        print('  key expansion =', d)
        print()
        print('trials:')
        for beta in betas:
            print()
            print('  key expansion(', beta, ') =', decompose_into_keys(q_key(beta)))
        print()
        print()
        print()
        print()
        print()
    print('discovered:')
    for beta in sorted(discovered):
        print(' *', beta, is_special(beta))
    print()
    print(len(set(discovered)))


def is_strict(c):
    c = [i for i in c if i != 0]
    return len(set(c)) == len(c)


def is_special(c):
    c = [i for i in c if i != 0]
    if len(set(c)) < len(c):
        return False
    if any(c[j] == c[i] - 1 for i in range(len(c)) for j in range(i + 1, len(c))):
        return False
    return True


def is_p_special(c):
    if not is_special(c):
        return False
    c = [i for i in c if i in [0, 1]]
    return (1 not in c) or (c[0] != 1)


def straighten(c):
    assert is_strict(c)
    while not is_special(c):
        inv = [(i, j) for i in range(len(c)) for j in range(i + 1, len(c)) if 0 < c[j] == c[i] - 1]
        i, j = inv[0]
        c = list(c)
        c[i], c[j] = c[j], c[i]
        c = tuple(c)
    return c


def p_straighten(c):
    if 0 not in c:
        c += (0,)
    c = [(i + 1) if i != 0 else i for i in c]
    i = [i for i in range(len(c)) if c[i] == 0][0]
    c[i] += 1
    return tuple((i - 1) if i > 0 else i for i in straighten(c))


# def test_special_atoms(m=5):
#     groups = {}
#     for n in range(m + 3):
#         for k in range(m):
#             for alpha in strict_weak_compositions(n, k, reduced=True):
#                 kappa = q_atom(alpha)
#                 rev = tuple(reversed(alpha))
#                 # print(alpha, kappa != 0, is_special(rev))
#                 assert (kappa != 0) == is_special(rev)
#                 kappa = p_atom(alpha)
#                 rep = p_straighten(alpha)
#                 if rep not in groups:
#                     groups[rep] = []
#                 alpha += (len(rep) - len(alpha)) * (0,)
#                 groups[rep].append((alpha, kappa != 0))
#     for rep in groups:
#         print(rep)
#         for alpha, boolean in groups[rep]:
#             print('  ', alpha, boolean)
#         print()
#         if 1 != len([a for a, b in groups[rep] if b]):
#             input('?')


# def special_composition_count(n, k):
#     a = 0
#     for c in strict_weak_compositions(n, k):
#         if not any(0 < c[j] == c[i] - 1 for i in range(len(c)) for j in range(i + 1, len(c))):
#             a += 1
#     return a


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
            for alpha in symmetric_weak_compositions(n, k):
                kappa = p_atom(alpha)
                assert kappa.is_positive()
                assert kappa.is_not_laurent_polynomial()
    print('success')


def test_q_atom():
    for n in range(8):
        for k in range(8):
            for alpha in symmetric_weak_compositions(n, k):
                kappa = q_atom(alpha)
                assert kappa.is_positive()
                assert kappa.is_not_laurent_polynomial()
    print('success')


def icode(w):
    diagram = w.involution_rothe_diagram()
    ans = w.rank * [0]
    for i, j in diagram:
        ans[j - 1] += 1
    while ans and ans[-1] == 0:
        ans = ans[:-1]
    return tuple(ans)


def fcode(w):
    diagram = w.fpf_rothe_diagram()
    ans = w.rank * [0]
    for i, j in diagram:
        ans[j - 1] += 1
    while ans and ans[-1] == 0:
        ans = ans[:-1]
    return tuple(ans)


# def test_leading_key(m=5):
#     for n in range(m):
#         for k in range(m):
#             for alpha in weak_compositions(n, k):
#                 while alpha and alpha[-1] == 0:
#                     alpha = alpha[:-1]
#                 kappa = key(alpha)
#                 betas = get_exponents(kappa)
#                 beta = min(betas)
#                 print(alpha, beta, betas, kappa)
#                 print()
#                 assert beta == alpha


def decompose_into_custom(kappa, dictionary):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        betas = [x for x in betas if x in dictionary]
        beta = betas[0]
        hd = dict_from_tuple(beta)
        coeff = kappa[hd]
        denom = dictionary[beta][hd]
        assert coeff % denom == 0
        coeff = coeff // denom
        kappa = kappa - coeff * dictionary[beta]
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


# def try_fpf_schubert_decompose(n, dictionary):
#     i = set(Permutation.fpf_involutions(n))
#     s = {w: FPFSchubert.get(w) for w in i}
#     for w in s:
#         print(w, decompose_into_keys(s[w]))
#         w.print_fpf_rothe_diagram()
#         try:
#             dec = decompose_into_custom(s[w], dictionary)
#             print('pkeys:', dec)
#         except:
#             print('failed')
#         print()


# def try_inv_schubert_decompose(n, dictionary):
#     i = set(Permutation.involutions(n))
#     s = {w: 2 ** w.number_two_cycles() * InvSchubert.get(w) for w in i}
#     for w in s:
#         print(w, decompose_into_keys(s[w]))
#         w.print_involution_rothe_diagram(sep='.')
#         try:
#             dec = decompose_into_custom(s[w], dictionary)
#             print('qkeys:', dec)
#         except:
#             print('failed')
#         print()
#         print()
#         print()


# def print_p_key_decomposition(m=5, excluded=False, verbose=True):
#     ans = {}
#     for n in range(m):
#         for k in range(m):
#             for alpha in symmetric_weak_compositions(n, k, reduced=True):
#                 kappa = p_key(alpha)
#                 dec = decompose_into_keys(kappa)
#                 if (alpha not in dec) == excluded:
#                     if verbose:
#                         print(alpha, dec)
#                         print()
#                     ans[alpha] = kappa
#     return ans


# def test_p_leading_filter(n):
#     ans = print_p_key_decomposition(n, False, verbose=False)
#     bns = print_p_key_decomposition(n, True, verbose=False)
#     assert set(bns.values()).issubset(set(ans.values()))
#     for alpha in ans:
#         coeff = decompose_into_keys(ans[alpha])[alpha]
#         assert coeff == 1
#     assert len(set(ans.values())) == len(ans)
#     for alpha in ans:
#         print(''.join([str(s) for s in alpha]), '\\\\')


# def print_q_key_decomposition(m=5, excluded=False, verbose=True):
#     ans = {}
#     for n in range(m):
#         for k in range(m):
#             for alpha in symmetric_weak_compositions(n, k, reduced=True):
#                 kappa = q_key(alpha)
#                 dec = decompose_into_keys(kappa)
#                 if (alpha not in dec) == excluded:
#                     if verbose:
#                         print(alpha, dec)
#                         print()
#                     ans[alpha] = kappa
#     return ans


# def test_q_leading_filter(n):
#     ans = print_q_key_decomposition(n, False, verbose=False)
#     bns = print_q_key_decomposition(n, True, verbose=False)
#     assert set(bns.values()).issubset(set(ans.values()))
#     for alpha in ans:
#         coeff = decompose_into_keys(ans[alpha])[alpha]
#         assert coeff == 2**len(set(alpha) - {0})
#     assert len(set(ans.values())) == len(ans)
#     for alpha in ans:
#         print(''.join([str(s) for s in alpha]), '\\\\')


def test_p_atom_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = p_atom(alpha)
                dec = decompose_into_atoms(kappa)
                print(alpha, kappa)
                print(set(dec.values()))
                print()
                assert min({0} | set(dec.values())) >= 0


def test_q_atom_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                e = 2 ** q_power(alpha)
                kappa = q_atom(alpha)
                dec = decompose_into_atoms(kappa)
                print(alpha, kappa)
                print(set(dec.values()))
                print()
                assert all(d % e == 0 for d in dec.values())
                assert min({0} | set(dec.values())) >= 0


def test_p_key_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = p_key(alpha)
                dec = decompose_into_keys(kappa)
                ex = min({0} | set(dec.values()))
                assert ex >= 0
                assert kappa != 0


def test_q_key_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = q_key(alpha)
                dec = decompose_into_keys(kappa)
                ex = min({0} | set(dec.values()))
                assert ex >= 0
                assert kappa != 0


# def test_leading_p_key(m=5):
#     for n in range(m + 3):
#         for k in range(m):
#             valuesdict = defaultdict(list)
#             for alpha in symmetric_weak_compositions(n, k, reduced=True):
#                 kappa = p_key(alpha)
#                 valuesdict[kappa].append(alpha)
#                 betas = decompose_into_keys(kappa)
#                 coeff = betas.get(alpha, 0)
#                 betas = sorted(betas, key=lambda x: (len(x), tuple(sorted(x))))
#                 if alpha != betas[0]:
#                     print(alpha, betas, coeff)
#                     gamma = betas[0]
#                     betas = decompose_into_keys(p_key(gamma))
#                     betas = sorted(betas, key=lambda x: (len(x), tuple(sorted(x))))
#                     print(gamma, betas)
#                     print()
#                 assert alpha == betas[0]
#                 assert coeff == 1
#             assert not any(len(v) > 1 for v in valuesdict.values())


def get_markers(decomp):
    return sorted(decomp, key=lambda x: (len(x), tuple(sorted(x))))


def test_leading_q_key(m=4):
    toprint = {}
    for n in range(m + 3):
        for k in range(m):
            valuesdict = defaultdict(list)
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = q_key(alpha)
                dec = decompose_into_keys(kappa)
                valuesdict[kappa].append(alpha)
                exponents = get_exponents(kappa)
                a, b = symmetric_halves(alpha)
                beta = exponents[0]
                assert beta == b
                assert a in exponents
                toprint[tuple(exponents)] = alpha, a, b, dec
            assert not any(len(v) > 1 for v in valuesdict.values())
    prev = None
    for betas in sorted(toprint):
        alpha, a, b, dec = toprint[betas]
        if prev is None or betas[0] != prev:
            print(2 * '\n')
        prev = betas[0]
        print(alpha, ':', b, a, '->', betas)
        print()
    print('success')
    #assert False


def test_inv_schubert(n=4):
    i = set(Permutation.involutions(n))
    s = {w: InvSchubert.get(w) * 2**w.number_two_cycles() for w in i}
    print('. . . s')
    # d = {w: decompose_into_keys(s[w]) for w in i}
    d = {}
    for count, w in enumerate(i):
        d[w] = decompose_into_keys(s[w])
        print('d  ', len(i) - count)
    print('. . . d')
    k = {w: [symmetric_double(a) for a in d[w] if is_special(a)] for w in i}
    print('. . . k')
    # x = {w for w in i if q_key(k[w]) == s[w]}
    x = set()
    for count, w in enumerate(i):
        print('x  ', len(i) - count, k[w])
        if any(q_key(a) == s[w] for a in k[w]):
            x.add(w)
    vex = {w for w in i if w.is_vexillary()}
    print('. . . vex')
    print()
    assert set(x).issubset(vex)
    for w in x:
        a = symmetric_double(icode(w))
        if q_key(a) != s[w]:
            print(w)
            print(decompose_into_keys(s[w]))
            print(decompose_into_keys(q_key(a)))
            print()
        assert q_key(a) == s[w]
    vex = list(vex - x)
    for w in list(x) + vex:
        betas = get_markers(d[w])
        e = list(reversed(get_exponents(s[w])))
        try:
            alpha = symmetric_composition_from_row_column_counts(e[0], e[-1])
        except:
            alpha = 'FAILED'
        # if alpha == w.code():
        #     continue
        y, z = symmetric_halves(w.code())
        print(w)
        print(alpha == w.code(), alpha, '->', e, '->', w.code(), ':', y, z)
        # print(len(str(alpha)) * ' ', '  ', sorted(e, key=lambda x: tuple(reversed(x))))
        # print()
        # print(decompose_into_keys(s[w]))
        # print()
        print()
    return i, s, d, k, x, vex


# def get_key(f):
#     d = decompose_into_keys(f)
#     return max([a for a in d if len(a) == min(map(len, d))])


# def test_fpf_schubert(n=4):
#     def isvex(w):
#         return True
#         f = FPFStanleyExpander(w).expand()
#         if len(f) > 1:
#             return False
#         if set(f.values()) != {1}:
#             return False
#         return True
#     i = set(Permutation.fpf_involutions(n))
#     s = {w: FPFSchubert.get(w) for w in i}
#     print('. . . s')
#     # d = {w: decompose_into_keys(s[w]) for w in i}
#     d = {}
#     for count, w in enumerate(i):
#         d[w] = decompose_into_keys(s[w])
#         print('d  ', len(i) - count)
#     print('. . . d')
#     k = {w: [a for a in d[w] if is_special(a)] for w in i}
#     print('. . . k')
#     # x = {w for w in i if p_key(k[w]) == s[w]}
#     x = set()
#     for count, w in enumerate(i):
#         print('x  ', len(i) - count, k[w])
#         if any(p_key(a) == s[w] for a in k[w]):
#             x.add(w)
#     print('. . . x')
#     vex = {w for w in i if isvex(w)}
#     print('. . . vex')
#     print()
#     for w in set(x) - vex:
#         f = FPFStanleyExpander(w)
#         v = f.expand()
#         print(w)
#         print(len(v))
#         print()
#         input('???')
#     y = set()
#     for w in x:
#         a = fcode(w)
#         if p_key(a) != s[w]:
#             print(decompose_into_keys(s[w]))
#             print(a, '!=', k[w])
#             w.print_fpf_rothe_diagram(sep='.')
#             print()
#         else:
#             y.add(w)
#     return i, s, k, x, y, vex
