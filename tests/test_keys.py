from keys import (
    weak_compositions,
    monomial,
    sorting_permutation,
    symmetric_weak_compositions,
    skew_symmetric_weak_compositions,
    key, atom,
    q_power,
    p_key, p_atom,
    q_key, q_atom,
    has_distinct_parts,
    get_exponents,
    decompose_into_keys, decompose_into_atoms, dict_from_tuple,
    symmetric_composition_from_row_column_counts, symmetric_halves,
    skew_symmetric_composition_from_row_column_counts, skew_symmetric_halves
)
from symmetric import FPFStanleyExpander, InvStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert, X
from permutations import Permutation
from collections import defaultdict
from words import Word
from partitions import Partition, StrictPartition
from tableaux import Tableau


q_alphas_cache = {}
q_halves_cache = {}

p_alphas_cache = {}
p_halves_cache = {}

q_insertion_cache = {}
p_insertion_cache = {}


def test_linear_dependence():
    assert q_key((1, 2, 3)) + q_key((0, 3, 2, 1)) == q_key((1, 3, 2)) + q_key((0, 2, 3, 1))
    assert q_key((3, 3, 4, 0, 0, 1)) + q_key((3, 4, 3, 0, 0, 0, 1)) == q_key((3, 4, 3, 0, 0, 1)) + q_key((3, 3, 4, 0, 0, 0, 1))


def test_symmetric_composition_from_row_column_counts(m=5):
    assert symmetric_composition_from_row_column_counts((), ()) == ()
    assert symmetric_composition_from_row_column_counts((5, 0, 0, 2, 0), (1, 1, 1, 2, 2)) == (5, 1, 1, 3, 2)
    assert symmetric_composition_from_row_column_counts((2,), (1, 1)) == (2, 1)
    assert symmetric_composition_from_row_column_counts((2, 2), (1, 2, 0, 1)) == (2, 3, 0, 1)
    assert symmetric_composition_from_row_column_counts((3, 1, 2), (1, 1, 3, 0, 0, 1)) == (3, 2, 4, 0, 0, 1)


def test_symmetric_halves(m=10, l=6):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                a, b = symmetric_halves(alpha)
                beta = symmetric_composition_from_row_column_counts(a, b)
                print(alpha)
                print(beta)
                print(a, b)
                print()
                assert alpha == beta


def test_skew_symmetric_halves(m=10, l=6):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
                a, b = skew_symmetric_halves(alpha)
                beta = skew_symmetric_composition_from_row_column_counts(a, b)
                print(alpha)
                print(beta)
                print(a, b)
                print()
                assert alpha == beta


# def partitions(n):
#     for mu in Partition.all(n):
#         yield tuple(mu.parts)


# def strict_partitions(n):
#     for mu in StrictPartition.all(n):
#         yield tuple(mu.parts)


# def symmetric_compositions(n):
#     for mu in strict_partitions(n):
#         m = mu[0] if mu else 0
#         shape = {(i, j) for i in range(1, len(mu) + 1) for j in range(i, i + mu[i - 1])}
#         shape |= {(j, i) for (i, j) in shape}
#         shape = tuple(sorted(shape))
#         #
#         seen = set()
#         toadd = {shape}
#         while toadd:
#             next_toadd = set()
#             for mu in toadd:
#                 yield Tableau({b: 1 for b in mu})
#                 seen.add(mu)
#                 for i in range(1, m):
#                     s = Permutation.s_i(i)
#                     nu = tuple(sorted((s(a), s(b)) for (a, b) in mu))
#                     if nu not in seen:
#                         next_toadd.add(nu)
#             toadd = next_toadd


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


def test_p_insertion_definition(n=2, positive=True, multiple=True):
    p_insertion_cache.clear()
    for w in Permutation.fpf_involutions(n):
        # for dream in w.get_fpf_involution_pipe_dreams():
        #     word = dream.word()
        #     p, q = sp_eg_insert(word)
        #     p_insertion_cache[(w, p)] = p_insertion_cache.get((w, p), 0) + dream.fpf_monomial()
        for a in w.get_fpf_atoms():
            for dream in a.get_pipe_dreams():
                word = dream.word()
                p, q = sp_eg_insert(word)
                p_insertion_cache[(w, p)] = p_insertion_cache.get((w, p), 0) + dream.monomial()
    for w, p in sorted(p_insertion_cache, key=lambda t: t[1].partition()):
        f = p_insertion_cache[(w, p)]
        d = try_to_decompose_p(f, p_halves_cache, p_alphas_cache, positive=positive, multiple=multiple)
        for decomp in d:
            print(w, '->', decomp)
        print()
        print(p)
        print()
        if len(d) == 0:
            print('keys:', decompose_into_keys(f), '\n***\n')
        assert len(d) >= 1
        assert all(len(decomp) == 1 for decomp in d)
        assert all(set(decomp.values()) == {1} for decomp in d)


def test_q_insertion_definition(n=2, positive=True, multiple=True):
    q_insertion_cache.clear()
    for w in Permutation.involutions(n):
        e = w.number_two_cycles()
        # for dream in w.get_involution_pipe_dreams():
        #     word = dream.word()
        #     p, q = o_eg_insert(word)
        #     q_insertion_cache[(w, p)] = q_insertion_cache.get((w, p), 0) + dream.inv_monomial() * 2**e
        for a in w.get_atoms():
            for dream in a.get_pipe_dreams():
                word = dream.word()
                p, q = o_eg_insert(word)
                q_insertion_cache[(w, p)] = q_insertion_cache.get((w, p), 0) + dream.monomial() * 2**e
    for w, p in sorted(q_insertion_cache, key=lambda t: t[1].partition()):
        f = q_insertion_cache[(w, p)]
        d = try_to_decompose_q(f, q_halves_cache, q_alphas_cache, positive=positive, multiple=multiple)
        for decomp in d:
            print(w, '->', decomp)
        print()
        print(p)
        print()
        w.print_essential_set()
        print()
        mu = sorted(list(d[0].keys())[0], reverse=True)
        print('half lambda =', symmetric_halves(mu)[0])
        print()
        if len(d) == 0:
            print('keys:', decompose_into_keys(f), '\n***\n')
        assert len(d) >= 1
        assert all(len(decomp) == 1 for decomp in d)


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
            for alpha in skew_symmetric_weak_compositions(n, k):
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


def test_leading_key(m=5):
    for n in range(m):
        for k in range(m):
            for alpha in weak_compositions(n, k, reduced=True):
                kappa = key(alpha)
                betas = get_exponents(kappa)
                beta = min(betas)
                print(alpha, beta, betas, kappa)
                print()
                assert beta == alpha


def test_p_atom_decomposition(m=5):
    for n in range(m + 3):
        for k in range(m):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
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
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
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


def test_leading_p_key(m=8, l=5):
    toprint = {}
    valuesdict = defaultdict(list)
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
                kappa = p_key(alpha)
                dec = decompose_into_keys(kappa)
                valuesdict[kappa].append(alpha)
                exponents = get_exponents(kappa)
                a, b = skew_symmetric_halves(alpha)
                beta = exponents[0]
                assert beta == b
                assert a in exponents
                toprint[tuple(exponents)] = alpha, a, b, dec
    if any(len(v) > 1 for v in valuesdict.values()):
        for kappa in sorted(valuesdict, key=lambda kappa: get_exponents(kappa)[0]):
            alphas = valuesdict[kappa]
            print(get_exponents(kappa)[0], '. . .')
            for a in alphas:
                print('  ', a)
            print()
    assert not any(len(v) > 1 for v in valuesdict.values())
    prev = None
    for betas in sorted(toprint):
        alpha, a, b, dec = toprint[betas]
        if prev is None or betas[0] != prev:
            print(2 * '\n')
        prev = betas[0]
        print(alpha, ':', b, a, '->', betas)
        print()


def test_leading_q_key(m=4):
    toprint = {}
    valuesdict = defaultdict(list)
    for n in range(m + 3):
        for k in range(m):
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


def q_update(targets, exponents, halves, alphas):
    for e in targets:
        for d in exponents:
            try:
                a = symmetric_composition_from_row_column_counts(d, e)
            except:
                continue
            assert symmetric_halves(a) == (d, e)
            if a not in alphas:
                alphas[a] = q_key(a)
                halves[e] = halves.get(e, []) + [(alphas[a], a)]


def test_decompose_q():
    f = q_key((0, 1))
    assert try_to_decompose_q(f, {}, {}, False, True) == [
        {(0, 1): 1}
    ]

    f = q_key((0, 1)) + 2 * q_key((1, 4, 2, 3))
    assert try_to_decompose_q(f, {}, {}, True, True) == [
        {(0, 1): 1, (1, 4, 2, 3): 2}
    ]


def try_to_decompose_q(f, halves={}, alphas={}, positive=True, multiple=False):
    if f == 0:
        return [{}]
    if positive and not f.is_positive():
        return []
    exponents = get_exponents(f)
    targets = [exponents[0]]
    q_update(targets, exponents, halves, alphas)
    answers = []
    for target in targets:
        dict_key = dict_from_tuple(target)
        for g, alpha in sorted(halves.get(target, [])):
            assert g == q_key(alpha)
            a = f[dict_key]
            b = g[dict_key]
            if a % b == 0:
                h = f - a // b * g
                assert h[dict_key] == 0
                for ans in try_to_decompose_q(h, halves, alphas, positive, multiple):
                    ans[alpha] = ans.get(alpha, 0) + a // b
                    if ans[alpha] == 0:
                        del ans[alpha]
                    if not multiple:
                        return [ans]
                    elif ans not in answers:
                        answers.append(ans)
    for ans in answers:
        g = 0
        for alpha, coeff in ans.items():
            g += coeff * q_key(alpha)
        assert f == g
    return answers


def test_inv_schubert(n=4, positive=True, multiple=True):
    i = list(Permutation.involutions(n))
    s = {w: InvSchubert.get(w) * 2**w.number_two_cycles() for w in i}
    print('. . . s')
    d = {}
    for t, w in enumerate(s):
        d[w] = try_to_decompose_q(s[w], q_halves_cache, q_alphas_cache, positive, multiple)
        for dec in d[w]:
            print(len(s) - t, ':', w, '->', dec)
        print()
        assert (not positive and multiple) or len(d[w]) == 1
        d[w] = d[w][0]
    qvex = {w: list(d[w])[0] for w in d if len(d[w]) == 1 and set(d[w].values()) == {1}}
    ivex = {w: w.code() for w in i if w.is_vexillary()}
    assert qvex == ivex
    print()
    print('caches =', len(q_halves_cache), len(q_alphas_cache))
    return i, s, d, qvex, ivex


def p_update(targets, exponents, halves, alphas):
    for e in targets:
        for d in exponents:
            try:
                a = skew_symmetric_composition_from_row_column_counts(d, e)
            except:
                continue
            assert skew_symmetric_halves(a) == (d, e)
            if a not in alphas:
                alphas[a] = p_key(a)
                halves[e] = halves.get(e, []) + [(alphas[a], a)]


def try_to_decompose_p(f, halves, alphas, positive=True, multiple=False):
    if f == 0:
        return [{}]
    if positive and not f.is_positive():
        return []
    exponents = get_exponents(f)
    targets = [exponents[0]]
    p_update(targets, exponents, halves, alphas)
    answers = []
    for target in targets:
        dict_key = dict_from_tuple(target)
        for g, alpha in sorted(halves.get(target, [])):
            a = f[dict_key]
            b = g[dict_key]
            if a % b == 0:
                for ans in try_to_decompose_p(f - a // b * g, halves, alphas, positive, multiple):
                    ans[alpha] = ans.get(alpha, 0) + a // b
                    if ans[alpha] == 0:
                        del ans[alpha]
                    if not multiple:
                        return [ans]
                    elif ans not in answers:
                        answers.append(ans)
    for ans in answers:
        g = 0
        for alpha, coeff in ans.items():
            g += coeff * p_key(alpha)
        assert f == g
    return answers


def is_fpf_vexillary(w):
    return True
    f = FPFStanleyExpander(w).expand()
    if len(f) > 1:
        return False
    if set(f.values()) != {1}:
        return False
    return True


def test_fpf_schubert(n=4, positive=True, multiple=True):
    i = list(Permutation.fpf_involutions(n))
    s = {w: FPFSchubert.get(w) for w in i}
    print('. . . s')
    d = {}
    for t, w in enumerate(s):
        d[w] = try_to_decompose_p(s[w], p_halves_cache, p_alphas_cache, positive, multiple)
        for dec in d[w]:
            print(len(s) - t, ':', w, '->', dec)
            assert set(dec.values()) == {1}
        print()
        # assert (not positive and multiple) or len(d[w]) == 1
        d[w] = sorted(d[w], key=lambda x: (len(x), sorted(x.values())))[0]
    pvex = {w: list(d[w])[0] for w in d if len(d[w]) == 1 and set(d[w].values()) == {1}}
    fvex = {w: w.code() for w in i if is_fpf_vexillary(w)}
    # assert pvex == fvex
    print()
    print('caches =', len(p_halves_cache), len(p_alphas_cache))
    return i, s, d, pvex, fvex
