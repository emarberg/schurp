from keys import (
    weak_compositions,
    monomial_from_composition,
    sorting_permutation,
    symmetric_weak_compositions,
    skew_symmetric_weak_compositions,
    key, atom,
    q_power,
    p_key, p_atom,
    q_key, q_atom,
    get_exponents,
    decompose_into_compositions,
    decompose_into_keys, decompose_into_atoms, dict_from_tuple,
    symmetric_composition_from_row_column_counts, symmetric_halves,
    skew_symmetric_composition_from_row_column_counts, skew_symmetric_halves,
    maximal_decreasing_factors,
    maximal_increasing_factors,
    compatible_sequences,
    weak_compatible_sequences,
    sp_knuth_class, o_knuth_class, coxeter_knuth_class, knuth_class,
    orthogonal_key_maps, symplectic_key_maps, nil_key_maps, key_maps,
    rsk_insert, key_tableau,
    shifted_knuth_class,
    shifted_key_maps,
    is_symmetric_composition,
    symmetric_partitions,
    skew_symmetric_partitions,
)
from symmetric import FPFStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert
from permutations import Permutation
from partitions import Partition, StrictPartition
from collections import defaultdict
from words import Word
from tableaux import Tableau
import pyperclip
import pytest


q_alphas_cache = {}
q_halves_cache = {}

p_alphas_cache = {}
p_halves_cache = {}

q_insertion_cache = {}
p_insertion_cache = {}


def _attempt_p_into_q(alpha, attempts=10):
    p = p_key(alpha)
    coeff = 1
    for i in range(attempts):
        try:
            return decompose_q(p), coeff
        except:
            p *= 2
            coeff *= 2


def test_p_key_into_q_key(m=4, l=4):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
                attempt = _attempt_p_into_q(alpha)
                print(alpha, '->', attempt)


def test_q_key_into_p_key(m=4, l=4):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                attempt = try_to_decompose_p(q_key(alpha))
                print(alpha, '->', attempt)


def altschurp(mu, n):
    from polynomials import X
    p = 1
    for i in range(len(mu)):
        p *= X(i)**mu[i]
    for i in range(len(mu)):
        for j in range(i + 1, n):
            p *= 1 + X(j) * X(i)**-1
    for i in range(n):
        for j in range(n):
            p = p.isobaric_divided_difference(j)
    return p


def schurp(mu, n):
    ans = 0
    for t in Tableau.get_semistandard_shifted(Partition(*mu), n):
        ans += monomial_from_composition(t.weight())
    return ans


def schurq(mu, n):
    return schurp(mu, n) * 2**len(mu)


def test_schurp(n=3, l=3):
    for m in range(n + 1):
        for mu in StrictPartition.all(m):
            start = mu(1) + 1
            for k in range(start, start + l):
                p = schurp(mu, k)
                try:
                    alpha = decompose_p(p)
                    print(mu, k, '->', alpha)
                except:
                    print(mu, k, '->', try_to_decompose_p(p, positive=True, multiple=True))


def test_schurq(n=3, l=3):
    for m in range(n + 1):
        for mu in StrictPartition.all(m):
            start = mu(1)
            for k in range(start, start + l):
                f = schurq(mu, k)
                try:
                    alpha = decompose_q(f)
                    print(mu, k, '->', alpha)
                except:
                    print(mu, k, '->', 'FAILED')


def test_brion_construction(n=4, m=4):
    def brion_key(z, mu):
        ans = 0
        mu += (z.rank - len(mu)) * (0,)
        alphas = set()
        for w in z.get_fpf_atoms():
            w = Permutation.longest_element(z.rank) * w
            alpha = tuple(mu[w(i + 1) - 1] for i in range(len(mu)))
            alphas.add(alpha)
        for alpha in alphas:
            ans += key(alpha)
        return ans

    failures, successes = 0, 0
    for z in Permutation.fpf_involutions(n):
        for r in range(m + 1):
            for mu in Partition.all(r):
                mu = tuple(mu.parts)
                f = brion_key(z, mu)
                d = try_to_decompose_p(f)
                if not d:
                    failures += 1
                    continue
                successes += 1
                print('z =', z, 'atoms:', list(z.get_fpf_atoms()))
                print('mu =', mu)
                print('brion key =', decompose_into_keys(f))
                print()
                print('decompositions:' + (' NONE' if not d else ''))
                print()
                for dec in d:
                    print('  ', d)
                print()
                print()
                print()
                print('failures:', failures, 'vs successes', successes)
                print()


def test_atom_operators_on_keys(k=3):
    for n in range(3 * k + 1):
        for alpha in weak_compositions(n, k, reduced=False):
            print('. . .', alpha)
            for w in Permutation.all(k):
                f = key(alpha)
                word = w.get_reduced_word()
                for i in reversed(word):
                    f = f.isobaric_divided_difference(i) - f
                dec = decompose_into_atoms(f)
                assert all(v > 0 for v in dec.values())


def test_linear_dependence():
    assert len({q_key((1, 2, 3)), q_key((0, 3, 2, 1)), q_key((1, 3, 2)), q_key((0, 2, 3, 1))}) == 4
    assert q_key((1, 2, 3)) + q_key((0, 3, 2, 1)) == q_key((1, 3, 2)) + q_key((0, 2, 3, 1))

    assert len({q_key((3, 3, 4, 0, 0, 1)), q_key((3, 4, 3, 0, 0, 0, 1)), q_key((3, 4, 3, 0, 0, 1)), q_key((3, 3, 4, 0, 0, 0, 1))})
    assert q_key((3, 3, 4, 0, 0, 1)) + q_key((3, 4, 3, 0, 0, 0, 1)) == q_key((3, 4, 3, 0, 0, 1)) + q_key((3, 3, 4, 0, 0, 0, 1))


def test_symmetric_composition_from_row_column_counts():
    assert symmetric_composition_from_row_column_counts((), ()) == ()
    assert symmetric_composition_from_row_column_counts((5, 0, 0, 2, 0), (1, 1, 1, 2, 2)) == (5, 1, 1, 3, 2)
    assert symmetric_composition_from_row_column_counts((2,), (1, 1)) == (2, 1)
    assert symmetric_composition_from_row_column_counts((2, 2), (1, 2, 0, 1)) == (2, 3, 0, 1)
    assert symmetric_composition_from_row_column_counts((3, 1, 2), (1, 1, 3, 0, 0, 1)) == (3, 2, 4, 0, 0, 1)


def test_symmetric_halves(m=6, l=6):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                a, b = symmetric_halves(alpha)
                beta = symmetric_composition_from_row_column_counts(a, b)
                if alpha != beta:
                    print(alpha)
                    print(beta)
                    print()
                    print(a)
                    print(b)
                    print()
                    print()
                assert alpha == beta


def test_skew_symmetric_halves(m=6, l=6):
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


def test_schubert_compatible_sequences(n=3):
    def schubert(w):
        ans = 0
        for seq in w.get_reduced_words():
            for a, x in compatible_sequences(seq):
                ans += x
        return ans

    for w in Permutation.all(n):
        assert Schubert.get(w) == schubert(w)


def words(n, k=None):
    for w in Word.all(n, k, packed=False):
        yield w.tuple()


def test_key_compatible_sequences(m=5, l=5):
    def test_key(p):
        ans = 0
        for seq in knuth_class(p):
            rseq = list(reversed(seq))
            for a, x in compatible_sequences(rseq):
                ans += x
        return ans

    seen = set()
    keys = {}
    for n in range(m + 1):
        for k in range(l + 1):
            for w in words(n, k):
                p = rsk_insert(w)[0]
                if p in seen:
                    continue
                seen.add(p)
                d = decompose_into_keys(test_key(p))
                assert len(d) == 1 and set(d.values()) == {1}
                alpha = list(d)[0]
                keys[alpha] = keys.get(alpha, []) + [p]
    return keys


def test_nil_key_compatible_sequences(m=3):
    def test_key(p):
        ans = 0
        for seq in coxeter_knuth_class(p):
            rseq = list(reversed(seq))
            for a, x in compatible_sequences(rseq):
                ans += x
        return ans

    seen = set()
    keys = {}
    for z in Permutation.all(m):
        for w in z.get_reduced_words():
            p = eg_insert(w)[0]
            if p in seen:
                continue
            seen.add(p)
            d = decompose_into_keys(test_key(p))
            assert len(d) == 1 and set(d.values()) == {1}
            alpha = list(d)[0]
            keys[alpha] = keys.get(alpha, []) + [p]
    return keys


def test_shifted_key_compatible_sequences(m=4, l=4):
    def test_shifted_key(p):
        ans = 0
        for seq in shifted_knuth_class(p):
            assert sagan_worley_insert(seq)[0] == p
            # seq = tuple(reversed(seq))
            for a, x in weak_compatible_sequences(seq):
                ans += x
        return ans
    #
    seen = set()
    keys = {}
    functions = set()
    successes, failures = 0, 0
    for n in range(m + 1):
        for k in range(l + 1):
            for w in words(n, k):
                p = sagan_worley_insert(w)[0]
                #
                # try:
                #     if o_eg_insert(w)[0] == p:
                #         seen.add(p)
                # except:
                #     pass
                #
                if p in seen:
                    continue
                seen.add(p)
                e = len(p.partition())
                kappa = test_shifted_key(p) * 2**e
                if kappa in functions:
                    continue
                functions.add(kappa)
                dec = decompose_into_keys(kappa)
                assert all(v > 0 for v in dec.values())
                d = try_to_decompose_q(kappa, positive=True, multiple=False)
                if not d:
                    failures += 1
                    continue
                else:
                    successes += 1
                print(p)
                print('kappa  =', dec)
                for decomp in d:
                    print('      ->', decomp)
                print()
                if len(d) == 1:
                    d = list(d)[0]
                    if len(d) == 1:
                        assert list(d.values()) == [1]
                        alpha = list(d)[0]
                        keys[alpha] = keys.get(alpha, []) + [p]
                print()
                print()
                print(sorted(keys))
                print()
                print('failures:', failures, 'vs successes', successes)
                print()
    return keys


def test_p_key_compatible_sequences(m=4):
    def test_p_key(p):
        ans = 0
        for seq in sp_knuth_class(p):
            assert sp_eg_insert(seq)[0] == p
            for a, x in compatible_sequences(seq):
                ans += x
        return ans
    #
    seen = set()
    keys = set()
    for z in Permutation.fpf_involutions(m):
        for w in z.get_fpf_involution_words():
            p = sp_eg_insert(w)[0]
            if p in seen:
                continue
            kappa = test_p_key(p)
            print(p)
            print('kappa  =', decompose_into_keys(kappa))
            d = try_to_decompose_p(kappa, p_halves_cache, p_alphas_cache, positive=True, multiple=False)
            for decomp in d:
                print('      ->', decomp)
            print()
            seen.add(p)
            assert len(d) > 0
            assert all(len(dec) == 1 for dec in d)
            assert all(v == 1 for dec in d for v in dec.values())
            keys.add(list(d[0].keys())[0])
    return keys


def test_q_key_compatible_sequences(m=4):
    def test_q_key(p):
        ans = 0
        for seq in o_knuth_class(p):
            assert o_eg_insert(seq)[0] == p
            for a, x in compatible_sequences(seq):
                ans += x
        return ans
    #
    seen = set()
    keys = set()
    for z in Permutation.involutions(m):
        for w in z.get_involution_words():
            p = o_eg_insert(w)[0]
            if p in seen:
                continue
            e = len(p.partition())
            kappa = test_q_key(p) * 2**e
            d = try_to_decompose_q(kappa, q_halves_cache, q_alphas_cache, positive=True, multiple=False)
            seen.add(p)
            assert len(d) > 0
            assert all(len(dec) == 1 for dec in d)
            assert all(v == 1 for dec in d for v in dec.values())
            keys.add(list(d[0].keys())[0])
    return keys


def sagan_worley_insert(sequence):
    return Word(*sequence).sagan_worley_insert()


def inverse_sagan_worley(p, q):
    w = Tableau.inverse_sagan_worley(p, q)
    return w


def test_inverse_sagan_worley(n=5):
    for w in words(n):
        p, q = sagan_worley_insert(w)
        assert w == inverse_sagan_worley(p, q)


def test_shifted_knuth_class(n=4):
    for w in words(n):
        p, q = sagan_worley_insert(w)
        for v in shifted_knuth_class(w):
            print(w, '~', v)
            print(p)
            t = sagan_worley_insert(v)[0]
            print(t)
            assert t == p


def eg_insert(sequence):
    return Word(*sequence).eg_insert()


def o_eg_insert(sequence):
    return Word(*sequence).involution_insert()


def sp_eg_insert(sequence):
    return Word(*sequence).fpf_insert()


def test_insertion_definition(n=4):
    dictionary = {}
    for w in Permutation.all(n):
        pipedreams = list(w.get_pipe_dreams())
        for dream in pipedreams:
            word = dream.word()
            p, q = eg_insert(reversed(word))
            dictionary[p] = dictionary.get(p, 0) + dream.monomial()
    keys = {}
    for p in dictionary:
        f = dictionary[p]
        d = decompose_into_keys(f)
        alpha = list(d)[0]
        assert len(d) == 1
        assert d[alpha] == 1
        keys[alpha] = keys.get(alpha, []) + [p]
    return keys


def test_p_insertion_definition(n=4, positive=True, multiple=True):
    p_insertion_cache.clear()
    invol = list(Permutation.fpf_involutions(n))
    for i, w in enumerate(invol):
        print('. . .', len(invol) - i)
        for a in w.get_fpf_atoms():
            for dream in a.get_pipe_dreams():
                word = dream.word()
                p, q = sp_eg_insert(word)
                p_insertion_cache[p] = p_insertion_cache.get(p, 0) + dream.monomial()
    keys = {}
    for p in sorted(p_insertion_cache, key=lambda t: t.partition()):
        f = p_insertion_cache[p]
        d = try_to_decompose_p(f, p_halves_cache, p_alphas_cache, positive=positive, multiple=multiple)
        assert len(d) >= 1
        assert all(len(decomp) == 1 for decomp in d)
        assert all(set(decomp.values()) == {1} for decomp in d)
        alpha = list(list(d)[0].keys())[0]
        p_insertion_cache[p] = alpha
        if alpha not in keys:
            keys[alpha] = []
        keys[alpha] += [p]
    return keys


def test_q_insertion_definition(n=5, positive=True, multiple=True):
    q_insertion_cache.clear()
    invol = list(Permutation.involutions(n))
    for i, w in enumerate(invol):
        print('. . .', len(invol) - i)
        for a in w.get_atoms():
            for dream in a.get_pipe_dreams():
                word = dream.word()
                p, q = o_eg_insert(word)
                e = len(p.partition())
                q_insertion_cache[p] = q_insertion_cache.get(p, 0) + dream.monomial() * 2**e
    keys = {}
    for p in sorted(q_insertion_cache, key=lambda t: t.partition()):
        f = q_insertion_cache[p]
        d = try_to_decompose_q(f, q_halves_cache, q_alphas_cache, positive=positive, multiple=multiple)
        assert len(d) >= 1
        assert all(len(decomp) == 1 for decomp in d)
        alpha = list(list(d)[0].keys())[0]
        q_insertion_cache[p] = alpha
        a, b = symmetric_halves(alpha)
        print(w)
        print(p)
        print('key:', alpha, '->', a, b)
        print()
        print()
        if alpha not in keys:
            keys[alpha] = []
        keys[alpha] += [p]
    return keys


# def print_shifted_keys(n=4):
#     keys = test_shifted_key_compatible_sequences(n, n)
#     for alpha in keys:
#         x, y = symmetric_halves(alpha)
#         print('alpha =', alpha, '->', x, y)
#         print()
#         for p in keys[alpha]:
#             a, b, c, d = shifted_key_maps(p)
#             dec = decompose_into_keys(q_key(alpha))
#             print(p)
#             print('increasing keys:', a, b)
#             print('decreasing keys:', c, d)
#             print('decomposition:', dec)
#             print(b == y)
#             if b != y:
#                 input('\n?')
#             print()
#         print()
#         assert is_symmetric_composition(alpha)


def print_keys(n=4):
    results = {}
    keys = test_key_compatible_sequences(n, n)
    for alpha in keys:
        print('alpha =', alpha)
        print()
        for p in keys[alpha]:
            a, b, c, d = key_maps(p)
            print(p)
            print('increasing keys:', a, b)
            print('decreasing keys:', c, d)
            print(c.weight() == alpha)
            if c.weight() != alpha:
                input('\n?')
            print()
            results[p] = (c, d)
        print()
    print_keys_table(results, unshifted=True)


def print_nil_keys(n=4):
    results = {}
    keys = test_insertion_definition(n)
    for alpha in keys:
        print('alpha =', alpha)
        print()
        for p in keys[alpha]:
            a, b, c, d = nil_key_maps(p)
            print(p)
            print('increasing keys:', a, b)
            print('decreasing keys:', c, d)
            print(c.weight() == alpha)
            if c.weight() != alpha:
                input('\n?')
            print()
            results[p] = (c, d)
        print()
    print_keys_table(results, unshifted=True)


def _print_composition(base, alpha):
    if base is None:
        base = alpha
    s = []
    for i, q in enumerate(alpha):
        s += [q * '* ' + ((base[i] - q) if i < len(base) else 0) * '. ']
    print('\n'.join(s))


def _summarize(alpha, p, x, y, a, b, c, d, dec, last=False):
    if last:
        print('alpha =', alpha, '->', x, y)
        print()
    print(p)
    print('increasing keys:')
    print(a)
    print(a.weight())
    print()
    print(b)
    print(b.weight())
    print()
    # if c and d:
    #     print('decreasing keys:')
    #     print(c)
    #     print(c.weight())
    #     print()
    #     print(d)
    #     print(d.weight())
    #     print()
    if dec:
        print('key decomposition:', dec)
    if y is not None and b.weight() != y:
        print()
        print(b.weight(), '!=', y)
    print()
    print()
    print()
    print()


def _discrepancies(discrep):
    if discrep:
        print('DISCREPANCIES:')
        print()
        for alpha, p in discrep:
            x, y, a, b, c, d, dec = discrep[(alpha, p)]
            _summarize(alpha, p, x, y, a, b, c, d, dec, True)
        print()


def print_o_keys(n=2, positive=True, multiple=True):
    keys = test_q_insertion_definition(n, positive, multiple)
    discrep = {}
    increasing_seen, decreasing_seen = set(), set()
    results = {}
    for alpha in keys:
        x, y = symmetric_halves(alpha)
        print('alpha =', alpha, '->', x, y)
        for p in keys[alpha]:
            a, b, c, d = orthogonal_key_maps(p)
            dec = decompose_into_keys(q_key(alpha))
            _summarize(alpha, p, x, y, a, b, c, d, dec)
            if b.weight() != y:
                discrep[(alpha, p)] = (x, y, a, b, c, d, dec)
            assert b.weight() in dec
            # assert (a, b) not in increasing_seen
            # assert (c, d) not in decreasing_seen
            increasing_seen.add((a, b))
            decreasing_seen.add((c, d))
            results[p] = (a, b)
        print()
    _discrepancies(discrep)
    print_keys_table(results)


def print_sp_keys(n=2, positive=True, multiple=True):
    keys = test_p_insertion_definition(n, positive, multiple)
    discrep = {}
    increasing_seen, decreasing_seen = set(), set()
    results = {}
    for alpha in keys:
        x, y = skew_symmetric_halves(alpha)
        print('alpha =', alpha, '->', x, y)
        for p in keys[alpha]:
            a, b, c, d = symplectic_key_maps(p)
            dec = decompose_into_keys(p_key(alpha))
            _summarize(alpha, p, x, y, a, b, c, d, dec)
            if b.weight() != y:
                discrep[(alpha, p)] = (x, y, a, b, c, d, dec)
            assert b.weight() in dec
            assert dec[b.weight()] == 1
            # assert (a, b) not in increasing_seen
            # assert (c, d) not in decreasing_seen
            increasing_seen.add((a, b))
            decreasing_seen.add((c, d))
            results[p] = (a, b)
        print()
    _discrepancies(discrep)
    print_keys_table(results)


def print_shifted_keys(n=2, positive=True, multiple=True):
    increasing_seen = set()
    results = {}
    for w in words(n):
        p = sagan_worley_insert(w)[0]
        if p not in results:
            a, b = shifted_key_maps(p)
            _summarize(None, p, None, None, a, b, None, None, None)
            assert (a, b) not in increasing_seen
            increasing_seen.add((a, b))
            results[p] = (a, b)


def print_keys_table(results, unshifted=False):
    s = []
    for p in sorted(results, key=lambda x: (len(x), x.partition(), x.row_reading_word())):
        if len(p) == 0:
            continue
        left_key, right_key = results[p]
        if sorted(left_key.mapping.values()) == sorted(p.mapping.values()):
            continue
        tex = tableau_tex if unshifted else shifted_tableau_tex
        s += [tex(p), '&']
        s += [tex(left_key), '&']
        s += [tex(right_key), '\\\\ & \\\\']
    s = '\n'.join(s[:-1])
    s = """
\\begin{figure}[h]
\\begin{center}
\\begin{tabular}{llllllll}
\\begin{tabular}[t]{l|l|l}
$T$ & $K_-(T)$ & $K_+(T)$ \\\\ \\hline & \\\\
""" + s + """
\\end{tabular}
\\end{tabular}
\\end{center}
\\caption{TODO}
\\end{figure}
"""
    pyperclip.copy(s)


def tableau_tex(p, shifted=False):
    s = []
    for i, j in sorted(p.mapping):
        if len(s) < i:
            s += [(i - 1) * ['\\none']] if shifted else [[]]
        s[i - 1] += [str(p.entry(i, j))]
    s = '\\\\\n'.join(' & '.join(row) for row in reversed(s))
    return '\\begin{ytableau}\n' + s + '\n\\end{ytableau}'


def shifted_tableau_tex(p):
    return tableau_tex(p, True)


def tableau_with_keywords_tex(p, get_class, shifted):
    def increasing_label(w):
        return sorted(map(len, maximal_increasing_factors(w)))

    def decreasing_label(w):
        return sorted(map(len, maximal_decreasing_factors(w)))

    def print_group(group):
        if any(len(f[0]) == len(g[0]) and f[0] != g[0] for f in group for g in group):
            return ''
        if any(len(f[-1]) == len(g[-1]) and f[-1] != g[-1] for f in group for g in group):
            return ''
        s = ''
        for f in group:
            s += '(' + ', '.join(map(cstr, f)) + ')\\\\\n'
        s += '\\\\ \n'
        return s

    s = '$\\begin{array}{l}'
    s += tableau_tex(p, shifted)
    s += '\n\\\\ \\\\ \n'
    klass = get_class(p)

    for tagline, get_label, get_factors in [
        ('increasing:', increasing_label, maximal_increasing_factors),
        ('decreasing:', decreasing_label, maximal_decreasing_factors),
    ]:
        s += '\\emph{%s}' % tagline + ' \\\\ \\\\ \n'
        klass = sorted(klass, key=get_label)
        group = []
        for i, w in enumerate(klass):
            if i > 0 and get_label(klass[i]) != get_label(klass[i - 1]):
                s += print_group(group)
                group = []
            group += [get_factors(w)]
        s += print_group(group)
    s += '\\end{array}$'
    return s


def cstr(x):
    return ''.join([str(_) for _ in x])


def print_nil_keywords_tex(n=4):
    keys = test_insertion_definition(n)
    _print_keywords_tex(keys, coxeter_knuth_class, None, shifted=False)


def print_o_keywords_tex(n=4, positive=True, multiple=True):
    keys = test_q_insertion_definition(n, positive, multiple)
    _print_keywords_tex(keys, o_knuth_class, symmetric_halves)


def print_sp_keywords_tex(n=6, positive=True, multiple=True):
    keys = test_p_insertion_definition(n, positive, multiple)
    _print_keywords_tex(keys, sp_knuth_class, skew_symmetric_halves)


def _print_keywords_tex(keys, get_class_fn, get_halves_fn, shifted=True):
    clip = []
    for alpha in sorted(keys, key=lambda t: (sorted(t, reverse=True), len(sorting_permutation(t)))):
        if sum(alpha) == 0:
            continue
        clip += ['\\newpage\n\\begin{tabular}[t]{ll}']
        clip += ['$T$ & $\\alpha$ \\\\ \\hline \\\\']
        for i, p in enumerate(sorted(keys[alpha], key=lambda t: t.row_reading_word())):
            clip += [tableau_with_keywords_tex(p, get_class_fn, shifted)]
            clip += ['&']
            if i == 0:
                clip += ['$\\begin{array}{l}']
                clip += [cstr(alpha)]
                if get_halves_fn is not None:
                    a, b = get_halves_fn(alpha)
                    clip += ['\\\\', cstr(a), '\\\\', cstr(b)]
                clip += ['\\end{array}$']
            clip += ['\\\\ \\\\']
        clip += ['\\end{tabular}', '\\newpage']
    clip = '\n'.join(clip)
    pyperclip.copy(clip)


def print_shifted_table_tex(keys, halves_fn):
    s = []
    for alpha in sorted(keys, key=lambda t: (sorted(t, reverse=True), len(sorting_permutation(t)))):
        if sum(alpha) == 0:
            continue
        a, b = halves_fn(alpha)
        for i, p in enumerate(sorted(keys[alpha], key=lambda t: t.row_reading_word())):
            s += [shifted_tableau_tex(p)]
            s += ['&']
            if i == 0:
                s += ['$\\begin{array}{l}']
                s += [cstr(alpha), '\\\\']
                s += [cstr(a), '\\\\']
                s += [cstr(b), '\\end{array}$']
            s += ['\\\\ \\\\']
    s = s[:-1]
    s = """
\\begin{figure}[h]
\\begin{center}
\\begin{tabular}{llllllll}
\\begin{tabular}[t]{ll}
$T$ & $\\alpha$ \\\\ \\hline \\\\
""" + '\n'.join(s) + """
\\end{tabular}
\\end{tabular}
\\end{center}
\\caption{TODO}
\\end{figure}
"""
    pyperclip.copy(s)


def print_o_table_tex(n=4, positive=True, multiple=True):
    keys = test_q_insertion_definition(n, positive, multiple)
    print_shifted_table_tex(keys, symmetric_halves)


def print_sp_table_tex(n=4, positive=True, multiple=True):
    keys = test_p_insertion_definition(n, positive, multiple)
    print_shifted_table_tex(keys, skew_symmetric_halves)


def print_table_tex(n=4):
    keys = test_insertion_definition(n)

    s = []
    for alpha in sorted(keys, key=lambda t: (len(keys[t][0]), keys[t][0].row_reading_word())):
        if len(keys[alpha]) <= 1:
            continue
        for i, p in enumerate(keys[alpha]):
            s += [tableau_tex(p)]
            s += ['&']
            if i == 0:
                s += [tableau_tex(key_tableau(alpha))]
            s += ['\\\\ \\\\']
    s = s[:-1]
    s = """
\\begin{figure}[h]
\\begin{center}
\\begin{tabular}{lcr}
\\begin{tabular}[t]{ll}
$T$ & $K_\\alpha$ \\\\ \\hline \\\\
""" + '\n'.join(s) + """
\\end{tabular}
\\end{tabular}
\\end{center}
\\caption{TODO}
\\end{figure}
"""
    pyperclip.copy(s)


def schur(partition):
    assert all(partition[i] >= partition[i + 1] for i in range(len(partition) - 1))
    w = Permutation.get_grassmannian(*partition)
    n = len(partition)
    return Schubert.get(w).truncate(n)


def test_schur():
    assert schur((3, 1)) == key((1, 3, 0, 0))


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
        monomial_from_composition((2, 1, 1, 0)) + \
        monomial_from_composition((1, 2, 1, 0)) + \
        monomial_from_composition((1, 1, 2, 0)) + \
        monomial_from_composition((2, 1, 0, 1)) + \
        monomial_from_composition((2, 0, 1, 1)) + \
        monomial_from_composition((1, 2, 0, 1)) + \
        monomial_from_composition((1, 1, 1, 1)) + \
        monomial_from_composition((1, 0, 2, 1))
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


def test_p_key_start_decomposition(m=5):
    # FAILS
    for n in range(m + 3):
        for alpha in symmetric_partitions(n):
            kappa = p_key(alpha)
            dec = decompose_into_keys(kappa)
            ex = min({0} | set(dec.values()))
            assert ex >= 0
            assert kappa != 0
            if set(dec.values()) != {1}:
                a, b = skew_symmetric_halves(alpha)
                print('alpha =', alpha, '->', a, b)
                print('decomposition =', dec)
                print('decomposition =', dec.values())
            assert set(dec.values()) == {1}


def test_q_key_start_decomposition(m=5):
    # FAILS
    for n in range(m + 1):
        for alpha in symmetric_partitions(n):
            kappa = q_key(alpha)
            dec = decompose_into_keys(kappa)
            ex = min({0} | set(dec.values()))
            assert ex >= 0
            assert kappa != 0
            if set(dec.values()) != {2**q_power(alpha)}:
                a, b = symmetric_halves(alpha)
                print('alpha =', alpha, '->', a, b)
                print('decomposition =', dec)
                print('decomposition =', dec.values())
            assert set(dec.values()) == {2**q_power(alpha)}


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
    for n in range(m + 1):
        for k in range(m + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = q_key(alpha)
                dec = decompose_into_keys(kappa)
                ex = min({0} | set(dec.values()))
                assert ex >= 0
                assert kappa != 0


# def _naive_next(w, base):
#     ell = len(base)
#     values = sorted(set(w) - {0})
#     y = 0
#     for v in values:
#         inds = [i for i in range(len(w)) if w[i] == v]
#         i, inds = inds[0], inds[1:]
#         if i < ell:
#             u = list(w)
#             if inds and any(u[k] == 0 for k in range(i + 1, inds[0])):
#                 continue
#             for j in inds:
#                 u[i], u[j] = u[i] + 1, u[j] - 1

#                 if u[i] > base[i]:
#                     break
#                 if j < ell and u[j] < base[j]:
#                     break

#                 tu = tuple(u)
#                 while tu and tu[-1] == 0:
#                     tu = tu[:-1]
#                 yield tu
#                 y += 1
#     if y == 0:
#         pass


# def _naive_next(w):
#     values = sorted(a[-1] for a in w if a)
#     for v in values:
#         inds = [i for i in range(len(w)) if w[i] and w[i][-1] == v]
#         u = [list(a) for a in w]
#         i = inds[0]
#         inds = inds[1:]
#         while inds:
#             j, inds = inds[0], inds[1:]
#             u[i], u[j] = u[i] + [u[j][-1] + 1], u[j][:-1]
#             for k in inds:
#                 u[k] = u[k][:-1] + [u[k][-1] + 1]
#             yield tuple(tuple(a) for a in u)


# def _naive_comp(w):
#     comp = tuple(len(a) for a in w)
#     while comp and comp[-1] == 0:
#         comp = comp[:-1]
#     return comp


# def _naive_string(w):
#     s = []
#     for a in w:
#         s += [((1 + max(a)) if a else 0) * ['  ']]
#         for i in a:
#             s[-1][i] = '* '
#         s[-1] = ''.join(s[-1])
#     s = '\n'.join(s)
#     return s


def _naive_next(w, base):
    for j in range(len(w)):
        for i in range(j + 1, len(w)):
            u = [list(a) for a in w]
            if w[i][j] == '*' and i == j + 1:
                k = i
                while w[j][k] == '*':
                    k += 1
                u[i][j], u[j][k] = '.', '*'
                yield tuple(tuple(a) for a in u)
            elif w[i][j] == '*' and i > j + 1 and w[i - 1][j] != '*':
                u[i][j], u[i - 1][j] = '.', '*'
                #u[i][j], u[j][i] = '.', '*'
                #u[i - 1][j], u[j][i - 1] = '*', '.'
                yield tuple(tuple(a) for a in u)


def _naive_comp(w):
    comp = tuple(len([i for i in a if i == '*']) for a in w)
    while comp and comp[-1] == 0:
        comp = comp[:-1]
    return comp


def _naive_string(w):
    def row(a):
        return ('',)
        # return (str(len([i for i in a if i == '*'])),)

    def col(w):
        ans = len(w) * [0]
        for i in range(len(w)):
            for j in range(len(w)):
                ans[j] += w[i][j] == '*'
        # return '\n' + ' '.join(map(str, ans))
        return ''

    # n = len(w)
    # v = [n * [' '] for _ in range(n)]
    # for i in range(n):
    #     for j in range(n):
    #         if w[i][j] == '*':
    #             if i == j:
    #                 v[i][0] = '*'
    #             elif i < j:
    #                 v[i][j] = '*'
    #             else:
    #                 v[i][j + 1] = '*'
    # w = tuple(tuple(a) for a in v)
    return '\n'.join([' '.join(a + row(a)) for a in w]) + col(w)


def p_naive_generator(alpha, base):
    _, x = symmetric_halves(alpha)
    n = max([0] + list(alpha))
    # start = tuple(tuple(j for j in range(a)) for a in x)
    start = [n * ['.'] for a in range(n)]
    for i in range(len(x)):
        for j in range(x[i]):
            if i != j:
                start[i][j] = '*'
    start = tuple(tuple(a) for a in start)
    seen = set()
    add = {start}
    while add:
        next_to_add = set()
        for w in add:
            yield w, _naive_comp(w), _naive_string(w)
            seen.add(w)
            for v in _naive_next(w, base):
                if v not in seen:
                    next_to_add.add(v)
        add = next_to_add


def naive_generator(alpha, base):
    _, x = symmetric_halves(alpha)
    n = max([0] + list(alpha))
    # start = tuple(tuple(j for j in range(a)) for a in x)
    start = [n * ['.'] for a in range(n)]
    for i in range(len(x)):
        for j in range(x[i]):
            start[i][j] = '*'
    start = tuple(tuple(a) for a in start)
    seen = set()
    add = {start}
    while add:
        next_to_add = set()
        for w in add:
            yield w, _naive_comp(w), _naive_string(w)
            seen.add(w)
            for v in _naive_next(w, base):
                if v not in seen:
                    next_to_add.add(v)
        add = next_to_add


def test_p_partition_key_expansion(m=5):
    success, failure = 0, 0
    for n in range(m + 1):
        for alpha in skew_symmetric_partitions(n):
            kappa = p_key(alpha)
            dec = decompose_into_keys(kappa)
            ex = min({0} | set(dec.values()))
            assert ex >= 0
            assert kappa != 0
            a, b = symmetric_halves(alpha)
            dec, val = set(dec), set(dec.values())
            trips = set(p_naive_generator(alpha, max(dec)))
            print('alpha =', alpha, '->', a, b)
            print()
            print('ACTUAL KEY DECOMPOSITION:', dec, val)
            # print('monomial decomposition:', dec)
            # for z in sorted(dec):
            #     print('\ncomposition =', z)
            #     print()
            #     _print_composition(max(dec), z)
            #     print()
            print()
            print()
            print('NAIVE GUESS:')
            dec = decompose_into_compositions(kappa)
            naive = {}
            for w, c, s in sorted(trips, key=lambda i: i[1]):
                naive[c] = naive.get(c, 0) + 1
                print('\ncomposition =', c)
                print()
                ss = [line + '  ->' for line in s.split('\n')]
                for v in _naive_next(w, None):
                    for i, line in enumerate(_naive_string(v).split('\n')):
                        ss[i] += '   ' + line
                s = '\n'.join(ss)
                print(s)
                print()
                # print('NAIVE GUESS:')
                # for z in sorted(naive):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(max(dec), z)
                #     print()
                # for _, c, s in sorted(trips, key=lambda t: t[1]):
                #     print('\ncomposition =', c)
                #     print()
                #     print(s)
                #     print()
                # print('WORKS?')
                # print()
                # print(dec == naive, dec.issubset(naive), naive.issubset(dec))
                # for z in sorted(naive - dec):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(a, z)
                #     print()
                # print()
                # print('---------')
                # print()
                # for z in sorted(dec - naive):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(max(dec), z)
                #     print()
                # print()
                # print()
                # print()
            success += dec == naive
            failure += dec != naive
            if dec != naive:
                print('dec =', dec)
                print('nai =', naive)
                print()
                input('')
            # assert naive.issubset(dec)
    print()
    print('success:', success)
    print('failure', failure)
    print()


def test_q_partition_key_expansion(m=5):
    success, failure = 0, 0
    for n in range(m + 1):
        for alpha in symmetric_partitions(n):
            kappa = q_key(alpha)
            dec = decompose_into_keys(kappa)
            ex = min({0} | set(dec.values()))
            assert ex >= 0
            assert kappa != 0
            a, b = symmetric_halves(alpha)
            dec, val = set(dec), set(dec.values())
            trips = set(naive_generator(alpha, max(dec)))
            print('alpha =', alpha, '->', a, b)
            print()
            print('ACTUAL KEY DECOMPOSITION:', sorted(dec), val)
            # print('monomial decomposition:', dec)
            # for z in sorted(dec):
            #     print('\ncomposition =', z)
            #     print()
            #     _print_composition(max(dec), z)
            #     print()
            print()
            print()
            print('NAIVE GUESS:')
            mon = decompose_into_compositions(kappa)
            naive = {}
            for w, c, s in sorted(trips, key=lambda i: i[1]):
                naive[c] = naive.get(c, 0) + 2 ** len([i for i in range(len(w)) if w[i][i] == '*'])
                if c not in dec:
                    continue
                print('\ncomposition =', c)
                print()
                ss = [line + '  ->' for line in s.split('\n')]
                for v in _naive_next(w, None):
                    for i, line in enumerate(_naive_string(v).split('\n')):
                        ss[i] += '   ' + line
                s = '\n'.join(ss)
                print(s)
                print()
                # print('NAIVE GUESS:')
                # for z in sorted(naive):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(max(dec), z)
                #     print()
                # for _, c, s in sorted(trips, key=lambda t: t[1]):
                #     print('\ncomposition =', c)
                #     print()
                #     print(s)
                #     print()
                # print('WORKS?')
                # print()
                # print(dec == naive, dec.issubset(naive), naive.issubset(dec))
                # for z in sorted(naive - dec):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(a, z)
                #     print()
                # print()
                # print('---------')
                # print()
                # for z in sorted(dec - naive):
                #     print('\ncomposition =', z)
                #     print()
                #     _print_composition(max(dec), z)
                #     print()
                # print()
                # print()
                # print()
            success += mon == naive
            failure += mon != naive
            if mon != naive:
                print('dec =', dec)
                print('nai =', naive)
                print()
                input('')
            # assert naive.issubset(dec)
    print()
    print('success:', success)
    print('failure', failure)
    print()


def is_power_of_two(x):
    assert type(x) == int
    if x == 1:
        return True
    if x <= 0 or x % 2 != 0:
        return False
    return is_power_of_two(x // 2)


@pytest.mark.slow
def test_leading_p_key(m=30, l=4):
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
                assert kappa[dict_from_tuple(a)] >= 1
                assert kappa[dict_from_tuple(b)] == 1
                toprint[tuple(exponents)] = alpha, a, b, dec
    if any(len(v) > 1 for v in valuesdict.values()):
        for kappa in sorted(valuesdict, key=lambda k: (len(valuesdict[k]), get_exponents(k)[0])):
            alphas = valuesdict[kappa]
            print(get_exponents(kappa)[0], '. . .')
            for a in alphas:
                b, c = skew_symmetric_halves(a)
                print('  ', a, '->', c, b)
            print()
    print()
    print('values:', len(valuesdict))
    assert not any(len(v) > 1 for v in valuesdict.values())
    prev = None
    for betas in sorted(toprint):
        alpha, a, b, dec = toprint[betas]
        if prev is None or betas[0] != prev:
            print(2 * '\n')
        prev = betas[0]
        print(alpha, ':', b, a, '->', betas)
        print()


def test_leading_q_key(m=4, l=4):
    toprint = {}
    valuesdict = defaultdict(list)
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                kappa = q_key(alpha)
                dec = decompose_into_keys(kappa)
                valuesdict[kappa].append(alpha)
                exponents = get_exponents(kappa)
                a, b = symmetric_halves(alpha)
                beta = exponents[0]
                assert beta == b
                assert a in exponents
                assert kappa[dict_from_tuple(a)] >= 2**q_power(alpha)
                assert kappa[dict_from_tuple(b)] == 2**q_power(alpha)
                toprint[tuple(exponents)] = alpha, a, b, dec
    prev = None
    for betas in sorted(toprint):
        alpha, a, b, dec = toprint[betas]
        if prev is None or betas[0] != prev:
            print(2 * '\n')
        prev = betas[0]
        print(alpha, ':', b, a, '->', betas)
        print()
    assert not any(len(v) > 1 for v in valuesdict.values())


def q_update(targets, exponents, halves, alphas):
    for e in targets:
        for d in exponents:
            try:
                a = symmetric_composition_from_row_column_counts(d, e)
                assert symmetric_halves(a) == (d, e)
            except:
                continue
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


def decompose_p(f):
    answers = try_to_decompose_p(f, positive=True)
    assert len(answers) == 1
    ans = answers[0]
    assert len(ans) == 1
    assert set(ans.values()) == {1}
    return list(ans.keys())[0]

def decompose_q(f):
    answers = try_to_decompose_q(f, positive=True, multiple=True)
    assert len(answers) == 1
    ans = answers[0]
    assert len(ans) == 1
    assert set(ans.values()) == {1}
    return list(ans.keys())[0]


def try_to_decompose_q(f, halves=None, alphas=None, positive=True, multiple=False):
    if halves is None:
        halves = q_halves_cache
    if alphas is None:
        alphas = q_alphas_cache
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
        w.print_rothe_diagram(sep='.')
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


def try_to_decompose_p(f, halves=None, alphas=None, positive=True, multiple=False):
    if halves is None:
        halves = p_halves_cache
    if alphas is None:
        alphas = p_alphas_cache
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
