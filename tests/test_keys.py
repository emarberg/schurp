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
    has_distinct_parts,
    get_exponents,
    decompose_into_keys, decompose_into_atoms, dict_from_tuple,
    symmetric_composition_from_row_column_counts, symmetric_halves,
    skew_symmetric_composition_from_row_column_counts, skew_symmetric_halves,
    symmetric_double
)
from symmetric import FPFStanleyExpander, InvStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert, X
from permutations import Permutation
from collections import defaultdict
from words import Word
from partitions import Partition, StrictPartition
from tableaux import Tableau
import pyperclip


q_alphas_cache = {}
q_halves_cache = {}

p_alphas_cache = {}
p_halves_cache = {}

q_insertion_cache = {}
p_insertion_cache = {}

o_reduced_tableau_cache = {}


def maximal_increasing_factors(w):
    factors = [[]]
    for a in w:
        if len(factors[-1]) == 0 or factors[-1][-1] < a:
            factors[-1].append(a)
        else:
            factors.append([a])
    return tuple(tuple(a) for a in factors)


def maximal_decreasing_factors(w):
    factors = [[]]
    for a in w:
        if len(factors[-1]) == 0 or factors[-1][-1] > a:
            factors[-1].append(a)
        else:
            factors.append([a])
    return tuple(tuple(a) for a in factors)


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


def colform(w):
    return Partition(*sorted([len(a) for a in maximal_decreasing_factors(w)], reverse=True)).tuple()


def is_key_word(p, w):
    return p == w.rsk_insert(w)[0] and colform(w) == p.partition().tuple()


def key_words(p):
    for w in knuth_class(p):
        if colform(w) == p.partition().transpose().tuple():
            yield w


def key_tableau(alpha):
    w = Permutation()
    for i in sorting_permutation(alpha):
        w *= Permutation.s_i(i)
    word = []
    for part in Partition(*sorted(alpha, reverse=True)).transpose().parts:
        word += sorted([w(i) for i in range(part, 0, -1)], reverse=True)
    return rsk_insert(word)[0]


# def symmetric_key_tableau(alpha):
#     mu = Partition(*sorted(alpha, reverse=True))
#     assert mu.is_symmetric()
#     w = Permutation()
#     for i in sorting_permutation(alpha):
#         w *= Permutation.s_i(i)
#     word = []
#     for part in mu.parts:
#         word += sorted([w(i) for i in range(part, 0, -1)], reverse=True)
#     return sagan_worley_insert(word)[0]


def get_keys(p):
    left_columns = {}
    right_columns = {}
    for w in key_words(p):
        f = maximal_decreasing_factors(w)
        a, b = f[0], f[-1]
        if len(a) in left_columns:
            assert left_columns[len(a)] == a
        if len(b) in right_columns:
            assert right_columns[len(b)] == b
        left_columns[len(a)] = a
        right_columns[len(b)] = b
    right_word = ()
    left_word = ()
    for part in p.partition().transpose().tuple():
        right_word += right_columns[part]
        left_word += left_columns[part]
    return rsk_insert(left_word)[0], rsk_insert(right_word)[0]


def compatible_sequences(seq, i_min=1):
    if len(seq) == 0:
        yield (), X(0)**0
    else:
        a, seq = seq[0], seq[1:]
        for i in range(i_min, a + 1):
            j_min = (i + 1) if (seq and a < seq[0]) else i
            for p, q in compatible_sequences(seq, j_min):
                yield (a,) + p, X(i) * q


def test_schubert_compatible_sequences(n=5):
    def schubert(w):
        ans = 0
        for seq in w.get_reduced_words():
            for a, x in compatible_sequences(seq):
                ans += x
        return ans

    for w in Permutation.all(n):
        assert Schubert.get(w) == schubert(w)


def test_key_compatible_sequences(m=5):
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
            d = decompose_into_keys(test_key(p))
            assert len(d) == 1 and set(d.values()) == {1}
            alpha = list(d)[0]
            seen.add(p)
            keys[alpha] = keys.get(alpha, []) + [p]
    print_tex(keys)


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
    keys = {}
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
            alpha = list(list(d)[0].keys())[0]
            if alpha not in keys:
                keys[alpha] = []
            keys[alpha] += [p]
    #
    for alpha in keys:
        print(alpha, '->', skew_symmetric_halves(alpha))
        for p in keys[alpha]:
            print(p)
        print()
        print()
    print('keys:', len(keys))


def tableau_tex(p):
    s = []
    for i, j in sorted(p.mapping):
        if len(s) < i:
            s += [[]]
        s[i - 1] += [str(p.entry(i, j))]
    s = '\\\\\n'.join(' & '.join(row) for row in reversed(s))
    return '\\begin{ytableau}\n' + s + '\n\\end{ytableau}'


def o_shifted_keys(p):
    nu, _ = symmetric_halves(symmetric_double(p.partition().tuple()))
    nu = tuple(sorted(nu, reverse=True))

    increasing_left_keys, increasing_right_keys = {}, {}
    for w in o_knuth_class(p):
        f = maximal_increasing_factors(w)
        mu = tuple(sorted(map(len, f), reverse=True))
        if mu == nu:
            j, k = len(f[0]), len(f[-1])
            assert j not in increasing_left_keys or increasing_left_keys[j] == f[0]
            assert k not in increasing_right_keys or increasing_right_keys[k] == f[-1]
            increasing_left_keys[j] = f[0]
            increasing_right_keys[k] = f[-1]

    incr_left_key = []
    for k, w in increasing_left_keys.items():
        for a in w:
            while a > len(incr_left_key):
                incr_left_key.append(0)
            incr_left_key[a - 1] += 1

    incr_right_key = []
    for k, w in increasing_right_keys.items():
        for a in w:
            while a > len(incr_right_key):
                incr_right_key.append(0)
            incr_right_key[a - 1] += 1

    decreasing_left_keys, decreasing_right_keys = {}, {}
    for w in o_knuth_class(p):
        f = maximal_decreasing_factors(w)
        mu = tuple(sorted(map(len, f), reverse=True))
        if mu == nu:
            j, k = len(f[0]), len(f[-1])
            assert j not in decreasing_left_keys or decreasing_left_keys[j] == f[0]
            assert k not in decreasing_right_keys or decreasing_right_keys[k] == f[-1]
            decreasing_left_keys[j] = f[0]
            decreasing_right_keys[k] = f[-1]

    decr_left_key = []
    for k, w in decreasing_left_keys.items():
        for a in w:
            while a > len(decr_left_key):
                decr_left_key.append(0)
            decr_left_key[a - 1] += 1

    decr_right_key = []
    for k, w in decreasing_right_keys.items():
        for a in w:
            while a > len(decr_right_key):
                decr_right_key.append(0)
            decr_right_key[a - 1] += 1

    return tuple(incr_left_key), tuple(incr_right_key), tuple(decr_left_key), tuple(decr_right_key)


def shifted_tableau_tex(p):
    def label(w):
        return sorted(map(len, maximal_increasing_factors(w)))

    def print_group(group):
        s = ''
        if not any(len(f[0]) == len(g[0]) and f[0] != g[0] for f in group for g in group) and not any(len(f[-1]) == len(g[-1]) and f[-1] != g[-1] for f in group for g in group):
            for f in group:
                s += '(' + ', '.join(map(cstr, f)) + ')\\\\\n'
            s += '\\\\ \n'
        return s

    s = []
    for i, j in sorted(p.mapping):
        if len(s) < i:
            s += [(i - 1) * ['\\none[\\cdot]']]
        s[i - 1] += [str(p.entry(i, j))]
    s = '\\\\\n'.join(' & '.join(row) for row in reversed(s))
    s = '\\begin{ytableau}\n' + s + '\n\\end{ytableau}'
    s = '$\\begin{array}{l}\n' + s + '\n\\\\ \\\\ \n'
    nu, qu = symmetric_halves(symmetric_double(p.partition().tuple()))
    nu = tuple(sorted(nu, reverse=True))
    qu = tuple(sorted(qu, reverse=True))
    while nu and nu[-1] == 0:
        nu = nu[:-1]
    while qu and qu[-1] == 0:
        qu = qu[:-1]
    klass = sorted(o_knuth_class(p), key=label)
    for w in klass:
        mu = tuple(sorted(map(len, maximal_increasing_factors(w)), reverse=True))
        if mu == nu:
            s += '(' + ', '.join(map(cstr, maximal_increasing_factors(w))) + ')\\\\\n'
    s += '\\\\ \n'
    group = []
    for i, w in enumerate(klass):
        if i > 0 and label(klass[i]) != label(klass[i - 1]):
            s += print_group(group)
            group = []
        group += [maximal_increasing_factors(w)]
        # mu = tuple(sorted(map(len, maximal_increasing_factors(w)), reverse=True))
        # if True:
        #     s += '(' + ', '.join(map(cstr, maximal_increasing_factors(w))) + ')\\\\\n'
    s += print_group(group)
    s += '\\\\ \n'
    for w in klass:
        mu = tuple(sorted(map(len, maximal_decreasing_factors(w)), reverse=True))
        if mu == nu:
            s += '(' + ', '.join(map(cstr, maximal_decreasing_factors(w))) + ')\\\\\n'
    s += '\\end{array}$'
    return s


def test_q_key_compatible_sequences(m=4):
    def test_q_key(p):
        ans = 0
        for seq in o_knuth_class(p):
            assert o_eg_insert(seq)[0] == p
            for a, x in compatible_sequences(seq):
                ans += x
        return ans

    seen = set()
    keys = {}
    for z in Permutation.involutions(m):
        for w in z.get_involution_words():
            p = o_eg_insert(w)[0]
            if p in seen:
                continue
            e = len(p.partition())
            kappa = test_q_key(p) * 2**e
            # print(p)
            d = try_to_decompose_q(kappa, q_halves_cache, q_alphas_cache, positive=True, multiple=False)
            # for decomp in d:
            #     print('kappa ->', decomp)
            # print()
            seen.add(p)
            assert len(d) > 0
            assert all(len(dec) == 1 for dec in d)
            assert all(v == 1 for dec in d for v in dec.values())
            alpha = list(list(d)[0].keys())[0]
            if alpha not in keys:
                keys[alpha] = []
            keys[alpha] += [p]
    for alpha in keys:
        x, y = symmetric_halves(alpha)
        print(alpha, '->', x, y)
        for p in keys[alpha]:
            a, b, c, d = o_shifted_keys(p)
            print(p)
            print('increasing keys:', a, b)
            print('decreasing keys:', c, d)
            print(b == y)
            if b != y:
                input('?')
            print()
        print()
    print('keys:', len(keys))
    # print_shifted_tex(keys)
    # return keys


def cstr(x):
    return ''.join([str(_) for _ in x])


def print_shifted_tex(keys):
    clip = ''

    def cprint(*args):
        return ' '.join(args) + '\n'

    for alpha in sorted(keys, key=lambda t: (t, len(sorting_permutation(t)))):
        if sum(alpha) == 0:
            continue
        #if len(keys[alpha]) <= 1:
        #    continue
        clip += cprint('\\begin{tabular}[t]{cc}')
        clip += cprint('$T$ & $\\alpha$ \\\\ \\hline \\\\')
        a, b = symmetric_halves(alpha)
        for i, p in enumerate(sorted(keys[alpha], key=lambda t: t.row_reading_word())):
            # if tuple(sorted([x for x in a if x], reverse=True)) == p.partition().tuple():
            #    continue
            clip += cprint(shifted_tableau_tex(p))
            clip += cprint('&')
            if i == 0:
                clip += cprint('$\\begin{array}{l}')
                clip += cprint(cstr(alpha), '\\\\')
                clip += cprint(cstr(a), '\\\\')
                clip += cprint(cstr(b), '\\end{array}$')
            clip += cprint('\\\\ \\\\')
        clip += cprint('\\end{tabular}')
        clip += cprint('\\newpage')
    pyperclip.copy(clip)


def print_tex(keys):
    for alpha in sorted(keys, key=lambda t: (len(keys[t][0]), keys[t][0].row_reading_word())):
        if len(keys[alpha]) <= 1:
            continue
        for i, p in enumerate(keys[alpha]):
            print(tableau_tex(p))
            print('&')
            if i == 0:
                print(cstr(alpha))
            print('\\\\ \\\\')


def sagan_worley_insert(sequence):
    return Word(*sequence).sagan_worley_insert()


def inverse_sagan_worley(p, q):
    w = Tableau.inverse_sagan_worley(p, q)
    return w


def test_inverse_sagan_worley(n=5):
    for word in Word.all(n):
        w = word.elements
        p, q = sagan_worley_insert(w)
        assert w == inverse_sagan_worley(p, q)


def rsk_insert(sequence):
    return Word(*sequence).rsk_insert()


def inverse_rsk(p, q):
    return Tableau.inverse_rsk(p, q)


def knuth_class(p):
    if type(p) != Tableau:
        p = rsk_insert(p)[0]
    mu = p.partition().tuple()
    for q in Tableau.standard(mu):
        yield inverse_rsk(p, q)


def shifted_knuth_class(p):
    if type(p) == Tableau:
        p = p.row_reading_word()
    seen = set()
    add = {p}
    while add:
        nextadd = set()
        for w in add:
            for v in knuth_class(w):
                if v not in seen:
                    seen.add(v)
                    if len(v) >= 2:
                        u = (v[1], v[0]) + v[2:]
                        if u not in seen:
                            nextadd.add(u)
                    yield v
        add = nextadd


def coxeter_knuth_class(p):
    def toggle(w):
        for i in range(len(w) - 2):
            if w[i] == w[i + 2]:
                yield w[:i] + (w[i + 1], w[i], w[i + 1]) + w[i + 3:]
            if w[i] < w[i + 2] < w[i + 1] or w[i + 1] < w[i + 2] < w[i]:
                yield w[:i] + (w[i + 1], w[i], w[i + 2]) + w[i + 3:]
            if w[i + 2] < w[i] < w[i + 1] or w[i + 1] < w[i] < w[i + 2]:
                yield w[:i] + (w[i], w[i + 2], w[i + 1]) + w[i + 3:]

    if type(p) == Tableau:
        p = p.row_reading_word()
    seen = set()
    add = {p}
    while add:
        nextadd = set()
        for v in add:
            if v not in seen:
                seen.add(v)
                for u in toggle(v):
                    if u not in seen:
                        nextadd.add(u)
                yield v
        add = nextadd


def sp_knuth_class(p):
    def toggle(w):
        if len(w) >= 2:
            if w[0] % 2 == w[1] % 2:
                yield (w[1], w[0]) + w[2:]
            if w[1] == w[0] - 1:
                yield (w[0], w[0] + 1) + w[2:]
            if w[1] == w[0] + 1:
                yield (w[0], w[0] - 1) + w[2:]
        for i in range(len(w) - 2):
            if w[i] == w[i + 2]:
                yield w[:i] + (w[i + 1], w[i], w[i + 1]) + w[i + 3:]
            if w[i] < w[i + 2] < w[i + 1] or w[i + 1] < w[i + 2] < w[i]:
                yield w[:i] + (w[i + 1], w[i], w[i + 2]) + w[i + 3:]
            if w[i + 2] < w[i] < w[i + 1] or w[i + 1] < w[i] < w[i + 2]:
                yield w[:i] + (w[i], w[i + 2], w[i + 1]) + w[i + 3:]

    if type(p) == Tableau:
        p = p.row_reading_word()
    seen = set()
    add = {p}
    while add:
        nextadd = set()
        for v in add:
            if v not in seen:
                seen.add(v)
                for u in toggle(v):
                    if u not in seen:
                        nextadd.add(u)
                yield v
        add = nextadd


def o_knuth_class(p):
    if type(p) == Tableau:
        p = p.row_reading_word()
    if p in o_reduced_tableau_cache:
        return o_reduced_tableau_cache[p]

    def toggle(w):
        if len(w) >= 2:
            yield (w[1], w[0]) + w[2:]
        for i in range(len(w) - 2):
            if w[i] == w[i + 2]:
                yield w[:i] + (w[i + 1], w[i], w[i + 1]) + w[i + 3:]
            if w[i] < w[i + 2] < w[i + 1] or w[i + 1] < w[i + 2] < w[i]:
                yield w[:i] + (w[i + 1], w[i], w[i + 2]) + w[i + 3:]
            if w[i + 2] < w[i] < w[i + 1] or w[i + 1] < w[i] < w[i + 2]:
                yield w[:i] + (w[i], w[i + 2], w[i + 1]) + w[i + 3:]

    ans = []
    if type(p) == Tableau:
        p = p.row_reading_word()
    seen = set()
    add = {p}
    while add:
        nextadd = set()
        for v in add:
            if v not in seen:
                seen.add(v)
                for u in toggle(v):
                    if u not in seen:
                        nextadd.add(u)
                ans.append(v)
        add = nextadd
    o_reduced_tableau_cache[p] = ans
    return ans


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
    invol = list(Permutation.involutions(n))
    for i, w in enumerate(invol):
        print('. . .', len(invol) - i)
        # for dream in w.get_involution_pipe_dreams():
        #     word = dream.word()
        #     p, q = o_eg_insert(word)
        #     q_insertion_cache[(w, p)] = q_insertion_cache.get((w, p), 0) + dream.inv_monomial() * 2**e
        for a in w.get_atoms():
            for dream in a.get_pipe_dreams():
                word = dream.word()
                p, q = o_eg_insert(word)
                e = len(p.partition())
                q_insertion_cache[(w, p)] = q_insertion_cache.get((w, p), 0) + dream.monomial() * 2**e
    keys = {}
    for w, p in sorted(q_insertion_cache, key=lambda t: t[1].partition()):
        f = q_insertion_cache[(w, p)]
        d = try_to_decompose_q(f, q_halves_cache, q_alphas_cache, positive=positive, multiple=multiple)
        for decomp in d:
            print(w, '->', decomp)
        print()
        print(p)
        print()
        if len(d) == 0:
            print('keys:', decompose_into_keys(f), '\n***\n')
        assert len(d) >= 1
        assert all(len(decomp) == 1 for decomp in d)
        alpha = list(list(d)[0].keys())[0]
        if alpha not in keys:
            keys[alpha] = []
        keys[alpha] += [p]
    print_shifted_tex(keys)


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
                assert kappa[dict_from_tuple(a)] >= 2**q_power(alpha)
                assert kappa[dict_from_tuple(b)] == 2**q_power(alpha)
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
