from keys import (
    weak_compositions,
    monomial_from_composition,
    sorting_permutation,
    symmetric_weak_compositions,
    skew_symmetric_weak_compositions,
    key, atom,
    q_power,
    p_key, p_atom, p_lascoux,
    q_key, q_atom, q_lascoux,
    get_exponents,
    decompose_key,
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
    symmetric_half,
    print_skew_symmetric_diagram,
    print_symmetric_diagram,
    symmetrize_strict_partition,
    skew_symmetrize_strict_partition,
)
from symmetric import FPFStanleyExpander
from schubert import Schubert, InvSchubert, FPFSchubert
from permutations import Permutation
from partitions import Partition, StrictPartition
from collections import defaultdict
from words import Word
from marked import MarkedNumber
from vectors import Vector
from tableaux import Tableau

from stable.utils import GP, GQ, Q, P, G_expansion, GP_expansion, SymmetricPolynomial

import pyperclip
import pytest
import time


q_alphas_cache = {}
q_lascoux_alphas_cache = {}
q_halves_cache = {}

p_alphas_cache = {}
p_lascoux_alphas_cache = {}
p_halves_cache = {}

q_insertion_cache = {}
p_insertion_cache = {}


def test_p_lascoux_to_gp(n=5):
    k = 0
    while n is None or k <= n:
        print(k)
        for mu in StrictPartition.all(k):
            mu = tuple(mu)
            lam = tuple(reversed(skew_symmetrize_strict_partition(mu)))
            nvars = max((0,) + mu) + 1
            print('  ', mu, '-->', lam, ':', nvars)
            f = p_lascoux(lam)
            g = GP(nvars, mu).polynomial()
            if f != g:
                print()
                print('FAIL,', 'L < GP :', f < g)
                print()
        k += 1


def test_q_lascoux_to_gq(n=5):
    nch = 150
    k = 0
    while n is None or k <= n:
        print(k)
        for mu in StrictPartition.all(k):
            mu = tuple(mu)
            lam = tuple(reversed(symmetrize_strict_partition(mu)))
            nvars = max((0,) + mu)
            print('  ', mu, '-->', lam, ':', nvars)
            f = q_lascoux(lam)
            gg = GQ(nvars, mu)
            g = gg.polynomial()
            if f != g:
                ff = SymmetricPolynomial.from_polynomial(f)
                print()
                print('FAIL,', 'L < GQ :', f < g, 'nvars =', nvars)
                print('    ', sorted(list(G_expansion(GQ(nvars, mu) - ff)), key=len))
                try:
                    print('    ', GP_expansion(ff))
                except:
                    print('    ', '(no GP expansion)')
                # print()
                # print(' L =', f)
                # print('GQ =', g)
                print()
        k += 1


def test_q_key_to_q(n=5):
    nch = 150
    k = 0
    while n is None or k <= n:
        print(k)
        for mu in StrictPartition.all(k):
            mu = tuple(mu)
            lam = tuple(reversed(symmetrize_strict_partition(mu)))
            nvars = max((0,) + mu)
            print('  ', mu, '-->', lam, ':', nvars)
            f = q_key(lam)
            g = Q(nvars, mu).polynomial()
            if f != g:
                print()
                print('FAIL,', 'k < Q :', f < g)
                print()
        k += 1


def test_distinct_atom(m=4, l=4):
    seen = {}
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in weak_compositions(n, k, reduced=True):
                f = atom(alpha)
                if f not in seen:
                    seen[f] = set()
                seen[f].add(alpha)
                assert len(seen[f]) == 1


def test_distinct_p_key(m=4, l=4):
    seen = {}
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
                f = p_key(alpha)
                if f not in seen:
                    seen[f] = set()
                seen[f].add(alpha)
                try:
                    assert len(seen[f]) == 1
                except:
                    print(seen[f], '-->', str(f)[:20])


def test_distinct_p_atom(m=4, l=4):
    seen = {}
    nonzero = 0
    zero = 0
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in skew_symmetric_weak_compositions(n, k, reduced=True):
                f = p_atom(alpha)
                if f == 0:
                    # print('zero:', alpha)
                    zero += 1
                    continue
                nonzero += 1
                if f not in seen:
                    seen[f] = set()
                seen[f].add(alpha)
                try:
                    assert len(seen[f]) == 1
                except:
                    print(seen[f], '-->', str(f)[:20])
    print(zero, 'nonzero:', nonzero)


def test_distinct_q_atom(m=4, l=4):
    seen = {}
    nonzero = 0
    zero = 0
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                f = q_atom(alpha)
                if f == 0:
                    # print('zero:', alpha)
                    zero += 1
                    continue
                nonzero += 1
                if f not in seen:
                    seen[f] = set()
                seen[f].add(alpha)
                try:
                    assert len(seen[f]) == 1
                except:
                    print(seen[f], '-->', str(f)[:20])
    print(zero, 'nonzero:', nonzero)


def test_distinct_q_key(m=4, l=4):
    seen = {}
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                f = q_key(alpha)
                if f not in seen:
                    seen[f] = set()
                seen[f].add(alpha)
                try:
                    assert len(seen[f]) == 1
                except:
                    print(seen[f], '-->', str(f)[:20])
                if len(seen) % 1000 == 0:
                    print('  seen:', len(seen))


def test_q_key_into_p_key(m=4, l=4):
    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                attempt = try_to_decompose_p(q_key(alpha))
                print(alpha, '->', attempt)

def flags(length, upperbound, lowerbound=0):
    if length == 0:
        yield ()
        return
    for i in range(lowerbound, upperbound):
        for f in flags(length - 1, upperbound, max(i, lowerbound + 1)):
            yield (i,) + f


def q_key_tableau(w, alpha):
    positions = [(i + 1, j) for i, v in enumerate(alpha) for j in range(1, v + 1)]
    sigma = Permutation.from_word(*sorting_permutation(alpha))
    positions = [(i, sigma(j)) for i, j in positions]
    positions = [(i, j) for i, j in positions if i <= j]
    positions = sorted(positions, key=lambda x: (x[1], -x[0]))
    mapping = {}
    for i, b in enumerate(positions):
        mapping[b] = w[i]
    return Tableau(mapping)


# def q_key_tableau(w, alpha):
#     _, cc = symmetric_halves(alpha)
#     positions = [(j, i + 1) for i, v in enumerate(cc) for j in range(v, 0, -1)]
#     mapping = {}
#     for i, b in enumerate(positions):
#         mapping[b] = w[i]
#     return Tableau(mapping)


def test_flagged_key(n=3):
    comps = set()
    seen = set()

    def add_comps(alpha):
        p = sum(alpha)
        q = len(alpha)
        if (p, q) not in seen:
            for beta in weak_compositions(p, q, reduced=True):
                mu = sorted(beta, reverse=True)
                if all(x == 0 or x <= n - i - 1 for i, x in enumerate(mu)):
                    comps.add(beta)
            seen.add((p, q))

    alphas = set()
    for i in Permutation.all(n):
        ck = defaultdict(list)
        for w in i.get_reduced_words():
            p, _ = eg_insert(w)
            ck[p].append(w)
        for p, ckclass in ck.items():
            print(p)
            m = max(ckclass[0]) if ckclass[0] else 0
            for f in flags(1 + m, 3 + n):
                ans = 0
                exp = defaultdict(list)
                for w in ckclass:
                    for _, x in compatible_sequences(w, flag=f):
                        exp[get_exponents(x)[0]].append(w)
                        ans += x

                alpha = decompose_key(ans)
                alphas.add(alpha)
                add_comps(alpha)

                print('flag:', f, '->', alpha, 'w =', w)
            print()
            print('missing:', comps - alphas)
            print()


def test_flagged_p_key(n=3):
    skew = set()
    seen = set()
    alphas = set()
    values = set()
    for i in Permutation.fpf_involutions(n):
        ck = defaultdict(list)
        for w in i.get_fpf_involution_words():
            p, _ = sp_eg_insert(w)
            ck[p].append(w)
        for p, ckclass in ck.items():
            print(p)
            m = max(ckclass[0]) if ckclass[0] else 0
            for f in flags(1 + m, 2 + n):
                ans = 0
                exp = defaultdict(list)
                for w in ckclass:
                    for _, x in compatible_sequences(w, flag=f):
                        exp[get_exponents(x)[0]].append(w)
                        ans += x
                alpha = decompose_p(ans)
                alphas.add(alpha)
                values.add(ans)

                p = sum(alpha)
                q = len(alpha)
                if (p, q) not in seen:
                    skew |= set(skew_symmetric_weak_compositions(p, q, reduced=True))
                    skew = {a for a in skew if p_key(a) not in values}
                    seen.add((p, q))

                rc, cc = skew_symmetric_halves(alpha)
                print('flag:', f, '->', alpha, ':', rc, cc, 'w =', w)
            print()
            print('missing:', skew - alphas)
            print()


def test_flagged_q_key(n=3):
    sym = set()
    seen = set()
    alphas = set()
    for i in Permutation.involutions(n):
        ck = defaultdict(list)
        va = {}
        for w in i.get_involution_words():
            p, _ = o_eg_insert(w)
            ck[p].append(w)
            va[p] = p.durfee()
        for p, ckclass in ck.items():
            m = max(ckclass[0]) if ckclass[0] else 0
            print(p)
            for f in flags(1 + m, 2 + n):
                ans = 0
                exp = defaultdict(list)
                for w in ckclass:
                    for _, x in compatible_sequences(w, flag=f):
                        exp[get_exponents(x)[0]].append(w)
                        ans += x
                ans *= 2**va[p]
                alpha = decompose_q(ans)
                alphas.add(alpha)

                u = sum(alpha)
                v = len(alpha)
                if (u, v) not in seen:
                    sym |= set(symmetric_weak_compositions(u, v, reduced=True))
                    seen.add((u, v))

                rc, cc = symmetric_halves(alpha)
                print('flag:', f, '->', alpha, ':', rc, cc, 'w =', w)
            print()
            print('missing:', sym - alphas)
            print()


def _flagged_q_key_insertion(n, p, va, ckclass):
    # if len(p.partition()) < 2:
    #    continue
    m = max(ckclass[0]) if ckclass[0] else 0
    for f in flags(1 + m, n):
        ans = 0
        exp = defaultdict(list)
        for w in ckclass:
            for _, x in compatible_sequences(w, flag=f):
                exp[get_exponents(x)[0]].append(w)
                ans += x
        ans *= 2**va
        alpha = decompose_q(ans)

        rc, cc = symmetric_halves(alpha)
        assert len(exp[cc]) == 1
        w = exp[cc][0]

        # scc = sorted(cc)
        # while scc and scc[0] == 0:
        #     scc = scc[1:]
        # if scc == sorted(symmetric_halves(sorted(alpha, reverse=True))[1]):
        #     return

        p, q = o_eg_insert(w)
        r = q_key_tableau(w, alpha)

        # if r != p:
        print('flag:', f, '->', alpha, ':', rc, cc)
        print()
        print('cc:', cc, '->', w)
        print(r)
        print(p)
        print(q)

        # print()
        # for e in exp:
        #     i = tuple(j + 1 for j, v in enumerate(e) for _ in range(v))
        #     for v in exp[e]:
        #         m = (1 + max(v)) if v else 0
        #         v = tuple(m - j for j in v)
        #         _, q = o_eg_insert(v, i)
        #         print(q)
        # print()
        break


def test_flagged_q_key_insertion(n=3):
    if type(n) == Tableau:
        word = n.row_reading_word()
        p, ckclass = n, o_knuth_class(word)
        n = (max(word) if word else 0) + 1
        _flagged_q_key_insertion(n, p, p.durfee(), ckclass)
        return

    for i in Permutation.involutions(n):
        ck = defaultdict(list)
        va = {}
        for w in i.get_involution_words():
            p, _ = o_eg_insert(w)
            ck[p].append(w)
            va[p] = p.durfee()
        for p, ckclass in ck.items():
            _flagged_q_key_insertion(n, p, va[p], ckclass)


def test_flagged_p_key_independence(n=4):
    subset = []
    progress = []
    alphas = defaultdict(list)
    seen = 0
    inv = list(Permutation.fpf_involutions(n))
    for e, i in enumerate(inv):
        ck = defaultdict(list)
        for w in i.get_fpf_involution_words():
            p, _ = sp_eg_insert(w)
            ck[p].append(w)
        for p, ckclass in ck.items():
            m = max(ckclass[0]) if ckclass[0] else 0
            for f in [None]: #flags(1 + m, 2 + n):
                ans = 0
                for w in ckclass:
                    for _, x in compatible_sequences(w, flag=f):
                        ans += x
                alpha = decompose_p(ans)
                if alpha not in alphas:
                    subset.append(ans)
                    t0 = time.time()
                    p0 = len(progress)
                    progress = Vector.reduce_linearly_independent_subset(subset, progress)
                    p1 = len(progress)
                    assert all(m is not None for m, v in progress[p0:p1])
                    t1 = time.time()
                    print('vectors:', len(subset), 'of', seen, 'seen, left:', len(inv) - e, '| independence check took %s milliseconds' % int(1000 * (t1 - t0)))
                seen += 1
                alphas[alpha].append((p, f))
                rc, cc = skew_symmetric_halves(alpha)
    assert Vector.is_linearly_independent_subset(subset)
    return subset, alphas


def test_flagged_q_key_independence(n=3):
    subset = []
    progress = []
    alphas = defaultdict(list)
    seen = 0
    inv = list(Permutation.involutions(n))
    for e, i in enumerate(inv):
        ck = defaultdict(list)
        va = {}
        for w in i.get_involution_words():
            p, _ = o_eg_insert(w)
            ck[p].append(w)
            va[p] = p.durfee()
        for p, ckclass in ck.items():
            m = max(ckclass[0]) if ckclass[0] else 0
            for f in [None]: #flags(1 + m, 2 + n):
                ans = 0
                for w in ckclass:
                    for _, x in compatible_sequences(w, flag=f):
                        ans += x
                ans *= 2**va[p]
                alpha = decompose_q(ans)
                if alpha not in alphas:
                    subset.append(ans)
                    t0 = time.time()
                    p0 = len(progress)
                    progress = Vector.reduce_linearly_independent_subset(subset, progress)
                    p1 = len(progress)
                    assert all(m is not None for m, v in progress[p0:p1])
                    t1 = time.time()
                    print('vectors:', len(subset), 'of', seen, 'seen, left:', len(inv) - e, '| independence check took %s milliseconds' % int(1000 * (t1 - t0)))
                seen += 1
                alphas[alpha].append((p, f))
                rc, cc = symmetric_halves(alpha)
    assert Vector.is_linearly_independent_subset(subset)
    return subset, alphas


def morse_schilling_f(decreasing_factorization, bounded=True):
    for i in range(len(decreasing_factorization)):
        if i == len(decreasing_factorization) - 1:
            decreasing_factorization += ((),)

        a = decreasing_factorization[i]
        b = decreasing_factorization[i + 1]
        unpaired = []

        while a:
            paired = False
            for j in range(len(b) - 1, -1, -1):
                if b[j] > a[0]:
                    b = b[:j] + b[j + 1:]
                    paired = True
                    break
            if not paired:
                unpaired += [a[0]]
            a = a[1:]

        if len(unpaired) == 0:
            continue

        x = unpaired[-1]
        y = x
        while y in decreasing_factorization[i + 1]:
            y -= 1
        df = list(decreasing_factorization)
        df[i] = tuple(e for e in df[i] if e != x)
        df[i + 1] = tuple(sorted(df[i + 1] + (y,), reverse=True))

        while df and df[-1] == ():
            df = df[:-1]

        if not bounded or all(df[i][-1] > i for i in range(len(df)) if df[i]):
            yield i, tuple(df)


def is_morse_schilling_highest_weight(decreasing_factorization):
    for i in range(len(decreasing_factorization) - 1):
        a = decreasing_factorization[i]
        b = decreasing_factorization[i + 1]
        while a and b:
            for j in range(len(b) - 1, -1, -1):
                if b[j] > a[0]:
                    b = b[:j] + b[j + 1:]
                    break
            a = a[1:]
        if len(b) != 0:
            return False
    return True


def get_morse_schilling_highest_weights(it):
    for word in it:
        df = maximal_decreasing_factors(word)
        if is_morse_schilling_highest_weight(df):
            yield df


def get_morse_schilling_demazure_weights(it):
    def weight(u):
        ans = tuple(map(len, u))
        while ans and ans[-1] == 0:
            ans = ans[:-1]
        return ans

    def sorted_weight(u):
        ans = tuple(sorted(weight(u), reverse=True))
        while ans and ans[-1] == 0:
            ans = ans[:-1]
        return ans

    for u, graph in get_morse_schilling_demazure(it):
        mu = weight(u)
        levels = []
        add = {u}
        while add:
            levels += [[]]
            next_add = set()
            for v in add:
                if sorted_weight(v) == mu:
                    levels[-1] += [v]
                for i in graph[v]:
                    w = graph[v][i]
                    next_add.add(w)
            add = next_add
        while levels and levels[-1] == []:
            levels = levels[:-1]
        assert len(levels[-1]) == 1
        yield levels[-1][0]


def get_morse_schilling_demazure(it):
    highest = sorted(get_morse_schilling_highest_weights(it))
    for u in highest:
        graph = {}
        add = {u}
        while add:
            next_add = set()
            for w in add:
                graph[w] = {}
                for i, df in morse_schilling_f(w):
                    graph[w][i] = df
                    next_add.add(df)
            add = next_add
        yield u, graph


def remove_hook(mu):
    h = (mu[0],) + (mu[0] - 1) * (1,)
    mu = list(mu)[1:]
    mu = [a - 1 for a in mu]
    while mu and mu[-1] == 0:
        mu = mu[:-1]
    mu = tuple(mu)
    return mu, h


def test_q_dominant(n0=5):
    weights = {}

    def get_weight(mu):
        if mu not in weights:
            rc, cc = symmetric_halves(mu)
            z = Permutation.from_involution_code(cc)
            it = z.get_involution_words()
            weights[mu] = sorted(get_morse_schilling_demazure_weights(it))
        return weights[mu]

    n = 0
    while True:
        n += 1
        if n > n0:
            break
        print('n =', n)
        print()
        for mu in symmetric_partitions(n):
            level = 0

            print('key =', decompose_into_keys(q_key(mu)))
            print()
            print('level', level, ': mu =', mu)
            print()
            for w in get_weight(mu):
                print('  ', w)
            print()

            mu, h = remove_hook(mu)
            level += 1

            print('level', level, ': mu =', mu, 'hook =', h)
            print()
            for w in get_weight(mu):
                print('  ', w)
            print()

            for w in get_weight(h):
                print('  ', w)
            print()

            # input('')


def test_p_multiplicity_free(n0=5):
    print()
    n = 0
    while True:
        n += 1
        if n > n0:
            break
        print(n)
        for mu in skew_symmetric_partitions(n):
            a, b = skew_symmetric_halves(mu)
            f = p_key(mu)
            dec = decompose_into_keys(f)
            if set(dec.values()) != {1}:
                print(mu, '->', a, b, ':', {x for x in dec if dec[x] != 1})
                print()
                return dec


def test_q_multiplicity_free(n0=5):
    print()
    n = 0
    while True:
        n += 1
        if n > n0:
            break
        print(n)
        for mu in symmetric_partitions(n):
            a, b = symmetric_halves(mu)
            f = q_key(mu)
            dec = decompose_into_keys(f)
            if len(set(dec.values())) > 1:
                print(mu, '->', a, b, ':', {x for x in dec if dec[x] != min(dec.values())})
                print()
                return dec


def test_leading_with_atoms(m=6):
    for n in range(m + 1):
        for mu in symmetric_partitions(n):
            b, c = symmetric_halves(mu)
            z = Permutation.from_involution_code(c)
            s = []
            for a in z.get_atoms():
                u = a.get_bottom_pipe_dream().weight()
                v = a.get_top_pipe_dream().weight()
                s += [(u, v, sorted(decompose_into_keys(Schubert.get(a))))]
            print()
            print('lambda =', mu, '->', c, b)
            for u, v, line in s:
                print(u, ':', line, ':', v)
            print()
            # input('?')


def test_leading(m=6, l=6, p=None):
    def toggle(a, i):
        if a[i - 1] <= a[i]:
            return a
        b = list(a)
        b[i - 1], b[i] = a[i], a[i - 1]
        return tuple(b)

    for n in range(m + 1):
        for k in range(l + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                word = tuple(reversed(sorting_permutation(alpha)))
                mu = tuple(reversed(sorted(alpha)))
                lam = {(i + 1, j + 1) for i in range(len(mu)) for j in range(mu[i])}
                alphas = sorted(decompose_into_keys(q_key(mu)))
                if p is not None and len(alphas) < p:
                    continue
                alphas = [a + (len(alpha) - len(a)) * (0,) for a in alphas]
                print('alpha =', alpha, '<-', mu, ':', word)
                print()
                print(Tableau({a: '0' for a in lam}))
                print(' '.join(map(str, alphas)))
                for i in word:
                    alphas = sorted(toggle(a, i) for a in alphas)
                    s = Permutation.s_i(i)
                    lam = {(s(a), s(b)) for (a, b) in lam}
                    print()
                    print(Tableau({a: '0' for a in lam}))
                    print(' '.join(map(str, alphas)))
                print()
                print()
                # input('')


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


def test_atom_operators_on_keys(k=2):
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


def test_key_compatible_sequences(m=4, l=4):
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
    # FAILS
    def test_shifted_key(p):
        ans = 0
        cl = list(shifted_knuth_class(p))
        for seq in cl:
            assert primed_sw_insert(seq)[0] == p
            seq = tuple(MarkedNumber(v) for v in seq)
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
                p = primed_sw_insert(w)[0]
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
                kappa = test_shifted_key(p) # * 2**len(p.partition())
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


def primed_sw_insert(sequence):
    return Word(*sequence).primed_sw_insert()


def inverse_sagan_worley(p, q):
    w = Tableau.inverse_sagan_worley(p, q)
    return w


def test_inverse_sagan_worley(n=5):
    for w in words(n):
        p, q = sagan_worley_insert(w)
        assert w == inverse_sagan_worley(p, q)


def test_shifted_knuth_class(n=4):
    # FAILS
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


def o_eg_insert(sequence, phi=None):
    return Word(*sequence).involution_insert(phi=phi)


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
            # if c.weight() != alpha:
            #    input('\n?')
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
            # if c.weight() != alpha:
            #    input('\n?')
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


def modify_row_counts(alpha, rc, up=True):
    rc += (len(alpha) - len(rc)) * (0,)
    mc = list(rc)
    repeat = True
    while repeat:
        repeat = False
        for i in range(len(rc) - 1):
            for j in range(i + 1, len(rc)):
                if alpha[i] == alpha[j] and ((up and mc[i] > mc[j]) or (not up and mc[i] < mc[j])):
                    mc[i], mc[j] = mc[j], mc[i]
                    repeat = True
    while mc and mc[-1] == 0:
        mc = mc[:-1]
    return tuple(mc)


def test_q_key_start_decomposition(m=5):
    for n in range(m + 1):
        for alpha in symmetric_partitions(n):
            kappa = q_key(alpha)
            dec = decompose_into_keys(kappa)
            ex = min({0} | set(dec.values()))
            assert ex >= 0
            assert kappa != 0

            rc, cc = symmetric_halves(alpha)
            mc = modify_row_counts(alpha, rc)
            tc = modify_row_counts(alpha, rc, False)

            print('alpha =', alpha, '->', cc, rc, 'modified =', mc, rc == tc)
            #print('decomposition =', sorted(dec))
            print()
            assert cc == min(dec)
            assert mc == max(dec)
            assert dec[cc] == 2**q_power(alpha)
            assert dec[mc] == 2**q_power(alpha)
            assert rc == tc
            # assert set(dec.values()) == {2**q_power(alpha)}


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
    # FAILS
    for n in range(m + 1):
        for k in range(m + 1):
            for alpha in symmetric_weak_compositions(n, k, reduced=True):
                if any(alpha[i] == alpha[j] for i in range(len(alpha)) for j in range(i + 1, len(alpha))):
                    continue
                kappa = q_key(alpha)

                rc, cc = symmetric_halves(alpha)
                mc = modify_row_counts(alpha, rc)

                dec = decompose_into_keys(kappa)
                ex = min({0} | set(dec.values()))

                is_dom = lambda x :all(x[i] >= x[i+1] for i in range(len(x) - 1)) # noqa
                lam = sorted(filter(is_dom, get_exponents(kappa)))[-1]
                assert lam == symmetric_half(tuple(reversed(sorted(alpha))))

                assert ex >= 0
                assert kappa != 0
                assert cc == min(dec)
                if mc not in dec:
                    print('alpha =', alpha, '->', cc, rc, 'modified =', mc)
                    print()
                    print(sorting_permutation(alpha))
                    print()
                    print(sorted(dec))
                    print()
                    # input('?')
                assert dec[cc] == 2**q_power(alpha)



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
                # input('')
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
                # input('')
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

                print(alpha)
                print_skew_symmetric_diagram(alpha)
                print()
                print(a, b)
                print(exponents[0], exponents[-1])
                print()
                print({e: kappa[dict_from_tuple(a)] for e in exponents})
                input('\n\n')

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
    try:
        assert not any(len(v) > 1 for v in valuesdict.values())
    except:
        for k, v in valuesdict.items():
            if len(v) > 1:
                print(k, '-->', v)
        assert False
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
                print(alpha)
                print_symmetric_diagram(alpha)
                print()
                print(a, b)
                print(exponents[0], exponents[-1])
                print()
                print({e: kappa[dict_from_tuple(a)] for e in exponents})
                # input('\n\n')
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


def q_update(targets, exponents, halves, alphas, functional):
    for e in targets:
        for d in exponents:
            try:
                a = symmetric_composition_from_row_column_counts(d, e)
                assert symmetric_halves(a) == (d, e)
            except:
                continue
            if a not in alphas:
                alphas[a] = functional(a)
                halves[e] = halves.get(e, []) + [(alphas[a], a)]


def _decompose(f, halves, alphas, positive, multiple, functional, update):
    if f == 0:
        return [{}]
    if positive and not f.is_positive():
        return []
    exponents = get_exponents(f)
    targets = [exponents[0]]
    update(targets, exponents, halves, alphas, functional)
    answers = []
    for target in targets:
        dict_key = dict_from_tuple(target)
        for g, alpha in sorted(halves.get(target, [])):
            assert g == functional(alpha)
            a = f[dict_key]
            b = g[dict_key]
            if a % b == 0:
                h = f - a // b * g
                assert h[dict_key] == 0
                for ans in _decompose(h, halves, alphas, positive, multiple, functional, update):
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
            g += coeff * functional(alpha)
        assert f == g
    return answers


def try_to_decompose_q(f, halves=None, alphas=None, positive=True, multiple=False):
    functional = q_key
    update = q_update
    halves = q_halves_cache if halves is None else halves
    alphas = q_alphas_cache if alphas is None else alphas
    return _decompose(f, halves, alphas, positive, multiple, functional, update)
    # if f == 0:
    #     return [{}]
    # if positive and not f.is_positive():
    #     return []
    # exponents = get_exponents(f)
    # targets = [exponents[0]]
    # update(targets, exponents, halves, alphas, functional)
    # answers = []
    # for target in targets:
    #     dict_key = dict_from_tuple(target)
    #     for g, alpha in sorted(halves.get(target, [])):
    #         assert g == functional(alpha)
    #         a = f[dict_key]
    #         b = g[dict_key]
    #         if a % b == 0:
    #             h = f - a // b * g
    #             assert h[dict_key] == 0
    #             for ans in try_to_decompose_q(h, halves, alphas, positive, multiple):
    #                 ans[alpha] = ans.get(alpha, 0) + a // b
    #                 if ans[alpha] == 0:
    #                     del ans[alpha]
    #                 if not multiple:
    #                     return [ans]
    #                 elif ans not in answers:
    #                     answers.append(ans)
    # for ans in answers:
    #     g = 0
    #     for alpha, coeff in ans.items():
    #         g += coeff * functional(alpha)
    #     assert f == g
    # return answers


def try_to_decompose_q_lascoux(f, halves=None, alphas=None, positive=True, multiple=False):
    functional = lambda a: q_lascoux(a).set(0, 1)
    update = q_update
    halves = q_halves_cache if halves is None else halves
    alphas = q_lascoux_alphas_cache if alphas is None else alphas
    return _decompose(f, halves, alphas, positive, multiple, functional, update)
    

def test_inv_schubert(n=4, positive=True, multiple=True):
    i = list(Permutation.involutions(n))
    s = {w: InvSchubert.get(w) * 2**w.number_two_cycles() for w in i}
    print('. . . s')
    d = {}
    for t, w in enumerate(s):
        isvex = w.is_vexillary()
        d[w] = try_to_decompose_q(s[w], q_halves_cache, q_alphas_cache, positive, multiple)
        for dec in d[w]:
            print(len(s) - t, ':', w, '->', dec, isvex, w.code())
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


def p_update(targets, exponents, halves, alphas, functional):
    for e in targets:
        for d in exponents:
            try:
                a = skew_symmetric_composition_from_row_column_counts(d, e)
            except:
                continue
            assert skew_symmetric_halves(a) == (d, e)
            if a not in alphas:
                alphas[a] = functional(a)
                halves[e] = halves.get(e, []) + [(alphas[a], a)]


def try_to_decompose_p(f, halves=None, alphas=None, positive=True, multiple=False):
    functional = p_key
    update = p_update
    halves = p_halves_cache if halves is None else halves
    alphas = p_alphas_cache if alphas is None else alphas
    return _decompose(f, halves, alphas, positive, multiple, functional, update)
    

def is_fpf_vexillary(w):
    f = FPFStanleyExpander(w).expand()
    if len(f) > 1:
        return False
    if set(f.values()) != {1}:
        return False
    return True


def vexify(w):
    for i in range(2, len(w.oneline) + 1):
        j = w(i)
        if j < i and all(w(k) < j for k in range(j + 1, i)):
            w *= Permutation.t_ij(j, i)
    return w


def special_fpf_code(w):
    ans = []
    for i, j in w.rothe_diagram():
        if i == j and w(i) == i + 1: # all(w(a) < i for a in range(i + 1, w(i))): # or all(w(i) < b for (a, b) in w.rothe_diagram() if i < a < w(i))):
            continue
        while i > len(ans):
            ans.append(0)
        ans[i - 1] += 1
    return tuple(ans)


def test_fpf_schubert(n=4, positive=True, multiple=True):
    i = list(Permutation.fpf_involutions(n))
    s = {w: FPFSchubert.get(w) for w in i}
    print('. . . s')
    d = {}
    for t, w in enumerate(s):
        isvex = is_fpf_vexillary(w)
        d[w] = try_to_decompose_p(s[w], p_halves_cache, p_alphas_cache, positive, multiple)
        for dec in d[w]:
            print(len(s) - t, ':', w.cycle_repr(), '->', dec, isvex, w.code())
            # if len(dec) == 1:
            #    print()
            #    print_skew_symmetric_diagram(next(iter(dec)))
            #    print()
            assert set(dec.values()) == {1}
        # w.print_rothe_diagram(sep='.')
        print()
        if any(len(dec) == 1 for dec in d[w]) and not isvex:
            input('\n!\n')
        # if isvex and not any(special_fpf_code(w) in dec for dec in d[w]):
        #    input('\n?\n')
        # assert (not positive and multiple) or len(d[w]) == 1
        d[w] = sorted(d[w], key=lambda x: (len(x), sorted(x.values())))[0]
        if vexify(w).is_vexillary():
            assert len(d[w]) == 1 and set(d[w].values()) == {1}
    pvex = {w: list(d[w])[0] for w in d if len(d[w]) == 1 and set(d[w].values()) == {1}}
    fvex = {w: w.code() for w in i if is_fpf_vexillary(w)}
    try:
        assert set(pvex) == set(fvex)
    except:
        print('p-vexillary not fpf vexillary:')
        for w in set(pvex) - set(fvex):
            print(w, pvex[w])
            print()
        print('fpf-vexillary not p-vexillary:')
        for w in set(fvex) - set(pvex):
            print(w, fvex[w])
            print()
        input('?')
    print()
    print('caches =', len(p_halves_cache), len(p_alphas_cache))
    return i, s, d, pvex, fvex
