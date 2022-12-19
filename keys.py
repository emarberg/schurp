# def tcode(diag):
#     ans = []
#     for (i, j) in diag:
#         if i != j:
#             while i - 1 >= len(ans):
#                 ans += [0]
#             ans[i - 1] += 1
#     return tuple(ans)


# def vex(w):
#     for i in range(2, len(w.oneline) + 1):
#         j = w(i)
#         if j < i and all(w(k) < j for k in range(j + 1, i)):
#             w *= Permutation.t_ij(j, i)
#     return w

# for w in pvex:
#     b = tcode(w.rothe_diagram())==tcode(skew_symmetric_diagram(pvex[w]))
#     print(w, w.cycle_repr(), pvex[w], pvex[w]==fvex[w], w.is_vexillary(), vex(w).is_vexillary())
#     if not vex(w).is_vexillary():
#         w.print_rothe_diagram(sep='.')
#         print()
#         print(vex(w))
#         vex(w).print_rothe_diagram(sep='.')
#         print()
#         print_skew_symmetric_diagram(pvex[w])
#         print()


from schubert import X, InvGrothendieck, InvSchubert
from partitions import Partition, StrictPartition
from tableaux import Tableau
from words import Word
from marked import MarkedNumber
from permutations import Permutation
import itertools


KEY_POLYNOMIAL_CACHE = {}
PKEY_POLYNOMIAL_CACHE = {}
QKEY_POLYNOMIAL_CACHE = {}

KEY_ATOM_CACHE = {}
PKEY_ATOM_CACHE = {}
QKEY_ATOM_CACHE = {}

LASCOUX_POLYNOMIAL_CACHE = {}
PLASCOUX_POLYNOMIAL_CACHE = {}
QLASCOUX_POLYNOMIAL_CACHE = {}

LASCOUX_ATOM_CACHE = {}
PLASCOUX_ATOM_CACHE = {}
QLASCOUX_ATOM_CACHE = {}

REDUCED_TABLEAU_CACHE = {}
O_REDUCED_TABLEAU_CACHE = {}
SP_REDUCED_TABLEAU_CACHE = {}
SHIFTED_REDUCED_TABLEAU_CACHE = {}


LASCOUX_BETA = X(0)


def _get_key_maps(tab, shape, get_class, get_factors):
    left_keys, right_keys = {}, {}
    for w in get_class(tab):
        f = get_factors(w)
        mu = tuple(sorted(map(len, f), reverse=True))
        if mu == shape:
            j, k = len(f[0]), len(f[-1])
            assert j not in left_keys or left_keys[j] == f[0]
            assert k not in right_keys or right_keys[k] == f[-1]
            left_keys[j] = f[0]
            right_keys[k] = f[-1]
    return left_keys, right_keys


def _map_to_shifted_key(shape, key_map):
    return _map_to_key(shape, key_map, True)


def _map_to_key(shape, key_map, shifted=False):
    ans = Tableau()
    for i in range(len(shape)):
        if i + 1 < len(shape):
            a = key_map[shape[i + 1]]
            b = key_map[shape[i]][1:] if shifted else key_map[shape[i]]
            assert set(a).issubset(set(b))
        if shifted:
            for j, a in enumerate(sorted(key_map[shape[i]])):
                ans = ans.add(i + 1, i + j + 1, a)
        else:
            for j, a in enumerate(sorted(key_map[shape[i]])):
                ans = ans.add(j + 1, i + 1, a)
    return ans


def orthogonal_key_maps(p):
    nu = p.partition().tuple()
    increasing_left_keys, increasing_right_keys = _get_key_maps(p, nu, o_knuth_class, maximal_increasing_factors)
    decreasing_left_keys, decreasing_right_keys = _get_key_maps(p, nu, o_knuth_class, maximal_decreasing_factors)
    return (_map_to_shifted_key(nu, increasing_left_keys),
            _map_to_shifted_key(nu, increasing_right_keys),
            _map_to_shifted_key(nu, decreasing_left_keys),
            _map_to_shifted_key(nu, decreasing_right_keys))


def symplectic_key_maps(p):
    nu = p.partition().tuple()
    increasing_left_keys, increasing_right_keys = _get_key_maps(p, nu, sp_knuth_class, maximal_increasing_factors)
    decreasing_left_keys, decreasing_right_keys = _get_key_maps(p, nu, sp_knuth_class, maximal_decreasing_factors)
    return (_map_to_shifted_key(nu, increasing_left_keys),
            _map_to_shifted_key(nu, increasing_right_keys),
            _map_to_shifted_key(nu, decreasing_left_keys),
            _map_to_shifted_key(nu, decreasing_right_keys))


def nil_key_maps(p):
    mu = p.partition().tuple()
    nu = p.partition().transpose().tuple()
    increasing_left_keys, increasing_right_keys = _get_key_maps(p, mu, coxeter_knuth_class, maximal_increasing_factors)
    decreasing_left_keys, decreasing_right_keys = _get_key_maps(p, nu, coxeter_knuth_class, maximal_decreasing_factors)
    return (_map_to_key(mu, increasing_left_keys),
            _map_to_key(mu, increasing_right_keys),
            _map_to_key(nu, decreasing_left_keys),
            _map_to_key(nu, decreasing_right_keys))


def key_maps(p):
    mu = p.partition().tuple()
    nu = p.partition().transpose().tuple()
    increasing_left_keys, increasing_right_keys = _get_key_maps(p, mu, knuth_class, maximal_weakly_increasing_factors)
    decreasing_left_keys, decreasing_right_keys = _get_key_maps(p, nu, knuth_class, maximal_decreasing_factors)
    return (_map_to_key(mu, increasing_left_keys),
            _map_to_key(mu, increasing_right_keys),
            _map_to_key(nu, decreasing_left_keys),
            _map_to_key(nu, decreasing_right_keys))


def shifted_key_maps(p):
    mu = p.partition().tuple()
    increasing_left_keys, increasing_right_keys = _get_key_maps(p, mu, shifted_knuth_class, maximal_weakly_increasing_factors)
    return (_map_to_shifted_key(mu, increasing_left_keys),
            _map_to_shifted_key(mu, increasing_right_keys))


def symmetric_composition_from_row_column_counts(row_counts, col_counts):
    shape = _symmetric_composition_from_row_column_counts(row_counts, col_counts)
    ans = []
    for i, j in shape:
        while not (i < len(ans)):
            ans.append(0)
        ans[i] += 1
    return tuple(ans)


def _symmetric_composition_from_row_column_counts(row_counts, col_counts):
    n = max(len(row_counts), len(col_counts))
    row_counts = tuple(row_counts) + (n - len(row_counts)) * (0,)
    col_counts = tuple(col_counts) + (n - len(col_counts)) * (0,)
    assert sum(row_counts) == sum(col_counts)

    if sum(row_counts) == 0:
        return set()

    c = tuple(row_counts[i] + col_counts[i] - 1 for i in range(n))
    m = [i for i in range(n) if c[i] == max(c)]

    a = list(row_counts)
    for i in m:
        a[i] = 0
    for i in range(n):
        if a[i] > 0:
            a[i] -= len([j for j in m if i < j])
    a = tuple(a)

    b = list(col_counts)
    for i in m:
        b[i] = 0
    for i in range(n):
        if b[i] > 0:
            b[i] -= len([j for j in m if j < i])
    b = tuple(b)

    mu = _symmetric_composition_from_row_column_counts(a, b)

    new_rows = n * [0]
    new_cols = n * [0]
    for i, j in mu:
        if i <= j:
            new_rows[i] += 1
            new_cols[j] += 1
    for i in range(n):
        if new_rows[i] < row_counts[i] or new_cols[i] < col_counts[i]:
            for j in m:
                mu.add((i, j))
                mu.add((j, i))
    return mu


# def symmetric_composition_from_row_column_counts(row_counts, col_counts):
#     def helper(rc, cc):
#         if len(rc) == 0 or max(rc) == 0:
#             yield set()
#             return
#         m = max(rc)
#         i = [i for i, a in enumerate(rc) if a == m][-1]
#         columns = [j for j, c in enumerate(cc) if c > 0 and j >= i]
#         for subset in itertools.combinations(columns, m):
#             new_rc = rc[:i] + (0,) + rc[i + 1:]
#             new_cc = tuple((a - 1) if j in subset else a for j, a in enumerate(cc))
#             for ans in helper(new_rc, new_cc):
#                 # print(new_rc, new_cc, '\n', Tableau({k: 1 for k in ans}))
#                 for j in subset:
#                     ans |= {(i + 1, j + 1), (j + 1, i + 1)}
#                 yield ans
#     #
#     n = max(len(row_counts), len(col_counts))
#     row_counts = tuple(row_counts) + (n - len(row_counts)) * (0,)
#     col_counts = tuple(col_counts) + (n - len(col_counts)) * (0,)
#     assert sum(row_counts) == sum(col_counts)
#     #
#     answers = []
#     for bns in helper(row_counts, col_counts):
#         ans = n * [0]
#         for i, j in bns:
#             ans[i - 1] += 1
#         ans = tuple(ans)
#         while ans and ans[-1] == 0:
#             ans = ans[:-1]
#         if Partition(*sorted(ans, reverse=True)).is_symmetric():
#             answers.append(ans)
#     if len(answers) > 1:
#         raise Exception('Failed uniqueness %s, %s: %s' % (str(row_counts), str(col_counts), str(answers)))
#     if len(answers) == 0:
#         raise Exception('Failed existence %s, %s' % (str(row_counts), str(col_counts)))
#     return answers[0]


def symmetrize_strict_partition(mu, n=None):
    shape = {(i, i + j - 1) for i in range(1, len(mu) + 1) for j in range(1, mu[i - 1] + 1)}
    shape |= {(j, i) for (i, j) in shape}
    ans = []
    for i, _ in shape:
        while i - 1 >= len(ans):
            ans.append(0)
        ans[i - 1] += 1

    if n is not None:
        ans += (n - len(ans)) * (0,)
    return tuple(ans)


def skew_symmetrize_strict_partition(mu, n=None):
    shape = {(i, i + j) for i in range(1, len(mu) + 1) for j in range(1, mu[i - 1] + 1)}
    shape |= {(j, i) for (i, j) in shape}
    k = 1
    while (k, k + 1) in shape:
        shape.add((k, k))
        k += 1
    k -= 1
    ans = []
    for i, _ in shape:
        while i - 1 >= len(ans):
            ans.append(0)
        ans[i - 1] += 1
    ans = tuple(ans)
    if not is_skew_symmetric_composition(ans):
        ans = list(ans)
        if (k, k + 1) not in shape:
            ans[k - 1] -= 1
        else:
            ans[k] += 1
        ans = tuple(ans)
    assert is_skew_symmetric_composition(ans)

    if n is not None:
        ans += (n - len(ans)) * (0,)
    return ans


def symmetric_double(alpha):
    word = sorting_permutation(alpha)
    mu = sorted(alpha, reverse=True)
    shape = {(i, j) for i in range(1, len(mu) + 1) for j in range(i, i + mu[i - 1])}
    shape |= {(j, i) for (i, j) in shape}
    t = Tableau({a: 1 for a in shape})
    mu = t.partition().parts
    for i in reversed(word):
        if i == len(mu):
            mu += [0]
        mu[i - 1], mu[i] = mu[i], mu[i - 1]
    ans = tuple(mu)
    while ans and ans[-1] == 0:
        ans = ans[:-1]
    return ans


def symmetric_half(alpha):
    standardize = sorted(
        [(i, a) for i, a in enumerate(alpha) if a != 0],
        key=lambda x: (-x[1], x[0])
    )
    ans = list(alpha)
    for diff, pair in enumerate(standardize):
        i, _ = pair
        ans[i] = max(ans[i] - diff, 0)
    ans = tuple(ans)
    while ans and ans[-1] == 0:
        ans = ans[:-1]
    return ans


def symmetric_halves(alpha, pad=False):
    def s(i):
        def f(x):
            return (x + 1) if x == i else (x - 1) if x == i + 1 else x
        return f

    word = sorting_permutation(alpha)
    mu = sorted(alpha, reverse=True)
    diagram = {(i, j) for i in range(1, 1 + len(mu)) for j in range(1, 1 + mu[i - 1])}
    for i in reversed(word):
        diagram = {(s(i)(a), s(i)(b)) for (a, b) in diagram}
    n = max([0] + [max(p) for p in diagram])
    rows, cols = n * [0], n * [0]
    for a, b in diagram:
        if a <= b:
            rows[a - 1] += 1
            cols[b - 1] += 1
    while rows and rows[-1] == 0:
        rows = rows[:-1]
    while cols and cols[-1] == 0:
        cols = cols[:-1]
    a, b = tuple(rows), tuple(cols)

    if pad:
        n = len(alpha)
        a += (n - len(a)) * (0,)
        b += (n - len(b)) * (0,)
    return a, b


def symmetric_diagram(alpha):
    return skew_symmetric_diagram(alpha)


def print_symmetric_diagram(alpha, sep='.'):
    print_skew_symmetric_diagram(alpha, sep)


def skew_symmetric_double(alpha):
    pass


def skew_symmetric_diagram(alpha):
    def s(i):
        def f(x):
            return (x + 1) if x == i else (x - 1) if x == i + 1 else x
        return f
    #
    word = sorting_permutation(alpha)
    mu = sorted(alpha, reverse=True)
    diagram = {(i, j) for i in range(1, 1 + len(mu)) for j in range(1, 1 + mu[i - 1])}
    for i in reversed(word):
        diagram = {(s(i)(a), s(i)(b)) for (a, b) in diagram}
    return diagram


def print_skew_symmetric_diagram(alpha, sep='.'):
    print(Permutation.print_diagram(skew_symmetric_diagram(alpha), sep=sep))


def skew_symmetric_halves(alpha):
    diagram = skew_symmetric_diagram(alpha)
    n = max([0] + [max(p) for p in diagram])
    rows, cols = n * [0], n * [0]
    for a, b in diagram:
        if a < b:
            rows[a - 1] += 1
            cols[b - 1] += 1
    while rows and rows[-1] == 0:
        rows = rows[:-1]
    while cols and cols[-1] == 0:
        cols = cols[:-1]
    return tuple(rows), tuple(cols)


def skew_symmetric_composition_from_row_column_counts(row_counts, col_counts):
    def helper(rc, cc):
        if len(rc) == 0 or max(rc) == 0:
            yield set()
            return
        m = max(rc)
        i = [i for i, a in enumerate(rc) if a == m][-1]
        columns = [j for j, c in enumerate(cc) if c > 0 and j > i]
        for subset in itertools.combinations(columns, m):
            new_rc = rc[:i] + (0,) + rc[i + 1:]
            new_cc = tuple((a - 1) if j in subset else a for j, a in enumerate(cc))
            for ans in helper(new_rc, new_cc):
                for j in subset:
                    ans |= {(i + 1, j + 1), (j + 1, i + 1)}
                yield ans
    #
    n = max(len(row_counts), len(col_counts))
    row_counts = tuple(row_counts) + (n - len(row_counts)) * (0,)
    col_counts = tuple(col_counts) + (n - len(col_counts)) * (0,)
    assert sum(row_counts) == sum(col_counts)
    #
    answers = []
    for cns in helper(row_counts, col_counts):
        s = list(range(1, n + 1))
        for k in range(n + 1):
            for diagonal in itertools.combinations(s, k):
                bns = cns.copy()
                for i in diagonal:
                    bns.add((i, i))
                # print(Tableau({box: 1 for box in bns}))
                ans = n * [0]
                for i, j in bns:
                    ans[i - 1] += 1
                ans = tuple(ans)
                while ans and ans[-1] == 0:
                    ans = ans[:-1]
                if is_skew_symmetric_composition(ans):
                    answers.append(ans)
    if len(answers) > 1:
        raise Exception('Failed uniqueness %s, %s: %s' % (str(row_counts), str(col_counts), str(answers)))
    if len(answers) == 0:
        raise Exception('Failed existence %s, %s' % (str(row_counts), str(col_counts)))
    return answers[0]


def weak_compositions(n, parts, allow_repeated_parts=True, reduced=False):
    assert n >= 0 and parts >= 0
    if parts == 0:
        if n == 0:
            yield ()
        return
    indices = tuple(range(1, n + parts))
    for c in itertools.combinations(indices, parts - 1):
        c = (0,) + c + (n + parts,)
        alpha = tuple(c[i] - c[i - 1] - 1 for i in range(1, len(c)))
        if not allow_repeated_parts and len({i for i in alpha if i != 0}) < len([i for i in alpha if i != 0]):
            continue
        if reduced and alpha and alpha[-1] == 0:
            continue
        yield alpha


def is_symmetric_composition(alpha):
    return Partition(*sorted(alpha, reverse=True)).is_symmetric()


def symmetric_partitions(n):
    for alpha in Partition.all(n):
        if alpha.is_symmetric():
            yield tuple(alpha.parts)


def skew_symmetric_partitions(n):
    for alpha in symmetric_partitions(n):
        if is_skew_symmetric_composition(alpha):
            yield alpha


def symmetric_weak_compositions(n, parts, reduced=False):
    seen = set()
    for mu in StrictPartition.all(n, parts):
        mu = StrictPartition.symmetric_double(mu)
        for gamma in itertools.permutations(mu):
            if gamma in seen:
                continue
            seen.add(gamma)
            if reduced:
                if len(gamma) == 0:
                    continue
                for indices in itertools.combinations(range(parts - 1),len(gamma) - 1):
                    alpha = parts * [0]
                    alpha[-1] = gamma[-1]
                    for i, a in enumerate(indices):
                        alpha[a] = gamma[i]
                    yield tuple(alpha)
            else:
                for indices in itertools.combinations(range(parts),len(gamma)):
                    alpha = parts * [0]
                    for i, a in enumerate(indices):
                        alpha[a] = gamma[i]
                    yield tuple(alpha)


def is_skew_symmetric_composition(alpha):
    if not is_symmetric_composition(alpha):
        return False
    mu = sorted(alpha, reverse=True) + [0]
    # if any(mu[i] == i and (i == 0 or mu[i - 1] > i) for i in range(len(mu))):
    #     return False
    if any(mu[i] == i and mu[i + 1] < i for i in range(len(mu) - 1)):
        return False
    if any(mu[i] == i + 1 and mu[i + 1] == i for i in range(len(mu) - 1)):
        return False
    # if any(
    #     alpha[i] > 0 == alpha[i + 1] and
    #     any(alpha[i] <= alpha[j] for j in range(i + 1, len(alpha))) and
    #     not any(0 < alpha[j] < alpha[i] for j in range(i + 1, len(alpha))
    # ) for i in range(len(alpha) - 1)):
    #     return False
    return True


def skew_symmetric_weak_compositions(n, parts, reduced=False):
    for alpha in symmetric_weak_compositions(n, parts, reduced):
        if is_skew_symmetric_composition(alpha):
            yield alpha

    # seen = set()
    # for mu in StrictPartition.all(n, parts - 1):
    #     mu = StrictPartition.skew_symmetric_double(mu)
    #     for gamma in itertools.permutations(mu):
    #         if gamma in seen:
    #             continue
    #         seen.add(gamma)
    #         if reduced:
    #             if len(gamma) == 0:
    #                 continue
    #             for indices in itertools.combinations(range(parts - 1),len(gamma) - 1):
    #                 alpha = parts * [0]
    #                 alpha[-1] = gamma[-1]
    #                 for i, a in enumerate(indices):
    #                     alpha[a] = gamma[i]
    #                 alpha = tuple(alpha)
    #                 if is_skew_symmetric_composition(alpha):
    #                     yield alpha
    #         else:
    #             for indices in itertools.combinations(range(parts),len(gamma)):
    #                 alpha = parts * [0]
    #                 for i, a in enumerate(indices):
    #                     alpha[a] = gamma[i]
    #                 alpha = tuple(alpha)
    #                 if is_skew_symmetric_composition(alpha):
    #                     yield alpha


def strict_weak_compositions(n, parts, reduced=False):
    for alpha in weak_compositions(n, parts, False, reduced):
        yield alpha


def q_power(alpha):
    mu = tuple(sorted(alpha, reverse=True))
    i = 0
    while i < len(mu) and mu[i] > i:
        i += 1
    return i


def monomial_from_composition(weak_composition):
    ans = X(0)**0
    for i, e in enumerate(weak_composition):
        ans *= X(i + 1) ** e
    return ans


def leading_monomial(weak_composition):
    return monomial_from_composition(reversed(sorted(weak_composition)))


def has_distinct_parts(mu):
    smu = [i for i in mu if i != 0]
    return len(smu) == len(set(smu))


def p_shifted_monomial(weak_composition):
    mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
    ans = X(0)**0
    for i in range(1, len(mu) + 1):
        for j in range(i + 1, mu[i - 1] + 1):
            ans *= X(i) + X(j)
    return ans


def q_shifted_monomial(weak_composition):
    mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
    ans = X(0)**0
    for i in range(1, len(mu) + 1):
        for j in range(i, mu[i - 1] + 1):
            ans *= X(i) + X(j)
    return ans
    # mu = symmetric_half(weak_composition)
    # w = Permutation.get_inv_grassmannian(*mu).shift(1)
    # w = Permutation.from_involution_shape(*mu).shift(1)
    # ans = InvSchubert.get(w).homogenize(sum(mu)) * 2**len(mu)
    # ans = restrict_variables(ans, max((1,) + mu))
    # return ans


def gp_shifted_monomial(beta=LASCOUX_BETA):
    def fun(weak_composition):
        mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
        ans = X(0)**0
        for i in range(1, len(mu) + 1):
            for j in range(i + 1, mu[i - 1] + 1):
                ans *= X(i) + X(j) + beta * X(i) * X(j)
        return ans
    return fun


def restrict_variables(ans, n):
    while any(i > n for i in ans.variables()):
        ans = ans.set(max(ans.variables()), 0)
    return ans


def gq_shifted_monomial(beta=LASCOUX_BETA):
    assert InvGrothendieck.beta == 1
    def fun(weak_composition):
        old = X(0)**0
        mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
        for i in range(1, len(mu) + 1):
            for j in range(i, mu[i - 1] + 1):
                old *= X(i) + X(j) + beta * X(i) * X(j)
        return old
        # lam = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
        # mu = symmetric_half(weak_composition)
        # w = Permutation.from_involution_shape(*mu)
        # w = Permutation.get_inv_grassmannian(*mu)
        # w = w.shift(1)
        # ans = InvGrothendieck.get(w).homogenize(sum(mu))
        # ans = ans.set(0, beta)
        # return ans
    return fun


def sorting_permutation(weak_comp):
    word = []
    n = len(weak_comp)
    weak_comp = list(weak_comp)
    for i in range(n):
        for j in range(i, 0, -1):
            if weak_comp[j] > weak_comp[j - 1]:
                word += [j]
                weak_comp[j - 1], weak_comp[j] = weak_comp[j], weak_comp[j - 1]
    return tuple(word)


def sorting_descent(weak_comp):
    for i in range(1, len(weak_comp)):
        if weak_comp[i] > weak_comp[i - 1]:
            new_comp = list(weak_comp)
            new_comp[i - 1], new_comp[i] = new_comp[i], new_comp[i - 1]
            return tuple(new_comp), i
    return weak_comp, None


def _generic_key(weak_comp, cache, name, monomial_fn, isobaric_divided_difference_fn):
    while weak_comp and weak_comp[-1] == 0:
        weak_comp = weak_comp[:-1]
    if weak_comp not in cache:
        new_comp, i = sorting_descent(weak_comp)
        if i is None:
            cache[weak_comp] = monomial_fn(weak_comp)
        else:
            f = _generic_key(new_comp, cache, name, monomial_fn, isobaric_divided_difference_fn)
            cache[weak_comp] = isobaric_divided_difference_fn(f, i)
        # if len(cache) % 100 == 0:
        #    print(' . . .', name, 'cache:', len(cache))
    return cache[weak_comp]


def key(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i)
    return _generic_key(weak_comp, KEY_POLYNOMIAL_CACHE, 'Key Polynomial', leading_monomial, isobaric_divided_difference_fn)


def atom(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, KEY_ATOM_CACHE, 'Key Atom', leading_monomial, isobaric_divided_difference_fn)


def p_key(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i)
    return _generic_key(weak_comp, PKEY_POLYNOMIAL_CACHE, 'PKey Polynomial', p_shifted_monomial, isobaric_divided_difference_fn)


def p_atom(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, PKEY_ATOM_CACHE, 'PKey Atom', p_shifted_monomial, isobaric_divided_difference_fn)


def q_key(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i)    
    return _generic_key(weak_comp, QKEY_POLYNOMIAL_CACHE, 'QKey Polynomial', q_shifted_monomial, isobaric_divided_difference_fn)


def q_atom(weak_comp):
    isobaric_divided_difference_fn = lambda f, i: f.isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, QKEY_ATOM_CACHE, 'QKey Atom', q_shifted_monomial, isobaric_divided_difference_fn)

def lascoux(weak_comp, beta=LASCOUX_BETA):
    isobaric_divided_difference_fn = lambda f, i: (f * (1 + beta * X(i + 1))).isobaric_divided_difference(i)
    return _generic_key(weak_comp, LASCOUX_POLYNOMIAL_CACHE[beta], 'Lascoux Polynomial', leading_monomial, isobaric_divided_difference_fn)


def lascoux_atom(weak_comp, beta=LASCOUX_BETA):
    LASCOUX_ATOM_CACHE[beta] = LASCOUX_ATOM_CACHE.get(beta, None) or {}
    isobaric_divided_difference_fn = lambda f, i: (f * (1 + beta * X(i + 1))).isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, LASCOUX_ATOM_CACHE[beta], 'Lascoux Atom', leading_monomial, isobaric_divided_difference_fn)


def p_lascoux(weak_comp, beta=LASCOUX_BETA):
    PLASCOUX_POLYNOMIAL_CACHE[beta] = PLASCOUX_POLYNOMIAL_CACHE.get(beta, None) or {}
    isobaric_divided_difference_fn = lambda f, i: (f * (1 + beta * X(i + 1))).isobaric_divided_difference(i)
    return _generic_key(weak_comp, PLASCOUX_POLYNOMIAL_CACHE[beta], 'PLascoux Polynomial', gp_shifted_monomial(beta), isobaric_divided_difference_fn)


def p_lascoux_atom(weak_comp, beta=LASCOUX_BETA):
    PLASCOUX_ATOM_CACHE[beta] = PLASCOUX_ATOM_CACHE.get(beta, None) or {}
    isobaric_divided_difference_fn = lambda f, i: (f * (1 + beta * X(i + 1))).isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, PLASCOUX_ATOM_CACHE[beta], 'PLascoux Atom', gp_shifted_monomial(beta), isobaric_divided_difference_fn)


def q_lascoux(weak_comp, beta=LASCOUX_BETA):
    QLASCOUX_POLYNOMIAL_CACHE[beta] = QLASCOUX_POLYNOMIAL_CACHE.get(beta, None) or {}
    term = lambda i: X(i) * (1 + beta * X(i + 1))
    isobaric_divided_difference_fn = lambda f, i: (f * term(i)).divided_difference(i)
    ans = _generic_key(weak_comp, QLASCOUX_POLYNOMIAL_CACHE[beta], 'QKey Polynomial', gq_shifted_monomial(beta), isobaric_divided_difference_fn)
    # n = len(weak_comp)
    # if n > 0 and min(weak_comp) > 0:
    #    ans = ans.set(n + 1, 0)
    return ans


def q_lascoux_atom(weak_comp, beta=LASCOUX_BETA):
    QLASCOUX_ATOM_CACHE[beta] = QLASCOUX_ATOM_CACHE.get(beta, None) or {}
    isobaric_divided_difference_fn = lambda f, i: (f * (1 + beta * X(i + 1))).isobaric_divided_difference(i) - f
    return _generic_key(weak_comp, QLASCOUX_ATOM_CACHE[beta], 'QLascoux Atom', gq_shifted_monomial(beta), isobaric_divided_difference_fn)


def composition_length(alpha):
    while alpha and alpha[-1] == 0:
        alpha = alpha[:-1]
    return len(alpha)

def tuplize(g):
    beta = []
    for i in g:
        while i > len(beta):
            beta += [0]
        beta[i - 1] = g[i]
    beta = tuple(beta)
    return beta


def get_exponents(kappa):
    # ans = sorted([tuplize(g) for g in kappa], key=lambda x: (sum(x), x))
    # return [a for a in ans if sum(a) == sum(ans[0])]
    return sorted([tuplize(g) for g in kappa])


def dict_from_tuple(beta):
    mon = X(0)**0
    for i, a in enumerate(beta):
        mon *= X(i + 1)**a
    return max(mon)


def decompose_into_compositions(kappa):
    def xbeta(beta):
        mon = X(0)**0
        for i, a in enumerate(beta):
            mon *= X(i + 1)**a
        return mon

    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * xbeta(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def decompose_into_lascoux(kappa):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * lascoux(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def decompose_lascoux(kappa):
    dec = decompose_into_lascoux(kappa)
    assert len(dec) == 1
    assert set(dec.values()) == {1}
    return list(dec)[0]


def decompose_into_keys(kappa):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * key(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def decompose_key(kappa):
    dec = decompose_into_keys(kappa)
    assert len(dec) == 1
    assert set(dec.values()) == {1}
    return list(dec)[0]


def decompose_into_atoms(kappa):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * atom(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def maximal_weakly_increasing_factors(w):
    factors = [[]]
    for a in w:
        if len(factors[-1]) == 0 or factors[-1][-1] <= a:
            factors[-1].append(a)
        else:
            factors.append([a])
    return tuple(tuple(a) for a in factors)


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


def colform(w):
    return Partition(*sorted([len(a) for a in maximal_decreasing_factors(w)], reverse=True)).tuple()


def is_key_word(p, w):
    return p == w.rsk_insert(w)[0] and colform(w) == p.partition().tuple()


def key_words(p):
    for w in knuth_class(p):
        if colform(w) == p.partition().transpose().tuple():
            yield w


def key_tableau(alpha):
    from permutations import Permutation
    w = Permutation()
    for i in sorting_permutation(alpha):
        w *= Permutation.s_i(i)
    word = []
    for part in Partition(*sorted(alpha, reverse=True)).transpose().parts:
        word += sorted([w(i) for i in range(part, 0, -1)], reverse=True)
    return rsk_insert(word)[0]


def weak_compatible_sequences(seq, i_min=1):
    if len(seq) == 0:
        yield (), X(0)**0
    else:
        a, seq = seq[0], seq[1:]
        for i in range(i_min, abs(a) + 1):
            j_min = (i + 1) if (seq and a <= seq[0]) else i
            for p, q in weak_compatible_sequences(seq, j_min):
                yield (a,) + p, X(i) * q


def compatible_sequences(seq, i_min=1, flag=None):
    def phi(a):
        if flag is None:
            return a
        if type(flag) in [list, tuple, dict]:
            return flag[a]
        return flag(a)

    if len(seq) == 0:
        yield (), X(0)**0
    else:
        a, seq = seq[0], seq[1:]
        for i in range(i_min, phi(abs(a)) + 1):
            j_min = (i + 1) if (seq and a < seq[0]) else i
            for p, q in compatible_sequences(seq, j_min, flag):
                yield (a,) + p, X(i) * q


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


def _equivalence_class(p, cache, toggle_fn):
    if p in cache:
        return cache[p]
    seen = set()
    add = {p}
    ans = []
    while add:
        nextadd = set()
        for v in add:
            if v not in seen:
                seen.add(v)
                for u in toggle_fn(v):
                    if u not in seen:
                        nextadd.add(u)
                ans.append(v)
        add = nextadd
    cache[p] = ans
    return ans


def shifted_knuth_class(p):
    def toggle(w):
        for i in range(len(w) - 2):
            a, c, b = w[i:i + 3]
            ma, mc, mb = MarkedNumber(a), MarkedNumber(c), MarkedNumber(b)
            if ma.ceil() <= mb <= mc.floor() or mc.ceil() <= mb <= ma.floor():
                yield w[:i] + (c, a, b) + w[i + 3:]

            b, c, a = w[i:i + 3]
            ma, mc, mb = MarkedNumber(a), MarkedNumber(c), MarkedNumber(b)
            if (ma <= mb.floor() and mb.ceil() <= mc) or (mc <= mb.floor() and mb.ceil() <= ma):
                yield w[:i] + (b, a, c) + w[i + 3:]

        if len(w) >= 2:
            a, b = w[:2]
            if abs(a) == abs(b):
                yield (a, -b) + w[2:]
            elif (a < 0 and b > 0) or (a > 0 and b < 0):
                yield (-b, -a) + w[2:]  # primed SW
                # yield (b, a) + w[2:]  # unprimed SW
            else:
                yield (b, a) + w[2:]

        if len(w) >= 1:
            yield (-w[0],) + w[1:]

        # for i in range(len(w) - 2):
        #     if w[i] <= w[i + 2] < w[i + 1] or w[i + 1] <= w[i + 2] < w[i]:
        #         yield w[:i] + (w[i + 1], w[i], w[i + 2]) + w[i + 3:]
        #     if w[i + 1] < w[i] <= w[i + 2] or w[i + 2] < w[i] <= w[i + 1]:
        #         yield w[:i] + (w[i], w[i + 2], w[i + 1]) + w[i + 3:]
        # if len(w) >= 2:
        #     yield (w[1], w[0]) + w[2:]

    p = p.row_reading_word() if type(p) == Tableau else p
    return _equivalence_class(p, SHIFTED_REDUCED_TABLEAU_CACHE, toggle)


def coxeter_knuth_class(p):
    def toggle(w):
        for i in range(len(w) - 2):
            if w[i] == w[i + 2]:
                yield w[:i] + (w[i + 1], w[i], w[i + 1]) + w[i + 3:]
            if w[i] < w[i + 2] < w[i + 1] or w[i + 1] < w[i + 2] < w[i]:
                yield w[:i] + (w[i + 1], w[i], w[i + 2]) + w[i + 3:]
            if w[i + 2] < w[i] < w[i + 1] or w[i + 1] < w[i] < w[i + 2]:
                yield w[:i] + (w[i], w[i + 2], w[i + 1]) + w[i + 3:]

    p = p.row_reading_word() if type(p) == Tableau else p
    return _equivalence_class(p, REDUCED_TABLEAU_CACHE, toggle)


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

    p = p.row_reading_word() if type(p) == Tableau else p
    return _equivalence_class(p, SP_REDUCED_TABLEAU_CACHE, toggle)


def o_knuth_class(p):
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

    p = p.row_reading_word() if type(p) == Tableau else p
    return _equivalence_class(p, O_REDUCED_TABLEAU_CACHE, toggle)
