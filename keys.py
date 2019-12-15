from schubert import X
from partitions import Partition
from tableaux import Tableau
import itertools


KEY_POLYNOMIAL_CACHE = {}
PKEY_POLYNOMIAL_CACHE = {}
QKEY_POLYNOMIAL_CACHE = {}

KEY_ATOM_CACHE = {}
PKEY_ATOM_CACHE = {}
QKEY_ATOM_CACHE = {}


def symmetric_composition_from_row_column_counts(row_counts, col_counts):
    def helper(rc, cc):
        if len(rc) == 0 or max(rc) == 0:
            yield set()
            return
        m = max(rc)
        i = [i for i, a in enumerate(rc) if a == m][-1]
        columns = [j for j, c in enumerate(cc) if c > 0 and j >= i]
        for subset in itertools.combinations(columns, m):
            new_rc = rc[:i] + (0,) + rc[i + 1:]
            new_cc = tuple((a - 1) if j in subset else a for j, a in enumerate(cc))
            for ans in helper(new_rc, new_cc):
                # print(new_rc, new_cc, '\n', Tableau({k: 1 for k in ans}))
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
    for bns in helper(row_counts, col_counts):
        ans = n * [0]
        for i, j in bns:
            ans[i - 1] += 1
        ans = tuple(ans)
        while ans and ans[-1] == 0:
            ans = ans[:-1]
        if Partition(*sorted(ans, reverse=True)).is_symmetric():
            answers.append(ans)
    if len(answers) > 1:
        raise Exception('Failed uniqueness %s, %s: %s' % (str(row_counts), str(col_counts), str(answers)))
    if len(answers) == 0:
        raise Exception('Failed existence %s, %s' % (str(row_counts), str(col_counts)))
    return answers[0]


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


def symmetric_halves(alpha):
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
    return tuple(rows), tuple(cols)


def skew_symmetric_double(alpha):
    pass


def skew_symmetric_halves(alpha):
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


def symmetric_weak_compositions(n, parts, reduced=False):
    for alpha in weak_compositions(n, parts, True, reduced):
        if Partition(*sorted(alpha, reverse=True)).is_symmetric():
            yield alpha


def is_skew_symmetric_composition(alpha):
    mu = sorted(alpha, reverse=True) + [0]
    if not Partition(*mu).is_symmetric():
        return False
    mu += [0]
    if any(mu[i] == i and mu[i + 1] < i for i in range(len(mu) - 1)):
        return False
    if any(mu[i] == i + 1 and mu[i + 1] == i for i in range(len(mu) - 1)):
        return False
    # if alpha and alpha[0] == 0:
    #     a = [a for a in alpha if a > 0]
    #     if a and a[0] == min(a):
    #         return False
    return True


def skew_symmetric_weak_compositions(n, parts, reduced=False):
    for alpha in symmetric_weak_compositions(n, parts, reduced):
        if is_skew_symmetric_composition(alpha):
            yield alpha


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


def _generic_key(weak_comp, cache, name, atomic, monomial_fn):
    while weak_comp and weak_comp[-1] == 0:
        weak_comp = weak_comp[:-1]
    if weak_comp not in cache:
        new_comp, i = sorting_descent(weak_comp)
        if i is None:
            cache[weak_comp] = monomial_fn(weak_comp)
        else:
            f = _generic_key(new_comp, cache, name, atomic, monomial_fn)
            cache[weak_comp] = f.isobaric_divided_difference(i) - (f if atomic else 0)
        if len(cache) % 100 == 0:
            print(' . . .', name, 'cache:', len(cache))
    return cache[weak_comp]


def key(weak_comp):
    return _generic_key(weak_comp, KEY_POLYNOMIAL_CACHE, 'Key Polynomial', False, leading_monomial)


def atom(weak_comp):
    return _generic_key(weak_comp, KEY_ATOM_CACHE, 'Key Atom', True, leading_monomial)


def p_key(weak_comp):
    return _generic_key(weak_comp, PKEY_POLYNOMIAL_CACHE, 'PKey Polynomial', False, p_shifted_monomial)


def p_atom(weak_comp):
    return _generic_key(weak_comp, PKEY_ATOM_CACHE, 'PKey Atom', True, p_shifted_monomial)


def q_key(weak_comp):
    return _generic_key(weak_comp, QKEY_POLYNOMIAL_CACHE, 'QKey Polynomial', False, q_shifted_monomial)


def q_atom(weak_comp):
    return _generic_key(weak_comp, QKEY_ATOM_CACHE, 'QKey Atom', True, q_shifted_monomial)


def tuplize(g):
    beta = []
    for i in g:
        while i > len(beta):
            beta += [0]
        beta[i - 1] = g[i]
    beta = tuple(beta)
    return beta


def get_exponents(kappa):
    return sorted([tuplize(g) for g in kappa])


def dict_from_tuple(beta):
    mon = X(0)**0
    for i, a in enumerate(beta):
        mon *= X(i + 1)**a
    return max(mon)


def decompose_into_keys(kappa):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * key(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}


def decompose_into_atoms(kappa):
    ans = {}
    while kappa != 0:
        betas = sorted(get_exponents(kappa), key=lambda x: (len(x), x))
        beta = betas[0]
        coeff = kappa[dict_from_tuple(beta)]
        kappa = kappa - coeff * atom(beta)
        ans[beta] = ans.get(beta, 0) + coeff
    return {k: v for k, v in ans.items() if v}
