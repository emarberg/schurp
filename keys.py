from schubert import X
import itertools


KEY_POLYNOMIAL_CACHE = {}
PKEY_POLYNOMIAL_CACHE = {}
QKEY_POLYNOMIAL_CACHE = {}

KEY_ATOM_CACHE = {}
PKEY_ATOM_CACHE = {}
QKEY_ATOM_CACHE = {}


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


def strict_weak_compositions(n, parts):
    for alpha in weak_compositions(n, parts, False):
        yield alpha


def monomial(weak_composition):
    ans = X(0)**0
    for i, e in enumerate(weak_composition):
        ans *= X(i + 1) ** e
    return ans


def leading_monomial(weak_composition):
    return monomial(reversed(sorted(weak_composition)))


def has_distinct_parts(mu):
    smu = [i for i in mu if i != 0]
    return len(smu) == len(set(smu))


def p_shifted_monomial(weak_composition):
    mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
    # has_distinct_parts(mu)
    ans = X(0)**0
    for i in range(1, len(mu) + 1):
        for j in range(1, mu[i - 1] + 1):
            ans *= X(i) + X(i + j)
    return ans


def q_shifted_monomial(weak_composition):
    mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
    # has_distinct_parts(mu)
    ans = X(0)**0
    for i in range(1, len(mu) + 1):
        for j in range(1, mu[i - 1] + 1):
            ans *= X(i) + X(i + j - 1)
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
