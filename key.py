from schubert import X
import itertools


def weak_compositions(n, parts, allow_repeated_parts=True):
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


def key(weak_composition):
    ans = leading_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i)
    return ans


def atom(weak_composition):
    ans = leading_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i) - ans
    return ans


def p_key(weak_composition):
    ans = p_shifted_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.divided_difference(i)
    return ans


def p_atom(weak_composition):
    ans = p_shifted_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i) - ans
    return ans


def q_key(weak_composition):
    ans = q_shifted_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.divided_difference(i)
    return ans


def q_atom(weak_composition):
    ans = q_shifted_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i) - ans
    return ans
