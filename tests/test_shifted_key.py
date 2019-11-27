from schubert import *


def monomial(weak_composition):
    ans = X(0)**0
    for i, e in enumerate(weak_composition):
        ans *= X(i + 1) ** e
    return ans


def leading_monomial(weak_composition):
    return monomial(reversed(sorted(weak_composition)))


def has_distinct_parts(mu):
    if any(mu[i] == mu[i + 1] for i in range(len(mu) - 1)):
        raise Exception('Weak composition has repeated parts')


def shifted_monomial(weak_composition):
    mu = tuple(sorted((i for i in weak_composition if i != 0), reverse=True))
    has_distinct_parts(mu)
    ans = X(0)**0
    for i in range(1, 1 + len(mu)):
        for j in range(1, mu[i - 1] + 1):
            ans *= (X(i) + X(i + j - 1)) if j > 1 else X(i)
    for i, e in enumerate(mu):
        ans *= X(i + 1) ** (e - i)
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


def test_sorting_permutation():
    assert sorting_permutation((1, 0, 2, 1)) == (2, 1, 3)
    assert sorting_permutation((0, 0, 0, 0)) == tuple()
    assert sorting_permutation((1, 2, 3)) == (1, 2, 1)


def key(weak_composition):
    ans = leading_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i)
    return ans


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


def shifted_key(weak_composition):
    ans = shifted_monomial(weak_composition)
    word = sorting_permutation(weak_composition)
    for i in reversed(word):
        ans = ans.isobaric_divided_difference(i)
    return ans


def schur(partition):
    assert all(partition[i] >= partition[i + 1] for i in range(len(partition) - 1))
    w = Permutation.get_grassmannian(*partition)
    n = len(partition)
    w = w.shift(w.rank)
    return Schubert.get(w).truncate(n)


def test_schur():
    assert schur((3, 1, 1)) == key((1, 1, 3, 0, 0))



def schurp(partition):
    assert all(partition[i] > partition[i + 1] for i in range(len(partition) - 1))
    w = Permutation.get_inv_grassmannian(*partition)
    n = len(partition)
    w = w.shift(w.rank)
    return InvSchubert.get(w)
