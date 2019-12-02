from key import (
    weak_compositions,
    monomial,
    sorting_permutation,
    strict_weak_compositions,
    key, atom,
    p_key, p_atom,
    q_key, q_atom,
    has_distinct_parts,
)
from schubert import Schubert, InvSchubert, FPFSchubert, X
from permutations import Permutation


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
    for n in range(5):
        for k in range(7):
            for alpha in weak_compositions(n, k):
                kappa = p_atom(alpha)
                if not kappa.is_positive():
                    assert not has_distinct_parts(alpha)
                elif not has_distinct_parts(alpha):
                    print(alpha, kappa)
                    print()
                assert kappa.is_not_laurent_polynomial()
    print('success')
#    assert False


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


def test_inv_schubert(n=4):
    i = set(Permutation.involutions(n))
    s = {w: InvSchubert.get(w) * 2**w.number_two_cycles() for w in i}
    m = max({w.involution_length() for w in i})
    e = {alpha: q_key(alpha) for l in range(m + 1) for alpha in strict_weak_compositions(l, n)}
    v = set(e.values())
    x = {w for w in s if s[w] in v}
    return i, s, x, e


def test_fpf_schubert(n=4):
    i = set(Permutation.fpf_involutions(n))
    s = {w: FPFSchubert.get(w) for w in i}
    m = max({w.fpf_involution_length() for w in i})
    e = {alpha: p_key(alpha) for l in range(m + 1) for alpha in strict_weak_compositions(l, n)}
    v = set(e.values())
    x = {w for w in s if s[w] in v}
    return i, s, x, e
