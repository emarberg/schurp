from schubert import GrothendieckB, GrothendieckC, GrothendieckD
from permutations import Permutation
from signed import SignedPermutation
from stable.utils import (
    G,
    GP,
    GP_expansion, 
    GQ, 
    GQ_expansion,
    SymmetricPolynomial,
)
from vectors import Vector
import itertools


def expand_reflection_chain(start, chain, length):
    ans = Vector()
    queue = [(1, start, chain)]
    while queue:
        (sgn, z, c) = queue[0]
        queue = queue[1:]
        if len(c) == 0:
            ans += Vector({z: sgn})
        else:
            queue.append((sgn, z, c[1:]))
            t = c[0]
            zt = z * t
            if length(zt) == length(z) + 1:
                queue.append((-sgn, zt, c[1:]))
    ans = Vector({start: 1}) - ans
    return ans


def test_a_symmetric_transitions(n=4, numvars=2):
    for w in Permutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if w(1) > 1 or len(r) == 0:
            continue
        r = max(r)
        s = max([(w(i), i) for i in range(r + 1, n + 1) if w(i) < w(r)])[1]
        v = w * Permutation.t_ij(r, s)

        chain = [Permutation.t_ij(i, r) for i in range(1, r)]
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        expected = sum([coeff * G(numvars, z).set_variable(0, -1) for (z, coeff) in ans.dictionary.items()])
        actual = G(numvars, w).set_variable(0, -1)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(chain)
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_b_symmetric_transitions(n=4, numvars=2):
    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = [SignedPermutation.reflection_t(i, r, n + 1) for i in range(r - 1, 0, -1)]
        chain += [SignedPermutation.reflection_s(r, r, n + 1), SignedPermutation.reflection_s(r, r, n + 1)]
        chain += [SignedPermutation.reflection_s(i, r, n + 1) for i in range(n + 1, 0, -1) if i != r]
        ans = expand_reflection_chain(v, chain, lambda x: x.length())
        expected = sum([coeff * GrothendieckB.symmetric(numvars, z).set(0, -1) for (z, coeff) in ans.dictionary.items()])
        actual = GrothendieckB.symmetric(numvars, w).set(0, -1)

        print('w =', w, 'r =', r, 's =', s, 'v =', v, chain)
        if expected != actual:
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_c_symmetric_transitions(n=4, numvars=2):
    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = [SignedPermutation.reflection_t(i, r, n + 1) for i in range(r - 1, 0, -1)]
        chain += [SignedPermutation.reflection_s(r, r, n + 1)]
        chain += [SignedPermutation.reflection_s(i, r, n + 1) for i in range(n + 1, 0, -1) if i != r]

        trials = []
        for ch in itertools.permutations(chain):
            ans = expand_reflection_chain(v, chain, lambda x: x.length())
            expected = sum([coeff * GrothendieckC.symmetric(numvars, z).set(0, -1) for (z, coeff) in ans.dictionary.items()])
            actual = GrothendieckC.symmetric(numvars, w).set(0, -1)
            if expected == actual:
                trials += [ch]


        print('w =', w, 'r =', r, 's =', s, 'v =', v, len(trials))
        print(' ', chain)
        #for ch in trials:
        #    print('  ', ch)
        # if expected != actual:
        #     print(ans)
        #     print(expected)
        #     print(actual)
        #     print()
        # assert expected == actual
        assert len(trials) > 0


def test_d_symmetric_transitions(n=4, numvars=2):
    for w in SignedPermutation.all(n, dtype=True):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = [SignedPermutation.reflection_t(i, r, n + 1) for i in range(r - 1, 0, -1)]
        chain += [SignedPermutation.reflection_s(i, r, n + 1) for i in range(n + 1, 0, -1) if i != r]

        trials = []
        for ch in itertools.permutations(chain):
            ans = expand_reflection_chain(v, ch, lambda x: x.dlength())
            expected = sum([coeff * GrothendieckD.symmetric(numvars, z).set(0, -1) for (z, coeff) in ans.dictionary.items()])
            actual = GrothendieckD.symmetric(numvars, w).set(0, -1)
            if expected == actual:
                trials += [ch]

        print('w =', w, 'r =', r, 's =', s, 'v =', v, len(trials))
        #for ch in trials:
        #    print('  ', ch)
        # if not any(trials):
        #     print(ans)
        #     print(expected)
        #     print(actual)
        #     print()
        assert len(trials) > 0


def test_b(rank=4, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckB.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_c(rank=4, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckC.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GQ_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_d(rank=4, numvars=2):
    for w in SignedPermutation.all(rank, dtype=True):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckD.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            print('  ', exp)
            assert all(v.is_positive() for v in exp.values())


def test_b_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_bc(rank):
        for n in range(numvars + 1):
            f = GrothendieckB.symmetric(n, w)
            g = GP(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g


def test_c_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_bc(rank):
        for n in range(numvars + 1):
            f = GrothendieckC.symmetric(n, w)
            g = GQ(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g


def test_d_grassmannian(rank=4, numvars=2):
    for w, mu in SignedPermutation.get_grassmannians_d(rank):
        for n in range(numvars + 1):
            f = GrothendieckD.symmetric(n, w)
            g = GP(n, mu).polynomial()
            print('w =', w, 'mu =', mu, 'n =', n)
            if f != g:
                print()
                print(f)
                print(g)
                print()
            assert f == g
