from schubert import (
    X, 
    Schubert, DoubleSchubert,
    Grothendieck, DoubleGrothendieck,
    GrothendieckB, DoubleGrothendieckB, 
    GrothendieckC, DoubleGrothendieckC,
    GrothendieckD, DoubleGrothendieckD,
)
from permutations import Permutation
from signed import SignedPermutation
from stable.utils import (
    G,
    GP,
    GP_expansion, 
    GQ, 
    GQ_expansion,
    SymmetricPolynomial,
    beta
)
from vectors import Vector
import itertools


def expand_reflection_chain(start, chain, length):
    return GrothendieckC.expand_reflection_chain(start, chain, length)


def test_general_a_grothendieck_transitions(n=4):
    for w in Permutation.all(n):
        for r in range(1, n + 1):
            chain = []
            chain += [(Permutation.t_ij(i, r), 1) for i in range(r - 1, 0, -1)]
            chain += [(Permutation.t_ij(r, j), -1) for j in range(n + 1, r, -1)]
            
            ans = expand_reflection_chain(w, chain, lambda x: x.length())
            expected = sum([coeff * Grothendieck.get(z).set_variable(0, -1) for (z, coeff) in ans.dictionary.items()])
            actual = Grothendieck.get(w).set_variable(0, -1) * (1 - X(r))

            print('n =', n, 'w =', w, 'r =', r)
            if expected != actual:
                print(chain)
                print(ans)
                print(expected)
                print(actual)
                print()
            assert expected == actual


def test_a_grothendieck_transitions(n=4):
    for w in Permutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if w(1) > 1 or len(r) == 0:
            continue
        r = max(r)
        s = max([(w(i), i) for i in range(r + 1, n + 1) if w(i) < w(r)])[1]
        v = w * Permutation.t_ij(r, s)

        chain = [(Permutation.t_ij(i, r), X(0)) for i in range(1, r)]
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        expected = Grothendieck.get(v) - (X(0) * X(r) + 1) * sum([coeff * Grothendieck.get(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = Grothendieck.get(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(chain)
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_double_a_grothendieck_transitions(n=4):
    for w in Permutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if w(1) > 1 or len(r) == 0:
            continue
        r = max(r)
        t, s = max([(w(i), i) for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = w * Permutation.t_ij(r, s)

        chain = [(Permutation.t_ij(i, r), X(0)) for i in range(1, r)]
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        expected = DoubleGrothendieck.get(v) - (X(0) * (X(r) + X(-t) + X(0)*X(r)*X(-t)) + 1) * sum([coeff * DoubleGrothendieck.get(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = DoubleGrothendieck.get(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(chain)
            print()
            print(ans)
            print()
            print('    :', expected)
            print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def subsets(coll):
    for k in range(1 + len(coll)):
        for sub in itertools.combinations(coll, k):
            yield sub


def signed(coll):
    for sub in subsets(coll):
        n = len(sub)
        for v in range(2**n):
            signs = []
            for _ in range(n):
                signs.append(1 if v % 2 == 0 else -1)
                v = v // 2
            yield tuple(zip(sub, signs))


def test_b_grothendieck_transitions(n=3):
    GROTH = GrothendieckB.get

    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain =  []
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        expected = GROTH(v) - (X(0) * X(r) + 1) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = GROTH(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print()
            print('base:', -X(0)**-1 * GROTH(v))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', -X(0)**-1 * ans[u] * GROTH(u) * (X(r) - 1))
                print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def test_double_b_grothendieck_transitions(n=3):
    GROTH = DoubleGrothendieckB.get

    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        adj_r = 2 * r - 1
        adj_s = 2 * abs(v(r))

        chain =  []
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), 1+X(0)*X(adj_s), X(0))]
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), 1, X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), 1, X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), 1, X(0)) for i in range(1, r)]

        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        if v(r) < 0:
            expected = X(0) * (X(adj_r) - X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (1 + X(0)*X(adj_s)) * ((1 + X(0)*X(adj_s)) * (X(0) * GROTH(w) + GROTH(v)) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]))
        else:
            expected = X(0) * (X(adj_r) + X(adj_s) + X(0)*X(adj_r)*X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (1 + X(0)*X(adj_s)) * (X(0) * GROTH(w) + GROTH(v)) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])

        print('w =', w, 'r =', r, 's =', s, 'v =', v, 'v(r) =', v(r))
        if expected != actual:
            print(ans)
            print()
            print('want:', actual.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print(' got:', expected.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print('diff:', (expected-actual).set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            input('\n\n\n\n')
        assert expected == actual


def test_c_grothendieck_transitions(n=3):
    GROTH = GrothendieckC.get

    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain =  []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        expected = GROTH(v) - (X(0) * X(r) + 1) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = GROTH(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print()
            print('base:', -X(0)**-1 * GROTH(v))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', -X(0)**-1 * ans[u] * GROTH(u) * (X(r) - 1))
                print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def test_double_c_grothendieck_transitions(n=3):
    GROTH = DoubleGrothendieckC.get

    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        adj_r = 2 * r - 1
        adj_s = 2 * abs(v(r))

        chain =  []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = expand_reflection_chain(v, chain, lambda x: x.length())

        if v(r) < 0:
            expected = X(0) * (X(adj_r) - X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (1 + X(0)*X(adj_s)) * (X(0) * GROTH(w) + GROTH(v) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]))
        else:
            expected = X(0) * (X(adj_r) + X(adj_s) + X(0)*X(adj_r)*X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (X(0) * GROTH(w) + GROTH(v) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]))

        print('w =', w, 'r =', r, 's =', s, 'v =', v, 'v(r) =', v(r))
        if expected != actual:
            print(ans)
            print()
            print(len(v), ':', v, ':', GROTH(v).set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', (ans[u] * GROTH(u)).set_variable(0,-1).tostring(Schubert.double_lettering()))
                print()
            print('want:', actual.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print(' got:', expected.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print('diff:', (expected-actual).set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print()
            print()
            print(sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]).set_variable(0,-1).tostring(Schubert.double_lettering()))
            input('\n\n\n\n')
        assert expected == actual


def test_d_grothendieck_transitions(n=4):
    GROTH = GrothendieckD.get

    for w in SignedPermutation.all(n, dtype=True):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain =  []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]

        ans = expand_reflection_chain(v, chain, lambda x: x.dlength())

        expected = GROTH(v) - (X(0) * X(r) + 1) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = GROTH(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print()
            print('base:', -X(0)**-1 * GROTH(v))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', -X(0)**-1 * ans[u] * GROTH(u) * -(X(0) * X(r) - 1))
                print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def test_double_d_grothendieck_transitions(n=3):
    GROTH = DoubleGrothendieckD.get

    for w in SignedPermutation.all(n, dtype=True):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        adj_r = 2 * r - 1
        adj_s = 2 * abs(v(r))

        chain =  []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]

        ans = expand_reflection_chain(v, chain, lambda x: x.dlength())

        if v(r) < 0:
            expected = X(0) * (X(adj_r) - X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (1 + X(0)*X(adj_s)) * (X(0) * GROTH(w) + GROTH(v) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]))
        else:
            expected = X(0) * (X(adj_r) + X(adj_s) + X(0)*X(adj_r)*X(adj_s)) * sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()])
            actual = (X(0) * GROTH(w) + GROTH(v) - sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]))

        print('w =', w, 'r =', r, 's =', s, 'v =', v, 'v(r) =', v(r))
        if expected != actual:
            print(ans)
            print()
            print(len(v), ':', v, ':', GROTH(v).set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', (ans[u] * GROTH(u)).set_variable(0,-1).tostring(Schubert.double_lettering()))
                print()
            print('want:', actual.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print(' got:', expected.set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print('diff:', (expected-actual).set_variable(0,-1).tostring(Schubert.double_lettering()))
            print()
            print()
            print()
            print(sum([coeff * GROTH(z) for (z, coeff) in ans.dictionary.items()]).set_variable(0,-1).tostring(Schubert.double_lettering()))
            input('\n\n\n\n')
        assert expected == actual


def test_a_symmetric_transitions(n=3, numvars=2):
    for w in Permutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if w(1) > 1 or len(r) == 0:
            continue
        r = max(r)
        s = max([(w(i), i) for i in range(r + 1, n + 1) if w(i) < w(r)])[1]
        v = w * Permutation.t_ij(r, s)

        chain = [(Permutation.t_ij(i, r), beta) for i in range(1, r)]
        ans = Vector({v: 1}) - expand_reflection_chain(v, chain, lambda x: x.length())
        ans = Vector({u: -beta**-1 * c for (u, c) in ans.dictionary.items()})

        expected = sum([coeff * G(numvars, z) for (z, coeff) in ans.dictionary.items()])
        actual = G(numvars, w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(chain)
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_b_symmetric_transitions(n=3, numvars=2):
    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = []
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = Vector({v: 1}) - expand_reflection_chain(v, chain, lambda x: x.length())
        ans *= -X(0)**-1

        expected = sum([coeff * GrothendieckB.symmetric(numvars, z) for (z, coeff) in ans.dictionary.items()])
        actual = GrothendieckB.symmetric(numvars, w)
        
        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_c_symmetric_transitions(n=3, numvars=2):
    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_s(r, r, n + 1), X(0))]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = Vector({v: 1}) - expand_reflection_chain(v, chain, lambda x: x.length())
        ans *= -X(0)**-1

        expected = sum([coeff * GrothendieckC.symmetric(numvars, z) for (z, coeff) in ans.dictionary.items()])
        actual = GrothendieckC.symmetric(numvars, w)
        
        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_d_symmetric_transitions(n=4, numvars=2):
    for w in SignedPermutation.all(n, dtype=True):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain =  []
        chain += [(SignedPermutation.reflection_s(i, r, n + 1), X(0)) for i in range(n + 1, 0, -1) if i != r]
        chain += [(SignedPermutation.reflection_t(i, r, n + 1), X(0)) for i in range(1, r)]
        
        ans = Vector({v: 1}) - expand_reflection_chain(v, chain, lambda x: x.dlength())
        ans *= -X(0)**-1

        expected = sum([coeff * GrothendieckD.symmetric(numvars, z) for (z, coeff) in ans.dictionary.items()])
        actual = GrothendieckD.symmetric(numvars, w)
        
        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print(ans)
            print(expected)
            print(actual)
            print()
        assert expected == actual


def test_b(rank=3, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckB.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            if exp != 0:
                print('  ', exp)
                print('  ', GrothendieckB.symmetric_simple(w))
                print()
            assert all(v.is_positive() for v in exp.values())


def test_c(rank=3, numvars=2):
    for w in SignedPermutation.all(rank):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckC.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GQ_expansion(f)
            if exp != 0:
                print('  ', exp)
                print('  ', GrothendieckC.symmetric_simple(w))
                print()
            assert all(v.is_positive() for v in exp.values())


def test_d(rank=4, numvars=2):
    for w in SignedPermutation.all(rank, dtype=True):
        print('w =', w)
        for n in range(numvars + 1):
            fp = GrothendieckD.symmetric(n, w)
            f = SymmetricPolynomial.from_polynomial(fp, n)
            exp = GP_expansion(f)
            if exp != 0:
                print('  ', exp)
                print('  ', GrothendieckD.symmetric_simple(w))
                print()
            assert all(v.is_positive() for v in exp.values())


def test_b_grassmannian(rank=3, numvars=2):
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


def test_c_grassmannian(rank=3, numvars=2):
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


def test_double_grothendieck(n=4):
    for w in Permutation.all(n):
        f = DoubleGrothendieck.get(w)
        g = DoubleSchubert.get(w)
        assert f.set_variable(0, 0).negate_vars(range(-1, -n, -1)) == g

