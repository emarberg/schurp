from schubert import X, Grothendieck, GrothendieckB, GrothendieckC, GrothendieckD
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
    newchain = [(c, -1) if type(c) in [Permutation, SignedPermutation] else c for c in chain]
    return GrothendieckC.expand_reflection_chain(start, newchain, length)


def test_a_grothendieck_transitions(n=4):
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

        expected = GrothendieckB.get(v) - (X(0) * X(r) + 1) * sum([coeff * GrothendieckB.get(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = GrothendieckB.get(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print()
            print('base:', -X(0)**-1 * GrothendieckB.get(v))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', -X(0)**-1 * ans[u] * GrothendieckB.get(u) * (X(r) - 1))
                print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def test_c_grothendieck_transitions(n=3):
    orders = {}
    for w in SignedPermutation.all(n):
        r = [i for i in range(1, n) if w(i) > w(i + 1)]
        if len(r) == 0:
            continue
        r = max(r)
        s = max([i for i in range(r + 1, n + 1) if w(i) < w(r)])
        v = (w * SignedPermutation.reflection_t(r, s, n)).inflate(n + 1)

        chain = []
        # chain += [SignedPermutation.reflection_t(r, j, n + 1) for j in range(r + 1, n + 2)]
        chain += [SignedPermutation.reflection_s(i, r, n + 1) for i in range(n + 1, 0, -1) if i != r]
        chain += [SignedPermutation.reflection_s(r, r, n + 1)]
        chain += [SignedPermutation.reflection_t(i, r, n + 1) for i in range(1, r)]

        worked = []
        for ch in [tuple(chain)]:
            ans = expand_reflection_chain(v, ch, lambda x: x.length())

            expected = GrothendieckC.get(v).set_variable(0, -1) + (X(r) - 1) * sum([coeff * GrothendieckC.get(z).set_variable(0, -1) for (z, coeff) in ans.dictionary.items()])
            actual = GrothendieckC.get(w).set_variable(0, -1)
            print('n =', n, 'w =', w, 'r =', r)
            print()
            print('base:', GrothendieckC.get(v).set_variable(0, -1))
            print()
            for v in sorted(ans, key=len):
                print(len(v), ':', v, ':', ans[v] * GrothendieckC.get(v).set_variable(0, -1) * (X(r) - 1))
                print()
            # print('***')
            # for v in bns:
            #     print(v, ':', GrothendieckC.get(v).set_variable(0, -1))
            #     print()
            print()
            #print(' got:', expected)
            #print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
            #input('')
            if expected == actual:
                worked.append(ch)

        if r not in orders:
            orders[r] = set(worked)
        else:
            orders[r] &= set(worked)

        print('n =', n, 'w =', w, 'r =', r, len(worked), len(orders[r]))
        input('')

        # for r in range(1, n + 1):
        #     chain = []
        #     chain += [SignedPermutation.reflection_t(r, j, n + 1) for j in range(r + 1, n + 2)]
        #     chain += [SignedPermutation.reflection_s(r, r, n + 1)]
        #     chain += [SignedPermutation.reflection_s(i, r, n + 1) for i in range(n + 1, 0, -1) if i != r]
        #     chain += [SignedPermutation.reflection_t(i, r, n + 1) for i in range(1, r)]

        #     worked = []
        #     for cha in signed(chain):
        #         for chb in signed(chain):
        #         # for v in [1]: #range(2**len(chain)):
        #         #     original_v = v
        #         #     signs = []
        #         #     for i in range(len(chain)):
        #         #         signs.append(1 if v % 2 == 0 else -1)
        #         #         v = v // 2
        #         #     sch = list(zip(ch, signs))
        #             #ans = expand_reflection_chain(w.inflate(n + 1), cha, lambda x: x.length())
        #             # for z in ans:
        #             #     print(z)
        #             # print()
                    
        #             ans = expand_reflection_chain(w.inflate(n + 1), cha, lambda x: x.length())
        #             bns = expand_reflection_chain(w.inflate(n + 1), chb, lambda x: x.length())
                    
        #             expected = sum([coeff * GrothendieckC.get(z).set_variable(0, -1) for (z, coeff) in ans.dictionary.items()])
        #             actual = sum([coeff * GrothendieckC.get(z).set_variable(0, -1) for (z, coeff) in bns.dictionary.items()]) * (1 - X(r))
        #             #actual = GrothendieckC.get(w.inflate(n + 1)).set_variable(0, -1) * (1 - X(r))

        #             if expected == actual:
        #                 worked.append((cha, chb))
        #             #else:
        #                 # print('v =', original_v)
        #             # print()
        #             # print(chain)
        #             # print((cha, chb))
        #             # print(ans)
        #             # print()
        #             # for v in ans:
        #             #     print(v, ':', GrothendieckC.get(v).set_variable(0, -1))
        #             #     print()
        #             # for v in bns:
        #             #     print(v, ':', GrothendieckC.get(v).set_variable(0, -1) * (1 - X(r)))
        #             #     print()
        #             # print()
        #             # print(expected)
        #             # print()
        #             # print(actual)
        #             # print()
        #             # print('diff:', expected - actual)
        #             # print()
        #             # if expected == actual and set(cha) | set(chb) == set(chain) and len(cha) + len(chb) == len(chain):
        #             #     input('')

        #     print('n =', n, 'w =', w, 'r =', r, len(worked))
        #     input('')

        #     # if expected != actual:
        #     #     print(chain)
        #     #     print(ans)
        #     #     print()
        #     #     for v in ans:
        #     #         print(v, ':', GrothendieckC.get(v).set_variable(0, -1))
        #     #         print()
        #     #     print()
        #     #     print(expected)
        #     #     print()
        #     #     print(actual)
        #     #     print()
        #     # input('\n' + str(expected == actual))
        #     #assert expected == actual


def test_d_grothendieck_transitions(n=4):
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

        expected = GrothendieckD.get(v) - (X(0) * X(r) + 1) * sum([coeff * GrothendieckD.get(z) for (z, coeff) in ans.dictionary.items()])
        expected *= -X(0)**-1

        actual = GrothendieckD.get(w)

        print('w =', w, 'r =', r, 's =', s, 'v =', v)
        if expected != actual:
            print()
            print('base:', -X(0)**-1 * GrothendieckD.get(v))
            print()
            for u in sorted(ans, key=len):
                print(len(u), ':', u, ':', -X(0)**-1 * ans[u] * GrothendieckD.get(u) * -(X(0) * X(r) - 1))
                print()
            print('want:', actual)
            print()
            print('diff:', expected - actual)
            print()
        assert expected == actual


def test_a_symmetric_transitions(n=4, numvars=2):
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


def test_b_symmetric_transitions(n=4, numvars=2):
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


def test_c_symmetric_transitions(n=4, numvars=2):
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


def test_b(rank=4, numvars=2):
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


def test_c(rank=4, numvars=2):
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
