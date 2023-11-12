from stable.kromatic import (
    kromatic,
    kmonomial,
    oriented_kromatic,
    oriented_kromatic_polynomial,
    reoriented_kromatic,
    reoriented_kromatic_polynomial,
    oriented_expander,
    strict_galois_number,
    posets,
    small_multipermutations,
    incomparability_graph,
    is_cluster_graph,
    is_connected_graph,
    kpowersum_expansion,
    Y_INDEX,
)
from stable.polynomials import Y, X
from stable.utils import mn_G_expansion, schur_expansion
from stable.symmetric import SymmetricPolynomial
from stable.partitions import Partition
from stable.tableaux import Tableau
from stable.quasisymmetric import Quasisymmetric
import itertools
import math
from sympy.ntheory import factorint


def test_alt_p_e_expansion(n=3, k=2):
    v = list(range(1, k + 1))
    x = kromatic(n, v, []) - kromatic(n, v, complete_graph(v))
    a = GE_expansion(x)
    ell = max([0] + [sum(mu) for mu in a])
    for i in range(ell, -1, -1):
        vec = Vector({mu: coeff for mu, coeff in a.items() if sum(mu) == i})
        if vec:
            print(vec, '\n\n+\n')


def test_p_expansion(n=3):
    x = kromatic(n, [1,2,3], [(1,2),(2,3),])
    a = p_expansion(x)

    fac = math.lcm(*[t.denominator for t in a.dictionary.values()])

    b = a * fac
    y = x.polynomial().set(0, 1) * fac

    assert all(coeff.denominator == 1 for (mu, coeff) in b.items())
    expected = sum(p(n, mu) * coeff.numerator for (mu, coeff) in b.items()).polynomial()
    if y != expected:
        print(y)
        print()
        print(expected)
        print()
    assert y == expected

    ell = max([sum(mu) for mu in a])
    for i in range(ell, -1, -1):
        subset = {coeff.denominator for mu, coeff in a.items() if sum(mu) == i}
        if subset:
            print(Vector({mu: coeff for mu, coeff in a.items() if sum(mu) == i}), '\n\n+\n')
            ell = math.lcm(*subset)
            print(i, ':', ell, '=', factorint(ell), '\n\n+\n')


def test_grothendieck_p_tableaux(nn=5, kk=5):
    for n in range(1 + nn):
        for p in posets(n):
            v, e = incomparability_graph(n, p)
            for k in range(1, 1 + kk):
                x = kromatic(k, v, e)
                exp = mn_G_expansion(x).set_variable(0, 1)
                print()
                print()
                print('new poset:', n, p, exp)
                print()
                r = max([sum(nu) for nu in exp], default=0)
                for nu in set(Partition.all(r)) | set(exp):
                    if len(nu) > k:
                        continue
                    tabs = list(Tableau.grothendieck_p_tableaux(n, p, nu))
                    print(nu, ':', len(tabs), exp[nu])
                    if len(tabs) != exp[nu]:
                        print(tabs)
                    assert len(tabs) == exp[nu]


def test_oriented_grothendieck_p_tableaux(nn=5, kk=5):
    for n in range(1 + nn):
        for p in posets(n):
            v, e = incomparability_graph(n, p)
            print()
            print()
            print('poset:', n, p, 'graph:', v, e)
            print()
            for k in range(1, 1 + kk):
                try:
                    x = oriented_kromatic(k, v, e)
                except:
                    print('nvars:', k, 'expansion:', '(not symmetric)', 'cluster graph:', is_cluster_graph(e))
                    print()
                    continue
                exp = oriented_expander(mn_G_expansion, x).set_variable(0, 1)
                print()
                print('nvars:', k, 'expansion:', exp)
                print()
                r = max([sum(nu) for nu in exp], default=0)
                for nu in set(Partition.all(r)) | set(exp):
                    if len(nu) > k:
                        continue
                    tabs = list(Tableau.grothendieck_p_tableaux(n, p, nu))
                    # print(nu, ':', len(tabs), '==', exp[nu])
                    counter = 0
                    for t in tabs:
                        inv = t.poset_inversions(p)
                        counter += Y()**inv
                    # print()
                    print('  ', nu, ':', exp[nu], '==', counter, exp[nu] == counter)
                    # print()
                    assert exp[nu] == counter
                    # assert len(tabs) == exp[nu]


def test_strict_galois_number(rr=8, nn=8):
    f = strict_galois_number
    for r in range(rr + 1):
        for n in range(nn + 1):
            expected = 0
            for i in range(n):
                qq = 1
                for j in range(r - i + 1, r + 1):
                    qq *= (1 - Y()**j)
                for j in range(n - 1 - i, n + 1):
                    expected += math.comb(n, j) * math.comb(j, n - 1 - i) * (-1)**i * qq * f(r - i, j)
            print(r, n)
            assert expected == f(r + 1, n)


def factorial(n):
    ans = 1
    for i in range(n):
        ans *= i + 1
    return ans


stirling_cache = {}


def stirling2(n, k):
    assert 0 <= k <= n
    if (n, k) not in stirling_cache:
        if n == k:
            ans = 1
        elif k == 0:
            ans = 0
        else:
            ans = stirling2(n - 1, k) * k + stirling2(n - 1, k - 1)
        stirling_cache[n, k] = ans
    return stirling_cache[n, k]


def test_kmonomial(nvars=3, k=10):
    def multiplicity(mu, m):
        return len([a for a in mu if a == m])

    def coefficient(nu, mu):
        ans = 1
        for a in set(mu):
            r = multiplicity(mu, a)
            s = multiplicity(nu, a)
            ans *= stirling2(s, r) * factorial(r)
        return ans

    def expand(mu, n):
        if len(mu) == 0:
            yield mu
        else:
            m = mu[0]
            r = multiplicity(mu, m)
            for a in range(n - len(mu) + 1):
                for nu in expand(mu[r:], n - r - a):
                    yield (m,) * (a + r) + nu

    for mu in Partition.all(k, max_part=3):
        if len(mu) == 0:
            continue
        mu = Partition.transpose(mu)
        print(mu)
        expected = 0
        for nu in expand(mu, nvars):
            expected += SymmetricPolynomial.monomial(nvars, nu) * coefficient(nu, mu)
        actual = kmonomial(nvars, mu).set_variable(0, 1)
        assert expected == actual


def natural_unit_interval_order_incompability_graph(n, r):
    v = list(range(1, n + 1))
    e = [(i, j) for i in v for j in v if 0 < j - i < r]
    return v, e


def graphs(n):
    v = list(range(1, n + 1))
    e = list(itertools.combinations(v, 2))
    for k in range(0, len(e) + 1):
        for edges in itertools.combinations(e, k):
            yield v, edges


def is_claw_free(vertices, edges):
    e = {v: set() for v in vertices}
    for a, b in edges:
        e[a].add(b)
        e[b].add(a)
    for v in vertices:
        for a in e[v]:
            for b in e[v]:
                if a == b:
                    continue
                for c in e[v]:
                    if a == c or b == c:
                        continue
                    if a not in e[b] and a not in e[c] and b not in e[c]:
                        return False
    return True


def test_G_expansion(nvars=3, nverts=3):
    for v, e in graphs(nverts):
        f = kromatic(nvars, v, e)
        exp = mn_G_expansion(f)
        icf = is_claw_free(v, e)
        pos = exp.is_nonnegative()
        print(v, e, icf, pos, exp)
        print()
        assert (not icf) or pos


def test_kpowersum_expansion(nvars=3, nverts=3, cutoff=7):
    for v, e in graphs(nverts):
        f = kromatic(nvars, v, e)
        exp = kpowersum_expansion(f, cutoff)
        icf = is_claw_free(v, e)
        pos = exp.is_nonnegative()
        print(v, e, icf, pos, exp)
        print()
        # assert (not icf) or pos


def test_G_oriented_expansion(nvars=3, nverts=3, re=False):
    # v1, e1 = [1, 2, 3, 4, 5], ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5))
    # v2, e2 = [1, 2, 3, 4, 5], ((1, 2), (1, 5), (2, 3), (3, 4), (4, 5))
    # for v, e in [(v1, e1), (v2, e2)]:
    for v, e in graphs(nverts):
        icg = is_cluster_graph(e)
        isc = is_connected_graph(v, e)
        icf = is_claw_free(v, e)
        if icg or not isc:
            continue
        try:
            f = reoriented_kromatic(nvars, v, e) if re else oriented_kromatic(nvars, v, e)
        except:
            ind = {i + 1: 1 for i in range(nvars)}

            f = reoriented_kromatic_polynomial(nvars, v, e) if re else oriented_kromatic_polynomial(nvars, v, e)
            f = f.set(0, 0).set(Y_INDEX, 0)

            coeffs = []
            for i in range(nvars):
                ind[i + 1] += 1
                coeffs.append(f[ind])
                ind[i + 1] -= 1

            print(v, e, '(not symmetric)', 'cluster graph:', icg)
            print()
            print('  ', coeffs)
            print()
            continue
        exp = oriented_expander(mn_G_expansion, f)
        pos = exp.is_nonnegative()
        print(v, e, 'claw-free:', icf, 'positive:', pos, 'cluster graph:', icg)
        print()
        assert (not icf) or pos


def test_G_reoriented_expansion(nvars=3, nverts=3):
    # expected: only symmetric for claster graphs, then alays G-positive
    test_G_oriented_expansion(nvars, nverts, True)


def test_multifundamental_oriented_expansion(nvars=3, nverts=3, re=False):
    # fails for v, e = [1, 2, 3, 4], ((1, 2), (1, 3))
    # fails for v, e = [1, 2, 3, 4], ((1, 2), (1, 3), (1, 4))
    for v, e in graphs(nverts):
        f = reoriented_kromatic_polynomial(nvars, v, e) if re else oriented_kromatic_polynomial(nvars, v, e)
        exp = Quasisymmetric.multifundamental_expansion(f)
        assert Quasisymmetric.from_expansion(nvars, exp, Quasisymmetric.multifundamental) == f
        icf = is_claw_free(v, e)
        pos = exp.is_nonnegative()
        print(v, e, icf, pos, exp)
        print()
        assert pos


def test_multifundamental_reoriented_expansion(nvars=3, nverts=3):
    test_multifundamental_oriented_expansion(nvars, nverts, True)


def composition_from_des(ell, des):
    if ell == 0:
        return ()
    # des = [0] + sorted(set(range(1, ell)) - set(des)) + [ell]
    des = [0] + sorted(des) + [ell]
    ans = []
    for i in range(1, len(des)):
        ans += [des[i] - des[i - 1]]
    return tuple(ans)


def leftreduce(w):
    seen = set()
    rw = []
    for a in w:
        if a not in seen:
            rw.append(a)
            seen.add(a)
    return tuple(rw)


def test_multifundamental_reoriented_posets(nn=5, kk=5):
    for n in range(1 + nn):
        for p in posets(n):
            v, e = incomparability_graph(n, p)
            print()
            print()
            print('poset:', n, p, 'graph:', v, e)
            print()
            for k in range(1, 1 + kk):
                x = reoriented_kromatic_polynomial(k, v, e)
                exp = oriented_expander(Quasisymmetric.multifundamental_expansion, x).set_variable(0, 1)
                print()
                print('nvars:', k, 'expansion:', exp)
                ted = exp._instantiate({})
                mm = max([sum(_) for _ in exp], default=0)
                for m in range(n, 1 + mm):
                    for w in small_multipermutations(n, m):
                        rw = leftreduce(w)
                        ginv = len([(i, j) for i in range(len(rw)) for j in range(i + 1, len(rw)) if rw[i] > rw[j] and ((rw[i], rw[j]) in e or (rw[j], rw[i]) in e)])
                        pdes = {m - (i + 1) for i in range(len(w) - 1) if (w[i + 1], w[i]) not in p}
                        alpha = composition_from_des(m, pdes)
                        if len(alpha) <= k:
                            # print(w, ginv, pdes, alpha)
                            ted += ted._instantiate({alpha: Y()**ginv})
                if exp != ted:
                    print('          expected:', ted)
                    print()
                assert exp == ted


def test_G_reoriented_expansion_nui(nvars=3, n=5, r=5):
    # expected: only symmetric for claster graphs
    test_G_oriented_expansion_nui(nvars, n, r, True)


def test_G_oriented_expansion_nui(nvars=3, n=5, r=5, re=False):
    seen = set()
    for nn in range(1, n + 1):
        for rr in range(1, r + 1):
            v, e = natural_unit_interval_order_incompability_graph(nn, rr)
            
            v, e = tuple(sorted(v)), tuple(sorted(e))
            if (v, e) not in seen:
                seen.add((v, e))
            else:
                continue
            
            icf = is_claw_free(v, e)
            try:
                f = reoriented_kromatic(nvars, v, e) if re else oriented_kromatic(nvars, v, e)
                exp = oriented_expander(mn_G_expansion, f)
                pos = exp.is_nonnegative()
                print(v, e, 'claw-free:', icf, 'positive:', pos, 'cluster graph:', is_cluster_graph(e))
                print()
                assert (not icf) or pos
            except:
                print(v, e, 'claw-free:', icf, '(not symmetric)', 'cluster graph:', is_cluster_graph(e))
                # assert False
