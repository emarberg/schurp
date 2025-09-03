from crystals import(
    AbstractGLCrystal,
    SuperGLCrystal,
)
from stable.tableaux import Tableau
from stable.symmetric import SymmetricPolynomial
from stable.polynomials import X, Y, Polynomial
from stable.utils import hs, HG, hs_expansion, HG_expansion, schur, G
from stable.vectors import Vector
from stable.partitions import Partition
from stable.permutations import Permutation
from tests.test_crystals import draw_graph


def test_tableau_raising(mm=2, nn=2, k=3):
    # assert mm == 1 or nn == 1

    def operator(b, w):
        m, n = b.rank
        ans = w
        for i in range(m - 1, -n, -1):
            ans = ans if ans is None else b.capital_e_operator(-i, ans)
        return ans

    for m in [mm]:
        for n in [nn]:
            for l in range(1, k + 1):
                b = SuperGLCrystal.sqrtcrystal_of_words(l, m, n)
                for comp in b.get_components():
                    ch = SymmetricPolynomial.from_super_polynomial(comp.character())
                    expand = HG_expansion(ch)
                    print('m =', m, 'n =', n, 'k =', l, ': ch( len', len(comp), ') =', expand)
                    assert min(expand.values()) > 0

                    vertices = []
                    edges = []
                    weights = {}

                    for w in comp:
                        p, q = b.bihecke(w)
                        if p not in vertices:
                            vertices.append(p)
                            weights[p] = b.weight(w)
                            fw = operator(b, w)
                            if fw is not None:
                                fp, fq = b.bihecke(fw)
                                edge = (0, p, fp)
                                if p != fp and edge not in edges:
                                    edges.append(edge)
                    
                    colors = lambda p: 'white' if (len(p[0]) + len(p[1])) % 2 else 'cyan'
                    draw_graph(vertices, edges, colors=colors)
                    input('\n')


def test_tableau_reduction(mm=2, nn=2, k=3):
    for m in [mm]:
        for n in [nn]:
            for l in range(1, k + 1):
                b = SuperGLCrystal.sqrtcrystal_of_words(l, m, n)
                for comp in b.get_components():
                    ch = SymmetricPolynomial.from_super_polynomial(comp.character())
                    expand = HG_expansion(ch)
                    print('m =', m, 'n =', n, 'k =', l, ': ch( len', len(comp), ') =', expand)
                    assert min(expand.values()) > 0

                    vertices = []
                    edges = []
                    weights = {}

                    diff = lambda p, q: abs(len(p[0]) + len(p[1]) - len(q[0]) - len(q[1])) == 1

                    # tab = {}
                    # for w in comp:
                    #     p, q = b.bihecke(w)
                    #     tab[p] = tab.get(p, set()) | {w}

                    # for p in tab:
                    #     vert_a = (tab[p], p)
                    #     vertices.append(vert_a)
                    #     for w in tab[p]:
                    #         for i in [0]: #b.indices:
                    #             fw = b.f_operator(i, w)
                    #             if fw is not None:
                    #                 fp, fq = b.bihecke(fw)
                    #                 vert_b = (tab[fp], fp)
                    #                 edge = (i, vert_a, vert_b)
                    #                 if p != fp and edge not in edges: # and diff(p, fp):
                    #                     edges.append(edge)
  
                    # colors = lambda p: 'white' if (len(p[1][0]) + len(p[1][1])) % 2 else 'cyan'
                    # draw_graph(vertices, edges, colors=colors)
                    # input('\n')

                    for w in comp:
                        p, q = b.bihecke(w)
                        if p not in vertices:
                            vertices.append(p)
                            weights[p] = b.weight(w)
                        for i in [0]: #b.indices:
                             fw = b.f_operator(i, w)
                             if fw is not None:
                                 fp, fq = b.bihecke(fw)
                                 edge = (i, p, fp)
                                 if p != fp and edge not in edges: # and diff(p, fp):
                                     edges.append(edge)
                    
                    colors = lambda p: 'white' if (len(p[0]) + len(p[1])) % 2 else 'cyan'
                    draw_graph(vertices, edges, colors=colors)
                    input('\n')


def test_sqrt_connected_lemma(k, *mu):
    def e(v, *args):
        ans = v
        for i in args:
            if ans is not None:
                ans = ans.sqrt_super_e_operator(i)
        return ans

    def f(v, *args):
        ans = v
        for i in args:
            if ans is not None:
                ans = ans.sqrt_super_f_operator(i)
        return ans

    def E(v, *args):
        ans = v
        for i in args:
            while True:
                bns = ans.sqrt_super_e_operator(i)
                if bns is None:
                    break
                ans = bns
        return ans

    def F(v, *args):
        ans = v
        for i in args:
            while True:
                bns = ans.sqrt_super_f_operator(i)
                if bns is None:
                    break
                ans = bns
        return ans

    assert mu[-1] > 0
    n = len(mu)
    v = {(1, i): (-1) for i in range(1, k + 1)}
    v[1, k + 1] = (-1, n)
    for col in range(len(mu)):
        for row in range(mu[col]):
            v[row + 2, col + 1] = col + 1
    v = Tableau(v)
    return v, e, f, E, F

def test_super_sergeev_pragacz(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            D = Polynomial.one()
            for a in range(1, i):
                for b in range(a + 1, i + 1):
                    D *= X(a) - X(b)
            for a in range(1, j):
                for b in range(a + 1, j + 1):
                    D *= Y(-a) - Y(-b)

            for kappa in Partition.subpartitions(i * (j,)):
                kappa_t = Partition.transpose(kappa)

                pi = Polynomial.one()
                for (a,b) in Partition.shape(kappa):
                    pi *= X(a) + Y(-b) + X(a) * Y(-b)

                for mu_t in Partition.all(k, max_part=len([a for a in kappa if a == j])):
                    mu = Partition.transpose(mu_t)
                    for nu_t in Partition.all(k, max_part=len([a for a in kappa_t if a == i])):
                        nu = Partition.transpose(nu_t)
                        lam = tuple(Partition.get(kappa, a) + Partition.get(mu, a) for a in range(1, i + 1)) + nu_t
                        print('m =', i, 'n =', j, 'lambda =', lam, 'kappa =', kappa, 'mu =', mu, 'nu =', nu)

                        expected = HG(i, j, lam).superpolynomial() * D

                        one = Polynomial.from_tuple([0] + [i - t + Partition.get(mu, t) for t in range(1, i + 1)])
                        for a in range(1, i + 1):
                            one *= (1 + X(a))**(a-1)

                        two = Polynomial.from_tuple([0] + [j - t + Partition.get(nu, t) for t in range(1, j + 1)]).swap_xy()
                        for b in range(1, j + 1):
                            one *= (1 + Y(-b))**(b-1)

                        product = pi * one * two

                        summation = Polynomial.zero()
                        for v in Permutation.all(i):
                            for w in Permutation.all(j):
                                summation += (-1)**(len(v) + len(w)) * product.permute_x(v).permute_y(w)

                        Partition.print_diagram(lam)
                        if expected != summation:
                            print()
                            print(expected)
                            print()
                            print(summation)
                        print()

                        assert expected == summation


def test_sergeev_pragacz(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            D = Polynomial.one()
            for a in range(1, i):
                for b in range(a + 1, i + 1):
                    D *= X(a) - X(b)
            for a in range(1, j):
                for b in range(a + 1, j + 1):
                    D *= Y(-a) - Y(-b)

            for kappa in Partition.subpartitions(i * (j,)):
                kappa_t = Partition.transpose(kappa)

                pi = Polynomial.one()
                for (a,b) in Partition.shape(kappa):
                    pi *= X(a) + Y(-b)

                for mu_t in Partition.all(k, max_part=len([a for a in kappa if a == j])):
                    mu = Partition.transpose(mu_t)
                    for nu_t in Partition.all(k, max_part=len([a for a in kappa_t if a == i])):
                        nu = Partition.transpose(nu_t)
                        lam = tuple(Partition.get(kappa, a) + Partition.get(mu, a) for a in range(1, i + 1)) + nu_t
                        print('m =', i, 'n =', j, 'lambda =', lam, 'kappa =', kappa, 'mu =', mu, 'nu =', nu)

                        expected = hs(i, j, lam).superpolynomial() * D

                        one = Polynomial.from_tuple([0] + [i - t + Partition.get(mu, t) for t in range(1, i + 1)])
                        two = Polynomial.from_tuple([0] + [j - t + Partition.get(nu, t) for t in range(1, j + 1)]).swap_xy()
                        product = pi * one * two

                        summation = Polynomial.zero()
                        for v in Permutation.all(i):
                            for w in Permutation.all(j):
                                summation += (-1)**(len(v) + len(w)) * product.permute_x(v).permute_y(w)

                        Partition.print_diagram(lam)
                        if expected != summation:
                            print()
                            print(expected)
                            print()
                            print(summation)
                        print()

                        assert expected == summation


def test_super_factorization(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            pi = Polynomial.one()
            for a in range(1, i + 1):
                for b in range(1, j + 1):
                    pi *= X(a) + Y(-b) + X(a) * Y(-b)

            for mu_t in Partition.all(k, max_part=i):
                mu = Partition.transpose(mu_t)
                for nu_t in Partition.all(k, max_part=j):
                    nu = Partition.transpose(nu_t)
                    lam = tuple(j + Partition.get(mu, a) for a in range(1, i + 1)) + nu_t
                    print('m =', i, 'n =', j, 'lambda =', lam, 'mu =', mu, 'nu =', nu)
                    expected = HG(i, j, lam).superpolynomial()
                    right = G(i, mu).polynomial().set_variable(0, 1)
                    bottom = G(j, nu).polynomial().swap_xy().set_variable(0, 1)
                    product = pi * right * bottom

                    # Partition.print_diagram(lam)
                    # print()
                    # print(expected)
                    # print(pi)
                    # print(right)
                    # print(bottom)
                    # print()

                    assert expected == product
                    cancel = expected.set_variable(-1, X(0)**-1 - 1).set_variable(1, X(0) - 1)
                    assert cancel.set(0, 0) == cancel


def test_factorization(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            pi = Polynomial.one()
            for ii in range(1, i + 1):
                for jj in range(1, j + 1):
                    pi *= X(ii) + Y(-jj)

            for mu_t in Partition.all(k, max_part=i):
                mu = Partition.transpose(mu_t)
                for nu_t in Partition.all(k, max_part=j):
                    nu = Partition.transpose(nu_t)
                    lam = tuple(j + Partition.get(mu, a) for a in range(1, i + 1)) + nu_t
                    print('m =', i, 'n =', j, 'lambda =', lam, 'mu =', mu, 'nu =', nu)
                    expected = hs(i, j, lam).superpolynomial()
                    right = schur(i, mu).polynomial()
                    bottom = schur(j, nu).polynomial().swap_xy()
                    product = pi * right * bottom

                    # Partition.print_diagram(lam)
                    # print()
                    # print(expected)
                    # print(pi)
                    # print(right)
                    # print(bottom)
                    # print()

                    assert expected == product
                    cancel = expected.set_variable(-1, -X(0)).set_variable(1, X(0))
                    assert cancel.set(0, 0) == cancel


def test_sqrt_super_crystals(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            b = SuperGLCrystal.standard_sqrtcrystal(i, j)
            for l in range(1, k + 1):
                for comp in b.get_components():
                    ch = SymmetricPolynomial.from_super_polynomial(comp.character())
                    expand = HG_expansion(ch)
                    print('m =', i, 'n =', j, 'k =', l, ': ch( component of len', len(comp), ') =', expand)
                    assert min(expand.values()) > 0
                if l < k:
                    b = b.tensor(b)
            print()


def test_super_crystals(m=2, n=2, k=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            b = SuperGLCrystal.standard_object(i, j)
            for l in range(1, k + 1):
                for comp in b.get_components():
                    ch = SymmetricPolynomial.from_super_polynomial(comp.character())
                    expand = hs_expansion(ch)
                    print('m =', i, 'n =', j, 'k =', l, ': ch( component of len', len(comp), ') =', expand)
                    assert min(expand.values()) > 0
                    assert len(expand) == 1
                    assert list(expand.values()) == [1]
                if l < k:
                    b = b.tensor(b)
            print()



def test_hook_schur_expansion(m=4, n=4):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            print('m =', i, 'n =', j)
            one = hs(i, j, (1,))
            two = hs(i, j, (1, 1)) + hs(i, j, (2,))
            hs_expansion(one) == Vector({(1,): 1})
            hs_expansion(two) == Vector({(1,1): 1, (2,): 1})


def test_hook_grothendieck_expansion(m=4, n=4):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            print('m =', i, 'n =', j)
            one = HG(i, j, (1,))
            two = HG(i, j, (1, 1)) + HG(i, j, (2,))
            HG_expansion(one) == Vector({(1,): 1})
            HG_expansion(two) == Vector({(1,1): 1, (2,): 1, (2, 1): 1})


def test_hook_schur(m=3, n=3):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            b = SuperGLCrystal.standard_object(i, j)
            c = b.tensor(b)
            
            print('m =', i, 'n =', j)
            
            f = b.character()
            symf = SymmetricPolynomial.from_super_polynomial(f)
            assert symf == hs(i, j, (1,))

            g = c.character()
            symg = SymmetricPolynomial.from_super_polynomial(g)
            assert symg == hs(i, j, (1, 1)) + hs(i, j, (2,))


def test_hook_grothendieck(m=2, n=2):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            b = SuperGLCrystal.standard_sqrtcrystal(i, j)
            c = b.tensor(b)
            
            print('m =', i, 'n =', j)
            
            f = b.character()
            symf = SymmetricPolynomial.from_super_polynomial(f)
            assert symf == HG(i, j, (1,))

            g = c.character()
            symg = SymmetricPolynomial.from_super_polynomial(g)
            assert symg == HG(i, j, (1, 1)) + HG(i, j, (2,)) + HG(i, j, (2, 1))

            tab11 = SuperGLCrystal.sqrtcrystal_from_partition((1, 1,), i, j)
            tab2 = SuperGLCrystal.sqrtcrystal_from_partition((2,), i, j)
            tab21 = SuperGLCrystal.sqrtcrystal_from_partition((2, 1,), i, j)
            for comp in c.get_components():
                print('* component of size', len(comp), 'compared to')
                print('  1,1', comp.character() == tab11.character(), ':', c.find_isomorphism(comp, tab11) is not None)
                print('  2,0', comp.character() == tab2.character(), ':', c.find_isomorphism(comp, tab2) is not None)
                print('  2,1', comp.character() == tab21.character(), ':', c.find_isomorphism(comp, tab21) is not None)
                print()

            
