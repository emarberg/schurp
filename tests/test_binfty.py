from keys import decompose_into_lascoux
from crystals import (
    AbstractGLCrystal, 
    AbstractQCrystal,
    AbstractPrimedQCrystal
)
from binfty import InfiniteCrystal
from permutations import Permutation
from stable.partitions import Partition
from stable.tableaux import Tableau
from stable.utils import G, G_expansion_no_beta, SymmetricPolynomial, schur_expansion, g_expansion, schur
from stable.vectors import Vector


dual_blambda_lookup = {}
dual_blambda_cache = {}
dual_blambda_character_cache = {}


def dual_blambda_character(n, mu):
    if (n, mu) not in dual_blambda_character_cache:
        b = dual_sqrt_b_lambda(n, mu)
        dual_blambda_character_cache[n, mu] = SymmetricPolynomial.from_polynomial(b.character())
    return dual_blambda_character_cache[n, mu]


def dual_expand(n, c):
    a = schur_expansion(c)
    ans = Vector()
    while a != 0:
        mu = max(a)
        c = a[mu]
        if (n, mu) not in dual_blambda_lookup:
            #print()
            #print('. . . computing X_%s%s)' % (n, str(mu)))
            #print()
            dual_blambda_lookup[(n, mu)] = schur_expansion(dual_blambda_character(n, mu))
        a = a - c * dual_blambda_lookup[(n, mu)]
        ans += Vector({mu: c})
    return ans


def test_dual_positivity(n, k):
    g = dual_sqrt_b_lambda(n, (1,))
    b = g
    for i in range(2, k + 1):
        print(i, 'tensor factors')
        print()
        b = g.tensor(b)
        for c in b.get_components():
            ch = SymmetricPolynomial.from_polynomial(c.character())
            actual = dual_expand(n, ch)

            expected = Vector()
            for e, wt in c.get_highest_weights():
                expected += Vector({Partition.trim(wt): 1})

            print('  ch =', actual)
            print('     =', expected)
            print()
            assert max(expected.values()) > 0
            assert actual == expected


def test_dual_pieri(n, k):
    from stable.utils import g

    for mu in Partition.all(k, max_row=n) if type(k) == int else k:
        if len(mu) <= 1 or max(mu) <= 1:
            continue
        for j in Partition.all(5, max_row=n):
            if len(j) <= 1 or max(j) <= 1:
                continue
            # base = g_expansion(g(n, mu) * g(n, j)).set_variable(0, 1)
            # print('    base:', base)
            for k in range(len(mu) + len(j), n + 1):
                c = dual_blambda_character(k, mu)
                o = dual_blambda_character(k, j)
                actual = dual_expand(k, c * o)
                expected = g_expansion(g(k, mu) * g(k, j)).set_variable(0, 1)
                print()
                print(mu, '*', j)
                print('    n =', k, ':', actual - expected)
                assert actual == expected


def test_lp(n, thresh=10):
    lp = InfiniteCrystal.sqrt_LP(n)
    elems = lp.vertices(thresh, thresh)

    def floor(x):
        return 2 * (x // 2)

    def ceil(x):
        return 2 * ((x + 1) // 2)

    def get(i, x):
        assert 1 <= i <= n
        return 0 if i == n else x[i - 1]

    def sub(i, x):
        x = list(x)
        x[i - 1] -= 1
        return tuple(x)

    def add(i, x):
        x = list(x)
        x[i - 1] += 1
        return tuple(x)

    for v in elems:
        weight = n * [0]
        for i in range(1, n):
            weight[i - 1] += ceil(get(i, v)) // 2
            weight[i] -= floor(get(i, v)) // 2
        weight = tuple(weight)
        assert weight == lp.weight(v)

        for i in range(1, n):
            assert get(i, v) <= 0
            assert get(i, v) <= ceil(get(i + 1, v))

            # test e
            ev = lp.e_operator(i, v)
            case = ''
            if get(i, v) == ceil(get(i + 1, v)):
                case += 'a'
                expected = None
            if (i > 1 and get(i - 1, v) > get(i, v)) and get(i, v) < ceil(get(i + 1, v)):
                case += 'b'
                expected = sub(i - 1, v)
            if (i == 1 or get(i - 1, v) <= get(i, v)) and get(i, v) < ceil(get(i + 1, v)):
                case += 'c'
                expected = add(i, v)
            if expected != ev:
                print('case', case, '| e_%s' % i, ':', v, '=', ev, '=?=', expected)
                #assert get(i, expected) > ceil(get(i + 1, expected))
                #assert ev is None
            assert len(case) == 1
            assert expected == ev
            
            # test f
            fv = lp.f_operator(i, v)
            case = ''
            if (i > 1 and get(i - 1, v) == ceil(get(i, v))):
                case += 'a'
                expected = None
            if (i > 1 and get(i - 1, v) == get(i, v) and get(i, v) % 2 != 0):
                case += 'b'
                expected = add(i - 1, v)
            if case == '':
                case += 'c'
                expected = sub(i, v)
            if expected != fv:
                print('case', case, '| f_%s' % i, ':', v, '=', fv, '=?=', expected)
            assert len(case) == 1
            assert expected == fv

            # test eps
            eps = lp.e_string(i, v)
            case = ''
            if i > 1 and get(i - 1, v) > get(i, v):
                case += 'a'
                expected = -get(i, v) + ceil(get(i + 1, v)) + 1
            if case == '':
                case += 'b'
                expected = -get(i, v) + ceil(get(i + 1, v))
            if expected != eps:
                print('case', case, '| eps_%s' % i, ':', v, '=', eps, '=?=', expected)
                input('\n')
            assert len(case) == 1
            assert expected == eps
            
            # test phi
            phi = lp.f_string(i, v)
            case = ''
            if i == 1:
                case += 'a'
                expected = get(i, v)
            if i > 1 and get(i - 1, v) <= get(i, v):
                case += 'b'
                expected = -floor(get(i - 1, v)) + get(i, v)
            if i > 1 and get(i - 1, v) > get(i, v):
                case += 'c'
                expected = -get(i - 1, v) + get(i, v) + 1
            if expected != phi:
                print('case', case, '| phi_%s' % i, ':', v, '=', phi, '=?=', expected)
                input('\n')
            assert len(case) == 1
            assert expected == phi


def test_elementary_squared(n=3, p=3, q=3):
    b = InfiniteCrystal.binfty(n)
    for w in Permutation.longest_element(n).get_reduced_words():
        c = InfiniteCrystal.binfty_squared(n, w)
        boolean = InfiniteCrystal.is_isomorphic(b, b.vertices(p,q), c, c.vertices(p,q))
        if boolean:
            print(w)


def test_r_lambda(n=3, k=10):
    for w in Permutation.longest_element(n).get_reduced_words():
        b = InfiniteCrystal.binfty(n, w)
        for mu in Partition.all(k, max_row=n):
            print(n, ':', 'mu =', mu)
            r = InfiniteCrystal.r_lambda(mu, n)
            t = InfiniteCrystal.tensor(r, b).finitize(0)
            u = AbstractGLCrystal.from_partition(mu, n)
            assert AbstractGLCrystal.find_isomorphism(t, u) is not None


def sqrt_r_tensor(mu, n, b):
    r = InfiniteCrystal.sqrt_r_lambda(mu, n)
    h = 0
    ans = []
    while True:
        h += 1
        # print('  h = ', h)
        bns = InfiniteCrystal.sqrt_tensor(r, b.temper(b.vertices(0, h))).finitize()
        if len(bns) > len(ans):
            ans = bns
        else:
            return ans


def test_sqrt_r_lambda(n=3, k=10, test_lasc=False):
    for x in _test_sqrt_r_lambda(n, k, test_lasc, False):
        yield x


def test_dual_sqrt_r_lambda(n=3, k=10, test_lasc=False):
    for x in _test_sqrt_r_lambda(n, k, test_lasc, True):
        yield x


def univar(f, n=10):
    for i in range(1, n + 1):
        f = f.set_variable(i, f.x(0))
    return f


def investigate_characters(n, k):
    from stable.utils import g_expansion, G_expansion_no_beta, schur_expansion, g, j, G, schur, g_expansion

    def character(b):
        ch = b.character()
        countval = sum(ch.coeffs.values())
        # return subs(ch, n, -1), countval
        return SymmetricPolynomial.from_polynomial(ch), countval

    def subs(f, n, a=1):
        f = f.polynomial() if type(f) == SymmetricPolynomial else f
        for i in range(1, n + 1):
            f = f.substitute(i, f.x(i) + a)
        return SymmetricPolynomial.from_polynomial(f)

    def count(ch):
        return sum(ch.coeffs.values())

    for mu in Partition.all(k, max_row=n) if type(k) == int else k:
        if max(mu, default=0) != 1:
            continue

        #mu = tuple(1 + a for a in mu)
        #while len(mu) < n:
        #    mu += (1,)

        print()
        print()
        print('mu =', mu)
        print()

        #b = sqrt_b_lambda(n, mu)
        #c = dual_sqrt_b_lambda(n, mu)
        
        #print('     ch =', schur_expansion(character(b)))
        #print('      G =', schur_expansion(G(n, mu).set_variable(0, 1)))
        #print()
        
        elems = []
        for k in range(max(1, len(mu)), n + 1):
            c = dual_sqrt_b_lambda(k, mu)
            countval = len(c)
            elems.append(countval)
            print(' n =', k, ' dual ch =', schur_expansion(dual_blambda_character(k, mu)))
            #print(k, '        =', G_expansion_no_beta(character(c)), ':', len(c))
        #print('        =', schur_expansion(subs(character(c), -1)))
        print()
        print(' |dual B_n(mu)| =', elems)
        print()

        # elems = []
        # for k in range(max(1, len(mu)), n + 1):
        #     ch = G(k, mu).set_variable(0, 1)
        #     countval = sum(ch.polynomial().coeffs.values())
        #     elems.append(countval)
        #     print(' n =', k, '  svt ch =', schur_expansion(ch))
        #     #print(k, '        =', G_expansion_no_beta(character(c)), ':', len(c))
        # #print('        =', schur_expansion(subs(character(c), -1)))
        # print()
        # print(' | svt B_n(mu)| =', elems)
        # print()
        # print()

        # for k in range(max(2, len(mu)), n + 1):
        #     #c = sqrt_b_lambda(k, mu)
        #     print(' n =', k, '      G =', schur_expansion(G(k, mu)).set_variable(0, 1))
        # print()

        # nn = n + 3
        # for k in range(1, nn):
        #     print(' n =', k, '      g =', schur_expansion(g(k, mu).set_variable(0, -1)).set_variable(0, -1))
        # print()
        # for k in range(1, nn):
        #     print(' n =', k, '      j =', schur_expansion(j(k, mu).set_variable(0, 1)).set_variable(0, -1))
        # print()
        nn = n + 3
        for k in range(1, nn):
            f = g(k, mu).set_variable(0, -1)
            f = subs(f, k)
            print(' n =', k, ' g(x+1) =', schur_expansion(f))
        print()
        # for k in range(1, nn):
        #     print(' n =', k, ' j(x+1) =', schur_expansion(subs(j(k, mu).set_variable(0, 1), k)).set_variable(0, -1))
        # print()
        # for k in range(1, nn):
        #     print(' n =', k, ' s(x+1) =', schur_expansion(subs(schur(k, mu), k)).set_variable(0, -1))
        # print()


def sqrt_b_lambda(n, mu):
    return _sqrt_b_lambda(n, mu, False)


def dual_sqrt_b_lambda(n, mu):
    if (n, mu) not in dual_blambda_cache:
        dual_blambda_cache[n, mu] = _sqrt_b_lambda(n, mu, True)
    return dual_blambda_cache[n, mu]


def _sqrt_b_lambda(n, mu, dual):
    b = InfiniteCrystal.dual_sqrt_binfty(n) if dual else InfiniteCrystal.sqrt_binfty(n)
    r = InfiniteCrystal.sqrt_r_lambda(mu, n)
    top = InfiniteCrystal.sqrt_tensor(r, b)
    tnaive = top.finitize(0)
    #t = sqrt_r_tensor(mu, n, b)
    #assert len(t) == len(tnaive)
    #return t
    return tnaive

def _test_sqrt_r_lambda(n, k, test_lasc, dual, words=None):
    gp = set()
    ngp = set()

    lp = set()
    nlp = set()

    if words is None:
        words = Permutation.longest_element(n).get_reduced_words()
    for w in words:
        gp.add(w)
        lp.add(w)

        b = InfiniteCrystal.dual_sqrt_binfty(n, w) if dual else InfiniteCrystal.sqrt_binfty(n, w)
        for mu in Partition.all(k, max_row=n):
            print(n, ':', 'mu =', mu)
            
            t = sqrt_r_tensor(mu, n, b)

            r = InfiniteCrystal.sqrt_r_lambda(mu, n)
            top = InfiniteCrystal.sqrt_tensor(r, b)
            tnaive = top.finitize(0)
            
            if len(t) != len(tnaive):
                #t.draw()
                #tnaive.draw()
                print('t =', len(t), 'tnaive =', len(tnaive))
                input('')
            
            # t = top.finitize(0)
            g = t.character()
            yield (w, mu, t)
            try:
                # u = AbstractGLCrystal.sqrtcrystal_from_partition(mu, n)
                # assert AbstractGLCrystal.find_isomorphism(t, u) is not None
                ch = G_expansion_no_beta(SymmetricPolynomial.from_polynomial(g))
            except:
                ch = None
                gp = gp - {w}
                ngp |= {(w, sum(mu))}

                print()
                print('  ', w, ':', len(t))
                print()
                print('  ', 'not symmetric')
                print()
                
                #yield (w, mu, t)

                if test_lasc:
                    lasc = decompose_into_lascoux(g)
                    # print('  ', lasc)
                    # print()
                    try:
                        assert min(set(lasc.values())) == 1
                    except:
                        # print('  ', 'not lascoux positive:', w)
                        if w in lp:
                            lp.remove(w)
                            nlp.add(w)
                        break
                else:
                    pass
                    break
            
            if ch is not None:
                try:
                    assert ch == Vector({mu: 1})
                except:
                    print()
                    print('ch =', ch)
                    print()

        print()
        print('groth works:')
        for a in gp:
            print('  ', a)
        print('groth fails:')
        for a in ngp:
            print('  ', a)
        print()

        if test_lasc:
            print()
            print('lasc works:')
            for a in lp:
                print('  ', a)
            print('lasc fails:')
            for a in nlp:
                print('  ', a)
            print()

        print('\n\n[continue]\n\n')


def test_odd_sqrt_r_lambda(n=3, k=10):
    b = InfiniteCrystal.odd_sqrt_binfty(n)
    for mu in Partition.all(k, max_row=n):
        print(n, ':', 'mu =', mu)
        r = InfiniteCrystal.sqrt_r_lambda(mu, n)
        top = InfiniteCrystal.sqrt_tensor(r, b)
        t = top.finitize(0)
        #u = AbstractGLCrystal.sqrtcrystal_from_partition(mu, n)
        #assert AbstractGLCrystal.find_isomorphism(t, u) is not None
        try:
            ch = G_expansion_no_beta(SymmetricPolynomial.from_polynomial(t.character()))
            assert ch == Vector({mu: 1})
        except:
            f = G(n, mu).polynomial().set_variable(0, 1)
            g = t.character()
            print('  ', f)
            print()
            print('  ', g)
            print()
            print('  ', g - f)
            print()
            top.draw_thresh(0)
            # return top
            input('')


def test_simple_sqrt_r_lambda(n=3, k=10):
    b = InfiniteCrystal.simple_sqrt_binfty(n)
    for mu in Partition.all(k, max_row=n):
        print(n, ':', 'mu =', mu)
        r = InfiniteCrystal.sqrt_r_lambda(mu, n)
        top = InfiniteCrystal.sqrt_tensor(r, b)
        t = top.finitize(0)
        try:
            ch = G_expansion_no_beta(SymmetricPolynomial.from_polynomial(t.character()))
            assert ch == Vector({mu: 1})
        except:
            f = G(n, mu).polynomial().set_variable(0, 1)
            g = t.character()
            print('  ', f)
            print()
            print('  ', g)
            print()
            print('  ', g - f)
            print()
            top.draw_thresh(0)
            # return top
            input('')


def test_demazure_braids(n=4, h=8):
    words = Permutation.longest_element(n).get_reduced_words()
    y = InfiniteCrystal.sqrt_tableaux(n)
    for thresh in range(h + 1):
        boolean = True
        for a in words:
            for b in words:
                if a < b:    
                    boolean &= y.is_demazure_isomorphic(thresh, a, b)
        print('rank =', n, 'threshhold =', thresh, boolean)


def embed_q(n, mu, tab, extended=False):
    if tab is None:
        return None
    try:
        rows = tab.get_rows(unpack=False)
    except TypeError:
        rows = tab.get_rows()
    while len(rows) < len(mu):
        rows.append([])
    for i, m in enumerate(mu):
        rows[i] = m * (n - i,) + tuple(rows[i])
        if extended and len(rows[i]) == m:
            rows[i] = (m - 1) * (n - i,) + (i - n,)
    return tab.from_rows(rows, shifted=True)


def embed_qplus(n, mu, tab):
    return embed_q(n, mu, tab, True)


def test_embed_q(n=3, k=10, cls=AbstractQCrystal):
    embed = embed_q if cls is AbstractQCrystal else embed_qplus
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    partitions = list(filter(Partition.is_strict, partitions))
    for index, lam in enumerate(partitions):
        b = cls.decomposition_tableaux_from_strict_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    try:
                        if u is not None:
                            assert embed(n, mu, u) == b.e_operator_on_decomposition_tableaux(i, embed(n, mu, t))
                    except:
                        print('mu =', mu, 'e_{%s}' % i)
                        print(t)
                        print(u)
                        print(embed(n, mu, t))
                        print(embed(n, mu, u))
                        print('yet:')
                        print(b.e_operator_on_decomposition_tableaux(i, embed(n, mu, t)))
                        return False

                    try:
                        u = b.f_operator(i, t)
                        if u is not None:
                            assert embed(n, mu, u) == b.f_operator_on_decomposition_tableaux(i, embed(n, mu, t))
                    except:
                        print('mu =', mu, 'f_{%s}' % i)
                        print(t)
                        print(u)
                        print(embed(n, mu, t))
                        print(embed(n, mu, u))
                        print('yet:')
                        print(b.f_operator_on_decomposition_tableaux(i, embed(n, mu, t)))
                        return False
    return True


def test_embed_qplus(n=2, k=10):
    # fails
    assert not test_embed_q(n, k, AbstractPrimedQCrystal)


def test_embed_sqrt_q(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    partitions = list(filter(Partition.is_strict, partitions))
    for index, lam in enumerate(partitions):
        b = AbstractQCrystal.sqrtcrystal_from_strict_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_q(n, mu, u) == embed_q(n, mu, t).sqrt_decomposition_e_operator(i)
                    
                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_q(n, mu, u) == embed_q(n, mu, t).sqrt_decomposition_f_operator(i)


def embed_gl(mu, tab):
    if tab is None:
        return None
    rows = tab.get_rows(unpack=False)
    while len(rows) < len(mu):
        rows.append([])
    for i, m in enumerate(mu):
        rows[i] = m * (i + 1,) + tuple(rows[i])
    return Tableau.from_rows(rows)


def test_embed_gl(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    for index, lam in enumerate(partitions):
        b = AbstractGLCrystal.from_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == b.e_operator_on_semistandard_tableaux(i, embed_gl(mu, t))

                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == b.f_operator_on_semistandard_tableaux(i, embed_gl(mu, t))


def test_embed_sqrt_gl(n=3, k=10):
    partitions = sorted({Partition.transpose(mu) for mu in Partition.all(k, max_part=n)})
    for index, lam in enumerate(partitions):
        b = AbstractGLCrystal.sqrtcrystal_from_partition(lam, n)
        print('n =', n, 'lam =', lam, b.indices, '#', index, 'of', len(partitions))
        for mu in partitions:
            # print('* mu =', mu)
            for i in b.indices:
                for t in b:
                    u = b.e_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == embed_gl(mu, t).sqrt_e_operator(i)
                    
                    u = b.f_operator(i, t)
                    if u is not None:
                        assert embed_gl(mu, u) == embed_gl(mu, t).sqrt_f_operator(i)


