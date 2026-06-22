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
from stable.utils import G, G_expansion_no_beta, SymmetricPolynomial
from stable.vectors import Vector


def test_lp(n, thresh=10):
    lp = InfiniteCrystal.sqrt_LP(n)
    elems = lp.vertices(thresh, thresh)

    def floor(x):
        return 2 * (x // 2)

    def ceil(x):
        return 2 * ((x + 1) // 2)

    def get(i, x):
        return 0 if i == n else x[i - 1]

    def sub(i, x):
        x = list(x)
        x[i - 1] -= 1
        return tuple(x)

    def add(i, x):
        x = list(x)
        x[i - 1] += 1
        return tuple(x)

    for i in range(1, n):
        for v in elems:
            # test e
            ev = lp.e_operator(i, v)
            case = ''
            if get(i, v) >= ceil(get(i + 1, v)):
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

            # test eps

            # test phi

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
    gp = set()
    ngp = set()

    lp = set()
    nlp = set()

    for w in Permutation.longest_element(n).get_reduced_words():
        gp.add(w)
        lp.add(w)

        b = InfiniteCrystal.sqrt_binfty(n, w)
        for mu in Partition.all(k, max_row=n):
            print(n, ':', 'mu =', mu)
            
            t = sqrt_r_tensor(mu, n, b)

            r = InfiniteCrystal.sqrt_r_lambda(mu, n)
            top = InfiniteCrystal.sqrt_tensor(r, b)
            tnaive = top.finitize(0)
            
            if len(t) != len(tnaive):
                t.draw()
                tnaive.draw()
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
                ngp |= {w}

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
                    #break
            
            if ch is not None:
                assert ch == Vector({mu: 1})

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


