from crystals import(
    AbstractGLCrystal,
    SuperGLCrystal,
)
from stable.symmetric import SymmetricPolynomial
from stable.utils import hs, HG, hs_expansion, HG_expansion
from stable.vectors import Vector


def test_sqrt_super_crystals(m=3, n=3, k=3):
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


def test_super_crystals(m=3, n=3, k=3):
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


def test_hook_schur(m=4, n=4):
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


def test_hook_grothendieck(m=3, n=3):
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

            