from crystals import(
    AbstractGLCrystal,
    SuperGLCrystal,
)
from stable.symmetric import SymmetricPolynomial
from stable.utils import hs, HG


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

            