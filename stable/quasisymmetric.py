from .polynomials import Polynomial, X, HashableDict
from .cached import cached_value
from .vectors import Vector
import itertools


MONOMIAL_CACHE = {}
FUNDAMENTAL_CACHE = {}
MULTIFUNDAMENTAL_CACHE = {}


class Quasisymmetric:

    @classmethod
    def monomial_expansion(cls, f):
        return cls.expansion(f, cls.monomial)

    @classmethod
    def fundamental_expansion(cls, f):
        return cls.expansion(f, cls.fundamental)

    @classmethod
    def multifundamental_expansion(cls, f):
        return cls.expansion(f, cls.multifundamental)

    @classmethod
    def from_expansion(cls, nvars, vec, creator):
        ans = Polynomial()
        for alpha, coeff in vec.items():
            ans += creator(nvars, alpha) * coeff
        return ans

    @classmethod
    def expansion(cls, f, creator):
        ans = Vector(sorter=sum)
        if f != 0:
            nvars = cls.nvars(f)
        while f != 0:
            alpha = cls.leading_term(f)
            coeff = cls.coefficient(f, alpha)
            g = creator(nvars, alpha)
            assert cls.coefficient(g, alpha) == 1
            ans += Vector({alpha: coeff}, sorter=sum)
            f = f - g * coeff
        return Vector(sorter=sum)

    @classmethod
    def nvars(cls, f):
        ans = 0
        for alpha in f.coeffs:
            ans = max(ans, max(alpha, default=0))
        return ans

    @classmethod
    def leading_term(cls, f):
        assert not f.is_zero()

        def key(alpha):
            tup = []
            for i in sorted(alpha):
                if i > 0:
                    tup += alpha[i] * (i,)
            return (len(tup), tuple(tup))

        alpha_dict = sorted(f.coeffs, key=key)[0]
        alpha = ()
        for i in sorted(alpha_dict):
            if i > 0:
                alpha += (alpha_dict[i],)
                assert i == len(alpha)
        return alpha

    @classmethod
    def coefficient(cls, f, alpha):
        ans = Polynomial()
        for i, c in f.coeffs.items():
            dictionary = {j: v for j, v in i.items() if j > 0}
            comparison = {i + 1: a for (i, a) in enumerate(alpha)}
            if dictionary == comparison:
                newdict = {j: v for j, v in i.items() if j <= 0}
                ans += Polynomial({HashableDict(newdict): c})
        return ans

    @cached_value(MONOMIAL_CACHE)
    def monomial(cls, nvars, alpha):
        assert all(type(a) == int and a > 0 for a in alpha)
        k = len(alpha)
        if k == 0:
            return Polynomial.one()
        inds = list(range(1, nvars + 1))
        ans = Polynomial()
        for subset in itertools.combinations(inds, k):
            term = Polynomial.one()
            for i, a in enumerate(subset):
                term *= X(a)**alpha[i]
            ans += term
        return ans

    @classmethod
    def composition_descents(cls, alpha):
        ans = []
        for a in alpha[:-1]:
            ans.append(a + (ans[-1] if ans else 0))   
        return ans   

    @cached_value(FUNDAMENTAL_CACHE)
    def fundamental(cls, nvars, alpha):
        assert all(type(a) == int and a > 0 for a in alpha)
        k = len(alpha)
        if k == 0:
            return Polynomial.one()

        ascents = cls.composition_descents(alpha)

        def indices(lowerbound, termsleft, init):
            if termsleft == 0:
                yield init
                return
            for i in range(lowerbound, nvars + 1):
                newbound = (i + 1) if (len(init) + 1) in ascents else i
                for ans in indices(newbound, termsleft - 1, init + (i,)):
                    yield ans

        ans = Polynomial()
        for seq in indices(1, sum(alpha), ()):
            term = Polynomial.one()
            for a in seq:
                term *= X(a)
            ans += term
        return ans

    @cached_value(MULTIFUNDAMENTAL_CACHE)
    def multifundamental(cls, nvars, alpha):
        assert all(type(a) == int and a > 0 for a in alpha)
        k = len(alpha)
        if k == 0:
            return Polynomial.one()

        ascents = cls.composition_descents(alpha)

        def indices(lowerbound, termsleft, init):
            if termsleft == 0:
                yield init
                return
            available = list(range(lowerbound, nvars + 1))
            for i in range(1, 1 + len(available)):
                for s in itertools.combinations(available, i):
                    newbound = max(s) + (1 if (len(init) + 1) in ascents else 0)
                    for ans in indices(newbound, termsleft - 1, init + (s,)):
                        yield ans

        ans = Polynomial()
        for seq in indices(1, sum(alpha), ()):
            term = Polynomial.one()
            for a in seq:
                for i in a:
                    term *= X(i)
            ans += term * X(0)**(term.total_degree() - sum(alpha))
        return ans