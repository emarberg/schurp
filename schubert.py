from vectors import Vector
from permutations import Permutation
from polynomials import (
    MPolynomial,
    one,
    x,
    y,
    X,
)
from tableaux import Tableau


SCHUBERT_CACHE = {}
DOUBLE_SCHUBERT_CACHE = {}
FPF_SCHUBERT_CACHE = {}
INV_SCHUBERT_CACHE = {}
GROTHENDIECK_CACHE = {}
FPF_GROTHENDIECK_CACHE = {}
INV_GROTHENDIECK_CACHE = {}
ALT_INV_GROTHENDIECK_CACHE = {}


class AbstractSchubert(object):

    @classmethod
    def cache(cls):
        raise NotImplementedError

    @classmethod
    def top(cls, w):
        raise NotImplementedError

    @classmethod
    def invalid_case(cls):
        raise Exception

    @classmethod
    def is_valid(cls, w):
        return True

    @classmethod
    def get_ascent(cls, w):
        raise NotImplementedError

    @classmethod
    def divided_difference(cls, f, i):
        return f.divided_difference(i)

    @classmethod
    def reduce(cls, oneline):
        return tuple(oneline)

    @classmethod
    def get(cls, w, verbose=False):
        assert cls.is_valid(w)
        oneline = cls.reduce(w.oneline)
        cache = cls.cache()
        if oneline not in cache:
            v, i = cls.get_ascent(Permutation(*oneline))
            if i is not None:
                s = cls.get(v)
                s = cls.divided_difference(s, i)
            elif v is not None:
                s = cls.top(v)
            else:
                s = cls.invalid_case()
            cache[oneline] = s
            if verbose:
                print(' . . .', cls.__name__, 'cache:', len(cache))
        return cache[oneline]

    @classmethod
    def from_code(cls, code):
        raise NotImplementedError

    @classmethod
    def least_term(cls, f):
        return min(f, key=lambda s: tuple(-i for k in sorted(s) for i in s[k] * [k]))

    @classmethod
    def decompose(cls, f):
        ans = Vector()
        while not f.is_zero():
            m = cls.least_term(f)
            code = max(m) * [0] if m else []
            for i, v in m.items():
                code[i - 1] = v
            code = tuple(code)
            w = cls.from_code(code)
            ans += Vector({w: f[m]})
            f = f - f[m] * cls.get(w)
        return ans


class Schubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return SCHUBERT_CACHE

    @classmethod
    def reduce(cls, oneline):
        while oneline and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        return tuple(oneline)

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if n == 0:
            return w.s_i(1), 1
        if w.is_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        return w * w.s_i(i), i

    @classmethod
    def from_code(cls, code):
        return Permutation.from_code(code)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.rothe_diagram():
            s *= x(i)
        return s


class Grothendieck(Schubert):

    beta = 1

    @classmethod
    def cache(cls):
        return GROTHENDIECK_CACHE

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)


class DoubleSchubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return DOUBLE_SCHUBERT_CACHE

    @classmethod
    def top(cls, w):
        s = one()
        for i, j in w.rothe_diagram():
                s *= (x(i) - y(j))
        return s

    @classmethod
    def get_ascent(cls, w):
        return Schubert.get_ascent(w)


class FPFSchubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return FPF_SCHUBERT_CACHE

    @classmethod
    def reduce(cls, oneline):
        while len(oneline) >= 2 and oneline[-2] == len(oneline) and oneline[-1] == len(oneline) - 1:
            oneline = oneline[:-2]
        return tuple(oneline)

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.fpf_rothe_diagram():
            s *= x(i) + x(j)
        return s

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_fpf_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        s = w.s_i(i)
        return s * w * s, i

    @classmethod
    def is_valid(cls, w):
        return all(len(c) == 2 for c in w.cycles)

    @classmethod
    def from_code(cls, code):
        return Permutation.from_fpf_involution_code(code)


class FPFGrothendieck(FPFSchubert):

    beta = 1

    @classmethod
    def cache(cls):
        return FPF_GROTHENDIECK_CACHE

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.fpf_rothe_diagram():
            s *= (x(i) + x(j) + cls.beta * x(i) * x(j))
        return s

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)


class InvSchubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return INV_SCHUBERT_CACHE

    @classmethod
    def top(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= (x(i) + x(j)) if i != j else x(i)
        return s


    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_dominant():
            return w, None
        i = min(set(range(1, n)) - w.right_descent_set)
        s = w.s_i(i)
        if s * w == w * s:
            return w * s, i
        else:
            return s * w * s, i

    @classmethod
    def is_valid(cls, w):
        return w.is_involution()

    @classmethod
    def from_code(cls, code):
        return Permutation.from_involution_code(code)


    @classmethod
    def vexillary_tableau_formula(cls, w):
        path = w.get_essential_path()
        
        t = [0, 0]
        for i in range(1, len(path)):
            a1, b1 = path[i - 1]
            a2, b2 = path[i]
            assert a1 == a2 or b1 == b2
            if a1 == a2:
                assert b2 == b1 - 1
                t.append(-X(b1))
            elif b1 == b2:
                assert a2 == a1 + 1
                t.append(X(a2))
        max_entry = path[0][0]
        mu = w.involution_shape()

        ans = 0
        for tab in Tableau.get_semistandard_shifted(mu, n=max_entry, diagonal_primes=True):
            summand = X(0)**0
            for (i, j) in tab.mapping:
                a = abs(tab.get(i, j))
                b = j - i + 1
                s = -1 if tab.get(i, j).is_primed() else 1
                summand *= X(a) + s * t[b]
            ans += summand
        return ans

class InvGrothendieck(InvSchubert):

    beta = 1

    @classmethod
    def grothendieck_involution_words(cls, w):
        assert w.is_vexillary()
        f = cls.get(w)
        dec = Grothendieck.decompose(f)
        for u in dec:
            for word in u.get_reduced_words():
                yield word, dec[u]

    @classmethod
    def cache(cls):
        return INV_GROTHENDIECK_CACHE

    @classmethod
    def invalid_case(cls):
        return MPolynomial.zero()

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if w.is_dominant():
            return w, None
        if not w.is_vexillary():
            return None, None
        for i in set(range(1, n)) - w.right_descent_set:
            s = w.s_i(i)
            ws = w * s if  s * w == w * s else s * w *s
            if ws.is_vexillary():
                return ws, i
        return None, None

    @classmethod
    def top(cls, w):
        return cls.product(w)

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= x(i) + x(j) + cls.beta * x(i) * x(j)
        return s

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 + cls.beta * MPolynomial.monomial(i + 1))).divided_difference(i)

    @classmethod
    def diagonal_product(cls, n):
        s = one()
        for i in range(1, n + 1):
            s *= (x(i) + x(i) + cls.beta * x(i) * x(i))
        return s


class AltInvGrothendieck(InvGrothendieck):

    @classmethod
    def cache(cls):
        return ALT_INV_GROTHENDIECK_CACHE

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.involution_rothe_diagram():
            s *= (x(i) + x(j) + cls.beta * x(i) * x(j)) if i != j else (x(i) + x(j))
        return s
