from vectors import Vector
from permutations import Permutation
from polynomials import (
    MPolynomial,
    one,
    x,
    y,
    X,
)


SCHUBERT_CACHE = {}
DOUBLE_SCHUBERT_CACHE = {}
FPF_SCHUBERT_CACHE = {}
INV_SCHUBERT_CACHE = {}
GROTHENDIECK_CACHE = {}
FPF_GROTHENDIECK_CACHE = {}


class AbstractSchubert(object):

    @classmethod
    def cache(cls):
        raise NotImplementedError

    @classmethod
    def top(cls, n):
        raise NotImplementedError

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
    def get(cls, w):
        assert cls.is_valid(w)
        oneline = cls.reduce(w.oneline)
        cache = cls.cache()
        if oneline not in cache:
            v, i = cls.get_ascent(Permutation(*oneline))
            if i is not None:
                s = cls.get(v)
                s = cls.divided_difference(s, i)
            else:
                s = cls.top(len(oneline))
            cache[oneline] = s
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
        if not f.is_zero():
            m = cls.least_term(f)
            code = max(m) * [0] if m else []
            for i, v in m.items():
                code[i - 1] = v
            code = tuple(code)
            w = cls.from_code(code)
            return Vector({w: f[m]}) + cls.decompose(f - f[m] * cls.get(w))
        else:
            return Vector()


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
    def top(cls, n):
        s = one()
        for i in range(1, n):
            s *= MPolynomial.monomial(i, n - i)
        return s

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if n == 0:
            return w.s_i(1), 1
        if len(w.right_descent_set) == n - 1:
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

    beta = X(0)

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
    def top(cls, n):
        s = one()
        for i in range(1, n):
            for j in range(1, n + 1 - i):
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
    def top(cls, n):
        s = one()
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1 - i):
                s *= (x(i) + x(j))
        return s

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if len(w.right_descent_set) in [0, n - 1]:
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

    beta = X(0)

    @classmethod
    def cache(cls):
        return FPF_GROTHENDIECK_CACHE

    @classmethod
    def top(cls, n):
        s = one()
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1 - i):
                s *= (x(i) + x(j) + cls.beta * x(i) * x(j))
        return s

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
    def top(cls, n):
        s = one()
        for i in range(1, n + 1):
            for j in range(i, n + 1 - i):
                if i == j:
                    s *= x(i)
                else:
                    s *= (x(i) + x(j))
        return s

    @classmethod
    def get_ascent(cls, w):
        n = len(w.oneline)
        if len(w.right_descent_set) in [0, n - 1]:
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
