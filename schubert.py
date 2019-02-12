

SCHUBERT_CACHE = {}
DOUBLE_SCHUBERT_CACHE = {}
FPF_SCHUBERT_CACHE = {}
INV_SCHUBERT_CACHE = {}
GROTHENDIECK_CACHE = {}
FPF_GROTHENDIECK_CACHE = {}


class HashableDict(dict):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))


class MPolynomial:

    """
    MPolynomial
    -----------

    Attributes:
     coeffs

    Methods:
     Constructor - takes dictionary whose keys are HashableDicts with integers key/val pairs
                   representing a monomial
     monomial - takes two ints an input: index then power
     is_zero
     nnz
     total_degree
     coefficient - returns MPolynomial which is coefficient of monomial of input (index,power)

    Overloaded Operators:
     + * ** [] () == !=


    """

    def __init__(self, coeffs={}):
        self.coeffs = coeffs

    @staticmethod
    def monomial(index, power=1):
        if power == 0:
            return MPolynomial({HashableDict({}): 1})

        ind = HashableDict({index: power})
        return MPolynomial({ind: 1})

    def coefficient(self, index, power):
        x = MPolynomial()
        for term in self.coeffs:
            if (index in term and term[index] == power) or (power == 0 and not (index in term)):
                new_term = HashableDict(term.copy())
                if power != 0:
                    del new_term[index]
                x = x + MPolynomial({new_term: self.coeffs[term]})
        return x

    def total_degree(self):
        ans = None
        for ind in self.coeffs:
            d = 0
            for i in ind:
                d = d + ind[i]

            if ans is None:
                ans = d
            else:
                ans = max(ans, d)

        return ans

    def __getitem__(self, i):
        i = HashableDict(i)
        if i in self.coeffs:
            return self.coeffs[i]
        return 0

    def __call__(self, x):
        ans = 0
        for ind in self.coeffs:
            factor = self.coeffs[ind]
            for j in ind:
                factor = factor * x.get(j, self.monomial(j, 1))**ind[j]
            ans = ans + factor
        return ans

    def __eq__(self, other):
        return (self - other).nnz() == 0

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if isinstance(other, int):
            other = other * MPolynomial.monomial(0, 0)

        newcoeffs = self.coeffs.copy()
        for i in other.coeffs:
            newcoeffs[i] = self[i] + other[i]
            if newcoeffs[i] == 0:
                del newcoeffs[i]

        return MPolynomial(newcoeffs)

    __radd__ = __add__

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -(self - other)

    def isobaric_divided_difference(self, i):
        return (self * MPolynomial.monomial(i, 1)).divided_difference(i)

    @classmethod
    def divided_difference_helper(cls, i, index, coeff):
        a = index.get(i, 0)
        b = index.get(i + 1, 0)
        d = max(a, b) - min(a, b)

        x = MPolynomial.monomial(i, 1)
        y = MPolynomial.monomial(i + 1, 1)
        ans = MPolynomial()

        new_index = HashableDict(index.copy())
        new_index[i] = min(a, b)
        new_index[i + 1] = min(a, b)
        tmp = MPolynomial({new_index: coeff})

        if a == max(a, b):
            sgn = 1
        else:
            sgn = -1
        for i in range(d):
            ans += sgn * tmp * x**i * y**(d - 1 - i)
        return ans

    def divided_difference(self, i):
        ans = MPolynomial()
        for index, coeff in self.coeffs.items():
            ans += self.divided_difference_helper(i, index, coeff)
        return ans

    def __mul__(self, f):
        if type(f) == int:
            return self * MPolynomial({HashableDict({}): f})
        newcoeffs = {}
        for i in self.coeffs:
            for j in f.coeffs:
                k = {t: i.get(t, 0) + j.get(t, 0) for t in i.keys() | j.keys()}
                k = HashableDict({t: k[t] for t in k if k[t] != 0})
                if k in newcoeffs:
                    newcoeffs[k] = newcoeffs[k] + self[i] * f[j]
                else:
                    newcoeffs[k] = self[i] * f[j]
                if newcoeffs[k] == 0:
                    del newcoeffs[k]
        return MPolynomial(newcoeffs)

    __rmul__ = __mul__

    def __pow__(self, i):
        if i == 0:
            return MPolynomial.monomial(0, 0)
        if i < 0:
            if len(self.coeffs) == 1:
                new_coeffs = {}
                for ind in self.coeffs:
                    if abs(self.coeffs[ind]) > 1:
                        return None
                    new_ind = ind.copy()
                    for key in new_ind:
                        new_ind[key] *= i
                    new_coeffs[HashableDict(new_ind)] = self.coeffs[ind]
                return MPolynomial(new_coeffs)
            return None
        return self * (self**(i - 1))

    def nnz(self):
        nonzeros = 0
        for i in self.coeffs:
            if self[i] != 0:
                nonzeros += 1
        return nonzeros

    def is_zero(self):
        return self.nnz() == 0

    def constant_term(self):
        return self[HashableDict({})]

    @staticmethod
    def letters(i):
        # if i == 0:
        #     return "x"
        # if i == 1:
        #     return "y"
        # if i == 2:
        #     return "z"
        # if i == 3:
        #     return "w"
        # if i == 4:
        #     return "v"
        # if i == 5:
        #     return "u"
        # if i == 6:
        #     return "t"
        # if i == 7:
        #     return "s"
        # if i == 8:
        #     return "r"
        if i >= 0:
            return "x_" + str(i)
        else:
            return "y_" + str(-i)

    @staticmethod
    def index_to_str(ind):

        s = ''
        for i in ind:
            if ind[i] != 0:
                s = s + ' ' + MPolynomial.letters(i)
                if ind[i] != 1:
                    s = s + "^" + str(ind[i])
        # s = '(' + s[1:] + ')'
        s = s[1:]
        if s == "()":
            s = ""
        return s

    def __repr__(self):
        if self.nnz() == 0:
            return '0'
        s = ''
        filtered = filter(lambda x: self[x] != 0, self.coeffs)

        def sorter(index):
            ans = []
            for i in sorted(index):
                ans += abs(index[i]) * [i]
            return ans

        for i in sorted(filtered, key=sorter):
            monomial = MPolynomial.index_to_str(i)
            coeff = str(abs(self[i]))

            if coeff == "1" and monomial != "":
                coeff = ""

            if self[i] < 0:
                coeff = " - " + coeff
            else:
                coeff = " + " + coeff

            s = s + coeff + monomial
        s = s[1:]
        if s[0] == "+":
            s = s[2:]
        else:
            s = "-" + s[2:]

        return s

    def __hash__(self):
        return hash(str(self))


def x(i):
    return MPolynomial.monomial(i)


def y(i):
    return MPolynomial.monomial(-i)


def one():
    return x(1)**0


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
    def get(cls, w):
        assert cls.is_valid(w)
        oneline = tuple(w.oneline)
        cache = cls.cache()
        if oneline not in cache:
            v, i = cls.get_ascent(w)
            if i is not None:
                s = cls.get(v)
                s = cls.divided_difference(s, i)
            else:
                s = cls.top(len(oneline))
            cache[oneline] = s
        return cache[oneline]


class Schubert(AbstractSchubert):

    @classmethod
    def cache(cls):
        return SCHUBERT_CACHE

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


class Grothendieck(Schubert):

    @classmethod
    def cache(cls):
        return GROTHENDIECK_CACHE

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 - MPolynomial.monomial(i + 1, 1))).divided_difference(i)


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


class FPFGrothendieck(FPFSchubert):
    @classmethod
    def cache(cls):
        return FPF_GROTHENDIECK_CACHE

    @classmethod
    def top(cls, n):
        s = one()
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1 - i):
                s *= (x(i) + x(j) - x(i) * x(j))
        return s

    @classmethod
    def product(cls, w):
        s = one()
        for i, j in w.fpf_rothe_diagram():
            s *= (x(i) + x(j) - x(i) * x(j))
        return s

    @classmethod
    def divided_difference(cls, f, i):
        return (f * (1 - MPolynomial.monomial(i + 1, 1))).divided_difference(i)


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
