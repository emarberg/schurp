

class HashableDict(dict):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def degree(self):
        d = 0
        for i in self:
            d = d + self[i]
        return d


def X(i=None):
    if i is None:
        i = Polynomial.MAX_INT
    assert i >= 0
    return Polynomial.x(i)


def Y(i=None):
    if i is None:
        i = Polynomial.MIN_INT
    assert i <= 0
    return Polynomial.x(i)


class Polynomial:

    """
    Polynomial
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
     coefficient - returns Polynomial which is coefficient of monomial of input (index,power)

    Overloaded Operators:
     + * ** [] () == !=


    """

    MAX_INT = 2**64 - 1
    MIN_INT = -MAX_INT

    def factor(self):
        pass

    def __bool__(self):
        return not self.is_zero()

    def __xor__(self, other):
        return self + other + X(0) * self * other

    def __init__(self, coeffs={}):
        self.coeffs = coeffs

    def __mod__(self, o):
        assert type(o) == int
        return Polynomial({k: v % o for k, v in self.coeffs.items() if v % o != 0})

    def __truediv__(self, o):
        assert type(o) == int
        assert all(v % o == 0 for v in self.coeffs.values())
        return Polynomial({k: v // o for k, v in self.coeffs.items()})

    def __floordiv__(self, o):
        assert type(o) == int
        return Polynomial({k: v // o for k, v in self.coeffs.items() if v // o != 0})

    def truncate(self, degree_bound):
        return self.__class__({
            mon: coeff
            for mon, coeff in self.coeffs.items()
            if mon.degree() <= degree_bound
        })

    @classmethod
    def monomial(cls, index, power=1):
        if power == 0:
            return cls({HashableDict({}): 1})

        ind = HashableDict({index: power})
        return cls({ind: 1})

    def coefficient(self, index, power):
        x = Polynomial()
        for term in self.coeffs:
            if (index in term and term[index] == power) or (power == 0 and not (index in term)):
                new_term = HashableDict(term.copy())
                if power != 0:
                    del new_term[index]
                x = x + Polynomial({new_term: self.coeffs[term]})
        return x

    def total_degree(self):
        ans = None
        for ind in self.coeffs:
            d = ind.degree()
            if ans is None:
                ans = d
            else:
                ans = max(ans, d)
        return ans

    def __iter__(self):
        return self.coeffs.__iter__()

    @classmethod
    def zero(cls):
        return cls.one() * 0

    @classmethod
    def one(cls):
        return cls.monomial(1, 0)

    def substitute(self, i, e):
        ans = 0
        for ind in self.coeffs:
            term = self.one() * self.coeffs[ind]
            for j in ind:
                if i != j:
                    term *= self.monomial(j, ind[j])
                else:
                    assert ind[j] >= 0
                    term *= e ** ind[j]
            ans = ans + term
        return ans

    def set(self, i, e):
        return self.substitute(i, e)

    def divide_linear(self, i, c):
        # divide by x(i) + c
        ans = self.substitute(i, x(i) - c) * self.monomial(i, -1)
        return ans.substitute(i, x(i) + c)

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
            other = other * Polynomial.monomial(0, 0)

        newcoeffs = self.coeffs.copy()
        for i in other.coeffs:
            newcoeffs[i] = self[i] + other[i]
            if newcoeffs[i] == 0:
                del newcoeffs[i]

        return Polynomial(newcoeffs)

    __radd__ = __add__

    def __len__(self):
        return len(self.coeffs)

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -(self - other)

    def is_positive(self):
        return self > 0

    def __lt__(self, other):
        other = Polynomial.one() * other if type(other) == int else other
        return all(v > 0 for v in (other - self).coeffs.values())

    def __gt__(self, other):
        other = Polynomial.one() * other if type(other) == int else other
        return all(v > 0 for v in (self - other).coeffs.values())

    def __mul__(self, f):
        if type(f) == int:
            return self * Polynomial({HashableDict({}): f})
        if type(f) != Polynomial:
            return f.__rmul__(self)
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
        return Polynomial(newcoeffs)

    def __rmul__(self, other):
        if type(other) not in [Polynomial, int]:
            return other.__mul__(self)
        else:
            return self.__mul__(other)

    def __pow__(self, i):
        if i == 0:
            return Polynomial.monomial(0, 0)
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
                return Polynomial(new_coeffs)
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

    def is_integer(self):
        return self == self.constant_term()

    @classmethod
    def letters(cls, i):
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
        if i == 0:
            return 'β'
        elif i == cls.MAX_INT:
            return 'x'
        elif i == cls.MIN_INT:
            return 'y'
        elif i > 0:
            return "x_" + str(i)
        else:
            return "y_" + str(-i)

    @classmethod
    def index_to_str(cls, ind):
        s = ''
        for i in ind:
            if ind[i] != 0:
                s = s + ' ' + cls.letters(i)
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

        def sorter(index):
            ans = []
            c = 0
            for i in sorted(index):
                c += index[i]
                ans += abs(index[i]) * [-i]
            return (c,) + tuple(ans)

        s = ''
        filtered = filter(lambda x: self[x] != 0, self.coeffs)
        for i in sorted(filtered, key=sorter):
            monomial = Polynomial.index_to_str(i)
            coeff = str(abs(self[i]))
            if coeff == "1" and monomial != "":
                coeff = ""
            coeff = (" - " if self[i] < 0 else " + ") + coeff
            s = s + coeff + monomial
        s = s[1:]
        s = s[2:] if s[0] == "+" else "-" + s[2:]
        return s

    def __hash__(self):
        return hash(str(self))

    @classmethod
    def x(cls, i):
        return cls.monomial(i)

    @classmethod
    def y(cls, i):
        return cls.monomial(-i)

    @classmethod
    def from_tuple(cls, tup):
        ans = X(0)**0
        for i, a in enumerate(tup):
            ans *= X(i)**a
        return ans


beta = X(0)
