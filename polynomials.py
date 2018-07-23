

def gcd(*args):
    if len(args) == 0:
        return 1
    if len(args) == 1:
        return args[0]
    if len(args) == 2:
        a, b = tuple(args)
        a = abs(a)
        b = abs(b)
        if b == 0:
            return a
        return gcd(b, a % b)
    else:
        return gcd(args[0], gcd(args[1:]))


class Polynomial:
    """
    Polynomial
    ----------

    Attributes:
     coeffs

    Methods:
     Constructor - takes dictionary of integer key/value pairs
     monomial
     divide - self and d must be polynomials;
              returns [q=pol,r=pol,c=nonzero int] such that (c*self) = q*d + r
     is_invertible
     is_zero
     nnz
     degree
     leading_coefficient


    Overloaded Operators:
     + * ** [] () == != %

    """

    def __hash__(self):
        return hash(str(self))

    def __init__(self, coeffs={}):
        self.coeffs = coeffs

    @staticmethod
    def monomial(i):
        return Polynomial({i: 1})

    @staticmethod
    def gcd(f, g):
        if g == 0:
            return f
        return Polynomial.gcd(g, f.divide(g)[1])

    def divide(self, pol):
        if self.is_zero():
            return (Polynomial(), Polynomial(), 1)
        if type(pol) == int:
            pol = Polynomial.monomial(0) * pol
        if pol == 0 or min(pol.coeffs) < 0 or min(self.coeffs) < 0:
            return None
        if self == 0:
            return (Polynomial(), Polynomial(), 1)
        remainder = self
        quotient = Polynomial()

        b = pol.leading_coefficient()
        c = 1

        x = Polynomial.monomial(1)
        e = pol.degree()
        while remainder.degree() >= e and remainder != 0:
            deg = remainder.degree()
            a = remainder[deg]
            d = gcd(a, b)
            c *= b // d
            quotient *= b // d
            remainder *= b // d
            quotient += a // d * x**(deg - e)
            remainder += -a // d * x**(deg - e) * pol
        if c < 0:
            return (-quotient, -remainder, -c)
        else:
            return (quotient, remainder, c)

    def __mod__(self, m):
        newcoeffs = {}
        for i in self.coeffs:
            c = self.coeffs[i] % m
            if c != 0:
                newcoeffs[i] = c
        return Polynomial(newcoeffs)

    def __getitem__(self, i):
        if i in self.coeffs:
            return self.coeffs[i]
        return 0

    def __call__(self, x):
        return sum(map(lambda i: x**i * self[i], self.coeffs))

    def __eq__(self, other):
        if other is None:
            return False
        return (self - other).nnz() == 0

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if isinstance(other, int):
            other = other * Polynomial.monomial(0)

        newcoeffs = self.coeffs.copy()
        for i in other.coeffs:
            newcoeffs[i] = self[i] + other[i]
            if newcoeffs[i] == 0:
                del newcoeffs[i]

        return Polynomial(newcoeffs)

    __radd__ = __add__

    def __neg__(self):
        return self * (-1)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -(self - other)

    def __div__(self, int):
        newcoeffs = {}
        for i in self.coeffs:
            c = self.coeffs[i] / int
            if c != 0:
                newcoeffs[i] = c
        return Polynomial(newcoeffs)

    def __mul__(self, f):
        if type(f) == int:
            return self * Polynomial({0: f})

        newcoeffs = {}
        for i in self.coeffs:
            for j in f.coeffs:
                if i + j in newcoeffs:
                    newcoeffs[i + j] = newcoeffs[i + j] + self[i] * f[j]
                else:
                    newcoeffs[i + j] = self[i] * f[j]

                if newcoeffs[i + j] == 0:
                    del newcoeffs[i + j]
        return Polynomial(newcoeffs)

    __rmul__ = __mul__

    def __pow__(self, i):
        if i == 0:
            return Polynomial.monomial(0)
        if i < 0:
            if self.is_invertible():
                return Polynomial({-self.coeffs.keys()[0]: int(self.coeffs.values()[0]**-1)})**-i
            return None

        return self * (self**(i - 1))

    def degree(self):
        if self.is_zero():
            return 0
        return max(max(self.coeffs), -min(self.coeffs))

    def leading_coefficient(self):
        if self.is_zero():
            return 0
        return self.coeffs[max(self.coeffs)]

    def nnz(self):
        nonzeros = 0
        for i in self.coeffs:
            if self[i] != 0:
                nonzeros += 1
        return nonzeros

    def is_zero(self):
        return self.nnz() == 0

    def repr_helper(self):
        coeffs = lambda x: (' + ')*(x > 0) + (' - ')*(x < 0) + (str(abs(x)) + '*')*(x != 0 and abs(x) !=1)
        exps = lambda x: ('q' + ('**' + str(x))*(x != 1))*(x != 0)
        next = lambda i: (coeffs(self[i]) + exps(i)) * (self[i] != 0)

        s = ''
        for i in filter(lambda x: self[x] != 0, self.coeffs):
            if i == 0:
                x = self[0]
                t = str(abs(x))
                s = s + (' + ') * (x > 0) + (' - ') * (x < 0) + t
            else:
                s = s + next(i)
        return s[1] * (s[1] == '-') + s[3:]

    def __repr__(self):
        if self.nnz() == 0:
            return '0'
        return self.repr_helper()


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
