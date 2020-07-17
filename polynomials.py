from permutations import Permutation


DIVIDED_DIFFERENCE_CACHE = {}


def q(i):
    return MPolynomial.monomial(0xFFFFFFFFFFFFFF)**i


def q_coeff(f, degree):
    return f.coefficient(0xFFFFFFFFFFFFFF, degree)


def X(i):
    return MPolynomial.monomial(i)


def D(i):
    return Operator.create(i) * (1 + 10 * X(i + 1))


class HashableDict(dict):

    def __hash__(self):
        return hash(tuple(sorted(self.items())))

    def __lt__(self, other):
        return tuple(sorted(self.items())) < tuple(sorted(other.items()))


class Operator:

    class Monomial:

        def __init__(self, index):
            if type(index) == int:
                index = (index,)
            assert type(index) == tuple
            assert all(type(i) == int for i in index)
            self.index = Permutation.from_word(index).get_reduced_word()

        def __repr__(self):
            return ' '.join(['D_%s' % i for i in self.index]) if self.index else 'I'

        def __hash__(self):
            return hash(self.index)

        def __eq__(self, other):
            assert type(other) == type(self)
            return self.index == other.index

        def __lt__(self, other):
            assert type(other) == type(self)
            return self.index < other.index

        def matches(self, other):
            return self.index and other.index and self.index[-1] == other.index[0]

    def __init__(self):
        self.dictionary = {}

    def __repr__(self):
        s = [
            '%s' % k if v == 1 else
            '(%s) * %s' % (v, k)
            for k, v in self.items()
        ]
        return ' + '.join(s) if s else '0'

    def __eq__(self, other):
        if other == 0:
            return len(self.dictionary) == 0
        assert type(other) == Operator
        return (self - other) == 0

    @classmethod
    def create(cls, *args):
        if len(args) == 1 and type(args[0]) in [tuple, list]:
            args = tuple(args[0])
        if len(args) == 1 and type(args[0]) == Operator.Monomial:
            args = args[0].index
        assert all(type(i) == int for i in args)
        ans = Operator()
        if not any(args[i] == args[i + 1] for i in range(len(args) - 1)):
            ans.dictionary[Operator.Monomial(tuple(args))] = 1
        return ans

    def __getitem__(self, item):
        return self.dictionary.get(item, 0)

    def keys(self):
        return self.dictionary.keys()

    def items(self):
        return self.dictionary.items()

    def __iter__(self):
        return self.dictionary.__iter__()

    def __add__(self, other):
        assert type(other) == Operator
        dictionary = {i: self[i] + other[i] for i in self.keys() | other.keys()}
        dictionary = {i: v for i, v in dictionary.items() if v}
        ans = Operator()
        ans.dictionary = dictionary
        return ans

    def __sub__(self, other):
        assert type(other) == Operator
        dictionary = {i: self[i] - other[i] for i in self.keys() | other.keys()}
        dictionary = {i: v for i, v in dictionary.items() if v}
        ans = Operator()
        ans.dictionary = dictionary
        return ans

    def __mul__(self, other):
        if type(other) == Operator:
            ans = Operator()
            for i in self:
                for j in other:
                    term = self[i] * Operator.create(i) * other[j]
                    dictionary = {Operator.Monomial(k.index + j.index): term[k] for k in term if not k.matches(j)}
                    term = Operator()
                    term.dictionary = dictionary
                    ans += term
            return ans

        if type(other) == int:
            other = MPolynomial.one() * other

        if type(other) in [MPolynomial]:
            queue = [(len(i.index) + 1, self[i],) + i.index + (other,) for i in self]
            ans = Operator()
            while queue:
                tup, queue = queue[0], queue[1:]
                i, tup = tup[0], tup[1:]
                if i == 0:
                    ans += tup[0] * Operator.create(tup[1:])
                if i == 1:
                    ans += (tup[0] * tup[1]) * Operator.create(tup[2:])
                else:
                    j = tup[i - 1]
                    one = (i - 1,) + tup[:i - 1] + (tup[i].divided_difference(j),) + tup[i + 1:]
                    two = (i - 1,) + tup[:i - 1] + (tup[i].toggle(j), j) + tup[i + 1:]
                    queue += [one, two]
            return ans

        raise Exception

    def __rmul__(self, other):
        if type(other) == Operator:
            other.__mul__(self)
        if type(other) in [int, MPolynomial]:
            ans = Operator()
            ans.dictionary = {i: other * v for i, v in self.items() if other * v}
            return ans
        raise Exception


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
    def any(self):
        assert len(self.coeffs) > 0
        return next(iter(self.coeffs))

    def variables(self):
        return {i for key in self.coeffs.keys() for i in key}

    def max_variable(self):
        v = self.variables()
        if v:
            return max(v)
        return 0

    def truncate(self, nvar):
        return MPolynomial({m: v for m, v in self.coeffs.items() if not any(i > nvar for i in m)})

    def __bool__(self):
        return not self.is_zero()

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

    def degree(self):
        return self.total_degree()

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

    def __iter__(self):
        return self.coeffs.__iter__()

    @classmethod
    def one(cls):
        return cls.monomial(1, 0)

    def set(self, i, e):
        return self.substitute(i, e)

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

    def __lt__(self, other):
        other = MPolynomial.one() * other if type(other) == int else other
        return all(v > 0 for v in (other - self).coeffs.values())

    def __gt__(self, other):
        other = MPolynomial.one() * other if type(other) == int else other
        return all(v > 0 for v in (self - other).coeffs.values())

    def is_positive(self):
        return self > 0

    def is_not_laurent_polynomial(self):
        return not any(v < 0 for c in self.coeffs for v in c.values())

    def isobaric_divided_difference(self, i):
        return (self * MPolynomial.monomial(i, 1)).divided_difference(i)

    def toggle(self, i):
        ans = MPolynomial()
        for index, coeff in self.coeffs.items():
            new_index = HashableDict(index.copy())
            if i + 1 in index:
                new_index[i] = index[i + 1]
            elif i in new_index:
                del new_index[i]
            if i in index:
                new_index[i + 1] = index[i]
            elif i + 1 in new_index:
                del new_index[i + 1]
            ans += MPolynomial({new_index: coeff})
        return ans

    @classmethod
    def divided_difference_helper(cls, i, index):
        if (i, index) not in DIVIDED_DIFFERENCE_CACHE:
            a = index.get(i, 0)
            b = index.get(i + 1, 0)
            d = max(a, b) - min(a, b)

            x = MPolynomial.monomial(i, 1)
            y = MPolynomial.monomial(i + 1, 1)
            ans = MPolynomial()

            new_index = HashableDict(index.copy())
            new_index[i] = min(a, b)
            new_index[i + 1] = min(a, b)
            tmp = MPolynomial({new_index: 1})

            sgn = 1 if a == max(a, b) else -1
            for j in range(d):
                ans += sgn * tmp * x**j * y**(d - 1 - j)
            DIVIDED_DIFFERENCE_CACHE[(i, index)] = ans
            ell = len(DIVIDED_DIFFERENCE_CACHE)
            # if ell % 100 == 0:
            #    print(' . . . Divided Differences cache:', ell)
        return DIVIDED_DIFFERENCE_CACHE[(i, index)]

    def divided_difference(self, i):
        ans = MPolynomial()
        for index, coeff in self.coeffs.items():
            ans += self.divided_difference_helper(i, index) * coeff
        return ans

    def __mul__(self, f):
        if type(f) == int:
            return self * MPolynomial({HashableDict({}): f})
        if type(f) != MPolynomial:
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
        return MPolynomial(newcoeffs)

    def __rmul__(self, other):
        if type(other) not in [MPolynomial, int]:
            return other.__mul__(self)
        else:
            return self.__mul__(other)

    def __floordiv__(self, other):
        assert type(other) in [int]
        coeffs = {m: c // other for (m, c) in self.coeffs.items() if c / other}
        return MPolynomial(coeffs)

    def __pow__(self, i):
        if i == 0:
            return MPolynomial.monomial(0, 0)
        if i == 1:
            return self
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
        return self**(i // 2) * self**(i - i // 2)

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
        if i == 0xFFFFFFFFFFFFFF:
            return 'q'
        if i > 0:
            return "x_" + str(i)
        if i == 0:
            return "\u03B2"
        if i < 0:
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
            c = 0
            for i in sorted(index):
                c += index[i]
                ans += abs(index[i]) * [-i]
            return (c,) + tuple(ans)

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
