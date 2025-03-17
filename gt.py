from cached import cached_value
from schubert import X
from permutations import Permutation
from tableaux import Tableau
import itertools

GT_CACHE = {}
STRICT_GT_CACHE = {}


def get_monomial(weight):
    ans = X(0) ** 0
    for i, w in enumerate(weight):
        ans *= X(i + 1) ** w
    return ans 


def ordinary_weight(rows):
    w = []
    for i in range(len(rows)): 
        diff = sum(map(abs, rows[-i - 1]))
        if i > 0:
            diff -= sum(map(abs, rows[-i]))
        w.append(diff)
    ans = get_monomial(w)


def strict_weight(rows):
    w = []
    for i in range(len(rows)): 
        diff = sum(map(abs, rows[-i - 1]))
        if i > 0:
            diff -= sum(map(abs, rows[-i]))
        w.append(diff)
    ans = get_monomial(w)
    # for i in range(1, len(rows)):
    #     for j in range(len(rows[i])):
    #         a, b, c = rows[i - 1][j], rows[i][j], rows[i - 1][j + 1]
    #         if a < b < c:
    #             ans *= 2
    return ans


class GTPattern:

    def __init__(self, array=(), compute_weight=None):
        assert all(len(row) == len(array) - i for i, row in enumerate(array))
        
        self.rows = tuple(tuple(row) for row in array)
        self._string = None
        self._length = len(self.rows)
        self._dkp = None
        self._skp = None
        
        if compute_weight is None:
            compute_weight = ordinary_weight
        self._weight = compute_weight(self.rows)

    def __str__(self):
        if self._string is None:
            space = ' '
            pad = max([0] + [len(str(x)) for row in self.rows for x in row])
            gap = (pad + 2) * space
            ans = []
            for i in range(len(self)):
                s = i * gap + gap.join([space + (pad - len(str(x))) * space + str(x) + space for x in self.rows[i]])
                ans.append(s[1:])
            self._string = '\n' + '\n'.join(ans) + '\n'
        return self._string

    @classmethod
    def from_tableau(cls, t, n=None):
        n = len(t.weight()) if n is None else n
        rows = []
        for i in range(n):
            mu = list(t.find(*[j for j in range(1, i + 2)]).partition().tuple())
            assert len(mu) <= i + 1
            while len(mu) < i + 1:
                mu.append(0)
            mu = list(reversed(mu))
            rows.append(mu)
        return GTPattern(list(reversed(rows)))

    def tableau(self):
        def shape(mu):
            return {(i + 1, j + 1) for i in range(len(mu)) for j in range(mu[i])}

        def diff(i):
            ans = shape(list(reversed(self.rows[i])))
            if i + 1 < len(self):
                ans -= shape(list(reversed(self.rows[i + 1])))
            return ans

        boxes = {}
        n = len(self)
        for i in range(n - 1, -1, -1):
            for (a, b) in diff(i):
                boxes[a, b] = n - i
        return Tableau(boxes)

    def shifted_tableau(self):
        pass

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(self.rows)

    def __eq__(self, other):
        return self.weight == other.weight and self.rows == other.rows

    def __len__(self):
        return self._length

    @classmethod
    def constant(cls, n, val):
        return cls([i * [val] for i in range(n, 0, -1)])

    @property
    def length(self):
        return len(self)

    @property
    def weight(self):
        return self._weight

    def __getitem__(self, *args):
        (i, j) = args if len(args) == 2 else (args[0][0], args[0][1])
        return self.get(i, j)

    def get(self, i, j):
        return self.rows[i][j]

    def strict_kogan_permutation(self):
        def dbl(a):
            return 2 * a if a >= 0 else -1 - 2 * a

        def cmp(a, b, c, d):
            return dbl(a) >= dbl(b) + dbl(d) and (c is None or dbl(a) > dbl(c) + dbl(d))

        if self._skp is None:
            n = len(self)
            ans = Permutation()
            for i in range(1, n):
                for j in range(i - 1, -1, -1):
                    v = n - (i - j)
                    d = self.get(0, v) - self.get(0, v - 1) - 1
                    if cmp(self.get(-i, j), abs(self.get(-i - 1, j)), self.get(-i, j - 1) if j > 0 else None, d):
                        ans = ans % Permutation.s_i(n - 1 - j)
            self._skp = ans
        return self._skp

    def dual_kogan_permutation(self):
        if self._dkp is None:
            n = len(self)
            ans = Permutation()
            for i in range(1, n):
                for j in range(i - 1, -1, -1):
                    if self.get(-i, j) == self.get(-i - 1, j + 1):
                        ans = ans % Permutation.s_i(n - 1 - j)
            self._dkp = ans
        return self._dkp

    def is_strict(self):
        return all(self[i, j - 1] < self[i, j] for i in range(len(self)) for j in range(1, len(self.rows[i])))

    @classmethod
    def strict_from_partition(cls, mu, size=None):
        if size is not None:
            assert size >= len(mu)
            mu = list(mu) + (size - len(mu)) * [0]
        for i in range(len(mu)):
            mu[i] += len(mu) - i - 1
        return cls.all_strict(mu)

    @classmethod
    def from_partition(cls, mu, size=None):
        if size is not None:
            assert size >= len(mu)
            mu = list(mu) + (size - len(mu)) * [0]
        return cls.all(mu)

    @classmethod
    def all(cls, *args):
        if len(args) == 1 and type(args[0]) in [list, tuple]:
            args = args[0]
        
        assert all(type(args[i]) == int for i in range(len(args)))
        assert all(args[i] >= 0 for i in range(len(args)))
        
        args = tuple(sorted(args))
        return cls._all(args)

    @cached_value(GT_CACHE)
    def _all(cls, mu):
        if len(mu) == 0:
            return [GTPattern()]
        if len(mu) == 1:
            return [GTPattern([mu])]
        ans = []
        intervals = [range(mu[i - 1], mu[i] + 1) for i in range(1, len(mu))]
        for nu in itertools.product(*intervals):
            ans.extend([GTPattern((mu,) + x.rows) for x in cls._all(nu)])
        return ans

    @classmethod
    def all_strict(cls, *args):
        if len(args) == 1 and type(args[0]) in [list, tuple]:
            args = args[0]
        
        assert all(type(args[i]) == int for i in range(len(args)))
        assert all(args[i] >= 0 for i in range(len(args)))
        
        args = tuple(sorted(args))
        assert all(args[i] < args[i + 1] for i in range(len(args) - 1))
        return cls._all_strict(args)

    @cached_value(STRICT_GT_CACHE)
    def _all_strict(cls, mu):
        if len(mu) == 0:
            return [GTPattern()]
        if len(mu) == 1:
            return [GTPattern([mu])]

        intervals = []
        for i in range(1, len(mu)):
            a = abs(mu[i - 1])
            b = abs(mu[i])
            val = list(range(a, b + 1)) + list(range(-b + 1, -a))
            intervals.append(val)
        ans = []
        for nu in itertools.product(*intervals):
            if all(abs(nu[i]) < abs(nu[i + 1]) for i in range(len(nu) - 1)):
                ans.extend([GTPattern((mu,) + x.rows) for x in cls._all_strict(nu)])
        return ans

