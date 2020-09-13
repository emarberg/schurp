from vectors import Vector
from signed import SignedPermutation
from even import EvenSignedPermutation
from words import Word, get_reduced_words
import itertools


class HopfPermutation:
    def __init__(self, *args, **kwargs):
        assert set(args) == set(range(1, len(args) + 1))
        s = kwargs.get('subset', None)
        if s is None:
            s = set(range(1, len(args)))
        self._subset = s
        self._vector = None
        self.oneline = tuple(args)

    def is_fully_commutative(self):
        w = self.oneline
        n = len(w)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if w[i] > w[j] > w[k]:
                        return False
        return True

    @property
    def vector(self):
        if self._vector is None:
            self._vector = Vector({
                Word(*w, subset=self._subset): 1
                for w in get_reduced_words(self.oneline)
            })
        return self._vector

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield HopfPermutation(*args)

    @classmethod
    def involutions(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            w = HopfPermutation(*args)
            if all(w(w(i)) == i for i in range(1, n + 1)):
                yield w

    @classmethod
    def fpf_involutions(cls, n):
        for w in cls.involutions(n):
            if all(w(i) != i for i in range(1, n + 1)):
                yield w

    @classmethod
    def oneline_from_word(cls, word, n):
        line = list(range(1, n + 1))
        for i in word:
            temp = line[i]
            line[i] = line[i - 1]
            line[i - 1] = temp
        return tuple(line)

    def __eq__(self, other):
        assert type(other) == HopfPermutation
        return self.oneline == other.oneline

    def __hash__(self):
        return hash(self.oneline)

    def _right_shift(self, i):
        assert i >= 0
        subset = {t + i for t in range(1, self.size)}
        oneline = list(range(1, i + 1)) + [n + i for n in self.oneline]
        return HopfPermutation(*oneline, subset=subset)

    @property
    def size(self):
        return len(self.oneline)

    def __len__(self):
        return len(next(iter(self.vector)))

    def find(self, i):
        for j in range(1, 1 + len(self.oneline)):
            if self(j) == i:
                return j
        return i

    def __repr__(self):
        return ''.join(str(i) for i in self.oneline)

    def __add__(self, other):
        if type(other) == HopfPermutation:
            return Vector({self: 1}) + Vector({other: 1})
        elif type(other) == Vector:
            return Vector({self: 1}) + other
        assert type(other) in [HopfPermutation, Vector]

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == HopfPermutation:
            return Vector({self: 1}) + Vector({other: -1})
        elif type(other) == Vector:
            return Vector({self: 1}) - other
        assert type(other) in [HopfPermutation, Vector]

    def coproduct(self, *subsets):
        # check that subsets are ordered partition of self.subset
        assert set(range(1, len(self.oneline))) == {i for x in subsets for i in x}
        assert {len(x & y) for x in subsets for y in subsets if x != y}.issubset({0})

        def st(subset, word):
            n = len(subset)
            sort = sorted(subset)
            mapping = {sort[i]: i + 1 for i in range(n)}
            newword = tuple(mapping[w] for w in word)
            return self.oneline_from_word(newword, n + 1)

        answer = set()
        for word in self.vector:
            for tup in word.coproduct(*subsets):
                new = tuple(st(subsets[i], tup[i]) for i in range(len(tup)))
                answer.add(new)
        return Vector({tup: 1 for tup in answer})

    @classmethod
    def test_product(cls, u, v):
        answer = Vector()
        for w in cls.test_product_helper(u.oneline, v.oneline):
            answer += Vector({HopfPermutation(*w): 1})
        return answer

    @classmethod
    def test_product_helper(cls, u, v):
        m = len(u) - 1
        n = len(v) - 1
        w = tuple(i + m for i in v)

        if u[-1] == m + 1:
            return {u[:-1] + w}
        if w[0] == m + 1:
            return {u + w[1:]}

        j = [i for i in range(m + 1) if u[i] == m + 1][0]
        k = [i for i in range(n + 1) if w[i] == m + 1][0]
        u_ = flatten(u[j + 1:])
        v_ = flatten(v[:k])

        ans = {merge(u[:j + 1], b, True) for b in cls.test_product_helper(u_, v)}
        ans |= {merge(w[k:], b, False) for b in cls.test_product_helper(u, v_)}
        return ans

    def __mul__(self, other):
        if type(other) == HopfPermutation:
            assert self.size >= 1 and other.size >= 1
            n = self.size + other.size - 1
            result = self.vector * other._right_shift(self.size - 1).vector
            answer = Vector()
            while result:
                key, value = next(iter(result.items()))
                sigma = HopfPermutation(*self.oneline_from_word(key, n))
                answer += Vector({sigma: value})
                result -= sigma.vector * value
            return answer

    def exclude(self, *args):
        if len(args) == 0:
            return self
        args = set(args)
        i = args.pop()
        assert i > 0
        args = {j - (j > i) for j in args}
        newline = tuple(j - (j > i) for j in self.oneline if j != i)
        return HopfPermutation(*newline).exclude(*args)

    def startswith(self, sequence):
        return all(self(i + 1) == sequence[i] for i in range(len(sequence)))

    def endswith(self, sequence):
        n = len(sequence)
        return all(self(self.size - i) == sequence[n - 1 - i] for i in range(n))

    @classmethod
    def reverse(cls, n):
        oneline = list(reversed(range(1, n + 1)))
        return HopfPermutation(*oneline)

    def __call__(self, i):
        if i < 1 or i > len(self.oneline):
            return i
        return self.oneline[i - 1]


class HopfSignedPermutation:
    def __init__(self, *args):
        assert len(args) > 0
        n = max([abs(a) for a in args])
        s = set(range(n))
        assert {i + 1 for i in range(n)} | {-i - 1 for i in range(n)} == {a for a in args} | {-a for a in args}
        self._subset = s
        self._vector = None
        self.oneline = tuple(args)

    @property
    def vector(self):
        if self._vector is None:
            self._vector = Vector({
                Word(*w, subset=self._subset): 1
                for w in SignedPermutation(*self.oneline).get_reduced_words()
            })
        return self._vector

    @classmethod
    def all(cls, n):
        for args in EvenSignedPermutation.all(n):
            yield HopfSignedPermutation(*args.oneline)

    def __eq__(self, other):
        assert type(other) == HopfSignedPermutation
        return self.oneline == other.oneline

    def __hash__(self):
        return hash(self.oneline)

    @property
    def size(self):
        return len(self.oneline)

    def __len__(self):
        return len(next(iter(self.vector)))

    def __repr__(self):
        return ''.join(str(-i) + '\u0305' if i < 0 else str(i) for i in self.oneline)

    def __add__(self, other):
        if type(other) == HopfSignedPermutation:
            return Vector({self: 1}) + Vector({other: 1})
        elif type(other) == Vector:
            return Vector({self: 1}) + other
        assert type(other) in [HopfSignedPermutation, Vector]

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == HopfSignedPermutation:
            return Vector({self: 1}) + Vector({other: -1})
        elif type(other) == Vector:
            return Vector({self: 1}) - other
        assert type(other) in [HopfSignedPermutation, Vector]

    @classmethod
    def oneline_from_word(cls, word, n):
        line = list(range(1, n + 1))
        for i in word:
            if i == 0:
                line[0] *= -1
            else:
                temp = line[i]
                line[i] = line[i - 1]
                line[i - 1] = temp
        return tuple(line)

    def __mul__(self, other):
        if type(other) == HopfPermutation:
            assert self.size >= 1 and other.size >= 1
            n = self.size + other.size - 1
            result = self.vector * other._right_shift(self.size - 1).vector
            answer = Vector()
            while result:
                key, value = next(iter(result.items()))
                sigma = HopfSignedPermutation(*self.oneline_from_word(key, n))
                answer += Vector({sigma: value})
                result -= sigma.vector * value
            return answer

    def __call__(self, i):
        if i < 1 or i > len(self.oneline):
            return i
        return self.oneline[i - 1]


class HopfEvenSignedPermutation:
    def __init__(self, *args):
        assert len(args) > 0
        n = max([abs(a) for a in args])
        s = set(range(n))
        assert {i + 1 for i in range(n)} | {-i - 1 for i in range(n)} == {a for a in args} | {-a for a in args}
        assert len([a for a in args if a < 0]) % 2 == 0
        self._subset = s
        self._vector = None
        self.oneline = tuple(args)

    @property
    def vector(self):
        if self._vector is None:
            self._vector = Vector({
                Word(*w, subset=self._subset): 1
                for w in EvenSignedPermutation(*self.oneline).get_reduced_words()
            })
        return self._vector

    @classmethod
    def all(cls, n):
        for args in EvenSignedPermutation.all(n):
            yield HopfEvenSignedPermutation(*args.oneline)

    def __eq__(self, other):
        assert type(other) == HopfEvenSignedPermutation
        return self.oneline == other.oneline

    def __hash__(self):
        return hash(self.oneline)

    @property
    def size(self):
        return len(self.oneline)

    def __len__(self):
        return len(next(iter(self.vector)))

    def __repr__(self):
        return ''.join(str(-i) + '\u0305' if i < 0 else str(i) for i in self.oneline)

    def __add__(self, other):
        if type(other) == HopfEvenSignedPermutation:
            return Vector({self: 1}) + Vector({other: 1})
        elif type(other) == Vector:
            return Vector({self: 1}) + other
        assert type(other) in [HopfEvenSignedPermutation, Vector]

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == HopfEvenSignedPermutation:
            return Vector({self: 1}) + Vector({other: -1})
        elif type(other) == Vector:
            return Vector({self: 1}) - other
        assert type(other) in [HopfEvenSignedPermutation, Vector]

    @classmethod
    def oneline_from_word(cls, word, n):
        line = list(range(1, n + 1))
        for i in word:
            if i == 0:
                line[0:2] = [-line[1], -line[0]]
            else:
                temp = line[i]
                line[i] = line[i - 1]
                line[i - 1] = temp
        return tuple(line)

    def __mul__(self, other):
        if type(other) == HopfPermutation:
            assert self.size >= 2 and other.size >= 1
            n = self.size + other.size - 1
            result = self.vector * other._right_shift(self.size - 1).vector
            answer = Vector()
            while result:
                key, value = next(iter(result.items()))
                sigma = HopfEvenSignedPermutation(*self.oneline_from_word(key, n))
                answer += Vector({sigma: value})
                result -= sigma.vector * value
            return answer

    def __call__(self, i):
        if i < 1 or i > len(self.oneline):
            return i
        return self.oneline[i - 1]
