from vectors import Vector
from collections import defaultdict
import itertools


class Word:
    def __init__(self, subset=None, *args):
        self.subset = subset or set()
        self.elements = tuple(args)
        assert all(i in subset for i in args)

    def __iter__(self):
        return self.elements.__iter__()

    def __hash__(self):
        return hash((self.elements, tuple(sorted(self.subset))))

    def __getitem__(self, i):
        return self.elements[i]

    def __len__(self):
        return len(self.elements)

    def __add__(self, other):
        if type(other) == Word:
            return Vector({self: 1}) + Vector({other: 1})
        elif type(other) == Vector:
            return Vector({self: 1}) + other
        assert type(other) in [Word, Vector]

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == Word:
            return Vector({self: 1}) + Vector({other: -1})
        elif type(other) == Vector:
            return Vector({self: 1}) - other
        assert type(other) in [Word, Vector]

    def __or__(self, other):
        assert type(other) == Word
        return Word(self.subset | other.subset, *(self.elements + other.elements))

    # def __lshift__(self, i):
    #     return Word(*[e - i for e in self.elements])

    # def __rshift__(self, i):
    #     return Word(*[e + i for e in self.elements])

    def __eq__(self, other):
        assert type(other) == Word
        return self.subset == other.subset and self.elements == other.elements

    def __lt__(self, other):
        assert type(other) == Word
        return self.elements < other.elements

    def __repr__(self):
        return '(%s | %s)' % (', '.join(map(str, self.elements)), str(self.subset))

    def _shuffle(self, other, indices):
        a, b = list(reversed(self.elements)), list(reversed(other.elements))
        word = []
        for i in range(len(self) + len(other)):
            if i in indices:
                word.append(a.pop())
            else:
                word.append(b.pop())
        return Word(self.subset | other.subset, *word)

    def __mul__(self, other):
        if type(other) == Word:
            assert self.subset.isdisjoint(other.subset)
            dictionary = defaultdict(int)
            for indices in itertools.combinations(set(range(len(self) + len(other))), len(self)):
                word = self._shuffle(other, indices)
                dictionary[word] += 1
            return Vector(dictionary)
        elif type(other) == int:
            return Vector({self: other})
        elif type(other) == Vector:
            return Vector({self: 1}) * other
        assert type(other) in [int, Word, Vector]

    def __rmul__(self, other):
        assert type(other) in [int, Word, Vector]
        return self.__mul__(other)

    def coproduct(self, *subsets):
        # check that subsets are ordered partition of self.subset
        assert self.subset == {i for x in subsets for i in x}
        assert {len(x & y) for x in subsets for y in subsets if x != y}.issubset({0})

        subwords = [tuple(i for i in self if i in a) for a in subsets]
        if tuple(i for word in subwords for i in word) == self.elements:
            return Vector.base(tuple(Word(subsets[i], *subwords[i]) for i in range(len(subsets))))
        else:
            return Vector()


REDUCED_WORDS = {(): {()}}


def reduce_oneline(oneline):
    while oneline and oneline[-1] == len(oneline):
        oneline = oneline[:-1]
    return oneline


def get_reduced_words(oneline):
    oneline = reduce_oneline(oneline)
    if oneline not in REDUCED_WORDS:
        words = set()
        for i in range(len(oneline) - 1):
            if oneline[i] > oneline[i + 1]:
                a, b = oneline[i:i + 2]
                newline = oneline[:i] + (b, a) + oneline[i + 2:]
                words |= {w + (i + 1,) for w in get_reduced_words(newline)}
        REDUCED_WORDS[oneline] = words
    return REDUCED_WORDS[oneline]


class Permutation:
    def __init__(self, *args, **kwargs):
        assert set(args) == set(range(1, len(args) + 1))
        s = kwargs.get('subset', None)
        if s is None:
            s = set(range(1, len(args)))
        self.vector = Vector({Word(s, *w): 1 for w in get_reduced_words(args)})
        self.oneline = tuple(args)

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield Permutation(*args)

    @classmethod
    def oneline_from_word(cls, word, n):
        line = list(range(1, n + 1))
        for i in word:
            temp = line[i]
            line[i] = line[i - 1]
            line[i - 1] = temp
        return tuple(line)

    def __eq__(self, other):
        assert type(other) == Permutation
        return self.oneline == other.oneline

    def __hash__(self):
        return hash(self.oneline)

    def _right_shift(self, i):
        assert i >= 0
        subset = {t + i for t in range(1, self.size)}
        oneline = list(range(1, i + 1)) + [n + i for n in self.oneline]
        return Permutation(*oneline, subset=subset)

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
        if type(other) == Permutation:
            return Vector({self: 1}) + Vector({other: 1})
        elif type(other) == Vector:
            return Vector({self: 1}) + other
        assert type(other) in [Permutation, Vector]

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == Permutation:
            return Vector({self: 1}) + Vector({other: -1})
        elif type(other) == Vector:
            return Vector({self: 1}) - other
        assert type(other) in [Permutation, Vector]

    def __mul__(self, other):
        if type(other) == Permutation:
            assert self.size >= 1 and other.size >= 1
            n = self.size + other.size - 1
            result = self.vector * other._right_shift(self.size - 1).vector
            answer = Vector()
            while result:
                key, value = next(iter(result.items()))
                sigma = Permutation(*self.oneline_from_word(key, n))
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
        return Permutation(*newline).exclude(*args)

    def startswith(self, sequence):
        return all(self(i + 1) == sequence[i] for i in range(len(sequence)))

    def endswith(self, sequence):
        n = len(sequence)
        return all(self(self.size - i) == sequence[n - 1 - i] for i in range(n))

    @classmethod
    def reverse(cls, n):
        oneline = list(reversed(range(1, n + 1)))
        return Permutation(*oneline)

    def __call__(self, i):
        if i < 1 or i > len(self.oneline):
            return i
        return self.oneline[i - 1]
