from vectors import Vector
from collections import defaultdict
from itertools import combinations


class Word:
    def __init__(self, *args):
        self.elements = tuple(args)

    def __iter__(self):
        return self.elements.__iter__()

    def __hash__(self):
        return hash(self.elements)

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
        return Word(*(self.elements + other.elements))

    def __lshift__(self, i):
        return Word(*[e - i for e in self.elements])

    def __rshift__(self, i):
        return Word(*[e + i for e in self.elements])

    def __eq__(self, other):
        assert type(other) == Word
        return self.elements == other.elements

    def __lt__(self, other):
        assert type(other) == Word
        return self.elements < other.elements

    def __repr__(self):
        return str(self.elements)

    def _shuffle(self, other, subset):
        a, b = list(reversed(self.elements)), list(reversed(other.elements))
        word = []
        for i in range(len(self) + len(other)):
            if i in subset:
                word.append(a.pop())
            else:
                word.append(b.pop())
        return Word(*word)

    def __mul__(self, other):
        if type(other) == Word:
            dictionary = defaultdict(int)
            for subset in combinations(set(range(len(self) + len(other))), len(self)):
                word = self._shuffle(other, subset)
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


REDUCED_WORDS = {(): {Word()}}


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
                words |= {w | Word(i + 1) for w in get_reduced_words(newline)}
        REDUCED_WORDS[oneline] = words
    return REDUCED_WORDS[oneline]


class Permutation:
    def __init__(self, *args):
        assert set(args) == set(range(1, len(args) + 1))
        self.vector = Vector({w: 1 for w in get_reduced_words(args)})
        self._oneline = None

    @classmethod
    def get_oneline_from_word(cls, word):
        line = list(range(1, max(set(word.elements) | {0}) + 2))
        for i in word:
            temp = line[i]
            line[i] = line[i - 1]
            line[i - 1] = temp
        return tuple(line)

    @property
    def oneline(self):
        if self._oneline is None:
            word = next(iter(self.vector.keys()))
            self._oneline = self.get_oneline_from_word(word)
        return self._oneline

    def __rshift__(self, i):
        assert i >= 0
        return Permutation(*(list(range(1, i + 1)) + [n + i for n in self.oneline]))

    @property
    def size(self):
        return len(self.oneline)

    def __repr__(self):
        return ''.join(str(i) for i in self.oneline)

    def __mul__(self, other):
        if type(other) == Permutation:
            result = self.vector * (other >> (self.size - 1)).vector
            answer = Vector()
            while result:
                print('...%s' % len(result))
                key, value = next(iter(result.items()))
                sigma = Permutation(*self.get_oneline_from_word(key))
                answer += Vector({sigma: value})
                result -= sigma.vector * value
            return answer

