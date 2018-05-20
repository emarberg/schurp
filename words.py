from vectors import Vector
from tableaux import Tableau
from numbers import MarkedNumber
from collections import defaultdict
import itertools
import subprocess


class Word:
    def __init__(self, *args, **kwargs):
        self.subset = kwargs.get('subset', None) or set(args)
        self.elements = tuple(args)
        assert all(i in self.subset for i in args)

    def is_increasing(self):
        for i in range(1, len(self)):
            if self[i - 1] >= self[i]:
                return False
        return True

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

    def __lshift__(self, i):
        return Word({x - i for x in self.subset}, *[e - i for e in self.elements])

    def __rshift__(self, i):
        return Word({x + i for x in self.subset}, *[e + i for e in self.elements])

    def __eq__(self, other):
        assert type(other) == Word
        return self.subset == other.subset and self.elements == other.elements

    def __lt__(self, other):
        assert type(other) == Word
        return self.elements < other.elements

    def __repr__(self):
        return str(self.elements)
#        return '(%s | %s)' % (', '.join(map(str, self.elements)), str(self.subset))

    def _shuffle(self, other, indices):
        a, b = list(reversed(self.elements)), list(reversed(other.elements))
        word = []
        for i in range(len(self) + len(other)):
            if i in indices:
                word.append(a.pop())
            else:
                word.append(b.pop())
        return Word(*word, subset=self.subset | other.subset)

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
            return Vector.base(tuple(Word(*subwords[i], subset=subsets[i]) for i in range(len(subsets))))
        else:
            return Vector()

    def hecke_insert(self):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.hecke_insert(MarkedNumber(a))

            v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def shifted_hecke_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.shifted_hecke_insert(MarkedNumber(a), verbose=verbose)

            if column_dir:
                v = MarkedNumber(-i)
            else:
                v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def involution_insert(self, verbose=True):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.involution_insert(MarkedNumber(a), verbose=verbose)

            if column_dir:
                v = MarkedNumber(-i)
            else:
                v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def eg_insert(self):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.eg_insert(MarkedNumber(a))

            v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def fpf_insert(self, verbose=True):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.fpf_insert(MarkedNumber(a), verbose=verbose)

            if column_dir:
                v = MarkedNumber(-i)
            else:
                v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q


def get_insertion_mapping(words):
    assert all(w.is_increasing() for w in words)
    elements = []
    mapping = {}
    for i, w in enumerate(words):
        e = len(elements) + 1
        for v in range(e, e + len(w)):
            mapping[MarkedNumber(v)] = MarkedNumber(i + 1)
            mapping[MarkedNumber(-v)] = MarkedNumber(-i - 1)
        elements += list(w.elements)
    return Word(*elements), mapping


def shifted_hecke_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.shifted_hecke_insert(verbose=False)
    return p, Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})


def involution_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.involution_insert(verbose=False)
    return p, Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})


def fpf_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.fpf_insert(verbose=False)
    return p, Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})


REDUCED_WORDS = {(): {()}}
INVOLUTION_WORDS = {(): {()}}
FPF_INVOLUTION_WORDS = {(): {()}}


HECKE_DATA = {}


def hecke_product(expr):
    if len(expr) == 0:
        return tuple()
    n = max(expr) + 1
    oneline = list(range(1, n + 1))
    for i in expr:
        if oneline[i - 1] < oneline[i]:
            oneline[i - 1], oneline[i] = oneline[i], oneline[i - 1]
    return reduce_oneline(tuple(oneline))


def split_hecke_product(w, n):
    expr = get_reduced_word(w)
    if len(expr) == 0:
        return tuple(), tuple()
    top = [i - n + 1 for i in expr if i >= n]
    bot = [i for i in expr if i < n]
    return hecke_product(bot), hecke_product(top)


def hecke_word_product(u, v):
    m, n = len(u), len(v)
    u, v = reduce_oneline(u), reduce_oneline(v)
    ans = set()
    for w in itertools.permutations(range(1, m + n)):
        if split_hecke_product(w, m) == (u, v):
            ans.add(w)
    return ans


def reduce_oneline(oneline):
    while oneline and oneline[-1] == len(oneline):
        oneline = oneline[:-1]
    return oneline


def get_reduced_word(oneline):
    return next(iter(get_reduced_words(oneline)))


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


def reduce_fpf(oneline):
    oneline = reduce_oneline(oneline)
    n = len(oneline)
    assert all(oneline[i - 1] != i == oneline[oneline[i - 1] - 1] for i in range(1, n + 1))
    while oneline and oneline[-1] == len(oneline) - 1 and oneline[-2] == len(oneline):
        oneline = oneline[:-2]
    return oneline


def get_fpf_involution_words(oneline):
    oneline = reduce_fpf(oneline)
    if oneline not in FPF_INVOLUTION_WORDS:
        words = set()
        for i in range(len(oneline) - 1):
            a, b = oneline[i:i + 2]
            if a > b and (a, b) != (i + 2, i + 1):
                newline = list(oneline[:])
                newline[a - 1] = i + 2
                newline[b - 1] = i + 1
                newline = newline[:i] + [newline[i + 1], newline[i]] + newline[i + 2:]
                newline = tuple(newline)
                words |= {w + (i + 1,) for w in get_fpf_involution_words(newline)}
        FPF_INVOLUTION_WORDS[oneline] = words
    return FPF_INVOLUTION_WORDS[oneline]


def get_involution_words(oneline):
    oneline = reduce_oneline(oneline)
    if oneline not in INVOLUTION_WORDS:
        words = set()
        for i in range(len(oneline) - 1):
            if oneline[i] > oneline[i + 1]:
                a, b = oneline[i:i + 2]
                if (a, b) != (i + 2, i + 1):
                    newline = list(oneline[:])
                    newline[a - 1] = i + 2
                    newline[b - 1] = i + 1
                    newline = newline[:i] + [newline[i + 1], newline[i]] + newline[i + 2:]
                    newline = tuple(newline)
                else:
                    newline = oneline[:i] + (b, a) + oneline[i + 2:]
                words |= {w + (i + 1,) for w in get_involution_words(newline)}
        INVOLUTION_WORDS[oneline] = words
    return INVOLUTION_WORDS[oneline]


def flatten(word):
    n = len(word)
    letters = sorted(word)
    new = []
    for a in word:
        for i in range(n):
            if letters[i] == a:
                new.append(i + 1)
                letters[i] = None
                break
    return tuple(new)


def merge(a, b, a_is_first=True):
    assert len(set(a)) == len(a)
    assert len(set(b)) == len(b)
    n = len(a) + len(b)
    mapping = dict(zip(sorted(b), [i for i in range(1, n + 1) if i not in a]))
    c = tuple(mapping[i] for i in b)
    if a_is_first:
        return tuple(a + c)
    else:
        return tuple(c + a)


def fpf_from_word(*w):
    if len(w) == 0:
        return tuple()
    n = max(w) + 1
    if n % 2 != 0:
        n += 1

    def conj(oneline, i):
        a, b = oneline[i - 1:i + 1]
        newline = oneline[:i - 1] + (b, a) + oneline[i + 1:]
        j = [t for t in range(len(newline)) if newline[t] == i][0]
        k = [t for t in range(len(newline)) if newline[t] == i + 1][0]
        newline = list(newline)
        newline[j] = i + 1
        newline[k] = i
        return tuple(newline)

    start = list(range(1, n + 1))
    for i in range(n):
        if i % 2 == 0:
            start[i] += 1
        else:
            start[i] -= 1
    start = tuple(start)

    for i in w:
        start = conj(start, i)
    return start


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


class ShiftedCrystalGenerator:

    DIRECTORY = '/Users/emarberg/Dropbox/shifted_crystals/examples/inv/'

    @classmethod
    def test_insertion_tableaux(cls, n, k):
        for i, w in enumerate(HopfPermutation.involutions(n)):
            cg = ShiftedCrystalGenerator(w.oneline, k)
            shapes = [
                {cg.insertion_tableau(i) for i in comp}
                for comp in cg.components
            ]
            b1 = all(len(sh) == 1 for sh in shapes)
            b2 = all(len(s & t) == 0 for s in shapes for t in shapes if s != t)
            print(i, w, b1, b2, '(n = %s, k = %s)' % (n, k))
            assert b1 and b2

    def print_zero_weight_space(self):
        for f in self.zero_weight_space:
            p, q = involution_insert(*f)
            print(p)
            print('')
            print(q)
            print('')
            print(f)
            print('\n')

    @property
    def _filename(self):
        w = ''.join([str(i) for i in self.oneline])
        return 'size%s_%s_%s' % (str(len(self._factorizations)), w, str(self.num_factors))

    @property
    def dot_filename(self):
        return self.DIRECTORY + 'dot/' + '%s.dot' % self._filename

    @property
    def png_filename(self):
        return self.DIRECTORY + 'png/' + '%s.png' % self._filename

    def node_label(self, i):
        pre = str(self.insertion_tableau(i))
        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(j) for j in w.elements]) for w in self[i]])
        return pre + '\n\n' + top + '\n\n' + bottom

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        s += ['    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), i) for (x, y, i) in self.edges]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    @classmethod
    def all(cls, n, k):
        for w in HopfPermutation.involutions(n):
            if w.oneline != reduce_oneline(w.oneline):
                continue
            cg = ShiftedCrystalGenerator(w.oneline, k)
            if cg.edges:
                cg.generate()

    def __init__(self, oneline, k):
        self.oneline = oneline
        self.words = self._get_words()
        self.num_factors = k
        self._factorizations = None
        self._weights = None
        self._ranks = None
        self._edges = None

    def _get_words(self):
        return get_involution_words(self.oneline)

    def __repr__(self):
        return 'Crystal generator for w = %s with ell = %s' % (self.oneline, self.num_factors)

    @property
    def factorizations(self):
        if self._factorizations is None:
            fac = [
                tuple([Word(*i) for i in tup])
                for w in self.words
                for tup in self.get_increasing_factorizations(w, self.num_factors)
            ]
            self._factorizations = sorted(fac, key=lambda f: tuple(len(a) for a in f), reverse=True)
        return self._factorizations

    @property
    def highest_weight(self):
        return self.weights[0]

    @property
    def weights(self):
        if self._weights is None:
            self._weights = [tuple(len(a) for a in f) for f in self.factorizations]
        return self._weights

    @property
    def zero_weight_space(self):
        return [f for f in self.factorizations if all(len(p) == 1 for p in f)]

    @property
    def ranks(self):
        if self._ranks is None:
            self._ranks = [self.get_rank(i) for i in range(len(self))]
        return self._ranks

    def __len__(self):
        return len(self.factorizations)

    def __getitem__(self, i):
        return self.factorizations[i]

    def insertion_tableau(self, i):
        return involution_insert(*self[i])[0]

    def recording_tableau(self, i):
        return involution_insert(*self[i])[1]

    @property
    def components(self):
        def dfs(seed):
            seen = {seed}
            added = {seed}
            while added:
                connected = set()
                for s in added:
                    connected |= \
                        {y for (x, y, e) in self.edges if s == x} | \
                        {x for (x, y, e) in self.edges if s == y}
                added = (connected - seen)
                seen |= added
            for x in seen:
                yield x

        rest = set(range(len(self)))
        components = []
        while rest:
            seed = rest.pop()
            comp = set(dfs(seed))
            rest -= comp
            components.append(comp)
        return components

    @property
    def edges(self):
        if self._edges is None:
            self._edges = []
            for x in range(len(self)):
                for i in range(self.num_factors):
                    y = self.f_i(x, i)
                    if y is not None:
                        self._edges.append((x, y, i))
        return self._edges

    def f_i(self, element_index, operator_index):
        if operator_index > 0:
            a = self[element_index][operator_index - 1].elements
            b = self[element_index][operator_index].elements
            while a and b:
                for j in range(len(a)):
                    if a[j] > b[-1]:
                        a = a[:j] + a[j + 1:]
                        break
                b = b[:-1]
            if len(a) == 0:
                return None
            x = max(a)
            a = set(self[element_index][operator_index - 1].elements) - {x}
            b = set(self[element_index][operator_index].elements)
            while x in b:
                x += 1
            b |= {x}
            ans = self[element_index]
            words = (Word(*sorted(a)), Word(*sorted(b)))
            ans = ans[:operator_index - 1] + words + ans[operator_index + 1:]
        if operator_index == 0:
            if self.num_factors < 2:
                return None
            a = self[element_index][0].elements
            b = self[element_index][1].elements
            if len(a) == 0:
                return None
            if len(b) > 0 and min(a) > min(b):
                return None
            x = min(a)
            a = set(a) - {x}
            b = set(b) | {x}
            ans = self[element_index]
            words = (Word(*sorted(a)), Word(*sorted(b)))
            ans = words + ans[2:]
        js = [j for j in range(len(self)) if self[j] == ans]
        assert len(js) == 1
        return js[0]

    def get_weight(self, i):
        return self.weights[i]

    def get_rank(self, j):
        wt = self.get_weight(j)
        n = len(wt)
        assert n == self.num_factors
        wt = [self.highest_weight[i] - wt[i] for i in range(len(wt))]
        ans = 0
        for i in range(n):
            ans += (n - i) * wt[i]
        return -ans

    @classmethod
    def get_increasing_factorizations(cls, w, k):
        def is_increasing(x):
            return all(x[i] > x[i - 1] for i in range(1, len(x)))
        #
        if k == 0 and len(w) == 0:
            yield tuple()
        elif k == 0:
            return
        elif len(w) == 0:
            yield tuple(k * [tuple()])
        else:
            for i in range(len(w) + 1):
                if is_increasing(w[:i]):
                    for tup in cls.get_increasing_factorizations(w[i:], k - 1):
                        yield (tuple(w[:i]),) + tup
                else:
                    break


class FPFCrystalGenerator(ShiftedCrystalGenerator):

    DIRECTORY = '/Users/emarberg/Dropbox/shifted_crystals/examples/fpf/'

    @classmethod
    def test_insertion_tableaux(cls, n, k):
        for i, w in enumerate(HopfPermutation.fpf_involutions(n)):
            cg = FPFCrystalGenerator(w.oneline, k)
            shapes = [
                {cg.insertion_tableau(i) for i in comp}
                for comp in cg.components
            ]
            b1 = all(len(sh) == 1 for sh in shapes)
            b2 = all(len(s & t) == 0 for s in shapes for t in shapes if s != t)
            print(i, w, b1, b2, '(n = %s, k = %s)' % (n, k))
            assert b1 and b2

    def insertion_tableau(self, i):
        return fpf_insert(*self[i])[0]

    def recording_tableau(self, i):
        return fpf_insert(*self[i])[1]

    def node_label(self, i):
        pre = str(self.insertion_tableau(i))
        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(j) for j in w.elements]) for w in self[i]])
        return pre + '\n\n' + top + '\n\n' + bottom

    def print_zero_weight_space(self):
        for f in self.zero_weight_space:
            p, q = fpf_insert(*f)
            print(p)
            print('')
            print(q)
            print('')
            print(f)
            print('\n')

    def _get_words(self):
        return get_fpf_involution_words(self.oneline)

    def __repr__(self):
        return 'FPF crystal generator for w = %s with ell = %s' % (self.oneline, self.num_factors)

    @classmethod
    def all(cls, n, k):
        assert n % 2 == 0
        for w in HopfPermutation.fpf_involutions(n):
            if w.oneline != reduce_fpf(w.oneline):
                continue
            cg = FPFCrystalGenerator(w.oneline, k)
            if cg.edges:
                cg.generate()

    def f_i(self, element_index, operator_index):
        if operator_index > 0:
            return super(FPFCrystalGenerator, self).f_i(element_index, operator_index)
        if self.num_factors < 2:
            return None
        a = self[element_index][0].elements
        b = self[element_index][1].elements
        if len(a) == 0:
            return None
        if len(b) > 0 and min(a) > min(b):
            return None
        if len(a) == 1:
            x = a[0]
            a = set(a) - {x}
            b = set(b) | {x}
        else:
            x, y = a[:2]
            if y % 2 == 0:
                a = set(a) - {x}
                b = set(b) | {x}
            elif y > x:
                a = set(a) - {y}
                b = set(b) | {y - 2}
        ans = self[element_index]
        words = (Word(*sorted(a)), Word(*sorted(b)))
        ans = words + ans[2:]
        js = [j for j in range(len(self)) if self[j] == ans]
        assert len(js) == 1
        return js[0]
