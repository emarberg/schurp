from vectors import Vector
from tableaux import Tableau
from marked import MarkedNumber
from collections import defaultdict
from signed import SignedPermutation, EvenSignedPermutation
import itertools
import numpy
from PIL import Image, ImageDraw, ImageColor
import random


class Word:

    def __init__(self, *args, **kwargs):
        self.subset = kwargs.get('subset', None) or set(args)
        self.elements = tuple(args)
        self._permutations = None
        self._fpf_involutions = None
        assert all(i in self.subset for i in args)

    def reverse(self):
        return Word(subset=self.subset, *tuple(reversed(self.elements)))

    def tuple(self):
        return self.elements

    @property
    def permutation_sequence(self):
        if self._permutations is None:
            n = max(self.elements) + 1
            words = [tuple(range(1, n + 1))]
            for counter, i in enumerate(self.elements):
                if counter % 1000 == 0:
                    print(counter, ':', len(self.elements))
                prev = list(words[-1])
                prev[i - 1:i + 1] = prev[i], prev[i - 1]
                words.append(tuple(prev))
            self._permutations = words
        return self._permutations

    @property
    def fpf_sequence(self):
        if self._fpf_involutions is None:
            n = max(self.elements) + 1
            invol = [tuple(i - 1 if i % 2 == 0 else i + 1 for i in range(1, n + 1))]
            for counter, i in enumerate(self.elements):
                if counter % 1000 == 0:
                    print(counter, ':', len(self.elements))
                prev = list(invol[-1])
                prev[i - 1:i + 1] = prev[i], prev[i - 1]
                prev = tuple(j if j not in {i, i + 1} else (i + 1 if j == i else i) for j in prev)
                invol.append(prev)
            self._fpf_involutions = [{i + 1: w[i] for i in range(n)} for w in invol]
        return self._fpf_involutions

    def draw_trajectories(self, k=10, filename='test.png'):
        n = max(self.elements) + 1
        t = random.sample(range(1, n + 1), k)
        colors = random.sample(list(ImageColor.colormap), k)
        pixels = {c: [c] for c in t}
        for counter, i in enumerate(self.elements):
            for c in t:
                if pixels[c][-1] == i:
                    pixels[c] += [i + 1]
                elif pixels[c][-1] == i + 1:
                    pixels[c] += [i]
                else:
                    pixels[c] += [pixels[c][-1]]
            if counter % 1000 == 0:
                print(counter, ':', len(self))

        xy = defaultdict(list)
        for c in t:
            for i, j in enumerate(pixels[c]):
                x, y = int(i * 800.0 / (len(self) + 1)), int((j - 1) * 500.0 / n)
                if xy[c] and (x, y) == xy[c][-1]:
                    continue
                xy[c] += [(x, y)]

        image = Image.new('RGBA', (800, 500))
        draw = ImageDraw.Draw(image)
        for c, color in zip(t, colors):
            for i in range(1, len(xy[c])):
                x1, y1 = xy[c][i - 1]
                x2, y2 = xy[c][i]
                draw.line((x1, y1, x2, y2), fill=color)
        image.save('images/trajectories/' + filename)

    def draw_fpf_trajectories(self, k=10, filename='test.png'):
        n = max(self.elements) + 1
        t = random.sample(range(1, n // 2 + 1), k)
        t = [(2 * i - 1, 2 * i) for i in t]
        colors = random.sample(list(ImageColor.colormap), k)
        pixels = {c: [c] for c in t}

        def toggle(x, i):
            if x == i:
                return i + 1
            elif x == i + 1:
                return i
            else:
                return x

        for counter, i in enumerate(self.elements):
            for c in t:
                x, y = pixels[c][-1]
                pixels[c] += [(toggle(x, i), toggle(y, i))]
            if counter % 1000 == 0:
                print(counter, ':', len(self))

        xy = defaultdict(list)
        for c in t:
            for i, pair in enumerate(pixels[c]):
                x = int(i * 800.0 / (len(self) + 1))
                y = int((pair[0] - 1) * 500.0 / n)
                z = int((pair[1] - 1) * 500.0 / n)
                if xy[c] and (x, y, z) == xy[c][-1]:
                    continue
                xy[c] += [(x, y, z)]

        image = Image.new('RGBA', (800, 500))
        draw = ImageDraw.Draw(image)
        for c, color in zip(t, colors):
            for i in range(1, len(xy[c])):
                x1, y1, z1 = xy[c][i - 1]
                x2, y2, z2 = xy[c][i]
                draw.line((x1, y1, x2, y2), fill=color)
                draw.line((x1, z1, x2, z2), fill=color)
        image.save('images/trajectories/' + filename)

    def print_permutation(self, m, filename='test.png'):
        pi = self.permutation_sequence[m]
        n = len(pi)
        h = 1000
        image = Image.new('RGBA', (h + 8, h + 8))
        draw = ImageDraw.Draw(image)

        for i, a in enumerate(pi):
            a -= 1
            x1, y1 = int(i * h / n), int(a * h / n)
            x2, y2 = x1 + 8, y1 + 8
            draw.ellipse((x1, y1, x2, y2), fill='black', outline='black')
        image.save('images/' + filename)

    def print_all(self):
        seq = self.permutation_sequence
        filename = 'test/test%09d.png'
        incr = max(1, len(seq) // 8)
        indices = list(range(0, len(seq), incr)) + [len(seq) - 1]
        for m, i in enumerate(indices):
            self.print_permutation(i, filename % m)
            print(m, '. . .')

    def print_fpf(self, m, filename='test.png'):
        invol = self.fpf_sequence[m]
        n = len(invol)
        h = 5 * n + 10
        image = Image.new('RGBA', (h, h))
        draw = ImageDraw.Draw(image)

        xcoord = {i + 1: h / 2.0 + (h - 10) / 2 * numpy.cos(numpy.pi * (0.5 + (2 * i + 1.0) / n)) for i in range(n)}
        ycoord = {i + 1: h / 2.0 - (h - 10) / 2 * numpy.sin(numpy.pi * (0.5 + (2 * i + 1.0) / n)) for i in range(n)}
        for i, j in invol.items():
            x1, y1 = xcoord[i], ycoord[i]
            x2, y2 = xcoord[j], ycoord[j]
            draw.line((x1, y1, x2, y2), fill='black')
        image.save('images/' + filename)

    def print_all_fpf(self):
        seq = self.fpf_sequence
        filename = 'fpf/test%09d.png'
        incr = max(1, len(seq) // 100)
        for i in list(range(0, len(seq), incr)) + [len(seq) - 1]:
            self.print_fpf(i, filename % i)
            print(i, '. . .')

    @classmethod
    def all(cls, n, l=None, packed=True):
        l = n if l is None else l
        for i in range(n + 1):
            for j in range(l**i):
                args = []
                for k in range(i):
                    args += [(j % l) + 1]
                    j = j // l
                if not packed or all(t == 1 or (t - 1) in args for t in args):
                    yield Word(*args)

    def wiring_diagram_tikz(self, n=None):
        n = max({0} | set(self.elements)) + 1 if n is None else n
        #
        pi = [tuple(range(1, n + 1))]
        for j in reversed(self.elements):
            tup = []
            for i in pi[-1]:
                if i == j:
                    tup += [i + 1]
                elif i == j + 1:
                    tup += [i - 1]
                else:
                    tup += [i]
            pi += [tuple(tup)]
        #
        s = []
        s += ['\\begin{center}']
        s += ['\\begin{tikzpicture}[scale=0.5]']
        for i in range(n):
            line = '\\draw (0,%s)' % i
            for j, x in enumerate(pi):
                if j == 0:
                    continue
                line += ' -- (%s,%s)' % (j, x[i] - 1)
            line += ';'
            s += [line]
        s += ['\\end{tikzpicture}']
        s += ['\\end{center}']
        return '\n'.join(s)

    @classmethod
    def increasing_zeta(cls, word):
        return 1 if all(word[i] < word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def decreasing_zeta(cls, word):
        return 1 if all(word[i] > word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_increasing_zeta(cls, word):
        return 1 if all(word[i] <= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def weakly_decreasing_zeta(cls, word):
        return 1 if all(word[i] >= word[i + 1] for i in range(len(word) - 1)) else 0

    @classmethod
    def unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def right_weakly_unimodal_zeta(cls, word):
        return sum([cls.decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def left_weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    @classmethod
    def weakly_unimodal_zeta(cls, word):
        return sum([cls.weakly_decreasing_zeta(word[:i]) * cls.weakly_increasing_zeta(word[i:]) for i in range(len(word) + 1)])

    def quasisymmetrize(self, zeta):
        if len(self) == 0:
            return Vector({(): 1})
        ans = defaultdict(int)
        for i in range(1, len(self) + 1):
            pre = Word(*self[:i])
            sub = Word(*self[i:])
            z = zeta(pre)
            if z != 0:
                for k, v in sub.quasisymmetrize(zeta).items():
                    ans[(i,) + k] += z * v
        ans = Vector({k: v for k, v in ans.items() if v != 0})
        return ans

    @classmethod
    def signed_involution_stable_grothendieck(cls, n, sigma, length_bound):
        ans = Vector()
        for w in sigma.get_involution_hecke_words(length_bound):
            ans += cls(*w).quasisymmetrize(cls.right_weakly_unimodal_zeta)

        def sort(t):
            return tuple(reversed(sorted(t)))

        assert all(ans.dictionary[sort(alpha)] == ans.dictionary[alpha] for alpha in ans.dictionary)
        ans = Vector({
            alpha: ans.dictionary[alpha]
            for alpha in ans.dictionary if sort(alpha) == alpha and len(alpha) <= n
        }, printer=lambda s: 'm[%s;%s]' % (','.join(map(str, s)), n) if s else '1')
        return ans

    @classmethod
    def permutations(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield Word(*args)

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
        return Word(subset=(self.subset | other.subset), *(self.elements + other.elements))

    def __lshift__(self, i):
        return Word(subset={x - i for x in self.subset}, *[e - i for e in self.elements])

    def __rshift__(self, i):
        return Word(subset={x + i for x in self.subset}, *[e + i for e in self.elements])

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

    def is_multishuffle(self):
        """Returns true if no adjancent letters are equal."""
        return not any(self[i] == self[i + 1] for i in range(len(self) - 1))

    def strict_coxeter_knuth_class(self):
        seen = set()
        add = {self}
        while add:
            newadd = set()
            for w in add:
                yield w
                seen.add(w)
                for i in range(len(w) - 2):
                    a, c, b = w[i:i + 3]
                    if a < b < c:
                        tup = w.elements[:i] + (c, a, b) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
                    c, a, b = w[i:i + 3]
                    if a < b < c:
                        tup = w.elements[:i] + (a, c, b) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
                    b, a, c = w[i:i + 3]
                    if a < b < c:
                        tup = w.elements[:i] + (b, c, a) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
                    b, c, a = w[i:i + 3]
                    if a < b < c:
                        tup = w.elements[:i] + (b, a, c) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
                    a, b, c = w[i:i + 3]
                    if c == a < b:
                        tup = w.elements[:i] + (b, a, b) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
                    b, a, c = w[i:i + 3]
                    if a < b == c:
                        tup = w.elements[:i] + (a, b, a) + w.elements[i + 3:]
                        newadd.add(Word(*tup, subset=w.subset))
            add = newadd - seen

    def modified_hecke_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.modified_hecke_insert(MarkedNumber(a), verbose=verbose)

            v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def rsk_insert(self):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.rsk_insert(MarkedNumber(a))
            v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

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

    def sagan_worley_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.sagan_worley_insert(MarkedNumber(a), verbose=verbose)
            v = MarkedNumber(-i if column_dir else i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def involution_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.involution_insert(MarkedNumber(a), verbose=verbose)
            v = MarkedNumber(-i if column_dir else i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def alt_involution_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.alt_involution_insert(MarkedNumber(a), verbose=verbose)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    if l == j:
                        q = q.set(k, l, MarkedNumber(-i))
                    if k == j:
                        q = q.set(k, l, MarkedNumber(i))
            assert p.shape() == q.shape()
        return p, q

    def clan_insert(self, n, verbose=True):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.clan_insert(n, MarkedNumber(a), verbose=verbose)

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

    def mystery_insert(self, verbose=True):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, p = p.mystery_insert(MarkedNumber(a), verbose=verbose)

            v = MarkedNumber(i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def fpf_insert(self, verbose=False):
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

    def alt_fpf_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            j, column_dir, p = p.alt_fpf_insert(MarkedNumber(a), verbose=verbose)

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
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
    for i in range(1, 1 + p.max_row):
        if p.entry(i, i).is_primed():
            p = p.set(i, i, -p.entry(i, i))
            q = q.set(i, i, -q.entry(i, i))
    return p, q


def fpf_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.fpf_insert(verbose=False)
    return p, Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})


def alt_fpf_insert(*words):
    w, mapping = get_insertion_mapping(words)
    mapping[MarkedNumber(0)] = MarkedNumber(0)
    p, q = w.alt_fpf_insert(verbose=False)
    p = p.halve()
    q = q.halve()
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
    oneline = tuple(oneline)
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

