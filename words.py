from vectors import Vector
from tableaux import Tableau
from marked import MarkedNumber
from collections import defaultdict
import itertools
import numpy
import random


class Word:

    def __init__(self, *args, **kwargs):
        args = args[0] if len(args) == 1 and type(args[0]) != int else args
        self.subset = kwargs.get('subset', None) or set(args)
        self.elements = tuple(args)
        self._permutations = None
        self._fpf_involutions = None
        assert all(i in self.subset for i in args)

    @classmethod
    def run_decomposition(cls, w):
        ans = []
        for i in range(len(w)):
            if i == 0 or w[i - 1] > w[i]:
                ans.append([])
            ans[-1].append(w[i])
        return ans

    @classmethod
    def drop(cls, top, bot=None):
        if bot is None:
            rho = cls.run_decomposition(top)
            should_continue = True
            while should_continue:
                should_continue = False
                for k in range(len(rho) - 1):
                    top = rho[k]
                    bot = rho[k + 1]
                    if len(top) > len(bot) or any(top[i] <= bot[i] for i in range(len(top))):
                        should_continue = True
                        rho[k], rho[k + 1] = cls.drop(top, bot)
                        break
            rho = list(reversed(rho))
            dictionary = {(i + 1, j + 1): rho[i][j] for i in range(len(rho)) for j in range(len(rho[i]))}
            ans = Tableau(dictionary)
            # assert ans.is_increasing()
            return ans
        if bot is not None:
            top, bot = cls.drop_alignment(top, bot)

            n = len(top)
            assert n == len(bot)

            for i in range(n):
                if bot[i] is None:
                    x = top[i]
                    bot[i] = x
                    top[i] = None

                    if i + 1 < n and x + 1 == top[i + 1]:
                        # error in Definition 4.1.3 in Assaf 1903.05802v1 (not corrected in published version):
                        # should be b_j = max{b : x_j + i == tau_i^(j) == sigma_i^(j) + 1 for 1<=i<=b }
                        j = 1
                        while i + j < n and bot[i + j] is not None and bot[i + j] + 1 == top[i + j] == x + j:
                            bot[i + j] += 1
                            j += 1
            return tuple(a for a in top if a is not None), tuple(a for a in bot if a is not None)

    @classmethod
    def drop_alignment(cls, top, bot):
        ans = [[], []]
        i = j = 0
        while i < len(top) or j < len(bot):
            if i == len(top):
                ans[0].append(None)
                ans[1].append(bot[j])
                j += 1
            elif j == len(bot):
                ans[0].append(top[i])
                ans[1].append(None)
                i += 1
            elif top[i] > bot[j]:
                ans[0].append(top[i])
                ans[1].append(bot[j])
                i += 1
                j += 1
            else:
                ans[0].append(top[i])
                ans[1].append(None)
                i += 1
        return ans

    @classmethod
    def lift_sequence(cls, rho, i, j=None):
        if rho is None:
            return None
        if j is not None:
            assert i <= j
            for t in range(i, j + 1):
                rho = cls.lift_sequence(rho, t)
            return rho
        if j is None:
            if i > len(rho):
                return None
            ans = rho[:]
            if i == len(ans):
                ans.append([])
            top = ans[i]
            bot = ans[i - 1]
            ntop, nbot = cls.lift(top, bot)
            if not (top and ntop and bot and ntop[0] == bot[0]):
                ans[i], ans[i - 1] = ntop, nbot
            return ans if ans != rho else None

    @classmethod
    def lift(cls, top, bot=None):
        # error after Definition 4.2.6 in Assaf 1903.05802v1 (corrected in published version):
        # lifting first example in Fig. 12 does not give every in Fig. 10

        if bot is None:
            rho = list(reversed(cls.run_decomposition(top)))
            for k in range(len(rho) - 1, -1, -1):
                j = rho[k][0] - 1
                indices = [i for i in range(1, j + 1) if cls.lift_sequence(rho, i, j) is not None]
                if indices:
                    i = indices[0]
                    rho = cls.lift_sequence(rho, i, j)
            dictionary = {(i + 1, j + 1): rho[i][j] for i in range(len(rho)) for j in range(len(rho[i]))}
            ans = Tableau(dictionary)
            # assert ans.is_increasing()
            return ans

        # errors/omissions in Definitions 4.19 and 4.22 of Assaf 1903.05802v1 published version:
        #
        # * lift_i(T) and lift(rho) are not ever formally defined?
        # * lift_i(T) should move row i to row i + 1 if row i + 1 is empty
        # * if P has n - 1 rows then j_k in Def. 4.22 should be (first entry in row n - k of P) - 1
        #
        if bot is not None:
            top, bot = cls.lift_alignment(top, bot)

            n = len(top)
            assert n == len(bot)

            for i in range(n - 1, -1, -1):
                if top[i] is None:
                    x = bot[i]
                    top[i] = x
                    bot[i] = None

                    if i + 1 < n and x + 1 == bot[i + 1]:
                        j = 1
                        while i + j < n and top[i + j] is not None and top[i + j] == bot[i + j] == x + j:
                            bot[i + j] -= 1
                            j += 1
            return tuple(a for a in top if a is not None), tuple(a for a in bot if a is not None)

    @classmethod
    def lift_alignment(cls, top, bot):
        ans = [[], []]
        i = len(top) - 1
        j = len(bot) - 1
        while i >= 0 or j >= 0:
            if i < 0:
                ans[0].append(None)
                ans[1].append(bot[j])
                j -= 1
            elif j < 0:
                ans[0].append(top[i])
                ans[1].append(None)
                i -= 1
            elif top[i] >= bot[j]:
                ans[0].append(top[i])
                ans[1].append(bot[j])
                i -= 1
                j -= 1
            else:
                ans[0].append(None)
                ans[1].append(bot[j])
                j -= 1
        top, bot = ans
        top = list(reversed(top))
        bot = list(reversed(bot))
        return [top, bot]

    @classmethod
    def _incr_pairing(cls, x, y):
        pairing = [[], []]
        for b in reversed(y):
            upper = [a for a in x if abs(a) > abs(b) and a not in pairing[0]]
            if upper:
                a = upper[0]
                pairing[0] = [a] + pairing[0]
                pairing[1] = [b] + pairing[1]
        return pairing

    @classmethod
    def incr_crystal_e(cls, tup, index):
        if len(tup) == 0:
            return None

        if index == 0:
            if len(tup[0]) == 0 or tup[0][0] > 0:
                return None
            return ((-tup[0][0],) + tuple(tup[0][1:]),) + tuple(tup[1:])

        elif index == -1:
            if len(tup[1]) == 0 or (len(tup[0]) > 0 and abs(tup[0][0]) <= abs(tup[1][0])):
                return None
            if len(tup[0]) == 0 or tup[0][0] * tup[1][0] > 0:
                x, y = (tup[1][0],) + tuple(tup[0]), tuple(tup[1][1:])
            else:
                x, y = (-tup[1][0], -tup[0][0]) + tuple(tup[0][1:]), tuple(tup[1][1:])
            return (x, y) + tuple(tup[2:])

        elif index > 0:
            ans = [list(tup[i]) for i in [abs(index) - 1, abs(index)]]
            pairing = cls._incr_pairing(ans[0], ans[1])
            left = [b for b in ans[1] if b not in pairing[1]]
            if left:
                b = left[0]
                ans[1] = [x for x in ans[1] if x != b]

                a = b
                while a in ans[0] or -a in ans[0]:
                    a = (a - 1) if a > 0 else (a + 1)
                ans[0] = sorted(ans[0] + [a], key=abs)

                if a > 0:
                    for t in range(a, b):
                        if -t - 1 in ans[0]:
                            ans[0] = sorted([x for x in ans[0] if x != -t - 1] + [t + 1], key=abs)
                            ans[1] = sorted([x for x in ans[1] if x != t] + [-t], key=abs)

            else:
                return None
            return tuple(tup[:index - 1]) + tuple(tuple(a) for a in ans) + tuple(tup[index + 1:])

        elif index < -1:
            pass

        else:
            raise Exception

    @classmethod
    def incr_crystal_f(cls, tup, index):
        if len(tup) == 0:
            return None

        if index == 'FPF':
            if len(tup[0]) == 0 or (len(tup[1]) > 0 and tup[1][0] <= tup[0][0]):
                return None
            if len(tup[0]) <= 1 or tup[0][0] + 1 not in tup[0]:
                x, y = tuple(tup[0][1:]), (tup[0][0],) + tuple(tup[1])
            else:
                x, y = (tup[0][0],) + tuple(tup[0][2:]), (tup[0][0] - 1,) + tuple(tup[1])
            return (x, y) + tuple(tup[2:])

        if index == 0:
            if len(tup[0]) == 0 or tup[0][0] < 0:
                return None
            return ((-tup[0][0],) + tuple(tup[0][1:]),) + tuple(tup[1:])

        elif index == -1:
            if len(tup[0]) == 0 or (len(tup[1]) > 0 and abs(tup[1][0]) <= abs(tup[0][0])):
                return None
            if len(tup[0]) <= 1 or tup[0][0] * tup[0][1] > 0:
                x, y = tuple(tup[0][1:]), (tup[0][0],) + tuple(tup[1])
            else:
                x, y = (-tup[0][1],) + tuple(tup[0][2:]), (-tup[0][0],) + tuple(tup[1])
            return (x, y) + tuple(tup[2:])

        elif index > 0:
            ans = [list(tup[i]) for i in [index - 1, index]]
            pairing = cls._incr_pairing(ans[0], ans[1])
            left = [a for a in ans[0] if a not in pairing[0]]
            if left:
                a = left[-1]
                ans[0] = [x for x in ans[0] if x != a]

                b = a
                while b in ans[1] or -b in ans[1]:
                    b = (b + 1) if b > 0 else (b - 1)
                ans[1] = sorted(ans[1] + [b], key=abs)

                if a > 0:
                    for t in range(a, b):
                        if -t in ans[1]:
                            ans[0] = sorted([x for x in ans[0] if x != t + 1] + [-t - 1], key=abs)
                            ans[1] = sorted([x for x in ans[1] if x != -t] + [t], key=abs)
            else:
                return None
            return tuple(tup[:index - 1]) + tuple(tuple(a) for a in ans) + tuple(tup[index + 1:])

        elif index < -1:
            pass

        else:
            raise Exception

    @classmethod
    def increasing_factorizations(cls, w, k):
        if k < 0:
            return

        w = w.tuple() if type(w) == Word else w

        def is_increasing(x):
            return all(abs(x[i]) > abs(x[i - 1]) for i in range(1, len(x)))

        if k == 0 and len(w) == 0:
            yield tuple()
            return
        elif k == 0:
            return
        elif len(w) == 0:
            yield tuple(Word() for _ in range(k))
            return

        for i in range(len(w) + 1):
            if not is_increasing(w[:i]):
                break
            for tup in cls.increasing_factorizations(w[i:], k - 1):
                yield (Word(*w[:i]),) + tup

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
        from PIL import Image, ImageDraw, ImageColor
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
        from PIL import Image, ImageDraw, ImageColor
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
        from PIL import Image, ImageDraw, ImageColor
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
        from PIL import Image, ImageDraw, ImageColor
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

    @classmethod
    def coxeter_knuth_move(cls, word, i):
        if 0 <= i < len(word) - 2:
            a, c, b = word[i:i + 3]
            if a < b < c or c < b < a:
                return word[:i] + (c, a, b) + word[i + 3:]
            b, a, c = word[i:i + 3]
            if a < b < c or c < b < a:
                return word[:i] + (b, c, a) + word[i + 3:]
            if b == c:
                return word[:i] + (a, b, a) + word[i + 3:]
        return word

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

    def sp_mixed_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            a = MarkedNumber(2 * a)
            j, column_dir, p = p.sp_mixed_insert(a, verbose=verbose)

            for k, l in p.shape():
                if (k, l) not in q.shape():
                    v = MarkedNumber(i * (-1 if p.get(k, l).number < 0 else 1))
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        mapping = {}
        for i, j in p:
            v = p.mapping[(i, j)]
            if v.number % 2 == 0:
                v = MarkedNumber(abs(v) // 2)
            else:
                v = MarkedNumber(-(abs(v) + 1) // 2)
            mapping[(i, j)] = v
        return Tableau(mapping), q

    def mixed_insert(self, verbose=False):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = i_zerobased + 1
            a = MarkedNumber(2 * a)
            j, column_dir, p = p.mixed_insert(a, verbose=verbose)

            for k, l in p.shape():
                if (k, l) not in q.shape():
                    v = MarkedNumber(i * (-1 if p.get(k, l).number < 0 and k != l else 1))
                    q = q.set(k, l, v)
            assert p.shape() == q.shape()
        mapping = {}
        for i, j in p:
            v = p.mapping[(i, j)]
            if i == j:
                assert v.number % 2 == 0
                v = MarkedNumber(v.number // 2)
            elif v.number % 2 == 0:
                v = MarkedNumber(abs(v) // 2)
            else:
                v = MarkedNumber(-(abs(v) + 1) // 2)
            mapping[(i, j)] = v
        return Tableau(mapping), q

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

    def primed_sw_insert(self, verbose=False, phi=None):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = (i_zerobased + 1) if phi is None else phi[i_zerobased]
            j, column_dir, p = p.primed_sw_insert(MarkedNumber(a), verbose=verbose)
            v = MarkedNumber(-i if column_dir else i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    if k == l and p.get(k, l).is_primed():
                        p = p.set(k, l, -p.get(k, l))
                        q = q.set(k, l, -v)
                    else:
                        q = q.set(k, l, v)
            assert p.shape() == q.shape()
        return p, q

    def involution_insert(self, verbose=False, phi=None):
        p, q = Tableau(), Tableau()
        for i_zerobased, a in enumerate(self):
            i = (i_zerobased + 1) if phi is None else phi[i_zerobased]
            j, column_dir, p = p.involution_insert(MarkedNumber(a), verbose=verbose)
            v = MarkedNumber(-i if column_dir else i)
            for k, l in p.shape():
                if (k, l) not in q.shape():
                    if k == l and p.get(k, l).is_primed():
                        p = p.set(k, l, -p.get(k, l))
                        q = q.set(k, l, -v)
                    else:
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
    words = [Word(*a) if type(a) != Word else a for a in words]
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


def primed_sw_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.primed_sw_insert(verbose=False)
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
    return p, q


def sagan_worley_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.sagan_worley_insert(verbose=False)
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
    return p, q


def mixed_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.mixed_insert(verbose=False)
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
    return p, q


def sp_mixed_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.sp_mixed_insert(verbose=False)
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
    return p, q


def involution_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.involution_insert(verbose=False)
    q = Tableau({(i, j): mapping[q.entry(i, j)] for (i, j) in q})
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


def eg_insert(*words):
    w, mapping = get_insertion_mapping(words)
    p, q = w.eg_insert()
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
