import itertools
from .symmetric import SymmetricPolynomial, SymmetricMonomial
from .words import Word
from .vectors import Vector
from .tableaux import Partition, Tableau
from .polynomials import beta as BETA # noqa
from collections import defaultdict

HECKE_CACHE = defaultdict(list)
SYMPLECTIC_HECKE_CACHE = defaultdict(list)
ORTHOGONAL_HECKE_CACHE = defaultdict(list)

FPF_INVOLUTION_WORDS = {(): {()}}
INVOLUTION_WORDS = {(): {()}}

SIGNED_REDUCED_WORDS = {(): {()}}

EVEN_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_COUNTS = {(): 1}
EVEN_SIGNED_INVOLUTION_WORDS = {(): [()]}

ATOMS_FPF_CACHE = {}
atoms_b_cache = {}
atoms_d_cache = {}


class Permutation:

    """
    Simple class for (signed) permutatations
    """

    def __init__(self, *oneline):
        if len(oneline) == 1 and type(oneline[0]) != int:
            oneline = tuple(oneline[0])
        while oneline and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        assert set(range(1, len(oneline) + 1)) == set(abs(i) for i in oneline)

        self._oneline = tuple(oneline)
        self._rank = len(oneline)
        # cached fields
        self._rdes = None
        self._ldes = None
        self._len = None

    def star(self):
        n = self.rank
        return Permutation(*[n + 1 - i for i in reversed(self.oneline)])

    @property
    def oneline(self):
        return self._oneline

    @property
    def rank(self):
        return self._rank

    def __repr__(self):
        if self.is_unsigned():
            return self.cycle_repr()
        else:
            return self.oneline_repr()

    def oneline_repr(self):
        return (',' if self.rank > 9 else '').join([str(i) for i in self.oneline])

    @property
    def cycles(self):
        oneline = self.oneline
        cycles = []
        numbers = list(range(1, len(oneline) + 1))
        while len(numbers) > 0:
            cycle = [numbers[0]]
            nxt = oneline[numbers[0] - 1]
            numbers.remove(numbers[0])
            while nxt != cycle[0]:
                cycle.append(nxt)
                numbers.remove(nxt)
                nxt = oneline[cycle[len(cycle) - 1] - 1]
            cycles.append(cycle)
        return cycles

    def cycle_repr(self):
        if len(self) == 0:
            return '1'

        s = ''
        for c in self.cycles:
            if len(c) > 1:
                s += '(' + ', '.join([str(x) for x in c]) + ')'
        return s

    @classmethod
    def from_word(cls, word):
        ans = Permutation()
        for i in word:
            ans *= Permutation.s_i(i)
        return ans

    @classmethod
    def from_involution_word(cls, word, strict=True):
        w = cls.identity()
        for i in word:
            s = Permutation.s_i(i)
            if i in w.right_descent_set and strict:
                return None
            elif s * w == w * s:
                w = w * s
            else:
                w = s * w * s
        return w

    @classmethod
    def from_fpf_involution_word(cls, word, strict=True):
        n = (1 + max(word)) if word else 0
        n = n if n % 2 == 0 else n + 1
        w = cls.identity()
        for i in range(1, n, 2):
            w *= cls.s_i(i)
        for i in word:
            s = Permutation.s_i(i)
            if i in w.right_descent_set and strict:
                return None
            elif s * w == w * s and strict:
                return None
            else:
                w = s * w * s
        return w

    @classmethod
    def all(cls, n, signed=False):
        for args in itertools.permutations(range(1, n + 1)):
            if not signed:
                yield Permutation(*args)
            else:
                for v in range(2**n):
                    oneline = []
                    for i in range(n):
                        oneline.append(args[i] * (-1) ** (v % 2))
                        v = v // 2
                    yield Permutation(*oneline)

    @classmethod
    def involutions(cls, n, signed=False, fpf=False):
        if n % 2 != 0 and fpf:
            return

        def create(delta, fixed, numbers):
            w = Permutation()
            for i in range(len(delta)):
                w *= cls.transposition(numbers[0], numbers[delta[i]])
                numbers = numbers[1:delta[i]] + numbers[delta[i] + 1:]
            return w

        if not signed:
            for k in [0] if fpf else range(n, -1, -2):
                for fixed in itertools.combinations(range(1, n + 1), k):
                    numbers = [i for i in range(1, n + 1) if i not in fixed]
                    seen = set()
                    queue = [tuple(range(n - k - 1, 0, -2))]
                    while queue:
                        delta, queue = queue[0], queue[1:]
                        yield create(delta, fixed, numbers)

                        for i in range(len(delta)):
                            if delta[i] > 1:
                                newdelta = delta[:i] + (delta[i] - 1,) + delta[i + 1:]
                                if newdelta not in seen:
                                    queue.append(newdelta)
                                    seen.add(newdelta)
        else:
            for w in cls.involutions(n, False):
                a = [i for i in range(1, n + 1) if i <= w(i)]
                for k in range(len(a) + 1):
                    for conjugators in itertools.combinations(a, k):
                        v = w
                        for i in conjugators:
                            s = cls.reflection_s(i, i)
                            v = v * s if v(i) == i else s * v * s
                        yield v

    @classmethod
    def fpf_involutions(cls, n):
        if n % 2 != 0:
            return
        crop = {Permutation.minfpf(n)}
        while crop:
            for z in crop:
                yield z
            crop = {
                Permutation.s_i(i) * z * Permutation.s_i(i)
                for z in crop
                for i in range(1, n)
                if z(i) < z(i + 1)
            }

    def rothe_diagram(self):
        n = len(self.oneline)
        return sorted([
            (i, self(j)) for i in range(1, n + 1) for j in range(i + 1, n + 1) if self(i) > self(j)
        ])

    def fpf_rothe_diagram(self, fpf=False):
        return self.involution_rothe_diagram(True)

    def involution_rothe_diagram(self, fpf=False):
        return [(i, j) for (i, j) in self.rothe_diagram() if i > j or (not fpf and i == j)]

    def print_rothe_diagram(self, french=False, sep=' '):
        print(self.print_diagram(self.rothe_diagram(), french=french, sep=sep))

    def print_fpf_rothe_diagram(self, french=False, sep=' '):
        print(self.print_diagram(self.involution_rothe_diagram(True), french=french, sep=sep))

    def print_involution_rothe_diagram(self, french=False, sep=' '):
        print(self.print_diagram(self.involution_rothe_diagram(False), french=french, sep=sep))

    @classmethod
    def print_diagram(cls, diagram, french=False, sep=' '):
        if not diagram:
            return ''
        rows = max([a[0] for a in diagram])
        cols = max([a[1] for a in diagram])
        arr = [[sep for i in range(cols)] for j in range(rows)]
        for a in diagram:
            i, j = tuple(a[:2])
            arr[i - 1][j - 1] = '*'
        tojoin = [''.join(row) for row in arr]
        if french:
            return '\n'.join(reversed(tojoin))
        else:
            return '\n'.join([''.join(row) for row in arr])

    def code(self):
        assert self.is_unsigned()
        code = ()
        for i in range(1, 1 + self.rank):
            code += (len([j for j in range(i + 1, self.rank + i) if self(j) < self(i)]),)
        return code

    def involution_code(self):
        assert self.is_unsigned() and self.is_involution()
        code = ()
        for i in range(1, 1 + self.rank):
            code += (len([
                j for j in range(i + 1, self.rank + i)
                if self(j) <= i and self(j) < self(i)
            ]),)
        return code

    def fpf_involution_code(self):
        assert self.is_unsigned() and self.is_fpf_involution()
        code = ()
        for i in range(1, 1 + self.rank):
            code += (len([
                j for j in range(i + 1, self.rank + i)
                if self(j) < i and self(j) < self(i)
            ]),)
        return code

    def shape(self):
        return Partition.sort(self.inverse().code(), trim=True)

    def involution_shape(self):
        assert self.is_involution()
        mu = Partition.sort(self.involution_code(), trim=True)
        return Partition.transpose(mu)

    @classmethod
    def grassmannians(cls, rank):
        delta = tuple(range(rank - 1, 0, -1))
        for mu in Partition.subpartitions(delta):
            yield cls.get_grassmannian(*mu)

    @classmethod
    def inv_grassmannians(cls, rank):
        delta = tuple(range(rank - 1, 0, -2))
        for mu in Partition.subpartitions(delta, strict=True):
            yield cls.get_inv_grassmannian(*mu)

    @classmethod
    def fpf_grassmannians(cls, rank):
        assert rank % 2 == 0
        delta = tuple(range(rank - 2, 0, -2))
        for mu in Partition.subpartitions(delta, strict=True):
            yield cls.get_fpf_grassmannian(*mu)

    @classmethod
    def get_grassmannian(cls, *mu):
        oneline = tuple(i + 1 + a for i, a in enumerate(sorted(mu)))
        if oneline:
            missing = set(range(1, oneline[-1] + 1)) - set(oneline)
            oneline += tuple(sorted(missing))
        return Permutation(*oneline).inverse()

    @classmethod
    def get_inv_grassmannian(cls, *mu):
        assert Partition.is_strict_partition(mu)
        ans = Permutation()
        for i in range(len(mu)):
            ans *= Permutation.transposition(1 + mu[0] - mu[i], i + 1 + mu[0])
        w0 = Permutation.longest_element(len(mu) + (mu[0] if mu else 0))
        ans = w0 * ans * w0
        return ans

    @classmethod
    def get_fpf_grassmannian(cls, *mu):
        assert Partition.is_strict_partition(mu)
        ans = Permutation()
        o = 1 if (mu and (len(mu) + 1 + mu[0]) % 2 != 0) else 0
        for i in range(len(mu)):
            ans *= Permutation.transposition(o + 1 + mu[0] - mu[i], o + i + 2 + mu[0])
        while not ans.is_fpf_involution():
            f = [i for i in range(1, ans.rank + 2) if ans(i) == i]
            ans *= Permutation.transposition(f[0], f[1])
        w0 = Permutation.longest_element(o + 1 + len(mu) + (mu[0] if mu else 0))
        ans = w0 * ans * w0
        return ans

    def fpf_involution_shape(self):
        assert self.is_fpf_involution()
        mu = Partition.sort(self.fpf_involution_code(), trim=True)
        return Partition.transpose(mu)

    def _fpf_grassmannian_shape(self):
        cycles = [(a, self(a)) for a in range(1, self.rank + 1) if a < self(a)]
        filtered = sorted([(a, b) for a, b in cycles if any(i < b and a < j < b for i, j in cycles)])
        r = len(filtered)
        if r == 0:
            return ()
        phi = [a for a, b in filtered]
        top = [b for a, b in filtered]
        n = top[0] - 1
        if phi[-1] >= n + 1 or top != [n + i + 1 for i in range(r)]:
            return None
        return tuple(n - a for a in phi)

    def is_grassmannian(self):
        return len(self.right_descent_set) <= 1

    def is_vexillary(self):
        n = self.rank
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    for l in range(k + 1, n + 1):
                        b, a, d, c = self(i), self(j), self(k), self(l)
                        if a < b < c < d:
                            return False
        return True

    def is_quasi_dominant(self):
        if not self.is_vexillary():
            return False
        n = len(self.oneline)
        for a in range(1, n + 1):
            if a < self(a) and not all(b < self(b) for b in range(1, a)):
                return False
        return True

    def is_dominant(self):
        n = len(self.oneline)
        for i in range(1, n - 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    if self(i) < self(k) < self(j):
                        return False
        return True

    def is_inv_grassmannian(self):
        return self.get_inv_grassmannian(*self.involution_shape()) == self

    def is_fpf_grassmannian(self):
        if not self.is_unsigned() and self.is_fpf_involution():
            return False
        return self._fpf_grassmannian_shape() is not None

    def is_fpf_involution(self):
        return self == self.inverse() and all(self(i) != i for i in range(1, self.rank + 1))

    def is_involution(self):
        return self == self.inverse()

    def is_identity(self):
        return self == self.identity()

    def is_unsigned(self):
        return all(i > 0 for i in self.oneline)

    @classmethod
    def fpf_involution_little_push(cls, word, i):
        new = word[:i] + (word[i] + 1,) + word[i + 1:]
        w = cls.from_fpf_involution_word(new, strict=False)
        if w.fpf_involution_length == len(new):
            return new, None
        elif w.fpf_involution_length == len(new) - 1:
            # print('  ', new, i, Permutation.from_word(new).length == len(new))
            return new, i
        v = cls.from_fpf_involution_word(word[:i] + word[i + 1:], strict=False)
        for j in range(len(new)):
            if i != j and v == cls.from_fpf_involution_word(new[:j] + new[j + 1:], strict=False):
                # print('  ', new, j, Permutation.from_word(new).length == len(new))
                return new, j
        raise Exception

    @classmethod
    def fpf_involution_little_bump(cls, word, *args):
        assert len(args) in [1, 2]

        w = cls.from_fpf_involution_word(word, strict=False)
        assert w.fpf_involution_length == len(word)

        if len(args) == 1:
            y = args[0]
            assert type(y) == Permutation
            assert y.is_fpf_involution()
            for i in range(len(word)):
                subword = word[:i] + word[i + 1:]
                if subword in y.get_fpf_involution_words():
                    while i is not None:
                        word, i = cls.fpf_involution_little_push(word, i)
                    return word
            return word
        else:
            j, k = args[0], args[1]
            assert type(j) == type(k) == int
            t = cls.transposition(j, k)

            n = (max(word) + 1) if word else 0
            n = n + 1 if n % 2 != 0 else n
            minfpf = Permutation.minfpf(n)

            subatoms = [x * t for x in w.get_fpf_atoms() if len(x * t) == len(x) - 1]
            subatoms = [x for x in subatoms if (x.inverse() * minfpf * x).fpf_involution_length == len(x)]
            if len(subatoms) == 0:
                return word
            y = subatoms[0].inverse() * minfpf * subatoms[0]
            return cls.fpf_involution_little_bump(word, y)

    @classmethod
    def minfpf(cls, n):
        assert n % 2 == 0
        return Permutation(*[i + (1 if i % 2 != 0 else -1) for i in range(1, n + 1)])

    @classmethod
    def involution_little_push(cls, word, i):
        new = word[:i] + (word[i] + 1,) + word[i + 1:]
        v = cls.from_involution_word(word[:i] + word[i + 1:], strict=False)
        w = cls.from_involution_word(new, strict=False)
        if w.involution_length == len(new):
            return new, None
        for j in range(len(new)):
            if i != j and v == cls.from_involution_word(new[:j] + new[j + 1:], strict=False):
                return new, j
        raise Exception

    @classmethod
    def involution_little_bump(cls, word, *args):
        assert len(args) in [1, 2]

        w = cls.from_involution_word(word, strict=False)
        assert w.involution_length == len(word)

        if len(args) == 1:
            y = args[0]
            assert type(y) == Permutation
            for i in range(len(word)):
                subword = word[:i] + word[i + 1:]
                if Permutation.from_involution_word(subword, strict=False) == y:
                    while i is not None:
                        word, i = cls.involution_little_push(word, i)
                    return word
            return word
        else:
            j, k = args[0], args[1]
            assert type(j) == type(k) == int
            t = cls.transposition(j, k)
            subatoms = [x * t for x in w.get_atoms() if len(x * t) == len(x) - 1]
            subatoms = [x for x in subatoms if (x.inverse() % x).involution_length == len(x)]
            if len(subatoms) == 0:
                return word
            y = subatoms[0].inverse() % subatoms[0]
            return cls.involution_little_bump(word, y)

    @classmethod
    def little_push(cls, word, i):
        new = word[:i] + (word[i] + 1,) + word[i + 1:]
        v = cls.from_word(word[:i] + word[i + 1:])
        w = cls.from_word(new)
        if len(w) == len(new):
            return new, None
        for j in range(len(new)):
            if i != j and v == cls.from_word(new[:j] + new[j + 1:]):
                return new, j
        raise Exception

    @classmethod
    def little_bump(cls, word, *args):
        assert len(args) in [1, 2]

        w = cls.from_word(word)
        assert w.length == len(word)

        if len(args) == 1:
            y = args[0]
            assert type(y) == Permutation
            for i in range(len(word)):
                subword = word[:i] + word[i + 1:]
                if Permutation.from_word(subword) == y:
                    while i is not None:
                        word, i = cls.little_push(word, i)
                    return word
            return word
        else:
            j, k = args[0], args[1]
            assert type(j) == type(k) == int
            t = cls.transposition(j, k)
            if len(w * t) != len(w) - 1:
                return word
            return cls.little_bump(word, w * t)

    @classmethod
    def hecke_words(cls, n, length_bound=None, ascent_bound=None):
        for level in cls.hecke_levels(n, length_bound, ascent_bound):
            for pi, w in level:
                yield w

    @classmethod
    def hecke_levels(cls, n, length_bound=None, ascent_bound=None):
        length_bound = -1 if length_bound is None else length_bound
        key = (n, ascent_bound)

        def generate():
            for level in HECKE_CACHE[key][:-1]:
                yield level
            if HECKE_CACHE[key]:
                level = HECKE_CACHE[key][-1]
            else:
                start = (cls.identity(), ())
                level = {start}
            while level:
                next_level = set()
                yield level
                for pi, w in level:
                    for i in range(1, n):
                        s = Permutation.s_i(i)
                        sigma = pi % s
                        v = w + (i,)
                        if ascent_bound is None or Word.ascents(v) < ascent_bound:
                            next_level.add((sigma, v))
                level = next_level

        for i, level in enumerate(generate()):
            if i >= len(HECKE_CACHE[key]):
                HECKE_CACHE[key].append(level)
            yield level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_hecke_words(self, length_bound=None, ascent_bound=None):
        for level in self.hecke_levels(self.rank, length_bound, ascent_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    @classmethod
    def symplectic_hecke_words(cls, n, length_bound=None, ascent_bound=None, hecke=True):
        assert n % 2 == 0
        for level in cls.symplectic_hecke_levels(n, length_bound, ascent_bound, hecke):
            for pi, w in level:
                yield w

    @classmethod
    def symplectic_hecke_levels(cls, n, length_bound=None, ascent_bound=None, hecke=True):
        assert n % 2 == 0
        length_bound = -1 if length_bound is None else length_bound
        key = (n, ascent_bound, hecke)

        def generate():
            for level in SYMPLECTIC_HECKE_CACHE[key][:-1]:
                yield level
            if SYMPLECTIC_HECKE_CACHE[key]:
                level = SYMPLECTIC_HECKE_CACHE[key][-1]
            else:
                a = cls.identity()
                for i in range(1, n, 2):
                    a *= cls.s_i(i)
                start = (a, ())
                level = {start}
            while level:
                next_level = set()
                yield level
                for pi, w in level:
                    for i in range(1, n):
                        s = Permutation.s_i(i)
                        if pi(i) != i + 1:
                            sigma = s % pi % s
                            v = w + (i,)
                            if ascent_bound is None or Word.ascents(v) < ascent_bound:
                                if hecke or sigma != pi:
                                    next_level.add((sigma, v))
                level = next_level

        for i, level in enumerate(generate()):
            if i >= len(SYMPLECTIC_HECKE_CACHE[key]):
                SYMPLECTIC_HECKE_CACHE[key].append(level)
            yield level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_symplectic_hecke_words(self, length_bound=None, ascent_bound=None):
        for level in self.symplectic_hecke_levels(self.rank, length_bound, ascent_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    @classmethod
    def involution_hecke_words(cls, n, length_bound=None, peak_bound=None):
        for level in cls.involution_hecke_levels(n, length_bound, peak_bound):
            for pi, w in level:
                yield w

    @classmethod
    def involution_hecke_levels(cls, n, length_bound=None, peak_bound=None, hecke=True):
        length_bound = -1 if length_bound is None else length_bound
        key = (n, peak_bound, hecke)

        def generate():
            for level in ORTHOGONAL_HECKE_CACHE[key][:-1]:
                yield level
            if ORTHOGONAL_HECKE_CACHE[key]:
                level = ORTHOGONAL_HECKE_CACHE[key][-1]
            else:
                start = (cls.identity(), ())
                level = {start}
            while level:
                next_level = set()
                yield level
                for pi, w in level:
                    for i in range(n):
                        s = Permutation.s_i(i)
                        sigma = s % pi % s
                        v = w + (i,)
                        if peak_bound is None or Word.peaks(v) <= peak_bound:
                            if hecke or sigma != pi:
                                next_level.add((sigma, v))
                level = next_level

        for i, level in enumerate(generate()):
            if i >= len(ORTHOGONAL_HECKE_CACHE[key]):
                ORTHOGONAL_HECKE_CACHE[key].append(level)
            yield level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_involution_hecke_words(self, length_bound=None, peak_bound=None):
        for level in self.involution_hecke_levels(self.rank, length_bound, peak_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    @classmethod
    def _symmetrize(cls, vector, n, check=True):
        def sort(t):
            return tuple(reversed(sorted(t)))

        ans = SymmetricPolynomial({
            SymmetricMonomial(n, alpha): max(vector[a] for a in vector if sort(a) == sort(alpha))
            for alpha in vector if sort(alpha) == alpha and len(alpha) <= n
        })

        if check:
            assert all(
                vector.dictionary[sort(alpha)] == vector.dictionary[alpha]
                for alpha in vector.dictionary if len(alpha) <= n
            )
        return ans

    def stable_grothendieck(self, n, degree_bound=None):
        assert self.is_unsigned()
        ans = Vector()
        ell = self.length
        for w in self.get_hecke_words(degree_bound, n):
            ans += BETA**(len(w) - ell) * Word.quasisymmetrize(w, Word.decreasing_zeta)
        return self._symmetrize(ans, n)

    def signed_involution_stable_grothendieck(self, n, degree_bound=None):
        assert self.is_involution()
        ans = Vector()
        ell = self.involution_length
        for w in self.get_involution_hecke_words(degree_bound, n):
            ans += BETA**(len(w) - ell) * Word.quasisymmetrize(w, Word.unimodal_zeta)
        return self._symmetrize(ans, n)

    def symplectic_stable_grothendieck(self, n, degree_bound=None):
        assert self.is_unsigned() and self.is_fpf_involution()
        ans = Vector()
        ell = self.fpf_involution_length
        for w in self.get_symplectic_hecke_words(degree_bound, n):
            ans += BETA**(len(w) - ell) * Word.quasisymmetrize(w, Word.decreasing_zeta)
        return self._symmetrize(ans, n)

    def involution_stable_grothendieck(self, n, degree_bound=None):
        assert self.is_involution() and self.is_unsigned()
        ans = Vector()
        ell = self.involution_length
        for w in self.get_involution_hecke_words(degree_bound, n):
            ans += BETA**(len(w) - ell) * Word.quasisymmetrize(w, Word.decreasing_zeta)
        return self._symmetrize(ans, n)

    @property
    def reduced_word(self):
        return self.get_reduced_word()

    def get_reduced_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = Permutation.s_i(i)
            return (i,) + (s * self).get_reduced_word()
        else:
            return ()

    @property
    def involution_word(self):
        return self.get_involution_word()

    def get_involution_word(self):
        assert self.is_involution()
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = Permutation.s_i(i)
            w = self * s if self * s == s * self else s * self * s
            return w.get_involution_word() + (i,)
        else:
            return ()

    @property
    def fpf_involution_word(self):
        return self.get_fpf_involution_word()

    def get_fpf_involution_word(self):
        assert self.is_fpf_involution()
        des = [i for i in self.left_descent_set if self(i) != i + 1]
        if des:
            i = min(des)
            s = Permutation.s_i(i)
            return (s * self * s).get_fpf_involution_word() + (i,)
        else:
            return ()

    @property
    def reduced_words(self):
        return self.get_reduced_words()

    def get_reduced_words(self):
        w = self
        oneline = w.oneline
        if oneline not in SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = Permutation.s_i(i)
                for e in (w * s).get_reduced_words():
                    word = e + (i,)
                    yield word
                    words.add(word)
            SIGNED_REDUCED_WORDS[oneline] = words
        else:
            for word in SIGNED_REDUCED_WORDS[oneline]:
                yield word

    @property
    def involution_words(self):
        return self.get_involution_words()

    def get_involution_words(self):
        w = self
        assert w.inverse() == w
        oneline = w.oneline
        if oneline not in INVOLUTION_WORDS:
            ans = set()
            for a in w.get_atoms():
                for word in a.get_reduced_words():
                    yield word
                    ans.add(word)
            INVOLUTION_WORDS[oneline] = ans
        else:
            for word in INVOLUTION_WORDS[oneline]:
                yield word

    @property
    def fpf_involution_words(self):
        return self.get_fpf_involution_words()

    def get_fpf_involution_words(self):
        w = self
        assert w.is_fpf_involution()
        oneline = w.oneline
        while oneline and oneline[-2] == len(oneline) and oneline[-1] == len(oneline) - 1:
            oneline = oneline[:-2]
        w = Permutation(*oneline)
        if oneline not in FPF_INVOLUTION_WORDS:
            ans = set()
            for i in w.right_descent_set:
                if w(i) != i + 1:
                    s = Permutation.s_i(i)
                    ans |= {word + (i,) for word in (s * w * s).get_fpf_involution_words()}
            FPF_INVOLUTION_WORDS[oneline] = ans
        return FPF_INVOLUTION_WORDS[oneline]

    def __call__(self, i):
        if 0 < i <= self.rank:
            return self.oneline[i - 1]
        elif -self.rank <= i < 0:
            return -self.oneline[-i - 1]
        else:
            return i

    def __hash__(self):
        return hash(self.oneline)

    def pair(self):
        n = self.rank
        return [
            (a, self(a))
            for a in range(-n, n + 1)
            if 0 < abs(a) < self(a)
        ]

    def neg(self):
        n = self.rank
        return [(-a, -a) for a in range(1, n + 1) if self(a) == -a]

    def fix(self):
        n = self.rank
        return [(a, a) for a in range(1, n + 1) if self(a) == a]

    def cyc(self):
        return sorted(self.pair() + self.neg() + self.fix())

    @classmethod
    def identity(cls):
        return Permutation()

    @classmethod
    def longest_element(cls, n, signed=False):
        if signed:
            return Permutation(*[-i for i in range(1, n + 1)])
        else:
            return cls.longest_unsigned(n)

    @classmethod
    def longest_unsigned(cls, n):
        return Permutation(*[i for i in range(n, 0, -1)])

    @classmethod
    def s_i(cls, i):
        assert 0 <= i
        if i == 0:
            oneline = [-1]
        else:
            oneline = list(range(1, i)) + [i + 1, i]
        return Permutation(*oneline)

    @property
    def right_descent_set(self):
        if self._rdes is None:
            self._rdes = set()
            if self.rank >= 1 and self(1) < 0:
                self._rdes.add(0)
            for i in range(1, self.rank):
                if self(i) > self(i + 1):
                    self._rdes.add(i)
        return self._rdes

    @property
    def left_descent_set(self):
        if self._ldes is None:
            self._ldes = self.inverse().right_descent_set
        return self._ldes

    def __len__(self):
        if self._len is None:
            biline = [-i for i in reversed(self.oneline)] + list(self.oneline)
            n = self.rank * 2
            inv = len([(i, j) for i in range(n) for j in range(i + 1, n) if biline[i] > biline[j]])
            inv_zero = len([i for i in self.oneline if i < 0])
            assert inv % 2 == inv_zero % 2
            self._len = (inv + inv_zero) // 2
        return self._len

    @property
    def length(self):
        return len(self)

    @property
    def involution_length(self):
        return (len(self.neg()) + len(self.pair()) + len(self)) // 2

    @property
    def fpf_involution_length(self):
        return self.involution_length - self.rank // 2

    def __mod__(self, other):
        assert type(other) == Permutation
        if other.left_descent_set:
            i = next(iter(other.left_descent_set))
            s = Permutation.s_i(i)
            if i in self.right_descent_set:
                return self % (s * other)
            else:
                return (self * s) % (s * other)
        else:
            return self

    def __mul__(self, other):
        if type(other) == Tableau:
            return Tableau({(i, j): tuple(self(v) for v in value) for i, j, value in other})
        assert type(other) == Permutation
        newline = [self(other(i)) for i in range(1, max(self.rank, other.rank) + 1)]
        return Permutation(*newline)

    def inverse(self):
        newline = self.rank * [0]
        for i in range(1, self.rank + 1):
            j = self(i)
            if j > 0:
                newline[j - 1] = i
            else:
                newline[-j - 1] = -i
        return Permutation(*newline)

    def __lt__(self, other):
        assert type(other) == Permutation
        return self.oneline < other.oneline

    def __eq__(self, other):
        assert type(other) == Permutation
        return self.oneline == other.oneline

    def __pow__(self, n):
        if n < 0:
            return self.inverse().__pow__(-n)
        elif n == 0:
            return Permutation.identity(self.rank)
        elif n == 1:
            return Permutation(*self.oneline)
        else:
            p = n // 2
            q = n - p
            return self.__pow__(p) * self.__pow__(q)

    def is_even_signed(self):
        return len([i for i in self.oneline if i < 0]) % 2 == 0

    def last_descent(self):
        n = self.rank
        descents = [i for i in range(1, n) if self(i) > self(i + 1)]
        if descents:
            return max(descents)

    @classmethod
    def reflection_s(cls, i, j=None):
        if j is None:
            j = i
        assert i <= j
        caller = list(range(1, j + 1))
        caller[i - 1] = -j
        caller[j - 1] = -i
        return cls(*caller)

    @classmethod
    def transposition(cls, i, j):
        return cls.reflection_t(i, j)

    @classmethod
    def cycle(cls, *args):
        ans = Permutation()
        for i in range(len(args) - 1):
            a, b = args[i], args[i + 1]
            ans *= Permutation.transposition(a, b)
        return ans

    @classmethod
    def reflection_t(cls, i, j):
        if j < i:
            j, i = i, j
        caller = list(range(1, j + 1))
        caller[i - 1] = j
        caller[j - 1] = i
        return cls(*caller)

    def tex(self):
        s = '$'
        for i in self.oneline:
            if i > 0:
                s += str(i) + '\\hs '
            else:
                s += '\\bar' + str(-i) + '\\hs '
        s = s[:-len('\\hs ')]
        s += '$'
        return s

    def _min_inv_atom_oneline(self):
        tup = tuple(i for p in self.cyc() for i in reversed(p))
        minimum = []
        for i in tup:
            if minimum and minimum[-1] == i:
                continue
            minimum += [i]
        return tuple(minimum)

    def get_min_atom(self):
        assert self == self.inverse()
        return Permutation(*self._min_inv_atom_oneline()).inverse()

    def get_fpf_atoms(self):
        assert self.is_fpf_involution()
        assert self.is_unsigned()
        w = self
        if w not in ATOMS_FPF_CACHE:
            des = [i for i in w.right_descent_set if w(i) != i + 1]
            if des:
                ans = set()
                for i in des:
                    s = Permutation.s_i(i)
                    v = s * w * s
                    ans |= {x * s for x in v.get_fpf_atoms()}
            else:
                ans = [Permutation()]
            ATOMS_FPF_CACHE[w] = list(ans)
        return ATOMS_FPF_CACHE[w]

    def get_atoms(self):
        assert self == self.inverse()
        w = self
        if w not in atoms_b_cache:
            atoms_b_cache[w] = list(w._get_atoms())
        return atoms_b_cache[w]

    def _get_atoms(self):
        def involution(oneline):
            word = reversed(Permutation(*oneline).get_reduced_word())
            return self.from_involution_word(word)

        def next(oneline):
            y = involution(oneline)
            assert y is not None
            for i in range(len(oneline) - 2):
                c, a, b = oneline[i:i + 3]
                if a < b < c:
                    newline = oneline[:i] + (b, c, a) + oneline[i + 3:]
                    yield newline
            for i in range(len(oneline) - 1):
                b_, a_ = oneline[i:i + 2]
                a, b = -a_, -b_
                if 0 < a < b == min(map(abs, oneline[:i + 1])):
                    newline = oneline[:i] + (a, -b) + oneline[i + 2:]
                    z = involution(newline)
                    if z and y == z:
                        yield newline

        minimum = self._min_inv_atom_oneline()
        add = {minimum}
        while add:
            for w in add:
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def get_atoms_d(self):
        assert self.is_even_signed()
        assert self == self.inverse()
        w = self
        if w not in atoms_d_cache:
            atoms_d_cache[w] = list(w._get_atoms_d())
        return atoms_d_cache[w]

    def _get_atoms_d(self):
        def length(w):
            ans = 0
            for i in range(1, w.rank + 1):
                for j in range(i + 1, w.rank + 1):
                    if w(i) > w(j):
                        ans += 1
                    if -w(i) > w(j):
                        ans += 1
            return ans

        if length(self) == 0:
            yield self
            return

        def s_i(i):
            return self.s_i(i) if i != 0 else self.s_i(0) * self.s_i(1) * self.s_i(0)

        for i in range(self.rank):
            s = s_i(i)
            w = self * s
            if length(w) < length(self):
                if w == s * self:
                    for a in w.get_atoms_d():
                        yield a * s
                else:
                    for a in (s * w).get_atoms_d():
                        yield a * s

    def bruhat_covers(self):
        word = self.get_reduced_word()
        for i in range(len(word)):
            w = Permutation.from_word(word[:i] + word[i + 1:])
            if len(w) == len(self) - 1:
                yield w

    def involution_bruhat_covers(self):
        assert self.is_involution()
        word = self.get_involution_word()
        for i in range(len(word)):
            w = Permutation.from_involution_word(word[:i] + word[i + 1:])
            if w.involution_length == self.involution_length - 1:
                yield w

    def fpf_involution_bruhat_covers(self):
        assert self.is_fpf_involution()
        word = self.get_fpf_involution_word()
        for i in range(len(word)):
            w = Permutation.from_fpf_involution_word(word[:i] + word[i + 1:])
            if w.fpf_involution_length == self.fpf_involution_length - 1:
                yield w

    def tau(self, i, j):
        n = self.rank
        assert 1 <= i < j

        def t(k, l):
            if k == l:
                return Permutation()
            elif k > l:
                return Permutation.transposition(l, k)
            else:
                return Permutation.transposition(k, l)

        y = self
        z = y * t(i, y(i)) if i == y(j) else y * t(i, y(i)) * t(j, y(j))

        if y(i) <= i < j < y(j) or y(i) < i < j <= y(j) or j < y(i) < y(j) or y(i) < y(j) < i:
            s = t(i, j)
            return s * y * s

        elif y(i) == i < y(j) < j:
            s = t(i, y(j))
            return s * y * s

        elif i < y(i) < j == y(j):
            s = t(y(i), j)
            return s * y * s

        elif i < y(j) < y(i) < j:
            s = t(i, y(j))
            return s * y * s

        elif i < y(i) < y(j) < j:
            return t(i, j) * z

        elif i < y(i) < j < y(j):
            return t(i, y(j)) * z

        elif y(i) < i < y(j) < j:
            return t(y(i), j) * z

        elif i == y(i) < j == y(j):
            return t(i, j) * y

        return y

    def upper_transitions(self, r):
        """Yields j > r such that self * self.transposition(r, j) covers self in Bruhat order."""
        n = self.rank
        for i in range(r + 1, n + 1):
            t = Permutation.transposition(r, i)
            if (t * self * t).length == self.length + 1:
                yield i

    def lower_transitions(self, r):
        """Yields i < r such that self * self.transposition(i, r) covers self in Bruhat order."""
        for i in range(1, r):
            t = Permutation.transposition(i, r)
            if (t * self * t).length == self.length + 1:
                yield i

    def upper_involution_transitions(self, r):
        """Yields j > r such that self.tau(r, j) covers self in Bruhat order."""
        assert self.is_involution()
        n = self.rank
        for i in range(r + 1, n + 1):
            if self.tau(r, i).involution_length == self.involution_length + 1:
                yield i

    def lower_involution_transitions(self, r):
        """Yields i < r such that self.tau(i, r) covers self in Bruhat order."""
        assert self.is_involution()
        for i in range(1, r):
            if self.tau(i, r).involution_length == self.involution_length + 1:
                yield i

    def upper_fpf_involution_transitions(self, r):
        """
        Yields j > r such that self.transposition(r, j) * self * self.transposition(r, j)
        covers self in Bruhat order.
        """
        assert self.is_fpf_involution()
        n = self.rank
        y = self * Permutation.s_i(n + 1)
        for j in range(r + 1, n + 3):
            t = Permutation.transposition(r, j)
            if (t * y * t).fpf_involution_length == y.fpf_involution_length + 1:
                yield j

    def lower_fpf_involution_transitions(self, r):
        """
        Yields i < r such that self.transposition(i, r) * self * self.transposition(i, r)
        covers self in Bruhat order.
        """
        assert self.is_fpf_involution()
        for j in range(1, r):
            t = Permutation.transposition(j, r)
            if (t * self * t).fpf_involution_length == self.fpf_involution_length + 1:
                yield j
