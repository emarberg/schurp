from permutations import Permutation
from collections import defaultdict
from vectors import Vector
from symmetric import SchurP
from partitions import StrictPartition
import itertools


D_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_COUNTS = {(): 1}
EVEN_SIGNED_INVOLUTION_WORDS = {(): [()]}
ATOMS_D_CACHE = {}
SCHURD_STANSYM_CACHE = {}


class EvenSignedPermutation:

    def __init__(self, *oneline):
        self.oneline = tuple(oneline)
        self.rank = len(oneline)
        assert set(range(1, self.rank + 1)) == set(abs(i) for i in self.oneline)
        assert len([i for i in oneline if i < 0]) % 2 == 0
        # cached fields
        self._rdes = None
        self._ldes = None
        self._len = None
        self._involution_length = None

    def __repr__(self):
        # return 'EvenSignedPermutation(' + ', '.join([repr(i) for i in self.oneline]) + ')'
        return str(self)

    def __str__(self):
        s = []
        for i in self.oneline:
            s += [str(abs(i))]
            if i < 0:
                s += ['\u0305']
        if s:
            return ''.join(s)
        else:
            return '1'

    @classmethod
    def from_word(cls, n, *args):
        if len(args) == 1 and type(args[0]) == tuple:
            args = args[0]
        w = EvenSignedPermutation.identity(n)
        for i in args:
            w *= EvenSignedPermutation.s_i(i, n)
        return w

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            for v in range(2**n):
                oneline = []
                for i in range(n):
                    oneline.append(args[i] * (-1) ** (v % 2))
                    v = v // 2
                try:
                    yield EvenSignedPermutation(*oneline)
                except:
                    pass

    @classmethod
    def permutations(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield EvenSignedPermutation(*args)

    @classmethod
    def involutions(cls, n):
        for w in Permutation.involutions(n):
            oneline = w.oneline
            oneline += tuple(range(len(oneline) + 1, n + 1))
            cycles = [{i, oneline[i] - 1} for i in range(n) if i <= oneline[i] - 1]
            k = len(cycles)
            for v in range(2**k):
                newline = list(oneline)
                for i in range(k):
                    if v % 2 != 0:
                        for j in cycles[i]:
                            newline[j] *= -1
                    v = v // 2
                try:
                    yield EvenSignedPermutation(*newline)
                except:
                    pass

    def is_abs_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(abs(self(i)) != i for i in range(1, n + 1))

    def is_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(self(i) != i for i in range(1, n + 1))

    @classmethod
    def fpf_involutions(cls, n):
        raise NotImplementedError

    def get_reduced_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = EvenSignedPermutation.s_i(i, self.rank)
            return (i,) + (s * self).get_reduced_word()
        else:
            return ()

    def get_involution_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = EvenSignedPermutation.s_i(i, self.rank)
            if s * self == self * s:
                return (s * self).get_involution_word() + (i,)
            else:
                return (s * self * s).get_involution_word() + (i,)
        else:
            return ()

    def get_twisted_involution_word(self):
        if self.right_descent_set:
            i = min(self.right_descent_set)
            s = EvenSignedPermutation.s_i(i, self.rank)
            t = s.star()
            if t * self == self * s:
                return (self * t).get_twisted_involution_word() + (i,)
            else:
                return (t * self * s).get_twisted_involution_word() + (i,)
        else:
            return ()

    def star(self):
        oneline = [i if abs(i) != 1 else -i for i in self.oneline]
        oneline[0] *= -1
        return EvenSignedPermutation(*oneline)

    def is_twisted_involution(self):
        return self.star() == self.inverse()

    @classmethod
    def flatten(cls, word):
        return tuple(i + (i == 0) for i in word)

    def get_flattened_reduced_words(self):
        prev = None
        for word in self.get_reduced_words():
            word = self.flatten(word)
            if prev is None:
                yield word
                prev = word
            elif prev != word:
                yield word
                prev = word

    def count_flattened_involution_words(self):
        p = self.inv_stanley_schur_p_decomposition

    def get_flattened_involution_words(self):
        prev = None
        for word in self.get_involution_words():
            word = self.flatten(word)
            if prev is None:
                yield word
                prev = word
            elif prev != word:
                yield word
                prev = word

    def count_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in EVEN_SIGNED_REDUCED_COUNTS:
            words = 0
            for i in w.right_descent_set:
                s = EvenSignedPermutation.s_i(i, w.rank)
                words += (w * s).count_reduced_words()
            EVEN_SIGNED_REDUCED_COUNTS[oneline] = words
        return EVEN_SIGNED_REDUCED_COUNTS[oneline]

    def get_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in EVEN_SIGNED_REDUCED_WORDS:
            words = []
            for i in w.right_descent_set:
                s = EvenSignedPermutation.s_i(i, w.rank)
                words += [e + (i,) for e in (w * s).get_reduced_words()]
            EVEN_SIGNED_REDUCED_WORDS[oneline] = sorted(words, key=lambda x: self.flatten(x))
        return EVEN_SIGNED_REDUCED_WORDS[oneline]

    def get_signed_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in D_SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = EvenSignedPermutation.s_i(i, w.rank)
                letters = {i, -i} if i > 1 else {1} if i == 1 else {-1}
                words |= {e + (j,) for e in (w * s).get_signed_reduced_words() for j in letters}
            D_SIGNED_REDUCED_WORDS[oneline] = words
        return D_SIGNED_REDUCED_WORDS[oneline]

    def get_involution_words(self):
        w = self.reduce()
        assert w.inverse() == w
        oneline = w.oneline
        if oneline not in EVEN_SIGNED_INVOLUTION_WORDS:
            words = []
            for a in w.get_atoms():
                for word in a.get_reduced_words():
                    words.append(word)
            EVEN_SIGNED_INVOLUTION_WORDS[oneline] = sorted(words, key=lambda x: self.flatten(x))
        return EVEN_SIGNED_INVOLUTION_WORDS[oneline]

    def get_fpf_involution_words(self):
        pass

    def __call__(self, i):
        if i == 0:
            return 0
        assert 1 <= abs(i) <= self.rank
        if i > 0:
            return self.oneline[i - 1]
        else:
            return -self.oneline[abs(i) - 1]

    def __hash__(self):
        return hash(self.oneline)

    def reduce(self):
        newline = self.oneline
        while newline and newline[-1] == len(newline):
            newline = newline[:-1]
        return EvenSignedPermutation(*newline)

    def involution_length(self):
        if self._involution_length is None:
            self._involution_length = len(self.get_involution_word())
        return self._involution_length

    def __neg__(self):
        return self * EvenSignedPermutation.longest_element(self.rank)

    @classmethod
    def identity(cls, n):
        return EvenSignedPermutation(*list(range(1, n + 1)))

    @classmethod
    def longest_element(cls, n):
        if n % 2 == 0:
            return EvenSignedPermutation(*[-i for i in range(1, n + 1)])
        else:
            return EvenSignedPermutation(*([1] + [-i for i in range(2, n + 1)]))

    @classmethod
    def s_i(cls, i, n):
        assert 0 <= i < n
        if i == 0:
            oneline = [-2, -1] + list(range(3, n + 1))
        else:
            oneline = list(range(1, i)) + [i + 1, i] + list(range(i + 2, n + 1))
        return EvenSignedPermutation(*oneline)

    @property
    def right_descent_set(self):
        if self._rdes is None:
            self._rdes = set()
            for i in range(self.rank):
                s = EvenSignedPermutation.s_i(i, self.rank)
                if len(self * s) < len(self):
                    self._rdes.add(i)
        return self._rdes

    @property
    def left_descent_set(self):
        if self._ldes is None:
            self._ldes = self.inverse().right_descent_set
        return self._ldes

    def __len__(self):
        if self._len is None:
            self._len = 0
            n = self.rank
            for i in range(1, n + 1):
                for j in range(i + 1, n + 1):
                    if self(i) > self(j):
                        self._len += 1
                    if -self(i) > self(j):
                        self._len += 1
        return self._len

    def __mul__(self, other):
        assert type(other) == EvenSignedPermutation
        assert self.rank == other.rank
        newline = [self(other(i)) for i in range(1, self.rank + 1)]
        return EvenSignedPermutation(*newline)

    def inverse(self):
        newline = self.rank * [0]
        for i in range(1, self.rank + 1):
            j = self(i)
            if j > 0:
                newline[j - 1] = i
            else:
                newline[-j - 1] = -i
        return EvenSignedPermutation(*newline)

    def __lt__(self, other):
        assert type(other) == EvenSignedPermutation
        return self.oneline < other.oneline

    def __eq__(self, other):
        assert type(other) == EvenSignedPermutation
        return self.oneline == other.oneline

    def __pow__(self, n):
        if n < 0:
            return self.inverse().__pow__(-n)
        elif n == 0:
            return EvenSignedPermutation.identity(self.rank)
        elif n == 1:
            return EvenSignedPermutation(*self.oneline)
        else:
            p = n // 2
            q = n - p
            return self.__pow__(p) * self.__pow__(q)

    def inv_stanley_schur_s_decomposition(self):
        assert self == self.inverse()
        ans = Vector()
        for x in self.get_atoms():
            ans += x.stanley_schur_p_decomposition()
        return SchurP.decompose_s_lambda(ans)

    def inv_stanley_schur_p_decomposition(self):
        assert self == self.inverse()
        ans = Vector()
        for x in self.get_atoms():
            ans += x.stanley_schur_p_decomposition()
        return ans

    def stanley_schur_p_decomposition(self):
        ans = Vector()
        for sh, i in self.stanley_schur_decomposition().items():
            ans += Vector({SchurP(StrictPartition(*sh)): i})
        return ans

    def stanley_schur_s_decomposition(self):
        ans = self.stanley_schur_p_decomposition()
        return SchurP.decompose_s_lambda(ans)

    def stanley_schur_decomposition(self,):
        cache = SCHURD_STANSYM_CACHE
        w = self.reduce()
        n = w.rank

        if w in cache:
            return cache[w]

        sh = w.increasing_shape()
        if sh is not None:
            return {sh: 1}

        r = w.last_descent()
        s = w.last_inversion(r)
        v = w * self.reflection_t(r, s, n)

        v_len = len(v)
        assert v_len + 1 == len(w)
        indices = [v * self.reflection_t(i, r, n) for i in range(1, r)]

        newline = v.oneline + (n + 1,)
        v = EvenSignedPermutation(*newline)
        indices += [v * self.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
        indices = [x.reduce() for x in indices if v_len + 1 == len(x)]

        ans = defaultdict(int)
        for x in indices:
            for sh, i in x.stanley_schur_decomposition().items():
                ans[sh] += i
        ans = dict(ans)

        cache[self] = ans
        return ans

    def increasing_shape(self):
        if all(self(i) < self(i + 1) for i in range(1, self.rank)):
            return tuple(abs(i) for i in self.oneline if i < 0)

    def last_descent(self):
        n = self.rank
        descents = [i for i in range(1, n) if self(i) > self(i + 1)]
        if descents:
            return max(descents)

    def last_inversion(self, r):
        n = self.rank
        return max(s for s in range(r + 1, n + 1) if self(s) < self(r))

    @classmethod
    def reflection_s(cls, i, j, n):
        assert i != j
        caller = list(range(1, n + 1))
        caller[i - 1] = -j
        caller[j - 1] = -i
        return EvenSignedPermutation(*caller)

    @classmethod
    def reflection_t(cls, i, j, n):
        assert i != j
        caller = list(range(1, n + 1))
        caller[i - 1] = j
        caller[j - 1] = i
        return EvenSignedPermutation(*caller)

    def inflate(self, rank):
        newline = self.oneline + tuple(range(self.rank + 1, rank + 1))
        return EvenSignedPermutation(*newline)

    def get_atoms(self):
        assert self == self.inverse()
        w = self.reduce()
        if w not in ATOMS_D_CACHE:
            ATOMS_D_CACHE[w] = list(w._get_atoms())
        ans = ATOMS_D_CACHE[w]
        return [x.inflate(self.rank) for x in ans]

    def _get_atoms(self):
        if len(self) == 0:
            yield self
            return

        for i in range(self.rank):
            s = EvenSignedPermutation.s_i(i, self.rank)
            w = self * s
            if len(w) < len(self):
                if w == s * self:
                    for a in w.get_atoms():
                        yield a * s
                else:
                    for a in (s * w).get_atoms():
                        yield a * s
