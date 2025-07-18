import itertools
from collections import defaultdict
from vectors import Vector
from symmetric import SchurP, SchurQ
from partitions import StrictPartition
from permutations import Permutation

DEMAZURE_FACTORIZATIONS = {}

B_SIGNED_REDUCED_WORDS = {(): {()}}
C_SIGNED_REDUCED_WORDS = {(): {()}}

SIGNED_REDUCED_WORDS = {(): {()}}
DTYPE_REDUCED_WORDS = {(): {()}}

SIGNED_HECKE_WORDS = {}
DTYPE_HECKE_WORDS = {}

B_COMPATIBLE_SEQUENCES = {}
C_COMPATIBLE_SEQUENCES = {}
D_COMPATIBLE_SEQUENCES = {}

PEAK_CSEQ_CACHE = {}

SIGNED_INVOLUTION_WORDS = {}
SIGNED_FPF_INVOLUTION_WORDS = {}

schurp_stansym_cache = {}
schurq_stansym_cache = {}
schurd_stansym_cache = {}
atoms_b_cache = {}
atoms_d_cache = {}


class SignedMixin:

    def neg(self):
        n = self.rank
        return [(-a, -a) for a in range(1, n + 1) if self(a) == -a]

    def fix(self):
        n = self.rank
        return [(a, a) for a in range(1, n + 1) if self(a) == a]

    def ell_zero(self):
        return len([i for i in range(1, self.rank + 1) if self(i) < 0])

    @classmethod
    def one_fpf_d(cls, n):
        if n % 2 == 0:
            return cls(*[a for i in range(0, n, 2) for a in [i + 2, i + 1]])
        else:
            return cls(*([1] + [a for i in range(1, n, 2) for a in [i + 2, i + 1]]))

    def __abs__(self):
        return Permutation(*[abs(self(i)) for i in range(1, self.rank + 1)])

    @property
    def ranktable(self):
        if self._ranks is None:
            w = self.inflate(self.rank + 1)
            n = w.rank
            self._ranks = {(i, j): 0 for i in range(-n, n + 1) for j in range(-n, n + 1)}
            for a in range(-n, n + 1):
                b = w(a)
                for i in range(a, n + 1):
                    for j in range(-n, b + 1):
                        self._ranks[i, j] += 1
        return self._ranks

    @property
    def emptyranks(self):
        if self._empty is None:
            n = self.rank
            self._empty = set()
            c = None
            for a in range(1, n + 1):
                c = abs(self(a)) if c is None else min(c, abs(self(a)))
                for b in range(1, c):
                    self._empty.add((a, b))
        return self._empty

    def embedded_permutation(self):
        n = self.rank
        t = lambda x: (x + n + 1) if x < 0 else x + n
        oneline = [t(self(i)) for i in range(-n, n + 1) if i != 0]
        return Permutation(*oneline)

    def dbruhat_less_equal(self, other):
        assert self.is_even_signed()
        assert other.is_even_signed()
        n = self.rank
        assert n == other.rank
        v = self
        u = other
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                if not (v.ranktable[i, j] <= u.ranktable[i, j]):
                    return False
        for (a, b) in v.emptyranks & u.emptyranks:
            if u.ranktable[-a - 1, b + 1] == v.ranktable[-a - 1, b + 1]:
                if (u.ranktable[-1, b + 1] - v.ranktable[-1, b + 1]) % 2 != 0:
                    return False
        return True

    def dbruhat_less_than(self, other):
        return self != other and self.dbruhat_less_equal(other)

    def strong_bruhat_less_equal(self, other, dtype=False):
        return self == other or self.strong_bruhat_less_than(other, dtype)

    def strong_bruhat_less_than(self, other, dtype=False):
        if dtype:
            return self.dbruhat_less_than(other)
        else:
            return self.embedded_permutation().strong_bruhat_less_than(other.embedded_permutation())

    def pair(self):
        n = self.rank
        return sorted([
            (a, self(a))
            for a in range(-n, n + 1)
            if 0 < abs(a) < self(a)
        ])

    def length(self):
        return len(self)

    def dlength(self):
        return len(self) - self.ell_zero()

    @classmethod
    def identity(cls, n):
        return cls(*list(range(1, n + 1)))

    def __mod__(self, other):
        assert type(other) == type(self)
        assert self.rank == other.rank
        if other.left_descent_set:
            i = next(iter(other.left_descent_set))
            s = self.s_i(i, self.rank)
            if i in self.right_descent_set:
                return self % (s * other)
            else:
                return (self * s) % (s * other)
        else:
            return self

    def dtype_demazure(self, other):
        assert type(other) == type(self)
        assert self.rank == other.rank
        if other.dtype_left_descent_set:
            i = next(iter(other.dtype_left_descent_set))
            s = self.ds_i(i, self.rank)
            if i in self.dtype_descent_set:
                return self.dtype_demazure(s * other)
            else:
                return (self * s).dtype_demazure(s * other)
        else:
            return self

    @classmethod
    def fpf_class(cls, n, dtype=False):
        if n % 2 == 0:
            oneline = [(i + 2) if i % 2 == 0 else i for i in range(n)]
            w = SignedPermutation(*oneline)
            return w.conjugacy_class(dtype=dtype, twisted=False)
        else:
            oneline = [1 if dtype else -1] + [(i + 3) if i % 2 == 0 else (i + 1) for i in range(n - 1)]
            w = SignedPermutation(*oneline)
            return w.conjugacy_class(dtype=dtype, twisted=dtype)

    def twisted_conjugacy_class(self):
        return self.conjugacy_class(True, True)

    def conjugacy_class(self, dtype=False, twisted=False):
        n = self.rank
        if twisted:
            assert dtype
            t = self.s_i(0, n)
        else:
            t = self.s_i(0, n)**2
        if dtype:
            g = [self.ds_i(i, n) for i in range(1, n)]
            if n >= 2:
                g += [self.ds_i(-1, n)]
        else:
            g = [self.s_i(i, n) for i in range(n)]
        ans = set()
        add = {self}
        while add:
            newadd = set()
            ans |= add
            for s in g:
                for a in add:
                    sas = t * s * t * a * s
                    if sas not in ans:
                        newadd.add(sas)
            add = newadd
        return ans

    def __mul__(self, other):
        assert type(other) == type(self)
        assert self.rank == other.rank
        newline = [self(other(i)) for i in range(1, self.rank + 1)]
        return self.__class__(*newline)

    def inverse(self):
        newline = self.rank * [0]
        for i in range(1, self.rank + 1):
            j = self(i)
            if j > 0:
                newline[j - 1] = i
            else:
                newline[-j - 1] = -i
        return self.__class__(*newline)

    def __lt__(self, other):
        assert type(other) == type(self)
        return self.oneline < other.oneline

    def __eq__(self, other):
        if type(other) != type(self):
            return False
        return self.oneline == other.oneline

    def __pow__(self, n):
        if n < 0:
            return self.inverse().__pow__(-n)
        elif n == 0:
            return self.identity(self.rank)
        elif n == 1:
            return self.__class__(*self.oneline)
        else:
            p = n // 2
            q = n - p
            return self.__pow__(p) * self.__pow__(q)

    @classmethod
    def t_ij(cls, i, j, n):
        return cls.reflection(i, j, n)

    @classmethod
    def reflection(cls, i, j, n):
        caller = list(range(1, n + 1))
        if i > 0:
            caller[i - 1] = j
        else:
            caller[-i - 1] = -j
        if j > 0:
            caller[j - 1] = i
        else:
            caller[-j - 1] = -i
        return cls(*caller)

    @classmethod
    def reflection_s(cls, i, j, n):
        caller = list(range(1, n + 1))
        caller[i - 1] = -j
        caller[j - 1] = -i
        return cls(*caller)

    @classmethod
    def reflection_t(cls, i, j, n):
        assert i != j
        caller = list(range(1, n + 1))
        caller[i - 1] = j
        caller[j - 1] = i
        return cls(*caller)

    def inflate(self, rank):
        newline = self.oneline + tuple(range(self.rank + 1, rank + 1))
        return self.__class__(*newline)

    def fixed_points(self):
        n = self.rank
        return tuple(a for a in range(-n, n + 1) if 0 != a == self(a))

    def negated_points(self):
        n = self.rank
        return tuple(a for a in range(-n, n + 1) if 0 != a == -self(a))

    def __call__(self, i):
        if not (1 <= abs(i) <= self.rank):
            return i
        if i > 0:
            return self.oneline[i - 1]
        else:
            return -self.oneline[abs(i) - 1]

    def __hash__(self):
        return hash(self.oneline)

    def reduce(self, limit=0):
        newline = self.oneline
        while newline and len(newline) > limit and newline[-1] == len(newline):
            newline = newline[:-1]
        return self.__class__(*newline)

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

    @classmethod
    def from_word(cls, n, *args):
        if len(args) == 1 and type(args[0]) == tuple:
            args = args[0]
        w = cls.identity(n)
        for i in args:
            w *= cls.s_i(i, n)
        return w

    @classmethod
    def double_involution_word(cls, word):
        if len(word) == 0:
            return tuple()
        n = 1 + max([i for i in word])
        w = cls.identity(n)
        ans = []
        for i in word:
            assert i not in w.right_descent_set
            s = cls.s_i(i, n)
            w *= s
            if i not in w.left_descent_set:
                ans = [i] + ans + [i]
                w = s * w
            else:
                ans = ans + [i]
        return tuple(ans)

    def __neg__(self):
        return self.longest_element(self.rank) * self

    def count_nontrivial_positive_cycles(self):
        return len([c for c in self.cycles() if min(c) > 0 and len(c) > 1])

    def count_positive_ascents(self):
        return len([i for i in range(1, self.rank + 1) if i > self(i)])

    def cycles(self):
        ans = []
        a = list(range(1, 1 + self.rank))
        while a:
            cyc = [a[0]]
            while self(cyc[-1]) != cyc[0]:
                cyc.append(self(cyc[-1]))
            ans.append(tuple(cyc))
            ncyc = [-i for i in cyc]
            if set(ncyc) != set(cyc):
                ans.append(tuple(ncyc))
            a = sorted(set(a) - set(cyc) - set(ncyc))
        return ans

    def cycle_repr(self):
        if len(self) == 0:
            return '1'

        SPACE = ' '
        DELIM = ''  # ','

        s = ''
        for c in self.cycles():
            # if len(c) == 1 and c[0] > 0:
            #    s += '%s' % c[0]
            # elif len(c) == 2 and c[0] == -c[1]:
            #    s += '%s' % (str(abs(c[0])) + '\u0305')
            # elif c[0] > 0:
            #    s += '(' + (DELIM + SPACE).join([(str(-x) + '\u0305') if x < 0 else str(x) for x in c]) + ')'
            if len(c) > 1:
                s += '(' + (DELIM + SPACE).join([str(x) for x in c]) + ')'
        return s.strip()

    def __iter__(self):
        return self.oneline.__iter__()
        
    def __str__(self):
        s = []
        for i in self.oneline:
            #s += [str(abs(i)) + ('\u0305' if i < 0 else '')]
            s += [str(i) if i < 0 else ('+' + str(i))]
        if s:
            sep = '' if len(self.oneline) < 10 else ' '
            return '[' + sep.join(s) + ']'
        else:
            return '[]'

    @classmethod
    def permutations(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield cls(*args)

    def get_reduced_word(self, dtype=False):
        des = self.dtype_descent_set if dtype else self.right_descent_set
        if des:
            i = min(des)
            s = self.ds_i(i, self.rank) if dtype else self.s_i(i, self.rank)
            return (self * s).get_reduced_word(dtype) + (i,)
        else:
            return ()

    def _get_reduced_words(self, des, cache):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in cache:
            words = set()
            for i in des(w):
                s = w.s_i(i, w.rank)
                words |= {e + (i,) for e in (w * s)._get_reduced_words(des, cache)}
            cache[oneline] = words
        return cache[oneline]

    def is_identity(self):
        return len(self) == 0


class SignedPermutation(SignedMixin):

    def __init__(self, *oneline):
        self.oneline = tuple(oneline)
        self.rank = len(oneline)
        assert set(range(1, self.rank + 1)) == set(abs(i) for i in self.oneline)
        # cached fields
        self._rdes = None
        self._ldes = None
        self._ddes = None
        self._len = None
        self._ranks = None
        self._empty = None

    def __repr__(self):
        # return 'SignedPermutation(' + ', '.join([repr(i) for i in self.oneline]) + ')'
        return str(self)

    def get_demazure_factorizations(self):
        w = self.reduce()
        key = tuple(w.oneline)
        n = len(key)

        if key not in DEMAZURE_FACTORIZATIONS:
            dictionary = {}
            for u in SignedPermutation.all(n):
                for v in Permutation.all(n):
                    x = u % SignedPermutation(*[v(i) for i in range(1, 1 + n)])
                    xkey = tuple(x.reduce().oneline)
                    if xkey not in dictionary:
                        dictionary[xkey] = set()
                    dictionary[xkey].add((u.reduce(), v))
            for xkey in dictionary:
                if xkey not in DEMAZURE_FACTORIZATIONS:
                    DEMAZURE_FACTORIZATIONS[xkey] = dictionary[xkey]

        return DEMAZURE_FACTORIZATIONS[key]

    @classmethod
    def get_grassmannians_bc(cls, n):
        for k in range(n + 1):
            for s in itertools.combinations(list(range(1, n + 1)), k):
                t = sorted(set(range(1, n + 1)) - set(s))
                mu = tuple(reversed(s))
                oneline = [-m for m in mu] + t
                yield (cls(*oneline), mu)

    @classmethod
    def get_grassmannians_d(cls, n):
        for k in range(0, n + 1, 2):
            for s in itertools.combinations(list(range(1, n + 1)), k):
                t = sorted(set(range(1, n + 1)) - set(s))
                mu = tuple(reversed([i - 1 for i in s]))
                if mu and mu[-1] == 0:
                    mu = mu[:-1]
                oneline = [-m for m in reversed(s)] + t
                yield (cls(*oneline), mu)

    @classmethod
    def reflections(cls, n):
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                yield cls.reflection_t(i, j, n)
                yield cls.reflection_s(i, j, n)
            yield cls.reflection_s(i, i, n)

    @classmethod
    def all(cls, n, dtype=False):
        for args in itertools.permutations(range(1, n + 1)):
            for v in range(2**n):
                oneline = []
                for i in range(n):
                    oneline.append(args[i] * (-1) ** (v % 2))
                    v = v // 2
                w = cls(*oneline)
                if not dtype or w.is_even_signed():
                    yield w

    @classmethod
    def involutions(cls, n, dtype=False, twisted=False):
        if twisted:
            assert dtype
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
                w = cls(*newline)
                if not dtype:                 
                    yield w
                elif dtype and not twisted and w.is_even_signed():
                    yield w
                elif dtype and twisted and not w.is_even_signed():
                    yield cls.s_i(0, n) * w

    @classmethod
    def abs_fpf_involutions(cls, n):
        for w in cls.involutions(n):
            if w.is_abs_fpf_involution():
                yield w

    @classmethod
    def fpf_involutions(cls, n):
        for w in cls.involutions(n):
            if w.is_fpf_involution():
                yield w

    def is_involution(self):
        return self == self.inverse()

    def is_abs_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(abs(self(i)) != i for i in range(1, n + 1))

    def is_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(self(i) != i for i in range(1, n + 1))

    def get_reduced_words(self):
        des = lambda w: w.right_descent_set
        return self._get_reduced_words(des, SIGNED_REDUCED_WORDS)

    def min_peaks(self):
        ans = None
        for w in self.get_reduced_words():
            a = len([i for i in range(1, len(w) - 1) if w[i - 1] < w[i] > w[i + 1]])
            ans = a if (ans is None or a < ans) else ans
        return 0 if ans is None else ans

    def get_dtype_reduced_words(self):
        des = lambda w: w.dtype_descent_set
        return self._get_reduced_words(des, DTYPE_REDUCED_WORDS)

    def dtype_min_peaks(self):
        ans = None
        for w in self.get_dtype_reduced_words():
            a = len([i for i in range(1, len(w) - 1) if abs(w[i - 1]) <= abs(w[i]) >= abs(w[i + 1])])
            ans = a if (ans is None or a < ans) else ans
        return 0 if ans is None else ans

    def _get_hecke_words(self, length, cache, des):
        w = self.reduce()
        oneline = tuple(w.oneline)
        key = (oneline, length)
        if key not in cache:
            if length < 0:
                words = set()
            elif length == 0:
                words = {()} if w.is_identity() else set()
            else:
                words = {e + (i,) for e in w._get_hecke_words(length - 1, cache, des) for i in des(w)}
                for i in des(w):
                    s = w.s_i(i, w.rank)
                    words |= {e + (i,) for e in (w * s)._get_hecke_words(length - 1, cache, des)}
            cache[key] = words
        return cache[key]

    def get_hecke_words(self, length):
        des = lambda w: w.right_descent_set
        return self._get_hecke_words(length, SIGNED_HECKE_WORDS, des)

    def get_dtype_hecke_words(self, length):
        des = lambda w: w.dtype_descent_set
        return self._get_hecke_words(length, DTYPE_HECKE_WORDS, des)

    @classmethod
    def _get_peak_compatible_sequences(cls, n, ascents, peaks):
        key = (n, ascents, peaks)
        cache = PEAK_CSEQ_CACHE
        if key not in cache:
            if len(peaks) == 0:
                ans = {()}
            else:
                ans = set()
                bns = cls._get_peak_compatible_sequences(n, ascents[:-1], peaks[:-1])
                for b in bns:
                    strict = (ascents[-1] and len(b) >= 1) or (peaks[-1] and len(b) >= 2 and b[-2] == b[-1])
                    start = 1 if len(b) == 0 else (b[-1] + 1) if strict else b[-1]
                    for e in range(start, n + 1):
                        ans.add(b + (e,))
            cache[key] = ans
        return cache[key]

    def _get_hecke_compatible_sequences(self, n, length, cache, get_hecke_words, get_ascents, get_peaks):
        w = self.reduce()
        key = (w.oneline, n, length)
        if key not in cache:
            ans = {}
            for ell in range(length + 1):
                for a in get_hecke_words(w, ell):
                    ans[a] = self._get_peak_compatible_sequences(n, get_ascents(a), get_peaks(a))
            cache[key] = ans
        return cache[key]

    def get_hecke_compatible_sequences_b(self, n, length):
        get_hecke_words = lambda w, ell: w.get_hecke_words(ell)
        get_ascents = lambda a: tuple(0 < i < len(a) and a[i - 1] == a[i] == 0 for i in range(len(a)))
        get_peaks = lambda a: tuple(1 < i < len(a) and a[i - 2] <= a[i - 1] >= a[i] for i in range(len(a)))
        return self._get_hecke_compatible_sequences(n, length, B_COMPATIBLE_SEQUENCES, get_hecke_words, get_ascents, get_peaks)
    
    def get_hecke_compatible_sequences_c(self, n, length):
        get_hecke_words = lambda w, ell: w.get_hecke_words(ell)
        get_ascents = lambda a: len(a) * (False,)
        get_peaks = lambda a: tuple(1 < i < len(a) and a[i - 2] <= a[i - 1] >= a[i] for i in range(len(a)))
        return self._get_hecke_compatible_sequences(n, length, C_COMPATIBLE_SEQUENCES, get_hecke_words, get_ascents, get_peaks)

    def get_hecke_compatible_sequences_d(self, n, length):
        get_hecke_words = lambda w, ell: w.get_dtype_hecke_words(ell)
        get_ascents = lambda a: tuple(0 < i < len(a) and a[i - 1] == a[i] and abs(a[i]) == 1 for i in range(len(a)))
        get_peaks = lambda a: tuple(1 < i < len(a) and abs(a[i - 2]) <= abs(a[i - 1]) >= abs(a[i]) for i in range(len(a)))
        return self._get_hecke_compatible_sequences(n, length, D_COMPATIBLE_SEQUENCES, get_hecke_words, get_ascents, get_peaks)
        
    def get_signed_reduced_words(self):
        return self.get_type_b_reduced_words()

    def get_type_b_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in B_SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = SignedPermutation.s_i(i, w.rank)
                letters = {i, -i}
                words |= {e + (j,) for e in (w * s).get_type_b_reduced_words() for j in letters}
            B_SIGNED_REDUCED_WORDS[oneline] = words
        return B_SIGNED_REDUCED_WORDS[oneline]

    def get_type_c_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in C_SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = SignedPermutation.s_i(i, w.rank)
                letters = {i + 1, -i - 1}
                words |= {e + (j,) for e in (w * s).get_type_c_reduced_words() for j in letters}
            C_SIGNED_REDUCED_WORDS[oneline] = words
        return C_SIGNED_REDUCED_WORDS[oneline]

    def get_fpf_involution_words(self):
        assert self.is_fpf_involution()
        w = self.reduce()
        n = w.rank

        oneline = w.oneline
        if oneline not in SIGNED_FPF_INVOLUTION_WORDS:
            words = set()
            if len(w) == (n + 1) // 2:
                words.add(())
            else:
                for i in w.right_descent_set:
                    if i == 0 and w(1) == -1:
                        continue
                    if i > 0 and w(i) == i + 1:
                        continue
                    s = self.s_i(i, w.rank)
                    if s * w != w * s:
                        v = s * w * s
                    else:
                        v = w * s
                    for word in v.get_fpf_involution_words():
                        words.add(word + (i,))
            SIGNED_FPF_INVOLUTION_WORDS[oneline] = words
        return SIGNED_FPF_INVOLUTION_WORDS[oneline]

    def get_involution_words(self):
        assert self.is_involution()
        w = self.reduce()
        oneline = w.oneline
        if oneline not in SIGNED_INVOLUTION_WORDS:
            words = set()
            if len(self.left_descent_set) == 0:
                words.add(())
            else:
                for i in self.left_descent_set:
                    s = self.s_i(i, self.rank)
                    w = s * self
                    if w != self * s:
                        w = w * s
                    for e in w.get_involution_words():
                        words.add(e + (i,))
            SIGNED_INVOLUTION_WORDS[oneline] = words
        return SIGNED_INVOLUTION_WORDS[oneline]

    @classmethod
    def involution_hecke_words(cls, n, length_bound=-1):
        for level in cls.involution_hecke_levels(n, length_bound):
            for pi, w in level:
                yield w

    @classmethod
    def involution_hecke_levels(cls, n, length_bound=-1):
        start = (cls.identity(n), ())
        level = {start}
        while level:
            next_level = set()
            yield level
            for pi, w in level:
                for i in range(n):
                    s = cls.s_i(i, n)
                    sigma = s % pi % s
                    next_level.add((sigma, w + (i,)))
            level = next_level
            if length_bound == 0:
                break
            length_bound -= 1

    def get_involution_hecke_words(self, length_bound):
        for level in self.involution_hecke_levels(self.rank, length_bound):
            for pi, w in level:
                if self == pi:
                    yield w

    def get_involution_word(self):
        assert self.is_involution()
        if self.right_descent_set:
            i = min(self.right_descent_set)
            s = self.s_i(i, self.rank)
            w = self * s
            if w != s * self:
                w = s * w
            return w.get_involution_word() + (i,)
        else:
            return ()

    def involution_length(self, dtype=False, twisted=False):
        if twisted:
            assert dtype
            return next(iter(self.get_atoms_d(twisted=True))).dlength()
        elif dtype:
            return next(iter(self.get_atoms_d())).dlength()
        else:
            return (len(self.neg()) + len(self.pair()) + len(self)) // 2

    @classmethod
    def dtype_longest_element(cls, n):
        if n % 2 == 0:
            return cls.longest_element(n)
        else:
            return cls.s_i(0, n) * cls.longest_element(n)

    @classmethod
    def longest_element(cls, n, k=None):
        k = n if k is None else k
        return SignedPermutation(*[(-i if i <= k else i) for i in range(1, n + 1)])

    @classmethod
    def longest_fpf(cls, n):
        oneline = [-j for i in range(1, n, 2) for j in [i + 1, i]]
        if n % 2 != 0:
            oneline += [-n]
        return SignedPermutation(*oneline)

    @classmethod
    def grassmannian_element(cls, n, k=None):
        k = n if k is None else k
        return SignedPermutation(*[(-(k + 1 - i) if i <= k else i) for i in range(1, n + 1)])

    @classmethod
    def s_i(cls, i, n):
        assert -1 <= i < n
        if i == -1:
            assert n >= 2
            oneline = [-2, -1] + list(range(3, n + 1))
        elif i == 0:
            oneline = [-1] + list(range(2, n + 1))
        else:
            oneline = list(range(1, i)) + [i + 1, i] + list(range(i + 2, n + 1))
        return SignedPermutation(*oneline)

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
    def dtype_left_descent_set(self):
        return self.inverse().dtype_descent_set

    @property
    def dtype_descent_set(self):
        if self._ddes is None:
            self._ddes = set()
            if self.rank >= 2 and -self(2) > self(1):
                self._ddes.add(-1)
            for i in range(1, self.rank):
                if self(i) > self(i + 1):
                    self._ddes.add(i)
        return self._ddes

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

    def is_even_signed(self):
        return len([i for i in self.oneline if i < 0]) % 2 == 0

    def fpf_stanley_schur_s_decomposition(self):
        assert self.is_fpf_involution()
        ans = Vector()
        for x in self.get_fpf_atoms():
            ans += x.stanley_schur_q_decomposition()
        return SchurQ.decompose_s_lambda(ans)

    def fpf_stanley_schur_p_decomposition(self):
        assert self.is_fpf_involution()
        ans = Vector()
        for x in self.get_fpf_atoms():
            ans += x.stanley_schur_p_decomposition()
        return ans

    def fpf_stanley_schur_q_decomposition(self):
        assert self.is_fpf_involution()
        ans = Vector()
        for x in self.get_fpf_atoms():
            ans += x.stanley_schur_q_decomposition()
        return ans

    def inv_stanley_schur_s_decomposition(self):
        assert self == self.inverse()
        ans = Vector()
        for x in self.get_atoms():
            ans += x.stanley_schur_q_decomposition()
        return SchurQ.decompose_s_lambda(ans)

    def inv_stanley_schur_p_decomposition(self):
        assert self == self.inverse()
        ans = Vector()
        for x in self.get_atoms():
            ans += x.stanley_schur_p_decomposition()
        return ans

    def inv_stanley_schur_q_decomposition(self):
        assert self == self.inverse()
        ans = Vector()
        for x in self.get_atoms():
            ans += x.stanley_schur_q_decomposition()
        return ans

    def stanley_schur_p_decomposition(self):
        ans = Vector()
        for sh, i in self.stanley_schur_decomposition('B').items():
            ans += Vector({SchurP(StrictPartition(*sh)): i})
        return ans

    def stanley_schur_q_decomposition(self):
        ans = Vector()
        for sh, i in self.stanley_schur_decomposition('C').items():
            ans += Vector({SchurQ(StrictPartition(*sh)): i})
        return ans

    def stanley_schur_s_decomposition(self):
        ans = self.stanley_schur_q_decomposition()
        return SchurQ.decompose_s_lambda(ans)

    def stanley_schur_d_decomposition(self):
        ans = Vector()
        for sh, i in self.stanley_schur_decomposition('D').items():
            ans += Vector({SchurP(StrictPartition(*sh)): i})
        return ans

    def _get_cache(self, bcd_type):
        assert bcd_type in ['B', 'C', 'D']
        if bcd_type == 'D':
            assert self.is_even_signed()
        caches = {
            'B': schurp_stansym_cache,
            'C': schurq_stansym_cache,
            'D': schurd_stansym_cache
        }
        return caches[bcd_type]

    def stanley_schur_decomposition(self, bcd_type):
        cache = self._get_cache(bcd_type)
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

        indices = []
        if bcd_type == 'C':
            indices += [v * self.reflection_s(r, r, n)]
        elif bcd_type == 'B':
            indices += 2 * [v * self.reflection_s(r, r, n)]
        elif bcd_type == 'D':
            pass
        indices += [v * self.reflection_t(i, r, n) for i in range(1, r)]

        newline = v.oneline + (n + 1,)
        v = SignedPermutation(*newline)
        indices += [v * self.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
        indices = [x.reduce() for x in indices if v_len + 1 == len(x)]

        ans = defaultdict(int)
        for x in indices:
            for sh, i in x.stanley_schur_decomposition(bcd_type).items():
                ans[sh] += i
        ans = dict(ans)

        cache[self] = ans
        return ans

    @classmethod
    def standardize(cls, oneline):
        distinct = sorted([abs(j) for j in oneline])
        assert len(set(distinct)) == len(oneline)
        mapping = {v: i + 1 for i, v in enumerate(distinct)}
        newline = tuple(mapping[v] if v > 0 else -mapping[-v] for v in oneline)
        return cls(*newline)

    @classmethod
    def get_grassmannian(cls, n):
        for w in cls.all(n):
            if w.is_grassmannian():
                yield w

    @classmethod
    def get_inv_grassmannian(cls, n):
        ans = defaultdict(list)
        for w in cls.involutions(n):
            s = w.inv_stanley_schur_s_decomposition()
            keys = list(s.keys())
            if len(keys) == 1 and s[keys[0]] == 1:
                ans[keys[0].mu.compact()] += [w]
        for k, v in sorted(ans.items(), key=lambda x: (len(x[0]), str(x[0]))):
            if v:
                print('\n' + k)
                for w in v:
                    print('  ', w.reduce(), w.get_min_atom().reduce())
        return ans

    def is_inv_grassmannian(self):
        return self == self.inverse() and self.get_min_atom().is_grassmannian()

    def is_grassmannian(self):
        return self.increasing_shape() is not None

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

    def cyc(self, sh=None):
        pair = self.pair()
        if sh is not None:
            pair += [(-b, a) for (a, b) in sh if 0 < a < b]
        neg = [x for x in self.neg() if sh is None or (x[0], -x[0]) in sh]
        fix = self.fix()
        return sorted(pair + neg + fix)

    @classmethod
    def ncsp_matchings(cls, base, trivial_allowed=True):
        base = set(base)
        assert all(-i in base and i != 0 for i in base)

        if len(base) == 0:
            yield ()
            return

        x = min(base)
        neg = {y for y in base if x < y < 0}
        for y in neg:
            left = {z for z in base if x < z < y or x < -z < y}
            right = {z for z in base if y < z < -y}
            for a in cls.ncsp_matchings(left, False):
                for b in cls.ncsp_matchings(right, trivial_allowed):
                    yield tuple(sorted(set(a) | set(b) | {(x, y), (-y, -x)}))
        if trivial_allowed:
            for a in cls.ncsp_matchings(base - {x, -x}, True):
                yield tuple(sorted(set(a) | {(x, -x)}))

    def fpf_shape(self, offset=None):
        n = self.rank
        if offset is None:
            offset = 0 if n % 2 == 0 else 1
        assert 0 <= offset <= n and (n - offset) % 2 == 0
        oneline = [-k for k in range(offset, 0, -1)]
        for i in range(offset + 1, n, 2):
            oneline += [i + 1, i]
        ans = (SignedPermutation(*oneline) * self).shape()

        o = list(self.inverse().oneline)
        sh = {(-abs(a), abs(a)) for a in o[:offset]}
        
        o = o[offset:]
        while o:
            b, c = o[:2]
            o = o[2:]
            if 0 < c < -b:
                sh.add((c, -b))
                sh.add((b, -c))

        assert sh == ans
        return sh

    def dtype_demazure_conjugate(self, twisted):
        if not twisted:
            return self.inverse().dtype_demazure(self)
        else:
            t = SignedPermutation.s_i(0, self.rank)
            return (t * self.inverse() * t).dtype_demazure(self)
            
    def is_atom_d(self, twisted):
        y = self.dtype_demazure_conjugate(twisted)
        return y.involution_length(dtype=True, twisted=twisted) == self.dlength()

    def fpf_dshape(self, verbose=False):
        n = self.rank
        w = self.one_fpf_d(n) * self
        twisted = n % 2 != 0
        assert w.is_atom_d(twisted)

        offset = 1 if twisted else 0
        o = list(w.inverse().oneline)
        init = [abs(i) for i in o[:offset]]
        ndes, fix, neg = w._ndes(o[offset:])
        
        desd = [(b, a) for a, b in ndes if a >= -b]
        negd = init + [abs(a) for a, b in ndes if a < -b] + [-b for a, b in ndes if a < -b]
        
        sh = set()
        for a, b in ndes:
            if 0 < a < -b:
                sh.add((b, -a))
                sh.add((a, -b))
            elif a < 0 < -b:
                sh.add((b, -b))
                sh.add((a, -a))
        for a in init:
            sh.add((-a, a))

        y = w.dtype_demazure_conjugate(twisted)
        z = SignedPermutation.identity(n)
        if twisted:
            z = z * SignedPermutation.s_i(0, n)
        for i in negd:
            z *= SignedPermutation.t_ij(-i, i, n)
        for a, b in desd:
            z *= SignedPermutation.t_ij(a, b, n)

        if verbose:
            print()
            print('*', 'y =', y, 'w =', w.inverse())
            print(' ', 'z =', z)
            print()
            print(' ndes =', ndes)
            print('  fix =', fix)
            print('  neg =', neg)
            print()
            print(' desd =', desd)
            print(' negd =', negd)
            print()
            print('   sh =', sh)

        assert y == z
        ans = sh

        sh = set()
        o = list(self.inverse().oneline)
        if n % 2 != 0:
            a = abs(o[0])
            sh.add((-a, a))
            o = o[1:]
        while o:
            a, b = o[:2]
            o = o[2:]
            if 0 < b < -a:
                sh.add((a, -b))
                sh.add((b, -a))
            elif b < 0 < -a:
                sh.add((b, -b))
                sh.add((a, -a))
        assert ans == sh

        return sh

    def dshape(self, offset=0, verbose=False):
        assert self.is_atom_d(offset % 2 != 0)

        o = list(self.inverse().oneline)
        init = [abs(a) for a in o[:offset]]
        ndes, fix, neg = self._ndes(o[offset:])
        
        ndes = [(abs(a), b) for (a, b) in ndes]

        if len(fix) > 0 and len(neg) % 2 != 0:
            a, b = -neg[-1], fix[0]
            if a < b:
                neg = neg[:-1]
                fix = (a,) + fix
            else:
                neg = neg + (b,)
                fix = fix[:-1]

        desd = [(b, a) for a, b in ndes if a >= -b]
        desd += [(-neg[i], neg[i + 1]) for i in range(0, len(neg) - 1, 2)]
        
        negd = init + [a for a, b in ndes if a < -b] + [-b for a, b in ndes if a < -b]
        
        sh = set()
        for a, b in ndes:
            if a < -b:
                sh.add((b, -a))
                sh.add((a, -b))
        for a in init:
            sh.add((-a, a))

        y = self.dtype_demazure_conjugate(offset % 2 != 0)
        
        z = SignedPermutation.identity(self.rank)
        if offset % 2 != 0:
            z = SignedPermutation.s_i(0, self.rank)
        
        for i in negd:
            z *= SignedPermutation.t_ij(-i, i, self.rank)
        for a, b in desd:
            z *= SignedPermutation.t_ij(a, b, self.rank)

        if verbose:
            print()
            print('*', 'y =', y, offset, 'offset', self.inverse())
            print(' ', 'z =', z)
            print()
            print(' ndes =', ndes)
            print('  fix =', fix)
            print('  neg =', neg)
            print()
            print(' desd =', desd)
            print(' negd =', negd)
            print()
            print('   sh =', sh)

        assert y == z
        return sh

    def shape(self, offset=0):
        o = list(self.inverse().oneline)
        ndes, fix, neg = self._ndes(o[offset:])
        neg += tuple(-abs(a) for a in o[:offset])

        desb = [(b, a) for a, b in ndes if not (0 < a < -b)]
        negb = [(-i, -i) for i in neg] + [(-a, -a) for a, b in ndes if 0 < a < -b] + [(b, b) for a, b in ndes if 0 < a < -b]
        n = self.rank
        y = SignedPermutation.identity(n)
        for (i, i) in negb:
            y *= SignedPermutation.t_ij(-i, i, n)
        for a, b in desb:
            y *= SignedPermutation.t_ij(a, b, n)
        assert self.inverse() % self == y

        sh = set()
        for a, b in ndes:
            if 0 < a < -b:
                sh.add((b, -a))
                sh.add((a, -b))
        for e in neg:
            sh.add((e, -e))
        return sh

    def ndes(self):
        ndes, fix, neg = self._ndes()
        return ndes

    def nfix(self):
        ndes, fix, neg = self._ndes()
        return fix

    def nneg(self):
        ndes, fix, neg = self._ndes()
        return neg

    def _ndes(self, oneline=None):
        if oneline is None:
            oneline = list(self.inverse().oneline)
        ndes = []
        while True:
            i = [i for i in range(len(oneline) - 1) if oneline[i] > oneline[i + 1]]
            if len(i) == 0:
                break
            i = i[0]
            a, b = oneline[i:i + 2]
            ndes.append((a, b))
            oneline = oneline[:i] + oneline[i + 2:]
        fix = tuple(i for i in oneline if i > 0)
        neg = tuple(i for i in oneline if i < 0)
        return tuple(sorted(ndes)), fix, neg

    @classmethod
    def ds_i(cls, i, n):
        assert 0 < abs(i) < n
        if i > 0:
            oneline = list(range(1, i)) + [i + 1, i] + list(range(i + 2, n + 1))
        else:
            i = -i
            oneline = list(range(1, i)) + [-i - 1, -i] + list(range(i + 2, n + 1))
        return SignedPermutation(*oneline)

    def dstar(self):
        return SignedPermutation(*[-self(i) if abs(self(i)) == i > 2 else self(i) for i in range(1, self.rank + 1)])

    def half_signs(self):
        k = len([i for i in range(1, self.rank + 1) if -i != self(i) < 0])
        assert k % 2 == 0
        return k // 2

    def brion_length_b(self):
        ans = 0
        n = self.rank
        y = SignedPermutation.identity(n)
        for i in self.get_reduced_word():
            assert i not in y.right_descent_set
            t = SignedPermutation.s_i(i, n)
            if i > 0 and y(i) == i and y(i + 1) == i + 1:
                ans += 1
            elif i == 0 and y(1) == 1:
                ans += 1
            y = y * t if t * y == y * t else t * y * t
        return ans

    def brion_length_c(self):
        ans = 0
        n = self.rank
        y = SignedPermutation.identity(n)
        for i in self.get_reduced_word():
            assert i not in y.right_descent_set
            t = SignedPermutation.s_i(i, n)
            if i > 0 and y(i) == i and y(i + 1) == i + 1:
                ans += 1
            y = y * t if t * y == y * t else t * y * t
        return ans

    def brion_length_d(self, twisted=False):
        ans = 0
        n = self.rank
        y = SignedPermutation.identity(n)
        u = SignedPermutation.s_i(0, n) if twisted else SignedPermutation.identity(n)
        for i in self.get_reduced_word(dtype=True):
            assert i not in y.right_descent_set
            t = SignedPermutation.ds_i(i, n)
            if i > 0 and (u * y)(i) == i and (u * y)(i + 1) == i + 1:
                ans += 1
            elif i == -1 and (u * y)(1) == 1 and (u * y)(2) == 2:
                ans += 1
            y = y * t if u * t * u * y == y * t else u * t * u * y * t
        return ans

    def _min_abs_fpf_inv_atom_oneline(self):
        return tuple(i for p in self.cyc() for i in p)

    @classmethod
    def _deduplicate_oneline(cls, tup):
        minimum = []
        for i in tup:
            if minimum and minimum[-1] == i:
                continue
            minimum += [i]
        return tuple(minimum)

    def _min_inv_atom_oneline(self, sh):
        tup = tuple(i for p in self.cyc(sh) for i in reversed(p))
        return self._deduplicate_oneline(tup)

    def _max_inv_atom_oneline(self, sh):
        cyc = sorted(self.cyc(sh), key=lambda pair: pair[1])
        tup = tuple(i for p in cyc for i in reversed(p))
        return self._deduplicate_oneline(tup)

    def get_min_abs_fpf_atom(self):
        assert self.is_abs_fpf_involution()
        return SignedPermutation(*self._min_abs_fpf_inv_atom_oneline()).inverse()

    def get_min_atom(self, sh=None):
        assert self == self.inverse()
        return SignedPermutation(*self._min_inv_atom_oneline(sh)).inverse()

    def get_max_atom(self, sh=None):
        assert self == self.inverse()
        return SignedPermutation(*self._max_inv_atom_oneline(sh)).inverse()

    @classmethod
    def get_minimal_fpf_involution(cls, n):
        if n % 2 == 0:
            oneline = [i for j in range(2, n + 1, 2) for i in [j, j - 1]]
        else:
            oneline = [-1] + [i for j in range(3, n + 1, 2) for i in [j, j - 1]]
        return SignedPermutation(*oneline)

    def get_fpf_atoms(self, offset=None):
        assert self.is_fpf_involution()
        n = self.rank
        if offset is None:
            offset = 0 if n % 2 == 0 else 1
        assert 0 <= offset <= n and (n - offset) % 2 == 0
        for w in self.get_atoms():
            oneline = w.inverse().oneline
            if all(a < 0 for a in oneline[:offset]) and all(oneline[i] > oneline[i + 1] for i in range(offset, len(oneline) - 1, 2)):
                newline = tuple(reversed([-a for a in oneline[:offset]]))
                for i in range(offset, len(oneline) - 1, 2):
                    newline += (oneline[i + 1], oneline[i])
                yield SignedPermutation(*newline).inverse()

    def get_abs_fpf_atoms(self):
        assert self.is_abs_fpf_involution()

        def upnext(oneline):
            for i in range(len(oneline) - 3):
                a, d, b, c = oneline[i:i + 4]
                if a < b < c < d:
                    newline = oneline[:i] + (b, c, a, d) + oneline[i + 4:]
                    yield newline

        minimum = self._min_abs_fpf_inv_atom_oneline()
        add = {minimum}
        while add:
            for w in add:
                yield SignedPermutation(*w).inverse()
            add = {new for w in add for new in upnext(w)}

    @classmethod
    def relative_atoms(cls, y, z):
        assert y.rank == z.rank
        if y == z:
            yield cls.identity(y.rank)
        elif y.involution_length() < z.involution_length():
            for i in range(y.rank):
                s = cls.s_i(i, y.rank)
                v = s % y % s
                if v != y:
                    for a in cls.relative_atoms(v, z):
                        yield s * a

    def get_atoms_by_shape(self):
        shapes = defaultdict(set)
        for w in self.get_atoms():
            shapes[tuple(sorted(w.shape()))].add(w)
        return shapes

    def get_involution_hecke_atoms(self, dtype=False, twisted=False):
        n = self.rank
        if twisted:
            assert dtype
            t = self.s_i(0, n)
            return [w for w in SignedPermutation.all(n, dtype=True) if (t * w * t).inverse().dtype_demazure(w) == self]
        elif dtype:
            return [w for w in SignedPermutation.all(n, dtype=True) if w.inverse().dtype_demazure(w) == self]
        else:   
            return [w for w in SignedPermutation.all(n) if w.inverse() % w == self]

    def get_atoms(self, offset=None):
        n = self.rank
        if offset is None:
            offset = 0
        assert 0 <= offset <= n

        w = self.reduce()
        if w not in atoms_b_cache:
            atoms_b_cache[w] = list(w._get_atoms())

        ans = [x.inflate(self.rank) for x in atoms_b_cache[w]]
        if offset == 0:
            return ans

        bns = []
        for w in ans:
            oneline = w.inverse().oneline
            if all(a < 0 for a in oneline[:offset]):
                newline = tuple(reversed([-a for a in oneline[:offset]])) + oneline[offset:]
                bns.append(SignedPermutation(*newline).inverse())
        return bns

    def _get_atoms(self):
        def involution(oneline):
            word = SignedPermutation(*oneline).get_reduced_word()
            n = len(oneline)
            w = self.identity(n)
            for i in reversed(word):
                s = SignedPermutation.s_i(i, n)
                if i in w.right_descent_set:
                    return None
                elif s * w == w * s:
                    w = w * s
                else:
                    w = s * w * s
            return w

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
                yield SignedPermutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def get_fpf_atoms_d(self):
        n = self.rank
        twisted = n % 2 != 0
        atoms = self.get_atoms_d(twisted)
        offset = 1 if twisted else 0
        for w in atoms:
            oneline = w.inverse().oneline
            if all(oneline[i] > oneline[i + 1] for i in range(offset, len(oneline) - 1, 2)):
                newline = oneline[:offset]
                for i in range(offset, len(oneline) - 1, 2):
                    newline += (oneline[i + 1], oneline[i])
                yield SignedPermutation(*newline).inverse()

    def get_atoms_d(self, twisted=False, offset=0):
        w = self.reduce(limit=max(2, offset))
        assert w.is_even_signed()
        assert offset <= w.rank
        key = (w, twisted, offset)
        if key not in atoms_d_cache:
            atoms_d_cache[key] = list(set(w._get_atoms_d(twisted, offset)))
        ans = atoms_d_cache[key]
        return [x.inflate(self.rank) for x in ans]

    @classmethod
    def bbase(cls, n, offset):
        return cls.longest_element(n, offset)

    @classmethod
    def dbase(cls, n, offset):
        oneline = [-i for i in range(1, offset + 1)] + [i for i in range(offset + 1, n + 1)]
        if offset % 2 != 0 and len(oneline) > 0:
            oneline[0] = 1
        return cls(*oneline)

    @classmethod
    def bbase_atom(cls, n, k):
        assert k <= n
        oneline = list(range(1, n + 1))
        for i in range(k):
            #oneline[i] = -k + i
            oneline[i] *= -1 if i % 2 != k % 2 else 1
        return SignedPermutation(*oneline) 

    @classmethod
    def bbase_fpf_atom(cls, n, k):
        assert k <= n
        oneline = list(range(1, n + 1))
        for i in range(k):
            oneline[i] = (k - i) * (-1 if i % 2 != 0 else 1)
        return SignedPermutation(*oneline) 

    @classmethod
    def dbase_atom(cls, n, twisted, k):
        assert k <= n
        # oneline = list(range(1, n + 1))
        # if k % 2 == 0:
        #     for i in range(0, k, 2):
        #         oneline[i] = 1 - k + i
        #         oneline[i + 1] = -k + i
        # else:
        #     for i in range(0, k - 1, 2):
        #         oneline[i] = 1 - k + i
        #         oneline[i + 1] = -k + i
        #     oneline[k - 1] = 1 
        # return SignedPermutation(*oneline) 
        oneline = list(range(1, n + 1))
        if k % 2 == 0:
            for i in range(1, k, 2):
                oneline[i] *= -1
            if (k // 2) % 2 != 0:
                oneline[0] *= -1
        else:
            for i in range(2, k, 2):
                oneline[i] *= -1
            if (k // 2) % 2 != 0:
                oneline[0] *= -1  
        w = cls(*oneline)

        if twisted and k > 0 and k % 2 == 0:
            w *= cls.ds_i(-1 if k % 4 == 0 else 1, n)

        z = cls.dbase(n, k)
        assert w in z.get_atoms_d(twisted=twisted)

        return w

    def _get_atoms_d(self, twisted, offset):
        n = self.rank
        t = self.s_i(0, n) if twisted else self.identity(n)
        length = lambda w: w.dlength()

        base = self.dbase(n, offset)

        if length(self) <= length(base):
            if self == base:
                yield SignedPermutation.identity(n)
            return

        for i in range(n):
            s = self.ds_i(i if i != 0 else -1, n)
            w = self * s
            if length(w) < length(self):
                if w == t * s * t * self:
                    for a in w.get_atoms_d(twisted, offset):
                        yield a * s
                else:
                    for a in (t * s * t * w).get_atoms_d(twisted, offset):
                        yield a * s
