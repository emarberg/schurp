from permutations import Permutation
from signed import SignedPermutation, SignedMixin
from collections import defaultdict
from vectors import Vector
from symmetric import SchurP, SchurQ
from partitions import StrictPartition
import itertools
import operator


D_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_COUNTS = {(): 1}

EVEN_SIGNED_INVOLUTION_WORDS = {}
EVEN_SIGNED_FPF_INVOLUTION_WORDS = {}

ATOMS_D_CACHE = {}
SCHURD_STANSYM_CACHE = {}


class EvenSignedPermutation(SignedMixin):

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
        self._twisted_involution_length = None

    def __iter__(self):
        return self.oneline.__iter__()
        
    def __repr__(self):
        # return 'EvenSignedPermutation(' + ', '.join([repr(i) for i in self.oneline]) + ')'
        return str(self)

    @classmethod
    def reflections(cls, n):
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                yield cls.reflection_t(i, j, n)
                yield cls.reflection_s(i, j, n)

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
    def twisted_involutions(cls, n):
        return cls.involutions(n, True)

    @classmethod
    def involutions(cls, n, twist=False):
        level = {EvenSignedPermutation.identity(n)}
        while level:
            nextlevel = set()
            for w in level:
                yield w
                for i in range(n):
                    if i not in w.right_descent_set:
                        s = EvenSignedPermutation.s_i(i, n)
                        t = EvenSignedPermutation.s_i(1 - i, n) if i in [0, 1] and twist else s
                        nextlevel.add(t % w % s)
            level = nextlevel

    @classmethod
    def fpf_involutions(cls, n):
        for w in (cls.involutions(n) if n % 2 == 0 else cls.twisted_involutions(n)):
            if w.is_fpf_involution():
                yield w

    def is_involution(self, twist=False):
        if twist:
            return self.star() == self.inverse()
        else:
            return self == self.inverse()

    def is_twisted_involution(self):
        return self.is_involution(True)

    def is_abs_fpf_involution(self):
        n = self.rank
        return all(abs(self(i)) != i for i in range(1, n + 1))

    def is_fpf_involution(self):
        n = self.rank
        if n % 2 == 0:
            return self.is_involution() and all(self(i) != i for i in range(1, n + 1))
        if n % 2 != 0:
            w = SignedPermutation.s_i(0, n) * SignedPermutation(*self.oneline)
            return w.is_fpf_involution()

    def get_twisted_involution_word(self):
        return self.get_involution_word(True)

    def get_involution_word(self, twist=False):
        if self.right_descent_set:
            i = min(self.right_descent_set)
            s = EvenSignedPermutation.s_i(i, self.rank)
            t = s.star() if twist else s
            if t * self == self * s:
                return (self * s).get_involution_word(twist) + (i,)
            else:
                return (t * self * s).get_involution_word(twist) + (i,)
        else:
            return ()

    def star(self):
        oneline = [i if abs(i) != 1 else -i for i in self.oneline]
        if oneline:
            oneline[0] *= -1
        return EvenSignedPermutation(*oneline)

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
        des = lambda w: w.right_descent_set
        return self._get_reduced_words(des, EVEN_SIGNED_REDUCED_WORDS)

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

    def get_involution_words(self, twist=False):
        assert self.is_involution(twist)
        w = self.reduce()
        oneline = w.oneline
        key = (oneline, twist)
        if key not in EVEN_SIGNED_INVOLUTION_WORDS:
            words = []
            for a in w.get_atoms(twist):
                for word in a.get_reduced_words():
                    words.append(word)
            EVEN_SIGNED_INVOLUTION_WORDS[key] = sorted(words, key=lambda x: self.flatten(x))
        return EVEN_SIGNED_INVOLUTION_WORDS[key]

    def get_twisted_involution_words(self):
        return self.get_involution_words(True)

    def get_fpf_involution_words(self):
        assert self.is_fpf_involution()
        key = self.oneline
        if key not in EVEN_SIGNED_FPF_INVOLUTION_WORDS:
            words = []
            for a in self.get_fpf_atoms():
                for word in a.get_reduced_words():
                    words.append(word)
            EVEN_SIGNED_FPF_INVOLUTION_WORDS[key] = sorted(words, key=lambda x: self.flatten(x))
        return EVEN_SIGNED_FPF_INVOLUTION_WORDS[key]

    def fpf_shape(self):
        raise NotImplementedError

    def involution_fixed_points(self, twist=False):
        y = self
        n = y.rank
        yfixed = {i for i in range(-n, n + 1) if i not in [-1, 0, 1] and y(i) == i}
        if not twist and y(1) == 1:
            yfixed |= {-1, 1}
        if twist and y(1) == -1:
            yfixed |= {-1, 1}
        return yfixed

    def twisted_shape(self, verbose=False):
        def vprint(*args):
            if verbose:
                print(*args) # noqa

        n = self.rank
        y = self.inverse().star() % self
        assert y.twisted_involution_length() == self.length()

        twist = n % 2 == 0
        w0 = EvenSignedPermutation.longest_element(n)
        y = w0 * y
        aword = list(reversed(self.get_reduced_word()))
        vprint('word:', aword)

        yfixed = y.involution_fixed_points(twist)
        v = EvenSignedPermutation.identity(n)
        sh = set()
        for a in aword:
            vprint('  v =', v, 'y =', y.cycle_repr(), 'shape =', sh, 'a =', a, 'leads to')

            if a > 0 and {a, a + 1}.issubset(y.involution_fixed_points(twist)):
                e, f = tuple(sorted([abs(v(a)), abs(v(a + 1))]))
                sh |= {(e, f), (-f, -e)}
            elif a == 0 and {1, 2}.issubset(y.involution_fixed_points(twist)):
                e, f = tuple(sorted([abs(v(1)), abs(v(2))]))
                sh |= {(e, f), (-f, -e)}
            s = EvenSignedPermutation.s_i(a, n)
            t = s.star() if twist else s

            u = v * s
            assert u.length() == v.length() + 1
            v = u

            z = t % y % s
            assert z.involution_length(twist) == y.involution_length(twist) + 1
            y = z
        vprint('  v =', v, 'y =', y.cycle_repr(), 'shape =', sh)
        vprint()
        f = {i for p in sh for i in p}
        return sh | {(-i, i) for i in yfixed - f if i > 0}

    def shape(self, verbose=False):
        def vprint(*args):
            if verbose:
                print(*args) # noqa

        n = self.rank
        y = self.inverse() % self
        assert y.involution_length() == self.length()

        twist = n % 2 != 0
        w0 = EvenSignedPermutation.longest_element(n)
        y = w0 * y
        aword = list(reversed(self.get_reduced_word()))
        vprint('word:', aword)

        yfixed = y.involution_fixed_points(twist)
        v = EvenSignedPermutation.identity(n)
        sh = set()
        for a in aword:
            vprint('  v =', v, 'y =', y.cycle_repr(), 'shape =', sh, 'a =', a, 'leads to')

            if a > 0 and {a, a + 1}.issubset(y.involution_fixed_points(twist)):
                e, f = tuple(sorted([abs(v(a)), abs(v(a + 1))]))
                sh |= {(e, f), (-f, -e)}
            elif a == 0 and {1, 2}.issubset(y.involution_fixed_points(twist)):
                e, f = tuple(sorted([abs(v(1)), abs(v(2))]))
                sh |= {(e, f), (-f, -e)}
            s = EvenSignedPermutation.s_i(a, n)
            t = s.star() if twist else s

            u = v * s
            assert u.length() == v.length() + 1
            v = u

            z = t % y % s
            assert z.involution_length(twist) == y.involution_length(twist) + 1
            y = z
        vprint('  v =', v, 'y =', y.cycle_repr(), 'shape =', sh)
        vprint()
        f = {i for p in sh for i in p}
        return sh | {(-i, i) for i in yfixed - f if i > 0}

    def twisted_pair(self):
        def y(i):
            a = self(i)
            if a in [1, -1]:
                return -a
            return a

        n = self.rank
        return sorted([
            (a, y(a))
            for a in range(-n, n + 1)
            if 0 < abs(a) < y(a)
        ])

    def twisted_fixed_points(self):
        n = self.rank
        ans = []
        for a in range(-n, n + 1):
            if a == -1 and self(a) == 1:
                ans.append(a)
            elif a == 1 and self(a) == -1:
                ans.append(a)
            elif abs(a) > 1 and self(a) == a:
                ans.append(a)
        return tuple(ans)

    def twisted_negated_points(self):
        n = self.rank
        ans = []
        for a in range(-n, n + 1):
            if a == -1 and self(a) == -1:
                ans.append(a)
            elif a == 1 and self(a) == 1:
                ans.append(a)
            elif abs(a) > 1 and self(a) == -a:
                ans.append(a)
        return tuple(ans)

    def twisted_nfix(self):
        ndes, fix = self._twisted_ndes()
        return fix

    def twisted_ndes(self):
        ndes, fix = self._twisted_ndes()
        return ndes

    def _twisted_ndes(self):
        y = self.inverse().star() % self
        assert y.twisted_involution_length() == self.length()

        oneline = tuple(self.inverse().oneline)
        ndes = []

        if oneline:
            a = abs(oneline[0])
            ndes.append((a, -a))
            oneline = oneline[1:]

        while True:
            if oneline and oneline[0] < 0:
                oneline = (-oneline[0],) + oneline[1:]
            i = [i for i in range(len(oneline) - 1) if oneline[i] > oneline[i + 1]]
            if len(i) == 0:
                break
            i = i[0]
            a, b = oneline[i:i + 2]
            ndes.append((abs(a), b))
            oneline = oneline[:i] + oneline[i + 2:]
        return tuple(sorted(ndes)), oneline

    def nfix(self):
        ndes, fix = self._ndes()
        return fix

    def ndes(self):
        ndes, fix = self._ndes()
        return ndes

    def _ndes(self):
        y = self.inverse() % self
        assert y.involution_length() == self.length()

        oneline = tuple(self.inverse().oneline)
        ndes = []
        while True:
            if oneline and oneline[0] < 0:
                oneline = (-oneline[0],) + oneline[1:]
            i = [i for i in range(len(oneline) - 1) if oneline[i] > oneline[i + 1]]
            if len(i) == 0:
                break
            i = i[0]
            a, b = oneline[i:i + 2]
            ndes.append((abs(a), b))
            oneline = oneline[:i] + oneline[i + 2:]
        return tuple(sorted(ndes)), oneline

    def twisted_involution_length(self):
        return self.involution_length(True)

    def involution_length(self, twist=False):
        if twist:
            if self._twisted_involution_length is None:
                self._twisted_involution_length = len(self.get_twisted_involution_word())
            return self._twisted_involution_length
        else:
            if self._involution_length is None:
                self._involution_length = len(self.get_involution_word())
            return self._involution_length

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

    def half_signs(self):
        k = len([i for i in range(1, self.rank + 1) if -i != self(i) < 0])
        assert k % 2 == 0
        return k // 2

    def length(self):
        return len(self)

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

    @classmethod
    def get_minimal_fpf_involution(cls, n):
        if n % 2 == 0:
            oneline = [i for j in range(2, n + 1, 2) for i in [j, j - 1]]
        else:
            oneline = [1] + [i for j in range(3, n + 1, 2) for i in [j, j - 1]]
        return cls(*oneline)

    def get_fpf_atoms(self):
        y = self.get_minimal_fpf_involution(self.rank)
        twist = self.rank % 2 != 0
        return self.relative_atoms(y, self, twist)

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

    def get_atoms_by_shape(self):
        shapes = defaultdict(set)
        for w in self.get_atoms():
            shapes[tuple(sorted(w.shape()))].add(w)
        return shapes

    @classmethod
    def relative_twisted_atoms(cls, y, z):
        for w in cls.relative_atoms(y, z, True):
            yield w

    @classmethod
    def relative_atoms(cls, y, z, twist=False):
        assert y.is_involution(twist)
        assert z.is_involution(twist)

        if y == z:
            yield EvenSignedPermutation.identity(y.rank)
        elif y.involution_length(twist) < z.involution_length(twist):
            for i in range(y.rank):
                s = EvenSignedPermutation.s_i(i, y.rank)
                t = s.star() if twist else s
                v = t % y % s
                if v != y:
                    for a in cls.relative_atoms(v, z, twist):
                        yield s * a

    def get_twisted_atoms(self):
        return self.get_atoms(True)

    def get_atoms(self, twist=False):
        assert self.is_involution(twist)
        w = self.reduce()
        key = (w, twist)
        if key not in ATOMS_D_CACHE:
            ATOMS_D_CACHE[key] = list(set(w._get_atoms(twist)))
        ans = ATOMS_D_CACHE[key]
        return [x.inflate(self.rank) for x in ans]

    def _get_atoms(self, twist=False):
        if len(self) == 0:
            yield self
            return

        for i in range(self.rank):
            s = EvenSignedPermutation.s_i(i, self.rank)
            t = s.star() if twist else s
            w = self * s
            if len(w) < len(self):
                if w == t * self:
                    for a in w.get_atoms(twist):
                        yield a * s
                else:
                    for a in (t * w).get_atoms(twist):
                        yield a * s
    
    def get_min_fpf_atom(self, matching=None):
        raise NotImplementedError

    def get_min_atom(self, matching=None):
        assert self.is_involution()
        if matching is None:
            g = sorted(self.negated_points())
            matching = {(g[i], g[i + 1]) for i in range(0, len(g), 2)}
        x = {i for m in matching for i in m if i > 0}

        init = {i for i in x if (-i, i) in matching}
        init = sorted(init)
        k = len(init)
        for i in range(1 if k % 2 == 0 else 1, k, 2):
            init[i] *= -1

        matching = [(a, b) for a, b in matching if a != -b]
        fix = [(i,) for i in self.fixed_points() if i > 0]
        neg = [(-i,) for i in self.negated_points() if i not in x and i > 0]
        pair = [(b, a) for a, b in self.pair()]
        des = [(a, -b) for a, b in matching if a > 0]
        oneline = init + [
            i
            for m in sorted(fix + neg + pair + des, key=operator.itemgetter(0))
            for i in m
        ]
        if len([i for i in oneline if i < 0]) % 2 != 0:
            oneline[0] *= -1
        w = EvenSignedPermutation(*oneline).inverse()
        assert w.inverse() % w == self
        assert self.involution_length() == w.length()
        return w

    def get_min_twisted_atom(self, matching=None):
        assert self.is_twisted_involution()
        if matching is None:
            g = sorted(self.twisted_negated_points())
            matching = {(g[i], g[i + 1]) for i in range(0, len(g) - 1, 2)}
        x = {i for m in matching for i in m if i > 0}

        init = {i for i in x if (-i, i) in matching}
        init = sorted(init)
        k = len(init)
        for i in range(1 if k % 2 == 0 else 1, k, 2):
            init[i] *= -1

        matching = [(a, b) for a, b in matching if a != -b]
        fix = [(i,) for i in self.twisted_fixed_points() if i > 0]
        neg = [(-i,) for i in self.twisted_negated_points() if i not in x and i > 0]
        pair = [(b, a) for a, b in self.twisted_pair()]
        des = [(a, -b) for a, b in matching if a > 0]
        oneline = init + [
            i
            for m in sorted(fix + neg + pair + des, key=operator.itemgetter(0))
            for i in m
        ]
        if len([i for i in oneline if i < 0]) % 2 != 0:
            a, b = oneline[:2]
            if abs(a) < abs(b):
                oneline[0] *= -1
            else:
                oneline[1] *= -1
        w = EvenSignedPermutation(*oneline).inverse()
        assert w.inverse().star() % w == self
        assert self.twisted_involution_length() == w.length()
        return w
