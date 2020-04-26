import itertools
import subprocess
from collections import defaultdict
from vectors import Vector
from symmetric import SchurP, SchurQ
from partitions import StrictPartition
from permutations import Permutation


C_SIGNED_REDUCED_WORDS = {(): {()}}
B_SIGNED_REDUCED_WORDS = {(): {()}}
D_SIGNED_REDUCED_WORDS = {(): [()]}


SIGNED_REDUCED_WORDS = {(): {()}}
EVEN_SIGNED_REDUCED_WORDS = {(): [()]}
EVEN_SIGNED_REDUCED_COUNTS = {(): 1}
EVEN_SIGNED_INVOLUTION_WORDS = {(): [()]}

schurp_stansym_cache = {}
schurq_stansym_cache = {}
schurd_stansym_cache = {}
atoms_b_cache = {}
atoms_d_cache = {}


class SignedPermutation:

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

    def shape(self):
        ndes, fix, neg = self._ndes()

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

    def _ndes(self):
        oneline = tuple(self.inverse().oneline)
        ndes = []
        while True:
            i = [i for i in range(len(oneline) - 1) if oneline[i] > oneline[i + 1]]
            if len(i) == 0:
                break
            i = i[0]
            ndes.append(oneline[i:i + 2])
            oneline = oneline[:i] + oneline[i + 2:]
        fix = [i for i in oneline if i > 0]
        neg = [i for i in oneline if i < 0]
        return ndes, fix, neg

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

    def __init__(self, *oneline):
        self.oneline = tuple(oneline)
        self.rank = len(oneline)
        assert set(range(1, self.rank + 1)) == set(abs(i) for i in self.oneline)
        # cached fields
        self._rdes = None
        self._ldes = None
        self._len = None

    @classmethod
    def from_word(cls, n, *args):
        if len(args) == 1 and type(args[0]) == tuple:
            args = args[0]
        w = SignedPermutation.identity(n)
        for i in args:
            w *= SignedPermutation.s_i(i, n)
        return w

    def __neg__(self):
        return self * SignedPermutation.longest_element(self.rank)

    def __repr__(self):
        # return 'SignedPermutation(' + ', '.join([repr(i) for i in self.oneline]) + ')'
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
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            for v in range(2**n):
                oneline = []
                for i in range(n):
                    oneline.append(args[i] * (-1) ** (v % 2))
                    v = v // 2
                yield SignedPermutation(*oneline)

    @classmethod
    def permutations(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield SignedPermutation(*args)

    def is_abs_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(abs(self(i)) != i for i in range(1, n + 1))

    def is_fpf_involution(self):
        n = self.rank
        return self.is_involution() and all(self(i) != i for i in range(1, n + 1))

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
                yield SignedPermutation(*newline)

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
                    s = SignedPermutation.s_i(i, n)
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

    def get_reduced_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = SignedPermutation.s_i(i, self.rank)
            return (i,) + (s * self).get_reduced_word()
        else:
            return ()

    def get_involution_word(self):
        assert self.is_involution()
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = SignedPermutation.s_i(i, self.rank)
            w = s * self
            if w != self * s:
                w = w * s
            return (i,) + w.get_involution_word()
        else:
            return ()

    def get_reduced_words(self):
        w = self.reduce()
        oneline = w.oneline
        if oneline not in SIGNED_REDUCED_WORDS:
            words = set()
            for i in w.right_descent_set:
                s = SignedPermutation.s_i(i, w.rank)
                words |= {e + (i,) for e in (w * s).get_reduced_words()}
            SIGNED_REDUCED_WORDS[oneline] = words
        return SIGNED_REDUCED_WORDS[oneline]

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
        w = self
        n = w.rank
        assert w.is_fpf_involution()

        if len(w) == (n + 1) // 2:
            yield ()
        else:
            for i in w.right_descent_set:
                if i == 0 and w(1) == -1:
                    continue
                if i > 0 and w(i) == i + 1:
                    continue
                s = SignedPermutation.s_i(i, w.rank)
                if s * w != w * s:
                    v = s * w * s
                else:
                    v = w * s
                for word in v.get_fpf_involution_words():
                    yield word + (i,)

    def get_involution_words(self):
        w = self.reduce()
        assert w.inverse() == w
        for a in w.get_atoms():
            for word in a.get_reduced_words():
                yield word

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
        return SignedPermutation(*newline)

    def involution_length(self):
        return (len(self.neg()) + len(self.pair()) + len(self)) // 2

    @classmethod
    def identity(cls, n):
        return SignedPermutation(*list(range(1, n + 1)))

    @classmethod
    def longest_element(cls, n):
        return SignedPermutation(*[-i for i in range(1, n + 1)])

    @classmethod
    def s_i(cls, i, n):
        assert 0 <= i < n
        if i == 0:
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

    def __mod__(self, other):
        assert type(other) == SignedPermutation
        assert self.rank == other.rank
        if other.left_descent_set:
            i = next(iter(other.left_descent_set))
            s = SignedPermutation.s_i(i, self.rank)
            if i in self.right_descent_set:
                return self % (s * other)
            else:
                return (self * s) % (s * other)
        else:
            return self

    def __mul__(self, other):
        assert type(other) == SignedPermutation
        assert self.rank == other.rank
        newline = [self(other(i)) for i in range(1, self.rank + 1)]
        return SignedPermutation(*newline)

    def inverse(self):
        newline = self.rank * [0]
        for i in range(1, self.rank + 1):
            j = self(i)
            if j > 0:
                newline[j - 1] = i
            else:
                newline[-j - 1] = -i
        return SignedPermutation(*newline)

    def __lt__(self, other):
        assert type(other) == SignedPermutation
        return self.oneline < other.oneline

    def __eq__(self, other):
        assert type(other) == SignedPermutation
        return self.oneline == other.oneline

    def __pow__(self, n):
        if n < 0:
            return self.inverse().__pow__(-n)
        elif n == 0:
            return SignedPermutation.identity(self.rank)
        elif n == 1:
            return SignedPermutation(*self.oneline)
        else:
            p = n // 2
            q = n - p
            return self.__pow__(p) * self.__pow__(q)

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

    def is_even_signed(self):
        return len([i for i in self.oneline if i < 0]) % 2 == 0

    def inv_stanley_schur_d_decomposition(self):
        assert self.is_even_signed()
        ans = Vector()
        for x in self.get_atoms_d():
            ans += x.stanley_schur_d_decomposition()
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
    def queue_stanley_decomposition(cls, n, verbose=True):
        def get_shape(oneline):
            while oneline[-1] == len(oneline):
                oneline = oneline[:-1]
            oneline = oneline[1:]
            ans = []
            while oneline:
                for i in range(len(oneline)):
                    a = oneline[i]
                    if i == 0 and a < 0:
                        ans += [(a, -a)]
                        oneline = oneline[1:]
                        break
                    if i + 1 >= len(oneline):
                        continue
                    b = oneline[i + 1]
                    if 0 < a < -b:
                        ans += [(a, -b)]
                        oneline = oneline[:i] + oneline[i + 2:]
                        break
            return ans

        def get_rs(x):
            r = 0
            s = n
            while r < s - 1 and (x(r + 1) == n or x(r + 1) == x(r) - 1):
                r += 1
            try:
                z = SignedPermutation(*x.oneline[r:]).reduce()
                a = SignedPermutation.longest_element(n - r - 1).get_atoms()
                assert z in a
            except:
                r += 1
                shape = get_shape(tuple(x(i) for i in range(r, n + 1)))
                s = x.inverse()(
                    max([a for a, b in shape if 0 < a < x(r) < b])
                )
            return r, s

        ans = []
        atoms = SignedPermutation.longest_element(n).get_atoms()
        preatoms = SignedPermutation.longest_element(n - 1).get_atoms()
        queue = [SignedPermutation(*((n + 1,) + x.oneline + (n,))) for x in preatoms]
        while queue:
            print('queue: %s' % len(queue))
            x = queue.pop(0).reduce()
            if x in atoms:
                print('\n*', x, 'is atom\n')
                ans += [x]
                continue

            n = x.rank
            r, s = get_rs(x)
            v = x * cls.reflection_t(r, s, n)
            v_len = len(v)
            assert v_len + 1 == len(x)

            newline = v.oneline + (n + 1,)
            new_v = SignedPermutation(*newline)
            lhs = [new_v * cls.reflection_t(r, i, n + 1) for i in range(r + 1, n + 2) if i != s]
            lhs = sorted([u.reduce() for u in lhs if v_len + 1 == len(u)])

            if lhs != sorted([u.reduce() for u in queue if u.reduce() in lhs]):
                queue += [x]
                continue

            yield (x.reduce(), v.reduce())
            for u in lhs:
                yield (u, v.reduce())

            queue = [u for u in queue if u.reduce() not in lhs]
            rhs = [v * cls.reflection_s(r, r, n)]
            rhs += [v * cls.reflection_t(i, r, n) for i in range(1, r)]
            rhs += [new_v * cls.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
            rhs = [u.reduce() for u in rhs if v_len + 1 == len(u)]
            queue += rhs

            for u in rhs:
                yield(v.reduce(), u)

        assert sorted(ans) == sorted(atoms)

    @classmethod
    def standardize(cls, oneline):
        distinct = sorted([abs(j) for j in oneline])
        assert len(set(distinct)) == len(oneline)
        mapping = {v: i + 1 for i, v in enumerate(distinct)}
        newline = tuple(mapping[v] if v > 0 else -mapping[-v] for v in oneline)
        return SignedPermutation(*newline)

    def conjectural_stanley_schur_decomposition(self, starting_rank=None, verbose=True, step=0):
        w = self.reduce()
        n = w.rank

        if starting_rank is None:
            starting_rank = n

        space = ((4 + starting_rank) * step) * ' '
        if w in SignedPermutation.longest_element(starting_rank - 1).get_atoms():
            print(space, '*', w, 'is atom')
            print()
            verbose = False
            yield (w, 1)
            return

        sh = w.increasing_shape()
        if sh is not None:
            return

        def get_shape(oneline):
            while oneline[-1] == len(oneline):
                oneline = oneline[:-1]
            oneline = oneline[1:]
            ans = []
            while oneline:
                for i in range(len(oneline)):
                    a = oneline[i]
                    if i == 0 and a < 0:
                        ans += [(a, -a)]
                        oneline = oneline[1:]
                        break
                    if i + 1 >= len(oneline):
                        continue
                    b = oneline[i + 1]
                    if 0 < a < -b:
                        ans += [(a, -b)]
                        oneline = oneline[:i] + oneline[i + 2:]
                        break
            return ans

        if verbose:
            r = 0
            s = n
            while r < s - 1 and (w(r + 1) == n or w(r + 1) == w(r) - 1):
                r += 1
            try:
                z = SignedPermutation(*w.oneline[r:]).reduce()
                a = SignedPermutation.longest_element(n - r - 1).get_atoms()
                assert z in a
            except:
                try:
                    r += 1
                    shape = get_shape(tuple(w(i) for i in range(r, n + 1)))
                    s = w.inverse()(
                        max([a for a, b in shape if 0 < a < w(r) < b])
                    )
                    print('(r,s) = (%s,%s)' % (r, s))
                except:
                    print(space, 'no inversions:', w, '\n')
                    w, r, s = eval(input('w, r, s = '))
                    print('(r,s) = (%s,%s)' % (r, s))
        else:
            r = w.last_descent()
            s = w.last_inversion(r)

        v = w * self.reflection_t(r, s, n)

        v_len = len(v)
        assert v_len + 1 == len(w)

        indices = [v * self.reflection_s(r, r, n)]
        indices += [v * self.reflection_t(i, r, n) for i in range(1, r)]
        newline = v.oneline + (n + 1,)
        v = SignedPermutation(*newline)
        indices += [v * self.reflection_s(i, r, n + 1) for i in range(1, n + 2) if i != r]
        indices = [x.reduce() for x in indices if v_len + 1 == len(x)]

        subindices = [v * self.reflection_t(r, i, n + 1) for i in range(r + 1, n + 2) if i != s]
        subindices = [x.reduce() for x in subindices if v_len + 1 == len(x)]

        if verbose:
            print(space, self, '->', v.reduce(), '->', indices, ('- ' + str(subindices)) if subindices else '')
            print()

        for x in indices:
            for a, coeff in x.conjectural_stanley_schur_decomposition(starting_rank, verbose, step + 1):
                yield (a, coeff)

        for x in subindices:
            yield (x, -1)

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
        for k, v in sorted(ans.items(), key=lambda x : (len(x[0]), str(x[0]))):
            if v:
                print('\n' + k)
                for w in v:
                    print('  ', w.reduce(), w.get_min_atom().reduce())
        return ans

    def is_involution(self):
        return self == self.inverse()

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
        return SignedPermutation(*newline)

    def ell_zero(self):
        return len([i for i in range(1, self.rank + 1) if self(i) < 0])

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

    def fixed_points(self):
        n = self.rank
        return [a for a in range(-n, n + 1) if 0 != a == self(a)]

    def negated_points(self):
        n = self.rank
        return [a for a in range(-n, n + 1) if 0 != a == -self(a)]

    def fix(self):
        n = self.rank
        return [(a, a) for a in range(1, n + 1) if self(a) == a]

    def cyc(self):
        return sorted(self.pair() + self.neg() + self.fix())

    def _min_fpf_inv_atom_oneline(self):
        return tuple(i for p in self.cyc() for i in p)

    def _min_inv_atom_oneline(self):
        tup = tuple(i for p in self.cyc() for i in reversed(p))
        minimum = []
        for i in tup:
            if minimum and minimum[-1] == i:
                continue
            minimum += [i]
        return tuple(minimum)

    def get_min_fpf_atom(self):
        assert self.is_abs_fpf_involution()
        return SignedPermutation(*self._min_fpf_inv_atom_oneline()).inverse()

    def get_min_atom(self):
        assert self == self.inverse()
        return SignedPermutation(*self._min_inv_atom_oneline()).inverse()

    @classmethod
    def get_minimal_fpf_involution(cls, n):
        assert n % 2 == 0
        oneline = [i for j in range(2, n + 1, 2) for i in [j, j - 1]]
        return SignedPermutation(*oneline)

    def get_fpf_atoms(self):
        assert self.is_abs_fpf_involution()

        def upnext(oneline):
            for i in range(len(oneline) - 3):
                a, d, b, c = oneline[i:i + 4]
                if a < b < c < d:
                    newline = oneline[:i] + (b, c, a, d) + oneline[i + 4:]
                    yield newline

        minimum = self._min_fpf_inv_atom_oneline()
        add = {minimum}
        while add:
            for w in add:
                yield SignedPermutation(*w).inverse()
            add = {new for w in add for new in upnext(w)}

    def get_atoms_by_shape(self):
        shapes = defaultdict(set)
        for w in self.get_atoms():
            shapes[tuple(sorted(w.shape()))].add(w)
        return shapes

    def get_atoms(self):
        assert self == self.inverse()
        w = self.reduce()
        if w not in atoms_b_cache:
            atoms_b_cache[w] = list(w._get_atoms())
        ans = atoms_b_cache[w]
        return [x.inflate(self.rank) for x in ans]

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

    def get_atoms_d(self):
        assert self.is_even_signed()
        assert self == self.inverse()
        w = self.reduce()
        if w not in atoms_d_cache:
            atoms_d_cache[w] = list(w._get_atoms_d())
        ans = atoms_d_cache[w]
        return [x.inflate(self.rank) for x in ans]

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

        def s_i(i, n):
            return self.s_i(i, n) if i != 0 else self.s_i(0, n) * self.s_i(1, n) * self.s_i(0, n)

        for i in range(self.rank):
            s = s_i(i, self.rank)
            w = self * s
            if length(w) < length(self):
                if w == s * self:
                    for a in w.get_atoms_d():
                        yield a * s
                else:
                    for a in (s * w).get_atoms_d():
                        yield a * s

    @classmethod
    def double_involution_word(cls, word):
        if len(word) == 0:
            return tuple()
        n = 1 + max([i for i in word])
        w = SignedPermutation.identity(n)
        ans = []
        for i in word:
            assert i not in w.right_descent_set
            s = SignedPermutation.s_i(i, n)
            w *= s
            if i not in w.left_descent_set:
                ans = [i] + ans + [i]
                w = s * w
            else:
                ans = ans + [i]
        return tuple(ans)


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
        cache = schurd_stansym_cache
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
        if w not in atoms_d_cache:
            atoms_d_cache[w] = list(w._get_atoms())
            print('atoms_d_cache:', len(atoms_d_cache))
        ans = atoms_d_cache[w]
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


class SignedAtomsGraph:

    DIRECTORY = '/Users/emarberg/Dropbox/stanley-type-c/atoms/'

    def tikz(self):
        pre = [
            '\\begin{center}',
            '\\begin{tikzpicture}',
        ]

        vertices = sorted({u for u, _ in self.edges} | {v for _, v in self.edges})
        labels = {v: str(self.n) + 'node' + str(i) for i, v in enumerate(vertices)}

        nodes = [
            '\\node (%s) {%s};' % (labels[v], v.tex()) for v in labels
        ]

        edges = [
            '\\draw [->] (%s) -- (%s);' % (labels[u], labels[v]) for u, v in self.edges
        ]

        post = [
            '\\end{tikzpicture}',
            '\\end{center}',
        ]

        return '\n'.join(pre + nodes + edges + post)

    @property
    def _filename(self):
        return 'atoms_graph_b%s' % self.n

    @property
    def dot_filename(self):
        return self.DIRECTORY + 'dot/' + '%s.dot' % self._filename

    @property
    def png_filename(self):
        return self.DIRECTORY + 'png/' + '%s.png' % self._filename

    def node_label(self, i):
        if i.rank != self.n + 1:
            assert SignedPermutation.standardize(i.oneline[1:]) in self.sub_atoms
        return str(i)
        # if i.rank == self.n + 1:
        #     return str(i.inverse())

        # sh = [(a, b) for (a, b) in i.get_atom_shape(i.oneline) if a + b != 0]
        # base = ' '.join([str(a) for a in range(1, self.n + 1)])

        # def index(j):
        #     return (j - 1) * 2

        # lines = []
        # for a, b in sorted(sh, key=lambda x: x[1] - x[0]):
        #     j, k = index(a), index(b)
        #     if not (lines and all(s == ' ' for s in lines[0][j:k])):
        #         lines = [(2 * self.n) * [' ']] + lines
        #     lines[0][j] = '\u256D'
        #     lines[0][k] = '\u256E'
        #     for l in range(j + 1, k):
        #         lines[0][l] = '\u2500'
        #     for l in range(1, len(lines)):
        #         lines[l][j] = '\u2502'
        #         lines[l][k] = '\u2502'
        # base = '\n'.join([''.join(l) for l in lines] + [base])

        # return str(i.inverse()) + '\n\n' + base

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"];']
        s += ['    "%s" -> "%s";' % (self.node_label(x), self.node_label(y)) for (x, y) in self.edges]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    def __init__(self, n):
        self.n = n
        self.atoms = SignedPermutation.longest_element(n).get_atoms()
        self.sub_atoms = SignedPermutation.longest_element(n - 1).get_atoms()
        self._edges = None

    @property
    def edges(self):
        if self._edges is None:
            self._edges = [
                (a, b)
                for (a, b) in SignedPermutation.queue_stanley_decomposition(self.n)
                if a.rank == self.n
            ]
        return self._edges

    def test(self):
        """Tests for cgraphs.tex"""
        def a_shape(w):
            oneline = [w(i) for i in range(1, w.rank + 1)]
            m = set()
            go = True
            while go and oneline:
                for i in range(len(oneline)):
                    go = False
                    if i + 1 < len(oneline) and oneline[i] > oneline[i + 1]:
                        a = oneline[i]
                        b = -oneline[i + 1]
                        assert 0 < a < b
                        m |= {(a, b), (-b, -a)}
                        oneline = oneline[:i] + oneline[i + 2:]
                        go = True
                        break
            assert all(c < 0 for c in oneline)
            m |= {(c, -c) for c in oneline}
            return m

        def q_shape(w):
            oneline = [w(i) for i in range(2, w.rank + 1)]
            m = set()
            go = True
            while go and oneline:
                for i in range(len(oneline)):
                    go = False
                    if i + 1 < len(oneline) and oneline[i] > oneline[i + 1]:
                        a = oneline[i]
                        b = -oneline[i + 1]
                        assert 0 < a < b
                        m |= {(a, b), (-b, -a)}
                        oneline = oneline[:i] + oneline[i + 2:]
                        go = True
                        break
            assert all(c < 0 for c in oneline)
            m |= {(c, -c) for c in oneline}
            return m

        def is_even(w):
            return w(1) < 0 or (w(1) - self.n - 1) % 2 == 0

        def up(w, i):
            if i + 2 <= w.rank:
                c, a, b = w(i), w(i + 1), w(i + 2)
                if a < b < c:
                    s = SignedPermutation.s_i(i, w.rank)
                    t = SignedPermutation.s_i(i + 1, w.rank)
                    return w * t * s

        def is_maximal(w):
            return all(up(w, i) is None for i in range(1, w.rank))

        def down(w, i):
            if i + 2 <= w.rank:
                b, c, a = w(i), w(i + 1), w(i + 2)
                if a < b < c:
                    s = SignedPermutation.s_i(i, w.rank)
                    t = SignedPermutation.s_i(i + 1, w.rank)
                    return w * s * t

        failures = 0

        # lengths
        for u, v in self.edges:
            try:
                if is_even(u):
                    assert len(v) == len(u) - 1
                else:
                    assert len(u) == len(v) - 1
            except:
                failures += 1

        # lemma 4.2(a)
        for u, v in self.edges:
            if not is_even(v):
                continue
            for i in range(2, self.n):
                uu = up(u, i)
                vv = up(v, i)
                try:
                    assert (vv is None) or ((uu, vv) in self.edges)
                except:
                    failures += 1
                    print('4.2(a)', u, '->', v, ', ', uu, '->', vv, ', ', i)

        # lemma 4.2(b)
        for v, w in self.edges:
            if not is_even(v):
                continue
            for i in range(2, self.n):
                vv = down(v, i)
                ww = down(w, i)
                try:
                    assert (vv is None) or ((vv, ww) in self.edges)
                except:
                    failures += 1
                    print('4.2(b)', v, '->', w, ', ', vv, '->', ww, ', ', i)

        q = {u for u, _ in self.edges} | {v for _, v in self.edges}
        q_odd = {u for u in q if not is_even(u)}
        q_even = q - q_odd
        atoms = {u for u in q_even if u(1) < 0 or any(a < u(1) < b == -a for a, b in q_shape(u))}
        q_even = q_even - atoms

        print('atoms:', len(atoms), ' even:', len(q_even), ' odd:', len(q_odd))

        # lemma 4.3
        for v in atoms:
            m = a_shape(v)
            if v(1) < 0:
                p = m - {(v(1), -v(1))}
                b = -v(1)
            else:
                a = v(1)
                b = [y for x, y in m if x == a][0]
                p = (m - {(a, b), (-b, -a)}) | {(-a, a)}
            try:
                assert not any(x < b < y for x, y in p)
                assert not any(x < c < y < d for x, y in p for c, d in p)
            except:
                failures += 1
            if is_maximal(v):
                pairs = sorted([(x, -y) for (x, y) in p if 0 < x < y])
                pairs = [x for pr in pairs for x in pr]
                fixed = sorted([x for (x, y) in p if x + y == 0])
                alpha = SignedPermutation(*([b] + fixed + pairs))
                try:
                    assert [alpha] == [u for (u, z) in self.edges if z == v]
                except:
                    print(alpha, ' v =', v, ' M\' =', p)
                    failures += 1

        # lemma 4.5(a) and theorem 4.8(a)
        for w in q_odd:
            b = w(1)
            for i in range(2, self.n):
                ww = up(w, i)
                if ww:
                    for a in range(1, b):
                        t = SignedPermutation.reflection_t(a, b, self.n)
                        if len(w) + 1 == len(t * w):
                            try:
                                assert up(t * w, i) == t * ww
                                assert len(ww) + 1 == len(t * ww)
                                assert (w, t * w) in self.edges
                            except:
                                print('w =', w, ' t =', t, ' w\' =', ww, ' i =', i, ' ? ', up(t * w, i), '!=', t * ww)
                                failures += 1

        # lemma 4.5(b) and theorem 4.8(b)
        for w in q_odd:
            b = w(1)
            for i in range(2, self.n):
                ww = down(w, i)
                if ww:
                    for c in range(b + 1, self.n + 1):
                        t = SignedPermutation.reflection_t(b, c, self.n)
                        if len(w) + 1 == len(t * w):
                            try:
                                assert down(t * w, i) == t * ww
                                assert len(ww) + 1 == len(t * ww)
                                assert (t * w, w) in self.edges
                            except:
                                print('w =', w, ' t =', t, ' w\' =', ww, ' i =', i, ' ? ', down(t * w, i), '!=', t * ww)
                                failures += 1

        # lemma 4.6
        for u in q_odd:
            p = q_shape(u)
            assert u(1) > 0
            m = p | {(-u(1), u(1))}
            try:
                assert not any(x < c < y < d for x, y in m for c, d in m)
                v = u * SignedPermutation.s_i(0, u.rank)
                assert v in atoms
                assert (u, v) in self.edges
            except:
                print(u, v, p, m)
                failures += 1

        # lemma 4.7
        for u in q_odd:
            p = q_shape(u)
            b = u(1)
            for a in range(1, b):
                if (-a, a) not in p:
                    continue
                t = SignedPermutation.reflection_t(a, b, u.rank)
                v = t * u
                if len(v) != len(u) + 1:
                    continue
                m = (p - {(-a, a)}) | {(a, b), (-b, -a)}
                try:
                    assert not any(x < c < y < d for x, y in m for c, d in m)
                    assert v in atoms
                    assert (u, v) in self.edges
                except:
                    print(u, v, a, b, p, m)
                    failures += 1

        print('failures:', failures)
