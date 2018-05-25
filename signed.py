import itertools
import subprocess
from collections import defaultdict
from vectors import Vector
from symmetric import SchurP, SchurQ
from partitions import StrictPartition
from permutations import Permutation


SIGNED_REDUCED_WORDS = {(): {()}}

schurp_stansym_cache = {}
schurq_stansym_cache = {}
schurd_stansym_cache = {}
atoms_b_cache = {}
atoms_d_cache = {}


class SignedPermutation:

    def __init__(self, *oneline):
        self.oneline = tuple(oneline)
        self.rank = len(oneline)
        assert set(range(1, self.rank + 1)) == set(abs(i) for i in self.oneline)
        # cached fields
        self._rdes = None
        self._ldes = None
        self._len = None

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

    def get_reduced_word(self):
        if self.left_descent_set:
            i = min(self.left_descent_set)
            s = SignedPermutation.s_i(i, self.rank)
            return (i,) + (s * self).get_reduced_word()
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

    # def inv_stanley_schur_d_decomposition(self):
    #     assert self.is_even_signed()
    #     ans = Vector()
    #     for x in self.get_atoms(type_d=True):
    #         for sh, i in x.stanley_schur_d_decomposition().items():
    #             ans += Vector({SchurP(StrictPartition(*sh)): i})
    #     return ans

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
        # print('  ', bcd_type, ':', len(cache))
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

    # @classmethod
    # def get_atom_shape(cls, oneline):
    #     while oneline and oneline[-1] == len(oneline):
    #         oneline = oneline[:-1]
    #     oneline = oneline[1:]
    #     ans = []
    #     while oneline:
    #         for i in range(len(oneline)):
    #             a = oneline[i]
    #             if i == 0 and a < 0:
    #                 ans += [(a, -a)]
    #                 oneline = oneline[1:]
    #                 break
    #             if i + 1 >= len(oneline):
    #                 continue
    #             b = oneline[i + 1]
    #             if 0 < a < -b:
    #                 ans += [(a, -b)]
    #                 oneline = oneline[:i] + oneline[i + 2:]
    #                 break
    #     return ans

    @classmethod
    def standardize(cls, oneline):
        distinct = sorted([abs(j) for j in oneline])
        assert len(set(distinct)) == len(oneline)
        mapping = {v: i + 1 for i, v in enumerate(distinct)}
        newline = tuple(mapping[v] if v > 0 else -mapping[-v] for v in oneline)
        return SignedPermutation(*newline)

        # if i.rank == self.n:
        #     try:
        #         assert self.standardize(i.inverse()) in self.sub_atoms
        #     except:
        #         print('!', i.inverse())

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
        #    return {sh: 1}

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
                # print('Computing atoms (B%s). . .' % (n - r - 1))
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

        # if len(w * self.reflection_t(r, s, n)) != len(w) - 1:
        #     print('failed')
        #     r = w.last_descent()
        #     s = w.last_inversion(r)

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
        # ans = defaultdict(int)
        # for x in indices:
        #     for sh, i in x.stanley_schur_decomposition(bcd_type, starting_rank, verbose, step + 1).items():
        #         ans[sh] += i
        # ans = dict(ans)

        # cache[self] = ans
        # # print('  ', bcd_type, ':', len(cache))
        # return ans

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
        return SignedPermutation(*self._min_inv_atom_oneline())

    def get_atoms(self):
        w = self.reduce()
        if w not in atoms_b_cache:
            atoms_b_cache[w] = list(w._get_atoms())
        ans = atoms_b_cache[w]
        return [x.inflate(self.rank) for x in ans]

    def _get_atoms(self):
        def next(oneline):
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
                    yield newline

        minimum = self._min_inv_atom_oneline()
        add = {minimum}
        while add:
            for w in add:
                yield SignedPermutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def involution_words(self):
        for w in self.get_atoms():
            for e in w.get_reduced_words():
                yield e

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


class SignedAtomsGraph:

    DIRECTORY = '/Users/emarberg/Dropbox/schubert-type-b/atoms/'

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
            self._edges = list(SignedPermutation.queue_stanley_decomposition(self.n))
        return self._edges


