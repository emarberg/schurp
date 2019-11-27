import itertools
from tableaux import Tableau

REDUCED_WORDS = {(): {()}}
PIPE_DREAMS = {(): {((),)}}
ATOMS_CACHE = {}
FPF_ATOMS_CACHE = {}
SYMPLECTIC_HECKE_ATOMS_CACHE = {}


class Permutation:

    @property
    def rank(self):
        return len(self.oneline)

    @classmethod
    def get_grassmannian(cls, *mu):
        oneline = tuple(i + 1 + a for i, a in enumerate(sorted(mu)))
        if oneline:
            missing = set(range(1, oneline[-1] + 1)) - set(oneline)
            oneline += tuple(sorted(missing))
        return Permutation(*oneline).inverse()

    @classmethod
    def get_inv_grassmannian(cls, *mu):
        ans = Permutation()
        for i in range(len(mu)):
            ans *= Permutation.transposition(1 + mu[0] - mu[i], i + 1 + mu[0])
        return ans

    @classmethod
    def get_fpf_grassmannian(cls, *mu):
        ans = Permutation()
        o = 1 if mu and (len(mu) + 1 + mu[0]) % 2 != 0 else 0
        for i in range(len(mu)):
            ans *= Permutation.transposition(o + 1 + mu[0] - mu[i], o + i + 2 + mu[0])
        while not ans.is_fpf_involution():
            f = [i for i in range(1, ans.rank + 2) if ans(i) == i]
            ans *= Permutation.transposition(f[0], f[1])
        return ans

    @classmethod
    def random_reduced_word(cls, n):
        return Tableau.random_sorting_network(n)

    @classmethod
    def random_involution_word(cls, n):
        return Tableau.random_inv_network(n)

    @classmethod
    def random_fpf_involution_word(cls, n):
        return Tableau.random_fpf_network(n)

    @classmethod
    def longest_element(cls, n):
        return Permutation(tuple(range(n, 0, -1)))

    @classmethod
    def double_involution_word(cls, word):
        w = Permutation()
        ans = []
        for i in word:
            assert i not in w.right_descent_set
            s = Permutation.s_i(i)
            w *= s
            if i not in w.left_descent_set:
                ans = [i] + ans + [i]
                w = s * w
            else:
                ans = ans + [i]
        return tuple(ans)

    def reduced_word(self):
        return self.get_reduced_word()

    def get_reduced_word(self):
        w = self
        ans = ()
        while len(w) > 0:
            i = min(w.left_descent_set)
            ans += (i,)
            w = Permutation.s_i(i) * w
        return ans

    def get_reduced_words(self):
        oneline = tuple(self.oneline)
        if oneline not in REDUCED_WORDS:
            words = set()
            for i in self.right_descent_set:
                s = Permutation.s_i(i)
                words |= {e + (i,) for e in (self * s).get_reduced_words()}
            REDUCED_WORDS[oneline] = words
        return REDUCED_WORDS[oneline]

    def get_bottom_pipe_dream(self):
        from pipedreams import Pipedream
        code = self.code()
        crossings = set()
        for i, c in enumerate(code):
            crossings |= {(i + 1, j + 1) for j in range(c)}
        return Pipedream(crossings)

    def get_top_pipe_dream(self):
        return self.inverse().get_bottom_pipe_dream().transpose()

    def get_pipe_dreams(self):
        bottom = self.get_bottom_pipe_dream()
        return bottom.upper_ladder_interval()

    def count_pipe_dreams(self):
        ans = 0
        for p in self.get_bottom_pipe_dream().upper_ladder_interval():
            ans += 1
        return ans

    def get_bottom_involution_pipe_dream(self):
        assert self.is_involution()
        return self.get_min_atom().get_bottom_pipe_dream()

    def _get_involution_pipe_dreams_slow(self, extended=False):
        assert self.is_involution()
        for w in self.get_atoms():
            for dream in w.get_pipe_dreams():
                if extended or all(i >= j for (i, j) in dream.crossings):
                    yield dream

    def get_involution_pipe_dreams(self, extended=False):
        return self.get_bottom_involution_pipe_dream().upper_involution_ladder_interval(extended)

    def count_involution_pipe_dreams(self):
        ans = 0
        kappa = len([i for i in self.oneline if self(i) < i])
        for p in self.get_bottom_involution_pipe_dream().upper_involution_ladder_interval():
            diagonal_entries = len([a for a in p if a[0] == a[1]])
            ans += 2**(kappa - diagonal_entries)
        return ans

    def get_bottom_fpf_pipe_dream(self):
        assert self.is_fpf_involution()
        return self.get_min_fpf_atom().get_bottom_pipe_dream()

    def _get_fpf_involution_pipe_dreams_slow(self, extended=False):
        assert self.is_fpf_involution()
        for w in self.get_fpf_atoms():
            for dream in w.get_pipe_dreams():
                if extended or all(i > j for (i, j) in dream.crossings):
                    yield dream

    def get_fpf_involution_pipe_dreams(self, extended=False):
        return self.get_bottom_fpf_pipe_dream().upper_fpf_involution_ladder_interval(extended)

    def count_fpf_involution_pipe_dreams(self):
        ans = 0
        for p in self.get_bottom_fpf_pipe_dream().upper_fpf_involution_ladder_interval():
            ans += 1
        return ans

    @classmethod
    def _get_pipe_dreams_helper(cls, word, lowerbound=0, upperbound=None):
        if len(word) == 0:
            yield ((),)
            return
        if word[0] <= lowerbound:
            return
        for dream in cls._get_pipe_dreams_helper(word, lowerbound + 1):
            yield ((),) + dream
        if upperbound is None or word[0] < upperbound:
            for dream in cls._get_pipe_dreams_helper(word[1:], lowerbound, upperbound=word[0]):
                newdream = list(dream)
                newdream[0] = (word[0],) + newdream[0]
                yield tuple(newdream)

    def get_atoms(self):
        if self not in ATOMS_CACHE:
            ATOMS_CACHE[self] = list(self._get_atoms())
        return ATOMS_CACHE[self]

    def _get_atoms(self):
        def next(oneline):
            for i in range(len(oneline) - 2):
                c, a, b = oneline[i:i + 3]
                if a < b < c:
                    newline = oneline[:i] + (b, c, a) + oneline[i + 3:]
                    yield newline

        minimum = tuple(self.get_min_atom().inverse().oneline)
        add = {minimum}
        while add:
            for w in add:
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def get_fpf_atoms(self):
        if self not in FPF_ATOMS_CACHE:
            FPF_ATOMS_CACHE[self] = list(self._get_fpf_atoms())
        return FPF_ATOMS_CACHE[self]

    def _get_fpf_atoms(self):
        def next(oneline):
            for i in range(0, len(oneline) - 3, 2):
                a, d, b, c = oneline[i:i + 4]
                if a < b < c < d:
                    newline = oneline[:i] + (b, c, a, d) + oneline[i + 4:]
                    yield newline

        minimum = tuple(self.get_min_fpf_atom().inverse().oneline)
        add = {minimum}
        while add:
            for w in add:
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)}

    def get_symplectic_hecke_atoms(self):
        if self not in SYMPLECTIC_HECKE_ATOMS_CACHE:
            SYMPLECTIC_HECKE_ATOMS_CACHE[self] = list(self._get_symplectic_hecke_atoms())
        return SYMPLECTIC_HECKE_ATOMS_CACHE[self]

    def _get_symplectic_hecke_atoms(self):
        def next(w):
            for i in range(0, len(w) - 3, 2):
                a, d, b, c = w[i: i + 4]
                if a < b < c < d:
                    yield w[:i] + (b, c, a, d) + w[i + 4:]
                    yield w[:i] + (b, d, a, c) + w[i + 4:]
                b, c, a, d = w[i:i + 4]
                if a < b < c < d:
                    yield w[:i] + (a, d, b, c) + w[i + 4:]
                    yield w[:i] + (b, d, a, c) + w[i + 4:]
                b, d, a, c = w[i:i + 4]
                if a < b < c < d:
                    yield w[:i] + (a, d, b, c) + w[i + 4:]
                    yield w[:i] + (b, c, a, d) + w[i + 4:]

        minimum = tuple(self.get_min_fpf_atom().inverse().oneline)
        add = {minimum}
        seen = set()
        while add:
            for w in add:
                seen.add(w)
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)} - seen

    def get_involution_words(self):
        assert self.inverse() == self
        for a in self.get_atoms():
            for word in a.get_reduced_words():
                yield word

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield Permutation(args)

    @classmethod
    def involutions(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            w = Permutation(args)
            if all(w(w(i)) == i for i in range(1, n + 1)):
                yield w

    @classmethod
    def dominant_involutions(cls, n):
        for w in cls.involutions(n):
            if not any(w(i) < w(k) < w(j) for i in range(1, n - 1) for j in range(i + 1, n) for k in range(j + 1, n + 1)):
                yield w

    @classmethod
    def fpf_involutions(cls, n):
        s = {i: Permutation.s_i(i) for i in range(1, n)}
        if n % 2 == 0:
            start = Permutation()
            for i in range(1, n, 2):
                start *= s[i]
            level = {start}
            while level:
                next_level = set()
                for w in level:
                    yield w
                    for i in range(1, n):
                        if w(i) < w(i + 1):
                            next_level.add(s[i] * w * s[i])
                level = next_level

    def shift(self, n):
        assert n >= 0
        oneline = [i + 1 for i in range(n)] + [i + n for i in self.oneline]
        return Permutation(oneline)

    def standardize(self, e):
        index = {b: a + 1 for a, b in enumerate(sorted(map(self, e)))}
        oneline = [index[self(i)] for i in sorted(e)]
        return Permutation(oneline)

    def get_descentset_fpf_visible(self):
        return sorted(
            [ij[0] for ij in self.get_fpf_visible_inversions() if ij[1] == ij[0] + 1]
        )

    def get_descentset_visible(self):
        return sorted(
            [ij[0] for ij in self.get_visible_inversions() if ij[1] == ij[0] + 1]
        )

    def psi(self):
        """See transitions paper, Schur P-positivity section."""
        max_vis_inv = self.max_visible_inversion()
        if not max_vis_inv:
            return None
        q, r = max_vis_inv
        a = set([self(q), self(r), q, r])
        t = Permutation.cycle([q, r])
        if len(a) == 2:
            return self * t
        else:
            return t * self * t

    def get_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j):
                    ans.add((i, j))
        return ans

    def get_visible_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j) and i >= self(j):
                    ans.add((i, j))
        return ans

    def get_fpf_visible_inversions(self):
        n = len(self.oneline)
        ans = set()
        for i in range(1, n):
            for j in range(i + 1, n + 1):
                if self(i) > self(j) and i > self(j):
                    ans.add((i, j))
        return ans

    def max_inversion(self):
        return max(self.get_inversions())

    def max_visible_inversion(self):
        return max(self.get_visible_inversions())

    def max_fpf_visible_inversion(self):
        return max(self.get_fpf_visible_inversions())

    def support(self):
        return [i + 1 for i in range(len(self.oneline)) if self.oneline[i] != i + 1]

    def is_dominant(self):
        n = len(self.oneline)
        for i in range(1, n - 1):
            for j in range(i + 1, n + 1):
                for k in range(j + 1, n + 1):
                    if self(i) < self(k) < self(j):
                        return False
        return True

    def is_fpf_dominant(self):
        assert self.is_fpf_involution()
        diagram = set(self.fpf_rothe_diagram())
        columns = len({j for i, j in diagram})
        mu = tuple(len({i for i, j in diagram if j == k}) for k in range(1, columns + 1))
        return diagram == {
            (i + j + 1, i) for i in range(1, columns + 1) for j in range(mu[i - 1])
        }

    def rothe_diagram(self):
        n = len(self.oneline)
        return sorted([
            (i, self(j)) for i in range(1, n + 1) for j in range(i + 1, n + 1) if self(i) > self(j)
        ])

    def code_helper(self, diag):
        n = len(self.oneline)
        ans = n * [0]
        for i, j in diag:
            ans[i - 1] += 1
        while ans and ans[-1] == 0:
            ans = ans[:-1]
        return tuple(ans)

    def code(self):
        return self.code_helper(self.rothe_diagram())

    def involution_code(self):
        return self.code_helper(self.involution_rothe_diagram())

    def fpf_involution_code(self):
        ans = self.code_helper(self.involution_rothe_diagram(True))
        return ans + (len(self.oneline) - len(ans)) * (0,)

    def shape(self):
        from partitions import Partition
        return Partition(*list(reversed(sorted(self.code())))).transpose()

    def involution_shape(self):
        from partitions import Partition
        return Partition(*list(reversed(sorted(self.involution_code())))).transpose()

    def fpf_involution_shape(self):
        from partitions import Partition
        return Partition(*list(reversed(sorted(self.fpf_involution_code())))).transpose()

    @classmethod
    def from_code(cls, *code):
        if len(code) == 1 and type(code[0]) in [list, tuple]:
            code = code[0]
        if len(code) == 0 or code[-1] != 0:
            code = tuple(code) + (0,)
        n = len(code)
        indices = [i for i in range(n - 1) if code[i] != 0 and code[i + 1] == 0]
        if indices:
            i = indices[0]
            newcode = list(code)
            newcode[i + 1] = newcode[i] - 1
            newcode[i] = 0
            return cls.from_code(newcode) * Permutation.s_i(i + 1)
        else:
            return Permutation()

    @classmethod
    def from_involution_code(cls, *code):
        if len(code) == 1 and type(code[0]) in [list, tuple]:
            code = code[0]
        if len(code) == 0 or code[-1] != 0:
            code = tuple(code) + (0,)
        n = len(code)
        indices = [i for i in range(n - 1) if code[i] != 0 and code[i + 1] == 0]
        if indices:
            i = indices[0]
            newcode = list(code)
            newcode[i + 1] = newcode[i] - 1
            newcode[i] = 0
            s = Permutation.s_i(i + 1)
            w = cls.from_involution_code(newcode)
            return w * s if s * w == w * s else s * w * s
        else:
            return Permutation()

    @classmethod
    def from_fpf_involution_code(cls, *code):
        if len(code) == 1 and type(code[0]) in [list, tuple]:
            code = code[0]
        if len(code) % 2 != 0:
            code = tuple(code) + (0,)
        r = len(code)
        if len(code) == 0 or code[-1] != 0:
            code = tuple(code) + (0, 0)
        n = len(code)
        indices = [i for i in range(n - 1) if code[i] != 0 and code[i + 1] == 0]
        if indices:
            i = indices[0]
            newcode = list(code)
            newcode[i + 1] = newcode[i] - 1
            newcode[i] = 0
            s = Permutation.s_i(i + 1)
            ans = s * cls.from_fpf_involution_code(newcode) * s
            while len(ans.oneline) < r:
                ans *= Permutation.s_i(len(ans.oneline) + 1)
            return ans
        else:
            return Permutation.shortest_fpf_involution(r)

    @classmethod
    def shortest_fpf_involution(cls, rank):
        assert rank % 2 == 0
        oneline = [i + 2 * ((i + 1) % 2) for i in range(rank)]
        return Permutation(oneline)

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

    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) in [list, tuple]:
            oneline = [i for i in args[0]]
        else:
            oneline = list(args)
        while len(oneline) > 0 and oneline[-1] == len(oneline):
            oneline = oneline[:-1]
        self.oneline = oneline
        self.cycles = self.get_cycles(oneline)

    def is_right_descent(self, i):
        return self(i) > self(i + 1)

    @property
    def left_descent_set(self):
        return self.get_descentset_L()

    @property
    def right_descent_set(self):
        return self.get_descentset_R()

    def get_descentset_R(self):
        ans = []
        for i in range(1, len(self.oneline)):
            if(self(i) > self(i + 1)):
                ans.append(i)
        return set(ans)

    def get_descentset_L(self):
        return (self.inverse()).get_descentset_R()

    def number_two_cycles(self):
        ans = 0
        for c in self.cycles:
            if len(c) == 2:
                ans += 1
        return ans

    def contains_atom(self, w):
        expr = w.reduced_expr()
        return self.expr_to_involution(expr) == self

    def is_atom(self):
        expr = self.reduced_expr()
        return self.expr_to_involution(expr).involution_length() == len(expr)

    def get_min_atom(self):
        assert self.is_involution()
        cycles = sorted([list(reversed(sorted(c))) for c in self.cycles], key=lambda x: x[-1])
        return Permutation([i for cycle in cycles for i in cycle])**-1

    def get_max_atom(self):
        assert self.is_involution()
        cycles = sorted([list(reversed(sorted(c))) for c in self.cycles], key=lambda x: x[0])
        return Permutation([i for cycle in cycles for i in cycle])**-1

    def get_min_fpf_atom(self):
        assert self.is_fpf_involution()
        cycles = sorted([list(sorted(c)) for c in self.cycles], key=lambda x: x[0])
        return Permutation([i for cycle in cycles for i in cycle])**-1

    def get_max_fpf_atom(self):
        assert self.is_fpf_involution()
        cycles = sorted([list(sorted(c)) for c in self.cycles], key=lambda x: x[-1])
        return Permutation([i for cycle in cycles for i in cycle])**-1

    def is_involution(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 2

    def is_fpf_involution(self):
        return all(self(i) != i and self(self(i)) == i for i in self)

    def is_identity(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 1

    def __iter__(self):
        return self.oneline.__iter__()

    # Input is [i_1, i_2, ... , i_k], returns permutation (i_1 i_2 ... i_k)
    @staticmethod
    def cycle(cyc):
        return Permutation.cycles([cyc])

    @staticmethod
    def cycles(cyc, oneline=None):

        if len(cyc) == 0:
            return Permutation()

        if oneline is None:
            n = max(map(max, cyc))
            oneline = list(range(1, 1 + n))
            for c in cyc:
                oneline[c[len(c) - 1] - 1] = c[0]
                for i in range(len(c) - 1):
                    oneline[c[i] - 1] = c[i + 1]
        ans = Permutation()
        ans.oneline = oneline
        ans.cycles = cyc
        return ans

    @classmethod
    def s_i(cls, i):
        return cls.cycle([i, i + 1])

    @classmethod
    def t_ij(cls, i, j):
        return cls.transposition(i, j)

    @classmethod
    def transposition(cls, i, j):
        return cls.cycle([i, j])

    def tau_ij(self, i, j):
        assert self.is_involution() and i < j
        w = self
        a_tup = tuple(sorted(set([i, j, w(i), w(j)])))
        if len(a_tup) == 2 and w(i) == i:
            r = Permutation.t_ij(i, j)
            return r * w
        elif len(a_tup) == 3:
            a, b, c = a_tup
            if (i, j) in [(b, c), (a, c)] and w(a) == b and w(c) == c:
                r = Permutation.t_ij(b, c)
                return r * w * r
            elif (i, j) in [(a, b), (a, c)] and w(a) == a and w(b) == c:
                r = Permutation.t_ij(a, b)
                return r * w * r
        elif len(a_tup) == 4:
            a, b, c, d = a_tup
            if (i, j) == (b, c) and w(a) == b and w(c) == d:
                r = Permutation.t_ij(i, j)
                return r * w * r
            elif (i, j) in [(a, c), (b, d), (a, d)] and w(a) == b and w(c) == d:
                s = Permutation.t_ij(a, b)
                r = Permutation.t_ij(a, c)
                return r * s * w * r
            elif (i, j) in [(a, b), (c, d), (a, d)] and w(a) == c and w(b) == d:
                r = Permutation.t_ij(a, b)
                return r * w * r
        return w

    def inverse_tau_ij(self, i, j):
        assert self.is_involution() and i < j
        w = self
        a_tup = tuple(sorted(set([i, j, w(i), w(j)])))
        if len(a_tup) == 2 and w(i) == j:
            r = Permutation.t_ij(i, j)
            return r * w
        elif len(a_tup) == 3:
            a, b, c = a_tup
            if (i, j) in [(b, c), (a, c)] and w(a) == c and w(b) == b:
                r = Permutation.t_ij(b, c)
                return r * w * r
            elif (i, j) in [(a, b), (a, c)] and w(a) == c and w(b) == b:
                r = Permutation.t_ij(a, b)
                return r * w * r
        elif len(a_tup) == 4:
            a, b, c, d = a_tup
            if (i, j) == (b, c) and w(a) == c and w(b) == d:
                r = Permutation.t_ij(i, j)
                return r * w * r
            elif (i, j) in [(a, c), (b, d), (a, d)] and w(a) == d and w(b) == b and w(c) == c:
                s = Permutation.t_ij(a, b)
                r = Permutation.t_ij(a, c)
                return s * r * w * r
            elif (i, j) in [(a, b), (c, d), (a, d)] and w(a) == d and w(b) == c:
                r = Permutation.t_ij(a, b)
                return r * w * r
        return w

    def inverse(self):
        oneline = list(range(1, 1 + len(self.oneline)))
        for i in range(1, 1 + len(self.oneline)):
            oneline[self(i) - 1] = i
        return Permutation(oneline)

    def __pow__(self, e):
        if e < 0:
            return self.inverse()**(-e)
        m = 1
        for ell in list(map(len, self.cycles)):
                m = m * ell
        e = e % m
        x = Permutation()
        for i in range(e):
            x = x * self
        return x

    # other is a permutation
    def __mul__(self, other):
        assert type(other) == Permutation
        newline = []
        n = max(len(self.oneline), len(other.oneline))
        for i in range(1, 1 + n):
            newline.append(self(other(i)))
        while len(newline) > 0 and newline[-1] == len(newline):
            newline = newline[:-1]
        return Permutation(newline)

    def __len__(self):
        if(len(self.cycles) == 0):
            return 0
        return len(self.get_inversions())

    def length(self):
        return len(self)

    def __call__(self, i):
        if i < 1:
            return i
        if i > len(self.oneline):
            return i
        return self.oneline[i - 1]

    def get_cycles(self, oneline):
        cycles = []
        numbers = list(range(1, len(oneline) + 1))
        while len(numbers) > 0:
            cycle = [numbers[0]]
            next = oneline[numbers[0] - 1]
            numbers.remove(numbers[0])

            while next != cycle[0]:
                cycle.append(next)
                numbers.remove(next)
                next = oneline[cycle[len(cycle) - 1] - 1]
            cycles.append(cycle)
        return cycles

    # Lexicographic order on one-line representation
    def __lt__(self, other):
        for i in range(1, 1 + max(max(self.oneline), max(other.oneline))):
            if self(i) >= other(i):
                return False
        return True

    def __eq__(self, other):
        if not isinstance(other, Permutation):
            return False
        return len(self * (other**-1)) == 0

    def __ne__(self, other):
        return not (self == other)

    def strong_bruhat_less_than(self, other):
        if self.length() >= other.length():
            return False
        des = self.get_descentset_L()
        for k in des:
            x = sorted(map(self, range(1, k + 1)))
            y = sorted(map(other, range(1, k + 1)))
            for i in range(k):
                if x[i] > y[i]:
                    return False
        return True

    def __repr__(self):
        sep = '' if len(self.oneline) < 10 else ','
        ans = sep.join([str(i) for i in self.oneline])
        return ans if ans else '()'
        # return self.cycle_repr()

    def cycle_repr(self):
        if len(self) == 0:
            return '1'
        EXCLUDE_SINGLE_CYCLES = True
        SPACE = ' '
        DELIM = ''  # ','

        s = ''
        for c in self.cycles:
            if not EXCLUDE_SINGLE_CYCLES or len(c) > 1:
                s += '(' + (DELIM + SPACE).join([str(x) for x in c]) + ')'
        return s

    def __hash__(self):
        return hash(tuple(self.oneline))

    def involution_length(self):
        return (self.length() + len(list(filter(lambda i: len(i) > 1, self.cycles)))) // 2

    def fpf_involution_length(self):
        assert self.is_fpf_involution()
        return (self.length() - len(list(filter(lambda i: len(i) > 1, self.cycles)))) // 2

    def fpf_trim(self):
        assert self.is_fpf_involution()
        n = len(self.oneline)
        return self if (n == 0 or self(n - 1) != n) else (self * self.s_i(n - 1)).fpf_trim()

    @classmethod
    def from_word(cls, *args):
        if len(args) == 1 and type(args[0]) == tuple:
            args = args[0]
        w = Permutation()
        for i in args:
            w *= Permutation.s_i(i)
        return w

    @classmethod
    def from_involution_word(cls, *word):
        w = Permutation()
        for i in word:
            s = Permutation.s_i(i)
            if i in w.right_descent_set:
                raise Exception
            elif s * w == w * s:
                w = w * s
            else:
                w = s * w * s
        return w

    @classmethod
    def from_fpf_involution_word(cls, *word):
        n = (1 + max(word)) if word else 0
        n = n if n % 2 == 0 else n + 1
        w = Permutation()
        for i in range(1, n, 2):
            w *= cls.s_i(i)
        for i in word:
            s = Permutation.s_i(i)
            if i in w.right_descent_set:
                raise Exception
            elif s * w == w * s:
                raise Exception
            else:
                w = s * w * s
        return w
