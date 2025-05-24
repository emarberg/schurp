import itertools
from tableaux import Tableau
from partitions import Partition
from words import Word
import operator
import random
from collections import deque

REDUCED_WORDS = {(): {()}}
HECKE_WORDS = {}

INVOLUTION_WORDS = {(): {()}}
FPF_INVOLUTION_WORDS = {(): {()}}
PRIMED_INVOLUTION_WORDS = {(): {()}}
PIPE_DREAMS = {(): {((),)}}
ATOMS_CACHE = {}

TWISTED_ATOMS_CACHE = {}
TWISTED_HECKE_ATOMS_CACHE = {}
TWISTED_INVOLUTION_WORDS_CACHE = {}
TWISTED_PRIMED_INVOLUTION_WORDS_CACHE = {}
FPF_ATOMS_CACHE = {}

SYMPLECTIC_HECKE_ATOMS_CACHE = {}
INVOLUTION_HECKE_ATOMS_CACHE = {}
EXTENDED_HECKE_ATOMS_CACHE = {}

K_BRUHAT_COVERS = {}
K_PIERI_CHAINS = {}

class Permutation:

    def identity(self, n=None):
        return Permutation()

    def star(self, n):
        return Permutation(*[n + 1 - self(i) for i in range(n, 0, -1)])

    def rank_table(self, p, q):
        ans = 0
        for i in range(1, p + 1):
            if self(i) <= q:
                ans += 1
        return ans

    def print_rank_table(self):
        n = self.rank
        for p in range(1, n + 1):
            print(' '.join([str(self.rank_table(p, q)) for q in range(1, n + 1)]))

    @classmethod
    def inversions(cls, word):
        inv = []
        for index, i in enumerate(word):
            a = i
            b = i + 1
            for k in word[index + 1:]:
                a = Permutation.s_i(k)(a)
                b = Permutation.s_i(k)(b)
            assert a < b
            inv.append((a, b))
        return tuple(inv)

    @classmethod
    def northeast_chains(cls, p, q, size):
        ans = set()
        for rows in itertools.combinations(range(1, p + 1), size):
            for cols in itertools.combinations(range(1, q + 1), size):
                ans.add(tuple(zip(rows, reversed(cols))))
        return ans

    @classmethod
    def reflected_northeast_chains(cls, p, q, size):
        ans = set()
        for chain in cls.northeast_chains(p, q, size):
            if not any(a == b for (a, b) in chain):
                ans.add(tuple(sorted({(a, b) if a > b else (b, a) for (a, b) in chain})))
        return ans

    def all_northeast_chains(self):
        n = self.rank
        ans = set()
        for p in range(1, n + 1):
            for q in range(1, n + 1):
                ans |= self.northeast_chains(p, q, 1 + self.rank_table(p, q))
        return ans

    def essential_northeast_chains(self):
        diagram = self.rothe_diagram()
        ans = set()
        for (p, q) in diagram:
            if (p + 1, q) not in diagram and (p, q + 1) not in diagram:
                ans |= self.northeast_chains(p, q, 1 + self.rank_table(p, q))
        return ans

    def all_reflected_northeast_chains(self):
        n = self.rank
        ans = set()
        for p in range(1, n + 1):
            for q in range(1, p + 1):
                ans |= self.reflected_northeast_chains(p, q, 1 + self.rank_table(p, q))
        return ans

    def essential_reflected_northeast_chains(self):
        diagram = self.fpf_rothe_diagram()
        ans = set()
        for (p, q) in diagram:
            if (p + 1, q) not in diagram and (p, q + 1) not in diagram:
                ans |= self.reflected_northeast_chains(p, q, 1 + self.rank_table(p, q))
        return ans

    @property
    def rank(self):
        return len(self.oneline)

    @classmethod
    def reflections(cls, n):
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                yield cls.t_ij(i, j)

    @classmethod
    def get_grassmannian(cls, *mu):
        oneline = tuple(i + 1 + a for i, a in enumerate(sorted(mu)))
        if oneline:
            missing = set(range(1, oneline[-1] + 1)) - set(oneline)
            oneline += tuple(sorted(missing))
        return Permutation(*oneline)

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
    def find_commutations(cls, iword):
        ans = []
        w = Permutation()
        for index, i in enumerate(iword):
            i = abs(i)
            s = Permutation.s_i(i)
            if w(i) == i and w(i + 1) == i + 1:
                ans.append(index)
                w = w * s
            else:
                w = s * w * s
        return ans

    @classmethod
    def random_primed_involution_word(cls, n):
        w = list(cls.random_involution_word(n))
        for i in cls.find_commutations(w):
            if random.randint(0, 1):
                w[i] *= -1
        return tuple(w)

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

    def get_hecke_words(self, length):
        oneline = tuple(self.oneline)
        key = (oneline, length)
        if key not in HECKE_WORDS:
            if length < 0:
                words = set()
            elif length == 0:
                words = {()} if self.is_identity() else set()
            else:
                words = {e + (i,) for e in self.get_hecke_words(length - 1) for i in self.right_descent_set}
                for i in self.right_descent_set:
                    s = Permutation.s_i(i)
                    words |= {e + (i,) for e in (self * s).get_hecke_words(length - 1)}
            HECKE_WORDS[key] = words
        return HECKE_WORDS[key]

    def get_increasing_factorizations(self, k):
        for w in self.get_reduced_words():
            for f in Word.increasing_factorizations(w, k):
                yield tuple(_.elements for _ in f)

    def get_bounded_increasing_factorizations(self, k=None, flag=None):
        if k is None:
            k = self.rank - 1
        if flag is None:
            flag = list(range(1, k + 1))
        for f in self.get_increasing_factorizations(k):
            if all(a[0] >= flag[-i - 1] for i, a in enumerate(f) if a):
                yield f

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

    def get_twisted_primed_involution_words(self, n):
        assert self.is_twisted_involution(n)
        if self not in TWISTED_PRIMED_INVOLUTION_WORDS_CACHE:
            TWISTED_PRIMED_INVOLUTION_WORDS_CACHE[(self, n)] = list(self._get_twisted_primed_involution_words(n))
        return TWISTED_PRIMED_INVOLUTION_WORDS_CACHE[(self, n)]

    def _get_twisted_primed_involution_words(self, n):
        if self.is_identity():
            return [()]
        words = set()
        for i in self.right_descent_set:
            s = Permutation.s_i(i)
            t = Permutation.s_i(n - i)
            if t * self == self * s:
                sub = (self * s).get_twisted_primed_involution_words(n)
                words |= {e + (i,) for e in sub} | {e + (-i,) for e in sub}
            else:
                sub = (t * self * s).get_twisted_primed_involution_words(n)
                words |= {e + (i,) for e in sub}
        return words

    def get_twisted_involution_words(self, n):
        assert self.is_twisted_involution(n)
        if self not in TWISTED_INVOLUTION_WORDS_CACHE:
            TWISTED_INVOLUTION_WORDS_CACHE[(self, n)] = list(self._get_twisted_involution_words(n))
        return TWISTED_INVOLUTION_WORDS_CACHE[(self, n)]

    def _get_twisted_involution_words(self, n):
        if self.is_identity():
            return [()]
        ans = []
        for i in self.right_descent_set:
            s = Permutation.s_i(i)
            t = Permutation.s_i(n - i)
            w = t * self * s
            if self == w:
                w = self * s
            ans += [e + (i,) for e in w.get_twisted_involution_words(n)]
        return ans

    def get_twisted_atoms(self, n, offset=None):
        assert self.is_twisted_involution(n)
        key = (self, n, offset)
        if key not in TWISTED_ATOMS_CACHE:
            TWISTED_ATOMS_CACHE[key] = list(self._get_twisted_atoms(n, offset))
        return TWISTED_ATOMS_CACHE[key]

    def _get_twisted_atoms(self, n, offset):
        if offset is None:
            base = Permutation()
        else:
            cyc = []
            for l in range(offset + 1, n):
                if l >= n + 1 - l:
                    break
                cyc.append((l, n + 1 - l))
            base = Permutation.from_cycles(*cyc)
        if len(self) <= len(base) and self != base:
            return []
        if self == base:
            return [Permutation()]
        ans = set()
        for i in self.right_descent_set:
            s = Permutation.s_i(i)
            t = Permutation.s_i(n - i)
            w = t * self * s
            if self == w:
                w = self * s
            ans |= {e * s for e in w.get_twisted_atoms(n, offset)}
        return list(ans)

    def fixed(self, n):
        return {i for i in range(1, n + 1) if self(i) == i}

    def twisted_shape(self, n):
        y = self.star(n).inverse() % self
        assert y.twisted_involution_length(n) == self.length()

        w = self.inverse()
        line = [w(i) for i in range(1, n + 1)]
        pairs = []
        while line:
            pairs += [(line[0], line[-1])]
            line = line[1:-1]

        return {(a, b) for b, a in pairs if a < b}

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

    def get_fpf_involution_words(self):
        assert self.is_fpf_involution()
        oneline = tuple(self.fpf_trim().oneline)
        if oneline not in FPF_INVOLUTION_WORDS:
            words = set()
            for i in self.right_descent_set:
                if self(i) == i + 1:
                    continue
                s = Permutation.s_i(i)
                w = self
                sws = s * w * s
                words |= {e + (i,) for e in sws.get_fpf_involution_words()}
            FPF_INVOLUTION_WORDS[oneline] = words
        return FPF_INVOLUTION_WORDS[oneline]

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

    def get_fpf_involution_word(self):
        assert self.is_fpf_involution()
        for a in self.get_fpf_atoms():
            return a.get_reduced_word()

    def get_involution_word(self):
        assert self.inverse() == self
        for a in self.get_atoms():
            return a.get_reduced_word()

    def get_involution_words(self):
        assert self.inverse() == self
        oneline = tuple(self.oneline)
        if oneline not in INVOLUTION_WORDS:
            words = set()
            for i in self.right_descent_set:
                s = Permutation.s_i(i)
                w = self
                ws = w * s
                if s * ws != w:
                    ws = s * ws
                words |= {e + (i,) for e in ws.get_involution_words()}
            INVOLUTION_WORDS[oneline] = words
        return INVOLUTION_WORDS[oneline]

    def get_primed_involution_words(self):
        assert self.inverse() == self
        oneline = tuple(self.oneline)
        if oneline not in PRIMED_INVOLUTION_WORDS:
            words = set()
            for i in self.right_descent_set:
                s = Permutation.s_i(i)
                if s * self == self * s:
                    sub = (self * s).get_primed_involution_words()
                    words |= {e + (i,) for e in sub} | {e + (-i,) for e in sub}
                else:
                    sub = (s * self * s).get_primed_involution_words()
                    words |= {e + (i,) for e in sub}
            PRIMED_INVOLUTION_WORDS[oneline] = words
        return PRIMED_INVOLUTION_WORDS[oneline]

    def get_involution_hecke_words(self):
        ans = set()
        for v in self.get_involution_hecke_atoms():
            ans |= v.get_reduced_words()
        return ans

    def get_extended_hecke_atoms(self):
        assert self.is_involution()
        if self not in EXTENDED_HECKE_ATOMS_CACHE:
            h = set(self.get_involution_hecke_atoms())
            j = min([i for (i, j) in self.rothe_diagram() if i == j], default=1) - 1
            k = max([i for (i, j) in self.rothe_diagram() if i == j], default=1)
            s = {w for v in h for w, _, _, _ in v.k_pieri_chains(k, k, lowerbound=j, z=self)}
            EXTENDED_HECKE_ATOMS_CACHE[self] = list(s)
        return EXTENDED_HECKE_ATOMS_CACHE[self]

    def get_twisted_hecke_atoms(self, n):
        if (self, n) not in TWISTED_HECKE_ATOMS_CACHE:
            w0 = Permutation.longest_element(n)
            dictionary = {}
            for w in Permutation.all(n):
                z = (w0 * w * w0).inverse() % w
                if (z, n) not in dictionary:
                    dictionary[z, n] = set()
                dictionary[z, n].add(w)
            for (key, val) in dictionary.items():
                if key not in TWISTED_HECKE_ATOMS_CACHE:
                    TWISTED_HECKE_ATOMS_CACHE[key] = val
        return TWISTED_HECKE_ATOMS_CACHE[self, n]

    def get_involution_hecke_atoms(self):
        if self not in INVOLUTION_HECKE_ATOMS_CACHE:
            INVOLUTION_HECKE_ATOMS_CACHE[self] = list(self._get_involution_hecke_atoms())
        return INVOLUTION_HECKE_ATOMS_CACHE[self]

    def _get_involution_hecke_atoms(self):
        def next(w):
            for i in range(0, len(w) - 2):
                c, b, a = w[i: i + 3]
                if a < b < c:
                    yield w[:i] + (b, c, a) + w[i + 3:]
                    yield w[:i] + (c, a, b) + w[i + 3:]
                b, c, a = w[i:i + 3]
                if a < b < c:
                    yield w[:i] + (c, b, a) + w[i + 3:]
                    yield w[:i] + (c, a, b) + w[i + 3:]
                c, a, b = w[i:i + 3]
                if a < b < c:
                    yield w[:i] + (c, b, a) + w[i + 3:]
                    yield w[:i] + (b, c, a) + w[i + 3:]

        minimum = tuple(self.get_min_atom().inverse().oneline)
        add = {minimum}
        seen = set()
        while add:
            for w in add:
                seen.add(w)
                yield Permutation(*w).inverse()
            add = {new for w in add for new in next(w)} - seen

    def get_twisted_involution_hecke_words(self, n):
        ans = set()
        for v in self.get_twisted_involution_hecke_atoms(n):
            ans |= v.get_reduced_words()
        return ans

    def get_twisted_involution_hecke_atoms(self, n):
        ans = set()
        n = self.rank
        for v in self.all(n):
            if v.inverse().star(n) % v == self:
                ans.add(v)
        return ans

    def dominant_component(self):
        ans = set()
        i = 1
        while self(i) > 1:
            j = 1
            while not any(y == self(x) for x in range(1, i + 1) for y in range(1, j +1)):
                ans.add((i, j))
                j += 1
            i += 1
        return ans

    def outer_corners(self):
        n = self.rank + 1
        dom = self.dominant_component()
        return {
            (i, j)
            for i in range(1, n + 1)
            for j in range(1, n + 1)
            if (i, j) not in dom and
            (i == 1 or (i - 1, j) in dom) and
            (j == 1 or (i, j - 1) in dom)
        }

    @classmethod
    def all(cls, n):
        for args in itertools.permutations(range(1, n + 1)):
            yield Permutation(args)

    def twisted_conjugacy_class(self, n):
        return self.conjugacy_class(n, True)

    def conjugacy_class(self, n, twisted=False):
        t = self.longest_element(n) if twisted else self.identity()
        g = [self.s_i(i) for i in range(1, n)]
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

    @classmethod
    def twisted_involutions(cls, n):
        level = {Permutation()}
        while level:
            new = set()
            for w in level:
                yield w
                for i in range(1, n):
                    if w(i) < w(i + 1):
                        s = Permutation.s_i(i)
                        t = Permutation.s_i(n - i)
                        v = t * w * s
                        if w == v:
                            new.add(w * s)
                        else:
                            new.add(v)
            level = new

    @classmethod
    def partial_involutions(cls, n, macaulay_output=False, partial=True):
        if not partial:
            ans = list(cls.involutions(n))
        elif n == 0:
            ans = [Permutation()]
        elif n == 1:
            ans = [Permutation(1, 2), Permutation(2, 1)]
        else:
            ans = []
            c = Permutation.cycle(list(range(n, 2 * n)))
            d = Permutation.cycle(list(range(n - 1, 2 * n + 1)))**2
            for w in cls.partial_involutions(n - 1, False):
                w = c * w * c**-1
                ans.append(w)
                m = max([n + 1] + [w(i) + 1 for i in range(1, n)])
                ans.append(w * w.t_ij(n, m))
            for w in cls.partial_involutions(n - 2, False):
                w = d * w * d**-1
                for k in range(1, n):
                    e = Permutation.cycle(list(range(k, n)))
                    ans.append(e * w * e**-1 * w.t_ij(k, n))
        if macaulay_output:
            s = '{' + ', '.join(['{ ' +  ', '.join([str(w(i + 1)) for i in range(max(n, w.rank))] + [str(n)]) + ' }' for w in ans]) + '};'
            import pyperclip
            pyperclip.copy(s)
        return ans

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

    def fpf_shift(self, n):
        assert n >= 0 and n % 2 == 0
        oneline = [i + 1 + (-1)**i for i in range(n)] + [i + n for i in self.oneline]
        return Permutation(oneline)

    def standardize(self, e):
        index = {b: a + 1 for a, b in enumerate(sorted(map(self, sorted(e))))}
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

    def dearc_L(self):
        assert self.is_fpf_involution()
        z = self
        oneline = [z(i) if any(j < z(j) for j in range(min(i, z(i)) + 1,  max(i, z(i)))) else i for i in range(1, len(self.oneline) + 1)]
        return Permutation(*oneline)

    def dearc_R(self):
        assert self.is_fpf_involution()
        z = self
        oneline = [z(i) if any(j > z(j) for j in range(min(i, z(i)) + 1,  max(i, z(i)))) else i for i in range(1, len(self.oneline) + 1)]
        return Permutation(*oneline)

    def is_fully_commutative(self):
        for i in range(1, self.rank + 1):
            for j in range(i + 1, self.rank + 1):
                for k in range(j + 1, self.rank + 1):
                    if self(i) > self(j) > self(k):
                        return False
        return True

    def skew_shape(self):
        assert self.is_fully_commutative()
        c = self.code()
        k = [i for i in range(len(c)) if c[i] > 0]

        pairs = [(i, j) for i in range(1, len(k) + 1) for j in range(i + 1 - k[i - 1] - c[k[i - 1]], i + 1 - k[i - 1])]
        m = max([0] + [1 - j for (i, j) in pairs])
        pairs = [(i, j + m) for (i, j) in pairs]
        print(pairs)
        lam = set()
        for i, j in pairs:
            for a in range(1, i + 1):
                for b in range(1, j + 1):
                    lam.add((a, b))
        mu = []
        for i, j in lam:
            while len(mu) < i:
                mu.append(0)
            mu[i - 1] = max(mu[i - 1], j)
        mu = tuple(mu)

        nu = list(mu)
        for i, j in pairs:
            print(nu, i, j)
            nu[i - 1] = min(nu[i - 1], j - 1)
        nu = tuple(nu)

        return mu, nu


    def is_sp_vexillary(self):
        assert self.is_fpf_involution()
        raise Exception

    def is_fpf_vexillary(self):
        assert self.is_fpf_involution()
        raise Exception

    def is_fpf_grassmannian(self):
        assert self.is_fpf_involution()
        raise Exception
    
    def get_essential_path(self):
        assert self.is_involution()
        assert self.is_vexillary()

        n = self.rank
        rothe = self.involution_rothe_diagram()
        ess = [(i, j) for (i, j) in rothe if (i + 1, j) not in rothe and (i, j + 1) not in rothe]
        ess = sorted(ess, reverse=True)
        
        path = [(n, 0)]
        while ess:
            a0, b0 = path[-1]
            a1, b1 = ess[0]
            ess = ess[1:]
            for i in range(a0 - 1, a1 - 1, -1):
                path.append((i, b0))
            for j in range(b0 + 1, b1 + 1):
                path.append((a1, j))
        while path[-1][0] > path[-1][1]:
            a, b = path[-1]
            path.append((a - 1, b))
        assert len(path) == n + 1
        return list(reversed(path))

    def is_grassmannian(self):
        return len(self.right_descent_set) <= 1

    def is_inv_grassmannian(self):
        return len(self.get_descentset_visible()) <= 1

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

    def print_matrix(self, sep='.'):
        n = self.rank
        ans = []
        for i in self.oneline:
            row = n * [sep]
            row[i - 1] = '*'
            ans += [' '.join(row)]
        print('\n'.join(ans))

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
    def from_shape(cls, *mu):
        return cls.from_code(*mu)

    @classmethod
    def from_cycles(cls, *c):
        return cls.cycles(c)

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
            ans = cls.from_code(newcode) * Permutation.s_i(i + 1)
        else:
            ans = Permutation()
        if ans.code() != Partition.trim(code):
            print('         code:', ans.code())
            print('expected code:', Partition.trim(code))
        assert ans.code() == Partition.trim(code)
        return ans

    @classmethod
    def from_involution_shape(cls, *mu):
        if len(mu) == 0:
            return Permutation()
        shape = {(i, i + j - 1) for i in range(1, len(mu) + 1) for j in range(1, mu[i - 1] + 1)}
        code = mu[0] * [0]
        for i, j in shape:
            code[j - 1] += 1
        ans = cls.from_involution_code(*code)
        assert sorted(ans.involution_shape()) == sorted([m for m in mu if m != 0])
        return ans

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
            ans = w * s if s * w == w * s else s * w * s
        else:
            ans = Permutation()
        assert ans.is_involution()
        if ans.involution_code() != Partition.trim(code):
            print('involution code:', ans.involution_code())
            print('  expected code:', Partition.trim(code))
        assert ans.involution_code() == Partition.trim(code)
        return ans

    @classmethod
    def from_fpf_involution_shape(cls, *mu):
        if len(mu) == 0:
            return Permutation()
        shape = {(i, i + j - 1) for i in range(1, len(mu) + 1) for j in range(1, mu[i - 1] + 1)}
        code = (mu[0] + 1) * [0]
        for i, j in shape:
            code[j] += 1
        return cls.from_fpf_involution_code(*code)

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
        else:
            ans = Permutation.shortest_fpf_involution(r)
        assert ans.is_fpf_involution()
        ans_code = Partition.trim(ans.fpf_involution_code())
        exp_code = Partition.trim(code)
        assert ans_code == exp_code
        return ans

    @classmethod
    def shortest_fpf_involution(cls, rank):
        assert rank % 2 == 0
        oneline = [i + 2 * ((i + 1) % 2) for i in range(rank)]
        return Permutation(oneline)

    def fpf_rothe_diagram(self, fpf=False):
        return self.involution_rothe_diagram(True)

    def involution_rothe_diagram(self, fpf=False):
        return [(i, j) for (i, j) in self.rothe_diagram() if i > j or (not fpf and i == j)]

    def print_essential_set(self, french=False, sep='.'):
        rothe = self.rothe_diagram()
        ess = [(i, j) for (i, j) in rothe if (i + 1, j) not in rothe and (i, j + 1) not in rothe]
        print(self.print_diagram(ess, french=french, sep=sep))

    def print_rothe_diagram(self, french=False, sep='.'):
        print(self.print_diagram(self.rothe_diagram(), french=french, sep=sep))

    def print_fpf_rothe_diagram(self, french=False, sep='.'):
        print(self.print_diagram(self.involution_rothe_diagram(True), french=french, sep=sep))

    def print_involution_rothe_diagram(self, french=False, sep='.'):
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

    def sgn(self):
        return (-1)**len([c for c in self.cycles if len(c) % 2 == 0])

    def sign(self):
        return self.sgn()

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

    def get_two_cycles(self):
        ans = []
        for c in self.cycles:
            if len(c) == 2:
                ans += [tuple(c)]
        return ans

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

    def twisted_fixed_points(self, n):
        return [i for i in range(1, n + 1) if n + 1 - self(i) == i]

    def twisted_cycles(self, n):
        return [(i, n + 1 - self(i)) for i in range(1, n + 1) if n + 1 - self(i) > i]

    def _get_min_twisted_matching(self, n, matching):
        if matching is None:
            fix = self.twisted_fixed_points(n)
            matching = []
            if len(fix) % 2 != 0:
                fix = fix[:-1]
            while len(fix) >= 2:
                a, b = fix[0], fix[1]
                matching.append((a, b))
                fix = fix[2:]
        return matching

    def _get_max_twisted_matching(self, n, matching):
        if matching is None:
            fix = self.twisted_fixed_points(n)
            matching = []
            if len(fix) % 2 != 0:
                fix = fix[:-1]
            while len(fix) >= 2:
                a, b = fix[0], fix[1]
                matching.append((a, b))
                fix = fix[2:]
        return matching

    def get_min_twisted_atom(self, n, matching=None):
        assert self.is_twisted_involution(n)
        matching = self._get_min_twisted_matching(n, matching)
        base = list(set(self.twisted_fixed_points(n)) - {a for m in matching for a in m})
        assert len(base) <= 1

        itemgetter = operator.itemgetter(0)
        cycles = sorted([(b, a) for (a, b) in matching] + self.twisted_cycles(n), key=itemgetter)
        for x, y in reversed(cycles):
            base = [x] + base + [y]
        return Permutation(*base).inverse()

    def get_max_twisted_atom(self, n, matching=None):
        assert self.is_twisted_involution(n)
        matching = self._get_max_twisted_matching(n, matching)
        base = list(set(self.twisted_fixed_points(n)) - {a for m in matching for a in m})
        assert len(base) <= 1

        itemgetter = lambda p: -p[1]  # noqa
        cycles = sorted([(b, a) for (a, b) in matching] + self.twisted_cycles(n), key=itemgetter)
        for x, y in reversed(cycles):
            base = [x] + base + [y]
        return Permutation(*base).inverse()

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

    def is_twisted_involution(self, n):
        return n >= len(self.oneline) and self.inverse() == self.star(n)

    def is_involution(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 2

    def is_fpf_involution(self):
        return all(self(i) != i and self(self(i)) == i for i in self)

    def is_identity(self):
        return len(self.cycles) == 0 or max(map(len, self.cycles)) <= 1

    def is_fireworks(self):
        decr = []
        for a in self.oneline:
            if len(decr) == 0 or decr[-1][-1] < a:
                decr.append([a])
            else:
                decr[-1].append(a)
        return all(decr[i][0] < decr[i + 1][0] for i in range(len(decr) - 1))

    def __iter__(self):
        return self.oneline.__iter__()

    # Input is [i_1, i_2, ... , i_k], returns permutation (i_1 i_2 ... i_k)
    @classmethod
    def cycle(cls, cyc):
        return cls.cycles([cyc])

    @classmethod
    def cycles(cls, cyc, oneline=None):

        if len(cyc) == 0:
            return cls()

        if oneline is None:
            n = max(map(max, cyc))
            oneline = list(range(1, 1 + n))
            for c in cyc:
                oneline[c[len(c) - 1] - 1] = c[0]
                for i in range(len(c) - 1):
                    oneline[c[i] - 1] = c[i + 1]
        ans = cls()
        return cls(*oneline)

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
        for i in range(1, 1 + max(len(self.oneline), len(other.oneline))):
            if self(i) < other(i):
                return True
            if self(i) > other(i):
                return False
        return True

    def __eq__(self, other):
        if type(other) != Permutation:
            return False
        return len(self * (other**-1)) == 0

    def __ne__(self, other):
        return not (self == other)

    def strong_bruhat_less_equal(self, other):
        return self == other or self.strong_bruhat_less_than(other)

    def strong_bruhat_less_than(self, other):
        if self.length() >= other.length():
            return False
        des = self.get_descentset_R()
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

    def oneline_repr(self, n):
        assert n >= len(self.oneline)
        return ('' if n < 10 else ',').join([str(self(i)) for i in range(1, n + 1)])

    def cycle_repr(self):
        if len(self) == 0:
            return '1'
        EXCLUDE_SINGLE_CYCLES = True
        SPACE = ' '
        DELIM = ','

        s = ''
        for c in self.cycles:
            if not EXCLUDE_SINGLE_CYCLES or len(c) > 1:
                s += '(' + (DELIM + SPACE).join([str(x) for x in c]) + ')'
        return s

    def __hash__(self):
        return hash(tuple(self.oneline))

    def twisted_involution_length(self, n):
        w0 = Permutation.longest_element(n)
        y = w0 * self
        return w0.involution_length() - y.involution_length()

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
        if len(args) == 1 and type(args[0]) in [list, tuple]:
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

    def inverse_upward_k_pieri_tree(self, k=None):
        v, e = self.upward_k_pieri_tree(k)
        return {(xc[0].inverse(), xc[1]) for xc in v}, {((xc[0].inverse(), xc[1]), (yc[0].inverse(), yc[1]), t) for (xc, yc, t) in e}

    def upward_k_pieri_tree(self, k=None):
        k = (self.inverse() % self).number_two_cycles() if k is None else k
        vertices = set()
        edges = set()
        q = deque([
            ((self, 2**k), set(), {}, None, self.rank + 1)
        ])
        while q:
            vc, seen_a, a_counter, last_a, last_b = q.popleft()
            v, coeff = vc
            vertices.add(vc)
            for b in range(k + 1, last_b + 1):
                e = max([v(a) for a in range(k + 1, b) if v(a) < v(b)], default=None)
                for a in range(k, 0, -1):
                    if v(a) > v(b):
                        continue
                    if e is not None and v(a) < e:
                        continue
                    e = v(a)
                    if last_a is not None and b == last_b and last_a > a and a_counter[last_a] > 1:
                        continue

                    w = v * v.t_ij(a, b)
                    assert len(w) == len(v) + 1
                    
                    new_seen_a = seen_a.copy() | {a}
                    
                    new_coeff = 2**(k - len(new_seen_a)) * abs(coeff) // coeff
                    if last_a is not None and b == last_b and last_a > a:
                        new_coeff *= -1

                    new_a_counter = a_counter.copy()
                    new_a_counter[a] = new_a_counter.get(a, 0) + 1

                    wc = (w, new_coeff)
                    q.append((wc, new_seen_a, new_a_counter, a, b))
                    edges.add((vc, wc, (a, b)))
        return vertices, edges

    def inverse_downward_k_pieri_tree(self, k):
        v, e = self.downward_k_pieri_tree(k)
        return {(xc[0].inverse(), xc[1]) for xc in v}, {((xc[0].inverse(), xc[1]), (yc[0].inverse(), yc[1]), t) for (xc, yc, t) in e}

    def downward_k_pieri_tree(self, k):
        vertices = set()
        edges = set()
        q = deque([
            ((self, 2**k), set(), set(), None, k + 1)
        ])
        while q:
            vc, prohibited, seen_a, last_a, last_b = q.popleft()
            v, coeff = vc
            vertices.add(vc)
            for b in range(last_b, v.rank + 2):
                e = min([v(a) for a in range(k + 1, b) if v(a) > v(b)], default=None)
                for a in range(k, 0, -1):
                    if v(a) < v(b):
                        continue
                    if a in prohibited:
                        e = v(a) if e is None else min(e, v(a))
                        continue
                    if e is not None and v(a) > e:
                        continue
                    e = v(a)
                    
                    w = v * v.t_ij(a, b)
                    assert len(w) == len(v) - 1
                    
                    new_seen_a = seen_a.copy() | {a}
                    new_coeff = 2**(k - len(new_seen_a)) * abs(coeff) // coeff
                    new_prohibited = prohibited.copy()
                    if last_a is not None and b == last_b and a > last_a:
                        new_prohibited.add(a)
                        new_coeff *= -1

                    wc = (w, new_coeff)
                    q.append((wc, new_prohibited, new_seen_a, a, b))
                    edges.add((wc, vc, (a, b)))
        return vertices, edges

    def inverse_k_pieri_chains(self, k, p):
        for v, forced, prohibited, length in self.inverse().k_pieri_chains(k, p):
            yield v.inverse(), forced, prohibited, length

    def k_pieri_chains(self, k, p, lowerbound=1, z=None):
        if (self, k, lowerbound) not in K_PIERI_CHAINS:
            yield self, 0, 0, ()
            ans = [(self, 0, 0, ())]
            
            q = deque([
                (self * self.t_ij(a, b), ((a, b),), 1, 0, {a: 1})
                for a, b in self.k_bruhat_covers(k)
            ]) 
            while q:
                v, path, forced, prohibited, a_counter = q.popleft()
                a0, b0 = path[-1]

                if forced > p:
                    continue
                if a0 < lowerbound:
                    continue
                if z is not None and a0 >= z(a0) and b0 <= z(b0):
                    continue
                
                yield v, forced, prohibited, path
                ans.append((v, forced, prohibited, path))

                for a1, b1 in v.k_bruhat_covers(k):
                    if b0 < b1:
                        continue
                    if a_counter[a0] > 1 and b0 == b1 and a0 >= a1:
                        continue
                    
                    new_forced = forced + int(b0 == b1 and a0 >= a1)
                    new_prohibited = prohibited + int(a1 in a_counter)

                    new_a_counter = a_counter.copy()
                    new_a_counter[a1] = new_a_counter.get(a1, 0) + 1

                    new_path = v * v.t_ij(a1, b1), path + ((a1, b1),), new_forced, new_prohibited, new_a_counter
                    q.append(new_path)
            K_PIERI_CHAINS[self, k, lowerbound] = ans
        else:
            for a in K_PIERI_CHAINS[self, k, lowerbound]:
                yield a

    def k_bruhat_covers(self, k):
        if (self, k) not in K_BRUHAT_COVERS:
            ans = []
            for a in range(1, k + 1):
                for b in range(k + 1, self.rank + 2):
                    if self(a) < self(b) and not any(self(a) < self(e) < self(b) for e in range(a + 1, b)):
                        ans.append((a, b))
            K_BRUHAT_COVERS[self, k] = ans
        return K_BRUHAT_COVERS[self, k]
