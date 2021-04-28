from words import (
    involution_insert,
    fpf_insert,
    Word,
    get_fpf_involution_words
)
from permutations import Permutation
from keys import monomial_from_composition, symmetric_halves
import tests.test_keys as testkeys
import subprocess
import time


BASE_DIRECTORY = '/Users/emarberg/examples/crystals/'


class AbstractCrystalMixin:

    def __init__(self, rank, vertices, edges, weights, printer=str):
        self._rank = rank
        self._vertices = set(vertices)
        self.f_operators = {}
        self.e_operators = {}
        for (i, v, w) in edges:
            self.f_operators[(i, v)] = w
            self.e_operators[(i, w)] = v
        self.weights = weights.copy()
        self.printer = printer
        self.e_strings = {}
        self.f_strings = {}
        self.s_operators = {}

    def draw(self, extended=False):
        s = ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        #
        for x in self:
            s += ['    "%s";' % self.printer(x)]
        #
        for v in self:
            for i in self.extended_indices if extended else self.indices:
                w = self.f_operator(i, v)
                if w is not None:
                    s += ['    "%s" -> "%s" [label="%s"];' % (self.printer(v), self.printer(w), i)]
        s += ['}']
        s = '\n'.join(s)
        #
        filename = self.filename(time.time)
        dot_filename = BASE_DIRECTORY + 'abstract/' + 'dot/' + '%s.dot' % filename
        png_filename = BASE_DIRECTORY + 'abstract/' + 'png/' + '%s.png' % filename
        with open(dot_filename, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])

    def e_string(self, i, v):
        assert i in self.indices
        assert v in self.vertices
        if (i, v) not in self.e_strings:
            k = 0
            w = self.e_operator(i, v)
            while w is not None:
                k += 1
                w = self.e_operator(i, w)
            self.e_strings[(i, v)] = k
        return self.e_strings[(i, v)]

    def f_string(self, i, v):
        assert i in self.indices
        assert v in self.vertices
        if (i, v) not in self.f_strings:
            k = 0
            w = self.f_operator(i, v)
            while w is not None:
                k += 1
                w = self.f_operator(i, w)
            self.f_strings[(i, v)] = k
        return self.f_strings[(i, v)]

    def e_operator(self, i, v):
        raise NotImplementedError

    def f_operator(self, i, v):
        raise NotImplementedError

    def s_operator(self, i, b):
        assert 1 <= i < self.rank
        assert b in self.vertices
        if (i, b) not in self.s_operators:
            k = self.f_string(i, b) - self.e_string(i, b)
            ans = b
            if k >= 0:
                for _ in range(k):
                    ans = self.f_operator(i, ans)
            else:
                for _ in range(-k):
                    ans = self.e_operator(i, ans)
            self.s_operators[(i, b)] = ans
        return self.s_operators[(i, b)]

    def weight(self, v):
        return self.weights[v]

    def __iter__(self):
        return iter(self._vertices)

    def __len__(self):
        return len(self.vertices)

    @property
    def vertices(self):
        return self._vertices

    @property
    def indices(self):
        raise NotImplementedError

    @property
    def extended_indices(self):
        raise NotImplementedError

    @property
    def rank(self):
        return self._rank

    def is_highest_weight(self, v):
        return all(self.e_operator(i, v) is None for i in self.extended_indices)

    def get_highest_weights(self):
        return [(v, self.weight(v)) for v in self if self.is_highest_weight(v)]

    def group_highest_weights(self):
        ans = {}
        for v, mu in self.get_highest_weights():
            ans[mu] = ans.get(mu, []) + [v]
        return ans

    def get_highest_weight_multiplicities(self):
        return {k: len(v) for k, v in self.group_highest_weights().items()}

    @classmethod
    def tensor_rank(cls, b, c):
        assert b.rank == c.rank
        return b.rank

    @classmethod
    def tensor_vertices(cls, b, c):
        return [(x, y) for x in b for y in c]

    @classmethod
    def tensor_printer(cls, b, c):
        def printer(pair):
            return b.printer(pair[0]) + ' ' + c.printer(pair[1])
        return printer

    @classmethod
    def tensor_weights(cls, b, c):
        def add_weights(w1, w2):
            assert len(w1) == len(w2)
            return tuple(w1[i] + w2[i] for i in range(len(w1)))

        return {(x, y): add_weights(b.weight(x), c.weight(y)) for x in b for y in c}

    @classmethod
    def tensor_edges(cls, b, c):
        raise NotImplementedError

    def tensor(self, c):
        b = self
        cls = type(b)
        assert type(b) == type(c)
        assert b.rank == c.rank
        rank = cls.tensor_rank(b, c)
        vertices = cls.tensor_vertices(b, c)
        edges = cls.tensor_edges(b, c)
        weights = cls.tensor_weights(b, c)
        printer = cls.tensor_printer(b, c)
        return cls(rank, vertices, edges, weights, printer)

    def isomorphic_subcrystals(self, a, b):
        if a is None and b is None:
            return True
        if a is None or b is None:
            return False
        assert a in self.vertices
        assert b in self.vertices
        children_a = {}
        none_a = []
        children_b = {}
        none_b = []
        for i in self.extended_indices:
            x = self.f_operator(i, a)
            if x is None:
                none_a.append(i)
            else:
                children_a[x] = children_a.get(x, []) + [i]
            x = self.f_operator(i, b)
            if x is None:
                none_b.append(i)
            else:
                children_b[x] = children_b.get(x, []) + [i]
        if none_a != none_b:
            return False
        children_a_values = {tuple(v) for v in children_a.values()}
        children_b_values = {tuple(v) for v in children_b.values()}
        if children_a_values != children_b_values:
            return False
        seen = set()
        for i in self.extended_indices:
            x = self.f_operator(i, a)
            y = self.f_operator(i, b)
            if (x, y) not in seen:
                if not self.isomorphic_subcrystals(x, y):
                    return False
                seen.add((x, y))
        return True


class AbstractGLCrystal(AbstractCrystalMixin):

    @classmethod
    def f_operator_on_words(cls, i, word):
        cl, word = type(word), list(word)
        stack = []
        for index, a in reversed(list(enumerate(word))):
            if a == i:
                stack.append(index)
            elif a == i + 1 and stack:
                stack.pop()
        if stack:
            word[stack[0]] += 1
            return cl(word)

    @classmethod
    def e_operator_on_words(cls, i, word):
        cl, word = type(word), list(word)
        stack = []
        for index, a in enumerate(word):
            if a == i + 1:
                stack.append(index)
            elif a == i and stack:
                stack.pop()
        if stack:
            word[stack[0]] -= 1
            return cl(word)

    @property
    def indices(self):
        return list(range(1, self.rank))

    @property
    def extended_indices(self):
        return self.indices

    def e_operator(self, i, v):
        assert i in self.indices
        assert v in self.vertices
        return self.e_operators.get((i, v), None)

    def f_operator(self, i, v):
        assert i in self.indices
        assert v in self.vertices
        return self.f_operators.get((i, v), None)

    @classmethod
    def tensor_edges(cls, b, c):
        assert b.rank == c.rank
        edges = []
        for x in b:
            for y in c:
                for i in range(1, b.rank):
                    if b.e_string(i, x) < c.f_string(i, y):
                        yy = c.f_operator(i, y)
                        if yy is not None:
                            edges.append((i, (x, y), (x, yy)))
                    else:
                        xx = b.f_operator(i, x)
                        if xx is not None:
                            edges.append((i, (x, y), (xx, y)))
        return edges

    def filename(self, ts=None):
        return "gl(%s)_crystal.%s" % (self.rank, len(self))

    @classmethod
    def standard_object(cls, rank):
        vertices = list(range(1, rank + 1))
        edges = [(i, i, i + 1) for i in range(1, rank)]
        weights = {}
        for v in vertices:
            wt = rank * [0]
            wt[v - 1] = 1
            weights[v] = tuple(wt)
        return cls(rank, vertices, edges, weights)


class AbstractQCrystal(AbstractCrystalMixin):

    @classmethod
    def f_operator_on_words(cls, i, word):
        if i > 0:
            return AbstractGLCrystal.f_operator_on_words(i, word)
        elif i == -1:
            cl, word = type(word), list(word)
            ones = [j for j in range(len(word)) if word[j] == 1]
            twos = [j for j in range(len(word)) if word[j] == 2]
            if not ones or (ones and twos and twos[0] < ones[0]):
                return None
            else:
                word[ones[0]] = 2
                return cl(word)

    @classmethod
    def e_operator_on_words(cls, i, word):
        if i > 0:
            return AbstractGLCrystal.e_operator_on_words(i, word)
        elif i == -1:
            cl, word = type(word), list(word)
            ones = [j for j in range(len(word)) if word[j] == 1]
            twos = [j for j in range(len(word)) if word[j] == 2]
            if not twos or (ones and twos and ones[0] < twos[0]):
                return None
            else:
                word[twos[0]] = 1
                return cl(word)

    @property
    def indices(self):
        return [-1] + list(range(1, self.rank))

    @property
    def extended_indices(self):
        return [i for i in range(-self.rank + 1, self.rank) if i != 0]

    def e_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        if i < -1:
            x = self.e_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.e_operators[(i, v)] = x
        return self.e_operators.get((i, v), None)

    def f_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        if i < -1:
            x = self.f_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.f_operators[(i, v)] = x
        return self.f_operators.get((i, v), None)

    @classmethod
    def tensor_edges(cls, b, c):
        edges = AbstractGLCrystal.tensor_edges(b, c)
        for x in b:
            xweight = b.weight(x)
            for y in c:
                if xweight[0] == xweight[1] == 0:
                    yy = c.f_operator(-1, y)
                    if yy is not None:
                        edges.append((-1, (x, y), (x, yy)))
                else:
                    xx = b.f_operator(-1, x)
                    if xx is not None:
                        edges.append((-1, (x, y), (xx, y)))
        return edges

    def filename(self, ts=None):
        return "q(%s)_crystal.%s" % (self.rank, len(self))

    @classmethod
    def standard_object(cls, rank):
        assert rank >= 2
        vertices = list(range(1, rank + 1))
        edges = [(i, i, i + 1) for i in range(1, rank)] + [(-1, 1, 2)]
        weights = {}
        for v in vertices:
            wt = rank * [0]
            wt[v - 1] = 1
            weights[v] = tuple(wt)
        return cls(rank, vertices, edges, weights)


class AbstractPrimedQCrystal(AbstractCrystalMixin):

    @classmethod
    def f_operator_on_words(cls, i, word):
        cl, word = type(word), list(word)
        if i > 0:
            u = AbstractGLCrystal.f_operator_on_words(i, [abs(x) for x in word])
            if u:
                u = cl([u[j] if word[j] >= 0 else -u[j] for j in range(len(u))])
            return u
        elif i == 0:
            ones = [j for j in range(len(word)) if word[j] in [-1, 1]]
            if not ones or (ones and word[ones[0]] == -1):
                return None
            word[ones[0]] = -1
            return cl(word)
        elif i == -1:
            ones = [j for j in range(len(word)) if word[j] in [-1, 1]]
            twos = [j for j in range(len(word)) if word[j] in [-2, 2]]
            if not ones or (ones and twos and twos[0] < ones[0]):
                return None
            a = word[ones[0]]
            if len(ones) == 1 or word[ones[0]] == word[ones[1]]:
                word[ones[0]] = 2 if a > 0 else -2
            else:
                word[ones[0]] = -2 if a > 0 else 2
                word[ones[1]] = 1 if a > 0 else -1
            return cl(word)

    @classmethod
    def e_operator_on_words(cls, i, word):
        cl, word = type(word), list(word)
        if i > 0:
            u = AbstractGLCrystal.e_operator_on_words(i, [abs(x) for x in word])
            if u:
                u = cl([u[j] if word[j] >= 0 else -u[j] for j in range(len(u))])
            return u
        elif i == 0:
            ones = [j for j in range(len(word)) if word[j] in [-1, 1]]
            if not ones or (ones and word[ones[0]] == 1):
                return None
            word[ones[0]] = 1
            return cl(word)
        elif i == -1:
            ones = [j for j in range(len(word)) if word[j] in [-1, 1]]
            twos = [j for j in range(len(word)) if word[j] in [-2, 2]]
            if not twos or (ones and twos and ones[0] < twos[0]):
                return None
            a = word[twos[0]]
            if len(ones) == 0 or word[ones[0]] * word[twos[0]] > 0:
                word[twos[0]] = 1 if a > 0 else -1
            else:
                word[twos[0]] = 1 if a < 0 else -1
                word[ones[0]] = -1 if a < 0 else 1
            return cl(word)

    def group_highest_weights(self):
        ans = {}
        for v, mu in self.get_highest_weights():
            mu = mu[1:]
            ans[mu] = ans.get(mu, []) + [v]
        return ans

    @property
    def indices(self):
        return list(range(-1, self.rank))

    @property
    def extended_indices(self):
        n = self.rank
        return sorted(set(range(-n + 1, n)) | set(range(0, n * n, n)))

    def e_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        if i < -1:
            x = self.e_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.e_operators[(i, v)] = x
        elif i >= self.rank:
            j = i // self.rank
            x = self.e_operator(i - self.rank, self.s_operator(j, v))
            x = x if x is None else self.s_operator(j, x)
            self.e_operators[(i, v)] = x
        return self.e_operators.get((i, v), None)

    def f_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        if i < -1:
            x = self.f_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.f_operators[(i, v)] = x
        elif i >= self.rank:
            j = i // self.rank
            x = self.f_operator(i - self.rank, self.s_operator(j, v))
            x = x if x is None else self.s_operator(j, x)
            self.f_operators[(i, v)] = x
        return self.f_operators.get((i, v), None)

    @classmethod
    def tensor_edges(cls, b, c):
        edges = AbstractGLCrystal.tensor_edges(b, c)
        for x in b:
            if b.e_string(0, x) + b.f_string(0, x) == 0:
                for y in c:
                    yy = c.f_operator(0, y)
                    if yy is not None:
                        edges.append((0, (x, y), (x, yy)))
            else:
                for y in c:
                    xx = b.f_operator(0, x)
                    if xx is not None:
                        edges.append((0, (x, y), (xx, y)))

            xweight = b.weight(x)
            fx = b.f_operator(-1, x)
            for y in c:
                if xweight[1] == xweight[2] == 0:
                    xx = x
                    yy = c.f_operator(-1, y)
                elif fx is not None and b.e_string(0, fx) == b.f_string(0, fx) < b.f_string(0, x) == c.e_string(0, y):
                    xx = b.f_operator(-1, b.f_operator(0, x))
                    yy = c.e_operator(0, y)
                elif fx is not None and b.e_string(0, fx) == b.f_string(0, fx) < b.e_string(0, x) == c.f_string(0, y):
                    xx = b.f_operator(-1, b.e_operator(0, x))
                    yy = c.f_operator(0, y)
                else:
                    xx = b.f_operator(-1, x)
                    yy = y

                if xx is not None and yy is not None:
                    edges.append((-1, (x, y), (xx, yy)))

        return edges

    def filename(self, ts=None):
        return "primed_q(%s)_crystal.%s" % (self.rank, len(self))

    @classmethod
    def standard_object(cls, rank):
        assert rank >= 2
        vertices = [i for i in range(-rank, rank + 1) if i != 0]
        edges = \
            [(i, i, i + 1) for i in range(1, rank)] + [(-1, 1, 2)] + \
            [(i, -i, -i - 1) for i in range(1, rank)] + [(-1, -1, -2)] + \
            [(0, 1, -1)]
        weights = {}
        for v in vertices:
            wt = (rank + 1) * [0]
            wt[0] = 2 if v > 0 else 1
            wt[abs(v)] = 1
            weights[v] = tuple(wt)
        printer = lambda x: str(abs(x)) + ("'" if x < 0 else "")
        return cls(rank, vertices, edges, weights, printer)


class OrthogonalCrystalGenerator:

    DIRECTORY = BASE_DIRECTORY + 'orthogonal/'

    def print_zero_weight_space(self):
        for f in self.zero_weight_space:
            p, q = involution_insert(*f)
            print(p)
            print('')
            print(q)
            print('')
            print(f)
            print('\n')

    @property
    def _filename(self):
        mu = ''.join([str(i) for i in sorted(self.alpha) if i])
        w = ''.join([str(i) for i in self.word])
        rank = str(self.num_factors)
        size = str(len(self._factorizations))
        return 'rank%s_alpha%s_size%s_%s' % (rank, mu, size, w)

    @property
    def dot_filename(self):
        return self.DIRECTORY + 'dot/' + '%s.dot' % self._filename

    @property
    def png_filename(self):
        return self.DIRECTORY + 'png/' + '%s.png' % self._filename

    def node_label(self, i):
        # turning increasing factorizations into decreasing by conjugating by long element
        pre = ''
        if i == 0:
            a, b = symmetric_halves(self.alpha, pad=True)
            pre += 'a: ' + str(self.alpha) + '\nr: ' + str(a) + '\nc: ' + str(b) + '\n\n' + str(self.tableau) + '\n\n'

        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(self.rank - abs(j)) + ("'" if j < 0 else "") for j in w.elements]) for w in self[i]])
        return pre + top + '\n\n' + bottom

    @classmethod
    def is_factorization_compatible(cls, f):
        return all(len(a) == 0 or i < min([abs(j) for j in a]) for i, a in enumerate(f))

    def is_highlighted(self, x):
        f = tuple(tuple(self.rank - a for a in part) for part in self[x])
        return self.is_factorization_compatible(f)

    def highlighted_nodes(self):
        s = []
        for x in range(len(self)):
            if self.is_highlighted(x):
                s += ['    "%s" [fillcolor=white];' % self.node_label(x)]
        return s

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=spline;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        #
        s += self.highlighted_nodes()
        #
        s += ['    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), i) for (x, y, i) in self.edges]# if self.is_highlighted(x) or self.is_highlighted(y)]
        s += ['}']
        s = '\n'.join(s)
        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    def is_alpha_increasing(self):
        return all(self.alpha[i] >= self.alpha[i + 1] for i in range(len(self.alpha) - 1))

    @property
    def alpha(self):
        if self._alpha is None:
            qkey = 0
            for x in range(len(self)):
                if self.is_highlighted(x):
                    qkey += monomial_from_composition(self.weights[x])
            for i in range(10):
                try:
                    self._alpha = testkeys.decompose_q(qkey)
                    return self._alpha
                except:
                    qkey *= 2
            print(qkey // 1024)
            print(self)
            raise Exception
        else:
            return self._alpha

    @classmethod
    def all(cls, n, k, dominant=False):
        for w in Permutation.involutions(n):
            for cg in cls.from_permutation(n, w, k):
                print()
                print('orthogonal crystal generator for word', cg.word)
                print('edges:', len(cg.edges))
                if not cg.edges:
                    continue
                try:
                    #if dominant and not cg.is_alpha_increasing():
                    #    continue
                    if len(cg.alpha) > k:
                        continue
                    print('generating . . .')
                    cg.generate()
                    assert cg.is_connected()
                except:
                    pass
            print()
            yield cg

    @classmethod
    def from_permutation(cls, n, pi, num_factors):
        fac = [
            tuple([Word(*i) for i in tup])
            for w in pi.get_involution_words()
            for tup in cls.get_increasing_factorizations(w, num_factors)
        ]
        dictionary = {}
        for f in fac:
            t = involution_insert(*f)[0]
            dictionary[t] = dictionary.get(t, []) + [f]
        for t in dictionary:
            fac = sorted(dictionary[t], key=lambda f: tuple(len(a) for a in f), reverse=True)
            yield cls(tuple(pi(i) for i in range(1, n + 1)), fac, num_factors)

    def __init__(self, oneline, factorizations, num_factors):
        self.rank = len(oneline)
        self.num_factors = num_factors
        self._factorizations = factorizations
        self._weights = None
        self._ranks = None

        def adjust(a):
            return (self.rank - abs(a)) * (-1 if a < 0 else 1)

        self.word = tuple(adjust(a) for a in self.insertion_tableau(0).row_reading_word())
        self.tableau = testkeys.o_eg_insert(self.word)[0]
        self._alpha = None
        self.f = {}
        self.compute_edges()

    def __repr__(self):
        return 'Crystal generator for P^O(w) = %s with ell = %s' % (self.insertion_tableau(0), self.num_factors)

    @property
    def factorizations(self):
        return self._factorizations

    @property
    def highest_weight(self):
        return self.weights[0]

    @property
    def weights(self):
        if self._weights is None:
            self._weights = [tuple(len(a) for a in f) for f in self.factorizations]
        return self._weights

    @property
    def zero_weight_space(self):
        return [f for f in self.factorizations if all(len(p) == 1 for p in f)]

    @property
    def ranks(self):
        if self._ranks is None:
            self._ranks = [self.get_rank(i) for i in range(len(self))]
        return self._ranks

    def __len__(self):
        return len(self.factorizations)

    def __getitem__(self, i):
        return self.factorizations[i]

    def insertion_tableau(self, i):
        return involution_insert(*self[i])[0]

    def recording_tableau(self, i):
        return involution_insert(*self[i])[1]

    @property
    def components(self):
        def dfs(seed):
            seen = {seed}
            added = {seed}
            while added:
                connected = set()
                for s in added:
                    connected |= \
                        {y for (x, y, e) in self.edges if s == x} | \
                        {x for (x, y, e) in self.edges if s == y}
                added = (connected - seen)
                seen |= added
            for x in seen:
                yield x

        rest = set(range(len(self)))
        components = []
        while rest:
            seed = rest.pop()
            comp = set(dfs(seed))
            rest -= comp
            components.append(comp)
        return components

    def is_connected(self):
        return len(self.components) == 1

    @property
    def edges(self):
        return self._edges

    def compute_edges(self):
        def label(i):
            return (str(-i) + "'") if i < 0 else str(i)

        self.e = {}
        self._edges = []
        for x in range(len(self)):
            for i in range(1, self.num_factors):
                y = self.f_i(x, i)
                self.e[(y, i)] = x
                if y is not None:
                    self._edges.append((x, y, label(i)))

        for x in range(len(self)):
            for i in range(1, self.num_factors):
                y = self.f_i(x, -i)
                self.e[(y, -i)] = x
                if y is not None:
                    self._edges.append((x, y, label(-i)))

    def sigma(self, element_index, operator_index):
        if element_index is None:
            return None
        k = self.phi(element_index, operator_index) - self.epsilon(element_index, operator_index)
        if k > 0:
            for i in range(k):
                element_index = self.f_i(element_index, operator_index)
        elif k < 0:
            for i in range(-k):
                element_index = self.e_i(element_index, operator_index)
        return element_index

    def phi(self, element_index, operator_index):
        k = 0
        while True:
            element_index = self.f_i(element_index, operator_index)
            if element_index is None:
                break
            k += 1
        return k

    def epsilon(self, element_index, operator_index):
        k = 0
        while True:
            element_index = self.e_i(element_index, operator_index)
            if element_index is None:
                break
            k += 1
        return k

    def e_i(self, element_index, operator_index):
        return self.e.get((element_index, operator_index), None)

    def f_i(self, element_index, operator_index):
        if (element_index, operator_index) not in self.f:
            self.f[(element_index, operator_index)] = self.compute_f_i(element_index, operator_index)
        return self.f[(element_index, operator_index)]

    def compute_f_i(self, element_index, operator_index):
        if operator_index < -1:
            i = abs(operator_index) - 1
            element_index = self.sigma(element_index, i)
            element_index = self.sigma(element_index, i + 1)
            element_index = self.f_i(element_index, -i)
            element_index = self.sigma(element_index, i + 1)
            element_index = self.sigma(element_index, i)
            return element_index
        elif operator_index > 0:
            a = self[element_index][operator_index - 1].elements
            b = self[element_index][operator_index].elements
            while a and b:
                for j in range(len(a)):
                    if abs(a[j]) > abs(b[-1]):
                        a = a[:j] + a[j + 1:]
                        break
                b = b[:-1]
            if len(a) == 0:
                return None

            x = max([abs(i) for i in a])
            sgn = -1 if x not in a else 1

            a = set(self[element_index][operator_index - 1].elements) - {sgn * x}
            b = set(self[element_index][operator_index].elements)
            while x in b or -x in b:
                if -x in b:
                    a = (a - {x + 1}) | {-x - 1}
                    b = (b - {-x}) | {x}
                x += 1
            b |= {sgn * x}

            ans = self[element_index]
            words = (Word(*sorted(a, key=abs)), Word(*sorted(b, key=abs)))
            ans = ans[:operator_index - 1] + words + ans[operator_index + 1:]
        elif operator_index == -1:
            if self.num_factors < 2:
                return None
            a = self[element_index][0].elements
            b = self[element_index][1].elements
            if len(a) == 0:
                return None
            if len(b) > 0 and abs(a[0]) > abs(b[0]):
                return None
            x = a[0]
            a = sorted(set(a) - {x}, key=abs)
            if a and (x < 0 < a[0] or a[0] < 0 < x):
                a[0] *= -1
                x *= -1
            b = sorted(set(b) | {x}, key=abs)
            ans = self[element_index]
            words = (Word(*a), Word(*b))
            ans = words + ans[2:]
        elif operator_index == 0:
            if not any(self[element_index]):
                return None
            i = 0
            while not self[element_index][i]:
                i += 1
            a = list(self[element_index][i].elements)
            a[0] *= -1
            ans = self[element_index]
            ans = ans[:i] + (Word(*a),) + ans[i + 1:]
        # print(self[element_index], '--[', operator_index, ']-->', ans)
        js = [j for j in range(len(self)) if self[j] == ans]
        assert len(js) == 1
        return js[0]

    def get_weight(self, i):
        return self.weights[i]

    def get_rank(self, j):
        wt = self.get_weight(j)
        n = len(wt)
        assert n == self.num_factors
        wt = [self.highest_weight[i] - wt[i] for i in range(len(wt))]
        ans = 0
        for i in range(n):
            ans += (n - i) * wt[i]
        return -ans

    @classmethod
    def get_increasing_primed_factorizations(cls, w, k):
        sw = tuple((-2 * i - 1) if i < 0 else 2 * i for i in w)
        for t in cls.get_increasing_factorizations(sw, k):
            st = tuple(tuple((i + 1) // -2 if i % 2 != 0 else i // 2 for i in a) for a in t)
            yield st

    @classmethod
    def get_increasing_factorizations(cls, w, k):
        def is_increasing(x):
            return all(x[i] > x[i - 1] for i in range(1, len(x)))
        #
        if k == 0 and len(w) == 0:
            yield tuple()
        elif k == 0:
            return
        elif len(w) == 0:
            yield tuple(k * [tuple()])
        else:
            for i in range(len(w) + 1):
                if is_increasing(w[:i]):
                    for tup in cls.get_increasing_factorizations(w[i:], k - 1):
                        yield (tuple(w[:i]),) + tup
                else:
                    break


class SymplecticCrystalGenerator(OrthogonalCrystalGenerator):

    DIRECTORY = BASE_DIRECTORY + 'symplectic/'

    @property
    def alpha(self):
        if self._alpha is None:
            pkey = 0
            for x in range(len(self)):
                if self.is_highlighted(x):
                    pkey += monomial_from_composition(self.weights[x])
            self._alpha = testkeys.decompose_p(pkey)
        return self._alpha

    @classmethod
    def all(cls, n, k, dominant=False):
        for w in Permutation.fpf_involutions(n):
            for cg in cls.from_permutation(n, w, k):
                print()
                print('symplectic crystal generator for word', cg.word)
                print('edges:', len(cg.edges))
                if not cg.edges:
                    continue
                try:
                    if dominant and not cg.is_alpha_increasing():
                        continue
                    if len(cg.alpha) > k:
                        continue
                    print('generating . . .')
                    cg.generate()
                    assert cg.is_connected()
                except:
                    pass
            print()
            yield cg

    @classmethod
    def from_permutation(cls, n, pi, num_factors):
        fac = [
            tuple([Word(*i) for i in tup])
            for w in pi.get_fpf_involution_words()
            for tup in cls.get_increasing_factorizations(w, num_factors)
        ]
        dictionary = {}
        for f in fac:
            t = fpf_insert(*f)[0]
            dictionary[t] = dictionary.get(t, []) + [f]
        for t in dictionary:
            fac = sorted(dictionary[t], key=lambda f: tuple(len(a) for a in f), reverse=True)
            yield cls(tuple(pi(i) for i in range(1, n + 1)), fac, num_factors)

    def __init__(self, oneline, factorizations, num_factors):
        self.rank = len(oneline)
        self.num_factors = num_factors
        self._factorizations = factorizations
        self._weights = None
        self._ranks = None

        def adjust(a):
            return (self.rank - abs(a)) * (-1 if a < 0 else 1)

        self.word = tuple(adjust(a) for a in self.insertion_tableau(0).row_reading_word())
        self.tableau = testkeys.sp_eg_insert(self.word)[0]

        self._alpha = None
        self.f = {}
        self.compute_edges()

    def insertion_tableau(self, i):
        return fpf_insert(*self[i])[0]

    def recording_tableau(self, i):
        return fpf_insert(*self[i])[1]

    def print_zero_weight_space(self):
        for f in self.zero_weight_space:
            p, q = fpf_insert(*f)
            print(p)
            print('')
            print(q)
            print('')
            print(f)
            print('\n')

    def _get_words(self):
        return get_fpf_involution_words(self.oneline)

    def __repr__(self):
        return 'FPF crystal generator for w = %s with ell = %s' % (self.oneline, self.num_factors)

    def compute_f_i(self, element_index, operator_index):
        if operator_index > 0 or operator_index < -1:
            return super(SymplecticCrystalGenerator, self).compute_f_i(element_index, operator_index)
        assert operator_index == -1
        if self.num_factors < 2:
            return None
        a = self[element_index][0].elements
        b = self[element_index][1].elements
        if len(a) == 0:
            return None
        if len(b) > 0 and min(a) > min(b):
            return None
        if len(a) == 1:
            x = a[0]
            a = set(a) - {x}
            b = set(b) | {x}
        else:
            x, y = a[:2]
            if y % 2 == 0:
                a = set(a) - {x}
                b = set(b) | {x}
            elif y > x:
                a = set(a) - {y}
                b = set(b) | {y - 2}
        ans = self[element_index]
        words = (Word(*sorted(a)), Word(*sorted(b)))
        ans = words + ans[2:]
        js = [j for j in range(len(self)) if self[j] == ans]
        assert len(js) == 1
        return js[0]
