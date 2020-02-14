from words import (
    involution_insert,
    fpf_insert,
    HopfPermutation,
    reduce_oneline,
    reduce_fpf,
    Word,
    get_involution_words,
    get_fpf_involution_words
)
from permutations import Permutation
from keys import monomial_from_composition
import tests.test_keys as testkeys
import subprocess


BASE_DIRECTORY = '/Users/emarberg/Desktop/examples/crystals/'


class OrthogonalCrystalGenerator:

    DIRECTORY = BASE_DIRECTORY + 'orthogonal/'

    @classmethod
    def test_insertion_tableaux(cls, n, k):
        for i, w in enumerate(HopfPermutation.involutions(n)):
            cg = OrthogonalCrystalGenerator(w.oneline, k)
            shapes = [
                {cg.insertion_tableau(i) for i in comp}
                for comp in cg.components
            ]
            b1 = all(len(sh) == 1 for sh in shapes)
            b2 = all(len(s & t) == 0 for s in shapes for t in shapes if s != t)
            print(i, w, b1, b2, '(n = %s, k = %s)' % (n, k))
            assert b1 and b2

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
        pre = (str(self.alpha) + '\n\n' + str(self.tableau) + '\n\n') if i == 0 else ''

        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(self.rank - j) for j in w.elements]) for w in self[i]])
        return pre + top + '\n\n' + bottom

    @classmethod
    def is_factorization_compatible(cls, f):
        return all(len(a) == 0 or i < min(a) for i, a in enumerate(f))

    def is_highlighted(self, x):
        f = tuple(tuple(self.rank - a for a in part) for part in self[x])
        return self.is_factorization_compatible(f)

    def highlighted_nodes(self):
        s = []
        for x in range(len(self)):
            if self.is_highlighted(x):
                s += ['    "%s" [fillcolor=white];' % self.node_label(x)]
        return s

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

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        #
        s += self.highlighted_nodes()
        #
        s += ['    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), i) for (x, y, i) in self.edges if self.is_highlighted(x) or self.is_highlighted(y)]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    @classmethod
    def all(cls, n, k):
        for w in Permutation.involutions(n):
            #if w.involution_length() > 2:
            #    continue
            for cg in OrthogonalCrystalGenerator.from_permutation(n, w, k):
                print(cg.word)
                if not cg.edges:
                    continue
                try:
                    if not all(cg.alpha[i] >= cg.alpha[i + 1] for i in range(len(cg.alpha) - 1)):
                        continue
                    if len(cg.alpha) > k:
                        continue
                    cg.generate()
                except:
                    pass
            print()

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
            yield OrthogonalCrystalGenerator(tuple(pi(i) for i in range(1, n + 1)), fac, num_factors)

    def __init__(self, oneline, factorizations, num_factors):
        self.rank = len(oneline)
        self.num_factors = num_factors
        self._factorizations = factorizations
        self._weights = None
        self._ranks = None
        self._edges = None

        self.word = tuple(self.rank - a for a in self.insertion_tableau(0).row_reading_word())
        self.tableau = testkeys.o_eg_insert(self.word)[0]
        self._alpha = None

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

    @property
    def edges(self):
        if self._edges is None:
            self._edges = []
            for x in range(len(self)):
                for i in range(self.num_factors):
                    y = self.f_i(x, i)
                    if y is not None:
                        self._edges.append((x, y, i))
        return self._edges

    def f_i(self, element_index, operator_index):
        if operator_index > 0:
            a = self[element_index][operator_index - 1].elements
            b = self[element_index][operator_index].elements
            while a and b:
                for j in range(len(a)):
                    if a[j] > b[-1]:
                        a = a[:j] + a[j + 1:]
                        break
                b = b[:-1]
            if len(a) == 0:
                return None
            x = max(a)
            a = set(self[element_index][operator_index - 1].elements) - {x}
            b = set(self[element_index][operator_index].elements)
            while x in b:
                x += 1
            b |= {x}
            ans = self[element_index]
            words = (Word(*sorted(a)), Word(*sorted(b)))
            ans = ans[:operator_index - 1] + words + ans[operator_index + 1:]
        if operator_index == 0:
            if self.num_factors < 2:
                return None
            a = self[element_index][0].elements
            b = self[element_index][1].elements
            if len(a) == 0:
                return None
            if len(b) > 0 and min(a) > min(b):
                return None
            x = min(a)
            a = set(a) - {x}
            b = set(b) | {x}
            ans = self[element_index]
            words = (Word(*sorted(a)), Word(*sorted(b)))
            ans = words + ans[2:]
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

    @classmethod
    def test_insertion_tableaux(cls, n, k):
        for i, w in enumerate(HopfPermutation.fpf_involutions(n)):
            cg = cls(w.oneline, k)
            shapes = [
                {cg.insertion_tableau(i) for i in comp}
                for comp in cg.components
            ]
            b1 = all(len(sh) == 1 for sh in shapes)
            b2 = all(len(s & t) == 0 for s in shapes for t in shapes if s != t)
            print(i, w, b1, b2, '(n = %s, k = %s)' % (n, k))
            assert b1 and b2

    def insertion_tableau(self, i):
        return fpf_insert(*self[i])[0]

    def recording_tableau(self, i):
        return fpf_insert(*self[i])[1]

    def node_label(self, i):
        pre = str(self.insertion_tableau(i))
        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(j) for j in w.elements]) for w in self[i]])
        return pre + '\n\n' + top + '\n\n' + bottom

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

    @classmethod
    def all(cls, n, k):
        assert n % 2 == 0
        for w in HopfPermutation.fpf_involutions(n):
            if w.oneline != reduce_fpf(w.oneline):
                continue
            cg = cls(w.oneline, k)
            if cg.edges:
                cg.generate()

    def f_i(self, element_index, operator_index):
        if operator_index > 0:
            return super(SymplecticCrystalGenerator, self).f_i(element_index, operator_index)
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
