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
        s += ['    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), i) for (x, y, i) in self.edges if self.is_highlighted(x) or self.is_highlighted(y)]
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
