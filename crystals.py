from words import (
    involution_insert,
    fpf_insert,
    alt_fpf_insert,
    HopfPermutation,
    reduce_oneline,
    reduce_fpf,
    Word,
    get_involution_words,
    get_fpf_involution_words
)
import subprocess


class ShiftedCrystalGenerator:

    DIRECTORY = '/Users/emarberg/Dropbox/shifted_crystals/examples/inv/'

    @classmethod
    def test_insertion_tableaux(cls, n, k):
        for i, w in enumerate(HopfPermutation.involutions(n)):
            cg = ShiftedCrystalGenerator(w.oneline, k)
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
        w = ''.join([str(i) for i in self.oneline])
        return 'size%s_%s_%s' % (str(len(self._factorizations)), w, str(self.num_factors))

    @property
    def dot_filename(self):
        return self.DIRECTORY + 'dot/' + '%s.dot' % self._filename

    @property
    def png_filename(self):
        return self.DIRECTORY + 'png/' + '%s.png' % self._filename

    def node_label(self, i):
        pre = str(self.insertion_tableau(i))
        top = str(self.recording_tableau(i))
        bottom = '/'.join([''.join([str(j) for j in w.elements]) for w in self[i]])
        return pre + '\n\n' + top + '\n\n' + bottom

    def write_dotfile(self):
        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=line;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        s += ['    "%s" -> "%s" [label="%s"];' % (self.node_label(x), self.node_label(y), i) for (x, y, i) in self.edges]
        s += ['}']
        s = '\n'.join(s)

        with open(self.dot_filename, 'w') as f:
            f.write(s)

    def generate(self):
        self.write_dotfile()
        subprocess.run(["dot", "-Tpng", self.dot_filename, "-o", self.png_filename])

    @classmethod
    def all(cls, n, k):
        for w in HopfPermutation.involutions(n):
            if w.oneline != reduce_oneline(w.oneline):
                continue
            cg = ShiftedCrystalGenerator(w.oneline, k)
            if cg.edges:
                cg.generate()

    def __init__(self, oneline, k):
        self.oneline = oneline
        self.words = self._get_words()
        self.num_factors = k
        self._factorizations = None
        self._weights = None
        self._ranks = None
        self._edges = None

    def _get_words(self):
        return get_involution_words(self.oneline)

    def __repr__(self):
        return 'Crystal generator for w = %s with ell = %s' % (self.oneline, self.num_factors)

    @property
    def factorizations(self):
        if self._factorizations is None:
            fac = [
                tuple([Word(*i) for i in tup])
                for w in self.words
                for tup in self.get_increasing_factorizations(w, self.num_factors)
            ]
            self._factorizations = sorted(fac, key=lambda f: tuple(len(a) for a in f), reverse=True)
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


class FPFCrystalGenerator(ShiftedCrystalGenerator):

    DIRECTORY = '/Users/emarberg/Dropbox/shifted_crystals/examples/fpf/'

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
            return super(FPFCrystalGenerator, self).f_i(element_index, operator_index)
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


class AltCrystalGenerator(FPFCrystalGenerator):

    DIRECTORY = '/Users/emarberg/Dropbox/shifted_crystals/examples/alt/'

    def insertion_tableau(self, i):
        return alt_fpf_insert(*self[i])[0]

    def recording_tableau(self, i):
        return alt_fpf_insert(*self[i])[1]
