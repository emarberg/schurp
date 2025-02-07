from partitions import Partition
from words import (
    involution_insert,
    fpf_insert,
    Word,
    get_fpf_involution_words
)
from stable.tableaux import Tableau
import subprocess
import time
import collections
import itertools


BASE_DIRECTORY = '/Users/emarberg/examples/crystals/'


class AbstractCrystalMixin:

    def is_even(self, v):
        return all(self.e_string(i, v) % 2 == 0 for i in self.indices)

    def character(self):
        from stable.polynomials import Polynomial
        ans = Polynomial()
        for b in self:
            ans += Polynomial.from_tuple((0,) + self.weight(b))
        return ans

    def squared_crystal(self):
        edges = {(i, v, self.f_operator(i, self.f_operator(i, v))) for i in self.indices for v in self if self.f_operator(i, v) is not None and self.f_operator(i, self.f_operator(i, v)) is not None}
        return self.__class__(self.rank, self.vertices, edges, self.weights, self.printer)

    def as_gl_crystal(self):
        edges = {(i, v, self.f_operator(i, v)) for i in range(1, self.rank) for v in self if self.f_operator(i, v) is not None}
        return AbstractGLCrystal(self.rank, self.vertices, edges, self.weights, self.printer)

    def __init__(self, rank, vertices, edges, weights, printer=str, provided_operators=None):
        self._rank = rank
        self._vertices = set(vertices)
        self.f_operators = {}
        self.e_operators = {}
        for (i, v, w) in edges:
            assert (i, v) not in self.f_operators
            assert (i, w) not in self.e_operators
            self.f_operators[(i, v)] = w
            self.e_operators[(i, w)] = v
        self.weights = {v: weights[v] for v in self._vertices}
        self.printer = printer
        self.e_strings = {}
        self.f_strings = {}
        self.s_operators = {}
        self.provided_operators = self.indices if provided_operators is None else provided_operators

    def truncate(self, subset):
        edges = []
        for i in self.extended_indices:
            for v in subset:
                w = self.f_operator(i, v)
                if w is not None and w in subset:
                    edges.append((i, v, w))
        return self.__class__(self.rank, subset, edges, self.weights, printer=self.printer, provided_operators=self.extended_indices)

    def index_printer(self, i, primed=False):
        return str(i)

    def draw(self, extended=(), highlighted_nodes=(), tex=False, exclude=False, node_width=0.0, node_height=0.0, neato=False):
        if type(extended) == bool:
            extended_operators = self.extended_indices if extended else ()
        else:
            extended_operators = extended

        def tex_tuple(n):
            parts = []
            for tup in n:
                if len(tup) == 0:
                    parts += ['\\cdot']
                else:
                    letters = [str(i) if i > 0 else str(-i) + "'" for i in tup]
                    parts += [''.join(letters)]
            return '$' + '\\hspace{0.0mm} /\\hspace{0.0mm} '.join(parts) + '$'

        def tex_string(n):
            try:
                return n.tex()
            except:
                return tex_tuple(n)

        printmap = {}

        def printer(n):
            if tex:
                if n not in printmap:
                    printmap[n] = str(len(printmap))
                return printmap[n]
            return self.printer(n)

        s = ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=true;']
        if not tex:
            s += ['    node [shape=box; fontname="courier"; style=filled];']
            # s += ['    node [shape=point,width=0.15,color=gray];']
        if tex:
            # s += ['    node [shape=box,style=filled,color=gray92];']
            s += ['    node [shape=box];']
            for x in self:
                if exclude and x not in highlighted_nodes:
                    continue
                ts = tex_string(x)
                color = 'gray95'
                if x in highlighted_nodes:
                    ts = '$\\boxed{' + ts[1:-1] + '}$'
                    color = 'none'
                s += ['    "%s" [color="%s",margin="0.0",width="%s",height="%s",texlbl="%s"];' % (printer(x), color, node_width, node_height, ts)]
        else:
            for x in self:
                if exclude and x not in highlighted_nodes:
                    continue
                if x in highlighted_nodes:
                    s += ['    "%s" [fillcolor=white];' % printer(x)]
                else:
                    s += ['    "%s";' % printer(x)]
        #
        for v in self:
            for i in set(self.provided_operators) | set(extended_operators):
                w = self.f_operator(i, v)
                if w is None:
                    continue
                if not exclude or (v in highlighted_nodes and w in highlighted_nodes):
                    if tex:
                        cstr = "blue" if i in [-1, 1] else "red" if i == 2 else "teal" if i == 3 else "black"
                        style = "dotted" if i == 0 else "dashed" if i < 0 else "solid"
                        s += ['    "%s" -> "%s" [style="%s",color="%s"];' % (printer(v), printer(w), style, cstr)]
                    else:
                        cstr = "black"
                        # cstr = "blue" if i in [-1, 1] else "red" if i == 2 else "teal" if i == 3 else "black"
                        istr = self.index_printer(i)
                        s += ['    "%s" -> "%s" [label="%s",color="%s"];' % (printer(v), printer(w), istr, cstr)]
                #
                # if i < 0 and extended:
                #     w = self.fprime_operator(i, v)
                #     if w is not None and (not exclude or (v in highlighted_nodes and w in highlighted_nodes)):
                #         if tex:
                #             cstr = "blue" if i == -1 else "red" if i == -2 else "teal"
                #             style = "dashed"
                #             s += ['    "%s" -> "%s" [style="%s",color="%s"];' % (printer(v), printer(w), style, cstr)]
                #         else:
                #             istr = self.index_printer(i, primed=True)
                #             s += ['    "%s" -> "%s" [label="%s"];' % (printer(v), printer(w), istr)]
                #
        s += ['}']
        s = '\n'.join(s)
        #
        filename = self.filename(time.time)
        dot_filename = BASE_DIRECTORY + 'abstract/' + 'dot/' + '%s.dot' % filename
        png_filename = BASE_DIRECTORY + 'abstract/' + 'png/' + '%s.png' % filename
        with open(dot_filename, 'w') as f:
            f.write(s)
        subprocess.run(["neato" if neato else "dot", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])

        if tex:
            ps = subprocess.Popen(("neato" if neato else "dot") + " -Txdot " + dot_filename + " | dot2tex -tmath --nominsize --figonly --figpreamble=\"\\small\" > ~/test.tex", stdin=subprocess.PIPE, shell=True)
            ps.communicate()

    def e_string(self, i, v):
        assert i in self.extended_indices
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
        assert i in self.extended_indices
        assert v in self.vertices
        if (i, v) not in self.f_strings:
            k = 0
            w = self.f_operator(i, v)
            while w is not None:
                k += 1
                w = self.f_operator(i, w)
            self.f_strings[(i, v)] = k
        return self.f_strings[(i, v)]

    def capital_e_operator(self, i, v):
        return self.eplus_operator(i, v)

    def rectify(self, v, reduced_word=None):
        assert v in self.vertices
        n = self.rank
        if reduced_word is None:
            reduced_word = [i for j in range(1, n) for i in range(j, 0, -1)]
        ans = v
        for i in reduced_word:
            ans = self.capital_e_operator(i, ans)
        return ans

    def e_operator(self, i, v):
        assert i in self.indices
        if v is None:
            return None
        assert v in self.vertices
        return self.e_operators.get((i, v), None)

    def f_operator(self, i, v):
        assert i in self.indices
        if v is None:
            return None
        assert v in self.vertices
        return self.f_operators.get((i, v), None)

    def eprime_operator(self, i, v):
        if i >= 0:
            return self.e_operator(i, v)
        return self.star_operator(self.f_operator(-self.rank - i, self.star_operator(v)))

    def fprime_operator(self, i, v):
        if i >= 0:
            return self.f_operator(i, v)
        return self.star_operator(self.e_operator(-self.rank - i, self.star_operator(v)))

    def s_operator(self, i, b):
        assert 0 <= i < self.rank
        if b not in self.vertices:
            print(self.vertices)
            print(i)
            print(b)
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

    def eplus_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        while self.e_operator(i, v) is not None:
            v = self.e_operator(i, v)
        return v

    def fplus_operator(self, i, v):
        assert i in self.extended_indices
        assert v in self.vertices
        while self.f_operator(i, v) is not None:
            v = self.f_operator(i, v)
        return v

    def star_operator(self, b):
        if b is None:
            return None
        n = self.rank
        a = b
        for j in range(n):
            for i in range(1, n - j):
                a = self.s_operator(i, a)
        return a

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

    def is_lowest_weight(self, v):
        return all(self.fprime_operator(i, v) is None for i in self.extended_indices)

    def get_lowest_weights(self):
        return [(v, self.weight(v)) for v in self if self.is_lowest_weight(v)]

    def group_lowest_weights(self):
        ans = {}
        for v, mu in self.get_lowest_weights():
            ans[mu] = ans.get(mu, []) + [v]
        return ans

    def get_lowest_weight_multiplicities(self):
        return {k: len(v) for k, v in self.group_lowest_weights().items()}

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

    def naive_highest_weights(self):
        return [b for b in self if all(self.e_operator(i, b) is None for i in self.indices)]

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

    @classmethod
    def find_isomorphism(cls, crystal_a, crystal_b, highest_a=None, highest_b=None):
        assert cls == crystal_a.__class__ == crystal_b.__class__ and crystal_a.rank == crystal_b.rank and crystal_a.extended_indices == crystal_b.extended_indices

        if highest_a is None:
            highest_a = [a for a in crystal_a if crystal_a.is_highest_weight(a)]
        else:
            highest_a = [highest_a]

        if highest_b is None:
            highest_b = [b for b in crystal_b if crystal_b.is_highest_weight(b)]
        else:
            highest_b = [highest_b]

        assert len(highest_a) == 1 or len(highest_b) == 1
        if len(highest_a) != len(highest_b):
            return None

        highest_a = highest_a[0]
        highest_b = highest_b[0]

        if crystal_a.weight(highest_a) != crystal_b.weight(highest_b):
            return None

        forward = {}
        queue = collections.deque([(highest_a, highest_b)])
        while queue:
            a, b = queue.popleft()
            if a in forward:
                if forward[a] != b:
                    return None
            else:
                forward[a] = b
                for i in crystal_a.extended_indices:
                    fa = crystal_a.f_operator(i, a)
                    fb = crystal_b.f_operator(i, b)
                    if fa is None:
                        if fb is not None:
                            return None
                    else:
                        if fb is None:
                            return None
                        queue.append((fa, fb))
        if len(forward) == len(crystal_b):
            return forward


    @classmethod
    def isomorphic_highest_weight_crystals(cls, crystal_a, highest_a, crystal_b, highest_b=None):
        return cls.find_isomorphism(crystal_a, crystal_b, highest_a, highest_b) is not None

        # a, b = highest_a, highest_b
        # if crystal_a.weight(a) != crystal_b.weight(b):
        #     return False

        # children_a, children_b = {}, {}
        # none_a, none_b = [], []
        
        # for i in crystal_a.extended_indices:
        #     x = crystal_a.f_operator(i, a)
        #     if x is None:
        #         none_a.append(i)
        #     else:
        #         children_a[x] = children_a.get(x, []) + [i]

        #     x = crystal_b.f_operator(i, b)
        #     if x is None:
        #         none_b.append(i)
        #     else:
        #         children_b[x] = children_b.get(x, []) + [i]
        
        # if none_a != none_b:
        #     return False
        
        # children_a_values = {tuple(v) for v in children_a.values()}
        # children_b_values = {tuple(v) for v in children_b.values()}
        # if children_a_values != children_b_values:
        #     return False
        # seen = {(None, None)}
        # for i in crystal_a.extended_indices:
        #     x = crystal_a.f_operator(i, a)
        #     y = crystal_b.f_operator(i, b)
        #     if (x, y) not in seen:
        #         if not cls.isomorphic_highest_weight_crystals(crystal_a, x, crystal_b, y):
        #             return False
        #         seen.add((x, y))
        # return True

    def get_component(self, a):
        rank = self.rank
        indices = self.indices
        e = lambda x, i: self.e_operator(i, x)
        f = lambda x, i: self.f_operator(i, x)
        wt = lambda x: self.weight(x)
        return self.from_element(a, rank, indices, e, f, wt)

    def get_components(self):
        rank = self.rank
        indices = self.provided_operators
        e = lambda x, i: self.e_operator(i, x)
        f = lambda x, i: self.f_operator(i, x)
        wt = lambda x: self.weight(x)
        
        elements = set(self)
        ans = []
        while elements:
            a = elements.pop()
            comp = self.from_element(a, rank, indices, e, f, wt)
            ans.append(comp)
            elements -= set(ans[-1])
        return ans

    def is_connected(self):
        return len(self.get_component(next(iter(self)))) == len(self)

    @classmethod
    def from_element(cls, a, rank, indices, e, f, wt=None):
        if wt is None:
            assert len(a) == rank
            wt = lambda x: tuple(len(_) for _ in x)
        vertices = []
        edges = set()
        weights = {}
        add = {a}
        while add:
            new_add = set()
            for a in add:
                vertices.append(a)
                weights[a] = wt(a)
                new_add |= {(a, i, e(a, i), False) for i in indices} | {(a, i, f(a, i), True) for i in indices}
            add = set()
            for a, i, b, fdir in new_add:
                if b is not None:
                    if b not in weights:
                        add.add(b)
                    new_edge = (i, a, b) if fdir else (i, b, a)
                    edges.add(new_edge)
        edges = list(edges)
        return cls(rank, vertices, edges, weights)

    def is_stembridge(self, return_counterexample=False):
        def retvalue(x):
            return self.get_component(x) if return_counterexample else False

        n = self.rank
        for i in range(1, n):
            for j in range(1, n):
                if i == j:
                    continue

                e_i = lambda x: self.e_operator(i, x)
                e_j = lambda x: self.e_operator(j, x)
                estr_i = lambda x: self.e_string(i, x)
                estr_j = lambda x: self.e_string(j, x)

                f_i = lambda x: self.f_operator(i, x)
                f_j = lambda x: self.f_operator(j, x)
                fstr_i = lambda x: self.f_string(i, x)
                fstr_j = lambda x: self.f_string(j, x)

                for x in self:
                    # S0
                    if e_i(x) is None and estr_i(x) != 0:
                        print('S0'); return retvalue(x)

                    # S1
                    y = e_i(x)
                    if y is not None:
                        if estr_j(y) - estr_j(x) not in [0, 1]:
                            print('S1a'); return retvalue(x)
                        if abs(i - j) > 1 and estr_j(y) != estr_j(x):
                            print('S1b'); return retvalue(x)

                    # S2
                    if estr_i(x) > 0 and estr_j(e_i(x)) == estr_j(x) > 0:
                        if e_i(e_j(x)) is None or e_j(e_i(x)) is None:
                            print('S2a'); return retvalue(x)
                        if e_i(e_j(x)) != e_j(e_i(x)):
                            print('S2b'); return retvalue(x)
                        if fstr_i(e_j(x)) != fstr_i(x):
                            print('S2c'); return retvalue(x)

                    # S3
                    if e_i(x) is not None and e_j(x) is not None and estr_j(e_i(x)) == estr_j(x) + 1 > 1 and estr_i(e_j(x)) == estr_i(x) + 1 > 1:
                        if e_i(e_j(x)) is None or e_i(e_i(e_j(x))) is None or e_j(e_i(e_i(e_j(x)))) is None or e_j(e_i(x)) is None or e_j(e_j(e_i(x))) is None or e_i(e_j(e_j(e_i(x)))) is None:
                            print('S3a'); return retvalue(x)
                        if e_j(e_i(e_i(e_j(x)))) != e_i(e_j(e_j(e_i(x)))):
                            print('S3b'); return retvalue(x)
                        if fstr_i(x) != fstr_i(e_j(x)) or fstr_i(e_j(x)) != fstr_i(e_j(e_j(e_i(x)))):
                            print('S3c'); return retvalue(x)
                        if fstr_j(x) != fstr_j(e_i(x)) or fstr_j(e_i(x)) != fstr_j(e_i(e_i(e_j(x)))):
                            print('S3d'); return retvalue(x)

                    # S0'
                    if f_i(x) is None and fstr_i(x) != 0:
                        print('dS0'); return retvalue(x)

                    # S1'
                    y = f_i(x)
                    if y is not None:
                        if fstr_j(y) - fstr_j(x) not in [0, 1]:
                            print('dS1a'); return retvalue(x)
                        if abs(i - j) > 1 and fstr_j(y) != fstr_j(x):
                            print('dS1b'); return retvalue(x)

                    # S2'
                    if fstr_i(x) > 0 and fstr_j(f_i(x)) == fstr_j(x) > 0:
                        if f_i(f_j(x)) is None or f_j(f_i(x)) is None:
                            print('dS2a'); return retvalue(x)
                        if f_i(f_j(x)) != f_j(f_i(x)):
                            print('dS2b'); return retvalue(x)
                        if estr_i(f_j(x)) != estr_i(x):
                            print('dS2c'); return retvalue(x)

                    # S3'
                    if f_i(x) is not None and f_j(x) is not None and fstr_j(f_i(x)) == fstr_j(x) + 1 > 1 and fstr_i(f_j(x)) == fstr_i(x) + 1 > 1:
                        if f_i(f_j(x)) is None or f_i(f_i(f_j(x))) is None or f_j(f_i(f_i(f_j(x)))) is None or f_j(f_i(x)) is None or f_j(f_j(f_i(x))) is None or f_i(f_j(f_j(f_i(x)))) is None:
                            print('dS3a'); return retvalue(x)
                        if f_j(f_i(f_i(f_j(x)))) != f_i(f_j(f_j(f_i(x)))):
                            print('dS3b'); return retvalue(x)
                        if estr_i(x) != estr_i(f_j(x)) or estr_i(f_j(x)) != estr_i(f_j(f_j(f_i(x)))):
                            print('dS3c'); return retvalue(x)
                        if estr_j(x) != estr_j(f_i(x)) or estr_j(f_i(x)) != estr_j(f_i(f_i(f_j(x)))):
                            print('dS3d'); return retvalue(x)
        return True


class AbstractGLCrystal(AbstractCrystalMixin):

    def virtual_type_c(self):
        assert self.rank % 2 == 0
        n = (self.rank + 1) // 2
        vertices = self.vertices
        edges = []
        weights = self.weights
        for b in vertices:
            for i in range(1, n):
                fb = self.f_operator(i, b)
                fb = None if fb is None else self.f_operator(self.rank - i, fb)
                if fb is not None:
                    edges.append((i, b, fb))

            fb = self.f_operator(n, b)
            fb = None if fb is None else self.f_operator(n, fb)
            if fb is not None:
                edges.append((n, b, fb))

        return AbstractCCrystal(n, vertices, edges, weights)

    @classmethod
    def strict_polarizations(cls, mu, rank):
        def col(s):
            return tuple(i for (i, j) in s)

        def f_operator(index, s):
            iword = [(k, p[0]) for (k, p) in enumerate(s) if p[0] in [index, index + 1]]
            stack = []
            for k, i in iword:
                if not stack or not (stack[-1][1] == index + 1 and i == index):
                    stack.append((k, i))
                elif stack[-1][1] == index + 1 and i == index:
                    stack.pop()
            stack = [(k, i) for (k, i) in stack if i == index]
            if stack:
                (k, i) = stack[-1]
                j = s[k][1]
                s = set(s)
                if j == i + 1:
                    s.remove((i, j))
                    s.add((i + 1, i))
                elif (j, i + 1) not in s:
                    return
                else:
                    s.remove((i, j))
                    s.remove((j, i + 1))
                    s.add((j, i))
                    s.add((i + 1, j))
                return tuple(sorted(s, key=lambda p: (p[1], -p[0])))

        if type(mu) != Partition:
            mu = Partition(*mu)
        assert mu.is_symmetric()
        shape = mu.shape

        n = rank
        vertices = []
        weights = {}

        for v in range(2**len({(i, j) for (i, j) in shape if i > j})):
            s = set()
            for (i, j) in shape:
                if i > j:
                    s.add((i, j) if v % 2 else (j, i))
                    v = v // 2
            s = tuple(sorted(s, key=lambda p: (p[1], -p[0])))
            vertices.append(s)
            
            w = n*[0]
            for (i, j) in s:
                w[i - 1] += 1
            weights[s] = tuple(w)

        edges = []
        for s in vertices:
            for i in range(1, n):
                t = f_operator(i, s)
                if t is not None:
                    assert t in weights
                    edges += [(i, s, t)]
        return cls(rank, vertices, edges, weights)

    def dual(self):
        cls = type(self)
        n = self.rank
        vertices = self.vertices
        edges = [(n - i, v, self.e_operator(i, v)) for i in range(1, n) for v in self if self.e_operator(i, v) is not None]
        weights = {v: tuple(reversed(self.weights[v])) for v in self}
        printer = self.printer
        return cls(n, vertices, edges, weights, printer)

    def reversal(self):
        n = self.rank
        vertices = self.vertices
        edges = [(i, self.f_operator(i, v), v) for i in range(1, n) for v in vertices if self.f_operator(i, v) is not None]
        weights = {v: tuple(-w for w in self.weight(v)) for v in vertices}
        return self.__class__(n, vertices, edges, weights)

    @classmethod
    def sqrtcrystal_of_words(cls, length, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.all(n, (length,), setvalued=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(1, n):
                u = t.sqrt_f_operator(i)
                if u is not None:
                    assert u.sqrt_e_operator(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def sqrtcrystal_from_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.semistandard(n, mu, setvalued=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(1, n):
                u = t.sqrt_f_operator(i)
                if u is not None:
                    assert u.sqrt_e_operator(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def standard_sqrtcrystal(cls, rank):
        n = rank
        ints = tuple(range(1, n + 1))
        vertices = [t for k in range(1, n + 1) for t in itertools.combinations(ints, k)]
        edges = []
        weights = {}
        for t in vertices:
            weights[t] = tuple(1 if i in t else 0 for i in range(1, n +1))
            for i in range(1, n):
                if i in t and i + 1 not in t:
                    u = tuple(sorted(t + (i + 1,)))
                    edges += [(i, t, u)]
                elif i in t and i + 1 in t:
                    u = tuple(_ for _ in t if _ != i)
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_partition(cls, mu, rank):
        return cls.semistandard_tableaux_from_partition(mu, rank)

    @classmethod
    def semistandard_tableaux_from_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.semistandard(n, mu):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(1, n):
                u = cls.f_operator_on_semistandard_tableaux(i, t)
                if u is not None:
                    assert cls.e_operator_on_semistandard_tableaux(i, u) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def setvalued_semistandard_tableaux_from_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.semistandard(n, mu, setvalued=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(1, n):
                u = t.f_operator_on_setvalued_semistandard_tableaux(i)
                if u is not None:
                    assert u.e_operator_on_setvalued_semistandard_tableaux(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_permutation(cls, z, n, increasing):
        rank = n
        vertices = []
        edges = []
        weights = {}
        if increasing:
            for w in z.get_reduced_words():
                for f1 in Word.increasing_factorizations(w, n):
                    f1 = tuple(_.tuple() for _ in f1)
                    vertices += [f1]
                    weights[f1] = tuple(len(_) for _ in f1)
                    for i in list(range(1, n)):
                        f2 = Word.incr_crystal_f(f1, i)
                        if f2 is not None:
                            edges += [(i, f1, f2)]
        else:
            m = z.rank
            c = cls.from_permutation(z.star(m), n, True)

            def unstar(t):
                return tuple(tuple(m - i for i in _) for _ in t)

            vertices = [unstar(t) for t in c.vertices]
            edges = [(i, unstar(x), unstar(c.f_operators[(i, x)])) for (i, x,) in c.f_operators]
            weights = {unstar(t): c.weights[t] for t in c.weights}
        return cls(rank, vertices, edges, weights)

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

    @classmethod
    def f_operator_on_semistandard_tableaux(cls, i, tab):
        word = tuple(tab.row_reading_word())
        ans = cls.f_operator_on_words(i, word)
        return None if ans is None else tab.from_row_reading_word(tuple(ans))

    @classmethod
    def e_operator_on_semistandard_tableaux(cls, i, tab):
        word = tuple(tab.row_reading_word())
        ans = cls.e_operator_on_words(i, word)
        return None if ans is None else tab.from_row_reading_word(tuple(ans))

    @property
    def indices(self):
        return list(range(1, self.rank))

    @property
    def extended_indices(self):
        return self.indices

    def filename(self, ts=None):
        return "gl%s_crystal.%s" % (self.rank, len(self))

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


class AbstractCCrystal(AbstractCrystalMixin):

    @property
    def indices(self):
        return list(range(1, self.rank + 1))

    @property
    def extended_indices(self):
        return self.indices

    def filename(self, ts=None):
        return "c%s_crystal.%s" % (self.rank, len(self))

    @classmethod
    def standard_object(cls, rank):
        vertices = [i + 1 for i in range(rank)] + [-(rank - i) for i in range(rank)]
        edges = [(i, i, i + 1) for i in range(1, rank)] + [(rank, rank, -rank)] + [(rank - 1 - i, -(rank - i), -(rank - i) + 1) for i in range(rank - 1)]
        weights = {}
        for v in vertices:
            wt = rank * [0]
            wt[abs(v) - 1] = 1 if v > 0 else -1
            weights[v] = tuple(wt)
        return cls(rank, vertices, edges, weights)


class SuperGLCrystal(AbstractCrystalMixin):

    @property
    def rank_m(self):
        return self.rank[0]

    @property
    def rank_n(self):
        return self.rank[1]

    @property
    def indices(self):
        return list(range(1 - self.rank_m, self.rank_n))

    @property
    def extended_indices(self):
        return self.indices

    def filename(self, ts=None):
        return "gl%s|%s_crystal.%s" % (self.rank_m, self.rank_n, len(self))

    @classmethod
    def f_operator_on_words(cls, i, word):
        reverse = lambda w: None if w is None else tuple(reversed(w))
        if i > 0:
            return AbstractGLCrystal.f_operator_on_words(i, word)
        elif i < 0:
            return reverse(AbstractGLCrystal.f_operator_on_words(i - 1, reverse(word)))
        elif i == 0:
            cl, word = type(word), list(word)
            ones = [j for j in range(len(word)) if word[j] == -1]
            twos = [j for j in range(len(word)) if word[j] == 1]
            if not ones or (ones and twos and twos[0] < ones[0]):
                return None
            else:
                word[ones[0]] = 1
                return cl(word)
        raise Exception

    @classmethod
    def e_operator_on_words(cls, i, word):
        reverse = lambda w: None if w is None else tuple(reversed(w))
        if i > 0:
            return AbstractGLCrystal.e_operator_on_words(i, word)
        elif i < 0:
            return reverse(AbstractGLCrystal.e_operator_on_words(i - 1, reverse(word)))
        elif i == 0:
            cl, word = type(word), list(word)
            ones = [j for j in range(len(word)) if word[j] == -1]
            twos = [j for j in range(len(word)) if word[j] == 1]
            if not twos or (ones and twos and ones[0] < twos[0]):
                return None
            else:
                word[twos[0]] = -1
                return cl(word)
        raise Exception

    @classmethod
    def standard_object(cls, rank_m, rank_n=None):
        if rank_n is None and type(rank_m) in [tuple, list] and len(rank_m) == 2:
            rank_m, rank_n = rank_m
        assert rank_n is not None
        rank = (rank_m, rank_n)

        vertices = [i for i in range(-rank_m, rank_n + 1) if i != 0]
        edges = [(-i, -i - 1, -i) for i in range(1, rank_m)]
        edges += [(i, i, i + 1) for i in range(1, rank_n)]
        edges += [(0, -1, 1)]
        weights = {}
        for v in vertices:
            wt = (rank_m + rank_n) * [0]
            if v > 0:
                wt[rank_m + v - 1] = 1
            else:
                wt[rank_m + v] = 1
            weights[v] = tuple(wt)
        return cls(rank, vertices, edges, weights)

    @classmethod
    def sqrtcrystal_from_partition(cls, mu, rank_m, rank_n=None):
        if rank_n is None and type(rank_m) in [tuple, list] and len(rank_m) == 2:
            rank_m, rank_n = rank_m
        assert rank_n is not None
        rank = (rank_m, rank_n)

        vertices = []
        edges = []
        weights = {}
        for t in Tableau.semistandard_super(rank_m, rank_n, mu, setvalued=True):
            vertices.append(t)
            weights[t] = t.superweight(rank_m, rank_n)
            for i in range(-rank_m + 1, rank_n):
                u = t.sqrt_super_f_operator(i)
                if u is not None:
                    assert u.sqrt_super_e_operator(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def standard_sqrtcrystal(cls, rank_m, rank_n):
        def unadjust(s):
            return tuple(x if x > 0 else x + 1 for x in s)

        def adjust(s):
            return tuple(x if x > 0 else x - 1 for x in s)

        rank = (rank_m, rank_n)
        ints = tuple(range(-rank_m + 1, rank_n + 1))
        vertices = [adjust(t) for k in range(1, rank_m + rank_n + 1) for t in itertools.combinations(ints, k)]
        edges = []
        weights = {}
        for t in vertices:
            wt = (rank_m + rank_n) * [0,]
            for i in unadjust(t):
                wt[i + rank_m - 1] = 1
            weights[t] = tuple(wt)
            for i in range(-rank_m + 1, rank_n):
                if i < 0:
                    if i - 1 in t and i not in t:
                        u = tuple(sorted(t + (i,)))
                        edges += [(i, t, u)]
                    elif i - 1 in t and i in t:
                        u = tuple(_ for _ in t if _ != i - 1)
                        edges += [(i, t, u)]
                elif i > 0:
                    if i in t and i + 1 not in t:
                        u = tuple(sorted(t + (i + 1,)))
                        edges += [(i, t, u)]
                    elif i in t and i + 1 in t:
                        u = tuple(_ for _ in t if _ != i)
                        edges += [(i, t, u)]
                elif i == 0:
                    if -1 in t and 1 not in t:
                        u = tuple(sorted(t + (1,)))
                        edges += [(0, t, u)]
                    elif -1 in t and 1 in t:
                        u = tuple(_ for _ in t if _ != -1)
                        edges += [(0, t, u)]
        return cls(rank, vertices, edges, weights)

    def character(self):
        from stable.polynomials import Polynomial, X, Y
        ans = Polynomial()
        for b in self:
            w = self.weight(b)
            mon = X(0)**0
            for i in range(1, self.rank_m + 1):
                mon *= Y(-i)**w[self.rank_m - i]
            for i in range(1, self.rank_n + 1):
                mon *= X(i)**w[self.rank_m + i - 1]
            ans += mon
        return ans

    @classmethod
    def tensor_edges(cls, b, c):
        edges = []
        for x in b:
            for y in c:
                for i in b.indices:
                    xx = x
                    yy = y
                    if i < 0:
                        if b.f_string(i, x) > c.e_string(i, y):
                            xx = b.f_operator(i, x)
                        else:
                            yy = c.f_operator(i, y)
                    elif i > 0:
                        if b.e_string(i, x) < c.f_string(i, y):
                            yy = c.f_operator(i, y)
                        else:
                            xx = b.f_operator(i, x)
                    elif i == 0:
                        if b.weight(x)[b.rank_m - 1] == b.weight(x)[b.rank_m] == 0:
                            yy = c.f_operator(i, y)
                        else:
                            xx = b.f_operator(i, x)
                    else:
                        raise Exception
                    if xx is not None and yy is not None:
                        edges.append((i, (x, y), (xx, yy)))
        return edges


class AbstractQCrystal(AbstractCrystalMixin):

    @classmethod
    def sqrtcrystal_of_words(cls, length, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.all(n, (length,), setvalued=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else [])  + list(range(1, n)):
                u = t.sqrt_f_operator(i)
                if u is not None:
                    assert u.sqrt_e_operator(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def sqrtcrystal_from_strict_partition(cls, mu, rank):
        return cls.decomposition_sqrtcrystal_from_strict_partition(mu, rank)

    @classmethod
    def standard_sqrtcrystal(cls, rank):
        n = rank
        ints = tuple(range(1, n + 1))
        vertices = [t for k in range(1, n + 1) for t in itertools.combinations(ints, k)]
        edges = []
        weights = {}
        for t in vertices:
            weights[t] = tuple(1 if i in t else 0 for i in range(1, n +1))
            for i in range(-1, n):
                if i == 0:
                    continue
                j = i
                i = abs(i)
                if i in t and i + 1 not in t:
                    u = tuple(sorted(t + (i + 1,)))
                    edges += [(j, t, u)]
                elif i in t and i + 1 in t:
                    u = tuple(_ for _ in t if _ != i)
                    edges += [(j, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def decomposition_sqrtcrystal_from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.setvalued_decomposition_tableaux(n, mu):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                u = t.sqrt_decomposition_f_operator(i)
                if u is not None:
                    assert u.is_decomposition_tableau()
                    assert u.sqrt_decomposition_e_operator(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def squared_setvalued_decomposition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.setvalued_decomposition_tableaux(n, mu):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                u = t.sqrt_decomposition_f_operator(i)
                if u is not None:
                    u = u.sqrt_decomposition_f_operator(i)
                    if u is not None:
                        assert u.is_decomposition_tableau()
                        assert u.sqrt_decomposition_e_operator(i).sqrt_decomposition_e_operator(i) == t
                        edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def setvalued_decomposition_tableaux_from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        for t in Tableau.setvalued_decomposition_tableaux(n, mu):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                u = t.f_operator_on_setvalued_decomposition_tableaux(i)
                if u is not None:
                    assert u.is_decomposition_tableau()
                    assert u.e_operator_on_setvalued_decomposition_tableaux(i) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def decomposition_tableaux_from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        from tableaux import Tableau
        for t in Tableau.get_semistandard_decomposition(n, mu):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                u = cls.f_operator_on_decomposition_tableaux(i, t)
                if u is not None:
                    assert cls.e_operator_on_decomposition_tableaux(i, u) == t
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        from tableaux import Tableau
        for t in Tableau.get_semistandard_shifted(mu, n, diagonal_primes=False):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                u = t.shifted_crystal_f(i)
                if u is not None:
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_strict_partition_quotient(cls, mu, rank):
        def unprime(t):
            mapping = {}
            seen = set()
            for (i, j) in sorted(t, key=lambda ij: (-ij[0], ij[1])):
                a = abs(t[(i, j)])
                if a not in seen:
                    mapping[(i, j)] = a
                    seen.add(a)
                else:
                    mapping[(i, j)] = t[(i, j)]
            return t.__class__(mapping)

        n = rank
        vertices = []
        edges = set()
        weights = {}
        provided_operators = list(range(1, n))
        c = cls.from_strict_partition(mu, n)
        for t in c.vertices:
            vertices += [unprime(t)]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in provided_operators:
                u = c.f_operator(i, t)
                if u is not None:
                    edges.add((i, unprime(t), unprime(u)))
        return cls(rank, vertices, list(edges), weights, provided_operators=provided_operators)

    @classmethod
    def from_involution(cls, z, n, increasing):
        rank = n
        vertices = []
        edges = []
        weights = {}
        if increasing:
            for w in z.get_involution_words():
                for f1 in Word.increasing_factorizations(w, n):
                    f1 = tuple(_.tuple() for _ in f1)
                    vertices += [f1]
                    weights[f1] = tuple(len(_) for _ in f1)
                    for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                        f2 = Word.incr_crystal_f(f1, i)
                        if f2 is not None:
                            edges += [(i, f1, f2)]
        else:
            m = z.rank
            c = cls.from_involution(z.star(m), n, True)

            def unstar(t):
                return tuple(tuple(m - i for i in _) for _ in t)

            vertices = [unstar(t) for t in c.vertices]
            edges = [(i, unstar(x), unstar(c.f_operators[(i, x)])) for (i, x,) in c.f_operators]
            weights = {unstar(t): c.weights[t] for t in c.weights}
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_fpf_factorization(cls, a, n, increasing):
        m = max([0] + [i for _ in a for i in _]) + 1
        m += int(m % 2 != 0)

        def star(b):
            return tuple(tuple(m - i for i in _) for _ in b) if b is not None else b

        def f(a, i):
            if increasing:
                return Word.incr_crystal_f(a, i if i != -1 else 'FPF')
            else:
                return star(Word.incr_crystal_f(star(a), i if i != -1 else 'FPF'))

        def e(a, i):
            if increasing:
                return Word.incr_crystal_e(a, i if i != -1 else 'FPF')
            else:
                return star(Word.incr_crystal_e(star(a), i if i != -1 else 'FPF'))

        indices = ([-1] if n >= 2 else []) + list(range(1, n))
        return cls.from_element(a, n, indices, e, f)

    @classmethod
    def from_fpf_involution(cls, z, n, increasing):
        rank = n
        vertices = []
        edges = []
        weights = {}
        if increasing:
            for w in z.get_fpf_involution_words():
                for f1 in Word.increasing_factorizations(w, n):
                    f1 = tuple(_.tuple() for _ in f1)
                    vertices += [f1]
                    weights[f1] = tuple(len(_) for _ in f1)
                    for i in ([-1] if n >= 2 else []) + list(range(1, n)):
                        f2 = Word.incr_crystal_f(f1, i if i != -1 else 'FPF')
                        if f2 is not None:
                            edges += [(i, f1, f2)]
        else:
            m = z.rank
            c = cls.from_fpf_involution(z.star(m), n, True)

            def unstar(t):
                return tuple(tuple(m - i for i in _) for _ in t)

            vertices = [unstar(t) for t in c.vertices]
            edges = [(i, unstar(x), unstar(c.f_operators[(i, x)])) for (i, x,) in c.f_operators]
            weights = {unstar(t): c.weights[t] for t in c.weights}
        return cls(rank, vertices, edges, weights)

    def index_printer(self, i, primed=False):
        if not primed:
            return str(i)
        else:
            assert -(self.rank) < i < 0
            return str(i) + "\'"

    @classmethod
    def s_operator_on_words(cls, j, w):
        assert j > 0
        if w is None:
            return None
        k = len([a for a in w if a == j]) - len([a for a in w if a == j + 1])
        for _ in range(abs(k)):
            w = cls.f_operator_on_words(j, w) if k >= 0 else cls.e_operator_on_words(j, w)
        return w

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
        elif i < -1:
            return cls.s_operator_on_words(-i - 1, cls.s_operator_on_words(-i, cls.f_operator_on_words(i + 1, cls.s_operator_on_words(-i, cls.s_operator_on_words(-i - 1, word)))))

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
        elif i < -1:
            return cls.s_operator_on_words(-i - 1, cls.s_operator_on_words(-i, cls.e_operator_on_words(i + 1, cls.s_operator_on_words(-i, cls.s_operator_on_words(-i - 1, word)))))

    @classmethod
    def f_operator_on_decomposition_tableaux(cls, i, tab):
        word = tuple(reversed(tab.row_reading_word()))
        ans = cls.f_operator_on_words(i, word)
        return None if ans is None else tab.decomposition_tableau_from_row_reading_word(tuple(reversed(ans)))

    @classmethod
    def e_operator_on_decomposition_tableaux(cls, i, tab):
        word = tuple(reversed(tab.row_reading_word()))
        ans = cls.e_operator_on_words(i, word)
        return None if ans is None else tab.decomposition_tableau_from_row_reading_word(tuple(reversed(ans)))

    @property
    def indices(self):
        return ([-1] if self.rank >= 2 else []) + list(range(1, self.rank))

    @property
    def extended_indices(self):
        return [i for i in range(-self.rank + 1, self.rank) if i != 0]

    def e_operator(self, i, v):
        assert i in self.extended_indices
        if v is None:
            return None
        assert v in self.vertices
        if i in self.provided_operators or (i, v) in self.e_operators:
            return self.e_operators.get((i, v), None)
        if i < -1:
            x = self.e_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.e_operators[(i, v)] = x
            return x
        raise Exception

    def f_operator(self, i, v):
        assert i in self.extended_indices
        if v is None:
            return None
        assert v in self.vertices
        if i in self.provided_operators or (i, v) in self.f_operators:
            return self.f_operators.get((i, v), None)
        if i < -1:
            x = self.f_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.f_operators[(i, v)] = x
            return x
        raise Exception

    @classmethod
    def tensor_edges(cls, b, c):
        edges = AbstractCrystalMixin.tensor_edges(b, c)
        for x in b:
            ex = b.e_operator(-1, x)
            fx = b.f_operator(-1, x)
            for y in c:
                if ex is None and fx is None:
                    xx = x
                    yy = c.f_operator(-1, y)
                else:
                    xx = fx
                    yy = y

                if xx is not None and yy is not None:
                   edges.append((-1, (x, y), (xx, yy)))
            # xweight = b.weight(x)
            # for y in c:
            #     if xweight[0] == xweight[1] == 0:
            #         yy = c.f_operator(-1, y)
            #         if yy is not None:
            #             edges.append((-1, (x, y), (x, yy)))
            #     else:
            #         xx = b.f_operator(-1, x)
            #         if xx is not None:
            #             edges.append((-1, (x, y), (xx, y)))
        return edges

    def filename(self, ts=None):
        return "q%s_crystal.%s" % (self.rank, len(self))

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
    def standard_sqrtcrystal(cls, rank):
        n = rank
        ints = set(range(-n, n + 1)) - {0}
        vertices = [tuple(sorted(t)) for k in range(1, 2 * n + 1) for t in itertools.combinations(ints, k)]

        gp_one = AbstractQCrystal.sqrtcrystal_of_words(1, rank)
        gp_two = AbstractQCrystal.sqrtcrystal_of_words(2, rank)

        edges = []
        weights = {}
        for t in vertices:
            w = [0 for _ in range(n)]
            for i in t:
                w[abs(i) - 1] += 1
            weights[t] = tuple(w)

            if 1 in t and -1 not in t:
                u = tuple(sorted(t + (-1,)))
                edges += [(0, t, u)]
            elif 1 in t and -1 in t:
                u = tuple(_ for _ in t if _ != 1)
                edges += [(0, t, u)]

            a = tuple(-i for i in t if i < 0)
            b = tuple(i for i in t if i > 0)

            for i in range(-1, n):
                if i == 0:
                    continue

                if a and not b:
                    tab = Tableau.from_rows([[a]])
                    uab = gp_one.f_operator(i, tab)
                    u = tuple(sorted([-x for x in uab.get(1, 1, unpack=False)])) if uab else None
                elif b and not a:
                    tab = Tableau.from_rows([[b]])
                    uab = gp_one.f_operator(i, tab)
                    u = tuple(sorted([x for x in uab.get(1, 1, unpack=False)])) if uab else None
                else:
                    tab = Tableau.from_rows([[a, b]])
                    uab = gp_two.f_operator(i, tab)
                    u = tuple(sorted([-x for x in uab.get(1, 1, unpack=False)] + [x for x in uab.get(1, 2, unpack=False)])) if uab else None

                if u:
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def f_operator_on_decomposition_tableaux(cls, i, tab):
        word = tuple(reversed(tab.row_reading_word()))
        ans = cls.f_operator_on_words(i, word)
        return None if ans is None else tab.decomposition_tableau_from_row_reading_word(tuple(reversed(ans)))

    @classmethod
    def e_operator_on_decomposition_tableaux(cls, i, tab):
        word = tuple(reversed(tab.row_reading_word()))
        ans = cls.e_operator_on_words(i, word)
        return None if ans is None else tab.decomposition_tableau_from_row_reading_word(tuple(reversed(ans)))

    @classmethod
    def decomposition_tableaux_from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        from tableaux import Tableau
        for t in Tableau.get_semistandard_decomposition(n, mu, primed=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(-1 if n >= 2 else 0, n):
                u = cls.f_operator_on_decomposition_tableaux(i, t)
                if u is not None:
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_inv_factorization(cls, a, n, increasing):
        m = max([0] + [abs(i) for _ in a for i in _]) + 1

        def invert(i):
            return m - i if i > 0 else -(m + i)

        def star(b):
            return tuple(tuple(invert(i) for i in _) for _ in b) if b is not None else b

        def f(a, i):
            if increasing:
                return Word.incr_crystal_f(a, i )
            else:
                return star(Word.incr_crystal_f(star(a), i))

        def e(a, i):
            if increasing:
                return Word.incr_crystal_e(a, i)
            else:
                return star(Word.incr_crystal_e(star(a), i))

        indices = ([-1] if n >= 2 else []) + list(range(0, n))
        return cls.from_element(a, n, indices, e, f)

    @classmethod
    def from_involution(cls, z, n, increasing):
        rank = n
        vertices = []
        edges = []
        weights = {}
        if increasing:
            for w in z.get_primed_involution_words():
                for f1 in Word.increasing_factorizations(w, n):
                    f1 = tuple(_.tuple() for _ in f1)
                    vertices += [f1]
                    weights[f1] = tuple(len(_) for _ in f1)
                    for i in ([-1] if n >= 2 else []) + list(range(0, n)):
                        f2 = Word.incr_crystal_f(f1, i)
                        if f2 is not None:
                            edges += [(i, f1, f2)]
        else:
            m = z.rank
            c = cls.from_involution(z.star(m), n, True)

            def invert(i):
                return m - i if i > 0 else -(m + i)

            def unstar(t):
                return tuple(tuple(invert(i) for i in _) for _ in t)

            vertices = [unstar(t) for t in c.vertices]
            edges = [(i, unstar(x), unstar(c.f_operators[(i, x)])) for (i, x,) in c.f_operators]
            weights = {unstar(t): c.weights[t] for t in c.weights}
        return cls(rank, vertices, edges, weights)

    @classmethod
    def from_strict_partition(cls, mu, rank):
        n = rank
        vertices = []
        edges = []
        weights = {}
        from tableaux import Tableau
        for t in Tableau.get_semistandard_shifted(mu, n, diagonal_primes=True):
            vertices += [t]
            wt = t.weight()
            weights[t] = wt + (n - len(wt)) * (0,)
            for i in range(-1 if n >= 2 else 0, n):
                u = t.shifted_crystal_f(i)
                if u is not None:
                    edges += [(i, t, u)]
        return cls(rank, vertices, edges, weights)

    def starb_operator(self, b):
        if b is None:
            return None
        n = self.rank
        a = b
        for j in range(n - 1, -1, -1):
            for i in range(j, 0, -1):
                a = self.s_operator(i, a)
            a = self.s_operator(0, a)
            for i in range(1, j + 1):
                a = self.s_operator(i, a)
        return a

    def star_plus_operator(self, b):
        if b is None:
            return None
        n = self.rank
        a = b
        for j in range(n):
            for i in range(0, n - j):
                a = self.s_operator(i, a)
        return a

    def index_printer(self, i, primed=False):
        if not primed:
            return str(i) if i < self.rank else ('0[' + str(i // self.rank) + ']')
        else:
            assert -(self.rank) < i < 0
            return str(self.rank + i) + "\'"

    @classmethod
    def s_operator_on_words(cls, j, w):
        assert j > 0
        if w is None:
            return None
        k = len([a for a in w if abs(a) == j]) - len([a for a in w if abs(a) == j + 1])
        for _ in range(abs(k)):
            w = cls.f_operator_on_words(j, w) if k >= 0 else cls.e_operator_on_words(j, w)
        return w

    @classmethod
    def f_operator_on_words(cls, i, word):
        if word is None:
            return None
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
        elif i < -1:
            return cls.s_operator_on_words(-i - 1, cls.s_operator_on_words(-i, cls.f_operator_on_words(i + 1, cls.s_operator_on_words(-i, cls.s_operator_on_words(-i - 1, word)))))

    @classmethod
    def e_operator_on_words(cls, i, word):
        if word is None:
            return None
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
        elif i < -1:
            return cls.s_operator_on_words(-i - 1, cls.s_operator_on_words(-i, cls.e_operator_on_words(i + 1, cls.s_operator_on_words(-i, cls.s_operator_on_words(-i - 1, word)))))


    def group_highest_weights(self):
        ans = {}
        for v, mu in self.get_highest_weights():
            mu = mu[1:]
            ans[mu] = ans.get(mu, []) + [v]
        return ans

    def group_lowest_weights(self):
        ans = {}
        for v, mu in self.get_lowest_weights():
            mu = mu[1:]
            ans[mu] = ans.get(mu, []) + [v]
        return ans

    @property
    def indices(self):
        return ([-1] if self.rank >= 2 else []) + list(range(self.rank))

    @property
    def extended_indices(self):
        n = self.rank
        return sorted(set(range(-n + 1, n)) | set(range(0, n * n, n)))

    def e_operator(self, i, v):
        assert i in self.extended_indices
        if v is None:
            return None
        assert v in self.vertices
        if i in self.provided_operators or (i, v) in self.e_operators:
            return self.e_operators.get((i, v), None)
        if i < -1:
            x = self.e_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.e_operators[(i, v)] = x
            return x
        elif i >= self.rank:
            j = i // self.rank
            x = self.e_operator(i - self.rank, self.s_operator(j, v))
            x = x if x is None else self.s_operator(j, x)
            self.e_operators[(i, v)] = x
            return x
        raise Exception

    def f_operator(self, i, v):
        assert i in self.extended_indices
        if v is None:
            return None
        assert v in self.vertices
        if i in self.provided_operators or (i, v) in self.f_operators:
            return self.f_operators.get((i, v), None)
        if i < -1:
            x = self.f_operator(i + 1, self.s_operator(-i, self.s_operator(-i - 1, v)))
            x = x if x is None else self.s_operator(-i - 1, self.s_operator(-i, x))
            self.f_operators[(i, v)] = x
            return x
        elif i >= self.rank:
            j = i // self.rank
            x = self.f_operator(i - self.rank, self.s_operator(j, v))
            x = x if x is None else self.s_operator(j, x)
            self.f_operators[(i, v)] = x
            return x
        raise Exception

    @classmethod
    def tensor_edges(cls, b, c):
        edges = AbstractCrystalMixin.tensor_edges(b, c)
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

            # xweight = b.weight(x)
            ex = b.e_operator(-1, x)
            fx = b.f_operator(-1, x)

            ffx = b.f_operator(-1, b.f_operator(0, x))
            effx = b.e_operator(0, ffx)
            fffx = b.f_operator(0, ffx)

            fex = b.f_operator(-1, b.e_operator(0, x))
            efex = b.e_operator(0, fex)
            ffex = b.f_operator(0, fex)

            for y in c:
                ey = c.e_operator(0, y)
                fy = c.f_operator(0, y)

                # if xweight[1] == xweight[2] == 0:
                #     assert ex is None and fx is None
                #     xx = x
                #     yy = c.f_operator(-1, y)
                # elif fx is not None and b.e_string(0, fx) == b.f_string(0, fx) < b.f_string(0, x) == c.e_string(0, y):
                #     assert ffx is not None and ey is not None and effx is None and fffx is None
                #     xx = b.f_operator(-1, b.f_operator(0, x))
                #     yy = c.e_operator(0, y)
                # elif fx is not None and b.e_string(0, fx) == b.f_string(0, fx) < b.e_string(0, x) == c.f_string(0, y):
                #     assert fex is not None and fy is not None and efex is None and ffex is None
                #     xx = b.f_operator(-1, b.e_operator(0, x))
                #     yy = c.f_operator(0, y)
                # else:
                #     xx = b.f_operator(-1, x)
                #     yy = y
                if ex is None and fx is None:
                    xx = x
                    yy = c.f_operator(-1, y)
                elif ffx is not None and ey is not None and effx is None and fffx is None:
                    xx = ffx
                    yy = ey
                elif fex is not None and fy is not None and efex is None and ffex is None:
                    xx = fex
                    yy = fy
                else:
                    xx = fx
                    yy = y

                if xx is not None and yy is not None:
                   edges.append((-1, (x, y), (xx, yy)))

        return edges

    def filename(self, ts=None):
        return "primed_q%s_crystal.%s" % (self.rank, len(self))

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
            wt = rank * [0]
            wt[abs(v) - 1] = 1
            #wt = (rank + 1) * [0]
            #wt[0] = 2 if v > 0 else 1
            #wt[abs(v)] = 1
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
        from keys import symmetric_halves
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
        from keys import monomial_from_composition
        if self._alpha is None:
            qkey = 0
            for x in range(len(self)):
                if self.is_highlighted(x):
                    qkey += monomial_from_composition(self.weights[x])
            for i in range(10):
                try:
                    from tests.test_keys import decompose_q
                    self._alpha = decompose_q(qkey)
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
        from permutations import Permutation
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
        from tests.test_keys import o_eg_insert
        self.tableau = o_eg_insert(self.word)[0]
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
        from keys import monomial_from_composition
        from tests.test_keys import decompose_p
        if self._alpha is None:
            pkey = 0
            for x in range(len(self)):
                if self.is_highlighted(x):
                    pkey += monomial_from_composition(self.weights[x])
            self._alpha = decompose_p(pkey)
        return self._alpha

    @classmethod
    def all(cls, n, k, dominant=False):
        from permutations import Permutation
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
        from tests.test_keys import sp_eg_insert
        self.tableau = sp_eg_insert(self.word)[0]

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
