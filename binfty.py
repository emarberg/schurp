from crystals import AbstractGLCrystal
from stable.tableaux import Tableau
import subprocess
import time
import collections
import itertools


BASE_DIRECTORY = '/Users/emarberg/examples/crystals/'


class InfiniteCrystal:

    @classmethod
    def is_isomorphic(cls, b, b_vertices, c, c_vertices, ignore_strings=False, verbose=False):
        assert b.indices == c.indices
        assert b.generator in b_vertices
        assert c.generator in c_vertices
        if len(b_vertices) != len(c_vertices):
            if verbose:
                print('b size:', len(b_vertices))
                print('c size:', len(c_vertices))
            #return False
        q = collections.deque([(b.generator, c.generator)])
        seen = {None: None}
        while q:
            x, y = q.popleft()
            if verbose:
                print('b element:', x)
                print('  weight =', b.weight(x))
                print('c element:', y)
                print('  weight =', c.weight(y))
                print()
            if x is None and y is not None:
                return False
            if x is not None and y is None:
                return False
            if x in seen:
                if y != seen[x]:
                    return False
                else:
                    continue
            else:
                seen[x] = y
            if b.weight(x) != c.weight(y):
                return False
            for i in b.indices:
                if not ignore_strings:
                    if b.f_string(i, x) != c.f_string(i, y):
                        return False
                    if b.e_string(i, x) != c.e_string(i, y):
                        return False
                
                fx, fy = b.f_operator(i, x), c.f_operator(i, y)
                ex, ey = b.e_operator(i, x), c.e_operator(i, y)

                fx = fx if fx in b_vertices else None
                fy = fy if fy in c_vertices else None

                ex = ex if ex in b_vertices else None
                ey = ey if ey in c_vertices else None

                q.append((fx, fy))
                q.append((ex, ey))
        return True

    @classmethod
    def tensor(cls, *args):
        return cls._tensor(1, *args)

    @classmethod
    def sqrt_tensor(cls, *args):
        return cls._tensor(2, *args)

    @classmethod
    def _tensor(cls, multiplier, *args):
        assert len(args) >= 1
        assert all(b.indices == args[0].indices for b in args)

        if len(args) == 1:
            return args[0]
        
        if len(args) > 2:
            return cls._tensor(multiplier, cls._tensor(multiplier, *args[:-1]), args[-1])

        b, c = args[0], args[1]
        generator = (b.generator, c.generator)
        indices = args[0].indices

        def e_operators(i, pair):
            x, y = pair

            ex = b.e_operator(i, x)
            ey = c.e_operator(i, y)

            phi = c.f_string(i, y)
            eps = b.e_string(i, x)
            
            if eps is None or (phi is not None and phi >= eps):
                return None if ey is None else (x, ey)
            else:
                return None if ex is None else (ex, y)

        def f_operators(i, pair):
            x, y = pair

            fx = b.f_operator(i, x)
            fy = c.f_operator(i, y)

            phi = c.f_string(i, y)
            eps = b.e_string(i, x)
            
            if phi is None or (eps is not None and phi <= eps):
                return None if fx is None else (fx, y)
            else:
                return None if fy is None else (x, fy)
                
        def e_strings(i, pair):
            x, y = pair
            d = multiplier * (c.weight(y)[i - 1] - c.weight(y)[i])
            x, y = b.e_string(i, x), c.e_string(i, y)
            return y if x is None else (x - d) if y is None else max(y, x - d)

        def f_strings(i, pair):
            x, y = pair
            d = multiplier * (b.weight(x)[i - 1] - b.weight(x)[i])
            x, y = b.f_string(i, x), c.f_string(i, y)
            return x if y is None else (y + d) if x is None else max(x, y + d)

        def weight_map(pair):
            x, y = pair
            x, y = b.weight(x), c.weight(y)
            assert len(x) == len(y)
            return tuple(x[i] + y[i] for i in range(len(x)))

        printer = lambda pair: b.printer(pair[0]) + ' ⊗ ' + c.printer(pair[1])

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map, printer)

    @classmethod
    def row(cls, j, n):
        generator = n * (0,)
        indices = list(range(1, n))

        def e_operators(i, x):
            if i < j:
                return None
            y = list(x)
            if i == j and y[i] > 0:
                y[i] -= 1
            elif y[i] > 0:
                y[i] -= 1
                y[i - 1] += 1
            else:
                return None
            return tuple(y)

        def f_operators(i, x):
            if i < j:
                return None
            y = list(x)
            if i == j:
                y[i] += 1
            elif y[i - 1] > 0:
                y[i - 1] -= 1
                y[i] += 1
            else:
                return None
            return tuple(y)

        e_strings = lambda i, x: None if i < j else x[i]
        f_strings = lambda i, x: None if i < j else -sum(x) if i == j else x[i - 1]

        def weight_map(x):
            ans = n * [0]
            for i in range(j, n):
                ans[j - 1] -= x[i]
                ans[i] += x[i]
            return tuple(ans)

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)
 
    @classmethod
    def sqrt_row(cls, j, n):
        generator = n * (0,)
        indices = list(range(1, n))

        def e_operators(i, x):
            y = list(x)
            if y[i] == 0:
                return None
            elif i == j or y[i - 1] > 0 != y[i] % 2:
                y[i] -= 1
            elif y[i] % 2 == 0:
                y[i - 1] += 2
                y[i] -= 1
            else:
                y[i - 1] += 1
            return tuple(y)

        def f_operators(i, x):
            y = list(x)
            if i < j or (i > j and y[i - 1] == 0):
                return None
            elif i == j or y[i] % 2 == 0:
                y[i] += 1
            elif y[i - 1] == 1:
                y[i - 1] -= 1
            else:
                y[i - 1] -= 2
                y[i] += 1
            return tuple(y)

        e_strings = lambda i, x: None if i < j else x[i] + (1 if i > j and x[i] % 2 != 0 and x[i - 1] == 0 else 0)
        
        def weight_map(x):
            ans = n * [0]
            for i in range(j, n):
                ans[j - 1] -= x[i] // 2
                ans[i] += (x[i] + 1) // 2
            return tuple(ans)

        def odd(x):
            return sum(map(lambda i: i % 2, x))

        #f_strings = lambda i, x: None if i < j else e_strings(i, x) + 2 * weight_map(x)[i - 1] - 2 * weight_map(x)[i]
        f_strings = lambda i, x: None if i < j else (-sum(x) + odd(x) - odd([x[j]])) if i == j else (x[i - 1] + x[i - 1] % 2 - (x[i] % 2) * (x[i - 1] > 0))

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)

    @classmethod
    def tableaux(cls, n):
        m = lambda t: None if t is None else t.marginalize(n)

        generator = m(Tableau())
        indices = list(range(1, n))

        e_operators = lambda i, x: m(AbstractGLCrystal.e_operator_on_semistandard_tableaux(i, x))
        f_operators = lambda i, x: m(AbstractGLCrystal.f_operator_on_semistandard_tableaux(i, x))

        def e_strings(i, x):
            q = []
            for w in x.row_reading_word():
                if w == i + 1:
                    q.append(w)
                elif w == i and q and q[-1] == i + 1:
                    q.pop()
            return len(q)

        def weight_map(x):
            tup = x.shape()
            ans = list(x.weight(n))
            for i in range(n):
                ans[i] -= tup[i] if i < len(tup) else 0
            return tuple(ans)

        def f_strings(i, x):
            eps = e_strings(i, x)
            wt = weight_map(x)
            return eps + wt[i - 1] - wt[i]

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)
 
    @classmethod
    def sqrt_tableaux(cls, n):
        m = lambda t: None if t is None else t.marginalize(n)

        generator = m(Tableau())
        indices = list(range(1, n))

        e_operators = lambda i, x: m(x.sqrt_e_operator(i))
        f_operators = lambda i, x: m(x.sqrt_f_operator(i))

        def e_strings(i, x):
            q = []
            for w in x.row_reading_word(setwise=True):
                if i + 1 in w and i not in w:
                    q.append(i + 1)
                    q.append(i + 1)
                elif i in w and i + 1 in w and not q:
                    q.append(i - 1)
                elif i in w and i + 1 not in w and q:
                    if q[-1] == i + 1:
                        q.pop()
                        q.pop()
                    elif q[-1] == i - 1:
                        q.pop()
            return len(q)

        def weight_map(x):
            tup = x.shape()
            ans = list(x.weight(n))
            for i in range(n):
                ans[i] -= tup[i] if i < len(tup) else 0
            return tuple(ans)

        def f_strings(i, x):
            eps = e_strings(i, x)
            wt = weight_map(x)
            return eps + 2 * (wt[i - 1] - wt[i])

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)
 
    @classmethod
    def elementary(cls, j, n):
        generator = 0
        indices = list(range(1, n))

        e_operators = lambda i, x: None if i != j else x + 1
        f_operators = lambda i, x: None if i != j else x - 1

        e_strings = lambda i, x: None if i != j else -x
        f_strings = lambda i, x: None if i != j else x

        weight_map = lambda x: tuple(0 if i not in [j - 1, j] else x if i == j - 1 else -x for i in range(n))

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)
 
    @classmethod
    def sqrt_elementary(cls, j, n):
        generator = 0
        indices = list(range(1, n))

        e_operators = lambda i, x: None if i not in [j, j + 1] else (x + 1) if i == j else ((x - 1) if x % 2 == 0 else None)
        f_operators = lambda i, x: None if i not in [j, j + 1] else (x - 1) if i == j else ((x + 1) if x % 2 != 0 else None)

        e_strings = lambda i, x: None if i not in [j, j + 1] else -x if i == j else x
        f_strings = lambda i, x: None if i not in [j, j + 1] else x if i == j else (0 if x % 2 == 0 else 1)

        weight_map = lambda x: tuple(0 if i not in [j - 1, j] else (x // 2 + int(x % 2)) if i == j - 1 else -(x // 2) for i in range(n))

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)

    @classmethod
    def odd_sqrt_elementary(cls, j, n):
        generator = 0
        indices = list(range(1, n))

        e_operators = lambda i, x: None if i not in [j, j - 1] else (x + 1) if i == j else ((x - 1) if x % 2 != 0 else None)
        f_operators = lambda i, x: None if i not in [j, j - 1] else (x - 1) if i == j else ((x + 1) if x % 2 == 0 else None)

        e_strings = lambda i, x: None if i not in [j, j - 1] else -x if i == j else (x + 1)
        f_strings = lambda i, x: None if i not in [j, j - 1] else x if i == j else (0 if x % 2 != 0 else 1)

        weight_map = lambda x: tuple(0 if i not in [j - 1, j] else (x // 2 + int(x % 2)) if i == j - 1 else -(x // 2) for i in range(n))

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)

    def __init__(self, generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map, printer=str):
        self.generator = generator
        self.indices = sorted(tuple(indices))
        self.e_operators = e_operators
        self.f_operators = f_operators
        self.e_strings = e_strings
        self.f_strings = f_strings
        self.weight_map = weight_map
        self.printer = printer

    def e_operator(self, i, x):
        return self.e_operators(i, x)

    def f_operator(self, i, x):
        return self.f_operators(i, x)

    def capital_e_operator(self, i, x):
        while x is not None:
            y = self.e_operator(i, x)
            if y is None:
                return x
            x = y

    def e_string(self, i, x):
        return self.e_strings(i, x)

    def f_string(self, i, x):
        return self.f_strings(i, x)

    def weight(self, x):
        return self.weight_map(x) if x is not None else None

    def filename(self, vertices, ts=None):
        n = str(len(self.indices))
        ell = str(len(vertices))
        ans = "gl%s_crystal.%s" % (n, ell)
        if ts is not None:
            ans += ".%s" % str(ts)
        return ans

    def draw_demazure(self, thresh, *args):
        vertices = self.demazure(thresh, *args)
        self.draw(vertices)

    def is_demazure_isomorphic(self, thresh, wa, wb):
        for (word_a, word_b) in [(wa, wb), (wb, wa)]:
            vertices = {(self.generator, 0)}
            for i in reversed(word_a):
                add = set()
                for (v, h) in vertices:
                    w = v
                    k = h
                    while w is not None and k <= thresh:
                        if k > h:
                            add.add((w, k))
                        w = self.f_operator(i, w)
                        k += 1
                vertices |= add

            def op(w):
                for i in reversed(word_b):
                    w = self.capital_e_operator(i, w)
                return w

            if not all(op(w) == self.generator for (w, h) in vertices):
                return False

        return True
        

    def demazure(self, thresh, *args):
        vertices = {(self.generator, 0)}
        for i in reversed(args):
            add = set()
            for (v, h) in vertices:
                w = v
                k = h
                while w is not None and k <= thresh:
                    if k > h:
                        add.add((w, k))
                    w = self.f_operator(i, w)
                    k += 1
            vertices |= add
        return {v for (v, h) in vertices if h <= thresh}

    def vertices(self, e_thresh, f_thresh):
        vertices = set()
        q = collections.deque([(self.generator, 0, 0)])
        while q:
            g, i, j = q.popleft()
            if g is not None and g not in vertices and i <= e_thresh and j <= f_thresh:
                vertices.add(g)
                for a in self.indices:
                    q.append((self.e_operator(a, g), i + 1, j))
                    q.append((self.f_operator(a, g), i, j + 1))
        return vertices

    def draw_thresh(self, e_thresh, f_thresh):
        vertices = self.vertices(e_thresh, f_thresh)
        self.draw(vertices)

    def draw(self, vertices):
        #vertices = {t for t in vertices if all(i > 1 or v == (1,) for i,j,v in t)}
        #vertices = {x for x in vertices if all(self.e_string(i, x) >= 0 for i in self.indices)}
        edges = []

        def printer(x):
            ans = [self.printer(x)]
            ans += [str(self.weight(x))]
            for i in self.indices:
                ans += ['string ' + str(i) + ': ' + str((self.e_string(i, x), self.f_string(i, x)))]
            return '\n'.join(ans)

        s = []
        s += ['digraph G {']
        s += ['    overlap=false;']
        s += ['    splines=true;']
        s += ['    node [shape=box; fontname="courier"; style=filled];']
        for v in vertices:
            s += ['    "%s";' % printer(v)]
        for v in vertices:
            for i in self.indices:
                w = self.f_operator(i, v)
                if w is None or w not in vertices:
                    continue
                s += ['    "%s" -> "%s" [label="%s",color="%s"];' % (printer(v), printer(w), str(i), "black")]
        s += ['}']
        s = '\n'.join(s)
        
        filename = self.filename(vertices)
        dot_filename = BASE_DIRECTORY + 'infinite/' + 'dot/' + '%s.dot' % filename
        png_filename = BASE_DIRECTORY + 'infinite/' + 'png/' + '%s.png' % filename
        with open(dot_filename, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])
