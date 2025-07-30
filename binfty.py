from crystals import AbstractGLCrystal
from stable.tableaux import Tableau
import subprocess
import time
import collections
import itertools


BASE_DIRECTORY = '/Users/emarberg/examples/crystals/'


class InfiniteCrystal:

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
    def tableaux(cls, n):
        m = lambda t: None if t is None else t.marginalize(n)

        generator = m(Tableau())
        indices = list(range(1, n))

        e_operators = lambda i, x: m(AbstractGLCrystal.e_operator_on_semistandard_tableaux(i, x))
        f_operators = lambda i, x: m(AbstractGLCrystal.f_operator_on_semistandard_tableaux(i, x))

        e_strings = lambda i, x: None
        f_strings = lambda i, x: None

        weight_map = lambda x: tuple(a - (n - i) for (i, a) in enumerate(x.weight(n)))

        return cls(generator, indices, e_operators, f_operators, e_strings, f_strings, weight_map)
 
    @classmethod
    def sqrt_tableaux(cls, n):
        m = lambda t: None if t is None else t.marginalize(n)

        generator = m(Tableau())
        indices = list(range(1, n))

        e_operators = lambda i, x: m(x.sqrt_e_operator(i))
        f_operators = lambda i, x: m(x.sqrt_f_operator(i))

        e_strings = lambda i, x: None
        f_strings = lambda i, x: None

        def weight_map(x):
            tup = x.shape()
            ans = list(x.weight(n))
            for i in range(n):
                ans[i] -= tup[i] if i < len(tup) else 0
            return tuple(ans)

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

        e_operators = lambda i, x: None if i != j else x + 1
        f_operators = lambda i, x: None if i != j else x - 1

        e_strings = lambda i, x: None if i != j else -x
        f_strings = lambda i, x: None if i != j else x

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

    def e_string(self, i, x):
        return self.e_strings(i, x)

    def f_string(self, i, x):
        return self.f_strings(i, x)

    def weight(self, x):
        return self.weight_map(x)

    def filename(self, ts=None):
        ans = "gl%s_crystal" % str(len(self.indices))
        if ts is not None:
            ans += ".%s" % str(ts)
        return ans

    def draw(self, e_thresh, f_thresh):
        vertices = set()
        q = collections.deque([(self.generator, 0, 0)])
        while q:
            g, i, j = q.popleft()
            if g is not None and g not in vertices and i <= e_thresh and j <= f_thresh:
                vertices.add(g)
                for a in self.indices:
                    q.append((self.e_operator(a, g), i + 1, j))
                    q.append((self.f_operator(a, g), i, j + 1))

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
        
        filename = self.filename()
        dot_filename = BASE_DIRECTORY + 'infinite/' + 'dot/' + '%s.dot' % filename
        png_filename = BASE_DIRECTORY + 'infinite/' + 'png/' + '%s.png' % filename
        with open(dot_filename, 'w') as f:
            f.write(s)
        subprocess.run(["dot", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])
