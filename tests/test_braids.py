from permutations import Permutation
import subprocess


def write_graph(labeled, edges, dot_filename, png_filename):
    directory = "/Users/emarberg/examples/reduced-word-graphs/" + ("labeled/" if labeled else "unlabeled/") 
    dot_filename = directory + dot_filename
    png_filename = directory + png_filename

    s = ["graph G {"]
    if not labeled:
        s += ["    node [label=\"\", shape=point; width=0.1;];"]
        s += ["    edge [color=gray];"]
    else:
        s += ["    node [margin=0; shape=plaintext; fontname=courier];"]
        s += ["    edge [color=gray];"]
    for a, b, t in edges:
        opt = ""
        if t[0] in ["c"]:
            opt = "style=solid" + ("; color=orange" if t in ["cp", "cs"] else "")
        if t[0] in ["b"]:
            opt = "" + ("color=orange" if t in ["bp", "bs"] else "")
        # if labeled:
        #    opt = ""
        s += ["    %s -- %s [%s];" % (a, b, opt)]
    s += ["}"]
    s = "\n".join(s)

    with open(dot_filename, "w") as f:
        f.write(s)
    if labeled:
        subprocess.run(["neato", "-Goverlap=scale", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])
        ps = subprocess.Popen("neato -Gstart=rand -Txdot -Goverlap=scale " + dot_filename + " | dot2tex -tmath --figonly --figpreamble=\"\\small\" > " + png_filename + ".tex", stdin=subprocess.PIPE, shell=True)
        ps.communicate()
    else:
        subprocess.run(["neato", "-Tpng", dot_filename, "-o", png_filename])
        subprocess.run(["open", png_filename])
        ps = subprocess.Popen("neato -Txdot " + dot_filename + " | dot2tex -tmath --figonly --figpreamble=\"\\small\" > " + png_filename + ".tex", stdin=subprocess.PIPE, shell=True)
        ps.communicate()


def primed_atoms_graph(w):
    w = Permutation.longest_element(w) if type(w) == int else w

    def span(n, cyclemap):
        a = tuple(n[-1].inverse().oneline)
        if len(a) == 0:
            return
        for i in range(len(a) - 1):
            if a[i] > a[i + 1]:
                nn = list(n)
                nn[cyclemap[(a[i + 1], a[i])]] ^= 1
                yield tuple(nn)
        for i in range(len(a) - 2):
            z, x, y = a[i:i + 3]
            if x < y < z:
                b = a[:i] + (y, z, x) + a[i + 3:]
                nn = list(n)
                nn[-1] = Permutation(*b).inverse()
                yield tuple(nn)
            y, z, x = a[i:i + 3]
            if x < y < z:
                b = a[:i] + (z, x, y) + a[i + 3:]
                nn = list(n)
                nn[-1] = Permutation(*b).inverse()
                yield tuple(nn)

    cycles = w.get_two_cycles()
    cyclemap = {c: i for i, c in enumerate(cycles)}
    nodes = []
    for a in w.get_atoms():
        for v in range(2**len(cycles)):
            n = []
            for i in range(len(cycles)):
                n.append(0 if v % 2 else 1)
                v = v // 2
            n.append(a)
            nodes.append(tuple(n))
    nodemap = {n: i for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b]) for a in nodes for b in span(a, cyclemap)}

    s = ["graph G {"]
    s += ["    node [label=\"\", shape=point];"]
    for a, b in edges:
        if a < b:
            s += ["    %s -- %s;" % (a, b)]
    s += ["}"]
    s = "\n".join(s)

    directory = "/Users/emarberg/examples/reduced-word-graphs/"
    dot_filename = directory + "atoms%s.dot" % str(w)
    png_filename = directory + "atoms%s.png" % str(w)
    with open(dot_filename, "w") as f:
        f.write(s)
    subprocess.run(["neato", "-Tpng", dot_filename, "-o", png_filename])
    subprocess.run(["open", png_filename])


def primed_involution_word_graph(w, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w

    def span(a):
        if len(a) > 0:
            yield (-a[0],) + a[1:], "cp"
        if len(a) >= 2 and (a[0] > 0 and a[1] > 0 and abs(a[0] - a[1]) == 1):
            yield (a[1], a[0],) + a[2:], "bs"
        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(abs(x) - abs(y)) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if abs(x) == abs(z):
                yield a[:i] + (y * abs(z) // z, abs(x), y * abs(x) // x) + a[i + 3:], "b"

    nodes = list(w.get_primed_involution_words())
    nodemap = {n: "\"" + "".join(map(lambda x: str(x) if x > 0 else str(-x) + "'", n)) + "\"" for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a) if a < b}

    dot_filename = "primed_invol%s.dot" % str(w)
    png_filename = "primed_invol%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


def involution_word_graph(w, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w

    def span(a):
        if len(a) >= 2 and (abs(a[0] - a[1]) == 1):
            yield (a[1], a[0],) + a[2:], "bs"
        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(x - y) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if x == z:
                yield a[:i] + (y, x, y) + a[i + 3:], "b"

    nodes = list(w.get_involution_words())
    nodemap = {n: "".join(map(str, n)) for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a) if a < b}
    dot_filename = "invol%s.dot" % str(w)
    png_filename = "invol%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


def hecke_involution_word_graph(w, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w

    def span(a):
        if len(a) >= 2 and (abs(a[0] - a[1]) == 1):
            x, y = a[0], a[1]
            v = Permutation.from_word(a[2:]).inverse()
            if v(x) < v(x + 1) and v(y) < v(y + 1):
                yield (x, y, x,) + a[2:], "bs"
                # yield (y, x, y,) + a[2:]
                # yield (y, x,) + a[2:]
        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(x - y) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if x == z:
                yield a[:i] + (y, x, y) + a[i + 3:], "b"

    nodes = list(w.get_involution_hecke_words())
    nodemap = {n: "".join(map(str, n)) for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a)}
    edges = {(x, y, t) for (x, y, t) in edges if (y, x, t) not in edges or x < y}
    dot_filename = "hecke%s.dot" % str(w)
    png_filename = "hecke%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


def twisted_involution_word_graph(w, rank=None, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w
    rank = w.rank if rank is None else rank

    def span(a):
        if len(a) >= 1 and abs(a[0] - (rank - a[0])) > 1:
            yield (rank - a[0],) + a[1:], "cs"
        if len(a) >= 2 and rank - a[0] == a[1]:
            yield (a[1], a[0],) + a[2:], "bs"
        if len(a) >= 4 and rank % 2 == 0:
            x, y, z = rank // 2 - 1, rank // 2, rank // 2 + 1
            if a[:4] == (y, z, x, y):
                yield (y, z, y, x) + a[4:], "bs"
            if a[:4] == (y, z, y, x):
                yield (y, z, x, y) + a[4:], "bs"
        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(x - y) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if x == z:
                yield a[:i] + (y, x, y) + a[i + 3:], "b"

    nodes = list(w.get_twisted_involution_words(rank))
    nodemap = {n: "".join(map(str, n)) for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a) if a < b}
    dot_filename = "twisted%s.dot" % str(w)
    png_filename = "twisted%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


def twisted_hecke_involution_word_graph(w, rank=None, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w
    rank = w.rank if rank is None else rank

    def span(a):
        if len(a) >= 1 and abs(a[0] - (rank - a[0])) > 1:
            x, y = a[0], rank - a[0]
            v = Permutation.from_word(a[1:]).inverse()
            if v(x) < v(x + 1) and v(y) < v(y + 1):
                yield (x, y,) + a[1:], "cs"
        if len(a) >= 2 and rank - a[0] == a[1]:
            x, y = a[0], a[1]
            v = Permutation.from_word(a[2:]).inverse()
            if abs(x - y) == 1 and v(x) < v(x + 1) and v(y) < v(y + 1):
                yield (x, y, x,) + a[2:], "bs"

        if len(a) >= 4 and rank % 2 == 0:
            y, z, x = a[:3]
            if y == rank // 2 and x == y - 1 and z == y + 1 and a[3] == y:
                v = Permutation.from_word(a[4:]).inverse()
                if v(x) < v(x + 1) < v(x + 2) < v(x + 3):
                    yield (y, z, y, x) + a[4:], "bs"
                    yield (y, z, y, x, y) + a[4:], "bs"

        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(x - y) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if x == z:
                yield a[:i] + (y, x, y) + a[i + 3:], "b"

    nodes = list(w.get_twisted_involution_hecke_words(rank))
    nodemap = {n: "".join(map(str, n)) for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a)}
    edges = {(x, y, t) for (x, y, t) in edges if (y, x, t) not in edges or x < y}
    dot_filename = "hecke_twisted%s.dot" % str(w)
    png_filename = "hecke_twisted%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


def twisted_primed_involution_word_graph(w, rank=None, labeled=False):
    w = Permutation.longest_element(w) if type(w) == int else w
    rank = w.rank if rank is None else rank

    def span(a):
        if len(a) >= 1 and abs(a[0]) == rank - abs(a[0]):
            yield (-a[0],) + a[1:], "cp"
        if len(a) >= 1 and abs(abs(a[0]) - (rank - abs(a[0]))) > 1:
            assert a[0] > 0
            yield (rank - abs(a[0]),) + a[1:], "cs"
        if len(a) >= 2 and rank - abs(a[0]) == abs(a[1]) and abs(abs(a[0]) - abs(a[1])) == 1:
            if a[0] > 0 and a[1] > 0 and abs(a[0] - a[1]) == 1:
                yield (a[1], a[0],) + a[2:], "bs"
            yield (a[0], -a[1],) + a[2:], "bp"

        if len(a) >= 4 and rank % 2 == 0:
            x, y, z = rank // 2 - 1, rank // 2, rank // 2 + 1
            if a[:4] == (y, z, x, y):
                yield (y, z, y, x) + a[4:], "bs"
                yield (y, z, x, -y) + a[4:], "bp"
            if a[:4] == (y, z, y, x):
                yield (y, z, x, y) + a[4:], "bs"
                yield (y, z, y, -x) + a[4:], "bp"
            if a[:4] == (y, z, x, -y):
                yield (y, z, x, y) + a[4:], "bs"
            if a[:4] == (y, z, y, -x):
                yield (y, z, y, x) + a[4:], "bp"

        for i in range(len(a) - 1):
            x, y = a[i], a[i + 1]
            if abs(abs(x) - abs(y)) > 1:
                yield a[:i] + (y, x) + a[i + 2:], "c"
        for i in range(len(a) - 2):
            x, y, z = a[i], a[i + 1], a[i + 2]
            if abs(x) == abs(z):
                yield a[:i] + (y * abs(z) // z, abs(x), y * abs(x) // x) + a[i + 3:], "b"

    nodes = list(w.get_twisted_primed_involution_words(rank))
    nodemap = {n: "\"" + "".join(map(lambda x: str(x) if x > 0 else str(-x) + "'", n)) + "\"" for i, n in enumerate(nodes)}
    edges = {(nodemap[a], nodemap[b], t) for a in nodes for b, t in span(a) if a < b}
    dot_filename = "primed_twisted%s.dot" % str(w)
    png_filename = "primed_twisted%s.png" % str(w)
    write_graph(labeled, edges, dot_filename, png_filename)


class H:

    @classmethod
    def braids(cls, w):
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if {a, b} == {1, 3}:
                yield w[:i] + (b, a) + w[i + 2:]
        for i in range(len(w) - 2):
            a, b, c = w[i], w[i + 1], w[i + 2]
            if (a == c == 3 and b == 2) or (a == c == 2 and b == 3):
                yield w[:i] + (b, a, b) + w[i + 3:]
        for i in range(len(w) - 4):
            a, b, c, d, e = w[i], w[i + 1], w[i + 2], w[i + 3], w[i + 4]
            if (a == c == e == 2 and b == d == 1) or (a == c == e == 1 and b == d == 2):
                yield w[:i] + (b, a, b, a, b) + w[i + 5:]

    @classmethod
    def expand(cls, base):
        ans = set(base)
        while True:
            toadd = set()
            for w in ans:
                for v in cls.braids(w):
                    if v not in ans:
                        toadd.add(v)
            if len(toadd) == 0:
                return ans
            ans |= toadd

    def __init__(self, *seed):
        self.words = self.expand({tuple(seed)})

    def __repr__(self):
        return str(min(self.words))

    def right_descent_set(self):
        return {w[-1] for w in self.words if w}

    def __mul__(self, other):
        other = (other,) if type(other) == int else other
        return self.__class__(*(next(iter(self.words)) + other))

    def __rmul__(self, other):
        other = (other,) if type(other) == int else other
        return self.__class__(*(other + next(iter(self.words))))

    @classmethod
    def generators(cls):
        return {1, 2, 3}

    @classmethod
    def m(cls, i, j):
        (i, j) = (j, i) if i > j else (i, j)
        if i == j:
            return 1
        if i == 1 and j == 2:
            return 5
        if i == 2 and j == 3:
            return 3
        return 2

    @classmethod
    def mstar(cls, i, j):
        if {cls.star(i), cls.star(j)} != {i, j}:
            return cls.m(i, j)
        if cls.m(i, j) % 2 != 0:
            return (cls.m(i, j) + 1) // 2
        if cls.star(i) == i:
            return cls.m(i, j) // 2 + 1
        return cls.m(i, j) // 2

    def __hash__(self):
        return hash(min(self.words))

    def __eq__(self, other):
        return min(self.words) == min(other.words)

    @classmethod
    def star(cls, i):
        return i

    def __len__(self):
        return len(min(self.words))

    def halfbraids(self):
        for s in self.generators():
            for t in self.generators():
                if s == t:
                    continue
                if {self.star(s), self.star(t)} != {s, t}:
                    continue
                (bestw, bestm) = (None, -1)
                for w in self.words:
                    m = 0
                    while m < len(w) and w[m] == (s if m % 2 == 0 else t):
                        m += 1
                    if m > bestm:
                        (bestw, bestm) = (w, m)
                if bestm < self.mstar(s, t):
                    continue
                for m in range(bestm, self.m(s, t) + 1):
                    yield tuple(s if i % 2 == 0 else t for i in range(m)) + bestw[bestm:]
                    yield tuple(t if i % 2 == 0 else s for i in range(m)) + bestw[bestm:]

    @classmethod
    def hecke(cls):
        ans = set()
        queue = {(cls(), cls())}
        while True:
            nextqueue = set()
            for z, w in queue:
                if z.right_descent_set() == cls.generators():
                    ans.add(w)
                asc = cls.generators() - w.right_descent_set()
                for i in asc:
                    y = z
                    if i not in z.right_descent_set():
                        y = z * i
                        if cls.star(i) * z != y:
                            y = cls.star(i) * y
                    nextqueue.add((y, w * i))
            if len(nextqueue) == 0:
                break
            queue = nextqueue
        ans = sorted(ans, key=len)
        links = set()
        for i in range(len(ans)):
            for w in ans[i].halfbraids():
                for j in range(len(ans)):
                    if i != j and w in ans[j].words:
                        links.add((i, j))
        components = []
        base = set(range(len(ans)))
        while base:
            comp = set()
            toadd = {next(iter(base))}
            while toadd:
                next_toadd = {i for (i, j) in links if j in toadd} | {j for (i, j) in links if i in toadd}
                comp |= toadd
                toadd = next_toadd - comp
            base = base - comp
            components.append(comp)
        return ans, links, components


class B(H):

    @classmethod
    def m(cls, i, j):
        (i, j) = (j, i) if i > j else (i, j)
        if i == j:
            return 1
        if i == 1 and j == 2:
            return 4
        if i == 2 and j == 3:
            return 3
        return 2

    @classmethod
    def braids(cls, w):
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if {a, b} == {1, 3}:
                yield w[:i] + (b, a) + w[i + 2:]
        for i in range(len(w) - 2):
            a, b, c = w[i], w[i + 1], w[i + 2]
            if (a == c == 3 and b == 2) or (a == c == 2 and b == 3):
                yield w[:i] + (b, a, b) + w[i + 3:]
        for i in range(len(w) - 3):
            a, b, c, d = w[i], w[i + 1], w[i + 2], w[i + 3]
            if (a == c == 2 and b == d == 1) or (a == c == 1 and b == d == 2):
                yield w[:i] + (b, a, b, a) + w[i + 4:]


class D(H):

    @classmethod
    def m(cls, i, j):
        (i, j) = (j, i) if i > j else (i, j)
        if i == j:
            return 1
        if i == 1 and j == 3:
            return 3
        if i == 2 and j == 3:
            return 3
        if i == 3 and j == 4:
            return 3
        return 2

    @classmethod
    def generators(cls):
        return {1, 2, 3, 4}

    @classmethod
    def braids(cls, w):
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if {a, b} in [{1, 2}, {1, 4}, {2, 4}]:
                yield w[:i] + (b, a) + w[i + 2:]
        for i in range(len(w) - 2):
            a, b, c = w[i], w[i + 1], w[i + 2]
            if (a == c == 1 and b == 3) or (a == c == 3 and b == 1) or (a == c == 2 and b == 3) or (a == c == 3 and b == 2) or (a == c == 4 and b == 3) or (a == c == 3 and b == 4):
                yield w[:i] + (b, a, b) + w[i + 3:]


class A(H):

    @classmethod
    def m(cls, i, j):
        (i, j) = (j, i) if i > j else (i, j)
        if i == j:
            return 1
        if i == 1 and j == 2:
            return 3
        if i == 2 and j == 3:
            return 3
        return 2

    @classmethod
    def star(cls, i):
        return 4 - i

    @classmethod
    def braids(cls, w):
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if {a, b} in [{1, 3}]:
                yield w[:i] + (b, a) + w[i + 2:]
        for i in range(len(w) - 2):
            a, b, c = w[i], w[i + 1], w[i + 2]
            if (a == c == 1 and b == 2) or (a == c == 2 and b == 1) or (a == c == 2 and b == 3) or (a == c == 3 and b == 2):
                yield w[:i] + (b, a, b) + w[i + 3:]
