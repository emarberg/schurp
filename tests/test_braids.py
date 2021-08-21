

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
