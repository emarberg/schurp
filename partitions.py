from collections import defaultdict, deque


class Shape:
    def __init__(self, positions=None):
        if positions is None:
            positions = set()
        self.positions = {(i, j) for i, j in positions}
        assert not any(i <= 0 or j <= 0 for i, j in self.positions)
        self.rows = defaultdict(set)
        self.columns = defaultdict(set)
        for i, j in self.positions:
            self.rows[i].add((i, j))
            self.columns[j].add((i, j))
        self.max_row = max({i for i, j in self.positions} | {0})
        self.max_column = max({j for i, j in self.positions} | {0})

    def __iter__(self):
        return self.positions.__iter__()

    def __hash__(self):
        return hash(tuple(sorted(self.positions)))

    def __sub__(self, shape):
        assert type(shape) == Shape and shape.positions.issubset(self.positions)
        return Shape(self.positions - shape.positions)

    def __add__(self, shape):
        assert type(shape) == Shape and not (shape.positions & self.positions)
        return Shape(self.positions | shape.positions)

    def contains(self, other):
        assert type(other) == Shape
        return self.positions.issuperset(other.positions)

    def transpose(self):
        return Shape({(j, i) for (i, j) in self.positions})

    def row(self, i):
        return Shape({(j, k) for (j, k) in self.positions if i == j})

    def column(self, i):
        return Shape({(j, k) for (j, k) in self.positions if i == k})

    def justify(self):
        if len(self.positions) == 0:
            return Shape()
        m = min(i for i, j in self.positions) - 1
        n = min(j for i, j in self.positions) - 1
        return Shape({(i - m, j - n) for i, j in self.positions})

    def northeast_position(self):
        if self.positions:
            return min(self.positions, key=lambda x: (x[0], -x[1]))

    def corners(self):
        return Shape({
            (i, j) for i, j in self.positions
            if (i + 1, j) not in self.positions and (i, j + 1) not in self.positions
        })

    def vertical_border_strips(self, exclude_diagonal=False):
        shapes = {x.transpose() for x in self.transpose().horizontal_border_strips()}
        return {
            x for x in shapes
            if not exclude_diagonal or not any(i == j for i, j in x)
        }

    def horizontal_border_strips(self):
        borders = self._horizontal_border_strips_helper()
        return {x for x in borders if x}

    def _horizontal_border_strips_helper(self):
        if len(self) == 0:
            yield Shape()
            return

        i, k = min(self.positions, key=lambda x: (x[0], -x[1]))  # northeast position
        subshape = self - self.row(i)

        for border in subshape._horizontal_border_strips_helper():
            j = k
            while True:
                yield border + Shape({(i, l) for l in range(j + 1, k + 1)})
                if (i + 1, j) in self.positions or (i, j) not in self.positions:
                    break
                j -= 1

    def __eq__(self, other):
        assert type(other) == Shape
        return self.positions == other.positions

    def __len__(self):
        return len(self.positions)

    def __nonzero__(self):
        return len(self.positions) > 0

    def __repr__(self):
        base = [[' ' for i in range(self.max_column)] for i in range(self.max_row)]
        for i, j in self.positions:
            base[i - 1][j - 1] = '*'
        return '\n'.join(' '.join(row) for row in base)


class Partition:
    def __init__(self, *args):
        self.parts = sorted(args, reverse=True)
        assert self.parts == list(args)

        while self.parts and self.parts[-1] == 0:
            self.parts.pop()
        # self.parts = tuple(self.parts)

        self.shape = Shape({
            (i + 1, j + 1) for i in range(len(self.parts)) for j in range(self.parts[i])
        })

    def add_box_to_row(self, row):
        p = self
        parts = [p(i) + (1 if i == row else 0) for i in range(1, max(len(p), row) + 1)]
        return Partition(*parts)

    def add_box_to_column(self, col, shift):
        shape = {(i, i + j - 1) for (i, j) in self.shape} if shift else set(self.shape)
        a = [i for (i, j) in shape if j == col]
        row = max(a) + 1 if a else 1
        shape.add((row, col))
        if shift:
            assert all((i, j - 1) in shape or i == j for (i, j) in shape)
        else:
            assert all((i, j - 1) in shape or j == 1 for (i, j) in shape)
        parts = []
        for (i, j) in shape:
            while i - 1 >= len(parts):
                parts += [0]
            parts[i - 1] += 1
        return Partition(*parts)

    @classmethod
    def union(cls, p, q):
        parts = [max(p(i), q(i)) for i in range(1, max(len(p), len(q)) + 1)]
        return Partition(*parts)

    @classmethod
    def shifted_growth_diagram(cls, dictionary, m=None, n=None):
        def shdiff(nu, lam):
            return next(iter(Partition.skew(nu, lam, shifted=True)))

        dictionary = {(a, i + 1) for i, a in enumerate(dictionary)} if type(dictionary) in [list, tuple] else dictionary
        dictionary = {k: 1 for k in dictionary} if type(dictionary) == set else dictionary
        n = max([0] + [i for i, _ in dictionary]) if n is None else n
        m = max([0] + [j for _, j in dictionary]) if m is None else m

        g = [[Partition() for _ in range(m + 1)] for _ in range(n + 1)]
        edges = [[False for _ in range(m + 1)] for _ in range(n + 1)]
        corners = [[None for _ in range(m + 1)] for _ in range(n + 1)]

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                v = dictionary.get((i, j), 0)
                assert v in [0, 1]
                lam, nu, mu = g[i - 1][j - 1], g[i - 1][j], g[i][j - 1]
                if v == 1 and lam(1) == mu(1):
                    # case (1)
                    # print('case 1')
                    gamma = mu.add_box_to_row(row=1)
                elif v == 1:
                    # case (2)
                    # print('case 2')
                    assert lam(1) + 1 == mu(1)
                    gamma = mu
                    corners[i][j] = 1
                elif mu == lam:
                    # case (3a)
                    # print('case 3a')
                    gamma = nu
                    edges[i][j] = edges[i - 1][j]
                    corners[i][j] = corners[i - 1][j]
                elif nu == lam and not edges[i - 1][j] and corners[i - 1][j] is None:
                    # case (3b)
                    # print('case 3b')
                    gamma = mu
                elif not mu.contains(nu):
                    # case (4)
                    # print('case 4')
                    gamma = Partition.union(nu, mu)
                    edges[i][j] = edges[i - 1][j]
                elif lam != nu:
                    a, b = shdiff(nu, lam)
                    row = {(x, x + y - 1) for (x, y) in mu.shape - lam.shape if x == a + 1}
                    col = {(x, x + y - 1) for (x, y) in mu.shape - lam.shape if x + y - 1 == b + 1}

                    if not edges[i - 1][j] and a != b:
                        if len(row) == 0:
                            # case (5)
                            gamma = mu.add_box_to_row(a + 1)
                        else:
                            # case (6)
                            gamma = mu
                            corners[i][j] = a + 1
                    else:
                        if len(col) == 0:
                            # case (7)
                            gamma = mu.add_box_to_column(b + 1, shift=True)
                            edges[i][j] = True
                        else:
                            # case (8)
                            gamma = mu
                            edges[i][j] = True
                            corners[i][j] = b + 1
                elif lam == nu and not edges[i - 1][j]:
                    a = corners[i - 1][j]
                    b = nu(a) + a - 1
                    skew = {(x, x + y - 1) for (x, y) in mu.shape - lam.shape}
                    if (a, b + 1) not in skew and (a + 1, b) not in skew:
                        # case (9)
                        gamma = mu
                        corners[i][j] = a
                    elif (a + 1, b) in skew:
                        # case (10)
                        gamma = mu
                        corners[i][j] = a + 1
                    elif (a, b + 1) in skew:
                        if any(x == a + 1 for (x, y) in skew):
                            # case (??)
                            gamma = mu
                            corners[i][j] = a + 1
                        elif a == b:
                            # case (12)
                            gamma = mu
                            edges[i][j] = True
                            corners[i][j] = b + 1
                        else:
                            # case (11)
                            gamma = mu.add_box_to_row(a + 1)
                    else:
                        raise Exception

                elif lam == nu and edges[i - 1][j]:
                    b = corners[i - 1][j]
                    a = max([x for (x, y) in nu.shape if x + y - 1 == b])
                    skew = {(x, x + y - 1) for (x, y) in mu.shape - lam.shape}

                    if (a, b + 1) not in skew and (a + 1, b) not in skew:
                        # case (13)
                        gamma = mu
                        corners[i][j] = b
                        edges[i][j] = True
                    elif (a, b + 1) in skew:
                        # case (14)
                        gamma = mu
                        corners[i][j] = b + 1
                        edges[i][j] = True
                    elif (a + 1, b) in skew:
                        if any(y == b + 1 for (x, y) in skew):
                            # case (??)
                            gamma = mu
                            corners[i][j] = b + 1
                            edges[i][j] = True
                        else:
                            # case (15)
                            gamma = mu.add_box_to_column(b + 1, shift=True)
                            edges[i][j] = True
                    else:
                        raise Exception
                else:
                    raise Exception

                g[i][j] = gamma
        return g, edges, corners

    @classmethod
    def growth_diagram(cls, dictionary, m=None, n=None):
        dictionary = {k: 1 for k in dictionary} if type(dictionary) == set else dictionary
        n = max([0] + [i for i, _ in dictionary]) if n is None else n
        m = max([0] + [j for _, j in dictionary]) if m is None else m
        g = [[Partition() for _ in range(m + 1)] for _ in range(n + 1)]
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                carry, k = dictionary.get((i, j), 0), 1
                rho, mu, nu = g[i - 1][j - 1], g[i - 1][j], g[i][j - 1]
                lam = []
                while True:
                    v = max(mu(k), nu(k)) + carry
                    if v == 0:
                        break
                    lam.append(v)
                    carry, k = min(mu(k), nu(k)) - rho(k), k + 1
                g[i][j] = Partition(*lam)
        return g

    @classmethod
    def print_growth_diagram(cls, g):
        g = cls.growth_diagram(g) if type(g) != list else g
        g = [[str(mu) for mu in row] for row in g]
        m = max([6] + [len(mu) for row in g for mu in row])
        g = [[mu + (m - len(mu)) * ' ' for mu in row] for row in g]
        print()
        for row in reversed(g):
            print(' '.join(row))
        print()

    @classmethod
    def skew(cls, mu, nu, shifted=False):
        s = set(mu.shape - nu.shape)
        return {(i, i + j - 1) for (i, j) in s} if shifted else s

    def tuple(self):
        return tuple(self.parts)

    def is_strict(self):
        return self.is_strict_partition(self.parts)

    @classmethod
    def is_strict_partition(cls, mu):
        return all(mu[i - 1] > mu[i] for i in range(1, len(mu)))

    def is_symmetric(self):
        return self == self.transpose()

    def __iter__(self):
        return self.parts.__iter__()

    def decrement(self, i):
        assert self(i) > self(i + 1)
        self.parts[i - 1] -= 1
        if self.parts[i - 1] == 0:
            self.parts = self.parts[:i - 1]

    def transpose(self):
        if len(self.parts) == 0:
            return Partition()
        else:
            return Partition(*[
                len([a for a in self.parts if a > i])
                for i in range(max(self.parts))
            ])

    def pieri(self, i):
        ranges = self._pieri_ranges(i)

        return {
            self._pieri_key(delta): self._pieri_value(delta)
            for delta in self._pieri_helper(i, len(ranges) - 1, ranges)
        }

    def _pieri_key(self, delta):
        return self.__class__(*[self(j + 1) + a for j, a in enumerate(delta)])

    def _pieri_value(self, delta):
        return 1

    def _pieri_ranges(self, i):
        return [i] + [self(j) - self(j + 1) for j in range(1, len(self) + 1)]

    @classmethod
    def _pieri_helper(cls, i, j, ranges):
        if j == 0:
            yield (i,)
            return
        for t in range(min(i, ranges[j]) + 1):
            for delta in cls._pieri_helper(i - t, j - 1, ranges):
                yield delta + (t,)

    def __hash__(self):
        return hash(tuple(self.parts))

    def __lt__(self, other):
        assert type(self) == type(other)
        return self.parts < other.parts

    def __eq__(self, other):
        assert type(self) == type(other)
        return self.parts == other.parts

    def __call__(self, i):
        if i <= 0 or i > len(self.parts):
            return 0
        else:
            return self.parts[i - 1]

    def __repr__(self):
        return str(self.parts)

    def __len__(self):
        return len(self.parts)

    def __abs__(self):
        ans = 0
        for i in self.parts:
            ans += i
        return ans

    def __nonzero__(self):
        return len(self.parts) > 0

    def to_grassmannian(self):
        from permutations import Permutation
        if len(self.parts) == 0:
            return Permutation()
        oneline = []
        for i, a in enumerate(reversed(self.parts)):
            oneline += [a + i + 1]
        n = max(oneline)
        oneline += [i for i in range(1, n + 1) if i not in oneline]
        return Permutation(oneline)

    def contains(self, other):
        assert type(other) == type(self)
        return self.shape.contains(other.shape)

    @classmethod
    def _contains(cls, bigger, smaller):
        """Returns true if mu subseteq nu as partitions."""
        if len(smaller) > len(bigger):
            return False
        return all(0 <= smaller[i] <= bigger[i] for i in range(len(smaller)))

    def compact(self):
        if self.parts:
            return ','.join([str(i) for i in self.parts])
        return '0'

    @classmethod
    def all(cls, n, max_part=None):
        if max_part is None:
            max_part = n

        if n == 0:
            yield Partition()
        else:
            for i in range(1, 1 + min(max_part, n)):
                for p in cls.all(n - i, i):
                    parts = [i] + p.parts
                    yield Partition(*parts)

    @classmethod
    def skew_pairs(cls, n):
        for mu in cls.all(n):
            for m in range(n + 1):
                for nu in cls.all(m):
                    if mu.contains(nu):
                        yield mu, nu

    @classmethod
    def trim(cls, mu):
        while mu and mu[-1] == 0:
            mu = mu[:-1]
        return tuple(mu)

    @classmethod
    def is_partition(cls, mu):
        return all(mu[i - 1] >= mu[i] for i in range(1, len(mu))) and (mu == () or mu[-1] >= 0)

    @classmethod
    def subpartitions(cls, mu, strict=False):

        def _subpartitions(mu, strict):
            if mu:
                for nu in _subpartitions(mu[1:], strict):
                    lb = (nu[0] + (1 if strict else 0)) if (nu and nu[0] > 0) else 0
                    ub = mu[0]
                    for a in range(lb, ub + 1):
                        yield (a,) + nu
            else:
                yield ()

        for nu in _subpartitions(mu, strict):
            yield cls.trim(nu)

    @classmethod
    def flags(cls, n):
        q = deque([()])
        while q:
            flag = q.popleft()
            if len(flag) == n:
                yield flag
                continue
            minpart = 1 if len(flag) == 0 else max(len(flag) + 1, flag[-1])
            for a in range(minpart, n + 1):
                q.append(flag + (a,))


class StrictPartition(Partition):
    def __init__(self, *args):
        if len(args) == 1 and type(args[0]) == Partition:
            args = args[0].parts
        self.parts = list(sorted(args, reverse=True))
        assert self.parts == list(args)
        assert all(p >= 0 for p in self.parts)

        while self.parts and self.parts[-1] == 0:
            self.parts.pop()
        # self.parts = tuple(self.parts)

        assert len(set(self.parts)) == len(self.parts)

        self.shape = Shape({
            (i + 1, i + j + 1) for i in range(len(self.parts)) for j in range(self.parts[i])
        })

    def __call__(self, i):
        return self.parts[i - 1] if 0 <= i - 1 < len(self.parts) else 0

    def _pieri_value(self, delta):
        return 2**(sum(x != 0 for x in delta) - 1)

    def _pieri_ranges(self, i):
        return [i] + [self(j) - self(j + 1) - 1 for j in range(1, len(self) + 1)]

    def to_grassmannian(self):
        from permutations import Permutation
        n = self.parts[0] if self.parts else 0
        w = Permutation()
        for i in range(len(self.parts)):
            w *= Permutation.cycle([n + 1 + i, n + 1 - self.parts[i]])
        return w

    def __add__(self, other):
        assert type(other) in [Partition, StrictPartition]
        parts = [self(i) + other(i) for i in range(1, 1 + max(len(self), len(other)))]
        return StrictPartition(*parts)

    @classmethod
    def symmetric_double(cls, mu):
        if type(mu) in [StrictPartition, Partition]:
            mu = mu.tuple()
        shape = {(i + 1, i + j + 1) for i in range(len(mu)) for j in range(mu[i])}
        shape = shape | {(j, i) for (i, j) in shape}
        mu = []
        for i, j in shape:
            while i - 1 >= len(mu):
                mu.append(0)
            mu[i - 1] += 1
        return tuple(mu)

    @classmethod
    def skew_symmetric_double(cls, mu):
        if type(mu) in [StrictPartition, Partition]:
            mu = mu.tuple()
        shape = {(i + 1, i + j + 2) for i in range(len(mu)) for j in range(mu[i])}
        shape = shape | {(j, i) for (i, j) in shape}
        i = 1
        while (i, i + 1) in shape:
            shape.add((i, i))
            i += 1
        if i > 1 and (i - 1, i + 1) not in shape:
            shape.add((i, i))
        mu = []
        for i, j in shape:
            while i - 1 >= len(mu):
                mu.append(0)
            mu[i - 1] += 1
        return tuple(mu)

    @classmethod
    def all(cls, n, max_part=None):
        if max_part is None:
            max_part = n

        if n == 0:
            yield StrictPartition()
        else:
            for i in range(1, 1 + min(max_part, n)):
                for p in cls.all(n - i, i - 1):
                    parts = [i] + p.parts
                    yield StrictPartition(*parts)
