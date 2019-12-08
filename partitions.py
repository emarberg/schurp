from collections import defaultdict


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
