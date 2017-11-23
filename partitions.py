from collections import defaultdict


class Shape:
    def __init__(self, positions=None):
        if positions is None:
            positions = set()
        assert not any(i <= 0 or j <= 0 for i, j in positions)
        self.positions = {(i, j) for i, j in positions}
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

    def transpose(self):
        return Shape({(j, i) for (i, j) in self.positions})

    def row(self, i):
        return Shape({(j, k) for (j, k) in self.positions if i == j})

    def column(self, i):
        return Shape({(j, k) for (j, k) in self.positions if i == k})

    def northeast_position(self):
        if self.positions:
            return min(self.positions, key=lambda x: (x[0], -x[1]))

    def corners(self):
        return {
            (i, j) for i, j in self.positions
            if (i + 1, j) not in self.positions and (i, j + 1) not in self.positions
        }

    def vertical_border_strips(self, exclude_diagonal=False):
        shapes = {x.transpose() for x in self.transpose().horizontal_border_strips()}
        return {
            x for x in shapes
            if not exclude_diagonal or not any(i == j for i, j in x)
        }

    def horizontal_border_strips(self):
        borders = self._horizontal_border_strips_helper()
        return {x for x in borders if len(x) > 0}

    def _horizontal_border_strips_helper(self):
        if len(self) == 0:
            return {Shape()}

        i, k = self.northeast_position()
        subshape = self - self.row(i)
        borders = subshape._horizontal_border_strips_helper()

        ans = set()
        j = k
        while True:
            for border in borders:
                ans.add(border + Shape({(i, l) for l in range(j + 1, k + 1)}))
            if (i + 1, j) in self.positions or (i, j) not in self.positions:
                break
            j -= 1
        return ans

    def __eq__(self, other):
        assert type(other) == Shape
        return self.positions == other.positions

    def __len__(self):
        return len(self.positions)

    def __repr__(self):
        base = [[' ' for i in range(self.max_column)] for i in range(self.max_row)]
        for i, j in self.positions:
            base[i - 1][j - 1] = '*'
        return '\n'.join(' '.join(row) for row in base)


class Partition:
    def __init__(self, *args):
        self.parts = tuple(sorted(args, reverse=True))
        assert self.parts == tuple(args)
        self.shape = Shape({
            (i + 1, j + 1) for i in range(len(self.parts)) for j in range(self.parts[i])
        })

    def __hash__(self):
        return hash(self.parts)

    def __eq__(self, other):
        return type(self) == type(other) and self.parts == other.parts

    def __call__(self, i):
        if i <= 0 or i >= len(self.parts):
            return 0
        else:
            return self.parts[i - 1]

    def __repr__(self):
        return str(self.shape)

    def __len__(self):
        return len(self.parts)


class StrictPartition(Partition):
    def __init__(self, *args):
        self.parts = tuple(sorted(args, reverse=True))
        assert len(set(self.parts)) == len(self.parts)
        assert self.parts == tuple(args)
        self.shape = Shape({
            (i + 1, i + j + 1) for i in range(len(self.parts)) for j in range(self.parts[i])
        })
