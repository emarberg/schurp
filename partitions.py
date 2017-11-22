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
        self.max_row = max(set(self.rows) | {0})
        self.max_column = max(set(self.columns) | {0})

    def __iter__(self):
        return self.positions.__iter__()

    def __sub__(self, shape):
        assert type(shape) == Shape and shape.positions.issubset(self.positions)
        return Shape(self.positions - shape.positions)

    def transpose(self):
        return Shape({(j, i) for (i, j) in self.positions})

    def row(self, i):
        return self.rows.get(i, set())

    def column(self, j):
        return self.columns.get(j, set())

    def ordered_row(self, i):
        return sorted(self.row(i))

    def ordered_column(self, j):
        return sorted(self.column(j))

    def northeast_position(self):
        if self.positions:
            return min(self.positions, key=lambda x: (x[0], -x[1]))

    def corners(self):
        rowends = {max(row) for row in self.rows.values()}
        colends = {max(col) for col in self.columns.values()}
        return rowends & colends

    def horizontal_border_strips(self):
        borders = self._horizontal_border_strips_helper()
        return {x for x in borders if x}

    def _horizontal_border_strips_helper(self):
        if len(self) == 0:
            return {()}

        i, k = self.northeast_position()
        subshape = Shape(self.positions - self.row(i))
        borders = subshape._horizontal_border_strips_helper()

        ans = set()
        j = k
        while True:
            for border in borders:
                newborder = border + tuple((i, l) for l in range(j + 1, k + 1))
                ans.add(tuple(sorted(newborder)))

            if (i + 1, j) in self.positions or (i, j) not in self.positions:
                break

            j -= 1
        return ans

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
