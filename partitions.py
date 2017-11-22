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

    def row(self, i):
        return self.rows.get(i, set())

    def column(self, j):
        return self.columns.get(j, set())

    def ordered_row(self, i):
        return sorted(self.row(i))

    def ordered_column(self, j):
        return sorted(self.column(j))

    def __repr__(self):
        base = [[' ' for i in range(self.max_column)] for i in range(self.max_row)]
        for i, j in self.positions:
            base[i - 1][j - 1] = '*'
        return '\n'.join(' '.join(row) for row in base)


class StrictPartition:
    def __init__(self, *args):
        self.parts = tuple(sorted(args, reverse=True))
        assert len(set(self.parts)) == len(self.parts)
        assert self.parts == tuple(args)
        self.shape = Shape({
            (i + 1, i + j)
            for i in range(len(self.parts)) for j in range(1, self.parts[i] + 1)
        })

    def __hash__(self):
        return hash(self.parts)

    def __eq__(self, other):
        return self.parts == other.parts

    def __call__(self, i):
        if i <= 0 or i >= len(self.parts):
            return 0
        else:
            return self.parts[i - 1]

    def __repr__(self):
        return str(self.shape)

    def __len__(self):
        return len(self.parts)
