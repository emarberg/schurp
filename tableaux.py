from collections import defaultdict
from partitions import StrictPartition, Shape
from numbers import MarkedNumber


class Tableau:
    def __init__(self, dictionary=None):
        if dictionary is None:
            dictionary = {}
        assert all(type(v) == MarkedNumber for v in dictionary.values())
        self.mapping = dictionary.copy()
        self.shape = Shape(self.mapping.keys())
        self.rows = {}
        self.columns = {}

    @classmethod
    def from_string(cls, string):
        def mark(i):
            i = i.strip()
            if i.endswith("'"):
                return MarkedNumber(-int(i[:-1]))
            else:
                return MarkedNumber(int(i))
        rows = [[mark(i) for i in row.split(',')] for row in string.split(';')]
        dictionary = {
            (i + 1, j + 1): rows[i][j]
            for i in range(len(rows)) for j in range(len(rows[i]))
        }
        return Tableau(dictionary)

    def shift(self):
        return Tableau({(i, i + j - 1): self.cell(i, j) for (i, j) in self.mapping})

    def cell(self, i, j):
        return self.mapping[(i, j)]

    def row(self, i):
        if i not in self.rows:
            self.rows[i] = (self.mapping[p] for p in self.shape.ordered_row(i))
        return self.rows[i]

    def column(self, j):
        if j not in self.columns:
            self.columns[j] = (self.mapping[p] for p in self.shape.ordered_column(j))
        return self.columns[j]

    @property
    def max_row(self):
        return self.shape.max_row

    @property
    def max_column(self):
        return self.shape.max_column

    def is_shifted(self):
        if any(j < i for i, j in self.mapping):
            return False

    def is_diagonally_unmarked(self):
        return not any(self.cell(i, i).is_marked() for i in range(1, self.max_row + 1))

    def is_increasing(self):
        pass

    def is_p_semistandard(self):
        pass

    def is_q_semistandard(self):
        pass

    def __repr__(self):
        width = max({len(str(v)) for i in range(1, self.max_row + 1) for v in self.row(i)} | {0})
        base = [[width * ' ' for i in range(self.max_column)] for i in range(self.max_row)]
        for i, j in self.mapping:
            v = str(self.mapping[(i, j)])
            base[i - 1][j - 1] = v + (width - len(v)) * ' '
        return '\n'.join(' '.join(row) for row in base)
