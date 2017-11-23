from partitions import Shape, Partition
from numbers import MarkedNumber


class Tableau:
    def __init__(self, dictionary=None):
        if dictionary is None:
            dictionary = {}
        assert all(type(v) == MarkedNumber for v in dictionary.values())
        self.mapping = dictionary.copy()
        self.max_row = max({i for i, j in self.mapping} | {0})
        self.max_column = max({j for i, j in self.mapping} | {0})

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
        return Tableau({(i, i + j - 1): self.entry(i, j) for (i, j) in self.mapping})

    def transpose(self):
        return Tableau({(j, i): self.cell(i, j) for i, j in self.mapping})

    def maximum(self):
        if self.mapping:
            return max(self.entries())

    def cells(self):
        return self.mapping.keys()

    def entries(self):
        return self.mapping.values()

    def entry(self, i, j):
        return self.mapping.get((i, j), None)

    def is_shifted(self):
        if any(j < i for i, j in self.mapping):
            return False

    def is_contiguous(self):
        values = {v.weight() for v in self.mapping.values()}
        return values == set(range(1, len(values) + 1))

    def is_unmarked(self):
        return not any(v.is_marked() for v in self.mapping.values())

    def is_diagonally_unmarked(self):
        return not any(self.entry(i, i).is_marked() for i in range(1, self.max_row + 1))

    def is_weakly_increasing(self):
        for i, j in self.cells():
            u = self.entry(i, j)
            v = self.entry(i, j + 1)
            w = self.entry(i + 1, j)
            if v and u == v and u.is_marked():
                return False
            if w and u == w and not u.is_marked():
                return False
            if (v and v < u) or (w and w < u):
                return False
        return True

    def is_increasing(self):
        for i, j in self.cells():
            u = self.entry(i, j)
            v = self.entry(i, j + 1)
            w = self.entry(i + 1, j)
            if (v and v <= u) or (w and w <= u):
                return False
        return True

    def is_semistandard(self):
        return self.is_unmarked() and self.is_weakly_increasing() and self.is_contiguous()

    def is_standard(self):
        return self.is_unmarked() and self.is_increasing() and self.is_contiguous()

    def is_shifted_semistandard(self):
        return self.is_diagonally_unmarked() and self.is_weakly_increasing() and self.is_contiguous()

    def is_shifted_standard(self):
        return self.is_diagonally_unmarked() and self.is_increasing() and self.is_contiguous()

    @classmethod
    def get_standard(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for i, j in shape.corners():
            for t in cls.get_standard(shape - Shape({(i, j)})):
                mapping = t.mapping
                mapping[(i, j)] = MarkedNumber(len(shape))
                ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_semistandard(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for border in shape.horizontal_border_strips():
            for t in cls.get_semistandard(shape - Shape(border)):
                if t.mapping:
                    n, mapping = t.maximum().increment(), t.mapping
                else:
                    n, mapping = MarkedNumber(1), {}
                for i, j in border:
                    mapping[(i, j)] = n
                ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_standard_shifted(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        ans = set()
        for i, j in shape.corners():
            for t in cls.get_standard_shifted(shape - Shape({(i, j)})):
                mapping = t.mapping
                mapping[(i, j)] = MarkedNumber(len(shape))
                ans.add(Tableau(mapping))
                if i != j:
                    mapping[(i, j)] = MarkedNumber(-len(shape))
                    ans.add(Tableau(mapping))
        return ans

    @classmethod
    def get_semistandard_shifted(cls, shape):
        if isinstance(shape, Partition):
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        borders = {
            (a, b)
            for a in shape.horizontal_border_strips() | {Shape()}
            for b in (shape - a).vertical_border_strips(exclude_diagonal=True) | {Shape()}
            if len(a) > 0 or len(b) > 0
        }

        ans = set()
        for border_h, border_v in borders:
            for t in cls.get_semistandard_shifted(shape - border_h - border_v):
                if t.mapping:
                    n, mapping = t.maximum().weight() + 1, t.mapping
                else:
                    n, mapping = 1, {}

                for i, j in border_h:
                    mapping[(i, j)] = MarkedNumber(n)
                for i, j in border_v:
                    mapping[(i, j)] = MarkedNumber(-n)
                ans.add(Tableau(mapping))
        return ans

    def __eq__(self, other):
        assert type(other) == Tableau
        return self.mapping == other.mapping

    def __hash__(self):
        return hash(tuple(sorted(self.mapping.items())))

    def __repr__(self):
        width = max({len(str(v)) for v in self.mapping.values()} | {0})
        base = [[width * ' ' for i in range(self.max_column)] for j in range(self.max_row)]
        for i, j in self.mapping:
            v = str(self.mapping[(i, j)])
            base[i - 1][j - 1] = v + (width - len(v)) * ' '
        return '\n'.join(' '.join(row) for row in base)
