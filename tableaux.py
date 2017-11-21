from partitions import StrictPartition
from numbers import MarkedNumber


class ShiftedTableau:
    def __init__(self, **kwargs):
        if 'string' in kwargs:
            def mark(i):
                i = i.strip()
                if i.endswith("'"):
                    return MarkedNumber(-int(i[:-1]))
                else:
                    return MarkedNumber(int(i))
            self.rows = [[mark(i) for i in row.split(',')] for row in kwargs['string'].split(';')]
        elif 'rows' in kwargs:
            self.rows = kwargs['rows'][:]
        self.shape = StrictPartition(*[len(row) for row in self.rows])

    def cell(self, i, j):
        return self.rows[i - 1][j - i]

    def row(self, i):
        return self.rows[i - 1]

    def column(self, j):
        return self.columns[j - 1]

    @property
    def num_rows(self):
        return self.shape.num_rows

    @property
    def num_columns(self):
        return self.shape.num_columns

    def is_diagonally_unmarked(self):
        return not any(self.get(i, i).is_marked() for i in range(1, self.num_rows + 1))

    def is_increasing(self):
        pass

    @property
    def _cell_width(self):
        if self.rows:
            return 1 + max(len(i) for row in self.rows for i in row)

    def __repr__(self):
        width = self._cell_width
        repr_rows = [
            index * width * ' ' + ''.join(str(i) + (width - len(i)) * ' ' for i in row)
            for index, row in enumerate(self.rows)
        ]
        return '\n'.join(repr_rows)
