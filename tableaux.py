from partitions import Shape, Partition
from numbers import MarkedNumber


class Tableau:
    def __init__(self, dictionary=None):
        if dictionary is None:
            dictionary = {}

        self.mapping = {}
        for key, value in dictionary.items():
            if type(value) == int:
                value = MarkedNumber(value)
            self.mapping[key] = value
        self.max_row = max({i for i, j in self.mapping} | {0})
        self.max_column = max({j for i, j in self.mapping} | {0})

    def __iter__(self):
        return self.mapping.__iter__()

    def shape(self):
        return Shape(self.mapping.keys())

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

    def set(self, i, j, v):
        mapping = self.mapping.copy()
        mapping[(i, j)] = v
        return Tableau(mapping)

    def toggle(self):
        subtableau = self.find(MarkedNumber(-2))
        if subtableau:
            assert len(subtableau) == 1
            i, j = max(subtableau.cells())
            return self.set(i, j, MarkedNumber(1))

        subtableau = self.find(MarkedNumber(1))
        if subtableau:
            i, j = max(subtableau.cells())
            if i == j:
                return self.set(i, j, MarkedNumber(2))
            else:
                return self.set(i, j, MarkedNumber(-2))

        subtableau = self.find(MarkedNumber(2))
        if subtableau:
            i, j = min(subtableau.mapping)
            return self.set(i, j, MarkedNumber(1))

        return self

    def find(self, *args):
        return Tableau({key: value for key, value in self.mapping.items() if value in list(args)})

    def shift(self):
        return Tableau({(i, i + j - 1): self.entry(i, j) for (i, j) in self.mapping})

    def transpose(self):
        return Tableau({(j, i): self.entry(i, j) for i, j in self.mapping})

    def double(self):
        assert self.is_shifted()
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping}
        for i, j in self.mapping:
            if i != j:
                mapping[(j, i)] = self.entry(i, j)
        return Tableau(mapping)

    def maximum(self):
        if self.mapping:
            return max(self.entries())

    def cells(self):
        return self.mapping.keys()

    def entries(self):
        return self.mapping.values()

    def entry(self, i, j):
        return self.mapping.get((i, j), None)

    def get_row(self, i):
        columns = sorted([j for (i_, j) in self.mapping if i == i_])
        return tuple(self.entry(i, j) for j in columns)

    def get_column(self, j):
        rows = sorted([i for (i, j_) in self.mapping if j == j_])
        return tuple(self.entry(i, j) for i in rows)

    def is_shifted(self):
        return not any(j < i for i, j in self.mapping)

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

    def __len__(self):
        return len(self.mapping)

    def __nonzero__(self):
        return len(self.mapping) > 0

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
        rows = [' '.join(row) for row in base]
        # return '\n'.join(reversed(rows))  # French
        return '\n'.join(rows)            # English

    @classmethod
    def decreasing_part(cls, row):
        if row:
            i = 1
            while i < len(row) and row[i] < row[i - 1]:
                i += 1
            return row[:i]

    @classmethod
    def increasing_part(cls, row):
        if row:
            rev = cls.decreasing_part(list(reversed(row)))
            return tuple(reversed(rev))[1:]

    # def is_decomposition_tableau(self):
    #     def is_row_unimodal(row):
    #         return len(self.decreasing_part(row)) + len(self.increasing_part(row)) == len(row)

    def replace_row(self, j, newrow, shifted=False):
        dictionary = {(j_, k): self.entry(j_, k) for (j_, k) in self.mapping if j != j_}
        for k_zerobased, v in enumerate(newrow):
            k = k_zerobased + 1
            assert type(v) == MarkedNumber
            if shifted:
                k += j - 1
            dictionary[(j, k)] = v
        return Tableau(dictionary)

    def replace_column(self, j, newcol):
        dictionary = {(i, j_): self.entry(i, j_) for (i, j_) in self.mapping if j != j_}
        for i_zerobased, v in enumerate(newcol):
            i = i_zerobased + 1
            assert type(v) == MarkedNumber
            dictionary[(i, j)] = v
        return Tableau(dictionary)

    def add_to_column(self, j, v):
        return self.replace_column(j, self.get_column(j) + (v, ))

    def add_to_row(self, j, v, shifted=False):
        return self.replace_row(j, self.get_row(j) + (v, ), shifted)

    @classmethod
    def bump(cls, p, column_dir, sequence):
        assert all(sequence[i + 1] > sequence[i] for i in range(len(sequence) - 1))
        if len(sequence) == 0 or p > sequence[-1]:
            newseq = sequence + (p,)
            q = None
        elif p == sequence[-1]:
            newseq = sequence
            q = None
        else:
            if p <= sequence[0]:
                i = 0
                column_dir = True
            else:
                i = [j for j in range(1, len(sequence)) if sequence[j - 1] < p <= sequence[j]][0]
            if p == sequence[i]:
                newseq = sequence
                q = sequence[i + 1]
            else:
                newseq = sequence[:i] + (p,) + sequence[i + 1:]
                q = sequence[i]
        return q, column_dir, newseq

    def hecke_insert(self, p, j=0):
        if p is None:
            return (j, self)

        def hecke_bump(a, tup):
            if len(tup) == 0 or a > tup[-1]:
                newtup = tup + (a,)
                q = None
            elif p == tup[-1]:
                newtup = tup
                q = None
            else:
                i = [j for j in range(len(tup)) if a < tup[j]][0]
                if i > 0 and a == tup[i - 1]:
                    newtup = tuple(tup)
                    q = tup[i]
                else:
                    newtup = tup[:i] + (a,) + tup[i + 1:]
                    q = tup[i]
            return q, newtup

        j += 1
        row = self.get_row(j)
        p, row = hecke_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)
        if tab.is_increasing():
            return tab.hecke_insert(p, j)
        else:
            return self.hecke_insert(p, j)

    def shifted_hecke_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = self.bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = self.bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if tab.is_increasing():
            return tab.shifted_hecke_insert(p, j, column_dir, verbose=verbose)
        else:
            return self.shifted_hecke_insert(p, j, column_dir, verbose=verbose)

    def eg_insert(self, p, j=0):
        if p is None:
            return (j, self)

        def eg_bump(a, tup):
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if a == b:
                    b = tup[i + 1]
                    new = tup
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, new)
            return (None, tup + (a,))

        j += 1
        row = self.get_row(j)
        p, row = eg_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)

        assert tab.is_increasing()
        return tab.eg_insert(p, j)

    def involution_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def involution_bump(a, cdir, tup):
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if a == b:
                    b = tup[i + 1]
                    new = tup
                    cdir = cdir or (i == 0)
                elif not cdir and i == 0:
                    new = (a,) + tup[1:]
                    cdir = True
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = involution_bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = involution_bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if verbose:
            print(tab, '\n')
        assert tab.is_increasing()
        return tab.involution_insert(p, j, column_dir, verbose=verbose)

    def fpf_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def fpf_bump(a, cdir, tup):
            # inserting `a` (number) into `tup` (tuple) in row (`cdir=False`) or column direction
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if a == b:
                    b = tup[i + 1]
                    assert b == a.increment()
                    new = tup
                elif not cdir and i == 0:
                    cdir = True
                    if a.number % 2 == 0:
                        new = (a,) + tup[1:]
                    else:
                        b = a.increment().increment()
                        new = tup
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if verbose:
            if column_dir:
                print('Inserting %s into column %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))
            else:
                print('Inserting %s into row %s of \n%s\n' % (
                    str(p),
                    str(j),
                    self
                ))

        if column_dir:
            p, column_dir, col = fpf_bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = fpf_bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        assert tab.is_increasing()
        return tab.fpf_insert(p, j, column_dir, verbose=verbose)
