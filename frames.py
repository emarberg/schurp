from tableaux import Tableau
from permutations import Permutation
from words import get_fpf_involution_words, Word


class Frame:
    def __init__(self, t=None):
        self.tableau = Tableau() if t is None else t

    def __repr__(self):
        return str(self.tableau)

    @classmethod
    def test(cls, n):
        for w in Permutation.fpf_involutions(n):
            for e in get_fpf_involution_words(w.oneline):
                f = Frame()
                for v in e:
                    f = f.go(v)
                p, q = Word(*e).fpf_insert()
                assert p == f.compress()

    def compress(self):
        assert self.odd_cell() is None
        return Tableau({(i // 2, j // 2): self.tableau.entry(i, j) for (i, j) in self.tableau})

    def odd_cell(self):
        cells = [(i, j) for (i, j) in self.tableau if i % 2 != 0 or j % 2 != 0]
        assert len(cells) <= 1
        return cells[0] if cells else None

    def free_cell_in_row(self, row):
        if self.odd_cell():
            return
        columns = [j for (i, j) in self.tableau if i == row]
        return max(columns) + 1 if columns else row - 1

    def free_cell_in_column(self, col):
        if self.odd_cell():
            return
        rows = [i for (i, j) in self.tableau if j == col]
        return max(rows) + 1 if rows else 1

    def add(self, v):
        if self.odd_cell() is not None:
            raise Exception('Cannot add cell to \n\n%s\n' % self)
        i, j = 2, self.free_cell_in_row(2)
        t = self.tableau.set(i, j, v)
        return Frame(t)

    def go(self, v):
        ans = self.add(v)
        while ans.odd_cell() is not None:
            ans = ans.next()
        return ans

    def next(self, verbose=False):
        if self.odd_cell() is None:
            return None
        i, j = self.odd_cell()
        v, t = self.tableau.pop(i, j)
        if i % 2 == 0:
            a, b, c = t.entry(i, j - 3), t.entry(i, j - 1), t.entry(i, j + 1)
            if verbose:
                print('(row) a, b, c =', a, b, c)

            # . . . . v . . -> . . . . . v .
            if b is None and c is None:
                return Frame(t.set(i, j + 1, v))
            # . a . b v . . -> . a . b . v .
            if b is not None and c is None and b < v:
                return Frame(t.set(i, j + 1, v))

            # . . . . v c . -> (next row)
            if b is None and c is not None:
                k = Frame(t).free_cell_in_row(i + 2)
                return Frame(t.set(i + 2, k, v))

            # . . . b v c . -> (special)
            if a is None and b is not None and v < b:
                k = Frame(t).free_cell_in_column(j + 1)
                if v.number % 2 == 0:
                    return Frame(t.set(k, j + 1, b).set(i, j - 1, v))
                else:
                    u = v.increment().increment()
                    return Frame(t.set(k, j + 1, u))

            # . a . b v c . -> . a v b . c .
            if b is not None and c is not None and b < c < v:
                return Frame(t.set(i, j - 2, v))
            if a is not None and b is not None and v < a < b:
                return Frame(t.set(i, j - 2, v))

            # . a . b v c . -> . a b v . c .
            if a is not None and b is not None and a < v < b:
                return Frame(t.set(i, j - 2, b).set(i, j - 1, v))

            # . a . b a c . -> b a . b . c .
            if a is not None and b is not None and a == v < b:
                return Frame(t.set(i, j - 4, b))
        if j % 2 == 0:
            a, b, c = t.entry(i - 3, j), t.entry(i - 1, j), t.entry(i + 1, j)
            if verbose:
                print('(column) a, b, c =', a, b, c)

            # .     .
            # .     v
            # v     .
            # .     .
            # .     .
            # .     .
            # .     .
            if b is None and c is None:
                return Frame(t.set(i + 1, j, v))
            # .     .
            # .     v
            # v     .
            # b     b
            # .     .
            # a     a
            # .     .
            if b is not None and c is None and b < v:
                return Frame(t.set(i + 1, j, v))

            # .     .  (v in next column)
            # c     c
            # v     .
            # .     .
            # .     .
            # .     .
            # .     .
            if b is None and c is not None:
                k = Frame(t).free_cell_in_column(j + 2)
                return Frame(t.set(k, j + 2, v))

            # .     .  (b in next column)
            # c     c
            # v     .
            # b     v
            # .     .
            # .     .
            # .     .
            if a is None and b is not None and v < b:
                k = Frame(t).free_cell_in_column(j + 2)
                return Frame(t.set(k, j + 2, b).set(i - 1, j, v))

            # .     .
            # c     c
            # v     .
            # b     b
            # .     v
            # a     a
            # .     .
            if b is not None and c is not None and b < c < v:
                return Frame(t.set(i - 2, j, v))
            if a is not None and b is not None and v < a < b:
                return Frame(t.set(i - 2, j, v))

            # .     .
            # c     c
            # v     .
            # b     v
            # .     b
            # a     a
            # .     .
            if a is not None and b is not None and a < v < b:
                return Frame(t.set(i - 2, j, b).set(i - 1, j, v))

            # .     .
            # c     c
            # a     .
            # b     b
            # .     .
            # a     a
            # .     b
            if a is not None and b is not None and a == v < b:
                return Frame(t.set(i - 4, j, b))
