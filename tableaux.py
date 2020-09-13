from cached import cached_value
from partitions import Shape, Partition, StrictPartition
from marked import MarkedNumber
import random
from collections import defaultdict


STANDARD_CACHE = {}
SEMISTANDARD_CACHE = {}
STANDARD_SHIFTED_MARKED_CACHE = {}
HORIZONTAL_STRIPS_CACHE = {}


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

    def shifted_crystal_e(self, index):
        pass

    def shifted_crystal_f(self, index):
        pass

    def restrict(self, n):
        n = MarkedNumber(n) if type(n) == int else n
        assert type(n) == MarkedNumber
        return Tableau({k: v for k, v in self.mapping.items() if v <= n})

    def durfee(self):
        i = 0
        while (i + 1, i + 1) in self.mapping:
            i += 1
        return i

    @classmethod
    def from_composition(cls, alpha):
        ans = Tableau()
        for i, a in enumerate(alpha):
            for col in range(1, a + 1):
                row = 1
                while (row, col) in ans:
                    row += 1
                ans = ans.set(row, col, i + 1)
        return ans

    def weight(self, as_dict=False):
        ans = {}
        for v in self.mapping.values():
            v = v.weight()
            if v not in ans:
                ans[v] = 0
            ans[v] += 1
        if as_dict:
            return ans
        assert all(type(m) == int and m > 0 for m in ans)
        m = max(ans) if ans else 0
        alpha = m * [0]
        for v in ans:
            alpha[v - 1] = ans[v]
        return tuple(alpha)

    def shape(self):
        return Shape(self.mapping.keys())

    def partition(self):
        rows = defaultdict(int)
        for i, j in self.mapping:
            rows[i] += 1
        return Partition(*sorted(rows.values(), reverse=True))

    def count_diagonal_cells(self):
        return len([(i, j) for i, j in self.mapping if i == j])

    @classmethod
    def from_string(cls, string):
        def mark(i):
            i = i.strip()
            if i == "":
                return None
            if i.endswith("'"):
                return MarkedNumber(-int(i[:-1]))
            else:
                return MarkedNumber(int(i))
        rows = [[mark(i) for i in row.split(',')] for row in string.split(';')]
        dictionary = {
            (i + 1, j + 1): rows[i][j]
            for i in range(len(rows)) for j in range(len(rows[i]))
            if rows[i][j]
        }
        return Tableau(dictionary)

    @classmethod
    def from_partition(cls, mu):
        mapping = {(i, j): 1 for (i, j) in mu.shape}
        return Tableau(mapping)

    def add(self, i, j, v):
        return self.set(i, j, v)

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
        args = [MarkedNumber(v) if type(v) == int else v for v in args]
        return Tableau({key: value for key, value in self.mapping.items() if value in args})

    def shift(self):
        return Tableau({(i, i + j - 1): self.entry(i, j) for (i, j) in self.mapping})

    def translate_left(self, steps=1):
        return Tableau({(i, j - steps): self.entry(i, j) for i, j in self.mapping})

    def transpose(self):
        return Tableau({(j, i): self.entry(i, j) for i, j in self.mapping})

    def double(self, shift=False):
        assert self.is_shifted()
        mapping = {(i, j + shift): self.entry(i, j) for i, j in self.mapping}
        for i, j in self.mapping:
            if (j, j) not in mapping or shift:
                mapping[(j, i)] = self.entry(i, j)
        return Tableau(mapping)

    def fpf_double(self):
        assert self.is_shifted()
        mapping = {(i, j + 1): self.entry(i, j) for i, j in self.mapping}
        for i, j in self.mapping:
            mapping[(j + 1, i)] = self.entry(i, j)
        offset = 0 if (1, 1) not in self.mapping else (self.entry(1, 1).number - 2) // 2
        for i, j in list(mapping.keys()):
            if (i, i + 1) in mapping or (i, i - 1) in mapping:
                mapping[(i, i)] = 2 * (i + offset) - 1
        return Tableau(mapping)

    def halve(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i <= j}
        return Tableau(mapping)

    def lower_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i <= j}
        return Tableau(mapping)

    def strict_lower_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i < j}
        return Tableau(mapping)

    def upper_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i >= j}
        return Tableau(mapping)

    def strict_upper_half(self):
        mapping = {(i, j): self.entry(i, j) for i, j in self.mapping if i > j}
        return Tableau(mapping)

    def maximum(self):
        if self.mapping:
            return max(self.entries())

    def cells(self):
        return self.mapping.keys()

    def entries(self):
        return self.mapping.values()

    def get(self, i, j):
        return self.entry(i, j)

    def entry(self, i, j):
        return self.mapping.get((i, j), None)

    def pop(self, i, j):
        assert (i, j) in self
        mapping = self.mapping.copy()
        del mapping[(i, j)]
        return self.entry(i, j), Tableau(mapping)

    def get_main_diagonal(self):
        return tuple(self.entry(i, i) for i in range(1, self.max_row + 1) if (i, i) in self.mapping)

    def get_row(self, i):
        columns = sorted([j for (i_, j) in self.mapping if i == i_])
        return tuple(self.entry(i, j) for j in columns)

    def get_column(self, j):
        rows = sorted([i for (i, j_) in self.mapping if j == j_])
        return tuple(self.entry(i, j) for i in rows)

    def get_hook(self, j, k):
        ans = {(x, y) for (x, y) in self if x == j and k <= y}
        ans |= {(x, y) for (x, y) in self if y == k and j < x}
        return ans

    def get_shifted_hook(self, j, k):
        assert self.is_shifted()
        ans = self.get_hook(j, k)
        if (k, k) in ans:
            ans |= {(x, y) for (x, y) in self if x == k + 1}
        return ans

    def hooks(self):
        return Tableau({(i, j): len(self.get_hook(i, j)) for (i, j) in self})

    def shifted_hooks(self):
        assert self.is_shifted()
        return Tableau({(i, j): len(self.get_shifted_hook(i, j)) for (i, j) in self})

    @classmethod
    def random_shifted(cls, mu):
        mu = StrictPartition(*mu.parts[:])
        nu = defaultdict(int)
        for i, j in mu.shape:
            nu[j] += 1
        nu = dict(nu)
        t = Tableau.from_partition(mu)

        for n in range(abs(mu), 0, -1):
            m = random.randint(1, abs(mu))
            i = 1
            while m > mu(i) and i <= len(mu):
                m -= mu(i)
                i += 1
            j = m + i - 1

            while True:
                row = mu(i) - j + i - 1
                col = nu[j] - i
                res = mu(j + 1)
                if row == 0 and col == 0 and res == 0:
                    break
                m = random.randint(1, row + col + res)
                if m <= row:
                    j += m
                elif row + col < m:
                    i = j + 1
                    j = j + m - row - col
                else:
                    i += m - row

            t.mapping[(i, j)] = MarkedNumber(-n if i != j and random.random() < 0.5 else n)
            mu.decrement(i)
            nu[j] -= 1
            if n % 1000 == 0:
                print(n, '. . .')
        return t

    @classmethod
    def inverse_fpf(cls, p, q):
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}

        order = len(q) * [0]
        for (i, j) in q:
            order[abs(q[(i, j)]) - 1] = (i, j, q[(i, j)] < 0)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j, signed = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]
            if signed:
                for col in range(j - 1, 0, -1):
                    row = i
                    while p.get((row + 1, col), a) < a:
                        row += 1
                    if row == col:
                        i, j = col, col
                        if a % 2 == 0:
                            a, p[(i, j)] = p[(i, j)], a
                        else:
                            a -= 2
                        break
                    elif (row + 1, col) in p and p[(row + 1, col)] == a:
                        a = p[(row, col)]
                    else:
                        a, p[(row, col)] = p[(row, col)], a
            for row in range(i - 1, 0, -1):
                col = j
                while p.get((row, col + 1), a + 1) <= a:
                    col += 1
                if p[(row, col)] == a:
                    a = p[(row, col - 1)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def inverse_inv(cls, p, q):
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}

        order = len(q) * [0]
        for (i, j) in q:
            order[abs(q[(i, j)]) - 1] = (i, j, q[(i, j)] < 0)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j, signed = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]
            if not signed:
                i = i - 1
            if signed:
                for col in range(j - 1, 0, -1):
                    row = i
                    while p.get((row + 1, col), a) < a:
                        row += 1
                    if row == col:
                        i, j = col, col
                        break
                    if (row + 1, col) in p and p[(row + 1, col)] == a:
                        a = p[(row, col)]
                    else:
                        a, p[(row, col)] = p[(row, col)], a
            for row in range(i, 0, -1):
                col = j
                while p.get((row, col + 1), a) < a:
                    col += 1
                if (row, col + 1) in p and p[(row, col + 1)] == a:
                    a = p[(row, col)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def random(cls, mu):
        mu = Partition(*mu.parts[:])
        nu = mu.transpose()
        t = Tableau.from_partition(mu)

        for n in range(abs(mu), 0, -1):
            m = random.randint(1, abs(mu))
            i = 1
            while m > mu(i) and i <= len(mu):
                m -= mu(i)
                i += 1
            j = m

            while True:
                row = mu(i) - j
                col = nu(j) - i

                if row == 0 and col == 0:
                    break
                m = random.randint(1, row + col)

                if m <= row:
                    j += m
                elif row < m:
                    i += m - row
            t.mapping[(i, j)] = MarkedNumber(n)
            mu.decrement(i)
            nu.decrement(j)
            if n % 1000 == 0:
                print(n, '. . .')
        return t

    @classmethod
    def inverse_eg(cls, p, q):
        ans = len(q) * [0]

        p = {k: v.number for k, v in p.mapping.items()}
        q = {k: v.number for k, v in q.mapping.items()}
        order = len(q) * [0]
        for (i, j) in q:
            order[q[(i, j)] - 1] = (i, j)

        for n in range(len(q), 0, -1):
            if n % 1000 == 0:
                print(n, '. . .')

            i, j = order[n - 1]
            a = p[(i, j)]
            del p[(i, j)]

            for row in range(i - 1, 0, -1):
                col = j
                while p.get((row, col + 1), a + 1) <= a:
                    col += 1
                if p[(row, col)] == a:
                    a = p[(row, col - 1)]
                else:
                    a, p[(row, col)] = p[(row, col)], a
            ans[n - 1] = a
        return tuple(ans)

    @classmethod
    def longest_eg_insertion_tableau(cls, n):
        return Tableau.from_string(';'.join([','.join([str(j) for j in range(i, n)]) for i in range(1, n)]))

    @classmethod
    def longest_inv_insertion_tableau(cls, n):
        partition = StrictPartition(*list(range(n - 1, 0, -2)))
        t = Tableau.from_partition(partition)
        for i in range(1, len(partition) + 1):
            for j in range(i, n - i + 1):
                t.mapping[(i, j)] = MarkedNumber(i + j - 1)
        return t

    @classmethod
    def longest_fpf_insertion_tableau(cls, n):
        partition = StrictPartition(*list(range(n - 2, 0, -2)))
        t = Tableau.from_partition(partition)
        for i in range(1, len(partition) + 1):
            for j in range(i, n - i):
                t.mapping[(i, j)] = MarkedNumber(i + j)
        return t

    @classmethod
    def random_sorting_network(cls, n):
        p = Tableau.longest_eg_insertion_tableau(n)
        q = cls.random(p.partition())
        return cls.inverse_eg(p, q)

    @classmethod
    def random_inv_network(cls, n):
        p = Tableau.longest_inv_insertion_tableau(n)
        q = cls.random_shifted(StrictPartition(p.partition()))
        return cls.inverse_inv(p, q)

    @classmethod
    def random_fpf_network(cls, n):
        p = Tableau.longest_fpf_insertion_tableau(n)
        q = cls.random_shifted(StrictPartition(p.partition()))
        return cls.inverse_fpf(p, q)

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
        if type(shape) in [tuple, list]:
            shape = Partition(*shape)
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
    def get_semistandard_shifted(cls, shape, n=None):
        if type(shape) == Partition:
            shape = StrictPartition(shape)
        if type(shape) == StrictPartition:
            shape = shape.shape

        if len(shape) == 0:
            return {Tableau()}

        n = len(shape) if n is None else n

        borders = {
            (a, b)
            for a in shape.horizontal_border_strips() | {Shape()}
            for b in (shape - a).vertical_border_strips(exclude_diagonal=True) | {Shape()}
            if len(a) > 0 or len(b) > 0
        }

        ans = set()
        for border_h, border_v in borders:
            for k in range(n):
                for t in cls.get_semistandard_shifted(shape - border_h - border_v, k):
                    mapping = t.mapping if t.mapping else {}
                    for i, j in border_h:
                        mapping[(i, j)] = MarkedNumber(k + 1)
                    for i, j in border_v:
                        mapping[(i, j)] = MarkedNumber(-k - 1)
                    ans.add(Tableau(mapping))
        return ans

    def __le__(self, other):
        assert type(other) == Tableau
        assert set(other.mapping) == set(self.mapping)
        return all(a <= other.mapping[x] for x, a in self.mapping.items())

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
        base = [['.' + (width - 1) * ' ' for i in range(self.max_column)] for j in range(self.max_row)]
        for i, j in self.mapping:
            v = str(self.mapping[(i, j)])
            base[i - 1][j - 1] = v + (width - len(v)) * ' '
        rows = [' '.join(row) for row in base]
        return '\n' + '\n'.join(reversed(rows)) + '\n'   # French
        #return '\n'.join(rows) + '\n'            # English

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
    def from_shifted_growth_diagram(cls, growth, edges, corners):
        def shdiff(nu, lam):
            return next(iter(Partition.skew(nu, lam, shifted=True)))

        p, q = Tableau(), Tableau()
        n, m = len(growth) - 1, len(growth[0]) - 1 if growth else 0

        for i in range(1, n + 1):
            mu, nu = growth[i][m], growth[i - 1][m]
            for a, b in Partition.skew(mu, nu, shifted=True):
                p = p.add(a, b, i)

        for i in range(1, m + 1):
            mu, nu = growth[n][i], growth[n][i - 1]
            v = -i if edges[n][i] else i
            j = corners[n][i]
            assert mu != nu
            a, b = shdiff(mu, nu)
            q = q.add(a, b, v)

        return p, q

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

    def remove(self, i, j):
        assert (i, j) in self
        return Tableau({a: b for a, b in self.mapping.items() if a != (i, j)})

    def row_reading_word(self):
        return tuple(
            self.mapping[key].number for key in sorted(self.mapping, key=lambda x: (-x[0], x[1]))
        )

    @classmethod
    def inverse_sagan_worley(cls, p, q):
        n = len(p)
        if n == 0:
            return ()
        if q.find(n):
            i, j = next(iter(q.find(n).mapping))
            a, p = p.entry(i, j), p.remove(i, j)
            cdir = False
        else:
            i, j = next(iter(q.find(-n).mapping))
            a, p = p.entry(i, j), p.remove(i, j)
            cdir = True
        while i > 1 or cdir:
            if cdir:
                j = j - 1
                i = max([k for (k, l) in p if l == j and p.entry(k, l) <= a])
                a, p = p.entry(i, j), p.set(i, j, a)
                cdir = (i != j)
            else:
                i = i - 1
                j = max([l for (k, l) in p if k == i and p.entry(k, l) < a])
                a, p = p.entry(i, j), p.set(i, j, a)
        return cls.inverse_sagan_worley(p, q) + (a.number,)

    @classmethod
    def inverse_rsk(cls, p, q):
        n = len(p)
        if n == 0:
            return ()
        i, j = next(iter(q.find(n).mapping))
        a = p.entry(i, j)
        p = p.remove(i, j)
        while i > 1:
            i = i - 1
            j = max([l for (k, l) in p if k == i and p.entry(k, l) < a])
            a, p = p.entry(i, j), p.set(i, j, a)
        return cls.inverse_rsk(p, q) + (a.number,)

    def rsk_insert(self, p, j=0):
        p = MarkedNumber(p) if type(p) == int else p
        if p is None:
            return (j, self)

        def rsk_bump(a, tup):
            if len(tup) == 0 or a >= tup[-1]:
                newtup = tup + (a,)
                q = None
            else:
                i = [j for j in range(len(tup)) if a < tup[j]][0]
                newtup = tup[:i] + (a,) + tup[i + 1:]
                q = tup[i]
            return q, newtup

        j += 1
        row = self.get_row(j)
        p, row = rsk_bump(p, row)
        tab = self.replace_row(j, row, shifted=False)
        return tab.rsk_insert(p, j)

    def hecke_insert(self, p, j=0):
        if p is None:
            return (j, self)

        def hecke_bump(a, tup):
            if len(tup) == 0 or a > tup[-1]:
                newtup = tup + (a,)
                q = None
            elif a == tup[-1]:
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

    def mixed_insert(self, p, j=1, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def sgn(x):
            return -1 if x.number < 0 else 1

        def bump(a, cdir, tup):
            for i, b in enumerate(tup):
                if abs(a) >= abs(b):
                    continue
                if not cdir and i == 0:
                    new = (a,) + tup[1:]
                    b = MarkedNumber((abs(b) - 1) * sgn(b))
                    cdir = True
                else:
                    new = tup[:i] + (MarkedNumber(abs(a) * sgn(b)),) + tup[i + 1:]
                    b = MarkedNumber(abs(b) * sgn(a))
                    cdir = abs(b) % 2 != 0
                return (b, cdir, new, i)
            return (None, cdir, tup + (a,), None)

        row, col = self.get_row(j), self.get_column(j)

        if column_dir:
            p, column_dir, col, i = bump(p, column_dir, col)
            tab = self.replace_column(j, col)

            if column_dir:
                j += 1
            else:
                j = i + 2
        else:
            p, column_dir, row, i = bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

            if column_dir:
                j = i + j + 1
            else:
                j += 1
        return tab.mixed_insert(p, j, column_dir, verbose=verbose)

    def sagan_worley_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def bump(a, cdir, tup):
            for i, b in enumerate(tup):
                if a > b:
                    continue
                if not cdir and a == b:
                    continue
                if not cdir and i == 0:
                    new = (a,) + tup[1:]
                    cdir = True
                else:
                    new = tup[:i] + (a,) + tup[i + 1:]
                return (b, cdir, new)
            return (None, cdir, tup + (a,))

        j += 1
        row, col = self.get_row(j), self.get_column(j)

        if column_dir:
            p, column_dir, col = bump(p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = bump(p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)
        return tab.sagan_worley_insert(p, j, column_dir, verbose=verbose)

    def involution_insert(self, p, j=0, column_dir=False, verbose=True):
        if p is None:
            return (j, column_dir, self)

        def involution_bump(j, a, cdir, tup):
            for i, b in enumerate(tup):
                if abs(a) > abs(b):
                    continue
                # print('j =', j, 'a =', a, 'cdir =', cdir, 'tup =', tup, 'i =', i + 1, 'b =', b)
                if a == b:
                    assert not a.is_primed()
                    b = MarkedNumber(abs(tup[i + 1]))
                    new = tup
                    cdir = cdir or (i == 0)

                elif a.is_primed() and a == -b:
                    b = -tup[i + 1]
                    assert b.is_primed()
                    new = tup
                elif b.is_primed() and a == -b and not (i == 0 and not cdir):
                    b = tup[i + 1]
                    assert not b.is_primed()
                    new = tup[:i] + (-tup[i], -tup[i + 1]) + tup[i + 2:]
                elif b.is_primed() and a == -b and i == 0 and not cdir:
                    b = tup[i + 1]
                    assert not b.is_primed()
                    new = tup
                    cdir = True

                elif not cdir and i == 0 and not b.is_primed() and not a.is_primed():
                    new = (a,) + tup[1:]
                    cdir = True
                elif not cdir and i == 0 and not b.is_primed() and a.is_primed():
                    b = -b
                    new = (-a,) + tup[1:]
                    cdir = True
                elif not cdir and i == 0 and b.is_primed() and not a.is_primed():
                    new = (-a,) + tup[1:]
                    cdir = True
                    b = -b
                elif not cdir and i == 0 and b.is_primed() and a.is_primed():
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
            p, column_dir, col = involution_bump(j, p, column_dir, col)
            tab = self.replace_column(j, col)
        else:
            p, column_dir, row = involution_bump(j, p, column_dir, row)
            tab = self.replace_row(j, row, shifted=True)

        if verbose:
            print(tab, '\n')
        assert tab.is_increasing()
        return tab.involution_insert(p, j, column_dir, verbose=verbose)

    def fpf_insert(self, p, j=0, column_dir=False, verbose=False):
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

    @cached_value(HORIZONTAL_STRIPS_CACHE)
    def _horizontal_strips(cls, mu, lam):  # noqa
        if not Partition._contains(mu, lam):
            return []

        core = [mu[i + 1] if i + 1 < len(mu) else 0 for i in range(len(mu))]
        for i in range(len(lam)):
            core[i] = max(core[i], lam[i])
        core = tuple(core)

        ans = []
        level = {core}
        while level:
            for nu in level:
                diff = {(i + 1, j + 1) for i in range(len(mu)) for j in range(nu[i], mu[i])}
                nu = nu if nu and nu[-1] > 0 else nu[:-1]
                corners = [(i + 1, nu[i]) for i in range(len(nu)) if core[i] < nu[i]]
                ans.append((nu, diff, corners))
            level = {
                nu[:i] + (nu[i] + 1,) + nu[i + 1:]
                for i in range(len(mu))
                for nu in level
                if nu[i] < mu[i]
            }
        return ans

    @classmethod
    def _subsets(cls, diff, corners, setvalued):
        if setvalued:
            for v in range(2**len(corners)):
                thisdiff = diff
                for i in range(len(corners)):
                    thisdiff = thisdiff if v % 2 == 0 else thisdiff | {corners[i]}
                    v = v // 2
                yield thisdiff
        else:
            yield diff

    @classmethod
    def semistandard(cls, max_entry, mu, nu=(), setvalued=False):  # noqa
        return cls._semistandard(max_entry, mu, nu, setvalued)

    @cached_value(SEMISTANDARD_CACHE)
    def _semistandard(cls, max_entry, mu, lam, setvalued):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam) and max_entry > 0:
            for nu, diff, corners in cls._horizontal_strips(mu, lam):
                for aug in cls._subsets(diff, corners, setvalued):
                    for tab in cls._semistandard(max_entry - 1, nu, lam, setvalued):
                        for (i, j) in aug:
                            tab = tab.add(i, j, max_entry)
                        ans.add(tab)
        return ans

    @classmethod
    def standard(cls, mu, nu=()):  # noqa
        return cls._standard(mu, nu)

    @cached_value(STANDARD_CACHE)
    def _standard(cls, mu, lam):  # noqa
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), mu[i]
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_partition(nu):
                    for tab in cls._standard(nu, lam):
                        ans.add(tab.add(row, col, n))
        return ans

    @classmethod
    def standard_shifted_marked(cls, mu, nu=(), diagonal_primes=False):  # noqa
        return cls._standard_shifted_marked(mu, nu, diagonal_primes)

    @cached_value(STANDARD_SHIFTED_MARKED_CACHE)
    def _standard_shifted_marked(cls, mu, lam, diagonal_primes):  # noqa
        assert Partition.is_strict_partition(mu)
        ans = set()
        if mu == lam:
            ans = {Tableau()}
        elif Partition._contains(mu, lam):
            n = sum(mu) - sum(lam)
            for i in range(len(mu)):
                row, col = (i + 1), (i + mu[i])
                nu = list(mu)
                nu[i] -= 1
                nu = Partition.trim(nu)
                if Partition.is_strict_partition(nu):
                    for tab in cls._standard_shifted_marked(nu, lam, diagonal_primes):
                        ans.add(tab.add(row, col, n))
                        if diagonal_primes or row != col:
                            ans.add(tab.add(row, col, -n))
        return ans
